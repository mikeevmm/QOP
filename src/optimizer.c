#include "include/optimizer.h"

AdadeltaSettings optimizer_adadelta_get_default() {
  AdadeltaSettings ada_default;
  ada_default.rho = 0.95;
  ada_default.epsilon = 1.0e-6;
  return ada_default;
}

Result optimizer_gate_param_init(GateParameterization *gate_param, Gate *gate,
                                 unsigned int param_count, double *params,
                                 double *deltas) {
  if (gate == NULL) return result_get_invalid_reason("gate is null");
  
  {
    void *mem = malloc(sizeof(double) * param_count);
    if (mem == NULL) return result_get_invalid_reason("could not malloc");
    memcpy(mem, params, param_count * sizeof(double));
    gate_param->params = (double *)mem;
  }
  
  {
    void *mem = malloc(sizeof(double) * param_count);
    if (mem == NULL) return result_get_invalid_reason("could not malloc");
    memcpy(mem, deltas, param_count * sizeof(double));
    gate_param->deltas = (double *)mem;
  }

  gate_param->gate = gate;
  gate_param->param_count = param_count;

  return result_get_valid_with_data(gate_param);
}

void optimizer_gate_param_free(GateParameterization *gate_param)
{
  free(gate_param->params);
  free(gate_param->deltas);
}

Result optimizer_settings_init(OptimizerSettings *opt_settings,
                               Circuit *circuit, double _Complex *hamiltonian,
                               double stop_at,
                               GateParameterization *parameterizations,
                               unsigned int parameterizations_count) {
  if (opt_settings == NULL)
    return result_get_invalid_reason("given OptmizerSettings* is null");
  if (circuit == NULL)
    return result_get_invalid_reason("given Circuit* is null");
  if (hamiltonian == NULL)
    return result_get_invalid_reason("given hamiltonian* is null");
  if (parameterizations == NULL)
    return result_get_invalid_reason("given parameterizations is null");

  double _Complex *zero_state;
  {
    unsigned int state_size = 1U << circuit->depth[0];
    void *mem = malloc(sizeof(double _Complex) * state_size);
    if (mem == NULL) {
      return result_get_invalid_reason("could not malloc");
    }
    zero_state = mem;
    memset(zero_state, 0, state_size * sizeof(double _Complex));
    zero_state[0] = (double _Complex)1.;
  }

  opt_settings->circuit = circuit;
  opt_settings->hamiltonian = hamiltonian;
  opt_settings->stop_at = stop_at;
  opt_settings->parameterizations = parameterizations;
  opt_settings->parameterization_count = parameterizations_count;
  opt_settings->zero_state = zero_state;

  return result_get_valid_with_data(opt_settings);
}

void optimizer_settings_free(OptimizerSettings *opt_settings) {
  free(opt_settings->zero_state);
}

Result optimizer_init(Optimizer *optimizer, OptimizerSettings opt_settings,
                      AdadeltaSettings ada_settings) {
  optimizer->opt_settings = opt_settings;
  optimizer->ada_settings = ada_settings;
  return result_get_empty_valid();
}

Result optimizer_optimize(Optimizer *optimizer) {
  const AdadeltaSettings ada_settings = optimizer->ada_settings;
  OptimizerSettings opt_settings = optimizer->opt_settings;
  Circuit *circuit = opt_settings.circuit;
  const unsigned int state_size = 1U << circuit->depth[0];

  if (optimizer == NULL) {
    return result_get_invalid_reason("optimizer is null");
  }

  // Count the total number of parameters
  //   This could be done (maybe) at the same time as another iteration
  //   through `opt_settings.parameterization`, but since it's only ran
  //   once (and not in the optimization cycle), the impact should be
  //   minimal.
  unsigned int param_count = 0;
  {
    Iter params_iter = iter_create(opt_settings.parameterizations,
                                   sizeof(GateParameterization),
                                   opt_settings.parameterization_count);
    Option next;
    while ((next = iter_next(&params_iter)).some) {
      param_count += ((GateParameterization *)next.data)->param_count;
    }
  }

  //  Gradient squared accumulation variable
  double acc_grad = 0;
  //  Parameter update accumulation variable
  double acc_dxsqr = 0;
  //  Tolerance to maximum component of update to parameters
  double max_grad = optimizer->opt_settings.stop_at + 1;

  // Buffers for simulation results at left and right
  double _Complex left_buffer[state_size];
  double _Complex right_buffer[state_size];

  // Optimization loop
  while (max_grad > opt_settings.stop_at) {
    // Offset all parameters to the left + reparameterize gates
    {
      Iter params_iter = iter_create(opt_settings.parameterizations,
                                     sizeof(GateParameterization),
                                     opt_settings.parameterization_count);
      Option next;
      while ((next = iter_next(&params_iter)).some) {
        GateParameterization *gate_param = (GateParameterization *)next.data;
        Gate *gate = gate_param->gate;

        for (unsigned int i = 0; i < gate_param->param_count; ++i) {
          *(gate_param->params + i) -= *(gate_param->deltas + i);
        }

        gate->reparamFn(&gate->matrix, gate_param->params);
      }
    }

    // Run simulation @ left
    {
      memcpy(left_buffer, opt_settings.zero_state, state_size);
      circuit_run(circuit, &left_buffer);
    }

    // Offset all parameters to the right
    // (de-offsetting the previous delta)
    {
      Iter params_iter = iter_create(opt_settings.parameterizations,
                                     sizeof(GateParameterization),
                                     opt_settings.parameterization_count);
      Option next;
      while ((next = iter_next(&params_iter)).some) {
        GateParameterization *gate_param = (GateParameterization *)next.data;
        Gate *gate = gate_param->gate;

        for (unsigned int i = 0; i < gate_param->param_count; ++i) {
          *(gate_param->params + i) += 2 * (*(gate_param->deltas + i));
        }

        gate->reparamFn(&gate->matrix, gate_param->params);
      }
    }

    // Run simulation @ right
    {
      memcpy(right_buffer, opt_settings.zero_state, state_size);
      circuit_run(circuit, &right_buffer);
    }

    // Calculate gradient
    double param_gradient[param_count];
    {
      double coef = 0.;
      const unsigned int state_size = 1U << opt_settings.circuit->depth[0];
      double _Complex *hamiltonian = opt_settings.hamiltonian;

      // Calculate the common (double sum) coefficient
      for (unsigned int i = 0; i < state_size; ++i) {
        double _Complex left_i = left_buffer[i];
        double _Complex right_i = right_buffer[i];

        for (unsigned int j = i + 1; j < state_size; ++j) {
          double _Complex hamilt_ij = *(hamiltonian + i * state_size + j);
          double _Complex left_j = left_buffer[j];
          double _Complex right_j = right_buffer[j];

          coef += 2. * creal(hamilt_ij *
                             (conj(right_i) * right_j - conj(left_i) * left_j));
        }

        double _Complex hamilt_ii = *(hamiltonian + i * state_size + i);
        coef += creal(hamilt_ii * (right_i + left_i) * conj(right_i - left_i));
      }

      coef /= 2.;

      // Set values of parameter gradient +
      // accumulate gradient
      acc_grad = ada_settings.rho * acc_grad;
      unsigned int param_offset = 0;

      for (unsigned int i = 0; i < param_count; ++i) {
        // Calculate ith element of param gradient
        GateParameterization *relevant_param =
            opt_settings.parameterizations + param_offset;
        double grad = coef / *((relevant_param)->deltas + i - param_offset);

        // Set values of parameter gradient
        param_gradient[i] = grad;

        // Accumulate gradient
        acc_grad += (1 - ada_settings.rho) * pow(grad, 2);

        // End of this parameter set?
        if (i - param_offset == relevant_param->param_count) {
          param_offset += 1;
        }
      }
    }

    // Update parameters
    // + find largest update
    // + accumulate updates
    {
      max_grad = 0;
      acc_dxsqr = ada_settings.rho * acc_dxsqr;

      unsigned int param_arr_index_offset = 0;
      Iter params_iter = iter_create(opt_settings.parameterizations,
                                     sizeof(GateParameterization),
                                     opt_settings.parameterization_count);
      Option next;
      while ((next = iter_next(&params_iter)).some) {
        GateParameterization *gate_param = (GateParameterization *)next.data;

        for (unsigned int i = 0; i < gate_param->param_count; ++i) {
          // The update to apply to the parameters
          double update = -sqrt(acc_dxsqr + ada_settings.epsilon) /
                          sqrt(acc_grad + ada_settings.epsilon) *
                          param_gradient[param_arr_index_offset + i];

          // Remember to de-offset right shift
          *(gate_param->params + i) += -(*(gate_param->deltas + i)) + update;

          // Save largest update
          if (max_grad < fabs(update)) max_grad = fabs(update);

          // Accumulate updates
          acc_dxsqr += (1. - ada_settings.rho) * pow(update, 2);
        }

        param_arr_index_offset += gate_param->param_count;
      }
    }

    // Perform updates according to ADADELTA algorithm
    acc_grad = ada_settings.rho * acc_grad;
  }

  // Update gates to match final parameters
  {
    Iter params_iter = iter_create(opt_settings.parameterizations,
                                   sizeof(GateParameterization),
                                   opt_settings.parameterization_count);
    Option next;
    while ((next = iter_next(&params_iter)).some) {
      GateParameterization *gate_param = (GateParameterization *)next.data;
      Gate *gate = gate_param->gate;
      gate->reparamFn(&gate->matrix, gate_param->params);
    }
  }

  return result_get_empty_valid();
}