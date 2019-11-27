#include "include/optimizer.h"

// Returns a default `AdadeltaSettings` struct.
// The default values here passed are the same used in section 4.1. of
// the ADADELTA paper (arXiv:1212.5701).
AdadeltaSettings optimizer_adadelta_get_default() {
  AdadeltaSettings ada_default;
  ada_default.rho = 0.95;
  ada_default.epsilon = 1.0e-6;
  return ada_default;
}
// Returns a default `LbfgsSettings` struct.
// The number 3 was chosen based on the remark by Nocedal (2006) that
// > Practical experience has shown that modest values of m
// > (between 3 and 20, say) often produce satisfactory results.
// The alpha of 1 is also based on the remark that
// > as a result [of the estimation of H] the step length alpha = 1 is
// accepted in most iterations
LbfgsSettings optimizer_lbfgs_get_default() {
  LbfgsSettings lbfgs_settings;
  lbfgs_settings.m = 3;
  lbfgs_settings.alpha = 1;
  return lbfgs_settings;
}

// TODO: Comment
Result optimizer_algo_settings_init(OptimizerAlgoSettings *algo_settings,
                                    OptimizerAlgorithm algorithm,
                                    void *specific_settings) {
  algo_settings->algorithm = algorithm;
  switch (algorithm) {
    case AlgoAdadelta: {
      memcpy(&algo_settings->specific.ada_settings, specific_settings,
             sizeof(AdadeltaSettings));
    } break;
    case AlgoLbfgs: {
      memcpy(&algo_settings->specific.lbfgs_settings, specific_settings,
             sizeof(LbfgsSettings));
    } break;
    default:
      return result_get_invalid_reason("unknown algorithm");
      break;
  }
  return result_get_empty_valid();
}

// Initializes a new `GateParameterization` object.
// There are heap memory allocations to copy `params` and `deltas` into
// which should later be freed with `optimizer_gate_param_free`. That
// call is the user's responsability.
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

// Frees the internal heap memory allocations made at initialization.
void optimizer_gate_param_free(GateParameterization *gate_param) {
  free(gate_param->params);
  free(gate_param->deltas);
}

// Initializes an `OptimizerSettings` object.
// There is an internal malloc to store a precomputed |000...> state.
// This should be later freed by the user using `optmizer_settings_free`.
// Otherwise, only that `*hamiltonian` is a row major flattened square
// matrix is of note, but the user is alerted in the documentation to
// make the appropriate cast. There is, unfortunately, no way I could
// find to sanity check it, although gcc does produce an invalid cast
// warning.
Result optimizer_settings_init(OptimizerSettings *opt_settings,
                               Circuit *circuit, double _Complex *hamiltonian,
                               double stop_at,
                               GateParameterization *parameterizations,
                               unsigned int parameterizations_count,
                               Option_Uint max_iterations) {
  if (opt_settings == NULL)
    return result_get_invalid_reason("given OptmizerSettings* is null");
  if (circuit == NULL)
    return result_get_invalid_reason("given Circuit* is null");
  if (hamiltonian == NULL)
    return result_get_invalid_reason("given hamiltonian* is null");
  if (parameterizations == NULL)
    return result_get_invalid_reason("given reparams is null");

  unsigned int state_size = 1U << circuit->depth[0];

  // Create |000...> state
  // malloc some memory, compute the state into that and save the pointer.
  double _Complex *zero_state;
  {
    void *mem = malloc(sizeof(double _Complex) * state_size);
    if (mem == NULL) {
      return result_get_invalid_reason("could not malloc");
    }
    zero_state = mem;
    memset(zero_state, 0, state_size * sizeof(double _Complex));
    zero_state[0] = 1.;
  }

  // Create a compacted representation of the hamiltonian matrix.
  // This should help with sparse matrices, which are common as physical
  // Hamiltonians
  Vector ham_rows;
  {
    Result init_r = vector_init(&ham_rows, sizeof(Vector), 0);
    if (!init_r.valid) {
      free(zero_state);
      return init_r;
    }
  }
  for (unsigned int u = 0; u < state_size; ++u) {
    Vector row;
    {
      Result init_r = vector_init(&row, sizeof(OptimizerDCPackedRowElem), 0);
      if (!init_r.valid) {
        vector_free(&ham_rows);
        free(zero_state);
        return init_r;
      }
    }

    for (unsigned int k = u; k < state_size; ++k) {
      // (double _Complex *)([][]) is row major
      double _Complex ham_uk = hamiltonian[u * state_size + k];
      if (ham_uk != 0) {
        OptimizerDCPackedRowElem elem;
        elem.value = ham_uk;
        elem.j = k;
        vector_push(&row, &elem);
      }
    }

    vector_push(&ham_rows, &row);
  }

  // Initialize fields of the new OptimizerSettings.
  opt_settings->circuit = circuit;
  opt_settings->hamiltonian = ham_rows;
  opt_settings->stop_at = stop_at;
  opt_settings->reparams = parameterizations;
  opt_settings->reparams_count = parameterizations_count;
  opt_settings->zero_state = zero_state;
  opt_settings->max_iterations = max_iterations;

  return result_get_valid_with_data(opt_settings);
}

// Frees the heap memory internally allocated during
// `optimizer_settings_init`.
void optimizer_settings_free(OptimizerSettings *opt_settings) {
  // Free the compact hamiltonian representation
  {
    Iter row_iter = vector_iter_create(&opt_settings->hamiltonian);
    Option next_row;
    while ((next_row = iter_next(&row_iter)).some) {
      Vector *row = (Vector *)next_row.data;
      vector_free(row);
    }
    iter_free(&row_iter);
    vector_free(&opt_settings->hamiltonian);
  }
  // Free |000>
  free(opt_settings->zero_state);
}

// Initializes a new `Optimizer` object.
// This is just a glorified (`OptimizerSettings`, `OptimizerAlgoSettings`)
// pair.
Result optimizer_init(Optimizer *optimizer, OptimizerSettings opt_settings,
                      OptimizerAlgoSettings algo_settings) {
  optimizer->opt_settings = opt_settings;
  optimizer->algo_settings = algo_settings;
  return result_get_empty_valid();
}

Result optimization_result_as_result(OptimizationResult opt_result) {
  Result cast;
  cast.valid = opt_result.valid;
  cast.content = opt_result.content;
  return cast;
}

OptimizationResult result_as_optimization_result(Result result) {
  OptimizationResult cast;
  cast.content = result.content;
  cast.valid = result.valid;
  return cast;
}

// Performs the actual optimization of the parameters pointed to in the
// `OptimizationSettings` object.
// The gradient calculation routine is as follows:
//
// 1. Offset each parameter $\theta_i$ by $+\Delta_i$, making the
//   appropriate reparameterizations.
// 2. Compute the output of the circuit with each new parameter (preserving
//   the others as originals), when given |000...> as input.
//   Store the output in `buffer_right`.
// 3. Offset each (original) parameter $\theta_i$ by $-\Delta_i$, making
//   the appropriate reparameterizations. Store the output in
//   `buffer_left`.
// 4. Compute the expectation value gradient, given these two states.
//
// The expectation value gradient is calculated based on an algebraically
// simplified expression, whose deduction follows:
//
//  Let $\vec{\Delta}: \vec{\Delta}_i = \Delta_i$,
//  Let $\ket{\phi}$ be the output state given some parameters,
//  Let $\Theta: \Theta_i = \theta_i$,
//  Let the hamiltonian be $H$, and the action of the quantum circuit be
//   the unitary matrix $U$,
//  Let $C$ be the cost function, which is identical to the expectation
//   value of the Hamiltonian
//  Let $\hat{e}_p$ denote a basis versor,
//  Let $q$ be the number of qubits in the circuit,
//  Following the implicit sum convention (sometimes made explicit),
//
// Then:
//  $$ \nabla_\Theta C = \hat{e}_p \pdv{\theta_p}
//                                         \qty( \expval{H}{\phi} ) = $$
//  $$                 = \hat{e}_p \pdv{\theta_p}
//                                   \qty( \expval{U^\dag H U}{0} ) = $$
//  $$ = \hat{e}_p \pdv{\theta_p} \qty( \delta_{i 1} U^\dag_{i u}
//                                        H_{uk} U_{kj} \delta_{j1} ) $$
//  $$ = \hat{e}_p \pdv{\theta_p} \qrt( U_{u1}^* H_{uk} U_{k1} )      $$
//  $$ = \hat{e}_p \qty( H_{uk} \pdv{U_{u1}^*}{\theta_p} U_{k1} +
//                           H_{uk} U_{u1}^* \pdv{U_{k1}}{\theta_p} ) $$
//  $$ = \hat{e}_p H_{uk} \qty(
//         \sum_{u=1}^{2^q} \sum_{k=1}^{2^q}
//                               \qty(U_{k1} \pdv{U_{u1}^*}{\theta_p}) +
//         \sum_{u=1}^{2^q} \sum_{k=1}^{2^q}
//                           \qty(U_{u1} \pdv{U_{k1}^*}{\theta_p})^*) $$
//  $$ = \hat{e}_p \qty(H_{uk} U_{k1} \pdv{U_{u1}^*}{\theta_p} +
//                  + \qty(H_{uk} U_{k1} \pdv{U_{u1}^*}{\theta_p})^*) $$
//  $$ = \hat{e}_p 2 \Re{ H_{uk} U_{k1} \pdv{U_{u1}^*}{\theta_p}}     $$
//
// Letting $A^p$ and $B^p$ be the results to the right and left of
// parameter $p$, respectively, we now approximate
//        $U_{i1} \approx (A_i^p + B_i^p)/2$ and
//        $\pdv{U_{u1}}{\theta_p} \approx  (A_u^p - B_u^p)/(2 \Delta_i)$
// and also note that, for a function f : (real --> complex),
//            $\pdv{f^*}{x} = (\pdv{f}{x})^*$
// so that
//
//  $$ = \hat{e}_p 2 \Re{ H_{uk} \qty(\frac{\vec{A} + \vec{B}}{2})_k
//              \times \qty(\frac{\vec{A} - \vec{B}}{2 \Delta})_u^* } $$
//  $$ = \Delta^{-1} \odot \hat{e}_p/2 \cdot
//               \Re{ H_{uk} (A_k^p + B_k^p) ({A_u^p}^* - {B_u^p}^*)} $$
//
// Where we defined $ \Delta^{-1}: \Delta^{-1}_k = 1/\Delta_k$
// (Hereafter suppressing the $p$ supraindex in $A$ and $B$)
//
//  $$ = \Delta^{-1} \odot \hat{e}_p/2 \Re{ \sum_{u=1}^{2^q} \qty[
//          \sum{k=1}^{u-1} \qty(H_{uk} (A_k + B_k)(A_u^* - B_u^*)) +
//          \sum{k=u+1}^{2^q} \qty(H_{uk} (A_k + B_k) (A_u^* - B_u^*)) +
//          H_{uu} (A_u + B_u) (A_u^* - B_u^*)]}                      $$
//  $$ = \Delta^{-1} \odot \hat{e}_p/2 \cdot \Re{
//        \sum_{k=1}^{2^q} \sum{u>k}{2^q} \qty(
//                                 H_{uk} (A_k + B_k)(A_u^* - B_u^*) ) +
//        \sum_{k=1}^{2^q} \sum{u>k}{2^q} \qty(
//                                 H_{uk} (A_k + B_k)(A_u^* - B_u^*) ) +
//        \sum{u=1}^{2^q} \qty( H_{uu} (A_u + B_u)(A_u^* - B_u^*) )}  $$
//  $$ = \Delta^{-1} \odot \hat{e}_p/2 \cdot \Re{
//        \sum_{k=1}^{2^q} \qty[ \sum{u>k}{2^q} \qty(
//              H_{ku} (A_u + B_u} (A_k^* - B_k^*) +
//              H_{uk} (A_k + B_k} (A_u^* - B_u^*)) +
//        H_{uu} (A_u + B_u)(A_u^* - B_u^*) ]}                        $$
//  $$ = \Delta^{-1} \odot \hat{e}_p/2 \cdot
//        \Re{ \sum_u \qty[ \sum_{k>u} \qty(
//          (H_{uk}^* A_u A_k^* + (H_{uk}^* A_u A_k^*)^*) -
//          (H_{uk}^* B_u B_k^* + (H_{uk}^* B_u B_k^*)^*) +
//          (H_{uk}^* B_u A_k^* - (H_{uk}^* B_u A_k^*)^*) -
//          (H_{uk}^* A_u B_k^* - (H_{uk}^* A_u B_k^*)^*) +
//          (H_{uu} (A_u + B_u)(A_u^* - B_u^*)) ]}                    $$
//  $$ = \Delta^{-1} \odot \hat{e}_p/2  \sum_u \qty[ \sum_{k>u} \qty(
//          2 \Re{H_{uk} A_u^* A_k} - 2 \Re{H_{uk} B_u^* B_k} +
//          2i \Im{H_{uk}^* B_u A_k^*} - 2i \Im{H_{uk}^* A_u B_k^*} ) +
//       + H_{uu} (A_u + B_u)(A_u^* - B_u^*) ]                        $$
//
// Yielding, finally,
//
// $$ \nabla_\Theta C = \Delta^{-1} \odot \hat{e}_p/2
//  \sum_{u=1}^{2^q} \qty[
//  \sum_{k = u+1}^{2^q} \qty( 2 \Re{H_{uk} ({A_u^p}^* {A_k^p} - {B_u^p}^*
//  {B_k^p})}) + \Re{ H_{uu} ({A_u^p} + {B_u^p})({A_u^p}^* - {B_u^p}^*) }] $$
OptimizationResult optimizer_optimize(Optimizer *optimizer,
                                      OptimizerParamCallback param_callback,
                                      OptimizerEnergyCallback energy_callback,
                                      void *callback_context) {
  switch (optimizer->algo_settings.algorithm) {
    case AlgoAdadelta:
      return optimizer_optimize_adadelta(optimizer, param_callback,
                                         energy_callback, callback_context);
      break;
    case AlgoLbfgs:
      return optimizer_optimize_lbfgs(optimizer, param_callback,
                                      energy_callback, callback_context);
      break;
    default:
      return result_as_optimization_result(
          result_get_invalid_reason("Unknown algorithm."));
      break;
  }
}

// Performs a circuit optimization using the ADADELTA algorithm, very
// neatly presented in the paper and copied here for convenience:
//
//  =====================================================
//  **Algorithm 1** Computing ADADELTA update at time $t$
//  -----------------------------------------------------
//  **Require**: Decay rate $\rho$, Constant $\epsilon$
//  **Require**: Initial parameter $x_1$
//  1: Initialize accumulation variables $E[g^2]_0 = 0$,
//     $E[\Delta x^2]_0 = 0$
//  2: **for** $t=1$ **do** %% Loop over # of updates
//  3:    Compute Gradient: $g_t$
//  4:    Accumulate Gradient:
//       $E[g^2]_t = \rho E[g^2]_{t-1} + (1 - \rho) g_t^2
//  5:    Compute Update:
//       $\Delta x_t = -\frac{RMS[\Delta x]_{t-1}}{RMS[g]_t} g_t$
//  6:    Accumulate Updates:
//       $E[\Delta x^2]_t = \rho E[\Delta x^2]_{t-1} +
//                                       (1 - \rho) \Delta x_t^2$
//  7:    Apply Update: $x_{t+1} = x_t + \Delta x_t$
//  =====================================================
//
// In this particular case, the "position" vector is considered to be
// the vector of all concatenated parameters, and the computation of the
// gradient is based on a simple finite difference approximation, with
// further algebraic optimization (see below). Further, it is also of
// notice that some of the above steps were compacted into a single
// pass over the gates (since the above algorithm is phrased "vectorially",
// while the programmatic approach is component-wise). Such compactions
// are commented appropriately inline.
OptimizationResult optimizer_optimize_adadelta(
    Optimizer *optimizer, OptimizerParamCallback param_callback,
    OptimizerEnergyCallback energy_callback, void *callback_context) {
  const AdadeltaSettings ada_settings =
      optimizer->algo_settings.specific.ada_settings;
  OptimizerSettings opt_settings = optimizer->opt_settings;
  Circuit *circuit = opt_settings.circuit;
  const unsigned int state_size = 1U << circuit->depth[0];

  if (optimizer == NULL) {
    Result invalid_result = result_get_invalid_reason("optimizer is null");
    return result_as_optimization_result(invalid_result);
  }

  // Count the number of parameters
  // It isn't great to have a whole cycle just for this, but it's just
  //  once in the optimization cycle...
  unsigned int abs_param_count = 0;
  {
    unsigned int reparams_count = opt_settings.reparams_count;
    for (unsigned int reparam_index = 0; reparam_index < reparams_count;
         ++reparam_index) {
      abs_param_count += opt_settings.reparams[reparam_index].param_count;
    }
  }

  // (Absolute) Max known component of gradient
  // used to terminate optimization cycle
  double max_grad = opt_settings.stop_at + 1;

  // Adadelta accumulation variables:
  double grad_sqr_acc = 0;
  double dx_sqr_acc = 0;

  // Iteration count; this will only be useful if there is a max
  // iteration count allowed
  unsigned long int iter_count = 0;

  // Optimization cycle:
  while (max_grad > opt_settings.stop_at) {
    // Max iteration?
    if (opt_settings.max_iterations.some) {
      iter_count++;
      if (iter_count >= opt_settings.max_iterations.data) break;
    }

    // Compute an adadelta update

    // Initialize the gradient array
    // The parameters here are "flattened"
    double param_gradient[abs_param_count];
    memset(param_gradient, 0, sizeof(double) * abs_param_count);

    // Buffers to write the result of the simulation at left and right
    // We don't need to initialize these right now because they're
    //  always initialized to |0> just before simulation
    double _Complex buff_left[state_size], buff_right[state_size];

    // First loop, calculate gradient and gradient accumulation
    {
      // Gradient accumulation will be done component wise, so we factor
      //  out the first term pre-component loop
      grad_sqr_acc = ada_settings.rho * grad_sqr_acc;

      // Flat index, to know what index of the gradient array to access
      unsigned int flat_param_index = 0;

      // Calculate also max grad
      max_grad = 0;

      // Loop through all parameters...
      for (unsigned int reparam_index = 0;
           reparam_index < opt_settings.reparams_count; ++reparam_index) {
        GateParameterization *reparam = opt_settings.reparams + reparam_index;
        Gate *reparam_gate_ptr = reparam->gate;

        for (unsigned int subparam_index = 0;
             subparam_index < reparam->param_count; ++subparam_index) {
          double *param_ptr = reparam->params + subparam_index;
          double delta = reparam->deltas[subparam_index];

          // Left shift
          {
            *param_ptr -= delta;
            reparam_gate_ptr->reparamFn(&reparam_gate_ptr->matrix,
                                        reparam->params);

            // Initialize left buffer to |0>
            memcpy(buff_left, opt_settings.zero_state,
                   sizeof(double _Complex) * state_size);

            // Simulate into left buffer
            {
              Result run_r = circuit_run(circuit, &buff_left);
              if (!run_r.valid) {
                return result_as_optimization_result(run_r);
              }
            }

            // Undo shift
            *param_ptr += delta;
          }

          // Right shift
          {
            *param_ptr += delta;
            reparam_gate_ptr->reparamFn(&reparam_gate_ptr->matrix,
                                        reparam->params);

            // Initialize right buffer to |0>
            memcpy(buff_right, opt_settings.zero_state,
                   sizeof(double _Complex) * state_size);

            // Simulate into right buffer
            {
              Result run_r = circuit_run(circuit, &buff_right);
              if (!run_r.valid) {
                return result_as_optimization_result(run_r);
              }
            }

            // Undo right shift; reparameterize
            *param_ptr -= delta;
            reparam_gate_ptr->reparamFn(&reparam_gate_ptr->matrix,
                                        reparam->params);
          }

          // Calculate this parameter's gradient component via the method
          //  presented above.
          double grad_component = 0;
          {
            for (unsigned int u = 0; u < state_size; ++u) {
              Vector *row = ((Vector *)opt_settings.hamiltonian.data) + u;
              double _Complex left_u = buff_left[u];
              double _Complex right_u = buff_right[u];

              for (unsigned int i = 0; i < row->size; ++i) {
                OptimizerDCPackedRowElem elem =
                    *((OptimizerDCPackedRowElem *)(row->data) + i);
                unsigned int k = elem.j;
                double _Complex left_k = buff_left[k];
                double _Complex right_k = buff_right[k];
                double _Complex ham_uk = elem.value;

                if (u == k) {
                  grad_component += creal(ham_uk * (right_u + left_u) *
                                          conj(right_u - left_u)) /
                                    2.;
                } else {
                  grad_component += creal(ham_uk * (conj(right_u) * right_k -
                                                    conj(left_u) * left_k));
                }
              }
            }

            grad_component /= delta;
          }

          // Accumulate the gradient (component wise)
          grad_sqr_acc += pow(grad_component, 2) * (1. - ada_settings.rho);

          // Save the gradient
          param_gradient[flat_param_index] = grad_component;
          ++flat_param_index;

          // Optionally report new parameter
          if (param_callback != NULL)
            param_callback(flat_param_index, *param_ptr, callback_context);

          // Max grad?
          if (fabs(grad_component) > max_grad) max_grad = fabs(grad_component);

          // We need a full loop before applying the update, as the update
          // depends on the complete value of grad_sqr_acc
        }  // Finished subparam
      }    // Finished reparam
    }      // Finished first loop

    // Store the dx_sqr_acc of the last full iteration loop, as
    // this is needed to compute the update of all the variables;
    // this allows to update the actual dx_sqr_acc as we update
    // the variables, without affecting the update calculation
    double old_dx_sqr_acc = dx_sqr_acc;

    // Second parameter loop; now calculating & applying the update
    {
      // Same as gradient accumulation; factor out non-component term
      dx_sqr_acc = ada_settings.rho * dx_sqr_acc;

      // Flat index, to know what index of the gradient array to access
      unsigned int flat_param_index = 0;

      // Loop over parameters...
      for (unsigned int reparam_index = 0;
           reparam_index < opt_settings.reparams_count; ++reparam_index) {
        GateParameterization reparam = opt_settings.reparams[reparam_index];

        for (unsigned int subparam_index = 0;
             subparam_index < reparam.param_count; ++subparam_index) {
          double *param_ptr = reparam.params + subparam_index;

          // Calculate update
          double update = -sqrt(old_dx_sqr_acc + ada_settings.epsilon) /
                          sqrt(grad_sqr_acc + ada_settings.epsilon) *
                          param_gradient[flat_param_index];

          // Accumulate dx component
          dx_sqr_acc += (1. - ada_settings.rho) * pow(update, 2);

          // Apply update!
          *param_ptr += update;

          // Flat index...
          ++flat_param_index;
        }  // Finished with subparam
      }    // Finished with reparams

      // Optionally calculate and report energy
      if (energy_callback != NULL) {
        double _Complex energy = 0;
        double _Complex current_phi_state[state_size];
        memcpy(current_phi_state, opt_settings.zero_state,
               sizeof(double _Complex) * (unsigned long int)state_size);
        circuit_run(circuit, &current_phi_state);

        for (unsigned int i = 0; i < opt_settings.hamiltonian.size; ++i) {
          Vector *row = (Vector *)opt_settings.hamiltonian.data + i;
          for (unsigned int u = 0; u < row->size; ++u) {
            OptimizerDCPackedRowElem elem =
                *((OptimizerDCPackedRowElem *)row->data + u);
            unsigned int j = elem.j;
            double _Complex value = elem.value;
            if (i == j)
              energy +=
                  conj(current_phi_state[i]) * value * current_phi_state[j];
            else
              energy +=
                  2 * conj(current_phi_state[i]) * value * current_phi_state[j];
          }
        }

        energy_callback(energy, callback_context);
      }
    }
  }

  // Done with optimization; reparameterize all gates so that they match
  // last calculated parameters
  for (unsigned int reparam_index = 0;
       reparam_index < opt_settings.reparams_count; ++reparam_index) {
    GateParameterization *param = opt_settings.reparams + reparam_index;
    param->gate->reparamFn(&param->gate->matrix, param->params);
  }

  OptimizationResult result;
  if (opt_settings.max_iterations.some &&
      iter_count >= opt_settings.max_iterations.data)
    result.quit_on_max_iter = true;
  else
    result.quit_on_max_iter = false;
  result.valid = true;
  result.content.data = circuit;
  return result;
}

// Helper function for the LBFGS optimization routine.
// Calculates the gradient over the parameters of the circuit w.r.t. the
// expectation value of the hamiltonian over the in-output of |0> This function
// is not used in the ADADELTA routine (which contains almost a copy of the
// code), because that version performs ADADELTA specific operations while
// calculating the gradient. (This is not the case for LBFGS.)
static Result calculate_gradient(double (*gradient)[],
                                 unsigned int abs_param_count,
                                 OptimizerSettings opt_settings,
                                 double _Complex (*buff_left)[],
                                 double _Complex (*buff_right)[],
                                 unsigned int state_size, Circuit *circuit) {
  unsigned int flat_param_index = 0;
  // For each gate reparameterization...
  Iter param_iter = iter_create_contiguous_memory(opt_settings.reparams,
                                                  sizeof(GateParameterization),
                                                  opt_settings.reparams_count);
  Option next;
  while ((next = iter_next(&param_iter)).some) {
    GateParameterization *param = (GateParameterization *)(next.data);
    Iter subparam_iter = iter_create_contiguous_memory(
        param->params, sizeof(double), param->param_count);
    Iter deltas_iter = iter_create_contiguous_memory(
        param->deltas, sizeof(double), param->param_count);
    // For each parameter in that gate parameterization...
    while ((next = iter_next(&subparam_iter)).some) {
      double *subparam = (double *)(next.data);
      double delta = *(double *)(iter_next(&deltas_iter).data);

      // Left shift
      {
        *subparam -= delta;
        param->gate->reparamFn(&param->gate->matrix, param->params);

        // Initialize left buffer to |0>
        memcpy(buff_left, opt_settings.zero_state,
               sizeof(double _Complex) * state_size);

        // Simulate into left buffer
        {
          Result run_r = circuit_run(circuit, buff_left);
          if (!run_r.valid) {
            iter_free(&subparam_iter);
            iter_free(&deltas_iter);
            iter_free(&param_iter);
            return run_r;
          }
        }

        // Undo shift
        *subparam += delta;
      }

      // Right shift
      {
        *subparam += delta;
        param->gate->reparamFn(&param->gate->matrix, param->params);

        // Initialize right buffer to |0>
        memcpy(buff_right, opt_settings.zero_state,
               sizeof(double _Complex) * state_size);

        // Simulate into right buffer
        {
          Result run_r = circuit_run(circuit, buff_right);
          if (!run_r.valid) {
            iter_free(&subparam_iter);
            iter_free(&deltas_iter);
            iter_free(&param_iter);
            return run_r;
          }
        }

        // Undo right shift; reparameterize
        *subparam -= delta;
        param->gate->reparamFn(&param->gate->matrix, param->params);
      }

      // Calculate this parameter's gradient component via the method
      //  presented above.
      double grad_component = 0;
      {
        for (unsigned int u = 0; u < state_size; ++u) {
          Vector *row = ((Vector *)opt_settings.hamiltonian.data) + u;
          double _Complex left_u = (*buff_left)[u];
          double _Complex right_u = (*buff_right)[u];

          for (unsigned int i = 0; i < row->size; ++i) {
            OptimizerDCPackedRowElem elem =
                *((OptimizerDCPackedRowElem *)(row->data) + i);
            unsigned int k = elem.j;
            double _Complex left_k = (*buff_left)[k];
            double _Complex right_k = (*buff_right)[k];
            double _Complex ham_uk = elem.value;

            if (u == k) {
              grad_component +=
                  creal(ham_uk * (right_u + left_u) * conj(right_u - left_u)) /
                  2.;
            } else {
              grad_component += creal(
                  ham_uk * (conj(right_u) * right_k - conj(left_u) * left_k));
            }
          }
        }

        grad_component /= delta;
        (*gradient)[flat_param_index] = grad_component;

        // Keep track of the "flattened" index of the subparameters.
        flat_param_index++;
      }
    }
    iter_free(&subparam_iter);
    iter_free(&deltas_iter);
  }
  iter_free(&param_iter);

  return result_get_empty_valid();
}

// Optimize a circuit's parameters w.r.t. the expectation value of the
// Hamiltonian using a Low-Memory Broyden-Fletcher-Goldfarb-Shanno algorithm
// TODO
OptimizationResult optimizer_optimize_lbfgs(
    Optimizer *optimizer, OptimizerParamCallback param_callback,
    OptimizerEnergyCallback energy_callback, void *callback_context) {
  // Aliases for commonly used values
  const LbfgsSettings lbfgs_settings =
      optimizer->algo_settings.specific.lbfgs_settings;
  OptimizerSettings opt_settings = optimizer->opt_settings;
  Circuit *circuit = opt_settings.circuit;
  const unsigned int state_size = 1U << circuit->depth[0];

  // Check validity of optimizer
  if (optimizer == NULL) {
    Result invalid_result = result_get_invalid_reason("optimizer is null");
    return result_as_optimization_result(invalid_result);
  }

  // Count the number of parameters
  unsigned int abs_param_count = 0;
  {
    unsigned int reparams_count = opt_settings.reparams_count;
    for (unsigned int reparam_index = 0; reparam_index < reparams_count;
         ++reparam_index) {
      abs_param_count += opt_settings.reparams[reparam_index].param_count;
    }
  }

  // Iteration count; this will only be useful if there is a max
  // iteration count allowed
  unsigned long int iter_count = 0;

  // 2-norm of the descent step
  // Used to break out of the optimization based on `stop_at`
  double step_two_norm = opt_settings.stop_at + 1;

  // Buffers for state at left and right of the parameters
  // This is needed to approximate the gradient
  double _Complex buff_left[abs_param_count];
  double _Complex buff_right[abs_param_count];

  // Gradient over the parameters, w.r.t. the expectation value of the
  //  Hamiltonian over the output
  double gradient[abs_param_count];

  // LBFGS "memory" variables
  Dequeue step_memory;
  Dequeue grad_delta_memory;
  {
    Result init_r;
    init_r = dequeue_init(&step_memory, sizeof(double) * abs_param_count);
    if (!init_r.valid) return result_as_optimization_result(init_r);
    init_r = dequeue_init(&grad_delta_memory, sizeof(double) * abs_param_count);
    if (!init_r.valid) {
      dequeue_free(&step_memory);
      return result_as_optimization_result(init_r);
    }
  }

  // An array to store the gradient of the last step; this is necessary to
  // compute the gradient delta
  double last_gradient[abs_param_count];
  memset(last_gradient, 0, sizeof(double) * abs_param_count);

  // Optimization cycle:
  while (step_two_norm > opt_settings.stop_at) {
    // If there's a maximum number of iterations, keep track of the
    // current iteration (and break as appropriate)
    if (opt_settings.max_iterations.some) {
      iter_count++;
      if (iter_count >= opt_settings.max_iterations.data) break;
    }

    // Choose H_k^0
    {  // Calculate gradient into `gradient`
      Result grad_r =
          calculate_gradient(&gradient, abs_param_count, opt_settings,
                             &buff_left, &buff_right, state_size, circuit);
      if (!grad_r.valid) {
        dequeue_free(&step_memory);
        dequeue_free(&grad_delta_memory);
        return result_as_optimization_result(grad_r);
      }
    }

    // Calculate the gradient delta
    // We need to do this now because `gradient` will be overwritten as a
    // consequence of the "recursion" over q
    double grad_delta[abs_param_count];
    for (unsigned int i = 0; i < abs_param_count; ++i) {
      grad_delta[i] = gradient[i] - last_gradient[i];
    }
    memcpy(last_gradient, gradient, sizeof(double) * abs_param_count);

    // Compute H⁰k
    // We use the approximation suggested by Nocedal (2006), p. 178
    // that H = \gamma I (eq. 7.20)
    // and so store only gamma
    double gamma;
    {
      // Calculate gamma
      // (if there is history to calculate it from, otherwise assume = 1)
      if (grad_delta_memory.size > 0) {
        // Safely peek the most recent gradient delta in memory
        double *last_grad_delta;
        {
          Result peek_r = dequeue_peek_from_front(&grad_delta_memory, 0);
          if (!peek_r.valid) {
            dequeue_free(&step_memory);
            dequeue_free(&grad_delta_memory);
            return result_as_optimization_result(peek_r);
          }
          last_grad_delta = (double *)peek_r.content.data;
        }

        if (step_memory.size == 0) {
          dequeue_free(&step_memory);
          dequeue_free(&grad_delta_memory);
          return result_as_optimization_result(result_get_invalid_reason(
              "size mismatch between memory dequeues"));
        }

        // Safely peek the most recent update delta in memory
        double *last_step;
        {
          Result peek_r = dequeue_peek_from_front(&step_memory, 0);
          if (!peek_r.valid) {
            dequeue_free(&step_memory);
            dequeue_free(&grad_delta_memory);
            return result_as_optimization_result(peek_r);
          }
          last_step = (double *)peek_r.content.data;
        }

        // Calculate the numerator and denominator sum in gamma
        double numerator = 0, denominator = 0;
        for (unsigned int i = 0; i < abs_param_count; ++i) {
          numerator += last_grad_delta[i] * last_step[i];
          denominator += last_grad_delta[i] * last_grad_delta[i];
        }

        // Finally, calculate gamma
        gamma = numerator / denominator;
      } else {
        gamma = 1;
      }
    }

    // Compute H_k Grad(f)_K
    // via two loop recursion
    // (i.e. calculate direction of descent)
    // ***This will overwrite the gradient!***
    // However, the gradient is still stored in `last_gradient`
    // for the rest of the iteration.

    double rho_cache[lbfgs_settings.m];    // Cache `rho`s on the first loop
    double alpha_cache[lbfgs_settings.m];  // Cache `alpha`s on the first loop
    // Note that the caches may have uninitialized elements,
    // if the dequeues have less than `m` elements still

    // First loop
    {
      Iter grad_delta_iter = dequeue_into_iterator_ftb(&grad_delta_memory);
      Iter step_delta_iter = dequeue_into_iterator_ftb(&step_memory);

      // Not the same as the index of the iteration! (i.e. "position" of the
      // iterators) Note that Because one loop is b.t.f. and the other is
      // f.t.b., the caching/access must account for this.
      unsigned int iter_index = grad_delta_memory.size - 1;

      Option next;
      while ((next = iter_next(&grad_delta_iter)).some) {
        double *grad_delta_i = (double *)next.data;

        next = iter_next(&step_delta_iter);
        if (!next.some) {
          dequeue_free(&step_memory);
          dequeue_free(&grad_delta_memory);
          return result_as_optimization_result(result_get_invalid_reason(
              "size mismatch between memory dequeues"));
        }
        double *step_i = (double *)next.data;

        // We now have s_i, y_i

        // Calculate rho
        double rho;
        {
          rho = 0;
          for (unsigned int i = 0; i < abs_param_count; ++i) {
            rho += grad_delta_i[i] * step_i[i];
          }
          rho = 1. / rho;
        }

        rho_cache[iter_index] = rho;

        // Calculate alpha
        double alpha;
        {
          alpha = 0;
          for (unsigned int i = 0; i < abs_param_count; ++i) {
            alpha += step_i[i] * gradient[i];
          }
          alpha = rho * alpha;
        }

        alpha_cache[iter_index] = alpha;

        // Update `q` (which is the gradient array, here)
        for (unsigned int i = 0; i < abs_param_count; ++i) {
          gradient[i] -= alpha * grad_delta_i[i];
        }

        iter_index--;
      }

      iter_free(&grad_delta_iter);
      iter_free(&step_delta_iter);
    }  // Finished first loop of two loop recursion

    // r <- H^0_k q
    // Take into account that we are taking H_0^k ~ gamma I
    {
      for (unsigned int i = 0; i < abs_param_count; ++i) {
        gradient[i] *= gamma;
      }
    }

    // Second loop
    {
      // Note that the iteration is now back to front
      Iter grad_delta_iter = dequeue_into_iterator_btf(&grad_delta_memory);
      Iter step_delta_iter = dequeue_into_iterator_btf(&step_memory);
      unsigned int iter_index = 0;

      Option next;
      while ((next = iter_next(&grad_delta_iter)).some) {
        double *grad_delta_i = (double *)next.data;

        next = iter_next(&step_delta_iter);
        if (!next.some) {
          dequeue_free(&step_memory);
          dequeue_free(&grad_delta_memory);
          return result_as_optimization_result(
              result_get_invalid_reason("size mismatch between dequeues"));
        }
        double *step_i = (double *)next.data;

        // We now have s_i, y_i

        // Calculate beta
        double beta;
        {
          beta = 0;
          for (unsigned int i = 0; i < abs_param_count; ++i) {
            beta += grad_delta_i[i] * gradient[i];
          }
          beta *= rho_cache[iter_index];
        }

        // Update r
        for (unsigned int i = 0; i < abs_param_count; ++i) {
          gradient[i] += step_i[i] * (alpha_cache[iter_index] - beta);
        }

        iter_index++;
      }

      iter_free(&grad_delta_iter);
      iter_free(&step_delta_iter);
    }  // Finished second loop of two loop recursion

    // `gradient` now holds H_k grad(f)_k
    // i.e. `-alpha*gradient` is the direction of descent
    // We will transform the components of gradient as we
    // reparameterize

    // Step all variables, perform the reparameterizations
    // Also calculate `step_two_norm`
    {
      step_two_norm = 0;
      unsigned int flat_param_index = 0;
      // For each gate reparameterization...
      Iter param_iter = iter_create_contiguous_memory(
          opt_settings.reparams, sizeof(GateParameterization),
          opt_settings.reparams_count);
      Option next;
      while ((next = iter_next(&param_iter)).some) {
        GateParameterization *param = (GateParameterization *)(next.data);
        Iter subparam_iter = iter_create_contiguous_memory(
            param->params, sizeof(double), param->param_count);
        Iter deltas_iter = iter_create_contiguous_memory(
            param->deltas, sizeof(double), param->param_count);
        // For each parameter in that gate parameterization...
        while ((next = iter_next(&subparam_iter)).some) {
          double *subparam = (double *)(next.data);
          *subparam -= lbfgs_settings.alpha * gradient[flat_param_index];

          // Also reverse gradient (+ multiply by alpha)
          // Careful that we're reversing *after* subtracting to subparam
          gradient[flat_param_index] =
              -lbfgs_settings.alpha * gradient[flat_param_index];

          // Also calculate two norm of descent
          step_two_norm += pow(gradient[flat_param_index], 2);

          flat_param_index++;
        }
      }
    }

    // If the dequeues are full, pop an item
    if (step_memory.size == lbfgs_settings.m) {
      Result pop_r = dequeue_pop_back(&step_memory, NULL);
      if (!pop_r.valid) {
        dequeue_free(&step_memory);
        dequeue_free(&grad_delta_memory);
        return result_as_optimization_result(pop_r);
      }
    }

    // Push the step to the memory; it's just `gradient`!
    // (Because we reversed and multiplied by alpha before)
    {
      Result push_r = dequeue_push_front(&step_memory, (void *)gradient);
      if (!push_r.valid) {
        dequeue_free(&step_memory);
        dequeue_free(&grad_delta_memory);
        return result_as_optimization_result(push_r);
      }
    }

    if (grad_delta_memory.size == lbfgs_settings.m) {
      Result pop_r = dequeue_pop_back(&grad_delta_memory, NULL);
      if (!pop_r.valid) {
        dequeue_free(&step_memory);
        dequeue_free(&grad_delta_memory);
        return result_as_optimization_result(pop_r);
      }
    }

    // Push the gradient delta to the memory
    {
      Result push_r =
          dequeue_push_front(&grad_delta_memory, (void *)grad_delta);
      if (!push_r.valid) {
        dequeue_free(&step_memory);
        dequeue_free(&grad_delta_memory);
        return result_as_optimization_result(push_r);
      }
    }  // done pushing
  }    // Done with optimize

  // Reparameterize all gates so that they match
  // last calculated parameters
  for (unsigned int reparam_index = 0;
       reparam_index < opt_settings.reparams_count; ++reparam_index) {
    GateParameterization *param = opt_settings.reparams + reparam_index;
    param->gate->reparamFn(&param->gate->matrix, param->params);
  }

  // Free local variables
  dequeue_free(&step_memory);
  dequeue_free(&grad_delta_memory);

  // Done!
  OptimizationResult result;
  if (opt_settings.max_iterations.some &&
      iter_count >= opt_settings.max_iterations.data)
    result.quit_on_max_iter = true;
  else
    result.quit_on_max_iter = false;
  result.valid = true;
  result.content.data = circuit;
  return result;
}