#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "include/circuit.h"
#include "include/gate.h"
#include "include/optimizer.h"
#include "include/option.h"
#include "include/vector.h"

#define CPARTS(C) creal(C), cimag(C)

static const double _Complex hamiltonian[4][4] = {
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}};
static const unsigned int qubit_count = 2;
static const unsigned int layer_count = 2;

int main(void) {
  Circuit circuit;
  result_unwrap(circuit_init(&circuit, qubit_count));

  Vector gates;
  result_unwrap(vector_init(&gates, sizeof(Gate), 0));
  Vector qubits;
  result_unwrap(vector_init(&qubits, sizeof(unsigned int), 0));
  Vector controls;
  result_unwrap(vector_init(&controls, sizeof(Option_Uint), 0));

  Vector reparams;
  result_unwrap(vector_init(&reparams, sizeof(GateParameterization), 0));

  // Fix the seed so that results are reproducible
  // Paul Dirac was born on 8th Aug 1902 :)
  //srand(8081902);
  srand(time(NULL));

  // First push all gates, to make sure their memory position doesn't
  // change
  for (int l = 0; l < layer_count; ++l) {
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < qubit_count; ++i) {
        Gate ry;
        double param[] = {((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926};
        result_unwrap(gate_init_from_identifier(&ry, GateRy, param));
        result_unwrap(vector_push(&gates, &ry));
        result_unwrap(vector_push(&qubits, &i));
        Option_Uint control = option_none_uint();
        result_unwrap(vector_push(&controls, &control));
      }

      for (unsigned int i = 0; i < qubit_count; i += 2) {
        Gate cz;
        result_unwrap(gate_init_from_identifier(&cz, GateZ, NULL));
        result_unwrap(vector_push(&gates, &cz));
        unsigned int qubit = (i + k) % qubit_count;
        result_unwrap(vector_push(&qubits, &qubit));
        Option_Uint control = option_from_uint((i + k + 1) % qubit_count);
        result_unwrap(vector_push(&controls, &control));
      }
    }
  }

  for (int i = 0; i < qubit_count; ++i) {
    Gate ry;
    double param[] = {0.01};
    result_unwrap(gate_init_from_identifier(&ry, GateRy, param));
    result_unwrap(vector_push(&gates, &ry));
    result_unwrap(vector_push(&qubits, &i));
    Option_Uint control = option_none_uint();
    result_unwrap(vector_push(&controls, &control));
  }

  // Push reparameterizations & gates into circuit
  {
    Iter gates_iter = vector_iter_create(&gates);
    Iter qubits_iter = vector_iter_create(&qubits);
    Iter controls_iter = vector_iter_create(&controls);
    Option next;
    double delta[] = {0.01};
    while ((next = iter_next(&gates_iter)).some) {
      double param[] = {((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926};
      Gate *gate = (Gate *)next.data;

      // Reparam
      if (gate->id == GateRy) {
        GateParameterization gparam;
        result_unwrap(
            optimizer_gate_param_init(&gparam, gate, 1, param, delta));
        result_unwrap(vector_push(&reparams, &gparam));
      }

      // Circuit
      unsigned int qubit = *(unsigned int *)iter_next(&qubits_iter).data;
      Option_Uint control = *(Option_Uint *)iter_next(&controls_iter).data;
      result_unwrap(circuit_add_gate(&circuit, gate, qubit, control));
    }
    iter_free(&gates_iter);
    iter_free(&qubits_iter);
    iter_free(&controls_iter);
  }

  result_unwrap(circuit_compact(&circuit));
  result_unwrap(circuit_harden(&circuit));

  OptimizerSettings opt_settings;
  result_unwrap(optimizer_settings_init(
      &opt_settings, &circuit, (double _Complex *)hamiltonian, 1e-2,
      reparams.data, reparams.size, option_none_uint()));
  OptimizerAlgoSettings algo_settings;
  {
    LbfgsSettings lbfgs_settings = optimizer_lbfgs_get_default();
    lbfgs_settings.alpha = 0.09;
    optimizer_algo_settings_init(&algo_settings, AlgoLbfgs,
                                 (void *)(&lbfgs_settings));
    /*
    AdadeltaSettings ada_settings = optimizer_adadelta_get_default();
    optimizer_algo_settings_init(&algo_settings, AlgoAdadelta,
                                 (void *)(&ada_settings));
    */
  }
  Optimizer opt;
  result_unwrap(optimizer_init(&opt, opt_settings, algo_settings));
  result_unwrap(optimization_result_as_result(
      optimizer_optimize(&opt, NULL, NULL, NULL)));

  double _Complex zero_state[1UL << qubit_count];
  zero_state[0] = 1;
  for (unsigned int i = 1; i < (1U << qubit_count); ++i) zero_state[i] = 0.;
  result_unwrap(circuit_run(&circuit, &zero_state));

  double norm = 0;
  for (int i = 0; i < (1UL << qubit_count); ++i)
    norm += (double)(zero_state[i] * conj(zero_state[i]));
  printf("NORM %e\n", norm);

  printf("PARAMETERS:\n");
  {
    Iter params_iter = vector_iter_create(&reparams);
    Option next;
    while ((next = iter_next(&params_iter)).some) {
      GateParameterization *param = (GateParameterization *)next.data;

      for (int i = 0; i < param->param_count; ++i) {
        printf("%e ", *(param->params + i));
      }
      printf("\n");
    }
    iter_free(&params_iter);
  }

  printf("FINAL EXPECTATION VALUE:\n");
  {
    double _Complex energy = 0;
    double _Complex current_phi_state[1U << qubit_count];
    memcpy(current_phi_state, opt_settings.zero_state,
           sizeof(double _Complex) * (1U << qubit_count));
    circuit_run(&circuit, &current_phi_state);

    for (unsigned int i = 0; i < opt_settings.hamiltonian.size; ++i) {
      Vector *row = (Vector *)opt_settings.hamiltonian.data + i;
      for (unsigned int u = 0; u < row->size; ++u) {
        OptimizerDCPackedRowElem elem =
            *((OptimizerDCPackedRowElem *)row->data + u);
        unsigned int j = elem.j;
        double _Complex value = elem.value;
        if (i == j)
          energy += conj(current_phi_state[i]) * value * current_phi_state[j];
        else
          energy +=
              2 * conj(current_phi_state[i]) * value * current_phi_state[j];
      }
    }
    printf("%e %e\n", creal(energy), cimag(energy));
  }

  optimizer_settings_free(&opt_settings);
  {
    Iter gates_iter = vector_iter_create(&gates);
    Option next;
    while ((next = iter_next(&gates_iter)).some) {
      gate_free((Gate *)next.data);
    }
    iter_free(&gates_iter);
  }
  vector_free(&gates);
  {
    Iter params_iter = vector_iter_create(&reparams);
    Option next;
    while ((next = iter_next(&params_iter)).some) {
      optimizer_gate_param_free((GateParameterization *)next.data);
    }
    iter_free(&params_iter);
  }
  vector_free(&reparams);
  circuit_free(&circuit);
}