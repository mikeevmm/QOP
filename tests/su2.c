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

static const unsigned int qubit_count = 4;
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

  // Random seed 
  srand(time(NULL));

  // First push all gates, to make sure their memory position doesn't
  // change
  for (int l = 0; l < layer_count; ++l) {
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < qubit_count; ++i) {
		double param[1];
		Option_Uint control;

		Gate rx;
		param[0] = ((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926;
		result_unwrap(gate_init_from_identifier(&rx, GateRx, param));
		result_unwrap(vector_push(&gates, &rx));
		result_unwrap(vector_push(&qubits, &i));
		control = option_none_uint();
		result_unwrap(vector_push(&controls, &control));

        Gate ry;
		param[0] = ((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926;
        result_unwrap(gate_init_from_identifier(&ry, GateRy, param));
        result_unwrap(vector_push(&gates, &ry));
        result_unwrap(vector_push(&qubits, &i));
        control = option_none_uint();
        result_unwrap(vector_push(&controls, &control));

        Gate rz;
		param[0] = ((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926;
        result_unwrap(gate_init_from_identifier(&rz, GateRz, param));
        result_unwrap(vector_push(&gates, &rz));
        result_unwrap(vector_push(&qubits, &i));
        control = option_none_uint();
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
	double param[1];
	Option_Uint control;

    Gate rx;
	param[0] = ((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926;
    result_unwrap(gate_init_from_identifier(&rx, GateRx, param));
    result_unwrap(vector_push(&gates, &rx));
    result_unwrap(vector_push(&qubits, &i));
    control = option_none_uint();
    result_unwrap(vector_push(&controls, &control));

    Gate ry;
	param[0] = ((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926;
    result_unwrap(gate_init_from_identifier(&ry, GateRy, param));
    result_unwrap(vector_push(&gates, &ry));
    result_unwrap(vector_push(&qubits, &i));
    control = option_none_uint();
    result_unwrap(vector_push(&controls, &control));

    Gate rz;
	param[0] = ((double)rand() / RAND_MAX * 2. - 1.) * 3.1415926;
    result_unwrap(gate_init_from_identifier(&rz, GateRz, param));
    result_unwrap(vector_push(&gates, &rz));
    result_unwrap(vector_push(&qubits, &i));
    control = option_none_uint();
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
      if (gate->id == GateRx || gate->id == GateRy || gate->id == GateRz) {
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

  // Perform a test run and check that it does not result in 0
  
  double _Complex state[1U << (unsigned int)qubit_count];
  state[0] = 1.;
  circuit_run(&circuit, &state);

  printf("[");
  for (unsigned int i = 0; i < 1U << (unsigned int)qubit_count; ++i)
  {
    printf("%f+i%f", CPARTS(state[i]));	
	if (i != (1U << (unsigned int)qubit_count)-1) printf(",");
  }
  printf("]\n");

  printf("Checking norm\n");
  double norm = 0.;
  for (unsigned int i = 0; i < 1U << (unsigned int)qubit_count; ++i)
  {
    norm += state[i] * conj(state[i]);
  }
  printf("%f\n", norm);

  // Free memory
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
