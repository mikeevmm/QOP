#include <assert.h>
#include <stdio.h>
#include "include/circuit.h"
#include "include/gate.h"
#include "include/option.h"
#include "include/vector.h"

#define CPARTS(C) creal(C), cimag(C)

int main() {
  Circuit circuit;
  result_unwrap(circuit_init(&circuit, 2));
  Gate x ;
  result_unwrap(gate_init_from_identifier(&x, GateX, NULL));
  Gate x2;
  result_unwrap(gate_init_from_identifier(&x2, GateX, NULL));
  Gate h ;
  result_unwrap(gate_init_from_identifier(&h, GateH, NULL));
  
  result_unwrap(circuit_add_gate(&circuit, &x, 1, option_none_uint()));
  result_unwrap(circuit_add_gate(&circuit, &x2, 0, option_none_uint()));
  result_unwrap(circuit_add_gate(&circuit, &h, 0, option_from_uint(1)));
  result_unwrap(circuit_compact(&circuit));

  double _Complex inout[1 << 2];
  memset(inout, 0, sizeof(inout));
  inout[0] = (double _Complex)1;

  result_unwrap(circuit_run(&circuit, &inout));

  for (int i = 0; i < (1 << 2); i++) {
    printf("%e+i*%e, ", CPARTS(inout[i]));
  }
  printf("\n");

  circuit_free(&circuit);
  gate_free(&x);
  gate_free(&x2);
  gate_free(&h);
}