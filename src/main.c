#include <assert.h>
#include <stdio.h>
#include <time.h>
#include "include/circuit.h"
#include "include/gate.h"
#include "include/optimizer.h"
#include "include/option.h"
#include "include/vector.h"

#define CPARTS(C) creal(C), cimag(C)

int main(void) {
  Circuit circuit;
  result_unwrap(circuit_init(&circuit, 2));
  Gate x;
  result_unwrap(gate_init_from_identifier(&x, GateX, NULL));
  Gate x2;
  result_unwrap(gate_init_from_identifier(&x2, GateX, NULL));
  Gate h;
  result_unwrap(gate_init_from_identifier(&h, GateH, NULL));
  Gate rx;
  {
    double params[] = {0};
    result_unwrap(gate_init_from_identifier(&rx, GateRx, params));
  }

  // result_unwrap(circuit_add_gate(&circuit, &x, 1, option_none_uint()));
  // result_unwrap(circuit_add_gate(&circuit, &x2, 0, option_none_uint()));
  // result_unwrap(circuit_add_gate(&circuit, &h, 0, option_from_uint(1)));
  result_unwrap(circuit_add_gate(&circuit, &rx, 0, option_none_uint()));
  result_unwrap(circuit_compact(&circuit));
  result_unwrap(circuit_harden(&circuit));

  double _Complex inout[1 << 2];
  memset(inout, 0, sizeof(inout));
  inout[0] = (double _Complex)1;

  result_unwrap(circuit_run(&circuit, &inout));

  Optimizer optimizer;
  AdadeltaSettings ada_default = optimizer_adadelta_get_default();

  GateParameterization rx_param;
  {
    double param[1] = {0.1};
    double delta[1] = {0.1};
    optimizer_gate_param_init(&rx_param, &rx, 1, param, delta);
  }

  OptimizerSettings opt_settings;
  {
    GateParameterization all_params[1] = {rx_param};
    double _Complex hamiltonian[4][4] = {{1., 0., 0., 0.},
                                         {0., -1., 0., 0.},
                                         {0., 0., 1., 0.},
                                         {0., 0., 0., -1.}};
    optimizer_settings_init(&opt_settings, &circuit,
                            (double _Complex *)hamiltonian, 1e-8, all_params,
                            1);
  }
  result_unwrap(optimizer_init(&optimizer, opt_settings, ada_default));

  {
    time_t start = time(NULL);
    result_unwrap(optimizer_optimize(&optimizer));
    time_t end = time(NULL);

    printf("OPTIMIZATION: \n");
    printf("%e\n", rx_param.params[0]);
    printf("in %d\n", end - start);
  }

  optimizer_settings_free(&opt_settings);
  optimizer_gate_param_free(&rx_param);
  circuit_free(&circuit);
  gate_free(&x);
  gate_free(&x2);
  gate_free(&h);
  gate_free(&rx);
}