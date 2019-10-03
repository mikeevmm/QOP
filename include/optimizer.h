#ifndef QOP_OPTIMIZER_H_
#define QOP_OPTIMIZER_H_

#include <complex.h>
#include "include/circuit.h"
#include "include/gate.h"
#include "include/vector.h"

typedef struct AdadeltaSettings {
  double rho;      // Decay rate
  double epsilon;  // Divide-by-zero prevention constant
} AdadeltaSettings;

AdadeltaSettings optimizer_adadelta_get_default();

typedef struct GateParameterization {
  Gate *gate;
  unsigned int param_count;
  double *params;
  double *deltas;
} GateParameterization;

Result optimizer_gate_param_init(GateParameterization *gate_param, Gate *gate,
                                 unsigned int param_count, double *params,
                                 double *deltas);
void optimizer_gate_param_free(GateParameterization *gate_param);

typedef struct OptimizerSettings {
  Circuit *circuit;
  double _Complex *hamiltonian;
  double stop_at;  // Tolerance to maximum component of update to parameters
  unsigned int parameterization_count;
  GateParameterization *parameterizations;
  double _Complex *zero_state;
} OptimizerSettings;

Result optimizer_settings_init(OptimizerSettings *opt_settings,
                               Circuit *circuit, double _Complex *hamiltonian,
                               double stop_at,
                               GateParameterization *parameterizations,
                               unsigned int parameterizations_count);

void optimizer_settings_free(OptimizerSettings *opt_settings);

typedef struct Optimizer {
  OptimizerSettings opt_settings;
  AdadeltaSettings ada_settings;
} Optimizer;

Result optimizer_init(Optimizer *optimizer, OptimizerSettings opt_settings,
                      AdadeltaSettings ada_settings);

Result optimizer_optimize(Optimizer *optimizer);

#endif  // QOP_OPTIMIZER_H_