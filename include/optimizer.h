/**
 * Defines all the structures needed to optimize a circuit.
 * The optimization is over the expectation value of a given Hamiltonian
 * over the output state of the circuit, given the |000...> state as
 * input.
 * The circuit and its optimization are two separate concepts, and
 * whereas gates are intrinsically parameterizable or not, not all
 * parameterizable gates need to be considered in an optimizer.
 * All methods declared in this file are prefixed with `optimizer_`.
 * The optimizer utilizes the ADADELTA method;
 *    see https://arxiv.org/abs/1212.5701
 **/

#ifndef QOP_OPTIMIZER_H_
#define QOP_OPTIMIZER_H_

#include <complex.h>
#include "include/circuit.h"
#include "include/gate.h"
#include "include/vector.h"

// Hyperparameters pertaining to the ADADELTA method.
typedef struct AdadeltaSettings {
  double rho;      // Decay rate
  double epsilon;  // Divide-by-zero prevention constant
} AdadeltaSettings;

// Constructs an `AdadeltaSettings` struct with some default parameters.
// This makes sense, because one of the key points of the ADADELTA
// method is that it is "robust to [...] selection of hyperparameters".
// The default hyperparameters given are those presented on section 4.1.
// of the ADADELTA paper (arXiv:1212.5701), namely
// `epsilon`   1.0E-6
// `rho`       0.95
AdadeltaSettings optimizer_adadelta_get_default();

// Description of a parameterized gate to be optmimized with an
// `Optimizer`.
// The preffered way to initialize a new `GateParameterization` object
// is with `optimizer_gate_param_init()`.
typedef struct GateParameterization {
  Gate *gate;  // the gate to be optimized (by reference)
  unsigned int param_count;
  double *params;
  double *deltas;  // deltas to consider in approximating a gradient in
                   // first order centered finite difference
} GateParameterization;

// Initializes a new `GateParameterization` object.
// This is the preffered way to create a `GateParameterization`, and
// will `malloc` memory` for and `memcpy` the `params` and `deltas`
// arrays, allowing the originals to be safely freed.
// This requires a call to `optimizer_gate_param_free()` when appropriate,
// otherwise resulting in a memory leak.
Result optimizer_gate_param_init(GateParameterization *gate_param, Gate *gate,
                                 unsigned int param_count, double *params,
                                 double *deltas);

// Frees the memory `malloc`ed for the `params` and `deltas` arrays,
// allocated at `optimizer_gate_param_init()`.
// If the `GateParameterization` was manually constructed (which is not
// recommended), it is best to manually free the appropriate memory.
void optimizer_gate_param_free(GateParameterization *gate_param);

// Object specifying the properties and scope of an `Optimizer` object.
// The preferred way of creating a new `OptimizerSettings` object is with
// the initialization function `optimizer_settings_init`.
// Note that `OptimizerSettings` objects are modified in place during
// optimization!
typedef struct OptimizerSettings {
  Circuit *circuit;
  double _Complex *hamiltonian;
  double stop_at;  // Tolerance to maximum component of gradient
  unsigned int reparams_count;
  GateParameterization *reparams;
  double _Complex *zero_state;  // The |000...> state of the circuit
} OptimizerSettings;

// Initializes a new `OptimizerSettings` object.
// This will `malloc` space for internal use, and a call to
// `optimizer_settings_free()` should be made when appropriate to free
// this memory.
// IMPORTANT: Note that the `hamiltonian` parameter is expected to be
// a `double _Complex *`, and not a 2D array passed by value (which
// would be a `double _Complex (*)[]`). This is no problem, since a
// 2D (row-major) Hamiltonian matrix can be cast directly into a single
// depth pointer, i.e.
//
// ```C
// double _Complex ham[2][2] = {{ 1., 0. },
//                              { 0., 1. }};
// optimizer_settings_init(..., (double _Complex *)ham, ...);
// ```
//
// Also note that the `reparams` array should be kept in scope/
// in memory, since `GateParameterizations` are modified in place during
// optimization.
Result optimizer_settings_init(OptimizerSettings *opt_settings,
                               Circuit *circuit, double _Complex *hamiltonian,
                               double stop_at,
                               GateParameterization *parameterizations,
                               unsigned int parameterizations_count);

// Frees the internal memory allocated at the initialization of the
// given `OptimizerSettings` object.
// This function should called for every `OptimizerSettings` object
// initialized with `optimizer_settings_init()`.
void optimizer_settings_free(OptimizerSettings *opt_settings);

// The optimizer object responsible to perform the optimization itself.
// It is just a glorified pair of `OptimizerSettings` and `AdadeltaSettings`.
typedef struct Optimizer {
  OptimizerSettings opt_settings;
  AdadeltaSettings ada_settings;
} Optimizer;

// Initializes an `Optimizer` object.
// There are no internal memory allocation concerns with the `Optimizer`
// object.
Result optimizer_init(Optimizer *optimizer, OptimizerSettings opt_settings,
                      AdadeltaSettings ada_settings);

// Performs the optimization itself.
// The final (optimized) parameters are stored in the
// `GateParameterization` objects that were passed to the
// `OptimizerSettings` at initialization.
Result optimizer_optimize(Optimizer *optimizer);

#endif  // QOP_OPTIMIZER_H_