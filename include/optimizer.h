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
#include "include/option.h"
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

typedef struct LbfgsSettings {
  // TODO:
} LbfgsSettings;

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
  Vector hamiltonian;
  double stop_at;  // Tolerance to maximum component of gradient
  unsigned int reparams_count;
  GateParameterization *reparams;
  double _Complex *zero_state;  // The |000...> state of the circuit
  Option_Uint max_iterations;
} OptimizerSettings;

// This is a purely internal use struct; it's used to identify the actual
// position of a collapsed element in a compacted sparse matrix.
typedef struct OptimizerDCPackedRowElem {
  double _Complex value;
  unsigned int j;
} OptimizerDCPackedRowElem;

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
// If no maximum number of iterations in the optimization cycle is to be
// considered, pass `-1` as the `max_iterations` argument.
Result optimizer_settings_init(OptimizerSettings *opt_settings,
                               Circuit *circuit, double _Complex *hamiltonian,
                               double stop_at,
                               GateParameterization *parameterizations,
                               unsigned int parameterizations_count,
                               Option_Uint max_iterations);

// Frees the internal memory allocated at the initialization of the
// given `OptimizerSettings` object.
// This function should called for every `OptimizerSettings` object
// initialized with `optimizer_settings_init()`.
void optimizer_settings_free(OptimizerSettings *opt_settings);

// The optimizer object responsible to perform the optimization itself.
// It is just a glorified pair of `OptimizerSettings` and some algorithm
// settings struct.
typedef struct Optimizer {
  OptimizerSettings opt_settings;
  union OptimizerAlgoSettings {
    AdadeltaSettings ada_settings;
    LbfgsSettings lbfgs_settings;
  } algo_settings;
} Optimizer;

// Initializes an `Optimizer` object.
// There are no internal memory allocation concerns with the `Optimizer`
// object.
Result optimizer_init(Optimizer *optimizer, OptimizerSettings opt_settings,
                      AdadeltaSettings ada_settings);

// Information about the results of an optimization.
// The structure is that of a simple `Result`, with extra fields specific
// to optimization at the end. This allows for casting to/from the `Result`
// struct, which may be useful for, e.g., returning an invalid result.
// The preferred way to perform this cast is via the
// `optimization_result_as_result()` function.
typedef struct OptimizationResult {
  bool valid;
  union ResultContent content;
  bool quit_on_max_iter;
} OptimizationResult;

Result optimization_result_as_result(OptimizationResult opt_result);
OptimizationResult result_as_optimization_result(Result result);

// Callback types for the optimization function.
// These callbacks are useful for inspecting the evolution of the
// optimization process, but are not a property of the optimization
// itself, and so are passed to the `optimizer_optimize` function,
// rather than included in the `Optimizer` object.
// The signatures for these callbacks are, respectively,
//  (flat parameter index, new value, context variables) -> void
//  (expectation value of hamiltonian, context variables) -> void
typedef void (*OptimizerParamCallback)(unsigned int, double, void *);
typedef void (*OptimizerEnergyCallback)(_Complex double, void *);

// Performs the optimization itself.
// The final (optimized) parameters are stored in the
// `GateParameterization` objects that were passed to the
// `OptimizerSettings` at initialization, and information
// about the optimization process is returned in the form
// of an `OpimizationResult` object.
// The callbacks are optional, and can be NULL.
// Passing a non-NULL callback may result in extra calculations.
OptimizationResult optimizer_optimize(Optimizer *optimizer,
                                      OptimizerParamCallback param_callback,
                                      OptimizerEnergyCallback energy_callback,
                                      void *callback_context);

// Perform an optimization using ADADELTA
// TODO: Comment
OptimizationResult optimizer_optimize_adadelta(
    Optimizer *optimizer, OptimizerParamCallback param_callback,
    OptimizerEnergyCallback energy_callback, void *callback_context);

// TODO: Comment
OptimizationResult optimizer_optimize_lbfgs(
    Optimizer *optimizer, OptimizerParamCallback param_callback,
    OptimizerEnergyCallback energy_callback, void *callback_context);

#endif  // QOP_OPTIMIZER_H_