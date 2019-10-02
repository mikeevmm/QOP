/**
 *  Declaration of all structures necessary to the definition of a gate
 * in an isolated sense. Also declared are macros and constants relevant
 * to operations with the `Gate` struct.
 * 
 * All methods that act upon a `Gate` are prefixed with `gate_`.
 **/

#pragma once
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "include/option.h"

// Macro to print the (complex) matrix of a `Gate` object, NumPy style.
#define PRINT_GATE_MATRIX(GATE)                                        \
  printf("[[%e+i*%e, %e+i*%e],\n [%e+i*%e, %e+i*%e]]\n",               \
         creal(GATE.gate.matrix[0][0]), cimag(GATE.gate.matrix[0][0]), \
         creal(GATE.gate.matrix[0][1]), cimag(GATE.gate.matrix[0][1]), \
         creal(GATE.gate.matrix[1][0]), cimag(GATE.gate.matrix[1][0]), \
         creal(GATE.gate.matrix[0][1]), cimag(GATE.gate.matrix[0][1]));

// Gate identifier for human-readable or filtering purposes.
typedef enum GateId {
  GateNoId,
  GateCustom,
  GateI,
  GateX,
  GateY,
  GateZ,
  GateH,
  GateSqrtX,
  GateT,
  GateS,
  GateRx,
  GateRy,
  GateRz,
} GateId;

// Constant gate size for use with malloc or memcpy
static const size_t GATE_SINGLE_QUBIT_SIZE = 4 * sizeof(double _Complex);

// Explicit static gate declarations for memcpy in `Gate` object creation
static double _Complex gate_static_i[2][2] = {{1., 0.}, {0., 1.}};
static double _Complex gate_static_x[2][2] = {{0., 1.}, {1., 0.}};
static double _Complex gate_static_y[2][2] = {{0., -1. * _Complex_I},
                                              {1. * _Complex_I, 0.}};
static double _Complex gate_static_z[2][2] = {{1., 0.}, {0., -1.}};
static double _Complex gate_static_h[2][2] = {
    {1.414213562373095, 1.414213562373095},
    {1.414213562373095, -1.414213562373095}};
static double _Complex gate_static_sqrt_x[2][2] = {
    {0.5 + 0.5 * _Complex_I, 0.5 - 0.5 * _Complex_I},
    {0.5 - 0.5 * _Complex_I, 0.5 + 0.5 * _Complex_I}};
static double _Complex gate_static_t[2][2] = {
    {1.0, 0.0},
    {0.0, 7.0710678118654752E-1 + 7.0710678118654752E-1 * _Complex_I}};
static double _Complex gate_static_s[2][2] = {{1.0, 0.0}, {0.0, _Complex_I}};

// `Gate` reparameterization function pointer type, of signature
//    <pointer to> `(matrix, parameter array) -> void`
typedef void (*ReparamFnPtr)(double _Complex[2][2], double[]);

// Representation of a `Gate` as an isolated object (i.e. outside of the
// context of a circuit).
// Note that the reparameterization function pointer `reparamFn` may be
// NULL if the gate is not parameterized/reparameterizable.  
typedef struct Gate {
  double _Complex matrix[2][2];
  ReparamFnPtr reparamFn;
  enum GateId id;
} Gate;

void gate_reparameterize_rx(double _Complex matrix[2][2], double params[]);
void gate_reparameterize_ry(double _Complex matrix[2][2], double params[]);
void gate_reparameterize_rz(double _Complex matrix[2][2], double params[]);

// Creates a new `Gate` object from a given matrix and reparameterization
// function pointer, storing it in the heap, and returns a pointer.
// Note that the `matrix` is copied into the `Gate` struct (which is
// possible as it is of known fixed size), and so the original may be
// safely freed or popped from stack.
// The same is not true of the reparameterization function, as it is
// passed by reference.
// The Gate object itself is stored in the heap, and should be freed
// with `gate_free` when appropriate!
Result gate_new_from_matrix(double _Complex matrix[2][2], ReparamFnPtr reparamFn);

// Creates a new gate from a `enum GateId` gate identifier, storing it
// in the heap, and returns a pointer.
// Some gates (such as the rotation gates) expect parameters, to be
// specified via `params`, but otherwise it is safe to pass `NULL` as
// the `params` argument.
Result gate_new_from_identifier(GateId identifier, double params[]);

// Unwraps the result of a `gate_new_...` call; attempts to unwrap
// the given `Result`; if successful, frees the heap memory of the data,
// returning a stack copy of it.
Gate gate_new_unwrap(Result result);

// Frees all of the heap memory internally allocated by the given `Gate`.
// Note that this function does not free `gate` itself, as it cannot
// ensure that the variable does not live in the stack.
void gate_free(Gate *gate);