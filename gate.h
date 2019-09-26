#pragma once
#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

enum GateId
{
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
};

// Constant sizes for use with malloc
static const int GATE_SINGLE_QUBIT_SIZE = 4 * sizeof(double _Complex);

// Explicit static gate declarations
static double _Complex static_gate_i[2][2] = {{1., 0.},
                                              {0., 1.}};
static double _Complex static_gate_x[2][2] = {{0., 1.},
                                              {1., 0.}};
static double _Complex static_gate_y[2][2] = {{0., -1. * _Complex_I},
                                              {1. * _Complex_I, 0.}};
static double _Complex static_gate_z[2][2] = {{1., 0.},
                                              {0., -1.}};
static double _Complex static_gate_h[2][2] = {{1.414213562373095, 1.414213562373095},
                                              {1.414213562373095, -1.414213562373095}};
static double _Complex static_gate_sqrt_x[2][2] = {{0.5 + 0.5 * _Complex_I, 0.5 - 0.5 * _Complex_I},
                                                   {0.5 - 0.5 * _Complex_I, 0.5 + 0.5 * _Complex_I}};
static double _Complex static_gate_t[2][2] = {{1.0, 0.0},
                                              {0.0, 7.0710678118654752E-1 + 7.0710678118654752E-1 * _Complex_I}};
static double _Complex static_gate_s[2][2] = {{1.0, 0.0},
                                              {0.0, _Complex_I}};

// Reparameterization function signature
//  (matrix head pointer, parameter) -> void
typedef void (*ReparamFnPtr)(double _Complex[2][2], double[]);

struct Gate
{
    // qubit matrix
    double _Complex matrix[2][2];
    // Reparameterization function pointer
    // (can be NULL)
    ReparamFnPtr reparamFn;
    // Gate human identifier
    enum GateId id;
};

// Reparameterization functions
void reparameterize_rx_gate(double _Complex matrix[2][2], double params[]);
void reparameterize_ry_gate(double _Complex matrix[2][2], double params[]);
void reparameterize_rz_gate(double _Complex matrix[2][2], double params[]);

// Create a gate from a matrix and an optional reparameterization function
struct Gate gate_new_from_matrix(double _Complex matrix[2][2], double params[], ReparamFnPtr reparamFn);