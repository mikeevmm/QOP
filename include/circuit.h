#pragma once
#include "include/gate.h"
#include "include/iter.h"
#include "include/option.h"
#include "include/vector.h"

typedef struct GatePosition {
  unsigned int slice;
  unsigned int qubit;
} GatePosition;

typedef struct SoftGate {
  GatePosition position;
  Option_Uint control;
  Gate gate;
} SoftGate;

typedef struct Circuit {
  Vector soft_gates;
  Vector slice_gate_count;
  unsigned int depth[2];
  SoftGate** hardened_gates;
} Circuit;

bool _circuit_filter_is_soft_gate(void* ptr);
Filter circuit_filter_slice_soft_gates(Circuit* circuit, unsigned int slice);

Result circuit_create(unsigned int qubit_count);
Result circuit_add_gate(Circuit* circuit, Gate* gate, unsigned int qubit,
                        Option_Uint control);
Result circuit_compact(Circuit* circuit);
Result circuit_run(Circuit* circuit, double _Complex (*inout)[]);
Result circuit_free(Circuit* circuit);