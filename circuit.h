#pragma once
#include "vector.h"
#include "gate.h"
#include "option.h"

typedef struct GatePosition
{
    unsigned int slice;
    unsigned int qubit;
} GatePosition;

typedef struct SoftGate
{
    GatePosition position;
    Option_Uint control;
    Gate gate;
} SoftGate;

typedef struct Circuit
{
    Vector soft_gates;
    unsigned int depth[2];
    Gate *hardened_gates;
} Circuit;

Result circuit_create(unsigned int qubit_count);
Result circuit_add_gate(Circuit *circuit, Gate gate, unsigned int qubit, Option_Uint control);
Result circuit_compact(Circuit *circuit);
Result circuit_harden(Circuit *circuit);
Result circuit_run(Circuit *circuit, double _Complex (*inout)[]);
Result circuit_free(Circuit *circuit);