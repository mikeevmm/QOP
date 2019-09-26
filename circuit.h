#pragma once
#include "gate.h"
#include "option.h"

typedef struct CircuitPosition
{
    unsigned int slice;
    unsigned int qubit;
} CircuitPosition;

typedef struct SoftGate
{
    CircuitPosition position;
    Option control;
    Gate *gate;
} SoftGate;