#include "main.h"

int main()
{
    Circuit circuit = *(Circuit *)result_unwrap(circuit_create(2));
    Gate x = *(Gate *)result_unwrap(gate_new_from_identifier(GateX, NULL));
    Gate h = *(Gate *)result_unwrap(gate_new_from_identifier(GateH, NULL));
    result_unwrap(circuit_add_gate(&circuit, x, 1, option_none_uint()));
    result_unwrap(circuit_add_gate(&circuit, h, 0, option_from_uint(1)));
    result_unwrap(circuit_compact(&circuit));
    result_unwrap(circuit_harden(&circuit));

    double _Complex inout[1 << 2];
    memset(inout, 0, sizeof(inout));
    inout[0] = 1;

    circuit_run(&circuit, &inout);
}