#include "main.h"

int main()
{
    Circuit circuit = *(Circuit *)result_unwrap(circuit_create(5));
    Gate x = *(Gate *)result_unwrap(gate_new_from_identifier(GateX, NULL));
    circuit_add_gate(&circuit, x, 1, option_none_uint());
}