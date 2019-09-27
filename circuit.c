#include "circuit.h"

Result circuit_create(unsigned int qubit_count)
{
    Vector soft_gates_vector;
    Result vector_init_r = vector_init(&soft_gates_vector, sizeof(SoftGate), qubit_count);
    if (!vector_init_r.valid)
    {
        return vector_init_r;
    }

    Circuit *new_circuit = malloc(sizeof(Circuit));
    if (new_circuit == NULL)
    {
        return result_get_invalid_reason("could not malloc circuit");
    }

    new_circuit->depth[0] = qubit_count;
    new_circuit->depth[1] = 0;
    new_circuit->soft_gates = soft_gates_vector;

    return result_get_valid_with_data(new_circuit);
}

Result circuit_add_gate(Circuit *circuit, Gate gate, unsigned int qubit, Option_Uint control)
{
    if (circuit == NULL)
    {
        return result_get_invalid_reason("circuit pointer is null");
    }

    if (qubit > circuit->depth[0] - 1)
        return result_get_invalid_reason("qubit is out of bounds");
    if (control.some)
    {
        unsigned int control_value = *(unsigned int *)control.data;
        if (control_value > circuit->depth[0] - 1)
            return result_get_invalid_reason("control is out of bounds");
        if (qubit == control_value)
            return result_get_invalid_reason("qubit cannot be the same as control");
    }

    GatePosition gate_position;
    gate_position.qubit = qubit;
    gate_position.slice = circuit->depth[1];
    circuit->depth[1] += 1;

    SoftGate new_soft_gate;
    new_soft_gate.gate = gate;
    new_soft_gate.control = control;
    new_soft_gate.position = gate_position;

    Result push_r = vector_push(&circuit->soft_gates, &new_soft_gate);
    if (!push_r.valid)
    {
        return push_r;
    }

    return result_get_valid_with_data(circuit);
}

Result circuit_compact(Circuit *circuit)
{
    unsigned int leftmost[circuit->depth[0]];
    Vector_Iter soft_gates_iter = vector_iter_create(&circuit->soft_gates);

    {
        Option next;
        while ((next = vector_iter_next(&soft_gates_iter)).some)
        {
            SoftGate gate = *(SoftGate *)next.data;
            unsigned int new_position;

            if (gate.control.some)
            {
                unsigned int control = *(unsigned int *)(gate.control.data);
                if (leftmost[gate.position.qubit] > leftmost[control])
                {
                    new_position = leftmost[gate.position.qubit];
                }
                else
                {
                    new_position = leftmost[control];
                }

                leftmost[gate.position.qubit] = new_position + 1;
                leftmost[control] = new_position + 1;
            }
            else
            {
                new_position = leftmost[gate.position.qubit];
                leftmost[gate.position.qubit] += 1;
            }

            gate.position.slice = new_position;
        }
    }

    {
        unsigned int max = 0;
        for (unsigned int i = 0; i < circuit->depth[0]; i++)
        {
            if (leftmost[i] > max)
                max = leftmost[i];
        }
        circuit->depth[1] = max;
    }

    return result_get_valid_with_data(circuit);
}

Result circuit_free(Circuit *circuit)
{
    if (circuit == NULL)
    {
        return result_get_invalid_reason("circuit pointer is null");
    }

    Result vector_free_r = vector_free(&circuit->soft_gates);
    if (!vector_free_r.valid)
    {
        return vector_free_r;
    }
    free(circuit);
    circuit = NULL;

    return result_get_empty_valid();
}