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
    new_circuit->hardened_gates = NULL;

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
        unsigned int control_value = control.data;
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
    if (circuit == NULL)
    {
        return result_get_invalid_reason("circuit pointer is null");
    }

    unsigned int leftmost[circuit->depth[0]];
    Vector_Iter soft_gates_iter = vector_iter_create(&circuit->soft_gates);

    {
        Option next;
        while (ITER_NEXT(next, soft_gates_iter))
        {
            SoftGate soft_gate = *(SoftGate *)next.data;
            unsigned int new_position;

            if (soft_gate.control.some)
            {
                unsigned int control = *(unsigned int *)(soft_gate.control.data);
                if (leftmost[soft_gate.position.qubit] > leftmost[control])
                {
                    new_position = leftmost[soft_gate.position.qubit];
                }
                else
                {
                    new_position = leftmost[control];
                }

                leftmost[soft_gate.position.qubit] = new_position + 1;
                leftmost[control] = new_position + 1;
            }
            else
            {
                new_position = leftmost[soft_gate.position.qubit];
                leftmost[soft_gate.position.qubit] += 1;
            }

            soft_gate.position.slice = new_position;
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

Result circuit_harden(Circuit *circuit)
{
    if (circuit == NULL)
    {
        return result_get_invalid_reason("circuit pointer is null");
    }

    if (circuit->hardened_gates != NULL)
    {
        free(circuit->hardened_gates);
    }

    Gate **hardened_gates = malloc(circuit->depth[0] * circuit->depth[1] * sizeof(Gate *));
    if (hardened_gates == NULL)
    {
        return result_get_invalid_reason("could not malloc hardened_gates");
    }

    Vector_Iter soft_gates_iter = vector_iter_create(&circuit->soft_gates);
    Option next;
    while (ITER_NEXT(next, soft_gates_iter))
    {
        SoftGate gate = *(SoftGate *)next.data;
        *(hardened_gates + gate.position.qubit + circuit->depth[0] * gate.position.slice) =
            &gate.gate;
    }

    circuit->hardened_gates = hardened_gates;

    return result_get_valid_with_data(circuit);
}

Result circuit_run(Circuit *circuit, double _Complex inout[])
{
    if (circuit == NULL)
    {
        return result_get_invalid_reason("circuit pointer is null");
    }

    if (circuit->hardened_gates == NULL)
    {
        return result_get_invalid_reason("circuit hardened_gates is null");
    }

    if (sizeof(inout) / sizeof(double _Complex) != (1 << (circuit->depth[0])))
    {
        return result_get_invalid_reason("inout array has wrong length");
    }

    double _Complex output[1 << (circuit->depth[0])];

    const size_t gates_count = circuit->soft_gates.size;
    unsigned int leftmost_bit = 1 << (circuit->depth[0] - 1);
    unsigned int mask = 0;
    unsigned int cmask = 0;
    unsigned int ncmask = 0;
    unsigned int relevant_count = 0;
    {
        Vector_Iter gates_iter = vector_iter_create(&circuit->soft_gates);
        Option next;
        while (ITER_NEXT(next, gates_iter))
        {
            SoftGate soft_gate = *(SoftGate *)next.data;
            mask |= 1 << (soft_gate.position.qubit);
            if ((mask ^ cmask) != 0)
                relevant_count += 1;

            cmask |= mask;
            if (soft_gate.control.some)
            {
                if ((cmask & (1 << soft_gate.control.data)) == 0)
                    relevant_count += 1;
                cmask |= (1 << soft_gate.control.data);
            }
        }
    }
    ncmask = (~cmask) & (leftmost_bit | (leftmost_bit - 1));

    unsigned int x = 0;
    for (int u = 0; u < (1 << relevant_count); ++u)
    {
        for (int y = 0; y < (1 << gates_count); ++y)
        {
            double _Complex coef = 0;
            unsigned int proj = 0;

            Vector_Iter gates_iter = vector_iter_create(&circuit->soft_gates);
            Option next;
            while (ITER_NEXT(next, gates_iter))
            {
                SoftGate gate = *(SoftGate *)next.data;
                unsigned int j = gates_iter.position;
                bool cset;

                proj |= (y >> j) << gate.position.qubit;
                cset = (!gate.control.some) ||
                       (((x >> gate.control.data) & 1) == 1);

                if (((x >> gate.position.qubit) & 1) == 1)
                {
                    if (cset)
                    {
                        coef *= gate.gate.matrix[1][(y >> j) & 1];
                    }
                    else
                    {
                        coef *= (y >> j) & 1;
                    }
                }
                else
                {
                    if (cset)
                    {
                        coef *= gate.gate.matrix[0][(y >> j) & 1];
                    }
                    else
                    {
                        coef *= 1 ^ (1 & (y >> j));
                    }
                }
            }

            unsigned int z = 0;
            for (unsigned int j = 0; j < (circuit->depth[0] - relevant_count); ++j)
            {
                unsigned int out_index = (proj & mask) | ((z | x) & (~mask));
                output[out_index] += coef * inout[z | x];

                z += 1;
                while ((z & ncmask) != z)
                {
                    unsigned int p = ~((z & -z) | ((z & -z) - 1));
                    z = (z & p) + ((ncmask & p) & (-(ncmask & p)));
                }
            }
        }

        x += 1;
        while ((x & cmask) != x)
        {
            unsigned int p = ~((x & -x) | ((x & -x) - 1));
            x = (x & p) + ((cmask & p) & (-(cmask & p)));
        }
    }

    void* copy = memcpy(&inout, &output, sizeof(output));
    if (copy == NULL)
        return result_get_invalid_reason("memcpy failed");
    return result_get_valid_with_data(copy);
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

    if (circuit->hardened_gates != NULL)
        free(circuit->hardened_gates);
    free(circuit);
    circuit = NULL;

    return result_get_empty_valid();
}