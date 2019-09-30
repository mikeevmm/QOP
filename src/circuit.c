#include "circuit.h"

unsigned int soft_gate_flat_position(Circuit *circuit, SoftGate *soft_gate)
{
    return soft_gate->position.qubit + soft_gate->position.slice * circuit->depth[0];
}

Filter circuit_filter_slice_soft_gates(Circuit *circuit, unsigned int slice)
{
    Iter slice_iter = iter_create(slice * circuit->depth[0],
                                  sizeof(SoftGate *),
                                  circuit->depth[1]);
    Filter slice_filter = filter_create(slice_iter, filter_generic_not_null);
    return slice_filter;
}

Result circuit_create(unsigned int qubit_count)
{
    if (qubit_count == 0)
    {
        return result_get_invalid_reason("cannot create a 0-qubit circuit");
    }

    Vector soft_gates_vector;
    {
        Result vector_create_r = vector_create(&soft_gates_vector, sizeof(SoftGate), qubit_count);
        if (!vector_create_r.valid)
        {
            return vector_create_r;
        }
    }

    Vector slice_gate_vector;
    {
        Result vector_create_r = vector_create(&slice_gate_vector, sizeof(unsigned int), 0);
        if (!vector_create_r.valid)
        {
            return vector_create_r;
        }
    }

    Circuit *new_circuit = (Circuit *)malloc(sizeof(Circuit));
    if (new_circuit == NULL)
    {
        return result_get_invalid_reason("could not malloc circuit");
    }

    new_circuit->depth[0] = qubit_count;
    new_circuit->depth[1] = 0;
    new_circuit->soft_gates = soft_gates_vector;
    new_circuit->slice_gate_count = slice_gate_vector;
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
    memset(leftmost, 0, sizeof(leftmost));

    {
        Iter soft_gates_iter = vector_iter_create(&circuit->soft_gates);
        Option next;
        while (ITER_NEXT(next, soft_gates_iter))
        {
            SoftGate *soft_gate = (SoftGate *)next.data;
            unsigned int new_position = 0;

            if (soft_gate->control.some)
            {
                unsigned int control = soft_gate->control.data;
                if (leftmost[soft_gate->position.qubit] > leftmost[control])
                {
                    new_position = leftmost[soft_gate->position.qubit];
                }
                else
                {
                    new_position = leftmost[control];
                }

                leftmost[soft_gate->position.qubit] = new_position + 1;
                leftmost[control] = new_position + 1;
            }
            else
            {
                new_position = leftmost[soft_gate->position.qubit];
                leftmost[soft_gate->position.qubit] += 1;
            }

            soft_gate->position.slice = new_position;
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

    if (circuit->hardened_gates != NULL)
    {
        free(circuit->hardened_gates);
    }

    {
        void *new_malloc = malloc(circuit->depth[0] * circuit->depth[1] * sizeof(SoftGate *));
        if (new_malloc == NULL)
        {
            return result_get_invalid_reason("could not malloc");
        }
        memset(new_malloc, 0, sizeof(new_malloc));
        circuit->hardened_gates = (SoftGate **)new_malloc;
    }

    {
        Iter gates_iter = vector_iter_create(&circuit->soft_gates);
        Option next;
        while (ITER_NEXT(next, gates_iter))
        {
            SoftGate *soft_gate = (SoftGate *)next.data;
            unsigned int flat_position = soft_gate_flat_position(circuit, soft_gate);
            SoftGate **ptr_pos = circuit->hardened_gates + flat_position * sizeof(SoftGate *);
            if (ptr_pos != NULL)
            {
                free(circuit->hardened_gates);
                return result_get_invalid_reason("soft gates memory position collision");
            }
            *(ptr_pos) = soft_gate;
        }
    }

    return result_get_valid_with_data(circuit);
}

Result circuit_run(Circuit *circuit, double _Complex (*inout)[])
{
    if (circuit == NULL)
    {
        return result_get_invalid_reason("circuit pointer is null");
    }

    if (circuit->hardened_gates == NULL)
    {
        return result_get_invalid_reason("circuit hardened_gates is null");
    }

    for (unsigned int slice = 0; slice < circuit->depth[1]; ++slice)
    {
        Filter slice_gates = circuit_filter_slice_soft_gates(circuit, slice);

        double _Complex output[1 << (circuit->depth[0])];
        memset(&output, 0, sizeof(output));

        const unsigned int qubits = (circuit->depth[0]);
        const unsigned int gates_len = (circuit->soft_gates.size);
        unsigned int x, y, k, z;
        unsigned int mask, proj, cmask, ncmask, relevant, out_index;
        double _Complex coef;
        bool cset;

        k = 1 << (qubits - 1);

        mask = 0;
        cmask = 0;
        relevant = 0;
        for (unsigned int j = 0; j < gates_len; ++j)
        {
            SoftGate sgj = *(SoftGate *)result_unwrap(vector_get_raw(&circuit->soft_gates, j));
            mask |= 1 << (sgj.position.qubit);
            if ((cmask ^ mask) != 0)
            {
                relevant += 1;
            }
            cmask |= mask;
            if (sgj.control.some)
            {
                if ((cmask & (1 << sgj.control.data)) == 0)
                {
                    relevant += 1;
                }
                cmask |= 1 << sgj.control.data;
            }
        }
        ncmask = (~cmask) & (k | (k - 1));

        x = 0;
        for (unsigned int u = 0; u < (1 << relevant); ++u)
        {
            for (unsigned int y = 0; y < (1 << gates_len); ++y)
            {
                coef = (double _Complex)1;
                proj = 0;

                for (unsigned int j = 0; j < gates_len; ++j)
                {
                    SoftGate sgj = *(SoftGate *)result_unwrap(vector_get_raw(&circuit->soft_gates, j));
                    proj |= (y >> j) << sgj.position.qubit;
                    cset = sgj.control.some && (((x >> sgj.control.data) & 1) == 1);
                    if (((x >> sgj.position.qubit) & 1) == 1)
                    {
                        if (cset)
                        {
                            coef *= sgj.gate.matrix[1][(y >> j) & 1];
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
                            coef *= sgj.gate.matrix[0][(y >> j) & 1];
                        }
                        else
                        {
                            coef *= ((y >> j) & 1) ^ 1;
                        }
                    }
                }

                z = 0;
                for (unsigned int j = 0; j < (1 << (qubits - relevant)); ++j)
                {
                    out_index = (proj & mask) | ((z | x) & (~mask));
                    output[out_index] += coef * (*inout)[z | x];

                    z += 1;
                    while ((z & ncmask) != z)
                    {
                        unsigned int p = ~((z & -z) | ((z & -z) - 1));
                        z = (z & p) + ((ncmask & p) & (-(ncmask & p)));
                    }
                }

                x += 1;
                while ((x & cmask) != x)
                {
                    unsigned int p = ~((x & -x) | ((x & -x) - 1));
                    x = x & p + ((cmask & p) & (-(cmask & p)));
                }
            }
        }

        void *copy = memcpy(inout, output, sizeof(output));
        if (copy == NULL)
            return result_get_invalid_reason("memcpy failed");
    }

    return result_get_valid_with_data(inout);
}

Result circuit_free(Circuit *circuit)
{
    if (circuit == NULL)
    {
        return result_get_invalid_reason("circuit pointer is null");
    }

    vector_free(&circuit->soft_gates);
    vector_free(&circuit->slice_gate_count);

    if (circuit->hardened_gates != NULL)
        free(circuit->hardened_gates);
    free(circuit);
    circuit = NULL;

    return result_get_empty_valid();
}