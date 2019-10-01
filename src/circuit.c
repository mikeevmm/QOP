#include <circuit.h>

unsigned int soft_gate_flat_position(Circuit *circuit, SoftGate *soft_gate)
{
    return soft_gate->position.qubit + soft_gate->position.slice * circuit->depth[0];
}

bool _circuit_filter_is_soft_gate(void *ptr)
{
    SoftGate *sg_ptr = *((SoftGate **)ptr);
    return sg_ptr != NULL;
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

Result circuit_add_gate(Circuit *circuit, Gate *gate, unsigned int qubit, Option_Uint control)
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
    new_soft_gate.gate = *gate;
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

    unsigned int slice_gate_count[circuit->depth[1]];
    memset(slice_gate_count, 0, sizeof(slice_gate_count));

    unsigned int leftmost[circuit->depth[0]];
    memset(leftmost, 0, sizeof(leftmost));

    {
        Iter soft_gates_iter = vector_iter_create(&circuit->soft_gates);
        Option next;
        while ((next = iter_next(&soft_gates_iter)).some)
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
            slice_gate_count[new_position] += 1;
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

    {
        Result result;
        result = vector_clean(&circuit->slice_gate_count);
        if (!result.valid)
        {
            return result;
        }

        result = vector_extend(&circuit->slice_gate_count, slice_gate_count, circuit->depth[1]);
        if (!result.valid)
        {
            return result;
        }
    }

    if (circuit->hardened_gates != NULL)
    {
        free(circuit->hardened_gates);
    }

    {
        size_t malloc_size = circuit->depth[0] * circuit->depth[1] * sizeof(SoftGate *);
        void *new_malloc = malloc(malloc_size);
        if (new_malloc == NULL)
        {
            return result_get_invalid_reason("could not malloc");
        }
        memset(new_malloc, (SoftGate *)NULL, malloc_size/sizeof(SoftGate *));
        circuit->hardened_gates = (SoftGate **)new_malloc;
    }

    {
        Iter gates_iter = vector_iter_create(&circuit->soft_gates);
        Option next;
        while ((next = iter_next(&gates_iter)).some)
        {
            SoftGate *soft_gate = (SoftGate *)next.data;
            unsigned int flat_position = soft_gate_flat_position(circuit, soft_gate);
            SoftGate **ptr_pos = circuit->hardened_gates + flat_position;
            if (*ptr_pos != NULL)
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
        unsigned int qubits = (circuit->depth[0]);
        unsigned int gates_len = *(unsigned int *)(vector_get_raw(&circuit->slice_gate_count, slice).data);

        {
            unsigned long int head_offset = slice * circuit->depth[0] * sizeof(SoftGate *);
            Iter slice_gates_iter = iter_create((void *)((char *)circuit->hardened_gates + head_offset),
                                                sizeof(SoftGate *),
                                                circuit->depth[0]);
            Filter slice_gates_filter = filter_create(slice_gates_iter, _circuit_filter_is_soft_gate);
            Option next;
            while ((next = filter_next(&slice_gates_filter)).some)
            {
                SoftGate *sg = *((SoftGate **)next.data);
                printf("(%u, %u)\n", sg->position.qubit, sg->position.slice);
            }
        }

        double _Complex output[1 << (circuit->depth[0])];
        memset(&output, 0, sizeof(output));

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