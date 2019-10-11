#include "include/circuit.h"

unsigned int soft_gate_flat_position(Circuit *circuit, SoftGate *soft_gate) {
  return soft_gate->position.qubit +
         soft_gate->position.slice * circuit->depth[0];
}

bool _circuit_filter_is_soft_gate(void *ptr) {
  SoftGate *sg_ptr = *((SoftGate **)ptr);
  return sg_ptr != NULL;
}

Result circuit_init(Circuit *circuit, unsigned int qubit_count) {
  if (qubit_count == 0) {
    return result_get_invalid_reason("cannot create a 0-qubit circuit");
  }

  // Initialize the vector to store `SoftGate`s belonging to the circuit,
  // pre-allocating space for a full slice.
  Vector soft_gates_vector;
  {
    Result vector_create_r =
        vector_init(&soft_gates_vector, sizeof(SoftGate), qubit_count);
    if (!vector_create_r.valid) {
      return vector_create_r;
    }
  }

  // Initialize the vector to store the gate count of each slice, which
  // is useful later, when compacting/running the circuit.
  // The vector is initialized with 0 capacity.
  Vector slice_gate_vector;
  {
    Result vector_create_r =
        vector_init(&slice_gate_vector, sizeof(unsigned int), 0);
    if (!vector_create_r.valid) {
      return vector_create_r;
    }
  }

  // Set circuit properties and return
  circuit->depth[0] = qubit_count;
  circuit->depth[1] = 0;
  circuit->soft_gates = soft_gates_vector;
  circuit->slice_gate_count = slice_gate_vector;
  circuit->hardened_gates = NULL;

  return result_get_valid_with_data(circuit);
}

// Makes a new `SoftGate` out of the given `Gate`, and adds that to the
// `Circuit`'s internal list. The circuit size is updated as needed.
// Returns `Result(circuit)`.
Result circuit_add_gate(Circuit *circuit, Gate *gate, unsigned int qubit,
                        Option_Uint control) {
  if (circuit == NULL) {
    return result_get_invalid_reason("circuit pointer is null");
  }

  if (qubit > circuit->depth[0] - 1)
    return result_get_invalid_reason("qubit is out of bounds");

  if (control.some) {
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
  if (!push_r.valid) {
    return push_r;
  }

  return result_get_valid_with_data(circuit);
}

// Compacts the circuit by pushing each gate as left as possible, while
// not allowing crossing of gates or gates with controls.
// This is done by keeping tabs on what the leftmost occupied position
// on each qubit line is, and then choosing for each gate the leftmost
// position available in the lines that it occupies.
// This is also the function responsible for checking how many gates are
// there in each slice and passing that to the `Circuit`'s
// `slice_gates_count`.
// This function does a single pass on the gates collection, with an
// extra cycle on the number of lines.
Result circuit_compact(Circuit *circuit) {
  if (circuit == NULL) {
    return result_get_invalid_reason("circuit pointer is null");
  }

  // Number of gates in each slice
  unsigned int slice_gate_count[circuit->depth[1]];
  memset(slice_gate_count, 0, sizeof(slice_gate_count));

  // Current leftmost position available
  unsigned int leftmost[circuit->depth[0]];
  memset(leftmost, 0, sizeof(leftmost));

  // Run over the `SoftGate`s and find the compacted position.
  // Update leftmost available positions accordingly.
  {
    Iter soft_gates_iter = vector_iter_create(&circuit->soft_gates);
    Option next;
    while ((next = iter_next(&soft_gates_iter)).some) {
      SoftGate *soft_gate = (SoftGate *)next.data;
      unsigned int new_position = 0;

      if (soft_gate->control.some) {  // Check gate and control
        unsigned int control = soft_gate->control.data;
        if (leftmost[soft_gate->position.qubit] > leftmost[control]) {
          new_position = leftmost[soft_gate->position.qubit];
        } else {
          new_position = leftmost[control];
        }

        leftmost[soft_gate->position.qubit] = new_position + 1;
        leftmost[control] = new_position + 1;
      } else {  // Just check gate
        new_position = leftmost[soft_gate->position.qubit];
        leftmost[soft_gate->position.qubit] += 1;
      }

      soft_gate->position.slice = new_position;
      slice_gate_count[new_position] += 1;
    }
  }

  // Update depth
  {
    unsigned int max = 0;
    for (unsigned int i = 0; i < circuit->depth[0]; i++) {
      if (leftmost[i] > max) max = leftmost[i];
    }
    circuit->depth[1] = max;
  }

  // Copy the new slice gate counts over to the `Circuit` object
  {
    Result result;
    result = vector_clean(&circuit->slice_gate_count);
    if (!result.valid) {
      return result;
    }

    result = vector_extend(&circuit->slice_gate_count, slice_gate_count,
                           circuit->depth[1]);
    if (!result.valid) {
      return result;
    }
  }

  return result_get_valid_with_data(circuit);
}

// Creates the "hardened representation" by allocating enough memory for
// it (known after compacting the circuit).
// The hardened representation is just a memory block of `SoftGate *`
// where the pointer to the gate at position (i,j) is located at
// `<hardened gate block head> + <flat position of (i, j)>`.
// This function does a single pass on the gates collection.
Result circuit_harden(Circuit *circuit) {
  // Create the hardened gate representation
  // Flat positions without a gate are left pointing to NULL; it's the
  // user's responsability to check whether they're pointing to a valid
  // location.
  if (circuit->hardened_gates != NULL) {
    free(circuit->hardened_gates);
  }

  // Allocate appropriate size...
  {
    size_t malloc_size =
        circuit->depth[0] * circuit->depth[1] * sizeof(SoftGate *);
    void *new_malloc = malloc(malloc_size);
    if (new_malloc == NULL) {
      return result_get_invalid_reason("could not malloc");
    }
    memset(new_malloc, 0, malloc_size);
    circuit->hardened_gates = (SoftGate **)new_malloc;
  }

  // Iterate over the gates and copy pointers to appropriate location
  {
    Iter gates_iter = vector_iter_create(&circuit->soft_gates);
    Option next;
    while ((next = iter_next(&gates_iter)).some) {
      SoftGate *soft_gate = (SoftGate *)next.data;
      unsigned int flat_position = soft_gate_flat_position(circuit, soft_gate);
      SoftGate **ptr_pos = circuit->hardened_gates + flat_position;
      if (*(int *)ptr_pos != 0) {
        free(circuit->hardened_gates);
        return result_get_invalid_reason(
            "soft gates memory position collision");
      }
      *(ptr_pos) = soft_gate;
    }
  }

  return result_get_valid_with_data(circuit);
}

// Performs a simulation of an input.
// The simulation can be thought of as a reduction over the slices;
//
// ```ASCII
//                           +-(when layers exhausted)->[Final output]
//                           |
//                     +-----+----+
//                     |  Slice   |   <[as input to next slice]
// [input array] --->--+          +-<-+
//   (*inout)          |Simulation|   |
//                     +-----+----+   |
//                           |        |
//                           +---->---+
//                   [slice output]>
// ```
//
// (This is also why `inout` is modified in place.)
// As such, only the slice simulation subroutine is discussed below.
// The slice simulation subroutine (henceforth referred to as the SSS,
// for short) operates on the superposition principle, whereas each
// component of the input can be taken separately, as long as the
// different wavefunction coefficients are added together.
// On the other hand, the SSS attempts to exploit the fact that qubits
// which are not affected by gates nor are control targets do not affect
// the wavefunction coefficient, and so a "batch write" of the coefficient
// can be made to all permutations of irrelevant qubits.
// Finally, the SSS exploits the "meaning" of matrix multiplication, in
// particular that the matrix corresponding to a gate can be thought of
// in terms of projectors:
//
//       <0| <1|
//     +-       -+
// |0> |  a   b  |
//     |         |
// |1> |  c   d  |
//     +-       -+
//
// Pseudocode of the SSS follows:
//
// 1. For all permutations of a string as big as the number of gates in
//    the slice
//                x := perm(string of len gates_len)
//    corresponding to the "input" in only the gates' qubits
// 2. For each x, consider also all permutations
//                y := perm(string of len gates_len)
//    corresponding to the "output" in the gates' qubits
// 3. For each y, calculate the coefficient c of the output by
// 3.1. Letting
//              c := 1.
// 3.2. For each gate in the slice,
// 3.2.1. Check if the gate is *not* controlled, or if its control is
//        set in the corresponding bit of x
// 3.2.2. If so, go to 3.2.2a., else go to 3.2.2b.
// 3.2.2a.  Multiply c by the element of the gate matrix M corresponding
//          to the relevant bits of x and y, x' and y' respectively:
//                c *= M[y', x']
// 3.2.2b.  Multiply c by the element of the 2-identity matrix
//          corresponding to the relevant bits of x and y, x' and y'
//          respectively:
//                c *= (x' XOR 1) XOR y'
// 3.3. Take all bitstrings of appropriate size (2**circuit qubit count)
//      whose relevant bits match x
//                z := bitstring matching x
// 3.4. Let the bitstring o be the bitstring whose irrelevant bits match
//      z, but whose relevant bits have been set to match y
//                o := z project y
// 3.5 Set the relevant output vector component as
//                output[o] += c * input[z]
//
//
// In order to consider all bitstrings varying in some bitmask, the
// following technique was employed:
//    bitstring := b
//    x := 0
//    x += 1
//    while ( x AND b ) != x    (meaning that there are digits outside
//                               the intended mask region)
//        Take the rightmost 1-bit of x and rightpropagate and invert
//        the value. Let this value be p.
//        Find what's the rightmost 1-bit of b  that is whithin p,
//        let this be delta.
//        Let
//            x := (x & p) + delta
//
// What this does, in a nutshell, is keep adding the rightmost offending
// bit to the value of x until it's within the intended mask.
Result circuit_run(Circuit *circuit, double _Complex (*inout)[]) {
  if (circuit == NULL) {
    return result_get_invalid_reason("circuit pointer is null");
  }

  if (circuit->hardened_gates == NULL) {
    return result_get_invalid_reason("circuit hardened_gates is null");
  }

  unsigned int qubits = (circuit->depth[0]);
  unsigned int leftmost = 1U << (qubits - 1U);
  Vector slice_gates;
  vector_init(&slice_gates, sizeof(SoftGate *), 0);

  // Find how many/what the gates for this slice are
  for (unsigned int slice = 0; slice < circuit->depth[1]; ++slice) {
    unsigned int gates_len =
        *(unsigned int *)(vector_get_raw(&circuit->slice_gate_count, slice));

    {
      unsigned long int head_offset =
          slice * circuit->depth[0] * sizeof(SoftGate *);
      Iter slice_gates_iter =
          iter_create((void *)((char *)circuit->hardened_gates + head_offset),
                      sizeof(SoftGate *), circuit->depth[0]);
      Filter slice_gates_filter =
          filter_create(slice_gates_iter, _circuit_filter_is_soft_gate);
      filter_into_vector(&slice_gates_filter, &slice_gates);
    }

    // Output vector to be written to during computation; its value
    // is transferred to inout after a slice cycle
    double _Complex output[1 << qubits];
    memset(output, 0, sizeof(output));

    // Bitmasks:
    //    Relevant bits (gates only)
    unsigned int mask = 0;
    //    Relevant bits (gates and controls)
    unsigned int cmask = 0;
    //    Irrelevant bits (not gates, not controls)
    unsigned int ncmask = 0;

    // Count of relevant bits
    unsigned int relevant = 0;

    // Determine the bitmasks + `relevant`
    for (unsigned int i = 0; i < gates_len; ++i) {
      unsigned int gate_i_qubit;
      Option_Uint gate_i_control;
      {
        SoftGate sg = **((SoftGate **)vector_get_raw(&slice_gates, i));
        gate_i_qubit = sg.position.qubit;
        gate_i_control = sg.control;
      }

      mask |= (1U << gate_i_qubit);
      if ((cmask ^ mask) != 0U) relevant += 1;
      cmask |= mask;
      if (gate_i_control.some) {
        if ((cmask & (1U << gate_i_control.data)) == 0) relevant += 1;
        cmask |= 1U << gate_i_control.data;
      }
    }
    ncmask = (~cmask) & (leftmost | (leftmost - 1));

    // Input string (see 1. above)
    unsigned int x = 0;
    for (unsigned int i = 0; i < (1U << relevant); ++i) {
      // Relevant output bits (see 2. above)
      for (unsigned int y = 0; y < (1U << gates_len); ++y) {
        double _Complex coef = 1.;
        unsigned int proj = 0;

        for (unsigned int j = 0; j < gates_len; ++j) {
          SoftGate gate_j = **(SoftGate **)(vector_get_raw(&slice_gates, j));

          proj |= (y >> j) << gate_j.position.qubit;
          // If there is no control, then the gate is always set,
          // therefore the ||
          bool cset =
              (!gate_j.control.some) || ((x >> gate_j.control.data) & 1U);

          unsigned int relevant_x_bit = (x >> gate_j.position.qubit) & 1U;
          unsigned int relevant_y_bit = (y >> j) & 1U;
          if (cset)
            coef *= gate_j.gate->matrix[relevant_y_bit][relevant_x_bit];
          else
            coef *= (double _Complex)(1U ^ relevant_x_bit ^ relevant_y_bit);
        }

        // Output string (see 3.4. above)
        unsigned int z = 0;
        for (unsigned int j = 0; j < (1U << (qubits - relevant)); ++j) {
          unsigned int out_index = (proj & mask) | ((z | x) & (~mask));
          output[out_index] += coef * (*inout)[z | x];

          // Move to the next number z in the "relevant qubit number space";
          // it's like counting using only some fingers!
          // See above for a detailed explanation
          z += 1;
          while ((z & ncmask) != z) {
            unsigned int p = ~((z & (-z)) | ((z & (-z)) - 1));
            z = (z & p) + ((ncmask & p) & (-(ncmask & p)));
          }
        }
      }

      // Move to the next relevant x
      // See above for a detailed explanation
      x = x + 1;
      while ((x & cmask) != x) {
        unsigned int p = ~((x & (-x)) | ((x & (-x)) - 1));
        x = (x & p) + ((cmask & p) & (-(cmask & p)));
      }
    }

    vector_clean(&slice_gates);

    void *copy =
        memcpy(*inout, output, sizeof(double _Complex) * (1 << qubits));
    if (copy == NULL) return result_get_invalid_reason("memcpy failed");
  }

  vector_free(&slice_gates);

  return result_get_valid_with_data(inout);
}

// Frees all the heap memory allocated by the circuit, but not the
// circuit itself.
Result circuit_free(Circuit *circuit) {
  if (circuit == NULL) {
    return result_get_invalid_reason("circuit pointer is null");
  }

  vector_free(&circuit->soft_gates);
  vector_free(&circuit->slice_gate_count);

  if (circuit->hardened_gates != NULL) free(circuit->hardened_gates);

  return result_get_empty_valid();
}