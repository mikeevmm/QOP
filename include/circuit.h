/**
 * Declaration of all structures and methods used to represent a
 * circuit, as defined by the gates it contains, their positions, and a
 * few reflexive properties. The `SoftGate` structure, here defined,
 * defines a gate in a circuit, as opposed to a gate as an isolated
 * entity.
 *  
 * All methods that operate on a circuit (here declared) are prefixed by
 * `circuit_`.
 **/

#pragma once
#include "include/gate.h"
#include "include/iter.h"
#include "include/option.h"
#include "include/vector.h"

// Position of a `SoftGate` in a circuit.
// `slice` is the horizontal position (left-to-right), and `qubit` is
// the vertical position.
typedef struct GatePosition {
  unsigned int slice;
  unsigned int qubit;
} GatePosition;

// Representation of a gate in a circuit.
typedef struct SoftGate {
  GatePosition position;
  Option_Uint control;
  Gate *gate;
} SoftGate;

typedef struct Circuit {
  Vector soft_gates;
  Vector slice_gate_count;
  unsigned int depth[2];
  SoftGate** hardened_gates;
} Circuit;

// Gives a 1D index of the `SoftGate`'s position, by counting along each
// slice sequentially.
// For example, a gate positioned on the second qubit of the third slice
// of a five qubit circuit has a flat index of `2*5+1=11`, since the
// third slice has index `2`, and the second qubit has index `1`.
unsigned int soft_gate_flat_position(Circuit *circuit, SoftGate *soft_gate);

// Function for internal use.
// Used with a `Filter` (see `Filter.h`) to skip over NULL pointers in
// the hardened representation.
// Takes a pointer and determines whether it is a pointer to a NULL
// SoftGate pointer.
bool _circuit_filter_is_soft_gate(void* ptr);

// Creates a `Circuit` object on the heap and returns a pointer to that
// memory block.
// `qubit_count` is the number of "wires" in the circuit.
Result circuit_create(unsigned int qubit_count);

// Unwraps the result of a `circuit_create` call; attempts to unwrap
// the given `Result`; if successful, frees the heap memory of the data,
// returning a stack copy of it.
Circuit circuit_unwrap_create(Result result);

// Adds a previously initialized `Gate` to the circuit, on line `qubit`,
// and optionally controlled by another line, as specified by `control`.
// Although a `SoftGate` is created internally from these arguments, the
// original gate is not copied into the `SoftGate`, but rather a pointer
// reference is kept. It is the user's responsability to ensure that the
// `Gate` object outlives the `Circuit` object, otherwise the behaviour
// is unspecified.
// The new gate's horizontal position is to the left of any previously
// added gates; a call to `circuit_compact` should be made afterwards,
// to reduce the horizontal depth of the circuit.
Result circuit_add_gate(Circuit* circuit, Gate* gate, unsigned int qubit,
                        Option_Uint control);

// Readjusts all `SoftGates`' positions in the circuit so that every
// gate is as leftmost as possible, while preserving order of addition
// of gates. For example:
//
// --[]-------o---[]--             --[]-o-[]--
// ------[]---|-------   becomes   --[]-|-----
// ----------[ ]------             ----[ ]----
//
// This function is also responsible for enumerating the gates in the
// circuit and creating the so-called "hardened representation",
// which is needed for simulating the circuit.
Result circuit_compact(Circuit* circuit);

// Calculates the simulated output state of the circuit, when given some
// input state as determined by `*inout`. The input state `|i>` is
// determined by `*inout` such that 
//    `|i> = sum( (*inout)[k] |k> for k in 0..2**q )`,
// with `q` the qubit count of the circuit, and `|k>` the `k`th basis
// state. 
// `*inout` is used as a buffer during computation, so that it should be
// passed in by reference, and the method is destructive under `*inout`.
// The final output is stored in `(*inout)[]`. 
Result circuit_run(Circuit* circuit, double _Complex (*inout)[]);

// Frees the members of the `Circuit` that were allocated at
// initialization or during a method (in particular the heap memory for
// `hardened_gates`, allocated in `circuit_harden`).
// This function DOES NOT free objects which it did not allocate (and
// therefore can't guarantee are not in the stack; in particular,
// `Gate`s added to this circuit are not freed) and does not free the
// `*circuit` pointer itself (because, again, it can be in the stack).
Result circuit_free(Circuit* circuit);