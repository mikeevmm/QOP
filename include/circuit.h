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

#ifndef QOP_CIRCUIT_H_
#define QOP_CIRCUIT_H_

#include <complex.h>
#include "include/bits.h"
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
  Gate* gate;
} SoftGate;

// Information on the relevant qubit lines in
// each slice of the circuit
typedef struct CircuitSliceInfo {
  Vector slice_sg_ptrs;
  unsigned int gate_mask;
  unsigned int ctrl_gate_mask;
  unsigned int relevant_count;
  unsigned int ctrl_relevant_count;
} CircuitSliceInfo;

typedef struct Circuit {
  Vector soft_gates;
  Vector slice_gate_count;
  unsigned int depth[2];
  Vector slice_info_vec;
} Circuit;

// Gives a 1D index of the `SoftGate`'s position, by counting along each
// slice sequentially.
// For example, a gate positioned on the second qubit of the third slice
// of a five qubit circuit has a flat index of `2*5+1=11`, since the
// third slice has index `2`, and the second qubit has index `1`.
unsigned int soft_gate_flat_position(Circuit* circuit, SoftGate* soft_gate);

// Function for internal use.
// Used with a `Filter` (see `Filter.h`) to skip over NULL pointers in
// the hardened representation.
// Takes a pointer and determines whether it is a pointer to a NULL
// SoftGate pointer.
bool _circuit_filter_is_soft_gate(void* ptr);

// Initializes a given `Circuit` object.
// `qubit_count` is the number of "wires" in the circuit.
Result circuit_init(Circuit* circuit, unsigned int qubit_count);

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
Result circuit_compact(Circuit* circuit);

// Enumerates the gates in the circuit and creates the so-called
// "hardened representation", which allows O(1) access to all gates
// given their position, at the expense of O(q²) memory use (with q the
// number of qubits).
// This is needed for simulating the circuit.
Result circuit_harden(Circuit* circuit);

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

#endif  // QOP_CIRCUIT_H_