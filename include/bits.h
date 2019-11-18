/**
 * Declarations of all helper functions pertaining to bit
 * operations.
 * 
 * All functions should be constant and pure.
**/

#ifndef QOP_BITS_H_
#define QOP_BITS_H_

#include <stdio.h>
#include <stdbool.h>

bool arch_is_32bit(void);
bool arch_is_64bit(void);

// Print the value in binary
void bit_print(unsigned int value);

// Right-propagates the rightmost 1 in the value
unsigned int right_propagate(unsigned int value);

// Isolates the rightmost 1 in the value
unsigned int isolate_rightmost(unsigned int value);

// Isolates the leftmost 1 in the value
unsigned int round_to_next_pow2(unsigned int value);

// Yields the ith number `x` such that `x & mask = x` 
unsigned int ith_under_mask(unsigned int i, unsigned int mask);

#endif // QOP_BITS_H_