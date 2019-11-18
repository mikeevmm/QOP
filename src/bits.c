#include "include/bits.h"

bool arch_is_32bit(void) {
  return sizeof(int *) == 4;
}

bool arch_is_64bit(void) {
  return sizeof(int *) == 8;
}

void bit_print(unsigned int value) {
  char bits[sizeof(unsigned int) * 8 + 1];
  for (unsigned long int i = 0; i < sizeof(unsigned int) * 8; ++i) {
    unsigned long int k = sizeof(unsigned int) * 8 - i - 1;
    bits[i] = (char)(48UL + (((1U << k) & value) >> k));
  }
  bits[sizeof(unsigned int) * 8] = 0;
  puts(bits);
}

unsigned int right_propagate(unsigned int value) {
  if (value == 0) return 0;
  return value | (value - 1);
}

unsigned int isolate_rightmost(unsigned int value) { return value & (-value); }

unsigned int round_to_next_pow2(unsigned int value) {
  value--;
  value |= value >> 1;
  value |= value >> 2;
  value |= value >> 4;
  value |= value >> 8;
  value |= value >> 16;
  if (arch_is_64bit()) value |= value >> 32;
  value++;
  return value;
}

unsigned int ith_under_mask(unsigned int i, unsigned int mask) {
  // Run through the bits of the mask; add offsets as appropriate
  unsigned int value = 0;
  unsigned int j = 0;  // Bit offset over i
  for (unsigned int k = 0; k < (sizeof(unsigned int)*8); ++k) {
    if ((1U << k) > mask) break;
    if ((1U << j) > i) break;
    if (!((mask >> k) & 1)) continue;
    value += ((i >> j) & 1U) << k;
    j += 1;
  }
  return value;
}
