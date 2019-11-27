#include <assert.h>
#include <stdio.h>
#include "include/iter.h"
#include "include/option.h"

int main(void) {
  unsigned int test_array[5] = {1, 2, 3, 4, 5};

  Iter array_iter = iter_create_contiguous_memory((void *)test_array,
                                                  sizeof(unsigned int), 5);
  Option next;
  unsigned int expected = 1;
  while ((next = iter_next(&array_iter)).some) {
    unsigned int *value = (unsigned int *)next.data;
    assert(*value == expected);
    expected++;
  }
}