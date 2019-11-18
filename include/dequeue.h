/**
 * Defines all the methods and structures needed for a double ended queue.
 **/
#ifndef QOP_DEQUEUE_H_
#define QOP_DEQUEUE_H_

#include <string.h>
#include "include/bits.h"
#include "include/option.h"

typedef struct Dequeue {
  Option head;
  unsigned int element_size;
  unsigned int head_index;
  unsigned int size;
  unsigned int capacity;
} Dequeue;

// Initializes a double ended queue.
// The initial capcity can be `0`, but this might result in a larger
// number of (re)allocations.
Result dequeue_init(Dequeue *dequeue, unsigned int element_size,
                    unsigned int initial_capacity);

// Reallocates memory internally so that dequeue can fit `fits`
// elements.
// Returns an empty Result if successful
Result dequeue_resize(Dequeue *dequeue, unsigned int fits);

// Pushes an element to the front of the double ended queue.
// The new element is `memcpy`'d into the dequeue's memory,
// so it can be freed from scope/memory safely after pushing.
// This may result in an internal reallocation of memory, if
// the queue's capacity needs to be increased.
// Returns Result(*element) if successful.
Result dequeue_push_front(Dequeue *dequeue, void *element);

// Pops an element from the back of the double ended queue.
// This may result in an internal reallocation of memory, if the
// capacity of the queue can be significantly decreased.
// The popped value is written into `into`.
// Returns `Result(*into)` if successful.
Result dequeue_pop_back(Dequeue *dequeue, void *into);

#endif  // QOP_DEQUEUE_H_