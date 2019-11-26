/**
 * Defines all the methods and structures needed for a double ended queue.
 **/
#ifndef QOP_DEQUEUE_H_
#define QOP_DEQUEUE_H_

#include <string.h>
#include "include/bits.h"
#include "include/iter.h"
#include "include/option.h"

typedef struct Dequeue {
  Option head;
  long unsigned int element_size;
  unsigned int head_index;
  unsigned int size;
  unsigned int capacity;
} Dequeue;

// Initializes a double ended queue.
Result dequeue_init(Dequeue *dequeue, long unsigned int element_size);

// Frees any internal memory that the dequeue may have allocated.
void dequeue_free(Dequeue *dequeue);

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
// The popped value is written into `into`. If `into` is NULL, the
// value is just discarded.
// Returns `Result(*into)` if successful.
Result dequeue_pop_back(Dequeue *dequeue, void *into);

// Creates an iterator over the elements of a double ended queue,
// from the oldest element to the most recent one.
Iter dequeue_into_iterator_btf(Dequeue *dequeue);

// Creates an iterator over the elements of a double ended queue,
// from the most recent element to the oldest one.
Iter dequeue_into_iterator_ftb(Dequeue *dequeue);

// Returns a pointer to the `back`th element of the dequeue, counting
// from the front (latest added) of the queue.
// This function performs boundary checks; if performance is critical,
// maybe consider working with the struct's internals directly.
Result dequeue_peek_from_front(Dequeue *dequeue, unsigned int back);

#endif  // QOP_DEQUEUE_H_