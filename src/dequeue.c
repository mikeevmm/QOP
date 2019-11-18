#include "include/dequeue.h"

Result dequeue_init(Dequeue *dequeue, unsigned int element_size,
                    unsigned int initial_capacity) {
  // Very basic sanitizing
  if (dequeue == NULL)
    return result_get_invalid_reason("given *dequeue is NULL");

  // Alocate memory to be the storage for the queue;
  // If 0 was specified, skip this step
  Option allocated;
  if (initial_capacity > 0) {
    void *alloc = malloc(initial_capacity * element_size);
    if (alloc == NULL) return result_get_invalid_reason("could not malloc");
    allocated = option_some_with_data(alloc);
  } else {
    allocated = option_none();
  }

  // Initialize the struct
  dequeue->head = allocated;
  dequeue->element_size = element_size;
  dequeue->head_index = 0;
  dequeue->size = 0;
  dequeue->capacity = initial_capacity;

  return result_get_empty_valid();
}

Result dequeue_resize(Dequeue *dequeue, unsigned int fits) {
  if (fits < dequeue->size) {
    return result_get_invalid_reason("cannot resize to smaller than contents");
  }

  // Resize the capacity to the smallest power of two bigger
  // than `fits`
  unsigned int new_capacity = round_to_next_pow2(fits);
  if (new_capacity == dequeue->capacity) return result_get_empty_valid();

  if (new_capacity > 0) {
    void *new_mem = malloc(dequeue->element_size * new_capacity);
    if (new_mem == NULL) return result_get_invalid_reason("could not malloc");

    // Allocate or reallocate memory
    if (dequeue->head.some) {
      // Reallocate
      // Since we are reallocating, sort the queue
      // Copy the queue in order
      char *actual_head = (char *)(dequeue->head.data) +
                          dequeue->head_index * dequeue->element_size;
      // Copy from `actual_head` to end
      {
        void *cpy_result = memcpy(
            new_mem, (void *)actual_head,
            (dequeue->size - dequeue->head_index) * dequeue->element_size);
        if (cpy_result == NULL) {
          free(new_mem);
          return result_get_invalid_reason("could not memcpy");
        }
      }
      // Copy from end to `actual_head`
      {
        void *cpy_result = memcpy(new_mem, dequeue->head.data,
                                  dequeue->head_index * dequeue->element_size);
        if (cpy_result == NULL) {
          free(new_mem);
          return result_get_invalid_reason("could not memcpy");
        }
      }
    }  // else There was no allocated memory; no sorting/copying

    dequeue->head = option_some_with_data(new_mem);
    dequeue->head_index = 0;
    dequeue->capacity = new_capacity;
  } else {  // Capacity == 0
    if (dequeue->head.some) free(dequeue->head.data);

    dequeue->head = option_none();
    dequeue->head_index = 0;
    dequeue->capacity = 0;
  }

  return result_get_empty_valid();
}

Result dequeue_push_front(Dequeue *dequeue, void *element) {
  // Very basic sanitizing
  if (dequeue == NULL) return result_get_invalid_reason("*dequeue is NULL");
  if (element == NULL) return result_get_invalid_reason("*element is NULL");

  // Check the size of the dequeue; resize if needed
  if (dequeue->size + 1 > dequeue->capacity) {
    Result resize_result = dequeue_resize(dequeue, dequeue->size + 1);
    if (!resize_result.valid) return resize_result;
  }

  // Push element; need to determine actual "physical" position
  unsigned int push_index =
      (dequeue->head_index + dequeue->size) % dequeue->capacity;
  char *push_pos =
      (char *)(dequeue->head.data) + push_index * dequeue->element_size;
  memcpy((void *)push_pos, element, dequeue->element_size);

  dequeue->size += 1;

  return result_get_empty_valid();
}

Result dequeue_pop_back(Dequeue *dequeue, void *into) {
  // Basic checks
  if (dequeue == NULL) return result_get_invalid_reason("dequeue is NULL");
  if (dequeue->size == 0) return result_get_invalid_reason("dequeue is empty");

  // Memcpy relevant element
  char *elem_pos = (char *)(dequeue->head.data) +
                   dequeue->head_index * dequeue->element_size;
  memcpy(into, (void *)elem_pos, dequeue->element_size);

  // Offset head_index; decrease size
  // There's no need to "delete" the popped element; it will be freed or
  // overwritten when appropriate
  dequeue->head_index = (dequeue->head_index + 1) % dequeue->capacity;
  dequeue->size -= 1;

  // Shrink if possible
  Result resize_result = dequeue_resize(dequeue, dequeue->size);
  if (!resize_result.valid) return resize_result;

  return result_get_empty_valid();
}
