#pragma once
#include "include/iter.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct Vector {
  bool init;
  void *data;
  size_t capacity;
  size_t size;
  size_t obj_size;
} Vector;

Result vector_create(Vector *v, size_t object_size, size_t init_capacity);
Result vector_get_raw(Vector *v, size_t index);
Result vector_resize(Vector *v, size_t size);
Result vector_push(Vector *v, void *object);
Result vector_raw_push(Vector *v, void *object);
Result vector_extend(Vector *v, void *object, size_t obj_count);
Result vector_extend_raw(Vector *v, void *object, size_t obj_count);
Result vector_pop(Vector *v);
Result vector_clean(Vector *v);
void vector_free(Vector *v);

Iter vector_iter_create(Vector *v);
Result filter_into_vector(Filter *filter, Vector *vector);
