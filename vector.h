#pragma once
#include "option.h"
#include <stdlib.h>
#include <string.h>

typedef struct Vector
{
    void *data;
    size_t capacity;
    size_t size;
    size_t obj_size;
} Vector;

typedef struct Vector_Iter
{
    Vector *vector;
    size_t position;
} Vector_Iter;

Result vector_init(Vector *v, size_t object_size, size_t init_capacity);
Result vector_resize(Vector *v, size_t size);
Result vector_push(Vector *v, void *object);
Result vector_raw_push(Vector *v, void *object);
Result vector_extend(Vector *v, void *object, size_t obj_count);
Result vector_extend_raw(Vector *v, void *object, size_t obj_count);
Result vector_pop(Vector *v, void *object);
Result vector_free(Vector *v);

Vector_Iter vector_iter_create(Vector *v);
Option vector_iter_next(Vector_Iter *vi);