#include "vector.h"

Result vector_create(Vector *v, size_t object_size, size_t init_capacity)
{
    if (init_capacity > 0)
    {
        v->data = (void *)malloc(init_capacity * object_size);
        if (!v->data)
            return result_get_invalid_reason("could not malloc");
    }
    else
    {
        v->data = NULL;
    }

    v->size = 0;
    v->capacity = init_capacity;
    v->obj_size = object_size;

    return result_get_valid_with_data(v);
}

Result vector_get_raw(Vector *v, size_t index)
{
    void *object = v->data + index * v->obj_size;
    return result_get_valid_with_data(object);
}

Result vector_resize(Vector *v, size_t size)
{
    if (size == v->capacity)
    {
        return result_get_valid_with_data(v);
    }
    else if (size > v->capacity)
    {
        // Resize to next power of 2
        size_t new_capacity = (v->capacity == 0 ? 1 : v->capacity);
        while (new_capacity < size)
            new_capacity <<= 1;

        void *resized = realloc(v->data, new_capacity * v->obj_size);
        if (!resized)
        {
            return result_get_invalid_reason("could not resize vector; realloc failed");
        }
        v->data = resized;
        v->capacity = new_capacity;
        return result_get_valid_with_data(v);
    }
    else if (size <= v->capacity / 2)
    {
        size_t new_capacity = (v->capacity == 0 ? 1 : v->capacity);
        while (new_capacity / 2 > size)
            new_capacity >>= 1;

        void *resized = realloc(v->data, new_capacity * v->obj_size);
        if (!resized)
        {
            return result_get_invalid_reason("could not resize vector; realloc failed");
        }
        v->data = resized;
        v->capacity = new_capacity;
        return result_get_valid_with_data(v);
    }
}

Result vector_push(Vector *v, void *object)
{
    if (v->size + 1 > v->capacity)
    {
        Result resize = vector_resize(v, v->size + 1);
        if (!resize.valid)
        {
            return resize;
        }
    }
    return vector_raw_push(v, object);
}

Result vector_raw_push(Vector *v, void *object)
{
    void *moved = memcpy(v->data + v->obj_size * v->size, object, v->obj_size);
    if (!moved)
    {
        return result_get_invalid_reason("could not memcpy");
    }
    v->size += 1;
    return result_get_valid_with_data(v);
}

Result vector_extend(Vector *v, void *object, size_t obj_count)
{
    if (v->size + obj_count > v->capacity)
    {
        Result resize = vector_resize(v, v->size + obj_count);
        if (!resize.valid)
        {
            return resize;
        }
    }
    return vector_extend_raw(v, object, obj_count);
}

Result vector_extend_raw(Vector *v, void *object, size_t obj_count)
{
    void *moved = memcpy(v->data + v->obj_size * v->size, object, v->obj_size * obj_count);
    if (!moved)
    {
        return result_get_invalid_reason("failed to memcpy");
    }
    v->size += obj_count;
    return result_get_valid_with_data(v);
}

Result vector_pop(Vector *v, void *object)
{
    void *poped_loc = (void *)malloc(v->obj_size);
    if (!poped_loc)
    {
        return result_get_invalid_reason("could not malloc");
    }
    void *moved = memcpy(poped_loc, v->data + v->size * v->obj_size, v->obj_size);
    if (!moved)
    {
        return result_get_invalid_reason("could not memmove");
    }
    free(v->data + v->size * v->obj_size);
    v->size -= 1;
    vector_resize(v, v->size);

    return result_get_valid_with_data(v);
}

void vector_free(Vector *v)
{
    free(v->data);
}

Iter vector_iter_create(Vector *v)
{
    Iter new_iter;
    new_iter.head = v->data;
    new_iter.length = v->size;
    new_iter.stride = v->obj_size;
    new_iter.position = 0;
    return new_iter;
}