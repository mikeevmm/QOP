#include "vector.h"

Result vector_init(Vector *v, size_t object_size, size_t init_capacity)
{
    v->data = malloc(init_capacity * object_size);
    if (!v->data)
        return result_get_invalid_reason("could not malloc");

    v->size = 0;
    v->capacity = init_capacity;
    v->obj_size = object_size;

    return result_get_empty_valid();
}

Result vector_resize(Vector *v, size_t size)
{
    if (size == v->capacity)
    {
        return result_get_empty_valid();
    }
    else if (size > v->capacity)
    {
        // Resize to next power of 2
        size_t new_capacity = 1;
        {
            size_t work = size;
            while (work >>= 1)
                new_capacity <<= 1;
            new_capacity <<= 1;
        }

        void *resized = realloc(v->data, new_capacity * v->obj_size);
        if (!resized)
        {
            return result_get_invalid_reason("could not resize vector; realloc failed");
        }
        v->data = resized;
        v->capacity = new_capacity;
        return result_get_empty_valid();
    }
    else if (size <= v->capacity / 2)
    {
        size_t new_capacity = v->capacity;
        while (new_capacity / 2 > size)
            new_capacity >>= 1;

        void *resized = realloc(v->data, new_capacity * v->obj_size);
        if (!resized)
        {
            return result_get_invalid_reason("could not resize vector; realloc failed");
        }
        v->data = resized;
        v->capacity = new_capacity;
        return result_get_empty_valid();
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
    void *copied = memcpy(v->data + v->obj_size * v->size, object, v->obj_size);
    if (!copied)
    {
        return result_get_invalid_reason("could not memcpy");
    }
    v->size += 1;
    return result_get_empty_valid();
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
    void *copied = memcpy(v->data + v->obj_size * v->size, object, v->obj_size * obj_count);
    if (!copied)
    {
        return result_get_invalid_reason("failed to memcpy");
    }
    v->size += obj_count;
    return result_get_empty_valid();
}

Result vector_pop(Vector *v, void *object)
{
    void *poped_loc = malloc(v->obj_size);
    if (!poped_loc)
    {
        return result_get_invalid_reason("could not malloc");
    }
    void *moved = memmove(poped_loc, v->data + v->size * v->obj_size, v->obj_size);
    if (!moved)
    {
        return result_get_invalid_reason("could not memmove");
    }
    v->size -= 1;
    vector_resize(v, v->size);
}
