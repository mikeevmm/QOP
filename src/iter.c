#include "iter.h"

Iter iter_create(void *head, unsigned int stride, unsigned int length)
{
    Iter new_iter;
    new_iter.position = 0;
    new_iter.head = head;
    new_iter.stride = stride;
    new_iter.length = length;
    return new_iter;
}

Option iter_next(Iter *iter)
{
    if (iter->position == iter->length)
    {
        return option_none();
    }

    void *next = iter->head + iter->stride * iter->position;
    iter->position += 1;
    return option_some_with_data(next);
}

Filter filter_create(Iter iter, FilterFn filter_fn)
{
    Filter new_filter;
    new_filter.iter = iter;
    new_filter.filter_fn = filter_fn;
    new_filter.position = 0;
    return new_filter;
}

Option filter_next(Filter *filter)
{
    Option next = iter_next(&filter->iter);
    while (next.some && !filter->filter_fn(next.data))
    {
        next = iter_next(&filter->iter);
    }
    filter->position += 1;
    return next;
}

bool filter_generic_not_null(void *ptr)
{
    return ptr != NULL;
}