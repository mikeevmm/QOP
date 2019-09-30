#include "iter.h"

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