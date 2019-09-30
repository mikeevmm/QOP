#pragma once
#include "option.h"

#define ITER_NEXT(NEXT, ITER) (NEXT = iter_next(&ITER)).some

typedef struct Iter
{
    void *head;
    unsigned int stride;
    unsigned int position;
    unsigned int length;
} Iter;

Option iter_next(Iter *iter);