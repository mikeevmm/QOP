#pragma once
#include <option.h>
#include <stdbool.h>
#include <stdlib.h>

#define ITER_NEXT(NEXT, ITER) (NEXT = iter_next(&ITER)).some
#define FILTER_NEXT(NEXT, FILTER) (NEXT = filter_next(&FILTER)).some

typedef struct Iter
{
    void *head;
    unsigned int stride;
    unsigned int length;
    unsigned int position;
} Iter;

Iter iter_create(void *head, unsigned int stride, unsigned int length);
Option iter_next(Iter *iter);

typedef bool (*FilterFn)(void *);
typedef struct Filter
{
    Iter iter;
    FilterFn filter_fn;
    unsigned int position;
} Filter;

Filter filter_create(Iter iter, FilterFn filter_fn);
Option filter_next(Filter *filter);

bool filter_generic_not_null(void *ptr);