#include <stdio.h>
#include "include/dequeue.h"
#include "include/option.h"

int main(void) {
    Dequeue dequeue;
    result_unwrap(dequeue_init(&dequeue, sizeof(int), 2));
    int element = 3;
    result_unwrap(dequeue_push_front(&dequeue, &element));
    result_unwrap(dequeue_pop_back(&dequeue, &element));
    if (element == 3) printf("Popped value matches pushed value.\n");
    else printf("Popped value does not match pushed value.\n");
}