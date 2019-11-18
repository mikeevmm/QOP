#include <stdio.h>
#include "include/dequeue.h"
#include "include/option.h"

int main(void) {
    printf("Pushing/popping single value...\n");
    {
        Dequeue dequeue;
        result_unwrap(dequeue_init(&dequeue, sizeof(int)));
        int element = 3;
        result_unwrap(dequeue_push_front(&dequeue, &element));
        result_unwrap(dequeue_pop_back(&dequeue, &element));
        assert(element == 3);
    }

    printf("Reversing sequence with dequeue...\n");
    {
        Dequeue dequeue;
        result_unwrap(dequeue_init(&dequeue, sizeof(int)));
        int elements[] = {1, 2, 3};
        for(unsigned int i = 0; i < sizeof(elements)/sizeof(int); ++i) {
            result_unwrap(dequeue_push_front(&dequeue, elements + i));
        }
        for(unsigned int i = 0; i < sizeof(elements)/sizeof(int); ++i) {
            int popped;
            result_unwrap(dequeue_pop_back(&dequeue, &popped));
            assert(popped == i + 1);
        }
    }
}