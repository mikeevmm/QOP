#include <stdio.h>
#include "include/dequeue.h"
#include "include/iter.h"
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
    dequeue_free(&dequeue);
  }

  printf("Reversing sequence with dequeue...\n");
  {
    Dequeue dequeue;
    result_unwrap(dequeue_init(&dequeue, sizeof(int)));
    int elements[] = {1, 2, 3};
    for (unsigned int i = 0; i < sizeof(elements) / sizeof(int); ++i) {
      result_unwrap(dequeue_push_front(&dequeue, elements + i));
    }
    for (unsigned int i = 0; i < sizeof(elements) / sizeof(int); ++i) {
      int peek = *(int *)result_unwrap(dequeue_peek_from_front(&dequeue, i));
      assert(peek == 3 - i);
    }
    for (unsigned int i = 0; i < sizeof(elements) / sizeof(int); ++i) {
      int popped;
      result_unwrap(dequeue_pop_back(&dequeue, &popped));
      assert(popped == i + 1);
    }
    dequeue_free(&dequeue);
  }

  printf("Iterating back to front and reverse...\n");
  {
    Dequeue dequeue;
    result_unwrap(dequeue_init(&dequeue, sizeof(int)));
    int elements[] = {1, 2, 3};
    for (unsigned int i = 0; i < sizeof(elements) / sizeof(int); ++i) {
      result_unwrap(dequeue_push_front(&dequeue, elements + i));
    }

    {
      Iter btf = dequeue_into_iterator_btf(&dequeue);
      Option next;
      int expected = 1;
      while ((next = iter_next(&btf)).some) {
        int value = *(int *)next.data;
        assert(value == expected);
        expected++;
      }
      iter_free(&btf);
    }

    {
      Iter ftb = dequeue_into_iterator_ftb(&dequeue);
      Option next;
      int expected = 3;
      while ((next = iter_next(&ftb)).some) {
        int value = *(int *)next.data;
        assert(value == expected);
        expected--;
      }
      iter_free(&ftb);
    }

    dequeue_free(&dequeue);
  }
}