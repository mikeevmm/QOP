#pragma once
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Result {
  bool valid;
  const char *reason;
  const char *file;
  unsigned int line;
  void *data;
} Result;

Result result_get_empty_valid();
Result result_get_valid_with_data(void *data);
Result result_get_invalid_reason_raw(const char *reason, const char *file,
                                     unsigned int line);
#define result_get_invalid_reason(reason) \
  result_get_invalid_reason_raw(reason, __FILE__, __LINE__)
void *result_unwrap(Result result);

typedef struct Option {
  bool some;
  void *data;
} Option;

typedef struct Option_Uint {
  bool some;
  unsigned int data;
} Option_Uint;

typedef struct Option_Int {
  bool some;
  int data;
} Option_Int;

typedef struct Option_Double {
  bool some;
  double data;
} Option_Double;

Option option_none();
Option_Int option_none_int();
Option_Uint option_none_uint();
Option_Double option_none_double();
Option option_some_with_data(void *data);
Option_Int option_from_int(int data);
Option_Uint option_from_uint(unsigned int data);
Option_Double option_from_double(double data);