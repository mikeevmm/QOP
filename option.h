#pragma once
#include <stdbool.h>
#include <stdlib.h>

typedef struct Result
{
    bool valid;
    char *reason;
    void *data;
} Result;

Result result_get_empty_valid();
Result result_get_valid_with_data(void *data);
Result result_get_invalid_reason(char *reason);

typedef struct Option
{
    bool some;
    void *data;
} Option;

typedef struct Option_Uint
{
    bool some;
    unsigned int data;
} Option_Uint;


typedef struct Option_Int
{
    bool some;
    int data;
} Option_Int;


typedef struct Option_Double
{
    bool some;
    double data;
} Option_Double;

Option option_none();
Option option_some_with_data(void *data);
Option_Int option_from_int(int data);
Option_Uint option_from_uint(unsigned int data);
Option_Double option_from_double(double data);