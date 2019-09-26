#pragma once
#include <stdbool.h>

typedef struct Result
{
    bool valid;
    char *reason;
    void *data;
} Result;

Result result_get_empty_valid();
Result result_get_invalid_reason(char *reason);

typedef struct Option
{
    bool some;
    void *data;
} Option;
