#pragma once
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdlib.h>

bool double_approx_eq(double lhs, double rhs);
bool float_approx_eq(float lhs, float rhs);

struct Result
{
    bool valid;
    char *reason;
    void *data;
};

void *get_data_from_result(struct Result result);