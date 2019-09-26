#include "util.h"

/**
 * Approximate equality for floating point types.
 *
 * @see https://floating-point-gui.de/errors/comparison/
 */
bool double_approx_eq(double lhs, double rhs)
{
    double abs_lhs = fabs(lhs);
    double abs_rhs = fabs(rhs);
    double abs_diff = fabs(lhs - rhs);

    if (lhs == rhs)
    {
        return true;
    }
    else if (lhs == 0 || rhs == 0 || (abs_lhs + abs_rhs < DBL_EPSILON))
    {
        return abs_diff < DBL_EPSILON;
    }
    else
    {
        return abs_diff / fmin(abs_lhs + abs_rhs, DBL_MAX) < DBL_EPSILON;
    }
}

bool float_approx_eq(float lhs, float rhs)
{
    float abs_lhs = fabsf(lhs);
    float abs_rhs = fabsf(rhs);
    float abs_diff = fabsf(lhs - rhs);

    if (lhs == rhs)
    {
        return true;
    }
    else if (lhs == 0 || rhs == 0 || (abs_lhs + abs_rhs < FLT_EPSILON))
    {
        return abs_diff < FLT_EPSILON;
    }
    else
    {
        return abs_diff / fminf(abs_lhs + abs_rhs, FLT_MAX) < FLT_EPSILON;
    }
}

void *get_data_from_result(struct Result result)
{
    if (result.valid) {
        return result.data;
    } else {
        return NULL;
    }
}