#include "option.h"

Result result_get_empty_valid()
{
    Result empty_valid;
    empty_valid.valid = true;
    return empty_valid;
}

Result result_get_invalid_reason(char *reason)
{
    Result invalid;
    invalid.valid = false;
    invalid.reason = reason;
    return invalid;
}
