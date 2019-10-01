#include <util.h>

/**
 * @see https://floating-point-gui.de/errors/comparison/
 */
bool double_approx_eq(double lhs, double rhs) {
  double abs_lhs = fabs(lhs);
  double abs_rhs = fabs(rhs);
  double abdiff = fabs(lhs - rhs);

  if (lhs == rhs) {
    return true;
  } else if (lhs == 0 || rhs == 0 || (abs_lhs + abs_rhs < DBL_EPSILON)) {
    return abdiff < DBL_EPSILON;
  } else {
    return abdiff / fmin(abs_lhs + abs_rhs, DBL_MAX) < DBL_EPSILON;
  }
}

bool float_approx_eq(float lhs, float rhs) {
  float abs_lhs = fabsf(lhs);
  float abs_rhs = fabsf(rhs);
  float abdiff = fabsf(lhs - rhs);

  if (lhs == rhs) {
    return true;
  } else if (lhs == 0 || rhs == 0 || (abs_lhs + abs_rhs < FLT_EPSILON)) {
    return abdiff < FLT_EPSILON;
  } else {
    return abdiff / fminf(abs_lhs + abs_rhs, FLT_MAX) < FLT_EPSILON;
  }
}