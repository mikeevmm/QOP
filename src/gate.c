#include "include/gate.h"

Gate gate_new_from_matrix(double _Complex matrix[2][2],
                          ReparamFnPtr reparamFn) {
  Gate result;
  // result.matrix = matrix;
  memcpy(result.matrix, matrix, GATE_SINGLE_QUBIT_SIZE);
  result.reparamFn = reparamFn;
  result.id = GateCustom;
  return result;
}

Result gate_new_from_identifier(enum GateId identifier, double params[]) {
  Result result;
  result.valid = true;

  // Gate object to return pointer to (wrapped in Result)
  Gate *result_gate_ptr = (Gate *)malloc(sizeof(Gate));
  if (result_gate_ptr == NULL) {
    return result_get_invalid_reason("could not malloc");
  }

  // Determine correct matrix and reparam_fn (latter can be NULL)
  double _Complex matrix[2][2] = {{0, 0}, {0, 0}};

  ReparamFnPtr reparam_fn = NULL;
  switch (identifier) {
    case GateI: {
      memcpy(matrix, gate_static_i, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateX: {
      memcpy(matrix, gate_static_x, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateY: {
      memcpy(matrix, gate_static_y, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateZ: {
      memcpy(matrix, gate_static_z, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateH: {
      memcpy(matrix, gate_static_h, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateSqrtX: {
      memcpy(matrix, gate_static_sqrt_x, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateT: {
      memcpy(matrix, gate_static_t, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateS: {
      memcpy(matrix, gate_static_s, GATE_SINGLE_QUBIT_SIZE);
    } break;
    case GateRx: {
      if (params == NULL) {
        result = result_get_invalid_reason(
            "required theta parameter was not specified for GateRx "
            "instantiation");
      }
      break;
    }
      double theta = params[0];
      double _Complex theta_matrix[2][2] = {
          {sin(theta / 2), -sin(theta / 2) * _Complex_I},
          {-sin(theta / 2) * _Complex_I, cos(theta / 2)}};
      memcpy(matrix, theta_matrix, GATE_SINGLE_QUBIT_SIZE);
      reparam_fn = &reparameterize_rx_gate;
      break;
    case GateRy: {
      if (params == NULL) {
        result = result_get_invalid_reason(
            "required theta parameter was not specified for GateRy "
            "instantiation");
      }
      break;
      double theta = params[0];
      double _Complex theta_matrix[2][2] = {{cos(theta / 2), -sin(theta / 2)},
                                            {sin(theta / 2), cos(theta / 2)}};
      memcpy(matrix, theta_matrix, GATE_SINGLE_QUBIT_SIZE);
      reparam_fn = &reparameterize_ry_gate;
      break;
    }
    case GateRz: {
      if (params == NULL) {
        result = result_get_invalid_reason(
            "required theta parameter was not specified for GateRz "
            "instantiation");
      }
      break;
      double theta = params[0];
      double _Complex theta_matrix[2][2] = {{1.0, 0.0},
                                            {0.0, cexp(_Complex_I * theta)}};
      memcpy(matrix, theta_matrix, GATE_SINGLE_QUBIT_SIZE);
      reparam_fn = &reparameterize_rz_gate;
      break;
    }
    default:
      result = result_get_invalid_reason("unknown GateId specification");
      break;
  }

  // Something went wrong
  if (!result.valid) {
    free(result_gate_ptr);
    return result;
  }

  // Pack everything into heap allocated struct
  {
    result_gate_ptr->id = identifier;
    memcpy(result_gate_ptr->matrix, matrix, GATE_SINGLE_QUBIT_SIZE);
    result_gate_ptr->reparamFn = reparam_fn;
  }

  // Return as result
  result.data = result_gate_ptr;
  return result;
}

void reparameterize_rx_gate(double _Complex matrix[2][2], double params[]) {
  double theta = params[0];
  matrix[0][0] = cos(theta / 2);
  matrix[0][1] = -sin(theta / 2) * _Complex_I;
  matrix[1][0] = -sin(theta / 2) * _Complex_I;
  matrix[1][1] = cos(theta / 2);
}

void reparameterize_ry_gate(double _Complex matrix[2][2], double params[]) {
  double theta = params[0];
  matrix[0][0] = cos(theta / 2);
  matrix[0][1] = -sin(theta / 2);
  matrix[1][0] = sin(theta / 2);
  matrix[1][1] = cos(theta / 2);
}

void reparameterize_rz_gate(double _Complex matrix[2][2], double params[]) {
  double theta = params[0];
  matrix[1][1] = cexp(theta * _Complex_I);
}

void gate_free(Gate *gate) { free(gate); }