#include "include/gate.h"

// Initialize a new gate from a matrix array.
// It is important that we `memcpy` the array, rather than assigning it,
// because we want to make sure that reparameterization doesn't affect
// other gates.
Result gate_init_from_matrix(Gate *gate, double _Complex matrix[2][2],
                             ReparamFnPtr reparamFn) {
  // result.matrix = matrix;
  memcpy(gate->matrix, matrix, GATE_SINGLE_QUBIT_SIZE);
  gate->reparamFn = reparamFn;
  gate->id = GateCustom;

  return result_get_valid_with_data(gate);
}

// Initializes a gate using an identifier from the GateId enum
// It's basically `gate_new_from_matrix` with a big switch statement.
Result gate_init_from_identifier(Gate *gate, GateId identifier,
                                 double params[]) {
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
        return result_get_invalid_reason(
            "required theta parameter was not specified for GateRx "
            "instantiation");
      }
      double theta = params[0];
      double _Complex theta_matrix[2][2] = {
          {cos(theta / 2), -sin(theta / 2) * _Complex_I},
          {-sin(theta / 2) * _Complex_I, cos(theta / 2)}};
      memcpy(matrix, theta_matrix, GATE_SINGLE_QUBIT_SIZE);
      reparam_fn = &gate_reparameterize_rx;
      break;
    }
    case GateRy: {
      if (params == NULL) {
        return result_get_invalid_reason(
            "required theta parameter was not specified for GateRy "
            "instantiation");
      }
      double theta = params[0];
      double _Complex theta_matrix[2][2] = {{cos(theta / 2), -sin(theta / 2)},
                                            {sin(theta / 2), cos(theta / 2)}};
      memcpy(matrix, theta_matrix, GATE_SINGLE_QUBIT_SIZE);
      reparam_fn = &gate_reparameterize_ry;
      break;
    }
    case GateRz: {
      if (params == NULL) {
        return result_get_invalid_reason(
            "required theta parameter was not specified for GateRz "
            "instantiation");
      }
      double theta = params[0];
      double _Complex theta_matrix[2][2] = {{1.0, 0.0},
                                            {0.0, cexp(_Complex_I * theta)}};
      memcpy(matrix, theta_matrix, GATE_SINGLE_QUBIT_SIZE);
      reparam_fn = &gate_reparameterize_rz;
      break;
    }
    default:
      return result_get_invalid_reason("unknown GateId specification");
      break;
  }

  // Pack everything into the given struct
  {
    gate->id = identifier;
    memcpy(gate->matrix, matrix, GATE_SINGLE_QUBIT_SIZE);
    gate->reparamFn = reparam_fn;
  }

  // Return pointer as result
  return result_get_valid_with_data(gate);
}

Result gate_clone(Gate *from, Gate *into)
{
  // There is currently no internal memory allocations in the implementation
  // of the Gate struct
  memcpy(into, from, sizeof(Gate));
  return result_get_empty_valid();
}

void gate_reparameterize_rx(double _Complex (*matrix)[2][2], double params[]) {
  double theta = params[0];
  (*matrix)[0][0] = (double _Complex)cos(theta / 2);
  (*matrix)[0][1] = -sin(theta / 2) * _Complex_I;
  (*matrix)[1][0] = -sin(theta / 2) * _Complex_I;
  (*matrix)[1][1] = (double _Complex)cos(theta / 2);
}

void gate_reparameterize_ry(double _Complex (*matrix)[2][2], double params[]) {
  double theta = params[0];
  (*matrix)[0][0] = (double _Complex)cos(theta / 2);
  (*matrix)[0][1] = (double _Complex)-sin(theta / 2);
  (*matrix)[1][0] = (double _Complex)sin(theta / 2);
  (*matrix)[1][1] = (double _Complex)cos(theta / 2);
}

void gate_reparameterize_rz(double _Complex (*matrix)[2][2], double params[]) {
  double theta = params[0];
  (*matrix)[1][1] = cexp(theta * _Complex_I);
}

void gate_free(Gate *gate) {
  /* No internal heap allocation in current implementation */
}