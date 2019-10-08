#include "python/qopmodule.h"

PyMODINIT_FUNC PyInit_qop(void) {
  PyObject *mod;

  if (PyType_Ready(&QopCircuitType) < 0) return NULL;

  if (PyType_Ready(&QopGateType) < 0) return NULL;

  mod = PyModule_Create(&qopmodule);
  if (mod == NULL) return NULL;

  import_array();
  if (PyErr_Occurred()) {
    return NULL;
  }

  QopError = PyErr_NewException("qop.error", NULL, NULL);
  Py_XINCREF(QopError);
  if (PyModule_AddObject(mod, "error", QopError) < 0) {
    Py_XDECREF(QopError);
    Py_CLEAR(QopError);
    Py_DECREF(mod);
    return NULL;
  }

  Py_INCREF(&QopCircuitType);
  if (PyModule_AddObject(mod, "Circuit", (PyObject *)&QopCircuitType) < 0) {
    Py_DECREF(&QopCircuitType);
    Py_XDECREF(QopError);
    Py_CLEAR(QopError);
    Py_DECREF(mod);
    return NULL;
  }

  Py_INCREF(&QopGateType);
  if (PyModule_AddObject(mod, "Gate", (PyObject *)&QopGateType) < 0) {
    Py_DECREF(&QopGateType);
    Py_DECREF(&QopCircuitType);
    Py_XDECREF(QopError);
    Py_CLEAR(QopError);
    Py_DECREF(mod);
    return NULL;
  }

  return mod;
}

static PyObject *qop_create_circuit(PyTypeObject *type, PyObject *args,
                                    PyObject *kwds) {
  unsigned int qubit_count;
  {
    int signed_qubit_count;
    char *kwarg_names[] = {"qubits", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i", kwarg_names,
                                     &signed_qubit_count))
      return NULL;
    if (signed_qubit_count < 0) {
      PyErr_SetString(PyExc_ValueError, "qubit count cannot be negative");
      return NULL;
    }
    qubit_count = (unsigned int)signed_qubit_count;
  }

  Circuit new_circuit;
  {
    Result init_r = circuit_init(&new_circuit, qubit_count);
    if (!init_r.valid) {
      PyErr_SetString(QopError, init_r.content.error_details.reason);
      return NULL;
    }
  }

  Vector gate_obj_refs;
  {
    Result init_r = vector_init(&gate_obj_refs, sizeof(QopGateObject *), 0);
    if (!init_r.valid) {
      PyErr_SetString(QopError, init_r.content.error_details.reason);
      return NULL;
    }
  }

  QopCircuitObject *self;
  self = (QopCircuitObject *)type->tp_alloc(type, 0);
  if (self == NULL) {
    circuit_free(&new_circuit);
    return NULL;
  }

  self->circuit = new_circuit;
  self->qubit_count = qubit_count;
  self->gate_obj_refs = gate_obj_refs;
  return (PyObject *)self;
}

static PyObject *qop_create_gate(PyTypeObject *type, PyObject *args,
                                 PyObject *kwds) {
  char *identifier_str;
  PyObject *matrix_optional = NULL;
  PyObject *params_optional = NULL;

  char *kwarg_names[] = {"identifier", "matrix", "parameters", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|$OO", kwarg_names,
                                   &identifier_str, &matrix_optional,
                                   &params_optional)) {
    PyErr_SetString(PyExc_ValueError,
                    "could not parse arguments; "
                    "expecting {identifier, matrix?, parameters?}");
    return NULL;
  }

  GateId id = GateNoId;
  {
    // strcmp returns the *difference* in strings (in char count)
    if (!strcmp(identifier_str, "c")) {
      id = GateCustom;
    } else if (!strcmp(identifier_str, "i")) {
      id = GateI;
    } else if (!strcmp(identifier_str, "x")) {
      id = GateX;
    } else if (!strcmp(identifier_str, "y")) {
      id = GateY;
    } else if (!strcmp(identifier_str, "z")) {
      id = GateZ;
    } else if (!strcmp(identifier_str, "h")) {
      id = GateH;
    } else if (!strcmp(identifier_str, "sqrtx")) {
      id = GateSqrtX;
    } else if (!strcmp(identifier_str, "t")) {
      id = GateT;
    } else if (!strcmp(identifier_str, "s")) {
      id = GateS;
    } else if (!strcmp(identifier_str, "rx")) {
      id = GateRx;
    } else if (!strcmp(identifier_str, "ry")) {
      id = GateRy;
    } else if (!strcmp(identifier_str, "rz")) {
      id = GateRz;
    } else {
      PyErr_SetString(PyExc_ValueError,
                      "unknown gate identifier; \
  use one of i,c,x,y,z,h,sqrtx,t,s,rx,ry,rz");
      return NULL;
    }
  }

  Gate gate;
  {
    Result init_result;

    switch (id) {
      case GateCustom: {
        PyArrayObject *matrix = NULL;

        // Parse matrix
        if (matrix_optional == NULL) {
          PyErr_SetString(PyExc_ValueError,
                          "specified custom gate, but did not specify matrix");
          return NULL;
        }

        matrix = (PyArrayObject *)PyArray_FROMANY(
            matrix_optional, NPY_CDOUBLE, 0, 0,
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);
        if (matrix == NULL) {
          return NULL;
        }

        if (PyArray_NDIM(matrix) != 2 || PyArray_DIMS(matrix)[0] != 2 ||
            PyArray_DIMS(matrix)[1] != 2) {
          PyErr_SetString(PyExc_ValueError,
                          "given matrix dimensions are not 2x2");
          return NULL;
        }

        double _Complex gate_matrix[2][2];
        memcpy(gate_matrix, PyArray_DATA(matrix), sizeof(double _Complex) * 4);

        init_result = gate_init_from_matrix(&gate, gate_matrix, NULL);
      } break;
      case GateRx:
      case GateRy:
      case GateRz: {
        if (params_optional == NULL) {
          PyErr_SetString(PyExc_ValueError,
                          "specified rotation gate but no parameter");
          return NULL;
        }

        PyArrayObject *params_array = (PyArrayObject *)PyArray_FROMANY(
            params_optional, NPY_DOUBLE, 0, 0,
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);
        if (PyArray_NDIM(params_array) != 1 ||
            PyArray_DIMS(params_array)[0] != 1) {
          PyErr_SetString(PyExc_ValueError,
                          "specified rotation gate but wrong number of "
                          "parameters (expecting 1)");
          return NULL;
        }

        double params[1];
        memcpy(params, PyArray_DATA(params_array), sizeof(double) * 1);

        init_result = gate_init_from_identifier(&gate, id, params);
      } break;
      default: {
        init_result = gate_init_from_identifier(&gate, id, NULL);
      } break;
    }

    if (!init_result.valid) {
      PyErr_SetString(QopError, init_result.content.error_details.reason);
      return NULL;
    }
  }

  QopGateObject *self;
  self = (QopGateObject *)type->tp_alloc(type, 0);
  if (self == NULL) {
    gate_free(&gate);
    return NULL;
  }

  self->gate = gate;
  return (PyObject *)self;
}

static PyObject *qop_circuit_add_gate(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwrds) {
  QopGateObject *gate_obj;
  int qubit = -1;
  int control = -1;

  char *kwarg_names[] = {"gate", "qubit", "control", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwrds, "O!i|i", kwarg_names,
                                   &QopGateType, &gate_obj, &qubit, &control)) {
    PyErr_SetString(PyExc_ValueError,
                    "could not parse arguments; "
                    "expecting {gate, qubit, control?}");
    return NULL;
  }

  if (qubit < 0 || (unsigned int)qubit >= self->qubit_count) {
    PyErr_SetString(PyExc_ValueError, "invalid qubit");
    return NULL;
  }

  if (control == qubit) {
    PyErr_SetString(PyExc_ValueError, "invalid control");
    return NULL;
  }

  Option_Uint control_opt;
  if (control < 0)
    control_opt = option_none_uint();
  else
    control_opt = option_from_uint((unsigned int)control);

  {
    // No not Py_INCREF gate! It seems like it's already done.
    Result add_r = circuit_add_gate(&self->circuit, &gate_obj->gate,
                                    (unsigned int)qubit, control_opt);
    if (!add_r.valid) {
      PyErr_SetString(QopError, add_r.content.error_details.reason);
      Py_DECREF(gate_obj);
      return NULL;
    }
  }

  {
    Result push_r =
        vector_push(&self->gate_obj_refs, (QopGateObject **)(&gate_obj));
    if (!push_r.valid) {
      PyErr_SetString(QopError, push_r.content.error_details.reason);
      Py_DECREF(gate_obj);
      return NULL;
    }
  }

  return (PyObject *)self;
}

static PyObject *qop_optimize_circuit(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds) {
  PyDictObject *settings = NULL;
  char *kwarg_names[] = {"settings", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!", kwarg_names, &PyDict_Type,
                                   &settings)) {
    PyErr_SetString(PyExc_ValueError,
                    "could not parse arguments to circuit optimization");
    return NULL;
  }

  AdadeltaSettings ada_settings = optimizer_adadelta_get_default();

  if (settings != NULL) {
    PyDictObject *ada_settings =
        (PyDictObject *)PyDict_GetItemString(settings, "ada");
    if (ada_settings != NULL) {

    }
  }
}

static void qop_circuit_obj_dealloc(QopCircuitObject *self) {
  printf("dealoc circuit\n");
  circuit_free(&self->circuit);
  {
    Iter refs_iter = vector_iter_create(&self->gate_obj_refs);
    Option next;
    while ((next = iter_next(&refs_iter)).some) {
      QopGateObject *gate_obj_ptr = *(QopGateObject **)next.data;
      Py_DECREF(gate_obj_ptr);
    }
  }
  vector_free(&self->gate_obj_refs);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static void qop_gate_obj_dealloc(QopGateObject *self) {
  printf("dealoc gate\n");
  gate_free(&self->gate);
  Py_TYPE(self)->tp_free((PyObject *)self);
}
