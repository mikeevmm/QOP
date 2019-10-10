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
                                     &signed_qubit_count)) {
      return NULL;
    }
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
  self->compacted_and_hardened = false;
  return (PyObject *)self;
}

static void qop_circuit_obj_dealloc(QopCircuitObject *self) {
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
  gate_free(&self->gate);
  free(self->params);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *qop_gate_create(PyTypeObject *type, PyObject *args,
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
                    "expecting (identifier, matrix?, parameters?)");
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
  double *params = NULL;
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

        {
          void *mem = malloc(sizeof(double) * 1);
          if (mem == NULL) {
            PyErr_NoMemory();
            return NULL;
          }
          params = (double *)mem;
        }
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
  self->params = params;
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

  Py_INCREF(gate_obj);

  {
    Result add_r = circuit_add_gate(&self->circuit, &gate_obj->gate,
                                    (unsigned int)qubit, control_opt);
    if (!add_r.valid) {
      PyErr_SetString(QopError, add_r.content.error_details.reason);
      Py_DECREF(gate_obj);
      return NULL;
    }
  }

  self->compacted_and_hardened = false;

  {
    Result push_r =
        vector_push(&self->gate_obj_refs, (QopGateObject **)(&gate_obj));
    if (!push_r.valid) {
      PyErr_SetString(QopError, push_r.content.error_details.reason);
      Py_DECREF(gate_obj);
      return NULL;
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *qop_circuit_optimize(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds) {
  PyDictObject *settings = NULL;
  PyObject *hamiltonian_arg;
  char *kwarg_names[] = {"hamiltonian", "settings", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O!", kwarg_names,
                                   &hamiltonian_arg, &PyDict_Type, &settings)) {
    PyErr_SetString(PyExc_ValueError,
                    "could not parse arguments to circuit optimization; "
                    "expected {hamiltonian, settings?}");
    return NULL;
  }

  PyArrayObject *hamiltonian_arr;
  {
    PyObject *from_obj =
        PyArray_FROMANY(hamiltonian_arg, NPY_CDOUBLE, 2, 2,
                        NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);
    if (from_obj == NULL) {
      PyErr_SetString(PyExc_ValueError,
                      "could not interpret hamiltonian as matrix");
      return NULL;
    }
    hamiltonian_arr = (PyArrayObject *)from_obj;
  }

  {
    int dims = PyArray_NDIM(hamiltonian_arr);
    npy_intp *size = PyArray_DIMS(hamiltonian_arr);
    if (dims != 2 || size[0] != size[1]) {
      PyErr_SetString(PyExc_ValueError,
                      "hamiltonian is malformed; must be square 2D array");
      Py_DECREF(hamiltonian_arr);
      return NULL;
    }

    if (size[0] != (1 << self->circuit.depth[0])) {
      PyErr_SetString(
          PyExc_ValueError,
          "hamiltonian is malformed; does not match circuit state size");
      Py_DECREF(hamiltonian_arr);
      return NULL;
    }
  }

  // Prepare the circuit
  if (!self->compacted_and_hardened) {
    {
      Result comp_r = circuit_compact(&self->circuit);
      if (!comp_r.valid) {
        PyErr_SetString(QopError, comp_r.content.error_details.reason);
        Py_DECREF(hamiltonian_arr);
        return NULL;
      }
    }
    {
      Result hard_r = circuit_harden(&self->circuit);
      if (!hard_r.valid) {
        PyErr_SetString(QopError, hard_r.content.error_details.reason);
        Py_DECREF(hamiltonian_arr);
        return NULL;
      }
    }
    self->compacted_and_hardened = true;
  }

  // Default settings; to be possibly changed using user-supplied
  // dictionary
  AdadeltaSettings ada_settings = optimizer_adadelta_get_default();
  double stop_at = 1e-8;
  Vector reparams_vec;
  vector_init(&reparams_vec, sizeof(GateParameterization), 0);
  Vector reparams_to_gate_obj_ptr;
  vector_init(&reparams_to_gate_obj_ptr, sizeof(QopGateObject *), 0);
  bool any_reparam_given = false;
  int max_iters = -1;

  /*
    {
      'ada': {
        'rho': val,
        'epsilon': val
      },
      'optimize':
      {
        'gates': [ val, val, ... ],
        'deltas': [ [val], [val, val], ... ],
        'stop_at': val,
        'max_iterations': val
      }
    }
  */
  // Parsing of options dictionary
  if (settings != NULL) {
    if (!PyDict_Check(settings)) {
      PyErr_SetString(PyExc_ValueError, "settings object must be a dictionary");
      vector_free(&reparams_vec);
      vector_free(&reparams_to_gate_obj_ptr);
      Py_DECREF(hamiltonian_arr);
      return NULL;
    }

    PyObject *ada_dict = PyDict_GetItemString((PyObject *)settings, "ada");
    if (ada_dict != NULL) {
      if (!PyDict_Check(ada_dict)) {
        PyErr_SetString(PyExc_ValueError,
                        "settings:ada object must be a dictionary");
        vector_free(&reparams_vec);
        vector_free(&reparams_to_gate_obj_ptr);
        Py_DECREF(hamiltonian_arr);
        return NULL;
      }

      PyObject *rho = PyDict_GetItemString(ada_dict, "rho");
      if (rho != NULL) {
        if (!PyFloat_Check(rho)) {
          PyErr_SetString(PyExc_ValueError,
                          "settings:ada:rho object must be a number");
          vector_free(&reparams_vec);
          vector_free(&reparams_to_gate_obj_ptr);
          Py_DECREF(hamiltonian_arr);
          return NULL;
        }
        ada_settings.rho = PyFloat_AsDouble(rho);
      }

      PyObject *epsilon = PyDict_GetItemString(ada_dict, "epsilon");
      if (epsilon != NULL) {
        if (!PyNumber_Check(epsilon)) {
          PyErr_SetString(PyExc_ValueError,
                          "settings:ada:epsilon object must be a number");
          vector_free(&reparams_vec);
          vector_free(&reparams_to_gate_obj_ptr);
          Py_DECREF(hamiltonian_arr);
          return NULL;
        }
        ada_settings.epsilon = PyFloat_AsDouble(epsilon);
      }
    }

    PyObject *opt_dict = PyDict_GetItemString((PyObject *)settings, "optimize");
    if (opt_dict != NULL) {
      if (!PyDict_Check(opt_dict)) {
        PyErr_SetString(PyExc_ValueError,
                        "settings:optimize object must be a dictionary");
        vector_free(&reparams_vec);
        vector_free(&reparams_to_gate_obj_ptr);
        Py_DECREF(hamiltonian_arr);
        return NULL;
      }

      PyObject *gates = PyDict_GetItemString(opt_dict, "gates");
      PyObject *deltas = PyDict_GetItemString(opt_dict, "deltas");

      if ((gates != NULL) ^ (deltas != NULL)) {
        PyErr_SetString(PyExc_ValueError,
                        "settings:optimize:deltas and settings:optimize:gates "
                        "must both be defined or undefined");
        vector_free(&reparams_vec);
        vector_free(&reparams_to_gate_obj_ptr);
        Py_DECREF(hamiltonian_arr);
        return NULL;
      }

      if (gates != NULL) {
        if (!(PyList_Check(gates) || PyTuple_Check(gates) ||
              PyIter_Check(gates))) {
          PyErr_SetString(PyExc_ValueError,
                          "settings:optimize:gates object must be iterable");
          vector_free(&reparams_vec);
          vector_free(&reparams_to_gate_obj_ptr);
          Py_DECREF(hamiltonian_arr);
          return NULL;
        }

        if (!(PyList_Check(deltas) || PyTuple_Check(deltas) ||
              PyIter_Check(deltas))) {
          PyErr_SetString(PyExc_ValueError,
                          "settings:optimize:deltas object must be iterable");
          vector_free(&reparams_vec);
          vector_free(&reparams_to_gate_obj_ptr);
          Py_DECREF(hamiltonian_arr);
        }

        PyObject *gates_iter = PyObject_GetIter(gates);
        PyObject *deltas_iter = PyObject_GetIter(deltas);

        if (gates_iter == NULL) {
          vector_free(&reparams_vec);
          vector_free(&reparams_to_gate_obj_ptr);
          Py_DECREF(hamiltonian_arr);
          return NULL;
        }

        if (deltas_iter == NULL) {
          vector_free(&reparams_vec);
          vector_free(&reparams_to_gate_obj_ptr);
          Py_DECREF(hamiltonian_arr);
          Py_DECREF(gates_iter);
          return NULL;
        }

        // Do not reparameterize any gate, even if the following iterators
        // are empty
        any_reparam_given = true;

        PyObject *next_gate;
        PyObject *next_deltas;
        while ((next_gate = PyIter_Next(gates_iter)) != NULL) {
          next_deltas = PyIter_Next(deltas_iter);
          if (next_deltas == NULL) {
            PyErr_SetString(PyExc_ValueError,
                            "size mismatch between settings:optimize:gates and "
                            "settings:optimize:deltas");
            vector_free(&reparams_vec);
            vector_free(&reparams_to_gate_obj_ptr);
            Py_DECREF(hamiltonian_arr);
            Py_XDECREF(next_gate);
            Py_XDECREF(next_deltas);
            Py_DECREF(gates_iter);
            Py_DECREF(deltas_iter);
            return NULL;
          }

          if (Py_TYPE(next_gate) != &QopGateType) {
            PyErr_SetString(PyExc_ValueError,
                            "settings:optimize:gates:element must be qop.Gate");
            vector_free(&reparams_vec);
            vector_free(&reparams_to_gate_obj_ptr);
            Py_DECREF(hamiltonian_arr);
            Py_XDECREF(next_gate);
            Py_XDECREF(next_deltas);
            Py_DECREF(gates_iter);
            Py_DECREF(deltas_iter);
            return NULL;
          }

          if (!(PyList_Check(next_deltas) || PyTuple_Check(next_deltas) ||
                PyIter_Check(next_deltas))) {
            PyErr_SetString(
                PyExc_ValueError,
                "settings:optimize:deltas:element must be iterable ");
            vector_free(&reparams_vec);
            vector_free(&reparams_to_gate_obj_ptr);
            Py_DECREF(hamiltonian_arr);
            Py_XDECREF(next_gate);
            Py_XDECREF(next_deltas);
            Py_DECREF(gates_iter);
            Py_DECREF(deltas_iter);
            return NULL;
          }

          PyObject *next_deltas_iter = PyObject_GetIter(next_deltas);
          if (next_deltas_iter == NULL) {
            vector_free(&reparams_vec);
            vector_free(&reparams_to_gate_obj_ptr);
            Py_DECREF(hamiltonian_arr);
            Py_XDECREF(next_gate);
            Py_XDECREF(next_deltas);
            Py_DECREF(gates_iter);
            Py_DECREF(deltas_iter);
            return NULL;
          }

          Vector deltas;
          vector_init(&deltas, sizeof(double), 1);
          {
            PyObject *next_delta;
            while ((next_delta = PyIter_Next(next_deltas_iter)) != NULL) {
              double delta = PyFloat_AsDouble(next_delta);
              Py_DECREF(next_delta);

              Result push_result = vector_push(&deltas, &delta);
              if (!push_result.valid) {
                PyErr_SetString(QopError,
                                push_result.content.error_details.reason);
                {
                  Iter reparam_iter = vector_iter_create(&reparams_vec);
                  Option next;
                  while ((next = iter_next(&reparam_iter)).some) {
                    optimizer_gate_param_free(
                        (GateParameterization *)next.data);
                  }
                  vector_free(&reparams_vec);
                  vector_free(&reparams_to_gate_obj_ptr);
                }
                vector_free(&deltas);
                Py_DECREF(hamiltonian_arr);
                Py_DECREF(next_delta);
                Py_DECREF(next_deltas_iter);
                Py_XDECREF(next_gate);
                Py_XDECREF(next_deltas);
                Py_DECREF(gates_iter);
                Py_DECREF(deltas_iter);
                return NULL;
              }
            }

            Py_DECREF(next_deltas_iter);
          }

          QopGateObject *gate = (QopGateObject *)next_gate;
          GateParameterization param;
          optimizer_gate_param_init(&param, &gate->gate, deltas.size,
                                    gate->params, deltas.data);
          vector_push(&reparams_vec, &param);
          vector_push(&reparams_to_gate_obj_ptr, &gate);
          vector_free(&deltas);

          Py_XDECREF(next_gate);
          Py_XDECREF(next_deltas);
        }
        Py_DECREF(gates_iter);
        Py_DECREF(deltas_iter);
      }

      PyObject *stop_at_obj = PyDict_GetItemString(opt_dict, "stop_at");
      if (stop_at_obj != NULL) {
        if (!PyNumber_Check(stop_at_obj)) {
          PyErr_SetString(PyExc_ValueError,
                          "settings:optimize:stop_at object must be a number");
          {
            Iter reparam_iter = vector_iter_create(&reparams_vec);
            Option next;
            while ((next = iter_next(&reparam_iter)).some) {
              optimizer_gate_param_free((GateParameterization *)next.data);
            }
            vector_free(&reparams_vec);
            vector_free(&reparams_to_gate_obj_ptr);
          }
          Py_DECREF(hamiltonian_arr);
          return NULL;
        }

        stop_at = PyFloat_AsDouble(stop_at_obj);
      }

      PyObject *max_iters_obj =
          PyDict_GetItemString(opt_dict, "max_iterations");
      if (max_iters_obj != NULL) {
        if (!PyNumber_Check(max_iters_obj)) {
          PyErr_SetString(
              PyExc_ValueError,
              "settings:optimize:max_iterations object must be a number");
          vector_free(&reparams_vec);
          vector_free(&reparams_to_gate_obj_ptr);
          Py_DECREF(hamiltonian_arr);
          return NULL;
        }

        max_iters = (int)PyFloat_AsDouble(max_iters_obj);
      }
    }
  }

  if (!any_reparam_given) {
    // Generate a reparam for all known parameterized gates in the
    // circuit
    Iter known_gates_iter = vector_iter_create(&self->gate_obj_refs);
    Option next;
    while ((next = iter_next(&known_gates_iter)).some) {
      QopGateObject *qop_gate = *(QopGateObject **)next.data;
      switch (qop_gate->gate.id) {
        case GateRx:
        case GateRy:
        case GateRz: {
          GateParameterization param;
          double delta[] = {0.1};
          {
            Result init_r = optimizer_gate_param_init(
                &param, &qop_gate->gate, 1, qop_gate->params, delta);
            if (!init_r.valid) {
              PyErr_SetString(QopError, init_r.content.error_details.reason);
              {
                Iter reparam_iter = vector_iter_create(&reparams_vec);
                Option next;
                while ((next = iter_next(&reparam_iter)).some) {
                  optimizer_gate_param_free((GateParameterization *)next.data);
                }
                vector_free(&reparams_vec);
                vector_free(&reparams_to_gate_obj_ptr);
              }
              vector_free(&reparams_vec);
              vector_free(&reparams_to_gate_obj_ptr);
              Py_DECREF(hamiltonian_arr);
            }
          }
          vector_push(&reparams_vec, &param);
          vector_push(&reparams_to_gate_obj_ptr, &qop_gate);
        } break;
        default:
          continue;
      }
    }
  }

  if (reparams_vec.size == 0) {
    PyErr_SetString(QopError, "nothing to optimize");
    vector_free(&reparams_vec);
    vector_free(&reparams_to_gate_obj_ptr);
    Py_DECREF(hamiltonian_arr);
    return NULL;
  }

  OptimizerSettings opt_settings;
  {
    Result init_result = optimizer_settings_init(
        &opt_settings, &self->circuit, PyArray_DATA(hamiltonian_arr), stop_at,
        reparams_vec.data, reparams_vec.size, max_iters);
    if (!init_result.valid) {
      PyErr_SetString(QopError, init_result.content.error_details.reason);
      {
        Iter reparam_iter = vector_iter_create(&reparams_vec);
        Option next;
        while ((next = iter_next(&reparam_iter)).some) {
          optimizer_gate_param_free((GateParameterization *)next.data);
        }
        vector_free(&reparams_vec);
        vector_free(&reparams_to_gate_obj_ptr);
      }
      Py_DECREF(hamiltonian_arr);
    }
  }

  // Perform optimization
  Optimizer optimizer;
  optimizer_init(&optimizer, opt_settings, ada_settings);
  OptimizationResult opt_result = optimizer_optimize(&optimizer);
  if (!opt_result.valid) {
    PyErr_SetString(QopError, opt_result.content.error_details.reason);
    {
      Iter reparam_iter = vector_iter_create(&reparams_vec);
      Option next;
      while ((next = iter_next(&reparam_iter)).some) {
        optimizer_gate_param_free((GateParameterization *)next.data);
      }
      vector_free(&reparams_vec);
      vector_free(&reparams_to_gate_obj_ptr);
    }
    Py_DECREF(hamiltonian_arr);
    optimizer_settings_free(&opt_settings);
  }

  // Free python objects
  Py_DECREF(hamiltonian_arr);

  // Free optimizer settings
  optimizer_settings_free(&opt_settings);

  // Free gate parameterizations, but save the results
  PyListObject *result_list = (PyListObject *)PyList_New(0);
  {
    Iter reparam_iter = vector_iter_create(&reparams_vec);
    Iter qop_gate_iter = vector_iter_create(&reparams_to_gate_obj_ptr);
    Option next_reparam;
    while ((next_reparam = iter_next(&reparam_iter)).some) {
      GateParameterization *param = (GateParameterization *)next_reparam.data;
      QopGateObject *qop_gate =
          *(QopGateObject **)(iter_next(&qop_gate_iter).data);

      PyObject *params_sublist = PyList_New(0);
      for (unsigned int i = 0; i < param->param_count; ++i) {
        PyFloatObject *param_value =
            (PyFloatObject *)PyFloat_FromDouble(*(param->params + i));
        PyList_Append((PyObject *)params_sublist, (PyObject *)param_value);

        memcpy(qop_gate->params, param->params,
               param->param_count * sizeof(double));
      }
      PyList_Append((PyObject *)result_list, (PyObject *)params_sublist);

      optimizer_gate_param_free(param);
    }
    vector_free(&reparams_vec);
    vector_free(&reparams_to_gate_obj_ptr);
  }

  return PyTuple_Pack(2, result_list,
                      PyBool_FromLong(opt_result.quit_on_max_iter));
}

static PyObject *qop_circuit_get_gates(QopCircuitObject *self) {
  PyListObject *result_list = (PyListObject *)PyList_New(0);
  if (result_list == NULL) {
    return NULL;
  }

  Iter ref_gates_iter = vector_iter_create(&self->gate_obj_refs);
  Option next;
  while ((next = iter_next(&ref_gates_iter)).some) {
    QopGateObject *gate_obj = *(QopGateObject **)next.data;
    Py_INCREF(gate_obj);
    int append_result =
        PyList_Append((PyObject *)result_list, (PyObject *)gate_obj);
    if (append_result < 0) {
      for (Py_ssize_t i = 0; i < PyList_Size((PyObject *)result_list); ++i) {
        Py_DECREF(PyList_GetItem((PyObject *)result_list, i));
      }
      Py_DECREF(result_list);
      return NULL;
    }
  }

  return (PyObject *)result_list;
}

static PyObject *qop_circuit_run(QopCircuitObject *self, PyObject *args,
                                 PyObject *kwds) {
  PyObject *state_in_obj;

  char *kwarg_names[] = {"state_in", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwarg_names,
                                   &state_in_obj)) {
    PyErr_SetString(PyExc_ValueError, "expecting state_in argument");
    return NULL;
  }

  double _Complex state_in[1U << self->qubit_count];

  {
    PyObject *as_arr =
        PyArray_FROMANY(state_in_obj, NPY_CDOUBLE, 1, 1,
                        NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);
    if (as_arr == NULL) {
      PyErr_SetString(PyExc_ValueError,
                      "cannot interpret state_in as 1D vector");
      return NULL;
    }
    unsigned int vec_size = PyArray_SIZE((PyObject *)as_arr);
    if (vec_size != (1U << self->qubit_count)) {
      PyErr_SetString(PyExc_ValueError,
                      "given input vector does not match size of circuit");
      return NULL;
    }
    memcpy(state_in, PyArray_DATA(as_arr), sizeof(double _Complex) * vec_size);
  }

  if (!self->compacted_and_hardened) {
    {
      Result comp_r = circuit_compact(&self->circuit);
      if (!comp_r.valid) {
        PyErr_SetString(QopError, comp_r.content.error_details.reason);
        return NULL;
      }
    }
    {
      Result hard_r = circuit_harden(&self->circuit);
      if (!hard_r.valid) {
        PyErr_SetString(QopError, hard_r.content.error_details.reason);
        return NULL;
      }
    }
    self->compacted_and_hardened = true;
  }

  {
    Result run_r = circuit_run(&self->circuit, &state_in);
    if (!run_r.valid) {
      PyErr_SetString(QopError, run_r.content.error_details.reason);
      return NULL;
    }
  }

  PyObject *return_list = PyList_New(0);
  if (return_list == NULL) {
    return NULL;
  }

  for (unsigned int i = 0; i < (1U << self->qubit_count); ++i) {
    double _Complex elem = state_in[i];
    PyList_Append(return_list, PyComplex_FromDoubles(creal(elem), cimag(elem)));
  }

  return return_list;
}

static PyObject *qop_gate_reparameterize(QopGateObject *self, PyObject *args,
                                         PyObject *kwds) {
  if (self->gate.reparamFn == NULL) {
    PyErr_SetString(QopError, "gate is not parameterized");
    return NULL;
  }

  PyObject *args_list;

  char *kwarg_names[] = {"parameters", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwarg_names, &args_list)) {
    PyErr_SetString(PyExc_ValueError, "expecting parameters iterable");
    return NULL;
  }

  if (!(PyList_Check(args_list) || PyTuple_Check(args_list) ||
        PyIter_Check(args_list))) {
    PyErr_SetString(PyExc_ValueError, "expecting iterable parameters");
    return NULL;
  }

  Vector params;
  vector_init(&params, sizeof(double), 1);
  {
    PyObject *args_iter = PyObject_GetIter(args_list);
    if (args_iter == NULL) {
      return NULL;
    }

    PyObject *next;
    while ((next = PyIter_Next(args_iter)) != NULL) {
      if (!PyNumber_Check(next)) {
        PyErr_SetString(PyExc_ValueError, "got unexpected non-numeric value");
        Py_DECREF(args_iter);
        return NULL;
      }
      double param = PyFloat_AsDouble(next);
      vector_push(&params, &param);

      Py_DECREF(next);
    }

    Py_DECREF(args_iter);
  }
  self->gate.reparamFn(&self->gate.matrix, params.data);
  memcpy(self->params, params.data, params.size * sizeof(double));
  vector_free(&params);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *qop_gate_get_parameters(QopGateObject *self) {
  PyListObject *result = (PyListObject *)PyList_New(0);
  unsigned int param_count;
  switch (self->gate.id) {
    case GateRx:
    case GateRy:
    case GateRz: {
      param_count = 1;
    } break;
    default: { param_count = 0; } break;
  }

  for (unsigned int i = 0; i < param_count; ++i) {
    PyList_Append((PyObject *)result, PyFloat_FromDouble(*(self->params + i)));
  }

  return (PyObject *)result;
}

static PyObject *qop_gate_get_matrix(QopGateObject *self) {
  PyListObject *result = (PyListObject *)PyList_New(0);
  for (unsigned int i = 0; i < 2; ++i) {
    PyListObject *row = (PyListObject *)PyList_New(0);
    for (unsigned int j = 0; j < 2; ++j) {
      double _Complex elem = self->gate.matrix[i][j];
      PyList_Append((PyObject *)row, (PyObject *)PyComplex_FromDoubles(
                                         creal(elem), cimag(elem)));
    }
    PyList_Append((PyObject *)result, (PyObject *)row);
  }
  return (PyObject *)result;
}