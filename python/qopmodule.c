#include "python/qopmodule.h"

// Initialization function of the module.
// Here we initialize not only the module, but also import the NumPy
// functions (via `import_array()`), and instantiate the internal error
// module exception type, `QopError`.
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

// Creates a new python object representing the circuit.
// This object is of type `QopCircuitObject`, and has internally stored
// a `Circuit` object, from the C qop library.
static PyObject *qop_create_circuit(PyTypeObject *type, PyObject *args,
                                    PyObject *kwds) {
  unsigned int qubit_count;

  // Parse the arguments
  // These should be (qubit_count)
  {
    int signed_qubit_count;
    char *kwarg_names[] = {"qubit_count", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i", kwarg_names,
                                     &signed_qubit_count)) {
      // An error in parsing the arguments occurred
      // Pass the error upstream
      return NULL;
    }

    // Ensure that given qubit count >= 0
    if (signed_qubit_count <= 0) {
      PyErr_SetString(PyExc_ValueError,
                      "qubit count must be greater than zero");
      return NULL;
    }

    // Save qubit count
    qubit_count = (unsigned int)signed_qubit_count;
  }

  // Create the C-qop Circuit object
  Circuit new_circuit;
  {
    Result init_r = circuit_init(&new_circuit, qubit_count);
    if (!init_r.valid) {
      // Could not initialize C-qop Circuit, pass the error upstream
      PyErr_SetString(QopError, init_r.content.error_details.reason);
      return NULL;
    }
  }

  // Create the vector that stores references to all *python* gate
  // objects this circuit references.
  // This is important because the circuit should be keeping a reference
  // to these gate objects, and these references should be freed when
  // the python circuit object is freed
  Vector gate_obj_refs;
  {
    Result init_r = vector_init(&gate_obj_refs, sizeof(PyObject *), 0);
    if (!init_r.valid) {
      // Could not initialize the `gate_obj_refs` vector;
      // the C-qop Circuit was already initialized at this stage and
      // must be freed
      circuit_free(&new_circuit);
      // Propagate the error upstream
      PyErr_SetString(QopError, init_r.content.error_details.reason);
      return NULL;
    }
  }

  // Create the actual python circuit object
  // This returns *a new reference*, and so it is safe to return
  // this new pointer at the end of the function without increasing
  // its reference.
  QopCircuitObject *self = (QopCircuitObject *)type->tp_alloc(type, 0);
  if (self == NULL) {
    // Could not create the python circuit object;
    // Release the already created associated circuit and refs vector
    circuit_free(&new_circuit);
    vector_free(&gate_obj_refs);
    // Propagate the error upstream
    return NULL;
  }

  // Set the properties of the new python object
  self->qubit_count = qubit_count;
  self->circuit = new_circuit;
  self->gate_obj_refs = gate_obj_refs;
  self->compacted_and_hardened = false;

  // (Reminder that there is no need to INCREF here, as tp_alloc
  // returns a new refrence)
  return (PyObject *)self;
}

// The function responsible for deallocating a python circuit object
// once it runs out of references.
// This function frees internally referenced memory that was allocated
// at initialization, and also DECREFs refrences to python gate objects
// that were added to it
static void qop_circuit_obj_dealloc(QopCircuitObject *obj) {
  // Free the C-qop circuit
  circuit_free(&obj->circuit);

  // DECREF gates that this circuit INCREFd
  {
    Iter refs_iter = vector_iter_create(&obj->gate_obj_refs);
    Option next;
    while ((next = iter_next(&refs_iter)).some) {
      PyObject *gate_pyobj_ptr = *(PyObject **)next.data;
      Py_DECREF(gate_pyobj_ptr);
    }
  }

  // Free the vector containing the references itself
  vector_free(&obj->gate_obj_refs);

  // Free the python object itself
  Py_TYPE(obj)->tp_free(obj);
}

// Creates a new python object representing a gate.
// This object is of type `QopGateType`, and contains a C-qop gate
// object inside.
static PyObject *qop_gate_create(PyTypeObject *type, PyObject *args,
                                 PyObject *kwds) {
  // Arguments coming in from python to be initialized
  char *identifier_str;
  PyObject *matrix_optional = NULL;
  PyObject *params_optional = NULL;

  // Parse arguments
  // We are expecting (identifier, $matrix?, $parameters?)
  char *kwarg_names[] = {"identifier", "matrix", "parameters", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|$OO", kwarg_names,
                                   &identifier_str, &matrix_optional,
                                   &params_optional)) {
    // An error occurred parsing the arguments; pass the error upstream
    return NULL;
  }

  // Match the gate string identifier to an enum equivalent
  GateId id = GateNoId;
  {
    // Note that strcmp returns the *difference* between strings,
    // so a perfect match returns `0`
    if (!strcmp(identifier_str, "c")) {  // Custom gate
      id = GateCustom;
    } else if (!strcmp(identifier_str, "i")) {  // Identity gate
      id = GateI;
    } else if (!strcmp(identifier_str, "x")) {  // X gate
      id = GateX;
    } else if (!strcmp(identifier_str, "y")) {  // Y gate
      id = GateY;
    } else if (!strcmp(identifier_str, "z")) {  // Z gate
      id = GateZ;
    } else if (!strcmp(identifier_str, "h")) {  // Hadammard gate
      id = GateH;
    } else if (!strcmp(identifier_str, "sqrtx")) {  // Sqrt(X) gate
      id = GateSqrtX;
    } else if (!strcmp(identifier_str, "t")) {  // T gate
      id = GateT;
    } else if (!strcmp(identifier_str, "s")) {  // S gate
      id = GateS;
    } else if (!strcmp(identifier_str, "rx")) {  // R_X gate
      id = GateRx;
    } else if (!strcmp(identifier_str, "ry")) {  // R_Y gate
      id = GateRy;
    } else if (!strcmp(identifier_str, "rz")) {  // R_Z gate
      id = GateRz;
    } else {  // Unrecognized string
      PyErr_SetString(PyExc_ValueError,
                      "unknown gate identifier; \
  use one of i,c,x,y,z,h,sqrtx,t,s,rx,ry,rz");
      return NULL;
    }
  }

  // Initialize the C-qop gate object
  Gate gate;
  unsigned int parameter_count = 0;
  double *params = NULL;
  {
    // The gate initialization result to be unwrap afterwards.
    // Because its value depends on the gate being initialized,
    // but its unwrapping is the same, it's pre-declared here
    Result gate_init_result;

    // Initialization will be slightly different for different gate
    // types
    switch (id) {
      case GateCustom: {
        // If a custom matrix was specified then the matrix argument
        // must be given as well
        if (matrix_optional == NULL) {
          // Nothing was allocated; we can just pass the error upstream
          PyErr_SetString(PyExc_ValueError,
                          "specified custom gate, but did not specify matrix");
          return NULL;
        }

        // Matrix was specified; for convenience we convert the given
        // object into a numpy matrix
        // The numpy matrix *steals a reference*! This means that the
        // `matrix_optional` reference is now a responsability of
        // `npy_matrix`, and since `matrix_optional` is a reference, then
        // `npy_matrix` *shouldn't* be dereferenced either.
        PyObject *npy_matrix =
            PyArray_FROMANY(matrix_optional, NPY_CDOUBLE, 0, 0,
                            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);

        if (npy_matrix == NULL) {
          // Something happened converting to a numpy matrix;
          // pass the error upstream
          return NULL;
        }

        // Validate the properties of the given matrix
        {
          int ndim = PyArray_NDIM((PyArrayObject *)npy_matrix);
          npy_intp *dims = PyArray_DIMS((PyArrayObject *)npy_matrix);

          if (ndim != 2 || dims[0] != 2 || dims[1] != 2) {
            // Single qubit gate matrix is malformed
            PyErr_SetString(
                PyExc_ValueError,
                "given custom matrix does not have correct 2x2 dimensions");
            return NULL;
          }

          // Copy the matrix into an array the C-qop gate object can own
          // Because we specified the alignment flags when initializing
          // the npy array, we can call memcpy on its data pointer
          double _Complex gate_matrix[2][2];
          memcpy(gate_matrix, PyArray_DATA((PyArrayObject *)npy_matrix),
                 GATE_SINGLE_QUBIT_SIZE);

          // Initialize the C-qop gate object, finally
          // Because its a custom matrix gate, it does not have a
          // reparameterization function.
          gate_init_result = gate_init_from_matrix(&gate, gate_matrix, NULL);
        }
      } break;
      case GateRx:
      case GateRy:
      case GateRz: {
        // All rotation gates have the same initialization
        // A rotation gate requires one parameter
        if (params_optional == NULL) {
          PyErr_SetString(
              PyExc_ValueError,
              "rotation gate requires a parameter, but none was given");
          return NULL;
        }

        // Convert the given parameters array into a numpy object, for
        // easier handling/validation
        // Reminder: this initialization *steals* the reference to the
        // argument, which is already borrowed, so the npy_params object
        // should not be DECREFd!
        PyObject *npy_params =
            PyArray_FROMANY(params_optional, NPY_DOUBLE, 0, 0,
                            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);

        // Validate npy_params
        {
          unsigned int ndims = PyArray_NDIM((PyArrayObject *)npy_params);
          npy_intp *dims = PyArray_DIMS((PyArrayObject *)npy_params);
          if (ndims != 1 || dims[0] != 1) {
            // The matrix does not have the correct shape; send error
            // upstream.
            PyErr_SetString(PyExc_ValueError,
                            "specified rotation gate but wrong number of "
                            "parameters (expecting 1)");
            return NULL;
          }
        }

        // "Move" the parameters to the heap, so that they become
        // responsibility of the C-qop gate object
        {
          void *mem = malloc(sizeof(double) * 1);
          if (mem == NULL) {
            PyErr_NoMemory();
            return NULL;
          }
          params = (double *)mem;
        }
        memcpy(params, PyArray_DATA((PyArrayObject *)npy_params),
               sizeof(double) * 1);
        parameter_count = 1;

        // Finally initialize the C-qop gate
        gate_init_result = gate_init_from_identifier(&gate, id, params);
      } break;
      default: {
        // All other gates are just initialized by their identifier
        // and have no parameters
        gate_init_result = gate_init_from_identifier(&gate, id, NULL);
      } break;
    }

    // "Unwrap" the gate initialization result
    if (!gate_init_result.valid) {
      // Free the parameters allocated memory; it might be NULL if
      // nothing was allocated, but that's ok
      free(params);

      // Send the error upstream
      PyErr_SetString(QopError, gate_init_result.content.error_details.reason);
      return NULL;
    }
  }

  // Finally, initialize the python gate object, being careful to free
  // the memory that was already initialized for the c-qop gate object
  // Again, `tp_alloc` creates a *new* reference, so there is no need to
  // INCREF before returning
  QopGateObject *self = (QopGateObject *)type->tp_alloc(type, 0);
  if (self == NULL) {
    // Something went wrong with instantiating the python gate object
    // Free the c-qop gate object & params in heap
    gate_free(&gate);
    free(params);

    // Send the error upstream
    return NULL;
  }

  // Set properties and return
  self->gate = gate;
  self->params = params;
  self->param_count = parameter_count;
  return (PyObject *)self;
}

// Dealocates a python gate object.
// This also releases the associated C-qop gate object, and the
// heap memory that was used to store the parameters
static void qop_gate_obj_dealloc(QopGateObject *obj) {
  gate_free(&obj->gate);
  free(obj->params);
  // Free the python object itself
  Py_TYPE(obj)->tp_free((PyObject *)obj);
}

// Add a new (python) gate object to a (python) circuit object.
// This handles the C-qop operations, as well as performing some
// sanity checks.
// Note that because the C-qop gate objects live in python gate
// objects, and these are in the heap, that there's not much concern
// with shifting addresses of the C-qop gates.
static PyObject *qop_circuit_add_gate(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds) {
  // Arguments to be parsed from python
  // Note that `gate_obj` is passed in by reference!
  // The appropriate INCREF should be made (so that it lives as long
  // as the circuit python object), and a DECREF on deallocation of
  // the circuit python object should be made too
  QopGateObject *gate_obj;
  unsigned int quibt;
  Option_Uint control = option_none_uint();

  // Parse the arguments
  {
    // "Raw" given values; may be invalid
    int arg_qubit = -1;
    int arg_control = -1;

    char *kwarg_names[] = {"gate", "qubit", "control", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!i|i", kwarg_names,
                                     &QopGateType, &gate_obj, &arg_qubit,
                                     &arg_control)) {
      // An error parsing the arguments occurred.
      // Send it upstream.
      // (No allocations were made so far)
      return NULL;
    }

    // Validate the values
    if (arg_qubit < 0 || (unsigned int)arg_qubit >= self->circuit.depth[0]) {
      PyErr_SetString(PyExc_ValueError, "invalid qubit");
      return NULL;
    }

    if (arg_qubit == arg_control ||
        (arg_control >= 0 &&
         (unsigned int)arg_control >= self->circuit.depth[0])) {
      PyErr_SetString(PyExc_ValueError, "invalid control");
      return NULL;
    }

    // Save valid values
    quibt = (unsigned int)arg_qubit;
    if (arg_control < 0)
      control = option_none_uint();
    else
      control = option_from_uint((unsigned int)arg_control);
  }

  // Got the arguments.

  // Add the C-qop gate to the C-qop circuit.
  // This is done before INCREF/DECREF so that if anything goes
  // wrong there's no need to handle those dangling references.
  {
    Result add_r =
        circuit_add_gate(&self->circuit, &gate_obj->gate, quibt, control);
    if (!add_r.valid) {
      // Something went wrong adding the C-qop gate.
      // Send error upstream
      PyErr_SetString(QopError, add_r.content.error_details.reason);
      return NULL;
    }
  }

  // INCREF the gate and save a reference to it, so that we can
  // DECREF it on deallocation
  // (We save the reference first, so that if an error occurs
  //  there's no need to DECREF)
  {
    // Push the object to the refs vector
    // Note that the vector is of **pointers**! so we'll be copying
    // `gate_obj` (and not `*gate_obj`) into the vector.
    // Since the vector memcopies data, this works out.
    Result push_r = vector_push(&self->gate_obj_refs, &gate_obj);
    if (!push_r.valid) {
      // Something went wrong pushing.
      PyErr_SetString(QopError, push_r.content.error_details.reason);
      return NULL;
    }

    // Successfully saved a reference to the gate, INCREF it
    Py_INCREF((PyObject *)gate_obj);
  }

  // Invalidate any compact/harden status
  self->compacted_and_hardened = false;

  // Return None (this requires an INCREF, because None is also a
  // PyObject)
  Py_INCREF(Py_None);
  return Py_None;
}

// Helper function for `parse_optimization_settings`
// Releases not only the vector, but all the `GateParameterization`s
// it contains.
static void free_reparams_vector(Vector *reparams_vector) {
  Iter param_iter = vector_iter_create(reparams_vector);
  Option next;
  while ((next = iter_next(&param_iter)).some) {
    GateParameterization *param = (GateParameterization *)next.data;
    optimizer_gate_param_free(param);
  }
}

// Helper function for `qop_circuit_optimize`.
// Turns a python dictionary object of structure
// ```
// {
//   {
//      'ada': {
//        'rho': value,
//        'epsilon': value
//      },
//      'optimize': {
//        'gates': [value, value, ...],
//        'deltas': [ [value], [value, value], ...],
//        'stop_at': value,
//        'max_iterations': value
//      }
//   }
// }
// ```
// into a complete C-qop `OptimizerSettings` object.
// The absence of some of these fields is handled.
// Exceptions that occur inside this function are immediately
// reported with `PyErr_SetString`, and `false` is returned,
// the caller being responsible to just return NULL.
// This is because the exception type may vary.
// If the function was successful, `true` is returned.
static bool parse_optimization_settings(
    PyObject *settings, QopCircuitObject *self, double _Complex *hamiltonian,
    Vector *reparams_vec, Vector *reparams_to_obj_pointer,
    OptimizerSettings *opt_settings, AdadeltaSettings *ada_settings) {
  // Default settings; the dictionary changes these properties.
  *ada_settings = optimizer_adadelta_get_default();
  double stop_at = 1e-4;
  Option_Uint max_iters = option_none_uint();
  bool reparam_given = false;  // If no reparam is given, create one for all
                               // known reparam'able gates

  // Initialize the vectors; these have to be freed on every exception.
  vector_init(reparams_vec, sizeof(GateParameterization), 0);
  vector_init(reparams_to_obj_pointer, sizeof(PyObject *), 0);

  // Very important note: PyDict_GetItemString returns a **borrowed**
  // reference! No need to Py_DECREF it!

  // Parse the dictionary
  if (settings != NULL) {
    // Is it a dictionary?
    if (!PyDict_Check(settings)) {
      PyErr_SetString(PyExc_ValueError, "settings object must be a dictionary");
      // Free the vectors
      free_reparams_vector(reparams_vec);
      vector_free(reparams_to_obj_pointer);

      // Report error
      return false;
    }

    // Try to get the ADADELTA subdictionary
    PyObject *ada_dict = PyDict_GetItemString(settings, "ada");
    if (ada_dict != NULL) {
      // Got an object for "ada"! Check it is dictionary
      if (!PyDict_Check(ada_dict)) {
        PyErr_SetString(PyExc_ValueError,
                        "settings.ada object must be a dictionary");
        // Free the vectors
        free_reparams_vector(reparams_vec);
        vector_free(reparams_to_obj_pointer);
        // Report error
        return false;
      }

      // Object *is* dictionary, grab properties
      // Rho property
      {
        PyObject *rho = PyDict_GetItemString(ada_dict, "rho");
        if (rho != NULL) {
          // rho property is specified, check that it's a number
          if (!PyNumber_Check(rho)) {
            // Got a non numeric rho
            PyErr_SetString(PyExc_ValueError,
                            "settings.ada.rho object must be a number");

            // Free the vectors
            free_reparams_vector(reparams_vec);
            vector_free(reparams_to_obj_pointer);

            // Report error
            return false;
          }
          // Confirmed rho is a number, parse to C-number
          ada_settings->rho = PyFloat_AsDouble(rho);
        }
      }
      // Epsilon property
      {
        PyObject *epsilon = PyDict_GetItemString(ada_dict, "epsilon");
        if (epsilon != NULL) {
          // epsilon property is specified, check that it's a number
          if (!PyNumber_Check(epsilon)) {
            // Got a non-numeric epsilon
            PyErr_SetString(PyExc_ValueError,
                            "settings.ada.epsilon object must be a number");

            // Free the vectors
            free_reparams_vector(reparams_vec);
            vector_free(reparams_to_obj_pointer);

            // Report error
            return false;
          }
          // Epsilon is valid, parse it
          ada_settings->epsilon = PyFloat_AsDouble(epsilon);
        }
      }
      // Done with ADADELTA settings dict
    }

    // Try to get "optimize" subdictionary
    PyObject *opt_dict = PyDict_GetItemString(settings, "optimize");
    if (opt_dict != NULL) {
      // Found a "optimize" object, confirm that it's a dictionary object
      if (!PyDict_Check(opt_dict)) {
        PyErr_SetString(PyExc_ValueError,
                        "settings.optimize object must be dictionary");

        // Free vectors
        free_reparams_vector(reparams_vec);
        vector_free(reparams_to_obj_pointer);

        // Report error
        return false;
      }

      // "optimize" *is* a dictionary
      // Try to get subfields

      // Try to parse `stop_at`
      {
        PyObject *stop_at_obj = PyDict_GetItemString(opt_dict, "stop_at");
        if (stop_at_obj != NULL) {
          // "stop_at" is defined; check type
          if (!PyNumber_Check(stop_at_obj)) {
            // Got non-numeric stop_at
            PyErr_SetString(PyExc_ValueError,
                            "settings.optimize.stop_at must be a number");

            // Free vectors
            free_reparams_vector(reparams_vec);
            vector_free(reparams_to_obj_pointer);
            return false;
          }

          // "stop_at" is a number
          stop_at = PyFloat_AsDouble(stop_at_obj);
        }
      }

      // Try to parse `max_iterations`
      {
        PyObject *max_iters_obj =
            PyDict_GetItemString(opt_dict, "max_iterations");
        if (max_iters_obj != NULL) {
          // Check that the value is numeric
          if (!PyNumber_Check(max_iters_obj)) {
            // Got a non-numeric max_iters
            PyErr_SetString(
                PyExc_ValueError,
                "settings.optimize.max_iterations must be a number");

            // Free vectors
            free_reparams_vector(reparams_vec);
            vector_free(reparams_to_obj_pointer);
            return false;
          }

          // "max_iterations" is a number
          // A value <= 0 signifies no maximum.
          {
            int max_iters_int = (int)PyFloat_AsDouble(max_iters_obj);
            if (max_iters_int <= 0)
              max_iters = option_none_uint();
            else
              max_iters = option_from_uint((unsigned int)max_iters_int);
          }
        }
      }

      // Try to parse `gates`/`deltas`
      PyObject *gates = PyDict_GetItemString(opt_dict, "gates");
      PyObject *deltas = PyDict_GetItemString(opt_dict, "deltas");

      // Both `gates` and `deltas` must be defined or undefined
      if ((gates != NULL) ^ (deltas != NULL)) {
        PyErr_SetString(PyExc_ValueError,
                        "settings.optimize.gates and settings.optimize.deltas "
                        "must both be defined or undefined");
        // Free vectors
        free_reparams_vector(reparams_vec);
        vector_free(reparams_to_obj_pointer);

        // Report error
        return false;
      }

      // Now, if `gates` is defined, we are guaranteed that `deltas`
      // is too
      if (gates != NULL) {
        // Check their types
        // Note that it's not enough to call `PyIter_Check`, as this
        // will fail for lists and tuples. Whether this is by design
        // or a bug is unknown, but the iterator API seems to work
        // fine with these two types as well.
        {
          // Check gates type
          if (!(PyList_Check(gates) || PyTuple_Check(gates) ||
                PyIter_Check(gates))) {
            // The gates object is not iterable
            PyErr_SetString(PyExc_ValueError,
                            "settings.optimize.gates object must be iterable");
            // DECREF appropriate objects
            // Use XDECREF for gates and deltas, since we don't know
            // which are NULL
            Py_XDECREF(gates);
            Py_XDECREF(deltas);

            // Free vectors
            free_reparams_vector(reparams_vec);
            vector_free(reparams_to_obj_pointer);

            // Report error
            return false;
          }

          // Check deltas type
          if (!(PyList_Check(deltas) || PyTuple_Check(deltas) ||
                PyIter_Check(deltas))) {
            // The gates object is not iterable
            PyErr_SetString(PyExc_ValueError,
                            "settings.optimize.deltas object must be iterable");
            // DECREF appropriate objects
            // Use XDECREF for gates and deltas, since we don't know
            // which are NULL
            Py_XDECREF(gates);
            Py_XDECREF(deltas);

            // Free vectors
            free_reparams_vector(reparams_vec);
            vector_free(reparams_to_obj_pointer);

            // Report error
            return false;
          }
        }

        // Iterate over `gates` and `deltas` simultaneously, using
        // the python iteration API.
        // To consider:
        //  The iterator itself is given by PyObject_GetIter, which
        //   returns a new reference, so the iterator must be DECREFed
        //   when done
        //  PyIter_Next returns a new reference, so the yielded object
        //   must be DECREFed when done
        {
          PyObject *gates_iter = PyObject_GetIter(gates);
          PyObject *delta_coll_iter = PyObject_GetIter(deltas);

          // Check if something went wrong creating either iterator
          if (gates_iter == NULL || delta_coll_iter == NULL) {
            // Dereference the iterators; because we don't know
            // which are NULL (they might both be, we XDECREF)
            Py_XDECREF(gates_iter);
            Py_XDECREF(delta_coll_iter);
            // DECREF other appropriate objects
            Py_XDECREF(gates);
            Py_XDECREF(deltas);

            // Free vectors
            free_reparams_vector(reparams_vec);
            vector_free(reparams_to_obj_pointer);

            // Just send the error upstream
            return false;
          }

          // Iterate over both iterators
          {
            PyObject *next_gate;
            PyObject *next_delta_collection;

            while ((next_gate = PyIter_Next(gates_iter)) != NULL) {
              next_delta_collection = PyIter_Next(delta_coll_iter);

              // Check if there's a mismatch between the two iterators
              if (next_delta_collection == NULL) {
                PyErr_SetString(PyExc_ValueError,
                                "size mismatch between settings.optimize.gates "
                                "and settings.optimize.deltas");

                // Decref appropriate objects
                Py_DECREF(next_gate);
                Py_DECREF(gates_iter);
                Py_DECREF(delta_coll_iter);

                // Free vectors
                free_reparams_vector(reparams_vec);
                vector_free(reparams_to_obj_pointer);

                // Report error
                return false;
              }

              // Check that the gate object is actually a python qop gate
              if (!PyObject_TypeCheck(next_gate, &QopGateType)) {
                // `next_gate` is not a python qop gate
                PyErr_SetString(
                    PyExc_ValueError,
                    "settings.optimize.gates element must be qop.Gate");

                // Decref all relevant objects
                Py_DECREF(next_gate);
                Py_DECREF(next_delta_collection);
                Py_DECREF(gates_iter);
                Py_DECREF(delta_coll_iter);

                // Free vectors
                free_reparams_vector(reparams_vec);
                vector_free(reparams_to_obj_pointer);

                // Report error
                return false;
              }

              // We have a delta collection object for this gate.
              // Check that it can be iterated over
              // Again, we must check for all of (list, tuple, iter)
              if (!(PyList_Check(next_delta_collection) ||
                    PyTuple_Check(next_delta_collection) ||
                    PyIter_Check(next_delta_collection))) {
                // The delta collection is not iterable
                PyErr_SetString(
                    PyExc_ValueError,
                    "settings.optimize.deltas element must be iterable");

                // Decref appropriate objects
                Py_DECREF(next_gate);
                Py_DECREF(next_delta_collection);
                Py_DECREF(gates_iter);
                Py_DECREF(delta_coll_iter);

                // Free vectors
                free_reparams_vector(reparams_vec);
                vector_free(reparams_to_obj_pointer);

                // Report error
                return false;
              }

              // Confirmed that we can iterate over the delta collection.
              // Iterate over the deltas; move these into a vector
              // so that it's easier to then plug them into a
              // GateParameterization.
              // The tradeoff is that this vector must then be freed
              // on every error
              Vector deltas_vec;
              {
                Result init_r = vector_init(&deltas_vec, sizeof(double), 0);
                if (!init_r.valid) {
                  // Internall error initializing the vector;
                  // send it upstream
                  PyErr_SetString(QopError,
                                  init_r.content.error_details.reason);

                  // Decref appropriate objects
                  Py_DECREF(next_gate);
                  Py_DECREF(next_delta_collection);
                  Py_DECREF(gates_iter);
                  Py_DECREF(delta_coll_iter);

                  // Free vectors
                  free_reparams_vector(reparams_vec);
                  vector_free(reparams_to_obj_pointer);

                  // Report error
                  return false;
                }
              }
              // Do actual iteration
              {
                PyObject *deltas_iter = PyObject_GetIter(next_delta_collection);
                if (deltas_iter == NULL) {
                  // Something went wrong creating the deltas iterator
                  // Just send the error upstream

                  // Decref appropriate objects
                  Py_DECREF(next_gate);
                  Py_DECREF(next_delta_collection);
                  Py_DECREF(gates_iter);
                  Py_DECREF(delta_coll_iter);

                  // Free vectors
                  free_reparams_vector(reparams_vec);
                  vector_free(reparams_to_obj_pointer);
                  vector_free(&deltas_vec);
                  return false;
                }

                // Perform actual iteration
                PyObject *delta;
                while ((delta = PyIter_Next(deltas_iter)) != NULL) {
                  // Check that the value is a numeric value
                  if (!PyNumber_Check(delta)) {
                    // Got non-numeric delta
                    PyErr_SetString(PyExc_ValueError, "delta must be a number");

                    // Decref appropriate objects
                    Py_DECREF(delta);
                    Py_DECREF(next_gate);
                    Py_DECREF(next_delta_collection);
                    Py_DECREF(gates_iter);
                    Py_DECREF(delta_coll_iter);
                    Py_DECREF(deltas_iter);

                    // Free vectors
                    free_reparams_vector(reparams_vec);
                    vector_free(reparams_to_obj_pointer);
                    vector_free(&deltas_vec);

                    // Send error upstream
                    return false;
                  }

                  // Confirmed that delta is numeric
                  double delta_num = PyFloat_AsDouble(delta);

                  // Push to deltas vector
                  {
                    Result push_r = vector_push(&deltas_vec, &delta_num);

                    if (!push_r.valid) {
                      // Something went wrong while pushing
                      // Got non-numeric delta
                      PyErr_SetString(QopError,
                                      push_r.content.error_details.reason);

                      // Decref appropriate objects
                      Py_DECREF(delta);
                      Py_DECREF(next_gate);
                      Py_DECREF(next_delta_collection);
                      Py_DECREF(gates_iter);
                      Py_DECREF(delta_coll_iter);
                      Py_DECREF(deltas_iter);

                      // Free vectors
                      free_reparams_vector(reparams_vec);
                      vector_free(reparams_to_obj_pointer);
                      vector_free(&deltas_vec);

                      // Send error upstream
                      return false;
                    }
                  }

                  // An iter next creates a new reference
                  Py_DECREF(delta);
                }

                // Done iterating on deltas collection; release the iterator
                Py_DECREF(deltas_iter);
              }

              // Pushed all of the deltas into the delta_vec vector
              // Create the C-qop GateParameterization
              GateParameterization gate_param;
              {
                QopGateObject *gate = (QopGateObject *)next_gate;
                // Note that according to documentation, both deltas
                // and params are memcpyd into the gate parameterization,
                // so there's no concerns about lifetime
                Result init_r = optimizer_gate_param_init(
                    &gate_param, &gate->gate, deltas_vec.size, gate->params,
                    deltas_vec.data);

                if (!init_r.valid) {
                  // Something went wrong
                  // Send the error upstream
                  PyErr_SetString(QopError,
                                  init_r.content.error_details.reason);

                  // Decref appropriate objects
                  Py_DECREF(next_gate);
                  Py_DECREF(next_delta_collection);
                  Py_DECREF(gates_iter);
                  Py_DECREF(delta_coll_iter);

                  // Free vectors
                  free_reparams_vector(reparams_vec);
                  vector_free(reparams_to_obj_pointer);
                  vector_free(&deltas_vec);

                  // Report error
                  return false;
                }
              }

              // Parameterization was defined; add it to the list!
              // Push reparams_vec
              {
                Result push_r = vector_push(reparams_vec, &gate_param);

                if (!push_r.valid) {
                  // Something went wrong
                  // Send the error upstream
                  PyErr_SetString(QopError,
                                  push_r.content.error_details.reason);

                  // Decref appropriate objects
                  Py_DECREF(next_gate);
                  Py_DECREF(next_delta_collection);
                  Py_DECREF(gates_iter);
                  Py_DECREF(delta_coll_iter);

                  // Free vectors
                  free_reparams_vector(reparams_vec);
                  vector_free(reparams_to_obj_pointer);
                  vector_free(&deltas_vec);

                  // Report error
                  return false;
                }
              }
              // Push reparams_to_obj_pointer
              {
                QopGateObject *gate = (QopGateObject *)next_gate;
                Result push_r = vector_push(reparams_to_obj_pointer, &gate);

                if (!push_r.valid) {
                  // Something went wrong
                  // Send the error upstream
                  PyErr_SetString(QopError,
                                  push_r.content.error_details.reason);

                  // Decref appropriate objects
                  Py_DECREF(next_gate);
                  Py_DECREF(next_delta_collection);
                  Py_DECREF(gates_iter);
                  Py_DECREF(delta_coll_iter);

                  // Free vectors
                  free_reparams_vector(reparams_vec);
                  vector_free(reparams_to_obj_pointer);
                  vector_free(&deltas_vec);

                  // Report error
                  return false;
                }
              }

              // A gate parameterization was fully defined by this point;
              // flag that
              reparam_given = true;

              // Decref the objects from this loop
              Py_DECREF(next_gate);
              Py_DECREF(next_delta_collection);
            }
          }

          // Done using the iterators
          Py_DECREF(gates_iter);
          Py_DECREF(delta_coll_iter);
        }  // Finished iterating over `gates` and `deltas`

      }  // Finished parsing `gates`

    }  // Finished parsing "optimize" dictionary
  }    // Finished parsing settings dictionary

  // Create parameterizations automatically if none was defined
  if (!reparam_given) {
    // Iterate over the python gate objects that this circuit stores
    // references to.
    // If we know how to parameterize them, create a suitable parameterization
    // and add that to the parameterizations vector
    Iter refd_gates_iter = vector_iter_create(&self->gate_obj_refs);
    Option next;
    while ((next = iter_next(&refd_gates_iter)).some) {
      QopGateObject *py_gate = *(QopGateObject **)next.data;

      switch (py_gate->gate.id) {
        case GateRx:
        case GateRy:
        case GateRz: {
          // We know how to parameterize a rotation gate.
          GateParameterization param;  // The new param. to add
          double delta[] = {0.1};
          {
            Result init_r = optimizer_gate_param_init(&param, &py_gate->gate, 1,
                                                      py_gate->params, delta);
            if (!init_r.valid) {
              // Something went wrong initializing the param.
              // Send error upstream
              PyErr_SetString(QopError, init_r.content.error_details.reason);

              // Free vectors
              free_reparams_vector(reparams_vec);
              vector_free(reparams_to_obj_pointer);
              return false;
            }
          }
          // Successfully initialized parameterization; push
          {
            Result push_r = vector_push(reparams_vec, &param);
            if (!push_r.valid) {
              // Something went wrong pushing
              // Send the error upstream
              PyErr_SetString(QopError, push_r.content.error_details.reason);

              // Free vectors
              free_reparams_vector(reparams_vec);
              vector_free(reparams_to_obj_pointer);
              return false;
            }
          }
          {
            Result push_r = vector_push(reparams_to_obj_pointer, &py_gate);
            if (!push_r.valid) {
              // Something went wrong pushing
              // Send the error upstream
              PyErr_SetString(QopError, push_r.content.error_details.reason);

              // Free vectors
              free_reparams_vector(reparams_vec);
              vector_free(reparams_to_obj_pointer);
              return false;
            }
          }
        } break;
        default:
          continue;
      }  // end switch(id)
    }    // done iterating over gates
  }      // done automatically creating reparams

  // Sanity check that there's something to optimize!
  if (reparams_vec->size == 0) {
    PyErr_SetString(QopError, "nothing to optimize");
    free_reparams_vector(reparams_vec);
    vector_free(reparams_to_obj_pointer);
    return false;
  }

  // Build the final OptimizerSettings
  {
    Result init_r = optimizer_settings_init(
        opt_settings, &self->circuit, hamiltonian, stop_at, reparams_vec->data,
        reparams_vec->size, max_iters);

    if (!init_r.valid) {
      // Something went wrong initializing the optimizer settings
      // Send the error upstream
      PyErr_SetString(QopError, init_r.content.error_details.reason);
      free_reparams_vector(reparams_vec);
      vector_free(reparams_to_obj_pointer);
      return false;
    }
  }

  // Success!
  return true;
}

// Perform the optimization routine (using the C-qop optimizer) on
// the C-qop circuit contained in the python circuit object.
// This abstracts a lot of steps from the end user, including
// compacting and hardening the C-qop circuit if needed, creating
// `GateParameterization`s if none were given, and generally parsing
// optimization settings from a more user friendly representation.
// Also, because the whole optimization runs under this function's
// scope, the memory scope is clear.
static PyObject *qop_circuit_optimize(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds) {
  // Arguments to be parsed
  PyObject *npy_hamiltonian;
  PyObject *settings = NULL;
  int writer_fd = -1;

  // Parse incoming arguments
  {
    PyObject *hamiltonian_obj = NULL;
    PyObject *given_writer = NULL;

    char *kwarg_names[] = {"hamiltonian", "settings", "writer", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O!O", kwarg_names,
                                     &hamiltonian_obj, &PyDict_Type, &settings,
                                     &given_writer)) {
      // Could not parse the arguments
      PyErr_SetString(PyExc_ValueError,
                      "could not parse arguments to circuit optimization; "
                      "expected (hamiltonian, settings?, writer?)");
      return NULL;
    }

    // Convert hamiltonian object to numpy array
    // Remember: creating a npy array steals a reference! and the
    //  argument is a borrow, so it should **not** be DECREFd!
    npy_hamiltonian =
        PyArray_FROMANY(hamiltonian_obj, NPY_CDOUBLE, 2, 2,
                        NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);
    if (npy_hamiltonian == NULL) {
      // Could not numpy parse the hamiltonian py object
      PyErr_SetString(PyExc_ValueError,
                      "could not interpret hamiltonian as matrix");
      return NULL;
    }

    // Sanity check the hamiltonian
    {
      int dims = PyArray_NDIM((PyArrayObject *)npy_hamiltonian);
      npy_intp *size = PyArray_DIMS((PyArrayObject *)npy_hamiltonian);
      if (dims != 2 || size[0] != size[1]) {
        PyErr_SetString(PyExc_ValueError,
                        "hamiltonian is malformed; must be square 2D array");
        return NULL;
      }

      if (size[0] != (1 << self->circuit.depth[0])) {
        PyErr_SetString(
            PyExc_ValueError,
            "hamiltonian is malformed; does not match circuit state size");
        return NULL;
      }
    }

    // Prepare the results writer if given
    if (given_writer != NULL) {
      // Ensure that the given object has a file descriptor
      int fdescriptor = PyObject_AsFileDescriptor(given_writer);
      if (fdescriptor < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid writer object");
        return NULL;
      }
      writer_fd = fdescriptor;
    }
  }

  // Prepare the circuit
  if (!self->compacted_and_hardened) {
    // Compact
    {
      Result comp_r = circuit_compact(&self->circuit);
      if (!comp_r.valid) {
        PyErr_SetString(QopError, comp_r.content.error_details.reason);
        return NULL;
      }
    }
    // Harden
    {
      Result hard_r = circuit_harden(&self->circuit);
      if (!hard_r.valid) {
        PyErr_SetString(QopError, hard_r.content.error_details.reason);
        return NULL;
      }
    }
    self->compacted_and_hardened = true;
  }

  // Create the optimizer settings
  Vector reparams_vec;
  Vector reparams_to_obj_pointer;  // (reparam index) -> (associated python
                                   //                     object gate)
  OptimizerSettings opt_settings;
  AdadeltaSettings ada_settings;
  {
    bool success = parse_optimization_settings(
        settings, self,
        (double _Complex *)PyArray_DATA((PyArrayObject *)npy_hamiltonian),
        &reparams_vec, &reparams_to_obj_pointer, &opt_settings, &ada_settings);
    if (!success) {
      // Had some error in creating the optimization settings.
      // The error has already been set, so we just send the error
      // upstream.
      return NULL;
    }
  }

  // Create optimizer
  Optimizer optimizer;
  {
    Result init_r = optimizer_init(&optimizer, opt_settings, ada_settings);
    if (!init_r.valid) {
      // Something went wrong with initializing the optimizer object
      optimizer_settings_free(&opt_settings);
      free_reparams_vector(&reparams_vec);
      vector_free(&reparams_to_obj_pointer);
      return NULL;
    }
  }

  // Perform actual optimization
  OptimizationResult optimization_result;
  {
    if (writer_fd > 0) {
      optimization_result =
          optimizer_optimize(&optimizer, qop_write_parameter, qop_write_energy,
                             (void *)&writer_fd);
    } else {
      optimization_result = optimizer_optimize(&optimizer, NULL, NULL, NULL);
    }
    if (!optimization_result.valid) {
      // Optimization failed for some reason
      PyErr_SetString(QopError,
                      optimization_result.content.error_details.reason);
      optimizer_settings_free(&opt_settings);
      free_reparams_vector(&reparams_vec);
      vector_free(&reparams_to_obj_pointer);
      return NULL;
    }
  }

  // Save results of optimization into tuple of tuples;
  // First make a C-list, then convert that to a python tuple of tuples

  // Vector initialization
  Vector results_vector;
  {
    Result init_r = vector_init(&results_vector, sizeof(PyTupleObject *), 0);
    if (!init_r.valid) {
      // Initializing the vector failed for some reason
      PyErr_SetString(QopError, init_r.content.error_details.reason);
      optimizer_settings_free(&opt_settings);
      free_reparams_vector(&reparams_vec);
      vector_free(&reparams_to_obj_pointer);
      return NULL;
    }
  }

  // Push parameters to vector
  {
    Iter reparams_iter = vector_iter_create(&reparams_vec);
    Option next;
    while ((next = iter_next(&reparams_iter)).some) {
      GateParameterization *reparam = (GateParameterization *)next.data;

      // Create a python tuple with enough size for all the reparam's
      // parameters
      PyObject *subtuple = PyTuple_New(reparam->param_count);

      // Copy the new parameters to the python object
      {
        // Grab the corresponding gate
        QopGateObject *py_gate =
            *(QopGateObject **)(reparams_to_obj_pointer.data +
                                (reparams_iter.position - 1) *
                                    sizeof(QopGateObject *));
        // Copy the parameters there
        memcpy(py_gate->params, reparam->params,
               sizeof(double) * reparam->param_count);
      }

      // Iterate over the subparameters; note that creating a python float
      // returns a *new* reference, and setting an element on a pytuple
      // *steals* that reference, so we don't need to DECREF or INCREF
      for (unsigned int i = 0; i < reparam->param_count; ++i) {
        double final_param = *(reparam->params + i);
        PyObject *py_param = PyFloat_FromDouble(final_param);
        PyTuple_SetItem(subtuple, i, py_param);
      }

      // Push **the subtuple pointer** to the results vector, to be assembled
      // to a tuple later
      {
        Result push_r = vector_push(&results_vector, &subtuple);
        if (!push_r.valid) {
          // Pushing failed for some reason
          // Send the error upstream
          PyErr_SetString(QopError, push_r.content.error_details.reason);
          // Free all the tuples that were already pushed
          {
            Iter results_iter = vector_iter_create(&results_vector);
            Option next;
            while ((next = iter_next(&results_iter)).some) {
              // next.data is a pointer to a pointer
              Py_DECREF(*(PyObject **)next.data);
            }
          }
          // Free stuff in general
          optimizer_settings_free(&opt_settings);
          free_reparams_vector(&reparams_vec);
          vector_free(&reparams_to_obj_pointer);
          vector_free(&results_vector);

          return NULL;
        }
      }
    }
  }

  // Assemble final results tuple
  PyObject *final_result;
  {
    final_result = PyTuple_New(results_vector.size);

    // Again, important note:
    // PyTuple_SetItem *steals* a reference, and when creating the
    // subtuple the reference is a *new* one, so there is no need
    // for any kind of INC/DECREFs here
    {
      Iter subtuple_iter = vector_iter_create(&results_vector);
      Option next;
      while ((next = iter_next(&subtuple_iter)).some) {
        PyObject *subtuple = *(PyObject **)next.data;
        PyTuple_SetItem(final_result, subtuple_iter.position - 1, subtuple);
      }
    }

    // Free the results_vector vector.
    // Note that this is safe, because it was a vector of pointers
    vector_free(&results_vector);
  }

  // Cleanup
  optimizer_settings_free(&opt_settings);
  free_reparams_vector(&reparams_vec);
  vector_free(&reparams_to_obj_pointer);

  // Return result
  // PyTuple_Pack creates a new reference, and so does PyBool_FromLong
  // so the references also work out here
  return PyTuple_Pack(2, final_result,
                      PyBool_FromLong(optimization_result.quit_on_max_iter));
}

// Return a tuple of the python gate objects this circuit references
static PyObject *qop_circuit_get_gates(QopCircuitObject *self) {
  PyObject *gates_tuple = PyTuple_New(self->gate_obj_refs.size);
  if (gates_tuple == NULL) {
    // Something went wrong initializing the tuple,
    // send the error upstream
    return NULL;
  }

  Iter ref_gates_iter = vector_iter_create(&self->gate_obj_refs);
  Option next;
  while ((next = iter_next(&ref_gates_iter)).some) {
    QopGateObject *gate_obj = *(QopGateObject **)next.data;
    // This gate object will now be referenced in the tuple,
    // as well as in the circuit, and PyTuple_SetItem steals a
    // reference, so we need to INCREF
    Py_INCREF(gate_obj);

    // Set the gate in the tuple
    PyTuple_SetItem(gates_tuple, ref_gates_iter.position - 1,
                    (PyObject *)gate_obj);
  }

  // The list was created with a new reference, so we can return
  // without further INCREFing it
  return gates_tuple;
}

// Run a state through the C-qop circuit associated with this python
// circuit object.
// The result is converted to a tuple.
static PyObject *qop_circuit_run(QopCircuitObject *self, PyObject *args,
                                 PyObject *kwds) {
  // Parse incoming state
  double _Complex state_in[1U << self->qubit_count];
  memset(state_in, 0, sizeof(double _Complex) * (1U << self->qubit_count));

  {
    // Parse argument tuple
    PyObject *state_in_obj;
    char *kwarg_names[] = {"state_in", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwarg_names,
                                     &state_in_obj)) {
      PyErr_SetString(PyExc_ValueError, "expecting state_in argument");
      return NULL;
    }

    // Convert to numpy array
    PyObject *as_arr =
        PyArray_FROMANY(state_in_obj, NPY_CDOUBLE, 1, 1,
                        NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED);
    if (as_arr == NULL) {
      PyErr_SetString(PyExc_ValueError,
                      "cannot interpret state_in as 1D vector");
      return NULL;
    }

    // Validate size of the array
    unsigned int vec_size = PyArray_SIZE((PyArrayObject *)as_arr);
    if (vec_size != (1U << self->qubit_count)) {
      PyErr_SetString(PyExc_ValueError,
                      "given input vector does not match size of circuit");
      return NULL;
    }

    // Validate 2-norm of the array
    double norm = 0;
    for (unsigned int i = 0; i < vec_size; ++i) {
      double _Complex elem =
          *(double _Complex *)PyArray_GETPTR1((PyArrayObject *)as_arr, i);
      norm += creal(elem * conj(elem));
    }
    if (fabs(norm - 1.) > 1e-6) {
      PyErr_SetString(
          PyExc_ValueError,
          "the given in state's 2-norm differs of 1 by more than 1E-6");
      return NULL;
    }

    // Copy contents to in_state
    memcpy(state_in, PyArray_DATA((PyArrayObject *)as_arr),
           sizeof(double _Complex) * vec_size);

    // Note that we **do not** DECREF `state_in_obj` or `as_arr` here,
    // since `state_in_obj` is passed in by reference, and `as_arr` steals
    // that reference
  }

  // Compact and harden the circuit if needed
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

  // Actually run the simulation!
  {
    Result run_r = circuit_run(&self->circuit, &state_in);
    if (!run_r.valid) {
      PyErr_SetString(QopError, run_r.content.error_details.reason);
      return NULL;
    }
  }

  // Put the results into a python tuple
  PyObject *final_result;
  {
    final_result = PyTuple_New(1U << self->qubit_count);

    for (unsigned int i = 0; i < (1U << self->qubit_count); ++i) {
      // References also check out here; the FromDoubles creates
      // a new reference, which is immediately stolen by the tuple
      double _Complex elem = state_in[i];
      PyTuple_SetItem(final_result, i,
                      PyComplex_FromDoubles(creal(elem), cimag(elem)));
    }
  }

  // final_result was created with a new reference, so there's no need
  // to INCREF before returning
  return final_result;
}

// Set new parameters for a gate.
// This is a very direct connection to the C-qop gate object, with
// a few sanity checks.
static PyObject *qop_gate_reparameterize(QopGateObject *self, PyObject *args,
                                         PyObject *kwds) {
  // Can we reparameterize the gate at all
  if (self->gate.reparamFn == NULL) {
    PyErr_SetString(QopError, "gate is not parameterized");
    return NULL;
  }

  // Parse the new given parameters
  Vector params;
  {
    PyObject *params_arg;

    char *kwarg_names[] = {"parameters", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwarg_names,
                                     &params_arg)) {
      PyErr_SetString(PyExc_ValueError, "expecting parameters iterable");
      return NULL;
    }

    // Validate that we can iterate over the parameters
    if (!(PyList_Check(params_arg) || PyTuple_Check(params_arg) ||
          PyIter_Check(params_arg))) {
      PyErr_SetString(PyExc_ValueError, "expecting iterable parameters");
      return NULL;
    }

    // Initialize the vector where we'll copy the new parameters
    {
      Result init_r = vector_init(&params, sizeof(double), 1);
      if (!init_r.valid) {
        // Could not initialize vector for some reason
        PyErr_SetString(QopError, init_r.content.error_details.reason);
        return NULL;
      }
    }

    // Iterate over the python object
    {
      PyObject *py_args_iter = PyObject_GetIter(params_arg);
      if (py_args_iter == NULL) {
        // Could not create a new iterator; send the error upstream
        vector_free(&params);
        return NULL;
      }

      PyObject *next;
      while ((next = PyIter_Next(py_args_iter)) != NULL) {
        // Validate that the argument is a number
        if (!PyNumber_Check(next)) {
          // Not a number!
          PyErr_SetString(PyExc_ValueError, "got unexpected non-numeric value");
          Py_DECREF(py_args_iter);
          vector_free(&params);
          return NULL;
        }

        // Parse value
        double param = PyFloat_AsDouble(next);

        // Push value
        {
          Result push_r = vector_push(&params, &param);
          if (!push_r.valid) {
            // Pushing failed for some reason,
            // send error upstream
            PyErr_SetString(QopError, push_r.content.error_details.reason);
            Py_DECREF(py_args_iter);
            vector_free(&params);
            return NULL;
          }
        }

        // Done with the iteration yield object
        Py_DECREF(next);
      }

      // Done with the iterator
      Py_DECREF(py_args_iter);
    }
  }

  // Validate parameter count
  if (params.size != self->param_count) {
    PyErr_SetString(PyExc_ValueError,
                    "got wrong number of parameters in reparameterization");
    vector_free(&params);
    return NULL;
  }

  // Actually perform reparameterization
  self->gate.reparamFn(&self->gate.matrix, params.data);

  // Save the parameters
  memcpy(self->params, params.data, params.size * sizeof(double));

  // Cleanup
  vector_free(&params);

  // Return None; still needs an INCREF!
  Py_INCREF(Py_None);
  return Py_None;
}

// Return a tuple with the current parameters of the given gate
static PyObject *qop_gate_get_parameters(QopGateObject *self) {
  PyObject *result = PyTuple_New(self->param_count);
  for (unsigned int i = 0; i < self->param_count; ++i) {
    // References work out, between new-ing and stealing
    PyTuple_SetItem(result, i, PyFloat_FromDouble(*(self->params + i)));
  }
  return result;
}

// Return a numpy style representation of the matrix; simply a
// nested (pseudo-row major) tuple of the matrix elements.
static PyObject *qop_gate_get_matrix(QopGateObject *self) {
  PyObject *result = PyTuple_New(2);
  for (unsigned int i = 0; i < 2; ++i) {
    PyObject *row = PyTuple_New(2);
    for (unsigned int j = 0; j < 2; ++j) {
      double _Complex elem = self->gate.matrix[i][j];
      PyTuple_SetItem(row, j, PyComplex_FromDoubles(creal(elem), cimag(elem)));
    }
    PyTuple_SetItem(result, i, row);
  }
  return result;
}

// Callback function for the optimizer, where the context given is the
// file descriptor. This should write the parameter to the file
static void qop_write_parameter(unsigned int flat_index, double param,
                                void *context) {
  int file_desc = *(int *)context;
  char *stringed;
  int byte_count = asprintf(&stringed, "%e ", param);
  if (byte_count > 0) {
    write(file_desc, stringed, byte_count);
  }
  free(stringed);
}

// Callback function for the optimizer, where the context given is the
// file descriptor. This should write the energy to the file.
static void qop_write_energy(_Complex double energy, void *context) {
  int file_desc = *(int *)context;
  char *stringed;
  int byte_count = asprintf(&stringed, "%e %e\n", creal(energy), cimag(energy));
  if (byte_count > 0) {
    write(file_desc, stringed, byte_count);
  }
  free(stringed);
}
