#ifndef QOP_QOPMODULE_H_
#define QOP_QOPMODULE_H_

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdbool.h>
#include <structmember.h>
#include "include/circuit.h"
#include "include/gate.h"
#include "include/iter.h"
#include "include/optimizer.h"
#include "include/option.h"
#include "include/vector.h"

static PyObject *QopError;

typedef struct QopCircuitObject {
  PyObject_HEAD Circuit circuit;
  unsigned int qubit_count;
  Vector gate_obj_refs;
  bool compacted_and_hardened;
} QopCircuitObject;

typedef struct QopGateObject {
  PyObject_HEAD Gate gate;
  double *params;
  unsigned int param_count;
} QopGateObject;

static PyObject *qop_create_circuit(PyTypeObject *type, PyObject *args,
                                    PyObject *kwds);
static void qop_circuit_obj_dealloc(QopCircuitObject *obj);
static PyObject *qop_gate_create(PyTypeObject *type, PyObject *args,
                                 PyObject *kwds);
static void qop_gate_obj_dealloc(QopGateObject *obj);

static PyObject *qop_circuit_add_gate(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds);
static bool parse_optimization_settings(
    PyObject *settings, QopCircuitObject *self, double _Complex *hamiltonian,
    Vector *reparams_vec, Vector *reparams_to_obj_pointer,
    OptimizerSettings *opt_settings, AdadeltaSettings *ada_settings);
static PyObject *qop_circuit_optimize(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds);
static PyObject *qop_circuit_get_gates(QopCircuitObject *self);
static PyObject *qop_circuit_run(QopCircuitObject *self, PyObject *args,
                                 PyObject *kwds);

static PyObject *qop_gate_reparameterize(QopGateObject *self, PyObject *args,
                                         PyObject *kwds);
static PyObject *qop_gate_get_parameters(QopGateObject *self);
static PyObject *qop_gate_get_matrix(QopGateObject *self);

static void qop_write_parameter(unsigned int flat_index, double param,
                                void *context);
static void qop_write_energy(_Complex double energy, void *context);

static PyMemberDef qop_circuit_obj_members[] = {
    {"qubits", T_INT, offsetof(QopCircuitObject, qubit_count), READONLY,
     "Number of qubits (lines) in the circuit."},
    {NULL}};

static PyMethodDef qop_circuit_obj_methods[] = {
    {"add_gate", (PyCFunction)qop_circuit_add_gate,
     METH_VARARGS | METH_KEYWORDS,
     "Add a previously created gate to this circuit."},
    {"get_gates", (PyCFunction)qop_circuit_get_gates, METH_NOARGS,
     "Returns a list of all gates added to this circuit."},
    {"optimize", (PyCFunction)qop_circuit_optimize,
     METH_VARARGS | METH_KEYWORDS, "Optimize all given gates"},
    {"run", (PyCFunction)qop_circuit_run, METH_VARARGS | METH_KEYWORDS,
     "Run a state through the circuit."},
    {NULL}};

static PyMethodDef qop_gate_obj_methods[] = {
    {"reparameterize", (PyCFunction)qop_gate_reparameterize,
     METH_VARARGS | METH_KEYWORDS,
     "Set new parameters on a reparameterized gate."},
    {"get_parameters", (PyCFunction)qop_gate_get_parameters, METH_NOARGS,
     "Get a tuple of the gate's current parameters."},
    {"get_matrix", (PyCFunction)qop_gate_get_matrix, METH_NOARGS,
     "Get the 2x2 matrix representation of this gate."},
    {NULL}};

static PyTypeObject QopCircuitType = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "qop.Circuit",
    .tp_doc = "A quantum circuit.",
    .tp_basicsize = sizeof(QopCircuitObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = qop_create_circuit,
    .tp_dealloc = (destructor)qop_circuit_obj_dealloc,
    .tp_members = qop_circuit_obj_members,
    .tp_methods = qop_circuit_obj_methods,
};

static PyTypeObject QopGateType = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_name = "qop.Gate",
    .tp_doc = "A quantum gate.",
    .tp_basicsize = sizeof(QopGateObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = qop_gate_create,
    .tp_dealloc = (destructor)qop_gate_obj_dealloc,
    .tp_methods = qop_gate_obj_methods};

static PyMethodDef qop_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef qopmodule = {
    PyModuleDef_HEAD_INIT, "qop",
    "Python interface for the qop C implementation", -1, qop_methods};

#endif  // QOP_QOPMODULE_H_