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
} QopCircuitObject;

typedef struct QopGateObject {
  PyObject_HEAD Gate gate;
  double *params;
} QopGateObject;

static PyObject *qop_create_circuit(PyTypeObject *type, PyObject *args,
                                    PyObject *kwds);
static PyObject *qop_create_gate(PyTypeObject *type, PyObject *args,
                                 PyObject *kwds);
static PyObject *qop_circuit_add_gate(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds);
static PyObject *qop_get_gates(QopCircuitObject *self);
static PyObject *qop_optimize_circuit(QopCircuitObject *self, PyObject *args,
                                      PyObject *kwds);
static void qop_circuit_obj_dealloc(QopCircuitObject *obj);
static void qop_gate_obj_dealloc(QopGateObject *obj);

static PyMemberDef qop_circuit_obj_members[] = {
    {"qubits", T_INT, offsetof(QopCircuitObject, qubit_count), READONLY,
     "Number of qubits (lines) in the circuit."},
    {NULL}};

static PyMethodDef qop_circuit_obj_methods[] = {
    {"add_gate", (PyCFunction)qop_circuit_add_gate,
     METH_VARARGS | METH_KEYWORDS,
     "Add a previously created gate to this circuit."},
    {"optimize", (PyCFunction)qop_optimize_circuit,
     METH_VARARGS | METH_KEYWORDS, "Optimize all given gates"},
    {NULL}};

static PyTypeObject QopCircuitType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "qop.Circuit",
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
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "qop.Gate",
    .tp_doc = "A quantum gate.",
    .tp_basicsize = sizeof(QopGateObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = qop_create_gate,
    .tp_dealloc = (destructor)qop_gate_obj_dealloc};

static PyMethodDef qop_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef qopmodule = {PyModuleDef_HEAD_INIT, "qop",
                                       NULL, /* module documentation */
                                       -1, qop_methods};

#endif  // QOP_QOPMODULE_H_