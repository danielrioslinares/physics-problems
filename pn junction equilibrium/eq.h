/* eq.h */

/* Standard bool (C99) */
#include <stdbool.h>

/* Python C API */
#include <Python.h>

/* Math */
#include <math.h>

/* eq.c functions */
//PyObject* equilibrium__Poisson_1D__BC_VonNeumann__direct_iterative__manual(PyObject* PyList_x, PyObject* PyList_e, PyObject* PyList_n, PyObject* PyList_D, PyObject* PyFloat_T, PyObject* PyLong_MAXCNT, PyObject* PyFloat_MINTOL);
//PyObject* equilibrium__Poisson_1D__BC_VonNeumann__newton_raphson__manual(PyObject* PyList_x, PyObject* PyList_e, PyObject* PyList_n, PyObject* PyList_D, PyObject* PyFloat_T, PyObject* PyLong_MAXCNT, PyObject* PyFloat_MINTOL);
//PyObject* equilibrium__Poisson_1D__BC_VonNeumann__newton_raphson_direct_iterative__manual(PyObject* PyList_x, PyObject* PyList_e, PyObject* PyList_n, PyObject* PyList_D, PyObject* PyFloat_T, PyObject* PyLong_MAXCNT, PyObject* PyFloat_MINTOL);
PyObject* equilibrium__Poisson_1D__BC_VonNeumann__direct_iterative__manual(
	// Mesh
	PyObject* PyList_x,
	// Semiconductor static features
	PyObject* PyList_e,
	PyObject* PyList_NC, PyObject* PyList_NV,
	PyObject* PyList_EC, PyObject* PyList_EV,
	PyObject* PyList_ND, PyObject* PyList_NA,
	PyObject* PyList_gD, PyObject* PyList_gA,
	PyObject* PyList_ED, PyObject* PyList_EA,
	// Temperature
	PyObject* PyFloat_T,
	// Loop properties
	PyObject* PyLong_MAXCNT, PyObject* PyFloat_MINTOL);
