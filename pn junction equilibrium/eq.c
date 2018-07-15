/* eq.c */

/* Header file */
#include "eq.h"





// <function equilibrium__Poisson_1D__BC_VonNeumann__direct_iterative__manual( TODO
//		PyObject* PyList_x, PyObject* PyList_e,
//		PyObject* PyList_NC, PyObject* PyList_NV,
//		PyObject* PyList_EC, PyObject* PyList_EV,
//		PyObject* PyList_D,
//		PyObject* PyFloat_T,
//		PyObject* PyLong_MAXCNT, PyObject* PyFloat_MINTOL)
//
// 	@argument <PyObject* PyList_x> : longitudinal mesh
// 	@argument <PyObject* PyList_e> : absolute permittivity [Farad/centimeter]
// 	@argument <PyObject* PyList_NC> : conductive band density of states vs x [centimeter^(-3)]
// 	@argument <PyObject* PyList_NV> : conductive band density of states vs x [centimeter^(-3)]
// 	@argument <PyObject* PyList_EC> : energy for conduction band vs x [centimeter^(-3)]
// 	@argument <PyObject* PyList_EV> : energy for valence band vs x [centimeter^(-3)]
// 	@argument <PyObject* PyList_D> : doping concentration vs x (ND⁻ - NA⁺) [centimeter^(-3)]
// 	@argument <PyObject* PyFloat_T> : temperature [Kelvin]
// 	@argument <PyObject* PyLong_MAXCNT> : maximum number of iterations [<integer>]
// 	@argument <PyObject* PyFloat_MINTOL> : minimum tolerance [<float>]
//
//	@returns <PyObject* PyList_v> : electric potential vs x [Volt]
//
// 	@description : solves the Maxwell's first equation (Poisson equation)
//		(ε·Ф')' = n(Ф,x) - p(Ф,x) - D(x) with Von Neumann boundary
//		conditions Ф'(0) = 0 Ф'(L) = 0 by direct substitution method
//
//	@name : equilibrium__Poisson_1D__BC_VonNeumann__direct_iterative__manual
//		    |            |       |   |   |          |                 |
//		    |            |       |   |   |          |                 manual mode
//		    |            |       |   |   |          Finite differences decoupled
//		    |            |       |   |   Von Neummann
//		    |            |       |   Boundary conditions...
//		    |            |       1D
//		    |            Poisson's equation (Maxwell's first equation)
//		    Equilibrium (np = n_i²)
//
// 	@author : Daniel Ríos Linares
//
//	@version : 0.1.0 July 10, 2018
//
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
	PyObject* PyLong_MAXCNT, PyObject* PyFloat_MINTOL) {

	/* Check list size */
	int N = PyList_Size(PyList_x);

	/* Temporal variables */
	int j;

	/* Declare array */
	double x[N+2]; // 1D mesh [nanometer]
	double e[N+2]; // absolute permittivity vs x [farad/centimeter]
	double NC[N+2]; // intrinsic carrier concentration (n_i) vs x [centimeter^(-3)]
	double NV[N+2]; // intrinsic carrier concentration (n_i) vs x [centimeter^(-3)]
	double EC[N+2]; // intrinsic carrier concentration (n_i) vs x [centimeter^(-3)]
	double EV[N+2]; // intrinsic carrier concentration (n_i) vs x [centimeter^(-3)]
	double ND[N+2]; // manufacturing extrinsic donor concentration [centimeter^(-3)]
	double NA[N+2]; // manufacturing extrinsic acceptor concentration [centimeter^(-3)]
	double gD[N+2]; // extrinsic donor ionization carrier coefficient [centimeter^(-3)]
	double gA[N+2]; // extrinsic acceptor ionization carrier coefficient [centimeter^(-3)]
	double ED[N+2]; // extrinsic donor ionization energy band [centimeter^(-3)]
	double EA[N+2]; // extrinsic acceptor ionization energy band [centimeter^(-3)]
	double D[N+2]; // doping concentration (N_D - N_A) vs x [centimeter^(-3)]

	/* Stop conditions */
	long MAXCNT = PyLong_AsLong(PyLong_MAXCNT); // Maximum number of
	double MINTOL = PyFloat_AsDouble(PyFloat_MINTOL); // Minimum tolerance before exit loop

	/* Get floats */
	static double k = 1.38064852e-023;//8.6173303e-005; // Boltmann constant [electronvolt/Kelvin]
	static double q = 1.6021766208e-019; // Elementary charge constant [Coulomb]
	double T = PyFloat_AsDouble(PyFloat_T); // Temperature [kelvin]
	double VT = k * T / q;

	/* Python list<float> to C double array */
	for (int i = 0 ; i < N ; i++) {
		x[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_x, i));
		e[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_e, i));
		NC[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_NC, i));
		NV[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_NV, i));
		EC[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_EC, i));
		EV[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_EV, i));
		gD[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_gD, i));
		gA[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_gA, i));
		ED[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_ED, i));
		EA[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_EA, i));
		ND[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_ND, i)) / (1 + gD[i+1] * exp(-ED[i+1] / VT));
		NA[i+1] = PyFloat_AsDouble(PyList_GetItem(PyList_NA, i)) / (1 + gA[i+1] * exp(+EA[i+1] / VT));
		D[i+1] = ND[i+1] - NA[i+1];
	}

	/* Von Neumann extension */
	x[0] = x[1] - (x[2] - x[1]); x[N+1] = x[N] + (x[N] - x[N-1]);
	e[0] = e[1]; e[N+1] = e[N];
	NC[0] = NC[1]; NC[N+1] = NC[N];
	NV[0] = NV[1]; NV[N+1] = NV[N];
	EC[0] = EC[1]; EC[N+1] = EC[N];
	EV[0] = EV[1]; EV[N+1] = EV[N];
	D[0] = D[1]; D[N+1] = D[N];

	////////////////////////////////////////////////////////////////////////////

	/** Solve the Poisson equation (ε·Ф')' = n(Ф,x) - p(Ф,x) - D(x) **/

	/* Fractions : only for convenience and reduce the computational cost */
	double v_pond_f_p1[N+2];
	double v_pond_f_m1[N+2];
	double p_m_n_pond_f_p0[N+2];
	for (int j = 0 ; j < N+2 ; j++) v_pond_f_p1[j] = ((e[j+1] + e[j]) / (x[j+1] - x[j])) /
		((e[j+1] + e[j]) / (x[j+1] - x[j]) + (e[j] + e[j-1]) / (x[j] - x[j-1]));
	for (int j = 0 ; j < N+2 ; j++) v_pond_f_m1[j] = ((e[j] + e[j-1]) / (x[j] - x[j-1])) /
		((e[j+1] + e[j]) / (x[j+1] - x[j]) + (e[j] + e[j-1]) / (x[j] - x[j-1]));
	for (int j = 0 ; j < N+2 ; j++) p_m_n_pond_f_p0[j] = q * (x[j+1] - x[j-1]) /
		((e[j+1] + e[j]) / (x[j+1] - x[j]) + (e[j] + e[j-1]) / (x[j] - x[j-1]));

	/* Declare target array */
	double v[N+2]; for (int j = 0 ; j < N+2 ; j++) v[j] = (EC[j] + EV[j]) / 2; // potential [Volts]

	/* Loop initialization */
	int cnt = 0; // Current iteration
	double tol = MINTOL * 1.01; // Current tolerance (current - last)

	/* Equilibrium */
	do {
		for (int j = 1 ; j < N+1 ; j++) {
			double last_vj = v[j];
			v[j] = v[j+1] * v_pond_f_p1[j] + v[j-1] * v_pond_f_m1[j] - p_m_n_pond_f_p0[j] *
				( NC[j] * exp((v[j] - EC[j]) / VT) - NV[j] * exp((EV[j] - v[j]) / VT) - D[j]);
			if (abs(v[j] - last_vj) > tol) tol = abs(v[j] - last_vj);
		}
		// Von Neumann conditions
		v[0] = v[1]; v[N+1] = v[N];
		// Loop stop number of cycles condition
		cnt += 1;
	} while (cnt <= MAXCNT && tol >= MINTOL);

	/* Store the important segment of the array v (the first and the last were
	used for applying the Von Neumann boundary conditions) */
	double t[N]; for (int j = 0 ; j < N ; j++) t[j] = v[j+1];

	////////////////////////////////////////////////////////////////////////////

	/* C double array to Python list<float> */
	// Electric potential
	PyObject* PyList_out_v = PyList_New(N);
    for (int i = 0; i < N; i++)
    	PyList_SetItem(PyList_out_v, i, PyFloat_FromDouble(t[i]));
	/* End of Python list fill */

    return PyList_out_v;
}













//
