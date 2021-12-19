#include <stdio.h> 
#include <stdlib.h>
#include <Python.h>

double lagrange(double X[], double Y[], double x, int len) {

    int i,j;
    double y = 0;

    double P[len];
    double N[len];
    double D[len];
    
    for (i=0;i<len;i++) {
        
        N[i] = 1.0;
        D[i] = 1.0;

        for (j=0;j<len-1;j++) {
            if (i != j) {
                N[i] *= x - X[j];
                D[i] *= X[i] - X[j];
            }
        }

        P[i] = Y[i]*N[i]/D[i];
        y += P[i];   
    }
    return y;
}

double newton(double X[], double Y[], double x, int len) {

    int i, j, k;

    double y = 0;

    double difference_table[len][len];

    for (k = 0; k < len; k++) {
        difference_table[0][k] = Y[k];
    }

    for (i = 1; i < len; i++) {
        for (j = 0; j < len-i; j++) {
            difference_table[i][j] = (difference_table[i-1][j+1] - difference_table[i-1][j])/(X[i+j] - X[j]);
        }     
    }

    for (k = len-1; k >= 0; k--) {
        y *= (x - X[k]);
        y += difference_table[k][0];
    }

    return y;

}

double monomial(double X[], double Y[], double x, int len) {

    int i, j, k;
    double y = 0;
    //Initialise Vandermonde
    double V[len][len];
    for (i = 0; i < len; i++) {
        for (j = 0; j < len; j++) {
            V[i][j] = 1;
            for (k = 0; k < j; k++) {
                V[i][j] *= X[i];
            }
        }
    }

    double matrix[len][len];
	for (i = 0; i < len; i++) {
		for (j = 0; j < 2 * len; j++) {

            if (j < len) {
                matrix[i][j] = V[i][j];
            }

			if (j == (i + len)) {
				matrix[i][j] = 1;
                continue;
            }

            matrix[i][j] = 0;
        }
    }

    double temp;

	for (i = len - 1; i > 0; i--) {

		if (matrix[i - 1][0] < matrix[i][0])
            for (int j = 0; j < 2 * len; j++) {
                temp = matrix[i][j];
                matrix[i][j] = matrix[i - 1][j];
                matrix[i - 1][j] = temp;
            }
	}


	for (i = 0; i < len; i++) {
        if(matrix[i][i] == 0.0)
        {
            printf("No inverse!\n");
            exit(0);
        }
		for (j = 0; j < len; j++) {

			if (j != i) {
				temp = matrix[j][i] / matrix[i][i];
				for (k = 0; k < 2 * len; k++) {
					matrix[j][k] -= matrix[i][k] * temp;
				}
			}
		}
	}

	for (i = 0; i < len; i++) {
		temp = matrix[i][i];
		for (j = 0; j < 2 * len; j++) {
			matrix[i][j] = matrix[i][j] / temp;
		}
	}

    double inv[len][len];

	for (i = 0; i < len; i++) {
		for (j = len; j < 2 * len; j++) {
			inv[i][j-len] = matrix[i][j];
		}
	}

    double coeffs[len];

    for (i = 0; i < len; i++) {
        coeffs[i] = 0;
        for (j = 0; j<len; j++) {
            coeffs[i] += inv[i][j]*Y[j];
        }
    }
    y = 0;
    float w = 0;
    for (i = len-1; i >=0; i--) {
        w = coeffs[i];
        for (j = 0; j<i; j++) {
            w *= x;
        }
        y += w;
    }

    return y;
}

double l_spline(double X[], double Y[], double x, int len) {

    int i;

    double y = 0;

    for (i = 0; i < len; i++) {
        if (X[i] > x) {
            break;
        }
    }

    double a = X[i-1];

    double dy = Y[i] - Y[i-1];
    double dx = X[i] - X[i-1];

    y = Y[i-1] + (dy/dx)*(x - a);
    return y;

}

static PyObject* py_lagrange(PyObject* self, PyObject* args) {

    PyObject* X, *Y, *x;
    int i;

    if(!PyArg_ParseTuple(args, "O|O|O", &X, &Y, &x)) {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y)) {
        return NULL;
    }

    int len = PyObject_Length(X);
    int len_samples = PyObject_Length(x);

    double c_X[len];
    double c_Y[len];
    double c_x[len_samples];
    PyObject* y = PyList_New(len_samples);

    for (i = 0; i < len; i++) {
        c_X[i] = (double)PyLong_AsLong(PyList_GetItem(X, i));
        c_Y[i] = (double)PyLong_AsLong(PyList_GetItem(Y, i));
    }

    for (i = 0; i < len_samples; i++) {
        c_x[i] = (double)PyLong_AsLong(PyList_GetItem(x, i));
    }

    for (i = 0; i <len_samples; i++) {
        PyList_SetItem(y, i, Py_BuildValue("f",lagrange(c_X, c_Y, c_x[i], len)));
    }

    return y;



}

static PyObject* py_newton(PyObject* self, PyObject* args) {

    PyObject* X, *Y, *x;
    int i;

    if(!PyArg_ParseTuple(args, "O|O|O", &X, &Y, &x)) {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y)) {
        return NULL;
    }

    int len = PyObject_Length(X);
    int len_samples = PyObject_Length(x);

    double c_X[len];
    double c_Y[len];
    double c_x[len_samples];
    PyObject* y = PyList_New(len_samples);

    for (i = 0; i < len; i++) {
        c_X[i] = (double)PyLong_AsLong(PyList_GetItem(X, i));
        c_Y[i] = (double)PyLong_AsLong(PyList_GetItem(Y, i));
    }

    for (i = 0; i < len_samples; i++) {
        c_x[i] = (double)PyLong_AsLong(PyList_GetItem(x, i));
    }

    for (i = 0; i <len_samples; i++) {
        PyList_SetItem(y, i, Py_BuildValue("f",newton(c_X, c_Y, c_x[i], len)));
    }

    return y;

}

static PyObject* py_monomial(PyObject* self, PyObject* args) {

    PyObject* X, *Y, *x;
    int i;

    if(!PyArg_ParseTuple(args, "O|O|O", &X, &Y, &x)) {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y)) {
        return NULL;
    }

    int len = PyObject_Length(X);
    int len_samples = PyObject_Length(x);

    double c_X[len];
    double c_Y[len];
    double c_x[len_samples];
    PyObject* y = PyList_New(len_samples);

    for (i = 0; i < len; i++) {
        c_X[i] = (double)PyLong_AsLong(PyList_GetItem(X, i));
        c_Y[i] = (double)PyLong_AsLong(PyList_GetItem(Y, i));
    }
    


    for (i = 0; i <len_samples; i++) {
        PyList_SetItem(y, i, Py_BuildValue("f",monomial(c_X, c_Y, c_x[i], len)));
    }

    return y;
 
}

static PyObject* py_l_spline(PyObject* self, PyObject* args) {

    PyObject* X, *Y, *x;
    int i;

    if(!PyArg_ParseTuple(args, "O|O|O", &X, &Y, &x)) {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y)) {
        return NULL;
    }

    int len = PyObject_Length(X);
    int len_samples = PyObject_Length(x);

    double c_X[len];
    double c_Y[len];
    double c_x[len_samples];
    PyObject* y = PyList_New(len_samples);

    for (i = 0; i < len; i++) {
        c_X[i] = (double)PyLong_AsLong(PyList_GetItem(X, i));
        c_Y[i] = (double)PyLong_AsLong(PyList_GetItem(Y, i));
    }

    for (i = 0; i < len_samples; i++) {
        c_x[i] = (double)PyLong_AsLong(PyList_GetItem(x, i));
    }

    for (i = 0; i <len_samples; i++) {
        PyList_SetItem(y, i, Py_BuildValue("f",l_spline(c_X, c_Y, c_x[i], len)));
    }

    return y;



}

static PyMethodDef interp_methods[] = {{"lagrange", py_lagrange, METH_VARARGS, NULL}, 
                                       {"newton", py_newton, METH_VARARGS, NULL},
                                       {"monomial", py_monomial, METH_VARARGS, NULL},
                                       {"l_spline", py_l_spline, METH_VARARGS, NULL},
                                       {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduleDef = {PyModuleDef_HEAD_INIT, "interp_lib", NULL, -1, interp_methods};

PyMODINIT_FUNC PyInit_interp_lib(void)
{
    PyObject *module = PyModule_Create(&moduleDef);
    return module;
}