#include <stdio.h> 
#include <stdlib.h>
#include <Python.h>

static PyObject* py_lagrange(PyObject* self, PyObject* args) {

    PyObject* X, *Y;
    double x;
    int i,j,k;
    double y = 0;

    if(!PyArg_ParseTuple(args, "O|O|d", &X, &Y, &x)) {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y)) {
        return NULL;
    }

    int len = PyObject_Length(X);

    double c_X[len];
    double c_Y[len];

    for (k = 0; k < len; k++) {
        c_X[k] = (double)PyLong_AsLong(PyList_GetItem(X, k));
        c_Y[k] = (double)PyLong_AsLong(PyList_GetItem(Y, k));
    }

    double P[len];
    double N[len];
    double D[len];
    
    for (i=0;i<len;i++) {
        
        N[i] = 1.0;
        D[i] = 1.0;

        for (j=0;j<len-1;j++) {
            if (i != j) {
                N[i] *= x - c_X[j];
                D[i] *= c_X[i] - c_X[j];
            }
        }

        P[i] = c_Y[i]*N[i]/D[i];
        y += P[i];   
    }
    return Py_BuildValue("d", y);

}

static PyObject* py_newton(PyObject* self, PyObject* args) {

    PyObject* X, *Y;
    double x;
    int i,j,k;
    double y = 0;

    if(!PyArg_ParseTuple(args, "O|O|d", &X, &Y, &x)) {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y)) {
        return NULL;
    }

    int len = PyObject_Length(X);

    double difference_table[len][len];

    double c_X[len];
    double c_Y[len];

    for (k = 0; k < len; k++) {
        c_X[k] = (double)PyLong_AsLong(PyList_GetItem(X, k));
        c_Y[k] = (double)PyLong_AsLong(PyList_GetItem(Y, k));
        difference_table[0][k] = c_Y[k];
    }

    for (i = 1; i < len; i++) {
        for (j = 0; j < len-i; j++) {
            difference_table[i][j] = (difference_table[i-1][j+1] - difference_table[i-1][j])/(c_X[i+j] - c_X[j]);
        }     
    }

    for (k = len-1; k >= 0; k--) {
        y *= (x - c_X[k]);
        y += difference_table[k][0];
    }

    return Py_BuildValue("d", y);

}

static PyMethodDef interp_methods[] = {{"lagrange", py_lagrange, METH_VARARGS, NULL}, 
                                       {"newton", py_newton, METH_VARARGS, NULL},
                                       {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduleDef = {PyModuleDef_HEAD_INIT, "interp_lib", NULL, -1, interp_methods};

PyMODINIT_FUNC PyInit_interp_lib(void)
{
    PyObject *module = PyModule_Create(&moduleDef);
    return module;
}