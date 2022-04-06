#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include <math.h>

double chebyshev(double a, double b, int i, int len)
{

    double arg = (M_PI / len) * (i + 0.5);

    return (a + b) / 2 + (b - a) * cos(arg) / 2;
}

double lagrange(double X[], double Y[], double x, int len)
{

    int i, j;
    double y = 0;

    double P[len];
    double N[len];
    double D[len];

    for (i = 0; i < len; i++)
    {

        N[i] = 1.0;
        D[i] = 1.0;

        for (j = 0; j < len - 1; j++)
        {
            if (i != j)
            {
                N[i] *= x - X[j];
                D[i] *= X[i] - X[j];
            }
        }

        P[i] = Y[i] * N[i] / D[i];
        y += P[i];
    }
    return y;
}

double newton(double X[], double Y[], double x, int len)
{

    int i, j, k;

    double y = 0;

    double v[len][len];

    for (k = 0; k < len; k++)
    {
        v[0][k] = Y[k];
    }

    for (i = 1; i < len; i++)
    {
        for (j = 0; j < len - i; j++)
        {
            v[i][j] = (v[i - 1][j + 1] - v[i - 1][j]) / (X[i + j] - X[j]);
        }
    }

    for (k = len - 1; k >= 0; k--)
    {
        y *= (x - X[k]);
        y += v[k][0];
    }

    return y;
}

double lin_spline(double X[], double Y[], double x, int len)
{

    int i = 0;
    int j;
    double y = 0;

    while (X[i] <= x)
    {
        i++;
    }

    j = i;

    float dy = Y[j] - Y[j - 1];
    float dx = X[j] - X[j - 1];

    y = Y[j - 1] + (x - X[j - 1]) * dy / dx;
    return y;
}

static PyObject *py_lagrange(PyObject *self, PyObject *args)
{

    PyObject *X, *Y, *x;
    int i;

    if (!PyArg_ParseTuple(args, "OOO", &X, &Y, &x))
    {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y))
    {
        return NULL;
    }

    int len = PyObject_Length(X);
    int len_samples = PyObject_Length(x);

    double c_X[len];
    double c_Y[len];
    double a;
    PyObject *y = PyList_New(len_samples);

    for (i = 0; i < len; i++)
    {
        c_X[i] = (double)PyFloat_AsDouble(PyList_GetItem(X, i));
        c_Y[i] = (double)PyFloat_AsDouble(PyList_GetItem(Y, i));
    }

    for (i = 0; i < len_samples; i++)
    {
        a = (double)PyFloat_AsDouble(PyList_GetItem(x, i));
        PyList_SetItem(y, i, Py_BuildValue("f", lagrange(c_X, c_Y, a, len)));
    }

    return y;
}

static PyObject *py_newton(PyObject *self, PyObject *args)
{

    PyObject *X, *Y, *x;
    int i;

    if (!PyArg_ParseTuple(args, "OOO", &X, &Y, &x))
    {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y))
    {
        return NULL;
    }

    int len = PyObject_Length(X);
    int len_samples = PyObject_Length(x);

    double c_X[len];
    double c_Y[len];
    double a;
    PyObject *y = PyList_New(len_samples);

    for (i = 0; i < len; i++)
    {
        c_X[i] = (double)PyFloat_AsDouble(PyList_GetItem(X, i));
        c_Y[i] = (double)PyFloat_AsDouble(PyList_GetItem(Y, i));
    }

    for (i = 0; i < len_samples; i++)
    {
        a = (double)PyFloat_AsDouble(PyList_GetItem(x, i));
        PyList_SetItem(y, i, Py_BuildValue("f", newton(c_X, c_Y, a, len)));
    }

    return y;
}

static PyObject *py_lin_spline(PyObject *self, PyObject *args)
{

    PyObject *X, *Y, *x;
    int i;

    if (!PyArg_ParseTuple(args, "OOO", &X, &Y, &x))
    {
        return NULL;
    }

    if (PyObject_Length(X) != PyObject_Length(Y))
    {
        return NULL;
    }

    int len = PyObject_Length(X);
    int len_samples = PyObject_Length(x);

    double c_X[len];
    double c_Y[len];
    double a;
    PyObject *y = PyList_New(len_samples);

    for (i = 0; i < len; i++)
    {
        c_X[i] = (double)PyFloat_AsDouble(PyList_GetItem(X, i));
        c_Y[i] = (double)PyFloat_AsDouble(PyList_GetItem(Y, i));
    }

    for (i = 0; i < len_samples; i++)
    {
        a = (double)PyFloat_AsDouble(PyList_GetItem(x, i));
        PyList_SetItem(y, i, Py_BuildValue("f", lin_spline(c_X, c_Y, a, len)));
    }

    return y;
}

static PyObject *py_chebyshev(PyObject *self, PyObject *args)
{

    PyObject *X;
    int i, L;

    if (!PyArg_ParseTuple(args, "Oi", &X, &L))
    {
        return NULL;
    }

    PyObject *y = PyList_New(L);

    double A = (double)PyFloat_AsDouble(PyList_GetItem(X, 0));
    double B = (double)PyFloat_AsDouble(PyList_GetItem(X, 1));

    for (i = 0; i < L; i++)
    {
        PyList_SetItem(y, i, Py_BuildValue("f", chebyshev(A, B, i, L)));
    }

    return y;
}

static PyMethodDef interp_methods[] = {{"lagrange", py_lagrange, METH_VARARGS, NULL},
                                       {"newton", py_newton, METH_VARARGS, NULL},
                                       {"lin_spline", py_lin_spline, METH_VARARGS, NULL},
                                       {"chebyshev", py_chebyshev, METH_VARARGS, NULL},
                                       {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduleDef = {PyModuleDef_HEAD_INIT, "interp_lib", NULL, -1, interp_methods};

PyMODINIT_FUNC PyInit_interp_lib(void)
{
    PyObject *module = PyModule_Create(&moduleDef);
    return module;
}