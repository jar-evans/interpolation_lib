#!/usr/bin/python3

import interp_lib as interp
import numpy as np 
import matplotlib.pyplot as plt 
from time import time

if __name__ == "__main__":

    X = [0.4, 1.1, 2.8, 3.4, 4.2, 5.5, 6.8, 7.8, 8.1, 9.6, 10]
    Y = [i*i*i for i in X]
    x = interp.chebyshev([0, 10], 10)

    y_lagrange = interp.lagrange(X, Y, x)
    y_newton = interp.newton(X, Y, x)
    
    y_lin_spline = interp.lin_spline(X, Y, x)

    plt.plot(x, y_lagrange, label='Lagrange', marker='D')
    plt.plot(x, y_newton, label='Newton', marker='x')

    plt.plot(x, y_lin_spline, label='Linear splines', marker='+')

    plt.legend()

    plt.show()
