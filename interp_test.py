#!/usr/bin/python3

import interp_lib as interp
import numpy as np 
import matplotlib.pyplot as plt 
from time import time

if __name__ == "__main__":
    X = [0, 2, 4, 6, 8]
    Y = [0, 4, 16, 36, 64]
    x = [0, 1, 3, 5, 7, 9]

    plt.scatter(X, Y)
    # plt.plot(x, interp.newton(X, Y, x))
    # plt.plot(x, interp.lagrange(X, Y, x))
    # plt.plot(x, interp.monomial(X, Y, x))
    plt.plot(x, interp.l_spline(X, Y, x))
    plt.show()