#!/usr/bin/python3

import interp_lib as interp
import numpy as np 
import matplotlib.pyplot as plt 
from time import time

if __name__ == "__main__":
    X = [0, 1]
    Y = [0, 1]
    x = [0, 1, 2, 3]
    y = []
    start = time()
    y = interp.newton(X, Y, x)
    end = time()

    print(end - start)

    plt.scatter(X, Y)
    plt.plot(x, y)
    plt.show()