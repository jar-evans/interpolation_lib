#!/usr/bin/python3

import interp_lib as interp
import numpy as np 
import matplotlib.pyplot as plt 
from time import time

if __name__ == "__main__":
    X = [0, np.pi/2, np.pi]
    Y = [0, 1, 0]
    x = np.arange(0, np.pi+0.1, 0.1).tolist()
    y = []
    for sample in x:
        y.append(interp.monomial(X, Y, sample))

    plt.scatter(X, Y)
    plt.plot(x, y)
    plt.show()