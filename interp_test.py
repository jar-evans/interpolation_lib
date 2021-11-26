#!/usr/bin/python3

import interp_lib as interp
import numpy as np 
import matplotlib.pyplot as plt 
from time import time

def f(x):
    return x*x

if __name__ == "__main__":
    X = np.array([-1,0,1])
    Y = f(X).tolist()
    X = X.tolist()
    xs = np.arange(-1, 1, 0.5)
    Z = []

    print(interp.newton())
   
    start = time()
    for x in xs:
        Z.append(interp.lagrange(X,Y,x+0.8))
    end = time()

    print(f"Runtime: {end - start}[s]")

    plt.plot(X,Y)
    plt.plot(xs,Z)
    plt.show()