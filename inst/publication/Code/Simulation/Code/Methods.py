#!/usr/bin/env python3
"""Helper functions for honeybee foraging simulation model"""
__appname__  = "Methods.py"
__author__   = "Joseph Palmer <joseph.palmer18@imperial.ac.uk>"
__version__  = "0.0.2"
__date__     = "02-2020"

import numpy as np
import matplotlib.pyplot as plt

def getdist(x2, y2):
    # return distance of resource to hive in 2D plane
    x1 = 0
    y1 = 0
    return np.sqrt(((x1 - x2)**2 + (y1 - y2)**2))

def get_ccdf(data):
    # Get CCDF for the data
    sorted_data = np.sort(data)[::-1]
    prob = np.linspace(0, 1, len(data))
    return np.column_stack((sorted_data, prob))

def plotccdf(data, title, log=True):
    # Plot cdf data
    cdf_dat = get_ccdf(data)
    if log:
        plt.plot(cdf_dat[:,0][1:], np.log10(cdf_dat[:,-1][1:]), "o", markersize=2)
    else:
        plt.plot(cdf_dat[:,0][1:], cdf_dat[:,-1][1:], "o", markersize=2)
    plt.ylabel("ln Probability")
    plt.xlabel("Foraging distance (km)")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.title(title)
    plt.show()
    return None


def qual2dance(q, x, alpha):
    res = q/(1+alpha*x)-1
    res[res<=0] = 0
    return res
