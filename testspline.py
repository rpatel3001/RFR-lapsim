"""Test plotting a spline."""

import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
import warnings


def plot_spline(x, y, lens):
    """Plot a spline."""
    lens = np.cumsum(lens)

    x = np.array(x)
    y = np.array(y)
    plt.plot(x, y, label='poly')

    # fit splines to x=f(u) and y=g(u), treating both as periodic. also note that s=0
    # is needed in order to force the spline fit to pass through all the input points.
    warnings.simplefilter("ignore")
    tck, u = interpolate.splprep([x, y], u=lens, s=0, per=True)
    warnings.simplefilter("default")
    # evaluate the spline fits for 1000 evenly spaced distance values
    for i in [0, 12.5, 25, 50, 100, 250, 500, 750, 1000]:
        xi, yi = interpolate.splev(i, tck)
        # plot the result
        plt.scatter(xi, yi, s=10, color='red', label='interp')
    plt.axis('equal')
    plt.legend(loc='best')
    plt.show()
