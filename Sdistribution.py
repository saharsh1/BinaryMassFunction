#               ---------------------------------------
#               The modified mass function distribution
#               ---------------------------------------
#
# This code utilizes the modified mass function, S, introduced
# by Shahaf, Mazeh and Faigler (2017, MNRAS). The main advantage of the
# modified mass function is that it mimics the shape of the mass-ratio
# distribution for samples of single-lined spectroscopic binaries.
# This routine calculates a model for the distribution of the modified
# mass function, given some input model for the mass-ratio distribution.
#
# The purpose of this routine is to demonstate the distribution of the
# modified mass function. Detailed calculation may require some modification
# to the integration tolerance and grid density.
#
#   Example:
#   import Sdistribution as SD
#   MRD = lambda x: x
#   fS, S = SD.calc_fS(MRD, PlotFlag=True)


import numpy as np
from scipy import integrate
from ModifiedMassFunction import BinaryMassFunction
import matplotlib.pyplot as plt


def calc_fS(MRD, PlotFlag=False):
    '''
    This function calculates the distribution of the modified mass function, S,
    for a given callable, MRD, that represents the mass-ratio distribution in
    the range 0-1. The function calculates the distribution of S over a pre
    determined grid.

    Input: MRD - callable, the mass ratio distirbution (0<q<1)
           PlotFlag (optional) - if true a plot with the given mass-ratio.
                 distribution and the calculated S distribution is generated.

    Output: fS - array with the calculated probability density
            S - the corresponding grid
    '''
    yGridVals = np.concatenate((10**np.linspace(-10, -1, num=500),
                                np.linspace(0.1001, 0.24999, num=250)),
                               axis=0)

    Grid = BinaryMassFunction(y=yGridVals).calc_q_minimum()
    Grid = BinaryMassFunction(y=yGridVals).calc_S()
    fS = np.zeros(Grid.y.shape)*np.NaN

    def kernelK(q, y):
        # This is the kernel of the integral used to calculate the
        # distirbution of the modified mass function. See equation (9)
        # in Sahahaf et al. (2017) for details.
        Knumer = (1+q)**(4/3)
        Kdenomer = 3*y**(1/3)*q*np.sqrt(q**2-y**(2/3)*Knumer)
        k = Knumer/Kdenomer
        return k

    def fqKernel(q, y, MRD):
        fqK = MRD(q)*kernelK(q, y)
        return fqK

    for I, val in enumerate(Grid.y):
            lowerBound = Grid.qmin[I]
            fy_of_fq = integrate.quadrature(fqKernel, lowerBound, 1,
                                            args=(val, MRD), tol=1e-04,
                                            rtol=1e-04, maxiter=5000)
            fy_of_1 = integrate.quadrature(kernelK, lowerBound, 1,
                                           args=(val), tol=1e-04,
                                           rtol=1e-04, maxiter=5000)
            fS[I] = fy_of_fq[0] / fy_of_1[0]

    if PlotFlag:
        q = np.linspace(0, 1, num=250)
        plt.plot(q, MRD(q), '--k')
        plt.plot(Grid.S, fS, 'r')
        plt.legend(['fq', 'fS'], loc='upper left')
        plt.grid()
        plt.xlabel('S / q')
        plt.ylabel('probability density')
        plt.show()

    return fS, Grid.S
