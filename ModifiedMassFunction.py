#                   ------------------------------
#                     The modified mass function
#                   ------------------------------
#
# This code utilizes the modified mass function, S, introduced
# by Shahaf, Mazeh and Faigler (2017, MNRAS). The main advantage of the
# modified mass function is that it mimics the shape of the mass-ratio
# distribution for samples of single-lined spectroscopic binaries.
#
# The defined class, BinaryMassFunction, contains the following methods:
#  1) calc_y - calculates the reduced mass function given the binary
#              orbital parameters.
#  2) calc_q_minimum - calculates the minimal mass-ratio possible for a
#              given value of reduced mass function.
#  3) calc_S - calculates the modified mass function for a list of given
#              list of reduced mass function values.
#   NOTE: A complete description appears in the header of each method.
#
#
#  Examples:
#    1)  Calculate the reduced mass function for three binaries with given P
#        and K, assuming a circular orbit and 1 Msun primary:
#        Y = BinaryMassFunction().calc_y(P=[10.,12.,1.17], K=[0.5, 1.17, 2.])
#        print(Y.y)
#
#    2)  Calculate the q minimum for a given array of reduced mass function:
#        Y = BinaryMassFunction(y=[0.1, 0.001, 0.25]).calc_q_minimum()
#        print(Y.qmin)
#
#    3)  Calculate the modified mass function for an array of y:
#        Y = BinaryMassFunction().calc_y(P=[10.,12.,1.17], K=[0.5, 1.17, 2.])
#        Y.calc_S()
#        print(Y.S)
#
# Dependencies: numpy, scipy.integrate

import numpy as np
from scipy import integrate


class BinaryMassFunction:

    # =============================================================================
    #                             Initialization
    # =============================================================================
    def __init__(self, **kwargs):
        '''
        Input: Input is optional, and needs to be called with its keyword.
               Below appears a list of the possible input variables.
        '''
        if 'y' in kwargs:
            self.y = np.array(kwargs['y'], dtype=float)
        else:
            self.y = []

    # =============================================================================
    #                                 calc_y
    # =============================================================================
    def calc_y(self, P, K, M1=1, e=0):
        '''
        Calculate the reduced mass function from the orbital parameters.
        Input: P (days), K(km/s), e (unitless), M1 (solar masses).
               Inputs must be arrays of the same size, or scalars.
               Default values: M1 = 1 solar mass and e = 0.

        Output: self.y, numpy array with the (unitless) reduced mass function
        The reduced mass function, y,

                               P * K^3
                    (1)  y = ----------- (1-e^2)^1.5  ,
                              2*pi*G*m1

         where P, K, e and m1 are the orbital period, RV semi-amplitude, eccentricity
         and primary mass, respectively. The reduced mass function can be expressed
         in terms of the mass ratio, q=m2/m1 and orbital inclination, i,

                               q^3
                    (2)  y = -------- * sin(i)^3  ,
                             (1+q)^2

        where q and i are the mass-ratio and orbital inclination, respectively.
        This reduced mass function links the observables from equation (1) to
        the unknown quantities in equation (2) for each single-lined system.
        '''
        # In the natural units (years, AU, solar mass)
        G = 4*np.pi**2
        d_to_yr = 0.00273790926
        kms_to_AUyr = 0.210945021

        # Convert the input to natural units
        K = np.multiply(kms_to_AUyr, K, dtype=float)
        P = np.multiply(d_to_yr, P, dtype=float)
        M1 = np.array(M1)
        e = np.array(e)

        # Calculate the reduced mass function
        self.y = np.array(P*(K**3)*(1-e**2)**(3/2)/(2*np.pi*G*M1))

        return self

    # =============================================================================
    #                            calc_q_minimum
    # =============================================================================
    def calc_q_minimum(self):
        '''
        This function calculates the minimal mass ratio possible for the
        given reduced mass function values, by assuming that the inclination
        angle, i, is 90 degrees and solving equation (2).

        Input: No additional input required.
        Output: self.qmin - an array of the same size as self.y, that
                            contains the calculated q minimum.

        Example:
        Y = BinaryMassFunction(y=[0.1, 0.15, 0.2, 0.25]).calc_q_minimum()
        print(Y.qmin)
        '''
        y = self.y
        h = (y/2 + (y**2)/3 + (y**3)/27
             + np.sqrt(3)/18*y*np.sqrt(4*y+27))**(1/3)

        self.qmin = np.array(h + (2*y/3 + (y**2)/9)/h + y/3)
        return self

    # =============================================================================
    #                                   calc_S
    # =============================================================================
    def calc_S(self):
        '''
        The modified mass function, S, is calculated by an integral expression,
        given in Shahaf et al. (2017). It is defined for mass ratios in the
        range of 0-1, i.e., 0 < y < 0.25. The possible range of values for
        the modified mass function is also between 0 to 1.

        Input: No additional input required.
        Output: self.S - an array of the same size as self.y, that
                         contains the calculated modified mass function.
                         If q minimum is larger than 1, S is not defined
                         and NaN is returned.

        Example:
        Y = BinaryMassFunction(y=[0.1, 0.15, 0.2, 0.25]).calc_S()
        print(Y.S)
        '''

        # Initialize:
        if 'qmin' not in self.__dict__:
            self.calc_q_minimum()
        S = np.zeros(self.y.shape)*np.NaN

        # Define the integrand of the modified mass function:
        def integrandS(x, y):
            intS = (1 - y**(2/3)*(1+x)**(4/3)*x**(-2))**(0.5)
            return intS

        # Calculate the modified mass function
        for I, val in enumerate(self.y):
            if val > 0 and val < 0.25:
                lowerBound = self.qmin[I]
                stmp = integrate.quad(integrandS, lowerBound, 1, args=(val))
                S[I] = 1 - stmp[0]
            elif val == 0:
                S[I] = 0
            elif val == 0.25:
                S[I] = 1

        self.S = S
        return self
