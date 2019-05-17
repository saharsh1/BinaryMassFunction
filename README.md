# BinaryMassFunction

The modified mass function was introduced by Shahaf, Mazeh and Faigler (2017, MNRAS). It's main advantage 
is that its distribution it mimics the shape of the mass-ratio distribution for samples of single-lined spectroscopic binaries.

The defined class, BinaryMassFunction, contains the following methods:
  1) calc_y - calculates the modified mass function given the binary
              orbital parameters.
  2) calc_q_minimum - calculates the minimal mass-ratio possible for a
              given value of reduced mass function.
  3) calc_S - calculates the modified mass function for a list of given
              list of reduced mass function values.
  A complete description appears in the header of each method.
