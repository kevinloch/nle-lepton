#include <stdio.h>
#include "nle-lepton.h"

void printUsage() {

printf("nle-lepton version %s\n", NLE_VERSION);
printf("\n\
NAME\n\
     nle-lepton -- Find polynomial-like non-linear equations that solve for the three charged lepton masses, or any three solution masses\n\
\n\
SYNOPSIS\n\
     nle-lepton [-c filename] [-h] [-s seed]\n\
 \n\
OPTIONS:\n\
     -c filename\n\
          Set configuration file name (default: nle-lepton.cfg)\n\
\n\
     -h\n\
          Show help\n\
\n\
     -s seed\n\
          Set external random seed.   This is combined with internal clock-based seeds to initialize srand48();\n\
\n\
DESCRIPTON\n\
 nle-lepton is a tool to search for polynomial-like non-linear equations (nle's) that generate the three charged lepton masses or any three \"solution masses\"\n\
 (sm1, sm2, sm3) specified by the user.  It does this by generating a random minimal formula structure (expnents and masses), solving it for real coefficients,\n\
 and then attempting to factor the coefficients.  When the factors are a close enough match an exact formula is contstructed from the implied factors and solved\n\
 for the solution masses and/or other variables.  The variables are then compared to their experimental uncertainties and the result is output if they match.\n\
 This multi-step process (solve for coefficients, factor, solve for masses) allows for billions of coefficient factors and trillions of coefficient, exponent\n\
 and mass ratio  combinations to be tested as efficiently as possible.\n\
\n\
 The nle's are designed to have three positive real roots corresponding to the three solution masses. This concept might explain why there are exactly three\n\
 generations of matter and an empirical formula might provide insight into the fermion mass generation process.  The solutions masses and their uncertainties\n\
 are user-configurable so they can be set to rest masses, running-masses or any desired values of interest.  The configuration file (by default nle-lepton.cfg)\n\
 also has extensive options to configure operating parameters and potential coefficient factors. As an open-source project the code can be modified by the user\n\
 and is intended to be a platform for quickly testing new formula types.\n\
\n\
 See README.txt in the source tree for the full description of operation and output format:\n\
\n\
   https://raw.githubusercontent.com/kevinloch/nle-lepton/master/README.txt\n\
\n\
 See nle-lepton-4.2.pdf for a scientific description operation and purpose:\n\
\n\
   https://raw.githubusercontent.com/kevinloch/nle-lepton/master/nle-lepton-4.2.pdf\n\
   http://kevinloch.com/papers/nle-lepton-4.2.pdf\n\
 \n");
}
