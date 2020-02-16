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
 nle-lepton is a user-configurable tool to search for polynomial-like non-linear equations (NLE's) that generate the three charged lepton masses or any three\n\
 \"solution masses\" (sm1, sm2, sm3) specified by the user.  The NLE's are designed to have three positive real roots corresponding to the three solution masses.\n\
 This concept might explain why there are exaclty three generations of matter and an empirical formula might provide insight into the fermion mass generation\n\
 process.  The solution masses and their uncertainties are user-configurable so they can be set to rest masses, running-masses or any desired values of\n\
 interest.  As an open source project the source code can be modified to add additional equation types or coefficient factors to be tested. nle-lepton uses a\n\
 configuration file (by default nle-lepton.cfg).   Edit that file to change search parameters and experimental reference values.\n\
\n\
 nle-lepton currently supports \"3-term\" mode, which uses three mass terms plus a constant term of -1, with a different exponent on each mass term. Exponents are\n\
 currently limited to fractional values between -1 and 1, excluding zero.  The use of fractional exponents allows the charged lepton mass spectrum to be\n\
 replicated with relatively small coefficients on each term, often in the range of 0.1 to 10.  With traditional integer exponent polynomials some coefficients\n\
 would be much larger or smaller and vary dramatically between terms. To ensure each term is dimensionless, the solution masses are combined with a reference\n\
 mass in  side each term as a \"solution mass ratio\".  The reference masses currently supported are the Planck mass \"mp\", Higgs vacuum expectation value \"v\", Z\n\
 boson mass \"mz\", W  boson mass \"mw\", and Higgs boson mass \"mH0\".  \n\
\n\
 Since fractional exponent terms might represent a relative length scale derived from a higher dimensional manifold, n-ball volume (nbv) and n-sphere surface\n\
 area (nss) geometric constant factors are included.  For example the geometric factor for the common 3-ball volume nbv[3]=4pi/3.  This appears in the formula\n\
 for the radius r of a 3-ball with volume V as r=(V/nbv[3])^(1/3).\n\
\n\
 Processing is broken into two phases.\n\
\n\
 Phase 1 - solve for coefficients using the known masses and then attempt to factor the coefficients\n\
   Randomly selected exponents put in the correct order with sign changes between terms to generate three real positive roots.  The three solution masses are\n\
   used as inputs and temporarily paired with reference mass \"v\" as a dimensionless solution mass ratio.  Inputs with experimental uncertainty are selected\n\
   randomly within their uncertainty range each time phase 1 is run.  \n\
\n\
   The formula is solved for the coefficients of each mass term.  These coefficients are then multiplied by a series of possible factor combinations and the\n\
   other reference masses.  Interesting matches to these factors as defined by 'phase1_int_match_max' and 'phase1_filter' in nle-lepton.cfg are saved for\n\
   processing in phase 2.\n\
\n\
   Up to this point factoring of the coefficients has been done on each term independently of the others.  When these terms are combined into complete formulas\n\
   in phase 2 each unique combination of potential factors for each term needs to be tried.  This often results in thousands of formulas to be processed.  To\n\
   avoid processing uninteresting formulas, symmetry and complexity scores are assigned to each combination of terms with higher symmetry and lower complexity\n\
   generally meaning a simpler and more interesting formula. The configuration options 'phase2_symmetry_min' and 'phase2_complexity_max' are provided to\n\
   filter the formulas allowed to be passed to phase 2.\n\
\n\
 Phase 2 - validate proposed formulas against experimental uncertainties of outputs:\n\
   In phase 2 processing is reversed - formulas are constructed with the factors found in phase 1 and then solved for the three solution masses and/or other\n\
   outputs.  Before solving each proposed formula, the experimentally known variables used in the formula are ranked by relative standard uncertainty and the\n\
   three (in 3-term mode) with the highest uncertainty are used as outputs (solved for) with the rest used as inputs.  This allows for the lowest possible\n\
   relative uncertainty in the outputs.  Only results with all outputs within 'phase2_results_window' (default=1.1) of experimental uncertainty are shown unless\n\
   'phase2_results_always' is set to 'yes'.\n\
\n\
 This two phase approach allows for rapid scanning of possible coefficient factors without having to solve a NLE each time a different coefficient factor is\n\
 tried.  This allows for billions of coefficient factors and trillions of coefficient, exponent and mass ratio  combinations to be tested as efficiently as\n\
 possible.\n\
\n\
 Each phase1 + coefficient scan + phase2 validation is stateless: exponents and variable inputs are randomly selected and can operate indepentenly of and \n\
 parallel to other threads. Validated results from phase 2 can optionally be uploaded to a remote server as soon as they are found.  This allows for massive \n\
 scaling using inexpensive \"spot\" instances on popular cloud services.  Sample cloud remote control scripts are provided in scripts/aws in the source tree.\n\
\n\
CAUTION\n\
 Ultimately the results from nle-lepton need to be reviewed by a physicist to determine any significance or usefulness (if any) of any result.  It is possible\n\
 that NLE's are the wrong approach to this problem or that the right ingredients (formula structure, factors) are not included yet. The author fully expects\n\
 that further development will be required before a useful formula offering clues to the underlying physics of fermion mass generation is found (if ever).\n\
 Regardless, nle-lepton  has been designed to be a flexible and extendable platform for conducting such searches.  See ToDo.txt in the source tree for planned\n\
 features.\n\
");
}
