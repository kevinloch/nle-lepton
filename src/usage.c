#include <stdio.h>
#include "nle-lepton.h"

void printUsage() {

printf("\n\
nle-lepton version %s\n", NLE_VERSION);
printf("\n\
usage: nle-lepton <seed> <exponentlimit> <phase1filter> <minsymmetry> <maxcomplexity>\n\
\n\
 example: nle-lepton 0 12 5 70 65 (no external seed, exponents from -1/12 to +1/12, phase1 filter of 1x10^-5,\n\
                             minimum phase 2 symmetry of 70, maximum phase 2 complexity of 65. (recommended values).\n\
\n\
 example: nle-lepton 1000000 26 6 30 150 (with an external seed, exponents from -1/26 to +1/26, phase1 fitler of 1x10^-6,\n\
                                       minimum phase 2 symmetry of 30, maximum phase 2 complexity of 150).\n\
\n\
 This program searches for polynomial-like non-linear equations that generate the observed charged lepton masses from simple coefficients. These three-term\n\
 (plus constant offset) formulas include rational exponents between -1 and 1 and are designed to give three positive real roots representing the charged\n\
 lepton masses. This is unlike conventional polynomials with integer exponents that have positive and negative roots and would require twice as many\n\
 non-constant terms.  Another interesting effect of this approach is that the charged lepton mass spectrum can be replicated with relatively small coefficients\n\
 on each term, often in the range of 0.1 to 10.  With traditional polynomials some coefficients would be quite large and vary dramatically between terms.\n\
 To ensure each term is dimensionless, the charged lepton masses are combined with a reference mass inside each term as a mass ratio.  The reference masses\n\
 currently supported are the Planck mass \"mp\", Higgs vacuum expectation value \"v\", Z boson mass \"mz\", W boson mass \"mw\", and Higgs boson mass \"mH0\".\n\
\n\
 Processing is broken into two phases.   Phase 1 starts with the known charged lepton masses, solves a proposed formula for the three coefficients and searches\n\
 for interesting multipliers (factors) for those coefficients. Phase 2 reverses the process by converting the coefficient factors from phase 1 to exact\n\
 formulas and solving them for the charged lepton masses and other outputs.  These results are then compared to their experimental values and uncertainties.\n\
\n\
 Phase 1:\n\
   A random combination of three exponents are selected and put in the correct order to generate three real positive roots.  The three charged lepton masses\n\
   are used as inputs with the Tau mass being varied randomly within it's experimental uncertainty range each time phase 1 is run. The electron\n\
   mass is temporarily used as the reference mass in each term.\n\
\n\
   This formula is then solved for the coefficients of each of the three exponent terms (left, middle, right).  Each coefficient is then multiplied by a\n\
   series of test multipliers (including the actual reference masses) to see if the result is close enough to a relatively small integer rational number. The\n\
   phase 1 limit determines which results are close enough to pass to phase 2.\n\
\n\
 Phase 2:\n\
   Formulas are constructed with the interesting ingredients found in phase 1 and then solved for the three charged lepton masses and other variables. In\n\
   phase 2 each exponent term is allowed to have a different mass ratio, as long as the product of the mass ratio and factors generates the correct total\n\
   value for that term.\n\
\n\
   Before solving each proposed formula, all of the experimentally known variables are ranked by relative standard uncertainty and the three with the highest\n\
   uncertainty are used as outputs (solved for) with the rest used as inputs.  Only results with all outputs within experimental uncertainty are shown unless\n\
   the ALWAYS_SHOW_RESULTS compile time directive is set in nle-lepton.h.\n\
\n\
   Symmetry and complexity scores are assigned to each result based on the numerical values in the multiplier, with higher symmetry and lower complexity\n\
   generally meaning a simpler formula. Minimum symmetry and maximum complexity filters are provided to restrict the formulas allowed to be solved.\n\
   This can greatly speed up finding the most interesting formulas as cpu time is not wasted solving the much more common higher complexity formulas.\n\
\n\
 seed            External integer used as part of a seed for srand48() along with second and microsecond clock data.\n\
                 This can help separate threads started at the same time on the same machine have different seeds.\n\
                 Recommended value is <cpu_index>*1000000 where 1 <= cpu_index <= 999.\n\
\n\
 exponent limit  Maximum inverse exponent (3 allows for exponents of -1, +1, -1/2, +1/2, -1/3, +1/3).\n\
\n\
 phase1 limit    Threshold for selecting (multiplier * coefficient) for further processing if the result is\n\
                 this close to an integer or simple rational number. (threshold = 1x10^-inlimit)\n\
\n\
 minsymmetry:    Minimum multiplier symmetry allowed to be processed in phase 2.\n\
\n\
 maxcomplexity:  Maximum multiplier complexity allowed to be processed in phase 2.\n\
");
}
