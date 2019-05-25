/*
 * Copyright (c) 2019, Kevin M. Loch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

###### Usage and Description of Operation ######

usage: lepton <seed> <exponentlimit> <phase1filter> <minsymmetry> <maxcomplexity>

 example: lepton 0 4 5 50 70 (no external seed, max exponent -1/4 or +1/4, phase1 filter of 1x10^-5,
                             minimum phase 2 symmetry of 50, maximum phase 2 complexity of 70. (recommended values).

 example: lepton 1000000 3 6 7 30 150 (with an external seed, max exponent -1/3 or +1/3, phase1 fitler of 1x10^-6,
                                       minimum phase 2 symmetry of 30, maximum phase 2 complexity of 150).

 *** Warning: very rough work in progress ***

 This program searches for polynomial-like formulas (polyforms) that might generate the charged lepton masses from simple coefficients.
 These three term (plus constant of -1) formulas use positive and negative rational exponents less than 1 and are designed to give exactly three positive real
 roots representing the lepton masses. Unlike integer exponents, rational exponents less than +/-1 generate positive-only roots with a spectrum similar to that
 of the observed masses using relatively small coefficients on each term.

 Processing is broken into two phases.   Phase 1 starts with the known charged lepton masses, solves a proposed polyform for the three coefficients and searches
 for interesting multipliers to those coefficients. Phase 2 reverses the process by converting the close multiplier matches from phase 1 to exact formulas
 and solving them for the charged lepton masses and other outputs.  These results are then filtered by comparing them to their experimental values and
 uncertainties.

 Phase 1:
   A random combination of three exponents are selected along with a mass ratio (like M/v or M/mZ etc.).   The three charged lepton masses are used as inputs
   along with the other mass used in the selected mass ratio. Any input with experimental uncertainty greater than 1x10^-7 is varied randomly within it's
   uncertainty range for each run to enable usefully tight filter limits in phase1.

   This formula is then solved for the coefficients for each of the three exponent terms (left, middle, right).  Each coefficient is then multiplied by a
   series of test multipliers to see if the result is a relatively small integer rational number.   These interesting (passing the phase 1 filter)
   coefficient * multiplier matches are then stored in a table for phase 2.  Each polyform (set of three exponents) is run through phase 1 for every mass ratio
   and all of the matches for a given polyform  are saved in the same table for phase 2 to enable mixed mass ratio formulas.

 Phase 2:
   formulas are constructed with the same exponents and various interesting multipliers from phase 1 and then solved for the three charged lepton masses and
   other variables.  Unlike in phase 1, in phase 2 each exponent term can have a different mass ratio, as long as the product of the mass ratio and multiplier
   generates the correct individual polynomial term.

   Before solving each proposed formula, all of the experimentally known variables are ranked by relative standard uncertainty and the three with the highest
   uncertainty are used as outputs (solved for) with the rest used as inputs. The resulting outputs are then checked against their experimental ranges.

   Symmetry and complexity scores are also
   assigned to each result based on the numerical values in the multiplier, with higher symmetry and lower complexity generally meaning a simpler formula.
   Minimum symmetry and maximum complexity
   filters are provided to restrict the formulas allowed to be solved.   This can greatly speed up finding the most interesting
   formulas as cpu time is not wasted solving the much more common higher complexity formulas.   The optimal maxcomplexity setting varies by exponent limit
   and the phase1 limit.

 seed            External integer used as part of a seed for srand48() along with second and microsecond clock data.
                 This can help separate threads started at the same time on the same machine have different seeds.
                 Recommended value is <cpu_index>*1000000 where 1 <= cpu_index <= 999.

 exponent limit  Maximum inverse exponent (3 allows for exponents of -1, +1, -1/2, +1/2, -1/3, +1/3).

 phase1 limit    Threshold for selecting (multiplier * coefficient) for further processing if the result is
                 this close to an integer or simple rational number. (threshold = 1x10^-inlimit)

 minsymmetry:    Minimum multiplier symmetry allowed to be processed in phase 2.

 maxcomplexity:  Maximum multiplier complexity allowed to be processed in phase 2.


###### Sample Output ######

"result" keyword for log extraction, complexity / symmetry score (lower is better), symmetry score (higher is better), complexity score (lower is better), Polyform ID, Mass ratio ID, formla hash, line number, data

result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 08, | mZ         || 9.118760000000000e+10 | 2.100e+06 | 2.303e-05 || 9.118760000000000e+10 | 2.100e+06 | 2.303e-05 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 09, | mW         || 8.037395792564896e+10 | 1.606e+06 | 1.998e-05 || 8.037900000000000e+10 | 1.200e+07 | 1.493e-04 || -5.0421e+06 | -6.2729e-05 |               |       *        |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 10, | mH0        || 1.251136342138832e+11 | 4.815e+05 | 3.849e-06 || 1.251800000000000e+11 | 1.600e+08 | 1.278e-03 || -6.6366e+07 | -5.3016e-04 |               |       *        |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 11, | sin2w      || 2.231107009828035e-01 | 4.749e-06 | 2.129e-05 || 2.231100000000000e-01 | 2.700e-04 | 1.210e-03 ||  7.0098e-07 |  3.1419e-06 |               |                |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 12, | Electron   || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 13, | Muon       || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 14, | Tau        || 1.776914965852280e+09 | 1.079e+05 | 6.073e-05 || 1.776860000000000e+09 | 1.200e+05 | 6.753e-05 ||  5.4966e+04 |  3.0934e-05 |               |       *        |
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 16,    16, '                 pi^(-1/1)                       cos2w^(-1/2) (                        pi^(-1/1)              mH0/M   )^(1/ 4)'
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 17,    19, '                 pi^(-1/1)                       cos2w^( 1/2) (( 1/ 3)                 pi^(-2/1)              mZ/M    )^(1/ 3)'
result, 0.7200,  75,  54, P-4-3+4, M422,   1805926051, 18,    19, '                 pi                              cos2w^( 1/2) (( 3/ 1) nbs^-1 2^(-1/2)                        M/mZ    )^(1/ 4)'

This particular formula would be written as:

                   left term                                                middle term                                                   right term                                     constant 
  ((1 / (pi * cos(thetaW))) * (mH0 / pi * M)^(1/4))  -  ((cos(thetaW) / pi) * (mZ / (3 * pi^2 * M))^(1/3))    +    ((pi * cos(thetaW)) * ((9 * M) / (8 * sqrt(2) pi^2 * mZ))^(1/4))    -    1       = 0


###### Abbreviations used in output ######

M     = unknown charged lepton mass (Electron, Muon, or Tau)
a     = alpha_em (fine structure constant)
sin2w = sin^2(thetaW) (W = weak mixing angle)
cos2w = cos^2(thetaW) (W = weak mixing angle)
mp    = Planck mass
v     = Higgs vacuum expecation value
mZ    = Z boson mass
mW    = W boson mass
mH0   = Higgs boson mass
G     = Netwon's gravitational constant

Polyform ID indicates the fracitonal exponent of each term and if the mass ratio is inverted (-) or not (+)

Mass ratio ID indicates the mass ratio used in each term using these indexes:

0  = M / (a^9 * mp)
1  = M / v
2  = M / mZ
3  = M / mW
4  = M / mH0

For exmpale: 422 = M / mH0 (left term), M / mZ (middle term), M / mZ (right term)

nbv = n-ball volume (index is term exponent, n-ball dimension):
  nbv[0]= 1.0;
  nbv[1]= 2.0;
  nbv[2]=            M_PI;
  nbv[3]= 4.0  *     M_PI       / 3.0;
  nbv[4]=        pow(M_PI, 2.0) / 2.0;
  nbv[5]= 8.0  * pow(M_PI, 2.0) / 15.0;
  nbv[6]=        pow(M_PI, 3.0) / 6.0;
  nbv[7]= 16.0 * pow(M_PI, 3.0) / 105.0;
  nbv[8]=        pow(M_PI, 4.0) / 24.0;
  nbv[9]= 32.0 * pow(M_PI, 4.0) / 945.0;
  nbv[10]=       pow(M_PI, 5.0) / 120.0;
  nbv[11]=64.0 * pow(M_PI, 5.0) / 10395.0;
  nbv[12]=       pow(M_PI, 6.0) / 720.0;

nbs = n-ball surface area (index is term exponent, n-ball dimension - 1);
  nbs[1]= 2.0  *      M_PI;
  nbs[2]= 4.0  *      M_PI;
  nbs[3]= 2.0  *  pow(M_PI, 2.0);
  nbs[4]= 8.0  *  pow(M_PI, 2.0) / 3.0;
  nbs[5]=         pow(M_PI, 3.0);
  nbs[6]= 16.0 *  pow(M_PI, 3.0) / 15.0;
  nbs[7]=         pow(M_PI, 4.0) / 3.0;
  nbs[8]= 32.0 *  pow(M_PI, 4.0) / 105.0;
  nbs[9]=         pow(M_PI, 5.0) / 12.0;
  nbs[10]=64.0 *  pow(M_PI, 5.0) / 945.0;
  nbs[11]=        pow(M_PI, 6.0) / 60.0;
  nbs[12]=128.0 * pow(M_PI, 6.0) / 10395.0;

