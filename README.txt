###### Usage and Description of Operation ######

usage: nle-lepton <seed> <exponentlimit> <phase1filter> <minsymmetry> <maxcomplexity>

 example: nle-lepton 0 12 5 70 65 (no external seed, exponents from -1/12 to +1/12, phase1 filter of 1x10^-5,
                             minimum phase 2 symmetry of 70, maximum phase 2 complexity of 65. (recommended values).

 example: nle-lepton 1000000 26 6 30 150 (with an external seed, exponents from -1/26 to +1/26, phase1 fitler of 1x10^-6,
                                       minimum phase 2 symmetry of 30, maximum phase 2 complexity of 150).

 This program searches for polynomial-like non-linear equations that generate the observed charged lepton masses from simple coefficients. These three-term
 (plus constant offset) formulas include rational exponents between -1 and 1 and are designed to give three positive real roots representing the charged
 lepton masses. This is unlike conventional polynomials with integer exponents that have positive and negative roots and would require twice as many
 non-constant terms.  Another interesting effect of this approach is that the charged lepton mass spectrum can be replicated with relatively small coefficients
 on each term, often in the range of 0.1 to 10.  With traditional polynomials some coefficients would be quite large and vary dramatically between terms.
 To ensure each term is dimensionless, the charged lepton masses are combined with a reference mass inside each term as a mass ratio.  The reference masses
 currently supported are the Planck mass "mp", Higgs vacuum expectation value "v", Z boson mass "mz", W boson mass "mw", and Higgs boson mass "mh0".

 Processing is broken into two phases.   Phase 1 starts with the known charged lepton masses, solves a proposed formula for the three coefficients and searches
 for interesting multipliers (factors) for those coefficients. Phase 2 reverses the process by converting the coefficient factors from phase 1 to exact
 formulas and solving them for the charged lepton masses and other outputs.  These results are then compared to their experimental values and uncertainties.

 Phase 1:
   A random combination of three exponents are selected and put in the correct order to generate three real positive roots.  The three charged lepton masses
   are used as inputs with the Tau mass being varied randomly within it's experimental uncertainty range each time phase 1 is run. The electron
   mass is temporarily used as the reference mass in each term.

   This formula is then solved for the coefficients of each of the three exponent terms (left, middle, right).  Each coefficient is then multiplied by a
   series of test multipliers (including the actual reference masses) to see if the result is close enough to a relatively small integer rational number. The
   phase 1 limit determines which results are close enough to pass to phase 2.

 Phase 2:
   Formulas are constructed with the interesting ingredients found in phase 1 and then solved for the three charged lepton masses and other variables. In
   phase 2 each exponent term is allowed to have a different mass ratio, as long as the product of the mass ratio and factors generates the correct total
   value for that term.

   Before solving each proposed formula, all of the experimentally known variables are ranked by relative standard uncertainty and the three with the highest
   uncertainty are used as outputs (solved for) with the rest used as inputs.  Only results with all outputs within experimental uncertainty are shown unless
   the ALWAYS_SHOW_RESULTS compile time directive is set in nle-lepton.h.

   Symmetry and complexity scores are assigned to each result based on the numerical values in the multiplier, with higher symmetry and lower complexity
   generally meaning a simpler formula. Minimum symmetry and maximum complexity filters are provided to restrict the formulas allowed to be solved.
   This can greatly speed up finding the most interesting formulas as cpu time is not wasted solving the much more common higher complexity formulas.

 seed            External integer used as part of a seed for srand48() along with second and microsecond clock data.
                 This can help separate threads started at the same time on the same machine have different seeds.
                 Recommended value is <cpu_index>*1000000 where 1 <= cpu_index <= 999.

 exponent limit  Maximum inverse exponent (3 allows for exponents of -1, +1, -1/2, +1/2, -1/3, +1/3).

 phase1 limit    Threshold for selecting (multiplier * coefficient) for further processing if the result is
                 this close to an integer or simple rational number. (threshold = 1x10^-inlimit)

 minsymmetry:    Minimum multiplier symmetry allowed to be processed in phase 2.

 maxcomplexity:  Maximum multiplier complexity allowed to be processed in phase 2.

###### Sample Output ######

"result" keyword for log extraction, complexity / symmetry score (lower is better), symmetry score (higher is better), complexity score (lower is better), Exponents ID, Mass ratio ID, formla hash, line number, data

result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 08, | mZ         || 9.118760000000000e+10 | 2.100e+06 | 2.303e-05 || 9.118760000000000e+10 | 2.100e+06 | 2.303e-05 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 09, | mW         || 8.037395792564896e+10 | 1.606e+06 | 1.998e-05 || 8.037900000000000e+10 | 1.200e+07 | 1.493e-04 || -5.0421e+06 | -6.2729e-05 |               |       *        |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 10, | mH0        || 1.251136342138832e+11 | 4.815e+05 | 3.849e-06 || 1.251800000000000e+11 | 1.600e+08 | 1.278e-03 || -6.6366e+07 | -5.3016e-04 |               |       *        |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 11, | sin2w      || 2.231107009828035e-01 | 4.749e-06 | 2.129e-05 || 2.231100000000000e-01 | 2.700e-04 | 1.210e-03 ||  7.0098e-07 |  3.1419e-06 |               |                |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 12, | Electron   || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 13, | Muon       || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 14, | Tau        || 1.776914965852280e+09 | 1.079e+05 | 6.073e-05 || 1.776860000000000e+09 | 1.200e+05 | 6.753e-05 ||  5.4966e+04 |  3.0934e-05 |               |       *        |
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 16, L-M+R-1=0
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 16, L='                 pi^(-1/1)                       cos2w^(-1/2) (                        pi^(-1/1)              mH0/M   )^(1/ 4)'
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 17, M='                 pi^(-1/1)                       cos2w^( 1/2) (( 1/ 3)                 pi^(-2/1)              mZ/M    )^(1/ 3)'
result, 0.7200,  75,  54, E-4-3+4, M422,   1805926051, 18, R='                 pi                              cos2w^( 1/2) (( 3/ 1) nbs^-1 2^(-1/2)                        M/mZ    )^(1/ 4)'

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

Exponents ID indicates the inverse exponent of each term and if the mass ratio is inverted (-) or not (+)

Mass ratio ID indicates the mass ratio used in each term using these indexes:

0  = M / mp
1  = M / v
2  = M / mZ
3  = M / mW
4  = M / mH0

For exmpale: 422 = M / mH0 (left term), M / mZ (middle term), M / mZ (right term)

nbv = geometric constant for n-ball volume (index is term exponent, n-ball dimension):
  nbv[0]=    1.0;
  nbv[1]=    2.0;
  nbv[2]=              M_PI;
  nbv[3]=    4.0 *     M_PI        / 3.0;
  nbv[4]=          pow(M_PI,  2.0) / 2.0;
  nbv[5]=    8.0 * pow(M_PI,  2.0) / 15.0;
  nbv[6]=          pow(M_PI,  3.0) / 6.0;
  nbv[7]=   16.0 * pow(M_PI,  3.0) / 105.0;
  nbv[8]=          pow(M_PI,  4.0) / 24.0;
  nbv[9]=   32.0 * pow(M_PI,  4.0) / 945.0;
  nbv[10]=         pow(M_PI,  5.0) / 120.0;
  nbv[11]=  64.0 * pow(M_PI,  5.0) / 10395.0;
  nbv[12]=         pow(M_PI,  6.0) / 720.0;
  nbv[13]= 128.0 * pow(M_PI,  6.0) / 135135.0;
  nbv[14]=         pow(M_PI,  7.0) / 5040.0;
  nbv[15]= 256.0 * pow(M_PI,  7.0) / 2027025.0;
  nbv[16]=         pow(M_PI,  8.0) / 40320.0;
  nbv[17]= 512.0 * pow(M_PI,  8.0) / 34459425.0;
  nbv[18]=         pow(M_PI,  9.0) / 362880.0;
  nbv[19]=1024.0 * pow(M_PI,  9.0) / 654729075.0;
  nbv[20]=         pow(M_PI, 10.0) / 3628800.0;
  nbv[21]=2048.0 * pow(M_PI, 10.0) / 13749310575.0;
  nbv[22]=         pow(M_PI, 11.0) / 39916800.0;
  nbv[23]=4096.0 * pow(M_PI, 11.0) / 316234143225.0;
  nbv[24]=         pow(M_PI, 12.0) / 479001600.0;
  nbv[25]=8192.0 * pow(M_PI, 12.0) / 7905853580625.0;
  nbv[26]=         pow(M_PI, 13.0) / 6227020800.0;

nbs = geometric constant for n-ball surface area (i.e. n-sphere, index is term exponent, n-ball dimension - 1);
  nbs[0]=     2.0;
  nbs[1]=     2.0 *     M_PI;
  nbs[2]=     4.0 *     M_PI;
  nbs[3]=     2.0 * pow(M_PI,  2.0);
  nbs[4]=     8.0 * pow(M_PI,  2.0) / 3.0;
  nbs[5]=           pow(M_PI,  3.0);
  nbs[6]=    16.0 * pow(M_PI,  3.0) / 15.0;
  nbs[7]=           pow(M_PI,  4.0) / 3.0;
  nbs[8]=    32.0 * pow(M_PI,  4.0) / 105.0;
  nbs[9]=           pow(M_PI,  5.0) / 12.0;
  nbs[10]=   64.0 * pow(M_PI,  5.0) / 945.0;
  nbs[11]=          pow(M_PI,  6.0) / 60.0;
  nbs[12]=  128.0 * pow(M_PI,  6.0) / 10395.0;
  nbs[13]=          pow(M_PI,  7.0) / 360.0;
  nbs[14]=  256.0 * pow(M_PI,  7.0) / 135135.0;
  nbs[15]=          pow(M_PI,  8.0) / 2520.0;
  nbs[16]=  512.0 * pow(M_PI,  8.0) / 2027025.0;
  nbs[17]=          pow(M_PI,  9.0) / 20160.0;
  nbs[18]= 1024.0 * pow(M_PI,  9.0) / 34459425.0;
  nbs[19]=          pow(M_PI, 10.0) / 181440.0;
  nbs[20]= 2048.0 * pow(M_PI, 10.0) / 654729075.0;
  nbs[21]=          pow(M_PI, 11.0) / 1814400.0;
  nbs[22]= 4096.0 * pow(M_PI, 11.0) / 13749310575.0;
  nbs[23]=          pow(M_PI, 12.0) / 19958400.0;
  nbs[24]= 8192.0 * pow(M_PI, 12.0) / 316234143225.0;
  nbs[25]=          pow(M_PI, 13.0) / 239500800.0;
  nbs[26]=16384.0 * pow(M_PI, 13.0) / 7905853580625.0;
