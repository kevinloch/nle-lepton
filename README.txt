NAME
     nle-lepton -- Find polynomial-like non-linear equations that solve for the three charged lepton masses, or any three solution masses

SYNOPSIS
     nle-lepton [-c filename] [-h] [-s seed]
 
OPTIONS:
     -c filename
          Set configuration file name (default: nle-lepton.cfg)

     -h
          Show help

     -s seed
          Set external random seed.   This is combined with internal clock-based seeds to initialize srand48();

DESCRIPTON
 nle-lepton is a tool to search for polynomial-like non-linear equations (nle's) that generate the three charged lepton masses or any three "solution masses"
 (sm1, sm2, sm3) specified by the user.  It does this by generating a random minimal formula structure (expnents and masses), solving it for real coefficients,
 and then attempting to factor the coefficients.  When the factors are a close enough match an exact formula is contstructed from the implied factors and solved
 for the solution masses and/or other variables.  The variables are then compared to their experimental uncertainties and the result is output if they match.
 This multi-step process (solve for coefficients, factor, solve for masses) allows for billions of coefficient factors and trillions of coefficient, exponent
 and mass ratio  combinations to be tested as efficiently as possible.

 The nle's are designed to have three positive real roots corresponding to the three solution masses. This concept might explain why there are exactly three
 generations of matter and an empirical formula might provide insight into the fermion mass generation process.  The solutions masses and their uncertainties
 are user-configurable so they can be set to rest masses, running-masses or any desired values of interest.  The configuration file (by default nle-lepton.cfg)
 also has extensive options to configure operating parameters and potential coefficient factors. As an open-source project the code can be modified by the user
 and is intended to be a platform for quickly testing new formula types.

 As of version 4.2, nle-lepton has three main operating modes: 3-term, 2-term mixed, and 2-term mixed with "1-minus". The original 3-term mode uses three mass
 terms plus a constant term of -1, with different exponents and alternating sign on each term to support three real roots. Exponents are currently limited to
 fractional values between -1 and 1, excluding zero.  The use of fractional exponents allows the charged lepton mass spectrum to be replicated with coefficients
 relatively close to 1, often in the range of 0.1 to 10.  With traditional integer exponent polynomials some coefficients would be much larger or smaller and
 vary dramatically between terms. To ensure each term is dimensionless, the solution masses are combined with a reference mass inside each term as a "solution
 mass ratio" (smr).  The reference masses currently supported are the Planck mass "mp", Higgs vacuum expectation value "v", Z boson mass "mz, W boson mass "mw",
 Higgs boson mass "mH0", and user-defined mass "m_user".  

 Since fractional exponent terms might represent a relative length scale derived from a higher dimensional manifold, n-ball volume (nbv) and n-sphere surface
 area (nss) geometric constant factors are tested.  For example the geometric factor for the common 3-ball volume nbv[3]=4pi/3.  This appears in the formula
 for the radius r of a 3-ball with volume V as r=(V/nbv[3])^(1/3).

 Processing is broken into two main phases.

 Phase 1 - solve for coefficients using the known masses and then attempt to factor the coefficients:
   Randomly selected exponents put in the correct order with sign changes between terms to generate three real positive roots.  The three solution masses are
   used as inputs and temporarily paired with reference mass "v" as a dimensionless solution mass ratio.  Inputs with experimental uncertainty are selected
   randomly within their uncertainty range each time phase 1 is run.  

   The formula is solved for the coefficients of each mass term.  These coefficients are then multiplied by a series of possible factor combinations and the
   other reference masses.  Interesting matches to these factors as defined by 'phase1_int_match_max' and 'phase1_filter' in nle-lepton.cfg are saved for
   processing in phase 2.

   Up to this point factoring of the coefficients has been done on each term independently of the others.  When these terms are combined into complete formulas
   in phase 2 each unique combination of potential factors for each term needs to be tried.  This often results in thousands of formulas to be processed.  To
   avoid processing uninteresting formulas, symmetry and complexity scores are assigned to each combination of terms with higher symmetry and lower complexity
   generally meaning a simpler and more interesting formula. The configuration options 'phase2_symmetry_min' and 'phase2_complexity_max' are provided to
   filter the formulas allowed to be passed to phase 2.

 Phase 2 - validate implied formulas against experimental uncertainties of outputs:
   In phase 2 processing is reversed - formulas are constructed with the factors found in phase 1 and then solved for the three solution masses and/or other
   outputs.  Before solving each proposed formula, the experimentally known variables used in the formula are ranked by relative standard uncertainty and the
   three (in 3-term mode) with the highest uncertainty are used as outputs (solved for) with the rest used as inputs.  This allows for the lowest possible
   relative uncertainty in the outputs.  Only results with all outputs within 'phase2_results_window' (default=1.1) of experimental uncertainty are shown unless
   'phase2_results_always' is set to 'yes'.

 Each phase1 + coefficient scan + phase2 validation is stateless: exponents and variable inputs are randomly selected and can operate indepentenly of and 
 parallel to other threads. Validated results from phase 2 can optionally be uploaded to a remote server as soon as they are found.  This allows for massive 
 scaling using inexpensive "spot" instances on popular cloud services.  Sample cloud remote control scripts are provided in scripts/aws in the source tree.

 Starting with version 4.1 2-term mixed mode is supported.  In this mode two mass terms are used with a constant term of -1 and the two mass terms are mixed to
 create a synthetic third (middle) term.  While this was inspired by two terms squared: (a-b)^2 = a^2 - 2ab + b^2, the way it is implemented supports other
 related mixing modes such as full mesh (sum of products) of two copies of a and b: a^2 - 4ab + b^2, or any mixing that can be approximated by a^2 - nab + b^2
 where n is an integer >= 1. One interesting property of this mode is that the relative mass spectrum is determined only by the exponents and mixing mode.
 Since the synthecitally generated middle coefficient is not independent of the other two there is no way to solve for a three mass spectrum unless it is
 inherrently supported by the exponents and mixing mode.  This provides a powerful check of the applicability of a formula without having to factor or
 interpret the factoring of the coefficients. This test is represented in the output as "two_term_test" which represents "n" in a^2 - nab + b^2. if
 two_term_test is not sufficiently close to an integer >= 1 the formula is considered invalid and coefficient factoring is bypassed.

 Starting with version 4.2 2-term mixed mode with (1-smr) is supported.  In 2-term mixed mode with smrfactor_1minus_enable=yes, processing is similar to 2-term
 except the first term uses 1-smr inside the radical, where smr=solution mass ratio and smr factors, while the second term uses just smr.  Due to the non-linear
 nature of 1-smr specific reference masses and other smr factors need to be explicitly tested each time phase 1 is run instead of being scanned for during
 coefficient factoring.  Configuration file options smrfactor_* control these tests.

 Starting with version 4.3 reference mass ratio (rmr) factors are supported outside the radical, for example: mh0/v.  As these would most likely belong on the constant term the option phase2_check_rmr can be set to yes to require all terms to have the same rmr.

CAUTION
 Ultimately the results from nle-lepton need to be reviewed by a physicist to determine any significance or usefulness (if any) of any result.  It is possible
 that these type of formulas are the wrong approach to this problem or that the right ingredients (formula structure, factors) are not included yet. The author
 fully expects that further development will be required before a useful formula offering clues to the underlying physics of fermion mass generation is found
 (if ever). Regardless, nle-lepton has been designed to be a flexible and extendable platform for conducting such searches.  See ToDo.txt in the source tree
 for planned features.

SAMPLE 3-TERM OUTPUT
 Key: "result" keyword for log extraction, complexity / symmetry score (lower is better), symmetry score (higher is better), complexity score (lower is better), Exponents ID, Mass ratio ID, formla hash, line number, data

 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 06, | mZ         || 9.118760000000000e+10 | 2.100e+06 | 2.303e-05 || 9.118760000000000e+10 | 2.100e+06 | 2.303e-05 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 08, | mW         || 8.038038511479678e+10 | 1.866e+06 | 2.321e-05 || 8.037900000000000e+10 | 1.200e+07 | 1.493e-04 ||  1.3851e+06 |  1.7232e-05 |               |       *        |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 09, | sin2w      || 2.229864465689631e-01 | 2.863e-07 | 1.284e-06 || 2.231100000000000e-01 | 2.700e-04 | 1.210e-03 || -1.2355e-04 | -5.5378e-04 |               |       D        |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 10, | mH0        || 1.252676782196380e+11 | 2.821e+06 | 2.252e-05 || 1.253500000000000e+11 | 1.500e+08 | 1.197e-03 || -8.2322e+07 | -6.5674e-04 |               |       *        |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 12, | sm1        || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 13, | sm2        || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 14, | sm3        || 1.776889906040733e+09 | 1.122e+05 | 6.317e-05 || 1.776860000000000e+09 | 1.200e+05 | 6.753e-05 ||  2.9906e+04 |  1.6831e-05 |               |       *        |
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 16, formula: term1 - term2 + term3 - 1 = 0,
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 17, term1='( 6/ 1) 2^( 1/2) pi^(-1/1)                                                              cos2w^(-1/1) (( 3/ 1) nss^-1 2^(-1/2)                                        mW/M        )^(1/11)'
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 18, term2='( 6/ 1) 2^( 1/2) pi^(-1/1)                                                              cos2w^(-1/1) (( 1/ 4) nss^-1                                                mH0/M        )^(1/10)'
 result, 0.3846, 104,  40, E-11-10+10, M342,   4107421516, 19, term3='( 1/ 2)                                                                                 cos2w^(-1/1) (( 4/ 3) nss^-1                                                 M/mZ        )^(1/10)'

SAMPLE 2-TERM WITH (1-SMR) OUTPUT
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 04, | alpha_em   || 7.297352569300000e-03 | 1.100e-12 | 1.507e-10 || 7.297352569300000e-03 | 1.100e-12 | 1.507e-10 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 05, | v          || 2.462196262074084e+11 | 4.635e+03 | 1.882e-08 || 2.462196510000000e+11 | 6.300e+04 | 2.559e-07 || -2.4793e+04 | -1.0069e-07 |               |       *        |
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 09, | sin2w      || 2.231398664015445e-01 | 3.153e-10 | 1.413e-09 || 2.231100000000000e-01 | 2.700e-04 | 1.210e-03 ||  2.9866e-05 |  1.3386e-04 |               |       *        |
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 12, | sm1        || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 || 5.109989500000000e+05 | 1.500e-04 | 2.935e-10 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 13, | sm2        || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 || 1.056583755000000e+08 | 2.300e+00 | 2.177e-08 ||  0.0000e+00 |  0.0000e+00 |       *       |                |
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 14, | sm3        || 1.776862320205652e+09 | 3.282e+01 | 1.847e-08 || 1.776862343312000e+09 | 1.200e+05 | 6.753e-05 || -2.3106e+01 | -1.3004e-08 |               |       *        |
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 16, formula: term1^2 - (term3 * term1 * term2) + term2^2 - 1 = 0,
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 17, term1='( 9/ 2)          pi^(-1/1)                                                              cos2w^(-1/2) (        nbv^-1 2^(-1/2)                                   (1-(smrf*M/v))   )^(1/ 3)'
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 18, term2='( 4/ 3)          pi^(-2/1)                                                              cos2w^(-1/2) (( 2/ 3) nbv^-1          pi^(-1/1)                             smrf*M/v     )^(1/ 9)'
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 19, term3='( 4/ 1)                                                                                                                                                                                  '
 result, 1.1613,  31,  36, E+3+9, M11,   3299600909, 20, smrf= '( 4/ 1)                    a^(-1/1)            '=5.481439963e+02

 Note: sample outputs are shown to illustrate output format only.   No special significance or usefulness of these specific formulas are intended.

ABBREVIATIONS USED IN OUTPUT
 M     = unknown solution mass (sm1, sm2 or sm3)
 a     = alpha_em (fine structure constant)
 sin2w = sin^2(thetaW) (thetaW = weak mixing angle)
 cos2w = cos^2(thetaW) (thetaW = weak mixing angle)
 mp    = Planck mass
 v     = Higgs vacuum expecation value
 mZ    = Z boson mass
 mW    = W boson mass
 mH0   = Higgs boson mass
 G     = Netwon's gravitational constant
 uout  = user-defined constant factor outside radical
 uin   = user-defined constant factor inside radical
 muser = user-defined reference mass for mass ratios
 D     = Appears in "Used as output" column when sin2w is derived from mw and mz instead of being solved for directly
 smrf  = solution mass ratio factors

 Exponents ID indicates the denomenator of the exponent of each term and if the mass ratio is inverted (-) or not (+)

 Mass ratio ID indicates the reference mass used in a mass ratio:
  0  = M / mp
  1  = M / v
  2  = M / mZ
  3  = M / mW
  4  = M / mH0
  5  = M / m_user
  For exmpale: 422 = M / mH0 (term 1), M / mZ (term 2), M / mZ (term 3)

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

 nss = geometric constant for n-sphere surface area (index is term exponent, n-sphere dimension embedded in n+1 euclidean space);
  nss[0]=     2.0;
  nss[1]=     2.0 *     M_PI;
  nss[2]=     4.0 *     M_PI;
  nss[3]=     2.0 * pow(M_PI,  2.0);
  nss[4]=     8.0 * pow(M_PI,  2.0) / 3.0;
  nss[5]=           pow(M_PI,  3.0);
  nss[6]=    16.0 * pow(M_PI,  3.0) / 15.0;
  nss[7]=           pow(M_PI,  4.0) / 3.0;
  nss[8]=    32.0 * pow(M_PI,  4.0) / 105.0;
  nss[9]=           pow(M_PI,  5.0) / 12.0;
  nss[10]=   64.0 * pow(M_PI,  5.0) / 945.0;
  nss[11]=          pow(M_PI,  6.0) / 60.0;
  nss[12]=  128.0 * pow(M_PI,  6.0) / 10395.0;
  nss[13]=          pow(M_PI,  7.0) / 360.0;
  nss[14]=  256.0 * pow(M_PI,  7.0) / 135135.0;
  nss[15]=          pow(M_PI,  8.0) / 2520.0;
  nss[16]=  512.0 * pow(M_PI,  8.0) / 2027025.0;
  nss[17]=          pow(M_PI,  9.0) / 20160.0;
  nss[18]= 1024.0 * pow(M_PI,  9.0) / 34459425.0;
  nss[19]=          pow(M_PI, 10.0) / 181440.0;
  nss[20]= 2048.0 * pow(M_PI, 10.0) / 654729075.0;
  nss[21]=          pow(M_PI, 11.0) / 1814400.0;
  nss[22]= 4096.0 * pow(M_PI, 11.0) / 13749310575.0;
  nss[23]=          pow(M_PI, 12.0) / 19958400.0;
  nss[24]= 8192.0 * pow(M_PI, 12.0) / 316234143225.0;
  nss[25]=          pow(M_PI, 13.0) / 239500800.0;
  nss[26]=16384.0 * pow(M_PI, 13.0) / 7905853580625.0;
