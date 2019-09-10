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

//
// *** Make sure you change each instance "localhost" in this file to your desired upload server url ***
// *** The url path does not need to exist (can return 404). See grep-logs.sh for sample script to   ***
// *** extract results from apache logs                                                              ***
//

const char version[20]="3.39";

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define UPLOAD                     // upload results to web server
#define SHOWSTATUS                 // print progress updates to stdout
#define NEGATIVEEXP                // allow negative exponents
#define IGNORE_SMALL_UNCERTAINTIES // use fixed reference values for electron mass and alpha_em in phase 2
//#define SIN2W                    // enable weak mixing angle terms (much slower)

// pick one G reference
#define CODATA_G
//#define ROSI_G
//#define WIDE_G // wide=larger uncertainty

//#define DEBUG10   // show periodic status and warnings if process is slow
//#define DEBUG11   // show status on each progress step
//#define DEBUG12   // show data on each sample tried.  Caution, extremely noisy

//#define DEBUG20   // show periodic status and warnings if process is slow
//#define DEBUG21   // show status on each progress step
//#define DEBUG22   // show data on each sample tried. Caution, extremely noisy
//#define DEBUG23   // show even more data on each sample tried. Caution, extremely noisy

// reference values (CODATA 2018/derived values unless otherwise stated)
// exact values
const double c_ref=             2.997924580000E+08;
const double h_ref=             6.62607015E-34;
const double hbar_ref=          1.05457181764616E-34;
//const double hbar_ref=h_ref/(2.0 * M_PI);
const double e_ref=             1.602176634E-19;
const double ev_to_kg=          1.78266192162790E-36;
const double kg_to_ev=          5.60958860380445E35;

// experimental values
const double alpha_ref=         7.2973525693E-3;
const double alpha_ref_error=   0.0000000011E-3;
const double me_ref=            0.51099895000E6;
const double me_ref_error=      0.00000000015E6;
const double mu_ref=          105.6583755E6;
const double mu_ref_error=      0.0000023E6;
const double tau_ref=        1776.86E6;
const double tau_ref_error=     0.12E6;
const double v_ref=           246.219651E9;      // derived from CODATA 2018 G0F
const double v_ref_error=       0.000063E9;      // derived from CODATA 2018 G0F
const double mz_ref=           91.1876E9;        // pdg
const double mz_ref_error=      0.0021E9;        // pdg
const double mw_ref=           80.379E9;         // pdg
const double mw_ref_error=      0.012E9;         // pdg
const double mh0_ref=         125.18E9;          // pdg
const double mh0_ref_error=     0.16E9;          // pdg
const double sin2w_ref=         0.22311;         // derived from pdg mw/mz
const double sin2w_ref_error=   0.00027;         // derived from pdg mw/mz

#ifdef CODATA_G
const double G_ref=             6.67430E-11; 
const double G_ref_error=       0.00015E-11; 
#elif defined ROSI_G
const double G_ref=             6.67191E-11;     // rosi
const double G_ref_error=       0.00099E-11;     // rosi
#elif defined WIDE_G
const double G_ref=             6.67400E-11;     // wide
const double G_ref_error=       0.00300E-11;     // wide
#endif

typedef struct {
  int alpha_em;
  int v;
  int G;
  int mz;
  int mw;
  int mh0;
  int sin2w;
} input_use;

typedef struct {
  int upin;
  int downin;
  int piupin;
  int pidownin;
  int aupin;
  int adownin;
  int e2upin;
  int e2downin;
  int nbvupin;
  int nbsupin;
  int multcomplexity;
  double mult[27];
  input_use uses;
} multipliers;

typedef struct {
  double tau_sample;
  double G_sample;
  double mp_sample;
  double mw_sample;
  double mz_sample;
  double mh0_sample;
  double sin2w_sample;
} random_input;


typedef struct {
  int invexp;
  int massratio;
  int matchup;
  int matchdown;
  int upin;
  int downin;
  int piupin;
  int pidownin;
  int aupin;
  int adownin;
  int e2upin;
  int e2downin;
  int nbvupin;
  int nbsupin;
  int upout;
  int downout;
  int piupout;
  int pidownout;
  int aupout;
  int adownout;
  int e2upout;
  int e2downout;
  int s2wupout;
  int s2wdownout;
  int c2wupout;
  int c2wdownout;
  int matchcomplexity;
  double matchmult;
  double match;
  long long matchhash;
  input_use uses;
} matches;

void initUses(input_use *uses) {
  uses->alpha_em=0;
  uses->v=0;
  uses->G=0;
  uses->mz=0;
  uses->mw=0;
  uses->mh0=0;
  uses->sin2w=0;
}

void addUses(input_use *dest, input_use *src) {
    dest->alpha_em= dest->alpha_em || src->alpha_em;
    dest->v=               dest->v || src->v;
    dest->G=               dest->G || src->G;
    dest->mz=             dest->mz || src->mz;
    dest->mw=             dest->mw || src->mw;
    dest->mh0=           dest->mh0 || src->mh0;
    dest->sin2w=       dest->sin2w || src->sin2w;
}

void printUses(input_use *uses) {
  printf("debug, uses:\n");
  printf("debug, ------------------------------\n");
  printf("debug, alpha_em: %d\n", uses->alpha_em);
  printf("debug, v:        %d\n", uses->v);
  printf("debug, G:        %d\n", uses->G);
  printf("debug, mz:       %d\n", uses->mz);
  printf("debug, mw:       %d\n", uses->mw);
  printf("debug, mh0:      %d\n", uses->mh0);
  printf("debug, sin2w:    %d\n", uses->sin2w);
  printf("debug, ------------------------------\n");
}

char *underscore(char *str, int len) {
  int i;

  for (i=0; i<len; i++) {
    if (str[i] == 32) {
      str[i]=95;
    }
  }
  return(str);
}

// from https://en.wikipedia.org/wiki/Binary_GCD_algorithm
unsigned int gcd(unsigned int u, unsigned int v)
{
    // simple cases (termination)
    if (u == v)
        return u;

    if (u == 0)
        return v;

    if (v == 0)
        return u;

/*
    if ((v == 5) || (v == 7) || (v == 9) || (v == 10) || (v == 11)) {
        return -1;
    }
*/

    // look for factors of 2
    if (~u & 1) // u is even
    {
        if (v & 1) // v is odd
            return gcd(u >> 1, v);
        else // both u and v are even
            return gcd(u >> 1, v >> 1) << 1;
    }

    if (~v & 1) // u is odd, v is even
        return gcd(u, v >> 1);

    // reduce larger argument
    if (u > v)
        return gcd((u - v) >> 1, v);

    return gcd((v - u) >> 1, u);
}

// from https://en.wikipedia.org/wiki/Binary_GCD_algorithm
unsigned int gcd2(unsigned int u, unsigned int v)
{
    // simple cases (termination)
    if (u == v)
        return u;

    if (u == 0)
        return v-1;

    if (v == 0)
        return u;

/*
    if ((v == 5) || (v == 7) || (v == 9) || (v == 10) || (v == 11)) {
        return -1;
    }
*/

    // look for factors of 2
    if (~u & 1) // u is even
    {
        if (v & 1) // v is odd
            return gcd(u >> 1, v);
        else // both u and v are even
            return gcd(u >> 1, v >> 1) << 1;
    }

    if (~v & 1) // u is odd, v is even
        return gcd(u, v >> 1);

    // reduce larger argument
    if (u > v)
        return gcd((u - v) >> 1, v);

    return gcd((v - u) >> 1, u);
}

void checkSymmetry(int *symmetry, int left, int middle, int right) {

  if (left == middle) {
    *symmetry+=2;
  }
  if (middle == right) {
    *symmetry+=2;
  }
  if (left == right) {
    *symmetry+=2;
  }
  if (left == -middle) {
    *symmetry+=1;
  }
  if (middle == -right) {
    *symmetry+=1;
  }
  if (left == -right) {
    *symmetry+=1;
  }
}


double solvePolyforMasses(char *poly, int leftinvexp, int middleinvexp, int rightinvexp, matches *leftmatchptr, matches *middlematchptr, matches *rightmatchptr, input_use *alluses, int maxcomplexity) {
  /***********/
  /* Phase 2 */
  /***********/
  // solve polynomial like function for particle masses using the supplied coefficients, exponents,
  // the electron mass and fine-structure constant as inputs.  Return the relative difference between
  // the computed and experimental muon mass (as an indication of accuracty)
  long int samples=0;
  int i;
  int arange, merange, murange, vrange, grange, mzrange, mwrange, mh0range, taurange, sin2wrange;
  struct timespec starttime, starttime2, endtime;
  double elapsedtime;
  int floatmu;
  int floatv;
  int floattau;
  int floatg;
  int floatmz;
  int floatmw;
  int floatmh0;
  int floatsin2w;
  double score;
  char massstr[20];
  char massstrinv[20];
  char execstr[352];
  char outstr01[320];
  char outstr02[320];
  char outstr03[320];
  char outstr04[320];
  //char outstr05[320];
  char outstr06[320];
  char outstr07[320];
  char outstr08[320];
  char outstr09[320];
  char outstr10[320];
  char outstr11[320];
  char outstr12[320];
  char outstr13[320];
  char outstr14[320];
  char outstr15[320];
  char outstr16[320];
  char outstr17[320];
  char outstr18[320];
  char usedasinput[5];
  char usedasoutput[5];
  long long resulthash;
  int complexity;
  double leftexp;
  double middleexp;
  double rightexp;
  matches *matchptr;
  char leftformulastr [288];
  char middleformulastr [288];
  char rightformulastr [288];
  int upout;
  int downout;
  char e2out[16];
  char updownout[16];
  char piout[16];
  char aout[16];
  char s2wout[16];
  char c2wout[16];
  int invexp;
  char updownin[16];
  char nbin[16];
  char e2in[16];
  char piin[16];
  char ain[16];
  int unknowns;
  int mwmzmode;
  int symmetry;
  float combinedscore;
 
  double alpha_ref_relerror=alpha_ref_error / alpha_ref;
  double me_ref_relerror=me_ref_error / me_ref;
  double mu_ref_relerror=mu_ref_error / mu_ref;  
  double tau_ref_relerror=tau_ref_error / tau_ref;
  double v_ref_relerror=v_ref_error / v_ref;
  double G_ref_relerror=G_ref_error / G_ref;
  double mz_ref_relerror=mz_ref_error / mz_ref;
  double mw_ref_relerror=mw_ref_error / mw_ref;
  double mh0_ref_relerror=mh0_ref_error / mh0_ref;
  double sin2w_ref_relerror=sin2w_ref_error / sin2w_ref;

  //  mc test vars
  double r;
  double e_test=0;
  double u_test=0;
  double t_test=0;
  double precision=0;
  double precision_last=0;
  double left, middle, right;
  double leftstatic, middlestatic, rightstatic;
  double lefts2w, middles2w, rights2w;
  double leftc2w, middlec2w, rightc2w;
  double leftmassterm, middlemassterm, rightmassterm;
  double leftmeterm, middlemeterm, rightmeterm;
  double leftmuterm, middlemuterm, rightmuterm;
  double leftmtterm, middlemtterm, rightmtterm;
  double mp;
  double worst_test;
  double rangefactor;
  double rangemultiplier=5.0;
  int stalled;
  int progress;

  // these tunings affect speed and reliability, adjust with extreme care
  double testratio=25.0;              // acceptable ratios of e_test/u_test/t_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
  int ratiograceperiod=15;            // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
  int stalledlimit=500000;            // most polyforms can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
  double defaultrangemultiplier=5.0;  // lowest practical range multiplier, fastest for most polyforms
  double stalledrangemultiplier=10.0; // this value works better for slow to solve polyforms and fast polyforms that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve polyforms this will continuously trigger
  int slowcheckpoint=1000000;         // progress point to check on slow processes
  double stuckprecision=1.0E-2;       // if precision is not past this level by slowcheckpoint, try resetting

  // mc outputs
  double alpha=0;
  double alpha_last=0;
  double alpha_center=0;
  double alpha_range=0;
  double me=0;
  double me_last=0;
  double me_center=0;
  double me_range=0;
  double mu=0;
  double mu_last=0;
  double mu_center=0;
  double mu_range=0;
  double mu_range_new=0;
  double v=0;
  double v_last=0;
  double v_center=0;
  double v_range=0;
  double v_range_new=0;
  double tau=0;
  double tau_last=0;
  double tau_center=0;
  double tau_range=0;
  double tau_range_new=0;
  double G=0;
  double G_last=0;
  double G_center=0;
  double G_range=0;
  double G_range_new=0;
  double mz=0;
  double mz_last=0;
  double mz_center=0;
  double mz_range=0;
  double mz_range_new=0;
  double mw=0;
  double mw_last=0;
  double mw_center=0;
  double mw_range=0;
  double mw_range_new=0;
  double mh0=0;
  double mh0_last=0;
  double mh0_center=0;
  double mh0_range=0;
  double mh0_range_new=0;
  double sin2w=0;
  double sin2w_last=0;
  double sin2w_center=0;
  double sin2w_range=0;
  double sin2w_range_new=0;
  double cos2w=0;

  // for reporting
  double alpha_out=0;
  double alpha_out_low=1.0E30;
  double alpha_out_high=-1.0E30;
  double alpha_out_c=0;
  double alpha_out_error=0;
  double alpha_out_relerror=0;
  double alpha_out_diff=0;
  double alpha_out_reldiff=0;
  double me_out=0;
  double me_out_low=1.0E30;
  double me_out_high=-1.0E30;
  double me_out_c=0;
  double me_out_error=0;
  double me_out_relerror=0;
  double me_out_diff=0;
  double me_out_reldiff=0;
  double mu_out=0;
  double mu_out_low=1.0E30;
  double mu_out_high=-1.0E30;
  double mu_out_c=0;
  double mu_out_error=0;
  double mu_out_relerror=0;
  double mu_out_diff=0;
  double mu_out_reldiff=0;
  double tau_out=0;
  double tau_out_low=1.0E30;
  double tau_out_high=-1.0E30;
  double tau_out_c=0;
  double tau_out_error=0;
  double tau_out_relerror=0;
  double tau_out_diff=0;
  double tau_out_reldiff=0;
  double v_out=0;
  double v_out_low=1.0E30;
  double v_out_high=-1.0E30;
  double v_out_c=0;
  double v_out_error=0;
  double v_out_relerror=0;
  double v_out_diff=0;
  double v_out_reldiff=0;
  double G_out=0;
  double G_out_low=1.0E30;
  double G_out_high=-1.0E30;
  double G_out_c=0;
  double G_out_error=0;
  double G_out_relerror=0;
  double G_out_diff=0;
  double G_out_reldiff=0;
  double mz_out=0;
  double mz_out_low=1.0E30;
  double mz_out_high=-1.0E30;
  double mz_out_c=0;
  double mz_out_error=0;
  double mz_out_relerror=0;
  double mz_out_diff=0;
  double mz_out_reldiff=0;
  double mw_out=0;
  double mw_out_low=1.0E30;
  double mw_out_high=-1.0E30;
  double mw_out_c=0;
  double mw_out_error=0;
  double mw_out_relerror=0;
  double mw_out_diff=0;
  double mw_out_reldiff=0;
  double mh0_out=0;
  double mh0_out_low=1.0E30;
  double mh0_out_high=-1.0E30;
  double mh0_out_c=0;
  double mh0_out_error=0;
  double mh0_out_relerror=0;
  double mh0_out_diff=0;
  double mh0_out_reldiff=0;
  double sin2w_out=0;
  double sin2w_out_low=1.0E30;
  double sin2w_out_high=-1.0E30;
  double sin2w_out_c=0;
  double sin2w_out_error=0;
  double sin2w_out_relerror=0;
  double sin2w_out_diff=0;
  double sin2w_out_reldiff=0;

  clock_gettime(CLOCK_REALTIME, &starttime);

  leftstatic=((double)leftmatchptr->matchup / (double)leftmatchptr->matchdown) / leftmatchptr->matchmult;
  middlestatic=((double)middlematchptr->matchup / (double)middlematchptr->matchdown) / middlematchptr->matchmult;
  rightstatic=((double)rightmatchptr->matchup / (double)rightmatchptr->matchdown) / rightmatchptr->matchmult;

  leftexp = 1.0 / (double)leftinvexp;
  middleexp = 1.0 / (double)middleinvexp;
  rightexp = 1.0 / (double)rightinvexp;

  // float the three used parameters with the highest experimental uncertainty, use the rest as inputs
  unknowns=0;
  // mH0 relative uncertainty 1.3E-3
  if (alluses->mh0 == 1) {
    floatmh0=1;
    unknowns++;
  } else {
    floatmh0=0;
    mh0_center=(double)mh0_ref;
    mh0_range=(double)mh0_ref_error;
  }

  // sin2w relative uncertainty 1.2E-3, sinw: 5.9E-4, cosw: 1.7E-4
  if (alluses->sin2w == 1) {
    if ((alluses->mw == 1) || (alluses->mz == 1)) {  // check if we are explicitly using mw or mz
      mwmzmode=1; // sin2w will be derived from mw and mz
      floatsin2w=0;
      sin2w_center=(double)sin2w_ref;
      sin2w_range=(double)sin2w_ref_error;
    } else {
      mwmzmode=0;
      floatsin2w=1;
      unknowns++;
    }
  } else {
    mwmzmode=0;
    floatsin2w=0;
    sin2w_center=(double)sin2w_ref;
    sin2w_range=(double)sin2w_ref_error;
  }

#ifdef WIDE_G // higher uncertainty puts G here when this is selected
  // G(wide) relative uncertainty 6.0E-4, mp(wide): 3.0E-4
  if (alluses->G == 1) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=(double)G_ref;
    G_range=(double)G_ref_error;
  }
#endif

  // after this point we must also check for three unknowns

  // mw relative uncertainty 1.5E-4
  if (((alluses->mw == 1) || (mwmzmode == 1)) && (unknowns < 3)) {
    floatmw=1;
    unknowns++;
  } else {
    floatmw=0;
    mw_center=(double)mw_ref;
    mw_range=(double)mw_ref_error;
  }

#ifdef ROSI_G // higher uncertainty puts G here when this is selected
  // G(rosi) relative uncertainty 1.5E-4, mp(rosi): 7.4E-5
  if ((alluses->G == 1) && (unknowns < 3)) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=(double)G_ref;
    G_range=(double)G_ref_error;
  }
#endif

  // tau relative uncertainty 6.8E-5
  if (unknowns < 3) { // tau is always used but only floated if we have less than 3 uknowns so far
    floattau=1;
    unknowns++;
  } else {
    floattau=0;
    tau_center=(double)mu_ref;
    tau_range=(double)mu_ref_error;
  }

#ifdef CODATA_G
  // G(codata) relative uncertainty 4.6E-5, mp(codata): 2.3E-5
  if ((alluses->G == 1) && (unknowns < 3)) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=(double)G_ref;
    G_range=(double)G_ref_error;
  }
#endif

  // mz relative uncertainty 2.3E-5
  if (((alluses->mz == 1) || (mwmzmode == 1)) && (unknowns < 3)) {
    floatmz=1;
    unknowns++;
  } else {
    floatmz=0;
    mz_center=(double)mz_ref;
    mz_range=(double)mz_ref_error;
  }

  // v relative uncertainty 2.6E-7
  if ((alluses->v == 1) && (unknowns < 3)) {
    floatv=1;
    unknowns++;
  } else {
    floatv=0;
    v_center=(double)v_ref;
    v_range=(double)v_ref_error;
  }

  // mu relative uncertainty 2.3E-8
  if (unknowns < 3) {  // mu is always used but only floated if we have less than 3 uknowns so far
    floatmu=1;
    unknowns++;
  } else {
    floatmu=0;
    mu_center=(double)mu_ref;
    mu_range=(double)mu_ref_error;
  }

  // we always have three unknowns at this point so the rest are never floated

  // me relative uncertainty 3.0E-10
  me_center=(double)me_ref;
  me_range=(double)me_ref_error;

  // alpha relative uncertainty 1.5E-10
  alpha_center=(double)alpha_ref;
  alpha_range=(double)alpha_ref_error;

  // systematically try all non-floated input extremes
#ifdef DEBUG21
  printf("debug, Begin phase 2 input loops, polyform: %s, unknowns: %d, floatmu: %d, floatv: %d, floatmz: %d, floatg: %d, floattau: %d, floatmw: %d, floatsin2w: %d, mwmzmode: %d, floatmh0: %d\n", poly, unknowns, floatmu, floatv, floatmz, floatg, floattau, floatmw, floatsin2w, mwmzmode, floatmh0);
  printUses(alluses);
  fflush(stdout);
#endif

#ifdef IGNORE_SMALL_UNCERTAINTIES
  // this speeds up phase 2 by up to 4x
  alpha=(double)alpha_ref;
  alpha_range=alpha_ref_error;
  arange=1;
  me=(double)me_ref;
  me_range=me_ref_error;
  merange=1;
#else
  for (arange=!alluses->alpha_em; arange <= 1; arange++) {
    // alpha is never floated, it is always used as high precision input when needed
    if (arange == 0) {
      alpha=(alpha_center - alpha_range);
    } else {
      alpha=(alpha_center + alpha_range);
    }  
    for (merange=0; merange <= 1; merange++) {
      // me is never floated it is always used as high precision input
      if (merange == 0) {
        me=(me_center - me_range);
      } else {
        me=(me_center + me_range);
      }
#endif
      for (murange=floatmu; murange <= 1; murange++) {
        // mu is always used but only floated if necessary
        if (floatmu == 0) {
          if (murange == 0) {
            mu=(mu_center - mu_range);
          } else {
            mu=(mu_center + mu_range);
          }
        }
        for (vrange=!(alluses->v || floatv); vrange <= 1; vrange++) {
          if (floatv == 0) {
            if (vrange == 0) {
              v=(v_center - v_range);
            } else {
              v=(v_center + v_range);
            }
          }
          for (mzrange=!((alluses->mz || mwmzmode) || floatmz); mzrange <= 1; mzrange++) {
            if (floatmz == 0) {
              if (mzrange == 0) {
                mz=(mz_center - mz_range);
              } else {
                mz=(mz_center + mz_range);
              }
            }
            for (grange=!(alluses->G || floatg); grange <= 1; grange++) {
              if (floatg == 0) { 
                if (grange == 0) {
                  G=(G_center - G_range);
                } else {
                  G=(G_center + G_range);
                } 
                mp=(double)kg_to_ev * (double)sqrt(hbar_ref * c_ref / G);
              } 
              for (taurange=floattau; taurange <= 1; taurange++) {
                // tau is always used but only floated if necessary
                if (floattau == 0) {
                  if (taurange == 0) {
                    tau=(tau_center - tau_range);
                  } else {
                    tau=(tau_center + tau_range);
                  }
                }
                for (mwrange=!((alluses->mw || mwmzmode) || floatmw); mwrange <= 1; mwrange++) {
                  if (floatmw == 0) {
                    if (mwrange == 0) {
                      mw=(mw_center - mw_range);
                    } else {
                      mw=(mw_center + mw_range);
                    }
                  }
                  for (sin2wrange=!(alluses->sin2w || floatsin2w); sin2wrange <= 1; sin2wrange++) {
                    if (floatsin2w == 0) {
                      if (sin2wrange == 0) {
                        sin2w=(sin2w_center - sin2w_range);
                      } else {
                        sin2w=(sin2w_center + sin2w_range);
                      }
                    }
                    for (mh0range=!(alluses->mh0 || floatmh0); mh0range <= 1; mh0range++) {
                      if (floatmh0 == 0) {
                        if (mh0range == 0) {
                          mh0=(mh0_center - mh0_range);
                        } else {
                          mh0=(mh0_center + mh0_range);
                        }
                      }
#ifdef DEBUG21
                      printf("debug, Begin    phase 2 samples loop, polyform: %s, arange: %d, merange: %d, murange: %d, vrange: %d, mzrange: %d, grange: %d, taurange: %d, mwrange: %d, sin2wrange: %d, mh0range: %d\n", poly, arange, merange, murange, vrange, mzrange, grange, taurange, mwrange, sin2wrange, mh0range);
                      fflush(stdout);
#endif
                      clock_gettime(CLOCK_REALTIME, &starttime2);
                      precision_last=1.0E99;
                      //  reset mc test vars and outputs
                      if (floatmu == 1) {
                        mu_last=(double)mu_ref;
                        mu_center=(double)mu_ref;
                        mu_range=(double)mu_ref * 0.1;
                      }
                      if (floatv == 1) {
                        v_last=(double)v_ref;
                        v_center=(double)v_ref; 
                        v_range=(double)v_ref * 0.1;
                      }
                      if (floatmz == 1) {
                        mz_last=(double)mz_ref;
                        mz_center=(double)mz_ref;
                        mz_range=(double)mz_ref * 0.1;
                      }
                      if (floatg == 1) {
                        G_last=(double)G_ref;
                        G_center=(double)G_ref; 
                        G_range=(double)G_ref * 0.1;
                      }
                      if (floattau == 1) {
                        tau_last=(double)tau_ref;
                        tau_center=(double)tau_ref; 
                        tau_range=(double)tau_ref * 0.1;
                      }
                      if (floatmw == 1) {
                        mw_last=(double)mw_ref;
                        mw_center=(double)mw_ref; 
                        mw_range=(double)mw_ref * 0.1;
                      }
                      if (floatsin2w == 1) {
                        sin2w_last=(double)sin2w_ref;
                        sin2w_center=(double)sin2w_ref;
                        sin2w_range=(double)sin2w_ref * 0.1;
                      }
                      if (floatmh0 == 1) {
                        mh0_last=(double)mh0_ref;
                        mh0_center=(double)mh0_ref;
                        mh0_range=(double)mh0_ref * 0.1;
                      }
                      precision_last=1E99;
                      stalled=0;
                      progress=0;
                      rangemultiplier=defaultrangemultiplier;
                      for (samples=0; precision_last > 1.0E-11; samples++) {
                        if ((samples > 1) && ((samples % slowcheckpoint) == 0)) { // check on slow processes
#ifdef DEBUG20
                          clock_gettime(CLOCK_REALTIME, &endtime);
                          elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                          if ((samples % 10000000) == 0) { // rate limit periodic debug prints
                            printf ("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last:  %.3e, precision: %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e, mz: %.9e, mz_range: %.4e, mw: %.9e, mw_range: %.4e, sin2w: %.9e, sin2w_range: %.4e, mh0: %.9e, mh0_range: %.4e\n", poly, samples, elapsedtime, progress, rangefactor, precision_last, precision, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range);
                            fflush(stdout);
                          }
#endif

                          if ((progress == ratiograceperiod) || (precision_last > stuckprecision)) { // it's stuck, try resetting
#ifdef DEBUG20
                            clock_gettime(CLOCK_REALTIME, &endtime);
                            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                            printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, resetting\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                            fflush(stdout);
#endif
                            //  reset mc test vars and outputs
                            if (floatmu == 1) {
                              mu_last=(double)mu_ref;
                              mu_center=(double)mu_ref;
                              mu_range=(double)mu_ref * 0.1;
                            }
                            if (floatv == 1) {
                              v_last=(double)v_ref;
                              v_center=(double)v_ref;
                              v_range=(double)v_ref * 0.1;
                            }
                            if (floatmz == 1) {
                              mz_last=(double)mz_ref;
                              mz_center=(double)mz_ref;
                              mz_range=(double)mz_ref * 0.1;
                            }
                            if (floatg == 1) {
                              G_last=(double)G_ref;
                              G_center=(double)G_ref;
                              G_range=(double)G_ref * 0.1;
                            }
                            if (floattau == 1) {
                              tau_last=(double)tau_ref;
                              tau_center=(double)tau_ref;
                              tau_range=(double)tau_ref * 0.1;
                            }
                            if (floatmw == 1) {
                              mw_last=(double)mw_ref;
                              mw_center=(double)mw_ref;
                              mw_range=(double)mw_ref * 0.1;
                            }
                            if (floatsin2w == 1) {
                              sin2w_last=(double)sin2w_ref;
                              sin2w_center=(double)sin2w_ref;
                              sin2w_range=(double)sin2w_ref * 0.1;
                            }
                            if (floatmh0 == 1) {
                              mh0_last=(double)mh0_ref;
                              mh0_center=(double)mh0_ref;
                              mh0_range=(double)mh0_ref * 0.1;
                            }
                            precision_last=1E99;
                            stalled=0;
                            progress=0;
                            rangemultiplier=defaultrangemultiplier;
                          }
                        } // end check on slow processes
                        if (stalled == stalledlimit) {
                          rangemultiplier=stalledrangemultiplier; // may be a slow solution, try bigger multiplier
                          rangefactor=worst_test * rangemultiplier;
#ifdef DEBUG20
                          clock_gettime(CLOCK_REALTIME, &endtime);
                          elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                          printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, stalled\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
#endif
                          if (floattau == 1) {
                            tau_range=tau_last * rangefactor;
                          }
                          if (floatmu == 1) {
                            mu_range=mu_last * rangefactor;
                          }
                          if (floatv == 1) {
                            v_range=v_last * rangefactor;
                          }
                          if (floatg == 1) {
                            G_range=G_last * rangefactor;
                          }
                          if (floatmz == 1) {
                            mz_range=mz_last * rangefactor;
                          }
                          if (floatmw == 1) {
                            mw_range=mw_last * rangefactor;
                          }
                          if (floatmh0 == 1) {
                            mh0_range=mh0_last * rangefactor;
                          }
                          if (floatsin2w == 1) {
                            sin2w_range=sin2w_last * rangefactor;
                          }
                        }
                        stalled++;
                        // guess random values for mc outputs
                        if (floatmu == 1) {
                          r=(double)drand48();
                          mu=((mu_center - mu_range) + (r * 2.0 * mu_range));
                          i=0;
                          while ((mu < mu_ref * 0.9) || (mu > (mu_ref * 1.1))) { // sanity check on mu mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mu range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mu_last=0;
                              mu_center=(double)mu_ref;
                              mu_range=(double)mu_ref * 0.1;
                            }
                            r=(double)drand48();
                            mu=((mu_center - mu_range) + (r * 2.0 * mu_range));
                            i++;
                          }
                        }
                        if (floatv == 1) {
                          r=(double)drand48();
                          v=((v_center - v_range) + (r * 2.0 * v_range));
                          i=0;
                          while ((v < (v_ref * 0.9)) || (v > (v_ref * 1.1))) { // sanity check to help convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, v range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              v_last=0;
                              v_center=(double)v_ref;
                              v_range=(double)v_ref * 0.1;
                            }
                            r=(double)drand48();
                            v=((v_center - v_range) + (r * 2.0 * v_range));
                            i++;
                          }
                        }
                        if (floatmz == 1) {
                          r=(double)drand48();
                          mz=((mz_center - mz_range) + (r * 2.0 * mz_range));
                          i=0;
                          while ((mz < mz_ref * 0.9) || (mz > (mz_ref * 1.1))) { // sanity check on mz mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mz range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mz_last=0;
                              mz_center=(double)mz_ref;
                              mz_range=(double)mz_ref * 0.1;
                            }
                            r=(double)drand48();
                            mz=((mz_center - mz_range) + (r * 2.0 * mz_range));
                            i++;
                          }
                        }
                        if (floatg == 1) {
                          r=(double)drand48();
                          G=((G_center - G_range) + (r * 2.0 * G_range));
                          i=0;
                          while ((G < (G_ref * 0.9)) || (G > (G_ref * 1.1))) {  // sanity check to help convergance
                            if (i > 50) {  // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, G range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              G_last=0;
                              G_center=(double)G_ref;
                              G_range=(double)G_ref * 0.1;
                            }
                            r=(double)drand48();
                            G=((G_center - G_range) + (r * 2.0 * G_range));
                            i++;
                          }
                          mp=(double)kg_to_ev * (double)sqrt(hbar_ref * c_ref / G);
                        }
                        if (floattau == 1) {
                          r=(double)drand48();
                          tau=((tau_center - tau_range) + (r * 2.0 * tau_range));
                          i=0;
                          while ((tau < tau_ref * 0.9) || (tau > (tau_ref * 1.1))) { // sanity check on tau mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, tau range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              tau_last=0;
                              tau_center=(double)tau_ref;
                              tau_range=(double)tau_ref * 0.1;
                            }
                            r=(double)drand48();
                            tau=((tau_center - tau_range) + (r * 2.0 * tau_range));
                            i++;
                          }
                        }
                        if (floatmw == 1) {
                          r=(double)drand48();
                          mw=((mw_center - mw_range) + (r * 2.0 * mw_range));
                          i=0;
                          while ((mw < mw_ref * 0.9) || (mw > (mw_ref * 1.1))) { // sanity check on mw mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mw range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mw_last=0;
                              mw_center=(double)mw_ref;
                              mw_range=(double)mw_ref * 0.1;
                            }
                            r=(double)drand48();
                            mw=((mw_center - mw_range) + (r * 2.0 * mw_range));
                            i++;
                          }
                        }
                        if (alluses->sin2w == 1) {  
                          if (mwmzmode == 0) {  // float sin2w
                            r=(double)drand48();
                            sin2w=((sin2w_center - sin2w_range) + (r * 2.0 * sin2w_range));
                            i=0;
                            while ((sin2w < sin2w_ref * 0.9) || (sin2w > (sin2w_ref * 1.1))) { // sanity check on sin2w to help solve correct root and speed convergance
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &endtime);
                                elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                                printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, sin2w range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                sin2w_last=0;
                                sin2w_center=(double)sin2w_ref; 
                                sin2w_range=(double)sin2w_ref * 0.1;
                              }
                              r=(double)drand48();
                              sin2w=((sin2w_center - sin2w_range) + (r * 2.0 * sin2w_range));
                              i++;
                            }
                            cos2w=1.0 - sin2w;
                          } else { // derive sin2w from mw/mz
                            cos2w=pow((mw/mz), 2.0);
                            sin2w=1.0 - cos2w;
                          } // end mwmzmode
                        } // end alluses sin2w
                        if (floatmh0 == 1) {
                          r=(double)drand48();
                          mh0=((mh0_center - mh0_range) + (r * 2.0 * mh0_range));
                          i=0;
                          while ((mh0 < mh0_ref * 0.9) || (mh0 > (mh0_ref * 1.1))) { // sanity check on mh0 mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mh0 range error\n", poly, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mh0_last=0;
                              mh0_center=(double)mh0_ref;
                              mh0_range=(double)mh0_ref * 0.1;
                            }
                            r=(double)drand48();
                            mh0=((mh0_center - mh0_range) + (r * 2.0 * mh0_range));
                            i++;
                          }
                        }

                        if (leftmatchptr->massratio == 0) {
                          leftmassterm=mp;
                        } else if (leftmatchptr->massratio == 1) {
                          leftmassterm=v;
                        } else if (leftmatchptr->massratio == 2) {
                          leftmassterm=mz;
                        } else if (leftmatchptr->massratio == 3) {
                          leftmassterm=mw;
                        } else if (leftmatchptr->massratio == 4) {
                          leftmassterm=mh0;
                        }
                        if (middlematchptr->massratio == 0) {
                          middlemassterm=mp;
                        } else if (middlematchptr->massratio == 1) {
                          middlemassterm=v;
                        } else if (middlematchptr->massratio == 2) {
                          middlemassterm=mz;
                        } else if (middlematchptr->massratio == 3) {
                          middlemassterm=mw;
                        } else if (middlematchptr->massratio == 4) {
                          middlemassterm=mh0;
                        }
                        if (rightmatchptr->massratio == 0) {
                          rightmassterm=mp;
                        } else if (rightmatchptr->massratio == 1) {
                          rightmassterm=v;
                        } else if (rightmatchptr->massratio == 2) {
                          rightmassterm=mz;
                        } else if (rightmatchptr->massratio == 3) {
                          rightmassterm=mw;
                        } else if (rightmatchptr->massratio == 4) {
                          rightmassterm=mh0;
                        }

                        leftmeterm=me / leftmassterm;
                        leftmuterm=mu / leftmassterm;
                        leftmtterm=tau / leftmassterm;
                        middlemeterm=me / middlemassterm;
                        middlemuterm=mu / middlemassterm;
                        middlemtterm=tau / middlemassterm;
                        rightmeterm=me / rightmassterm;
                        rightmuterm=mu / rightmassterm;
                        rightmtterm=tau / rightmassterm;

                        if (leftmatchptr->s2wupout != 0) {
                          lefts2w=pow(sin2w, ((double)leftmatchptr->s2wupout / (double)leftmatchptr->s2wdownout));
                        } else {
                          lefts2w=1.0;
                        }
                        if (leftmatchptr->c2wupout != 0) {
                          leftc2w=pow(cos2w, ((double)leftmatchptr->c2wupout / (double)leftmatchptr->c2wdownout));
                        } else {
                          leftc2w=1.0;
                        }

                        if (middlematchptr->s2wupout != 0) {
                          middles2w=pow(sin2w, ((double)middlematchptr->s2wupout / (double)middlematchptr->s2wdownout));
                        } else {
                          middles2w=1.0;
                        }
                        if (middlematchptr->c2wupout != 0) {
                          middlec2w=pow(cos2w, ((double)middlematchptr->c2wupout / (double)middlematchptr->c2wdownout));
                        } else {
                          middlec2w=1.0;
                        }

                        if (rightmatchptr->s2wupout != 0) {
                          rights2w=pow(sin2w, ((double)rightmatchptr->s2wupout / (double)rightmatchptr->s2wdownout));
                        } else {
                          rights2w=1.0;
                        }
                        if (rightmatchptr->c2wupout != 0) {
                          rightc2w=pow(cos2w, ((double)rightmatchptr->c2wupout / (double)rightmatchptr->c2wdownout));
                        } else {
                          rightc2w=1.0;
                        }

                        left=  (leftstatic   / lefts2w  ) / leftc2w;
                        middle=(middlestatic / middles2w) / middlec2w;
                        right= (rightstatic  / rights2w ) / rightc2w;

                        e_test=(left * pow(leftmeterm, leftexp)) - (middle * pow(middlemeterm, middleexp)) + (right * pow(rightmeterm, rightexp)) - 1.0;
                        u_test=(left * pow(leftmuterm, leftexp)) - (middle * pow(middlemuterm, middleexp)) + (right * pow(rightmuterm, rightexp)) - 1.0;
                        t_test=(left * pow(leftmtterm, leftexp)) - (middle * pow(middlemtterm, middleexp)) + (right * pow(rightmtterm, rightexp)) - 1.0;

#ifdef DEBUG23
                        printf("debug, polyform: %s, samples: %ld, left:   %.6e, leftstatic:   %.6e, lefts2w:   %.6e, leftc2w:   %.6e, lefts2wupout:   %d, lefts2wdownout:   %d, leftc2wupout:   %d, leftc2wdownout:   %d\n", poly, samples, left, leftstatic, lefts2w, leftc2w, leftmatchptr->s2wupout, leftmatchptr->s2wdownout, leftmatchptr->c2wupout, leftmatchptr->c2wdownout);
                        printf("debug, polyform: %s, samples: %ld, middle: %.6e, middlestatic: %.6e, middles2w: %.6e, middlec2w: %.6e, middles2wupout: %d, middles2wdownout: %d, middlec2wupout: %d, middlec2wdownout: %d\n", poly, samples, middle, middlestatic, middles2w, middlec2w, middlematchptr->s2wupout, middlematchptr->s2wdownout, middlematchptr->c2wupout, middlematchptr->c2wdownout);
                        printf("debug, polyform: %s, samples: %ld, right:  %.6e, rightstatic:  %.6e, rights2w:  %.6e, rightc2w:  %.6e, rights2wupout:  %d, rights2wdownout:  %d, rightc2wupout:  %d, rightc2wdownout:  %d\n", poly, samples, right, rightstatic, rights2w, rightc2w, rightmatchptr->s2wupout, rightmatchptr->s2wdownout, rightmatchptr->c2wupout, rightmatchptr->c2wdownout);
                        printf("debug, polyform: %s, samples: %ld, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, left: %.3e, middle: %.3e, right: %.3e\n", poly, samples, e_test, u_test, t_test, left, middle, right);
                        fflush(stdout);
#endif
                        if ((progress < ratiograceperiod) || (((fabs(e_test) / fabs(u_test)) < testratio) && ((fabs(u_test) / fabs(e_test)) < testratio))) {
                          if ((progress < ratiograceperiod) || (((fabs(e_test) / fabs(t_test)) < testratio) && ((fabs(t_test) / fabs(e_test)) < testratio) &&\
                             ((fabs(u_test) / fabs(t_test)) < testratio) && ((fabs(t_test) / fabs(u_test)) < testratio))) {
                            precision=fabs(e_test) + fabs(u_test) + fabs(t_test);
#ifdef DEBUG22
                            clock_gettime(CLOCK_REALTIME, &endtime);
                            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                            printf ("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last:  %.3e, precision: %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e, mz: %.9e, mz_range: %.4e, mw: %.9e, mw_range: %.4e, sin2w: %.9e, sin2w_range: %.4e, mh0: %.9e, mh0_range: %.4e\n", poly, samples, elapsedtime, progress, rangefactor, precision_last, precision, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range);
                            fflush(stdout);
#endif
                            if (precision < precision_last) {
                              progress++;
                              stalled=0;
                              precision_last=precision;
                              alpha_last=alpha;
                              me_last=me;
                              mu_last=mu;
                              tau_last=tau;
                              v_last=v;
                              G_last=G;
                              mz_last=mz;
                              mw_last=mw;
                              mh0_last=mh0;
                              sin2w_last=sin2w;
                              if (fabs(e_test) > fabs(u_test)) {
                                worst_test=fabs(e_test);
                              } else {
                                worst_test=fabs(u_test);
                              }
                              if (fabs(t_test) > worst_test) {
                                worst_test=fabs(t_test);
                              }
                              rangemultiplier=defaultrangemultiplier;
                              rangefactor=worst_test * rangemultiplier;
                              if (rangefactor > 0.1) {
                                //  use default ranges
                                if (floatmu == 1) {
                                  mu_center=mu_last;
                                  mu_range=(double)mu_ref * 0.1;
                                }
                                if (floatv == 1) {
                                  v_center=v_last;
                                  v_range=(double)v_ref * 0.1;
                                }
                                if (floatmz == 1) {
                                  mz_center=mz_last;
                                  mz_range=(double)mz_ref * 0.1;
                                }
                                if (floatg == 1) {
                                  G_center=G_last;
                                  G_range=(double)G_ref * 0.1;
                                }
                                if (floattau == 1) {
                                  tau_center=tau_last;
                                  tau_range=(double)tau_ref * 0.1;
                                }
                                if (floatmw == 1) {
                                  mw_center=mw_last;
                                  mw_range=(double)mw_ref * 0.1;
                                }
                                if (floatsin2w == 1) {
                                  sin2w_center=sin2w_last;
                                  sin2w_range=(double)sin2w_ref * 0.1;
                                }
                                if (floatmh0 == 1) {
                                  mh0_center=mh0_last;
                                  mh0_range=(double)mh0_ref * 0.1;
                                }
                              } else { 
                                if (floattau == 1) {
                                  tau_center=tau_last;
                                  tau_range_new=tau_last * rangefactor;
                                  tau_range=((tau_range + tau_range_new + tau_range_new) / 3.0);
                                }
                                if (floatmu == 1) {
                                  mu_center=mu_last;
                                  mu_range_new=mu_last * rangefactor;
                                  mu_range=((mu_range + mu_range_new + mu_range_new) / 3.0);
                                }
                                if (floatv == 1) {
                                  v_center=v_last;
                                  v_range_new=v_last * rangefactor;
                                  v_range=((v_range + v_range_new + v_range_new) / 3.0);
                                }
                                if (floatg == 1) {
                                  G_center=G_last;
                                  G_range_new=G_last * rangefactor;
                                  G_range=((G_range + G_range_new + G_range_new) / 3.0);
                                }
                                if (floatmz == 1) {
                                  mz_center=mz_last;
                                  mz_range_new=mz_last * rangefactor;
                                  mz_range=((mz_range + mz_range_new + mz_range_new) / 3.0);
                                }
                                if (floatmw == 1) {
                                  mw_center=mw_last;
                                  mw_range_new=mw_last * rangefactor;
                                  mw_range=((mw_range + mw_range_new + mw_range_new) / 3.0);
                                }
                                if (floatmh0 == 1) {
                                  mh0_center=mh0_last;
                                  mh0_range_new=mh0_last * rangefactor;
                                  mh0_range=((mh0_range + mh0_range_new + mh0_range_new) / 3.0);
                                }
                                if (floatsin2w == 1) {
                                  sin2w_center=sin2w_last;
                                  sin2w_range_new=sin2w_last * rangefactor;
                                  sin2w_range=((sin2w_range + sin2w_range_new + sin2w_range_new) / 3.0);
                                }
                              }
#ifdef DEBUG21
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf ("debug, polyform: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last:  %.3e, precision: %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e, mz: %.9e, mz_range: %.4e, mw: %.9e, mw_range: %.4e, sin2w: %.9e, sin2w_range: %.4e, mh0: %.9e, mh0_range: %.4e\n", poly, samples, elapsedtime, progress, rangefactor, precision_last, precision, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range);
                              fflush(stdout);
#endif
                            } // end if  precision < precision_last
                          } // end if t_test/u_test/e_test
                        } // end if e_test / u_test
                      } // end for samples 

#ifdef DEBUG21
                      clock_gettime(CLOCK_REALTIME, &endtime);
                      elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                      printf("debug, Finished phase 2 samples loop, polyform: %s, samples: %ld, mass mode: %d%d%d, precision: %.6e (%6.4fs)\n", poly, samples, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, precision_last, elapsedtime);
                      fflush(stdout);
#endif
                      alpha_out=alpha_last;
                      if (alpha_out < alpha_out_low) {
                        alpha_out_low=alpha_out;
                      }
                      if (alpha_out > alpha_out_high) {
                        alpha_out_high=alpha_out;
                      }
                      alpha_out_c=((alpha_out_high + alpha_out_low) / 2.0);
                      alpha_out_error=(alpha_out_high - alpha_out_c);
                      alpha_out_relerror=alpha_out_error / alpha_out_c;
                      alpha_out_diff=alpha_out_c - alpha_ref;
                      alpha_out_reldiff=alpha_out_diff / alpha_ref;

                      me_out=me_last;
                      if (me_out < me_out_low) {
                        me_out_low=me_out;
                      }
                      if (me_out > me_out_high) {
                        me_out_high=me_out;
                      }
                      me_out_c=((me_out_high + me_out_low) / 2.0);
                      me_out_error=(me_out_high - me_out_c);
                      me_out_relerror=me_out_error / me_out_c;
                      me_out_diff=me_out_c - me_ref;
                      me_out_reldiff=me_out_diff / me_ref;

                      mu_out=mu_last;
                      if (mu_out < mu_out_low) {
                        mu_out_low=mu_out;
                      }
                      if (mu_out > mu_out_high) {
                        mu_out_high=mu_out;
                      }
                      mu_out_c=((mu_out_high + mu_out_low) / 2.0);
                      mu_out_error=(mu_out_high - mu_out_c);
                      mu_out_relerror=mu_out_error / mu_out_c;
                      mu_out_diff=mu_out_c - mu_ref;
                      mu_out_reldiff=mu_out_diff / mu_ref;

                      tau_out=tau_last;
                      if (tau_out < tau_out_low) {
                        tau_out_low=tau_out;
                      }
                      if (tau_out > tau_out_high) {
                        tau_out_high=tau_out;
                      }
                      tau_out_c=((tau_out_high + tau_out_low) / 2.0);
                      tau_out_error=(tau_out_high - tau_out_c);
                      tau_out_relerror=tau_out_error / tau_out_c;
                      tau_out_diff=tau_out_c - tau_ref;
                      tau_out_reldiff=tau_out_diff / tau_ref;

                      v_out=v_last;
                      if (v_out < v_out_low) {
                        v_out_low=v_out;
                      }
                      if (v_out > v_out_high) {
                        v_out_high=v_out;
                      }
                      v_out_c=((v_out_high + v_out_low) / 2.0);
                      v_out_error=(v_out_high - v_out_c);
                      v_out_relerror=v_out_error / v_out_c;
                      v_out_diff=v_out_c - v_ref;
                      v_out_reldiff=v_out_diff / v_ref;

                      G_out=G_last;
                      if (G_out < G_out_low) {
                        G_out_low=G_out;
                      }
                      if (G_out > G_out_high) {
                        G_out_high=G_out;
                      }
                      G_out_c=((G_out_high + G_out_low) / 2.0);
                      G_out_error=(G_out_high - G_out_c);
                      G_out_relerror=G_out_error / G_out_c;
                      G_out_diff=G_out_c - G_ref;
                      G_out_reldiff=G_out_diff / G_ref;

                      mz_out=mz_last;
                      if (mz_out < mz_out_low) {
                        mz_out_low=mz_out;
                      }
                      if (mz_out > mz_out_high) {
                        mz_out_high=mz_out;
                      }
                      mz_out_c=((mz_out_high + mz_out_low) / 2.0);
                      mz_out_error=(mz_out_high - mz_out_c);
                      mz_out_relerror=mz_out_error / mz_out_c;
                      mz_out_diff=mz_out_c - mz_ref;
                      mz_out_reldiff=mz_out_diff / mz_ref;

                      mw_out=mw_last;
                      if (mw_out < mw_out_low) {
                        mw_out_low=mw_out;
                      }
                      if (mw_out > mw_out_high) {
                        mw_out_high=mw_out;
                      }
                      mw_out_c=((mw_out_high + mw_out_low) / 2.0);
                      mw_out_error=(mw_out_high - mw_out_c);
                      mw_out_relerror=mw_out_error / mw_out_c;
                      mw_out_diff=mw_out_c - mw_ref;
                      mw_out_reldiff=mw_out_diff / mw_ref;

                      mh0_out=mh0_last;
                      if (mh0_out < mh0_out_low) {
                        mh0_out_low=mh0_out;
                      }
                      if (mh0_out > mh0_out_high) {
                        mh0_out_high=mh0_out;
                      }
                      mh0_out_c=((mh0_out_high + mh0_out_low) / 2.0);
                      mh0_out_error=(mh0_out_high - mh0_out_c);
                      mh0_out_relerror=mh0_out_error / mh0_out_c;
                      mh0_out_diff=mh0_out_c - mh0_ref;
                      mh0_out_reldiff=mh0_out_diff / mh0_ref;

                      sin2w_out=sin2w_last;
                      if (sin2w_out < sin2w_out_low) {
                        sin2w_out_low=sin2w_out;
                      }
                      if (sin2w_out > sin2w_out_high) {
                        sin2w_out_high=sin2w_out;
                      }
                      sin2w_out_c=((sin2w_out_high + sin2w_out_low) / 2.0);
                      sin2w_out_error=(sin2w_out_high - sin2w_out_c);
                      sin2w_out_relerror=sin2w_out_error / sin2w_out_c;
                      sin2w_out_diff=sin2w_out_c - sin2w_ref;
                      sin2w_out_reldiff=sin2w_out_diff / sin2w_ref;
                    } // end mh0range
                  } // end sin2wrange
                } // end mwrange
              } // end taurange
            } // end grange
          } // end mzrange
        } // end vrange
      }  // end murange
#ifndef IGNORE_SMALL_UNCERTAINTIES
    } // end merange
  }  // end arange
#endif

  score=((fmax((fabs(mu_out_reldiff) / mu_ref_relerror), 1.0) + fmax((fabs(tau_out_reldiff) / tau_ref_relerror), 1.0) + (alluses->G * fmax((fabs(G_out_reldiff) / G_ref_relerror), 1.0)) + (alluses->v * fmax((fabs(v_out_reldiff) / v_ref_relerror), 1.0)) + ((alluses->mz || mwmzmode) * fmax((fabs(mz_out_reldiff) / mz_ref_relerror), 1.0))\
       + ((alluses->mw || mwmzmode) * fmax((fabs(mw_out_reldiff) / mw_ref_relerror), 1.0)) + (alluses->mh0 * fmax((fabs(mh0_out_reldiff) / mh0_ref_relerror), 1.0)) + (alluses->sin2w * fmax((fabs(sin2w_out_reldiff) / sin2w_ref_relerror), 1.0))) / (2.0 + (double)alluses->G + (double)alluses->v + (double)(alluses->mz || mwmzmode) + (double)(alluses->mw || mwmzmode) + (double)alluses->mh0 + (double)alluses->sin2w)) - 1.0; 
  if (score == 0.0) {
    complexity=leftmatchptr->matchcomplexity + middlematchptr->matchcomplexity + rightmatchptr->matchcomplexity;
    //resulthash=lrand48();
    resulthash=(leftmatchptr->matchhash ^ middlematchptr->matchhash) ^ rightmatchptr->matchhash;

    // calculate symmetry score.   This measures how many factors are identical or inverse identical between the left, middle and right terms
    symmetry=0;
    checkSymmetry(&symmetry, (leftmatchptr->matchup * leftmatchptr->downout), (middlematchptr->matchup * middlematchptr->downout), (rightmatchptr->matchup * rightmatchptr->downout));
    checkSymmetry(&symmetry, (leftmatchptr->matchdown * leftmatchptr->upout), (middlematchptr->matchdown * middlematchptr->upout), (rightmatchptr->matchdown * rightmatchptr->upout));
    checkSymmetry(&symmetry, leftmatchptr->e2upout, middlematchptr->e2upout, rightmatchptr->e2upout);
    checkSymmetry(&symmetry, leftmatchptr->piupout, middlematchptr->piupout, rightmatchptr->piupout);
    checkSymmetry(&symmetry, (leftmatchptr->aupout * leftmatchptr->adownout), (middlematchptr->aupout * middlematchptr->adownout), (rightmatchptr->aupout * rightmatchptr->adownout));
    checkSymmetry(&symmetry, (leftmatchptr->s2wupout * leftmatchptr->s2wdownout), (middlematchptr->s2wupout * middlematchptr->s2wdownout), (rightmatchptr->s2wupout * rightmatchptr->s2wdownout));
    checkSymmetry(&symmetry, (leftmatchptr->c2wupout * leftmatchptr->c2wdownout), (middlematchptr->c2wupout * middlematchptr->c2wdownout), (rightmatchptr->c2wupout * rightmatchptr->c2wdownout));
    checkSymmetry(&symmetry, leftmatchptr->upin, middlematchptr->upin, rightmatchptr->upin);
    checkSymmetry(&symmetry, leftmatchptr->downin, middlematchptr->downin, rightmatchptr->downin);
    checkSymmetry(&symmetry, leftmatchptr->nbvupin, middlematchptr->nbvupin, rightmatchptr->nbvupin);
    checkSymmetry(&symmetry, leftmatchptr->nbsupin, middlematchptr->nbsupin, rightmatchptr->nbsupin);
    checkSymmetry(&symmetry, leftmatchptr->e2upin, middlematchptr->e2upin, rightmatchptr->e2upin);
    checkSymmetry(&symmetry, leftmatchptr->piupin, middlematchptr->piupin, rightmatchptr->piupin);
    checkSymmetry(&symmetry, (leftmatchptr->aupin * leftmatchptr->adownin), (middlematchptr->aupin * middlematchptr->adownin), (rightmatchptr->aupin * rightmatchptr->adownin));

    combinedscore = (float)complexity / (float)symmetry;

/* create formula output */
    matchptr=leftmatchptr;
    invexp=leftinvexp;
    // note: all terms except matchup and matchdown are inverted here as they represent offsets to the real coefficient
    upout=matchptr->matchup * matchptr->downout;
    downout=matchptr->matchdown * matchptr->upout;
    if ((upout == 1) && (downout == 1)) {
      sprintf(updownout, "       ");
    } else {
      sprintf(updownout, "(%2d/%2d)", upout, downout);
    }
    if (matchptr->e2upout == 0) {
      sprintf(e2out, "        ");
    } else {
      sprintf(e2out, "2^(%2d/%d)", -(matchptr->e2upout), matchptr->e2downout);
    }
    if (matchptr->piupout == 0) {
      sprintf(piout, "         ");
    } else {
      if ((matchptr->piupout == -1) && (matchptr->pidownout == 1)) {
        sprintf(piout, "pi       ");
      } else {
        sprintf(piout, "pi^(%2d/%d)", -(matchptr->piupout), matchptr->pidownout);
      }
    }
    if (matchptr->aupout == 0) {
      sprintf(aout, "        ");
    } else {
      if ((matchptr->aupout == -1) && (matchptr->adownout == 1)) {
        sprintf(aout, "a       ");
      } else {
        sprintf(aout, "a^(%2d/%d)", -(matchptr->aupout), matchptr->adownout);
      }
    }
    if (matchptr->s2wupout == 0) {
      sprintf(s2wout, "            ");
    } else {
      if ((matchptr->s2wupout == -1) && (matchptr->s2wdownout == 1)) {
        sprintf(s2wout, "sin2w       ");
      } else {
        sprintf(s2wout, "sin2w^(%2d/%d)", -(matchptr->s2wupout), matchptr->s2wdownout);
      }
    }
    if (matchptr->c2wupout == 0) { 
      sprintf(c2wout, "            ");
    } else {
      if ((matchptr->c2wupout == -1) && (matchptr->c2wdownout == 1)) { 
        sprintf(c2wout, "cos2w       ");
      } else {
        sprintf(c2wout, "cos2w^(%2d/%d)", -(matchptr->c2wupout), matchptr->c2wdownout);
      }
    }
    if ((matchptr->upin == 1) && (matchptr->downin == 1)) {
      sprintf(updownin, "       ");
    } else {
      sprintf(updownin, "(%2d/%2d)", matchptr->downin, matchptr->upin);
    }
    if ((matchptr->nbvupin == 0) && (matchptr->nbsupin == 0)) {
      sprintf(nbin, "      ");
    } else {
      if (matchptr->nbvupin == -1) {
        sprintf(nbin, "nbv   ");
      } else if (matchptr->nbvupin == 1) {
        sprintf(nbin, "nbv^-1");
      } else if (matchptr->nbsupin == -1) {
        sprintf(nbin, "nbs   ");
      } else {
        sprintf(nbin, "nbs^-1");
      }
    }
    if (matchptr->e2upin == 0) {
      sprintf(e2in, "        ");
    } else {
      sprintf(e2in, "2^(%2d/%d)", -(matchptr->e2upin), matchptr->e2downin);
    }
    if (matchptr->piupin == 0) {
      sprintf(piin, "         ");
    } else if ((matchptr->piupin == -1) && (matchptr->pidownin == 1)) {
      sprintf(piin, "pi       ");
    } else {
      sprintf(piin, "pi^(%2d/%d)", -(matchptr->piupin), matchptr->pidownin);
    }
    if (matchptr->aupin == 0) {
      sprintf(ain, "        ");
    } else if ((matchptr->aupin == -1) && (matchptr->adownin == 1)) {
      sprintf(ain, "a       ");
    } else {
      sprintf(ain, "a^(%2d/%d)", -(matchptr->aupin), matchptr->adownin);
    }
    if (matchptr->massratio == 0) {
      sprintf(massstr,       "    M/mp    ");
      sprintf(massstrinv,    "    mp/M    ");
    } else if (matchptr->massratio == 1) {
      sprintf(massstr,       "    M/v     ");
      sprintf(massstrinv,    "    v/M     ");
    } else if (matchptr->massratio == 2) {
      sprintf(massstr,       "    M/mZ    ");
      sprintf(massstrinv,    "    mZ/M    ");
    } else if (matchptr->massratio == 3) {
      sprintf(massstr,       "    M/mW    ");
      sprintf(massstrinv,    "    mW/M    ");
    } else if (matchptr->massratio == 4) {
      sprintf(massstr,       "    M/mH0   ");
      sprintf(massstrinv,    "    mH0/M   ");
    }             
    if (invexp == 1) {
      sprintf(leftformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  massstr);
    } else if (invexp == -1) {
      sprintf(leftformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv);
    } else if (invexp > 1) {
      sprintf(leftformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstr, invexp);
    } else {
      sprintf(leftformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv, -invexp);
    }

    matchptr=middlematchptr;
    invexp=middleinvexp;
    // note: all terms except matchup and matchdown are inverted here as they represent offsets to the real coefficient
    upout=matchptr->matchup * matchptr->downout;
    downout=matchptr->matchdown * matchptr->upout;
    if ((upout == 1) && (downout == 1)) {
      sprintf(updownout, "       ");
    } else {
      sprintf(updownout, "(%2d/%2d)", upout, downout);
    }
    if (matchptr->e2upout == 0) {
      sprintf(e2out, "        ");
    } else {
      sprintf(e2out, "2^(%2d/%d)", -(matchptr->e2upout), matchptr->e2downout);
    }
    if (matchptr->piupout == 0) {
      sprintf(piout, "         ");
    } else {
      if ((matchptr->piupout == -1) && (matchptr->pidownout == 1)) {
        sprintf(piout, "pi       ");
      } else {
        sprintf(piout, "pi^(%2d/%d)", -(matchptr->piupout), matchptr->pidownout);
      }
    }
    if (matchptr->aupout == 0) {
      sprintf(aout, "        ");
    } else {
      if ((matchptr->aupout == -1) && (matchptr->adownout == 1)) {
        sprintf(aout, "a       ");
      } else {
        sprintf(aout, "a^(%2d/%d)", -(matchptr->aupout), matchptr->adownout);
      }
    }
    if (matchptr->s2wupout == 0) {
      sprintf(s2wout, "            ");
    } else {
      if ((matchptr->s2wupout == -1) && (matchptr->s2wdownout == 1)) {
        sprintf(s2wout, "sin2w       ");
      } else {
        sprintf(s2wout, "sin2w^(%2d/%d)", -(matchptr->s2wupout), matchptr->s2wdownout);
      }
    }
    if (matchptr->c2wupout == 0) { 
      sprintf(c2wout, "            ");
    } else {
      if ((matchptr->c2wupout == -1) && (matchptr->c2wdownout == 1)) { 
        sprintf(c2wout, "cos2w       ");
      } else {
        sprintf(c2wout, "cos2w^(%2d/%d)", -(matchptr->c2wupout), matchptr->c2wdownout);
      }
    }
    if ((matchptr->upin == 1) && (matchptr->downin == 1)) {
      sprintf(updownin, "       ");
    } else {
      sprintf(updownin, "(%2d/%2d)", matchptr->downin, matchptr->upin);
    }
    if ((matchptr->nbvupin == 0) && (matchptr->nbsupin == 0)) {
      sprintf(nbin, "      ");
    } else {
      if (matchptr->nbvupin == -1) {
        sprintf(nbin, "nbv   ");
      } else if (matchptr->nbvupin == 1) {
        sprintf(nbin, "nbv^-1");
      } else if (matchptr->nbsupin == -1) {
        sprintf(nbin, "nbs   ");
      } else {
        sprintf(nbin, "nbs^-1");
      }
    }
    if (matchptr->e2upin == 0) {
      sprintf(e2in, "        ");
    } else {
      sprintf(e2in, "2^(%2d/%d)", -(matchptr->e2upin), matchptr->e2downin);
    }
    if (matchptr->piupin == 0) {
      sprintf(piin, "         ");
    } else if ((matchptr->piupin == -1) && (matchptr->pidownin == 1)) {
      sprintf(piin, "pi       ");
    } else {
      sprintf(piin, "pi^(%2d/%d)", -(matchptr->piupin), matchptr->pidownin);
    }
    if (matchptr->aupin == 0) {
      sprintf(ain, "        ");
    } else if ((matchptr->aupin == -1) && (matchptr->adownin == 1)) {
      sprintf(ain, "a       ");
    } else {
      sprintf(ain, "a^(%2d/%d)", -(matchptr->aupin), matchptr->adownin);
    }
    if (matchptr->massratio == 0) {
      sprintf(massstr,       "    M/mp    ");
      sprintf(massstrinv,    "    mp/M    ");
    } else if (matchptr->massratio == 1) {
      sprintf(massstr,       "    M/v     ");
      sprintf(massstrinv,    "    v/M     ");
    } else if (matchptr->massratio == 2) {
      sprintf(massstr,       "    M/mZ    ");
      sprintf(massstrinv,    "    mZ/M    ");
    } else if (matchptr->massratio == 3) {
      sprintf(massstr,       "    M/mW    ");
      sprintf(massstrinv,    "    mW/M    ");
    } else if (matchptr->massratio == 4) {
      sprintf(massstr,       "    M/mH0   ");
      sprintf(massstrinv,    "    mH0/M   ");
    }
    if (invexp == 1) {
      sprintf(middleformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  massstr);
    } else if (invexp == -1) {
      sprintf(middleformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv);
    } else if (invexp > 1) {
      sprintf(middleformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstr, invexp);
    } else {
      sprintf(middleformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv, -invexp);
    }

    matchptr=rightmatchptr;
    invexp=rightinvexp;
    // note: all terms except matchup and matchdown are inverted here as they represent offsets to the real coefficient
    upout=matchptr->matchup * matchptr->downout;
    downout=matchptr->matchdown * matchptr->upout;
    if ((upout == 1) && (downout == 1)) {
      sprintf(updownout, "       ");
    } else {
      sprintf(updownout, "(%2d/%2d)", upout, downout);
    }
    if (matchptr->e2upout == 0) {
      sprintf(e2out, "        ");
    } else {
      sprintf(e2out, "2^(%2d/%d)", -(matchptr->e2upout), matchptr->e2downout);
    }
    if (matchptr->piupout == 0) {
      sprintf(piout, "         ");
    } else {
      if ((matchptr->piupout == -1) && (matchptr->pidownout == 1)) {
        sprintf(piout, "pi       ");
      } else {
        sprintf(piout, "pi^(%2d/%d)", -(matchptr->piupout), matchptr->pidownout);
      }
    }
    if (matchptr->aupout == 0) {
      sprintf(aout, "        ");
    } else {
      if ((matchptr->aupout == -1) && (matchptr->adownout == 1)) {
        sprintf(aout, "a       ");
      } else {
        sprintf(aout, "a^(%2d/%d)", -(matchptr->aupout), matchptr->adownout);
      }
    }
    if (matchptr->s2wupout == 0) {
      sprintf(s2wout, "            ");
    } else {
      if ((matchptr->s2wupout == -1) && (matchptr->s2wdownout == 1)) {
        sprintf(s2wout, "sin2w       ");
      } else {
        sprintf(s2wout, "sin2w^(%2d/%d)", -(matchptr->s2wupout), matchptr->s2wdownout);
      }
    }
    if (matchptr->c2wupout == 0) { 
      sprintf(c2wout, "            ");
    } else {
      if ((matchptr->c2wupout == -1) && (matchptr->c2wdownout == 1)) { 
        sprintf(c2wout, "cos2w       ");
      } else {
        sprintf(c2wout, "cos2w^(%2d/%d)", -(matchptr->c2wupout), matchptr->c2wdownout);
      }
    }
    if ((matchptr->upin == 1) && (matchptr->downin == 1)) {
      sprintf(updownin, "       ");
    } else {
      sprintf(updownin, "(%2d/%2d)", matchptr->downin, matchptr->upin);
    }
    if ((matchptr->nbvupin == 0) && (matchptr->nbsupin == 0)) {
      sprintf(nbin, "      ");
    } else {
      if (matchptr->nbvupin == -1) {
        sprintf(nbin, "nbv   ");
      } else if (matchptr->nbvupin == 1) {
        sprintf(nbin, "nbv^-1");
      } else if (matchptr->nbsupin == -1) {
        sprintf(nbin, "nbs   ");
      } else {
        sprintf(nbin, "nbs^-1");
      }
    }
    if (matchptr->e2upin == 0) {
      sprintf(e2in, "        ");
    } else {
      sprintf(e2in, "2^(%2d/%d)", -(matchptr->e2upin), matchptr->e2downin);
    }
    if (matchptr->piupin == 0) {
      sprintf(piin, "         ");
    } else if ((matchptr->piupin == -1) && (matchptr->pidownin == 1)) {
      sprintf(piin, "pi       ");
    } else {
      sprintf(piin, "pi^(%2d/%d)", -(matchptr->piupin), matchptr->pidownin);
    }
    if (matchptr->aupin == 0) {
      sprintf(ain, "        ");
    } else if ((matchptr->aupin == -1) && (matchptr->adownin == 1)) {
      sprintf(ain, "a       ");
    } else {
      sprintf(ain, "a^(%2d/%d)", -(matchptr->aupin), matchptr->adownin);
    }
    if (matchptr->massratio == 0) {
      sprintf(massstr,       "    M/mp    ");
      sprintf(massstrinv,    "    mp/M    ");
    } else if (matchptr->massratio == 1) {
      sprintf(massstr,       "    M/v     ");
      sprintf(massstrinv,    "    v/M     ");
    } else if (matchptr->massratio == 2) {
      sprintf(massstr,       "    M/mZ    ");
      sprintf(massstrinv,    "    mZ/M    ");
    } else if (matchptr->massratio == 3) {
      sprintf(massstr,       "    M/mW    ");
      sprintf(massstrinv,    "    mW/M    ");
    } else if (matchptr->massratio == 4) {
      sprintf(massstr,       "    M/mH0   ");
      sprintf(massstrinv,    "    mH0/M   ");
    }             
    if (invexp == 1) {
      sprintf(rightformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  massstr);
    } else if (invexp == -1) {
      sprintf(rightformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv);
    } else if (invexp > 1) {
      sprintf(rightformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstr, invexp);
    } else {
      sprintf(rightformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv, -invexp);
    }
/* end create formula output */

    sprintf(outstr01, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr01);
    sprintf(outstr02, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr02);
    sprintf(outstr03, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr03);
    if (alluses->alpha_em == 1) {
      sprintf(usedasinput, "*");
      sprintf(outstr04, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 04, | alpha_em   || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |                |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, alpha_out_c, alpha_out_error, alpha_out_relerror, alpha_ref, alpha_ref_error, alpha_ref_relerror, alpha_out_diff, alpha_out_reldiff, usedasinput);
      printf("%s\n", outstr04);
    } else {
      outstr04[0]=0;
    }
    if (alluses->v == 1) {
      if (floatv == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr06, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 06, | v          || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, v_out_c, v_out_error, v_out_relerror, v_ref, v_ref_error, v_ref_relerror, v_out_diff, v_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr06);
    } else {
      outstr06[0]=0;
    }
    if ((alluses->mz == 1) || (mwmzmode == 1)) {
      if (floatmz == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr07, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 08, | mZ         || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mz_out_c, mz_out_error, mz_out_relerror, mz_ref, mz_ref_error, mz_ref_relerror, mz_out_diff, mz_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr07);
    } else {
      outstr07[0]=0;
    }
    if (alluses->G == 1) {
      if (floatg == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr08, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 07, | G          || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, G_out_c, G_out_error, G_out_relerror, G_ref, G_ref_error, G_ref_relerror, G_out_diff, G_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr08);
    } else {
      outstr08[0]=0;
    }
    if ((alluses->mw == 1) || (mwmzmode)) {
      if (floatmw == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr09, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 09, | mW         || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mw_out_c, mw_out_error, mw_out_relerror, mw_ref, mw_ref_error, mw_ref_relerror, mw_out_diff, mw_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr09);
    } else {
      outstr09[0]=0;
    }
    if (alluses->sin2w == 1) {
      if (floatsin2w == 1 ) {
        sprintf(usedasoutput, "*");
      } else {
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr10, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 11, | sin2w      || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |               |       %s        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, sin2w_out_c, sin2w_out_error, sin2w_out_relerror, sin2w_ref, sin2w_ref_error, sin2w_ref_relerror, sin2w_out_diff, sin2w_out_reldiff, usedasoutput);
      printf("%s\n", outstr10);
    } else {
      outstr10[0]=0;
    }
    if (alluses->mh0 == 1) {
      sprintf(usedasoutput, "*");
      sprintf(outstr11, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 10, | mH0        || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |               |       %s        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mh0_out_c, mh0_out_error, mh0_out_relerror, mh0_ref, mh0_ref_error, mh0_ref_relerror, mh0_out_diff, mh0_out_reldiff, usedasoutput);
      printf("%s\n", outstr11);
    } else {
      outstr11[0]=0;
    }
    sprintf(outstr12, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 12, | Electron   || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       *       |                |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, me_out_c, me_out_error, me_out_relerror, me_ref, me_ref_error, me_ref_relerror, me_out_diff, me_out_reldiff);
    printf("%s\n", outstr12);
    if (floatmu == 1) {
      sprintf(usedasinput, " ");
      sprintf(usedasoutput, "*");
    } else {
      sprintf(usedasinput, "*");
      sprintf(usedasoutput, " ");
    }
    sprintf(outstr13, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 13, | Muon       || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mu_out_c, mu_out_error, mu_out_relerror, mu_ref, mu_ref_error, mu_ref_relerror, mu_out_diff, mu_out_reldiff, usedasinput, usedasoutput);
    printf("%s\n", outstr13);
    sprintf(outstr14, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 14, | Tau        || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |               |       *        |", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, tau_out_c, tau_out_error, tau_out_relerror, tau_ref, tau_ref_error, tau_ref_relerror, tau_out_diff, tau_out_reldiff);
    printf("%s\n", outstr14);
    sprintf(outstr15, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr15);
    sprintf(outstr16, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 16, %5d, %s", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, leftmatchptr->matchcomplexity, leftformulastr);
    printf("%s\n", outstr16);
    sprintf(outstr17, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 17, %5d, %s", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, middlematchptr->matchcomplexity, middleformulastr);
    printf("%s\n", outstr17);
    sprintf(outstr18, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 18, %5d, %s", combinedscore, symmetry, complexity, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, rightmatchptr->matchcomplexity, rightformulastr);
    printf("%s\n", outstr18);
    fflush(stdout);
#ifdef UPLOAD
    if (complexity <= maxcomplexity) {
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr01, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr02, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr03, 320));
      system(execstr);
      if (outstr04[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr04, 320));
        system(execstr);
      }
/*
      if (outstr05[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr05, 320));
        system(execstr);
      }
*/
      if (outstr06[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr06, 320));
        system(execstr);
      }
      if (outstr07[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr07, 320));
        system(execstr);
      }
      if (outstr08[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr08, 320));
        system(execstr);
      }
      if (outstr09[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr09, 320));
        system(execstr);
      }
      if (outstr10[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr10, 320));
        system(execstr);
      }
      if (outstr11[0] != 0) {
        sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr11, 320));
        system(execstr);
      }
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr12, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr13, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr14, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr15, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr16, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr17, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr18, 320));
      system(execstr);
    } // end if complexity
#endif
  } // end if score

  return(precision_last);
}

int verifyMatches(matches *matchstart, int *nummatches, char *poly, int leftinvexp, int middleinvexp, int rightinvexp, int random_input_count, int minsymmetry, int maxcomplexity) {
  //  For polynomials with interesting coefficients on all three exponent terms,
  //  separate the match table into a separate list for each exponent, then test all unique combinations of coefficients
  //  for accuracy by comparing the computed muon mass to it's experimental value.  Print results matching a minimum threshold of interest.
  int i,j;
  matches *match;
  matches *leftmatches;
  matches *middlematches;
  matches *rightmatches;
  matches *leftmatchptr;
  matches *middlematchptr;
  matches *rightmatchptr;
  matches *tmpmatchptr;
  int numleftmatches, nummiddlematches, numrightmatches;
  int l,m,r;
  int dupe;
  struct timespec starttime;
  struct timespec endtime;
  double elapsedtime;
  double precision;
  int tmpmatchup;
  int tmpmatchdown;
  int tmpmatchcomplexity;
  long long tmphash;
  long totalcombos;
  long combo;
  input_use leftuses;
  input_use middleuses;
  input_use rightuses;
  input_use alluses;
  int complexity;
  int symmetry;

  leftmatches = (matches *)malloc(100000 * sizeof(matches));
  middlematches = (matches *)malloc(100000 * sizeof(matches));
  rightmatches = (matches *)malloc(100000 * sizeof(matches));

  // extract all left term coefficients
  match=matchstart;
  leftmatchptr=leftmatches;
  numleftmatches=0;
  for (i=0; i<*nummatches; i++) {
    if (match->invexp == leftinvexp) {
     // determine integer/rational match value
      if (match->match > 1.0) {
       tmpmatchup=(int)(match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / match->match) + 0.5);
      }
      tmphash=(long long)match->massratio ^ ((long long)((((tmpmatchup / tmpmatchdown) * (1.0 / match->matchmult)) * (long long)1.0E9) + 0.5));
      tmpmatchcomplexity=(match->matchcomplexity + tmpmatchup + tmpmatchdown);
      // search existing match table for dupes and see if we have lower complexity
      tmpmatchptr=leftmatches;
      dupe=0;
      for (j=0; j< numleftmatches; j++) {
        if (tmphash == tmpmatchptr->matchhash) {
          if (tmpmatchcomplexity < tmpmatchptr->matchcomplexity) {
            // replace
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->matchup=tmpmatchup;
            tmpmatchptr->matchdown=tmpmatchdown;
            tmpmatchptr->invexp=match->invexp;
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->upin=match->upin;
            tmpmatchptr->downin=match->downin;
            tmpmatchptr->piupin=match->piupin;
            tmpmatchptr->pidownin=match->pidownin;
            tmpmatchptr->aupin=match->aupin;
            tmpmatchptr->adownin=match->adownin;
            tmpmatchptr->e2upin=match->e2upin;
            tmpmatchptr->e2downin=match->e2downin;
            tmpmatchptr->nbvupin=match->nbvupin;
            tmpmatchptr->nbsupin=match->nbsupin;
            tmpmatchptr->upout=match->upout;
            tmpmatchptr->downout=match->downout;
            tmpmatchptr->piupout=match->piupout;
            tmpmatchptr->pidownout=match->pidownout;
            tmpmatchptr->aupout=match->aupout;
            tmpmatchptr->adownout=match->adownout;
            tmpmatchptr->e2upout=match->e2upout;
            tmpmatchptr->e2downout=match->e2downout;
            tmpmatchptr->s2wupout=match->s2wupout;
            tmpmatchptr->s2wdownout=match->s2wdownout;
            tmpmatchptr->c2wupout=match->c2wupout;
            tmpmatchptr->c2wdownout=match->c2wdownout;
            tmpmatchptr->matchcomplexity=tmpmatchcomplexity;
            tmpmatchptr->matchmult=match->matchmult;
            tmpmatchptr->match=match->match;
            tmpmatchptr->matchhash=tmphash;
            initUses(&tmpmatchptr->uses);
            addUses(&tmpmatchptr->uses, &match->uses);
          }
          dupe=1;
          break;
        }
        tmpmatchptr++;
      } // end for j
      if (dupe ==0) {
        leftmatchptr->massratio=match->massratio;
        leftmatchptr->matchup=tmpmatchup;
        leftmatchptr->matchdown=tmpmatchdown;
        leftmatchptr->invexp=match->invexp;
        leftmatchptr->massratio=match->massratio;
        leftmatchptr->upin=match->upin;
        leftmatchptr->downin=match->downin;
        leftmatchptr->piupin=match->piupin;
        leftmatchptr->pidownin=match->pidownin;
        leftmatchptr->aupin=match->aupin;
        leftmatchptr->adownin=match->adownin;
        leftmatchptr->e2upin=match->e2upin;
        leftmatchptr->e2downin=match->e2downin;
        leftmatchptr->nbvupin=match->nbvupin;
        leftmatchptr->nbsupin=match->nbsupin;
        leftmatchptr->upout=match->upout;
        leftmatchptr->downout=match->downout;
        leftmatchptr->piupout=match->piupout;
        leftmatchptr->pidownout=match->pidownout;
        leftmatchptr->aupout=match->aupout;
        leftmatchptr->adownout=match->adownout;
        leftmatchptr->e2upout=match->e2upout;
        leftmatchptr->e2downout=match->e2downout;
        leftmatchptr->s2wupout=match->s2wupout;
        leftmatchptr->s2wdownout=match->s2wdownout;
        leftmatchptr->c2wupout=match->c2wupout;
        leftmatchptr->c2wdownout=match->c2wdownout;
        leftmatchptr->matchcomplexity=tmpmatchcomplexity;
        leftmatchptr->matchmult=match->matchmult;
        leftmatchptr->match=match->match;
        leftmatchptr->matchhash=tmphash;
        initUses(&tmpmatchptr->uses);
        addUses(&tmpmatchptr->uses, &match->uses);
        numleftmatches++;
        leftmatchptr++;
      } // end if not dupe
    }  // end if invexp
    match++;
  } // end for i

  // extract all middle term coefficients
  match=matchstart;
  middlematchptr=middlematches;
  nummiddlematches=0;
  for (i=0; i<*nummatches; i++) {
    if (match->invexp == middleinvexp) {
     // determine integer/rational match value
      if (match->match > 1.0) {
       tmpmatchup=(int)(match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / match->match) + 0.5);
      }
      tmphash=(long long)match->massratio ^ ((long long)((((tmpmatchup / tmpmatchdown) * (1.0 / match->matchmult)) * (long long)1.0E9) + 0.5));
      tmpmatchcomplexity=(match->matchcomplexity + tmpmatchup + tmpmatchdown);
      // search existing match table for dupes and see if we have lower complexity
      tmpmatchptr=middlematches;
      dupe=0;
      for (j=0; j< nummiddlematches; j++) {
        if (tmphash == tmpmatchptr->matchhash) {
          if (tmpmatchcomplexity < tmpmatchptr->matchcomplexity) {
            // replace
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->matchup=tmpmatchup;
            tmpmatchptr->matchdown=tmpmatchdown;
            tmpmatchptr->invexp=match->invexp;
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->upin=match->upin;
            tmpmatchptr->downin=match->downin;
            tmpmatchptr->piupin=match->piupin;
            tmpmatchptr->pidownin=match->pidownin;
            tmpmatchptr->aupin=match->aupin;
            tmpmatchptr->adownin=match->adownin;
            tmpmatchptr->e2upin=match->e2upin;
            tmpmatchptr->e2downin=match->e2downin;
            tmpmatchptr->nbvupin=match->nbvupin;
            tmpmatchptr->nbsupin=match->nbsupin;
            tmpmatchptr->upout=match->upout;
            tmpmatchptr->downout=match->downout;
            tmpmatchptr->piupout=match->piupout;
            tmpmatchptr->pidownout=match->pidownout;
            tmpmatchptr->aupout=match->aupout;
            tmpmatchptr->adownout=match->adownout;
            tmpmatchptr->e2upout=match->e2upout;
            tmpmatchptr->e2downout=match->e2downout;
            tmpmatchptr->s2wupout=match->s2wupout;
            tmpmatchptr->s2wdownout=match->s2wdownout;
            tmpmatchptr->c2wupout=match->c2wupout;
            tmpmatchptr->c2wdownout=match->c2wdownout;
            tmpmatchptr->matchcomplexity=tmpmatchcomplexity;
            tmpmatchptr->matchmult=match->matchmult;
            tmpmatchptr->match=match->match;
            tmpmatchptr->matchhash=tmphash;
            initUses(&tmpmatchptr->uses);
            addUses(&tmpmatchptr->uses, &match->uses);
          }
          dupe=1;
          break;
        }
        tmpmatchptr++;
      } // end for j
      if (dupe ==0) {
        middlematchptr->massratio=match->massratio;
        middlematchptr->matchup=tmpmatchup;
        middlematchptr->matchdown=tmpmatchdown;
        middlematchptr->invexp=match->invexp;
        middlematchptr->massratio=match->massratio;
        middlematchptr->upin=match->upin;
        middlematchptr->downin=match->downin;
        middlematchptr->piupin=match->piupin;
        middlematchptr->pidownin=match->pidownin;
        middlematchptr->aupin=match->aupin;
        middlematchptr->adownin=match->adownin;
        middlematchptr->e2upin=match->e2upin;
        middlematchptr->e2downin=match->e2downin;
        middlematchptr->nbvupin=match->nbvupin;
        middlematchptr->nbsupin=match->nbsupin;
        middlematchptr->upout=match->upout;
        middlematchptr->downout=match->downout;
        middlematchptr->piupout=match->piupout;
        middlematchptr->pidownout=match->pidownout;
        middlematchptr->aupout=match->aupout;
        middlematchptr->adownout=match->adownout;
        middlematchptr->e2upout=match->e2upout;
        middlematchptr->e2downout=match->e2downout;
        middlematchptr->s2wupout=match->s2wupout;
        middlematchptr->s2wdownout=match->s2wdownout;
        middlematchptr->c2wupout=match->c2wupout;
        middlematchptr->c2wdownout=match->c2wdownout;
        middlematchptr->matchcomplexity=tmpmatchcomplexity;
        middlematchptr->matchmult=match->matchmult;
        middlematchptr->match=match->match;
        middlematchptr->matchhash=tmphash;
        initUses(&tmpmatchptr->uses);
        addUses(&tmpmatchptr->uses, &match->uses);
        nummiddlematches++;
        middlematchptr++;
      } // end if not dupe
    }  // end if invexp
    match++;
  } // end for i

  // extract all right term coefficients
  match=matchstart;
  rightmatchptr=rightmatches;
  numrightmatches=0;
  for (i=0; i<*nummatches; i++) {
    if (match->invexp == rightinvexp) {
     // determine integer/rational match value
      if (match->match > 1.0) {
       tmpmatchup=(int)(match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / match->match) + 0.5);
      }
      tmphash=(long long)match->massratio ^ ((long long)((((tmpmatchup / tmpmatchdown) * (1.0 / match->matchmult)) * (long long)1.0E9) + 0.5));
      tmpmatchcomplexity=(match->matchcomplexity + tmpmatchup + tmpmatchdown);
      // search existing match table for dupes and see if we have lower complexity
      tmpmatchptr=rightmatches;
      dupe=0;
      for (j=0; j< numrightmatches; j++) {
        if (tmphash == tmpmatchptr->matchhash) {
          if (tmpmatchcomplexity < tmpmatchptr->matchcomplexity) {
            // replace
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->matchup=tmpmatchup;
            tmpmatchptr->matchdown=tmpmatchdown;
            tmpmatchptr->invexp=match->invexp;
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->upin=match->upin;
            tmpmatchptr->downin=match->downin;
            tmpmatchptr->piupin=match->piupin;
            tmpmatchptr->pidownin=match->pidownin;
            tmpmatchptr->aupin=match->aupin;
            tmpmatchptr->adownin=match->adownin;
            tmpmatchptr->e2upin=match->e2upin;
            tmpmatchptr->e2downin=match->e2downin;
            tmpmatchptr->nbvupin=match->nbvupin;
            tmpmatchptr->nbsupin=match->nbsupin;
            tmpmatchptr->upout=match->upout;
            tmpmatchptr->downout=match->downout;
            tmpmatchptr->piupout=match->piupout;
            tmpmatchptr->pidownout=match->pidownout;
            tmpmatchptr->aupout=match->aupout;
            tmpmatchptr->adownout=match->adownout;
            tmpmatchptr->e2upout=match->e2upout;
            tmpmatchptr->e2downout=match->e2downout;
            tmpmatchptr->s2wupout=match->s2wupout;
            tmpmatchptr->s2wdownout=match->s2wdownout;
            tmpmatchptr->c2wupout=match->c2wupout;
            tmpmatchptr->c2wdownout=match->c2wdownout;
            tmpmatchptr->matchcomplexity=tmpmatchcomplexity;
            tmpmatchptr->matchmult=match->matchmult;
            tmpmatchptr->match=match->match;
            tmpmatchptr->matchhash=tmphash;
            initUses(&tmpmatchptr->uses);
            addUses(&tmpmatchptr->uses, &match->uses);
          }
          dupe=1;
          break;
        }
        tmpmatchptr++;
      } // end for j
      if (dupe ==0) {
        rightmatchptr->massratio=match->massratio;
        rightmatchptr->matchup=tmpmatchup;
        rightmatchptr->matchdown=tmpmatchdown;
        rightmatchptr->invexp=match->invexp;
        rightmatchptr->massratio=match->massratio;
        rightmatchptr->upin=match->upin;
        rightmatchptr->downin=match->downin;
        rightmatchptr->piupin=match->piupin;
        rightmatchptr->pidownin=match->pidownin;
        rightmatchptr->aupin=match->aupin;
        rightmatchptr->adownin=match->adownin;
        rightmatchptr->e2upin=match->e2upin;
        rightmatchptr->e2downin=match->e2downin;
        rightmatchptr->nbvupin=match->nbvupin;
        rightmatchptr->nbsupin=match->nbsupin;
        rightmatchptr->upout=match->upout;
        rightmatchptr->downout=match->downout;
        rightmatchptr->piupout=match->piupout;
        rightmatchptr->pidownout=match->pidownout;
        rightmatchptr->aupout=match->aupout;
        rightmatchptr->adownout=match->adownout;
        rightmatchptr->e2upout=match->e2upout;
        rightmatchptr->e2downout=match->e2downout;
        rightmatchptr->s2wupout=match->s2wupout;
        rightmatchptr->s2wdownout=match->s2wdownout;
        rightmatchptr->c2wupout=match->c2wupout;
        rightmatchptr->c2wdownout=match->c2wdownout;
        rightmatchptr->matchcomplexity=tmpmatchcomplexity;
        rightmatchptr->matchmult=match->matchmult;
        rightmatchptr->match=match->match;
        rightmatchptr->matchhash=tmphash;
        initUses(&tmpmatchptr->uses);
        addUses(&tmpmatchptr->uses, &match->uses);
        numrightmatches++;
        rightmatchptr++;
      } // end if not dupe
    }  // end if invexp
    match++;
  } // end for i

  // send all combinations of left, middle, right to solution function
  leftmatchptr=leftmatches;
  combo=0;
#ifdef SHOWSTATUS
  totalcombos=(numleftmatches * nummiddlematches * numrightmatches);
#endif
#ifdef SHOWSTATUS
  printf("status, Solving phase 2 formulas for masses, random input: %d, polyform %s,                 progress: total (0/%ld) left (0/%d) middle (0/%d) right (0/%d)\n", random_input_count, poly, totalcombos, numleftmatches, nummiddlematches, numrightmatches);
  fflush(stdout);
#endif
  for (l=0; l<numleftmatches; l++) {
    initUses(&leftuses);
    addUses(&leftuses, &leftmatchptr->uses);
    middlematchptr=middlematches;
    for (m=0; m<nummiddlematches; m++) {
      initUses(&middleuses);
      addUses(&middleuses, &middlematchptr->uses);
      rightmatchptr=rightmatches;
      for (r=0; r<numrightmatches; r++) {
        combo++;
        complexity=leftmatchptr->matchcomplexity + middlematchptr->matchcomplexity + rightmatchptr->matchcomplexity;
        // calculate symmetry score.   This measures how many factors are identical or inverse identical between the left, middle and right terms
        symmetry=0;
        checkSymmetry(&symmetry, (leftmatchptr->matchup * leftmatchptr->downout), (middlematchptr->matchup * middlematchptr->downout), (rightmatchptr->matchup * rightmatchptr->downout));
        checkSymmetry(&symmetry, (leftmatchptr->matchdown * leftmatchptr->upout), (middlematchptr->matchdown * middlematchptr->upout), (rightmatchptr->matchdown * rightmatchptr->upout));
        checkSymmetry(&symmetry, leftmatchptr->e2upout, middlematchptr->e2upout, rightmatchptr->e2upout);
        checkSymmetry(&symmetry, leftmatchptr->piupout, middlematchptr->piupout, rightmatchptr->piupout);
        checkSymmetry(&symmetry, (leftmatchptr->aupout * leftmatchptr->adownout), (middlematchptr->aupout * middlematchptr->adownout), (rightmatchptr->aupout * rightmatchptr->adownout));
        checkSymmetry(&symmetry, (leftmatchptr->s2wupout * leftmatchptr->s2wdownout), (middlematchptr->s2wupout * middlematchptr->s2wdownout), (rightmatchptr->s2wupout * rightmatchptr->s2wdownout));
        checkSymmetry(&symmetry, (leftmatchptr->c2wupout * leftmatchptr->c2wdownout), (middlematchptr->c2wupout * middlematchptr->c2wdownout), (rightmatchptr->c2wupout * rightmatchptr->c2wdownout));
        checkSymmetry(&symmetry, leftmatchptr->upin, middlematchptr->upin, rightmatchptr->upin);
        checkSymmetry(&symmetry, leftmatchptr->downin, middlematchptr->downin, rightmatchptr->downin);
        checkSymmetry(&symmetry, leftmatchptr->nbvupin, middlematchptr->nbvupin, rightmatchptr->nbvupin);
        checkSymmetry(&symmetry, leftmatchptr->nbsupin, middlematchptr->nbsupin, rightmatchptr->nbsupin);
        checkSymmetry(&symmetry, leftmatchptr->e2upin, middlematchptr->e2upin, rightmatchptr->e2upin);
        checkSymmetry(&symmetry, leftmatchptr->piupin, middlematchptr->piupin, rightmatchptr->piupin);
        checkSymmetry(&symmetry, (leftmatchptr->aupin * leftmatchptr->adownin), (middlematchptr->aupin * middlematchptr->adownin), (rightmatchptr->aupin * rightmatchptr->adownin));
        if ((symmetry >= minsymmetry) && (complexity <= maxcomplexity)) {
          initUses(&rightuses);
          addUses(&rightuses, &rightmatchptr->uses);
          initUses(&alluses);
          addUses(&alluses, &leftuses);
          addUses(&alluses, &middleuses);
          addUses(&alluses, &rightuses);
          clock_gettime(CLOCK_REALTIME, &starttime);
          precision=solvePolyforMasses(poly, leftinvexp, middleinvexp, rightinvexp, leftmatchptr, middlematchptr, rightmatchptr, &alluses, maxcomplexity);
#ifdef SHOWSTATUS
          clock_gettime(CLOCK_REALTIME, &endtime);
          elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
          printf("status, Solved  phase 2 formula  for masses, random input: %d, polyform: %s, mass mode: %d%d%d, progress: total (%ld/%ld) left (%d/%d) middle (%d/%d) right (%d/%d), precision: %.3e, (%6.4fs)\n", random_input_count, poly, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, combo, totalcombos, l+1, numleftmatches, m+1, nummiddlematches, r+1, numrightmatches, precision, elapsedtime);
          fflush(stdout);
#endif
        }
        rightmatchptr++;
      } // for r
      middlematchptr++;
    } // for m
    leftmatchptr++;
  } // for l
  free(leftmatches);
  free(middlematches);
  free(rightmatches);
  return(0);
}

int interesting(int range, double inld) {
  double rangelow;
  double rangehigh;
  double testld;
  int testint;

/*
  if (range == 2) {
    rangehigh=1.01;
    rangelow=0.99;
  if (range == 3) {
    rangehigh=1.001;
    rangelow=0.999;
  if (range == 4) {
    rangehigh=1.0001;
    rangelow=0.9999;
*/
  if (range == 5) {
    rangehigh=1.00001;
    rangelow=0.99999;
  } else if (range == 6) {
    rangehigh=1.000001;
    rangelow=0.999999;
/*
  } else if (range == 7) {
    rangehigh=1.0000001;
    rangelow=0.9999999;
  } else if (range == 8) {
    rangehigh=1.00000001;
    rangelow=0.99999999;
*/
  }

  // test if too big or small
  if ((inld < 0.620) || (inld > 16.1)) {  // 27=0.037 / 64=0.0156 // 9=0.1 // 16=.0620 // 32=0.0311
    return(0);
  }

  // test if int > 4 are not divisible by 2 or 3
  testint=(int)(inld + 0.5);
  if ((testint > 4) && ((testint % 2) != 0) && ((testint %3) != 0)) {
    return(0);
  }

  // test if close to int
  testld=fmodl(inld, 1.0); 
  if ((testld >= rangelow) && (testld <=rangehigh)) {
    return(1);
  }

  if  (inld < 1.0) {
    //  test if 1/int > 4 are divisible by 2 or 3
    testint=(int)((1.0 / inld) + 0.5);
    if ((testint > 4) && ((testint % 2) != 0) && ((testint %3) != 0)) {
      return(0);
    }

    // test if small number matches 1/int
    testld=fmodl((1.0 / inld), 1.0);
    if ((testld >= rangelow) && (testld <=rangehigh)) {
      return(1);
    }
  }

  return(0);
}

void initMultiplierArray(multipliers *mult, int *nummult) {
  // store pre-computed multiplier terms.   Store enough terms to accelerate cscanner but not so many as to 
  // break the l2 cache benefits

  struct timespec starttime;
  struct timespec endtime;
  double elapsedtime;

  int upin=1, downin=1;
  int piupin=0, pidownin=1;
  int aupin=0, adownin=1;
  int e2upin=0, e2downin=2;
  int nbvupin=0;
  int nbsupin=0;
  double updownin;
  double piin;
  double ain;
  double e2in;
  double nbv[27];
  double nbs[27];

  unsigned int u, v;
  int i;
  multipliers *m;

  printf("init, Initializing pre-computed multiplier accelerator array (a subset of multiplier terms that easily fits in on-chip cache)\n");
  clock_gettime(CLOCK_REALTIME, &starttime);

  // nball volume (index is term exponent and n-ball dimension)
  // nbv[n] = nbs[n-1]/n (n > 0)
  // nbv[n] = 2pi * nbv[n-2]/n (n > 1)
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

  // nball surface area (index is term exponent and n-sphere dimension and n-ball dimension - 1)
  // nbs[n] = 2pi * nbv[n-1] (n > 0)
  // nbs[n] = 2pi * nbs[n-2]/(n-1) (n > 1)
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

  m=mult;
  for (upin=1; upin<=32; upin++) {
    for (downin=1; downin<=32; downin++) {
      u=upin;
      v=downin;
      if ((upin <= 4) || ((downin == 1) && ((upin == 8) || (upin == 9) || (upin == 12) || (upin == 16) || (upin == 18) || (upin == 24) || (upin == 27) || (upin == 32)))) {
        if ((downin <= 4) || ((upin == 1) && ((downin == 8) || (downin == 9) || (downin == 12) || (downin == 16) || (downin == 18) || (downin == 24) || (downin == 27) || (downin == 32)))) { 
          if (gcd(u, v) == 1) {
            updownin=(double)upin / (double)downin;
            for (nbvupin = -1; nbvupin <=1; nbvupin++) {
              for (nbsupin = -1; nbsupin <=1; nbsupin++) {
                if (((nbvupin == 0) && (nbsupin == 0)) || ((nbvupin != 0) && (nbsupin == 0)) || ((nbvupin == 0) && (nbsupin != 0))) {

                  for (piupin = -2; piupin <=2; piupin++) {
/*
                    for (pidownin=1; pidownin <= 1; pidownin++) {
                      u=abs(piupin);
                      v=pidownin;
                      if (gcd(u, v) == 1) {
*/
                        piin=updownin * pow(M_PI, ((float)piupin / 1.0));
                        for (aupin=-2; aupin<=2; aupin++) {
                          for (adownin=1; adownin <= 2; adownin++) {
                            u=abs(aupin);
                            v=adownin;
                            if (gcd(u, v) == 1) {
                              ain=piin * pow(alpha_ref, ((float)aupin / 1.0));
                              for (e2upin=-1; e2upin<=1; e2upin++) {
                                for (e2downin=2; e2downin <= 2; e2downin++) {
                                  u=abs(e2upin);
                                  v=e2downin;
                                  if (gcd2(u, v) == 1) {
                                    e2in=ain * pow(2.0, ((float)e2upin / (float)e2downin));
                                    mult->upin=upin;
                                    mult->downin=downin;
                                    mult->piupin=piupin;
                                    mult->pidownin=pidownin;
                                    mult->aupin=aupin;
                                    mult->adownin=adownin;
                                    mult->e2upin=e2upin;
                                    mult->e2downin=e2downin;
                                    mult->nbvupin=nbvupin;
                                    mult->nbsupin=nbsupin;
                                    mult->multcomplexity=\
                                             (upin + downin)\
                                           + abs(nbvupin)\
                                           + abs(nbsupin)\
                                           + (abs(piupin) + pidownin)\
                                           + (abs(aupin) + adownin)\
                                           + (abs(e2upin) + e2downin);
                                    for (i=1; i<=27; i++) {
                                      if (i==1) {
                                        mult->mult[1]=e2in * pow(nbv[1], (float)nbvupin) * pow(nbs[1], (float)nbsupin);
                                      } else if (i==2) {
                                        mult->mult[2]=pow(e2in * pow(nbv[2], (float)nbvupin) * pow(nbs[2], (float)nbsupin), (1.0 / 2.0));
                                      } else if (i==3) {
                                        mult->mult[3]=pow(e2in * pow(nbv[3], (float)nbvupin) * pow(nbs[3], (float)nbsupin), (1.0 / 3.0));
                                      } else if (i==4) {
                                        mult->mult[4]=pow(e2in * pow(nbv[4], (float)nbvupin) * pow(nbs[4], (float)nbsupin), (1.0 / 4.0));
                                      } else if (i==5) {
                                        mult->mult[5]=pow(e2in * pow(nbv[5], (float)nbvupin) * pow(nbs[5], (float)nbsupin), (1.0 / 5.0));
                                      } else if (i==6) {
                                        mult->mult[6]=pow(e2in * pow(nbv[6], (float)nbvupin) * pow(nbs[6], (float)nbsupin), (1.0 / 6.0));
                                      } else if (i==7) {
                                        mult->mult[7]=pow(e2in * pow(nbv[7], (float)nbvupin) * pow(nbs[7], (float)nbsupin), (1.0 / 7.0));
                                      } else if (i==8) {
                                        mult->mult[8]=pow(e2in * pow(nbv[8], (float)nbvupin) * pow(nbs[8], (float)nbsupin), (1.0 / 8.0));
                                      } else if (i==9) {
                                        mult->mult[9]=pow(e2in * pow(nbv[9], (float)nbvupin) * pow(nbs[9], (float)nbsupin), (1.0 / 9.0));
                                      } else if (i==10) {
                                        mult->mult[10]=pow(e2in * pow(nbv[10], (float)nbvupin) * pow(nbs[10], (float)nbsupin), (1.0 / 10.0));
                                      } else if (i==11) {
                                        mult->mult[11]=pow(e2in * pow(nbv[11], (float)nbvupin) * pow(nbs[11], (float)nbsupin), (1.0 / 11.0));
                                      } else if (i==12) {
                                        mult->mult[12]=pow(e2in * pow(nbv[12], (float)nbvupin) * pow(nbs[12], (float)nbsupin), (1.0 / 12.0));
                                      } else if (i==13) {
                                        mult->mult[13]=pow(e2in * pow(nbv[13], (float)nbvupin) * pow(nbs[13], (float)nbsupin), (1.0 / 13.0));
                                      } else if (i==14) {
                                        mult->mult[14]=pow(e2in * pow(nbv[14], (float)nbvupin) * pow(nbs[14], (float)nbsupin), (1.0 / 14.0));
                                      } else if (i==15) {
                                        mult->mult[15]=pow(e2in * pow(nbv[15], (float)nbvupin) * pow(nbs[15], (float)nbsupin), (1.0 / 15.0));
                                      } else if (i==16) {
                                        mult->mult[16]=pow(e2in * pow(nbv[16], (float)nbvupin) * pow(nbs[16], (float)nbsupin), (1.0 / 16.0));
                                      } else if (i==17) {
                                        mult->mult[17]=pow(e2in * pow(nbv[17], (float)nbvupin) * pow(nbs[17], (float)nbsupin), (1.0 / 17.0));
                                      } else if (i==18) {
                                        mult->mult[18]=pow(e2in * pow(nbv[18], (float)nbvupin) * pow(nbs[18], (float)nbsupin), (1.0 / 18.0));
                                      } else if (i==19) {
                                        mult->mult[19]=pow(e2in * pow(nbv[19], (float)nbvupin) * pow(nbs[19], (float)nbsupin), (1.0 / 19.0));
                                      } else if (i==20) {
                                        mult->mult[20]=pow(e2in * pow(nbv[20], (float)nbvupin) * pow(nbs[20], (float)nbsupin), (1.0 / 20.0));
                                      } else if (i==21) {
                                        mult->mult[21]=pow(e2in * pow(nbv[21], (float)nbvupin) * pow(nbs[21], (float)nbsupin), (1.0 / 21.0));
                                      } else if (i==22) {
                                        mult->mult[22]=pow(e2in * pow(nbv[22], (float)nbvupin) * pow(nbs[22], (float)nbsupin), (1.0 / 22.0));
                                      } else if (i==23) {
                                        mult->mult[23]=pow(e2in * pow(nbv[23], (float)nbvupin) * pow(nbs[23], (float)nbsupin), (1.0 / 23.0));
                                      } else if (i==24) {
                                        mult->mult[24]=pow(e2in * pow(nbv[24], (float)nbvupin) * pow(nbs[24], (float)nbsupin), (1.0 / 24.0));
                                      } else if (i==25) {
                                        mult->mult[25]=pow(e2in * pow(nbv[25], (float)nbvupin) * pow(nbs[25], (float)nbsupin), (1.0 / 25.0));
                                      } else if (i==26) {
                                        mult->mult[26]=pow(e2in * pow(nbv[26], (float)nbvupin) * pow(nbs[26], (float)nbsupin), (1.0 / 26.0));
                                      }  // if i
                                    } // for i
                                    initUses(&mult->uses);
                                    if (aupin != 0) {
                                      mult->uses.alpha_em=1;
                                    }
                                    *nummult=*nummult+1;
                                    mult++;
                                  } // e2up gcd
                                } // e2downin
                              } // e2upin
                            } // gcd ain
                          } // adownin
                        } // aupin
/*
                      } // gcd piin
                    } // pidownin
*/
                  } // piupin
                } // not both nbvupin and nbsupin
              } // nbsupin
            } // nbvupin
          } // gcd updownin
        }  // downmod
      }  // upmod
    } // downin
  } // upin

  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);

  printf("init, Initialized %d pre-computed multipliers (%6.4fs)\n", *nummult, elapsedtime);
}

void cscanner(multipliers *multstart, int *nummult, matches **matchptr, int *nummatches, int *coffhit, int range, char *poly, int leftinvexp, int middleinvexp, int rightinvexp, double cleft, double cmiddle, double cright, random_input random_inputs) {
  /*************/
  /* Phase 1.5 */
  /*************/
  //  Each coefficient is multiplied by various numbers and the resulting value is tested to see if it is
  // close to an interesting integer or simple rational number.  The results are stored in the match table.
  int i;
  double multiplierout;
  multipliers *mult;
  unsigned int u, v;
  matches *match;
  int complexity;

  int massratio;
  double massterm;
  double leftmassterm;
  double middlemassterm;
  double rightmassterm;
  int upout=1, downout=1;
  int piupout=0, pidownout=1;
  int aupout=1, adownout=1;
  int e2upout=0, e2downout=2;
  int s2wupout=0, s2wdownout=1;
  int c2wupout=0, c2wdownout=1;
  double updownout;
  double piout;
  double aout;
  double e2out;
  double s2wout=1.0;
  double c2wout=1.0;

  double sin2w;
  double cos2w;

  sin2w=random_inputs.sin2w_sample;
  cos2w=1.0 - sin2w;

  match=*matchptr;
  // here we substitute the M/me mass ratio used in phase 1 with the actual test mass ratio
  for (massratio=0; massratio<=4; massratio++) {
    if (massratio == 0) {         
      massterm=random_inputs.mp_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mp, ");             
      fflush(stdout);               
#endif
    } else if (massratio == 1) {  
      massterm=v_ref/me_ref;
#ifdef SHOWSTATUS
      printf("M/v, ");              
      fflush(stdout);               
#endif
    } else if (massratio == 2) {  
      massterm=random_inputs.mz_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mz, ");             
      fflush(stdout);               
#endif
    } else if (massratio == 3) {  
      massterm=random_inputs.mw_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mw, ");             
      fflush(stdout);               
#endif
    } else if (massratio == 4) {  
      massterm=random_inputs.mh0_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mh0");             
      fflush(stdout);               
#endif
    }                             
    leftmassterm=pow(massterm, (1.0 / (double)leftinvexp));
    middlemassterm=pow(massterm, (1.0 / (double)middleinvexp));
    rightmassterm=pow(massterm, (1.0 / (double)rightinvexp));

    for (upout=1; upout <= 3; upout++) {
      for (downout=1; downout <=3; downout++) {
        u=upout;
        v=downout;
        if (gcd(u, v) == 1) {
          updownout=(double)upout / (double)downout;
          for (piupout=-2; piupout <= 2; piupout++) {
/*
            for (pidownout=1; pidownout <= 1; pidownout++) {
              u=abs(piupout);
              v=pidownout;
              if (gcd(u, v) == 1) {
*/
                piout=updownout * pow(M_PI, ((float)piupout / (float)pidownout));
                for (aupout=-2; aupout <= 2; aupout++) {
                  for (adownout=1; adownout <= 2; adownout++) {
                    u=abs(aupout);
                    v=adownout;
                    if (gcd(u, v) == 1) {
                      aout=piout * pow(alpha_ref, ((float)aupout / (float)adownout));
                      for (e2upout=-1; e2upout <= 1; e2upout++) {
                        e2out=aout * pow(2.0, ((float)e2upout / (float)e2downout));
// s2w and c2w are calculated separately since they contain significant uncertainty and are reconstructed separately in phase 2
#ifdef SIN2W
                        for (s2wupout=-1; s2wupout <= 1; s2wupout++) {
                          for (s2wdownout=1; s2wdownout <= 2; s2wdownout++) {
                            u=abs(s2wupout);
                            v=s2wdownout;
                            if ((gcd2(u, v) == 1) && ((v % 2) == 0)) {
                              s2wout=pow(sin2w, ((float)s2wupout / (float)s2wdownout));
                              for (c2wupout=-1; c2wupout <= 1; c2wupout++) {
                                for (c2wdownout=1; c2wdownout <= 2; c2wdownout++) {
                                  u=abs(c2wupout);
                                  v=c2wdownout;
                                  if ((gcd2(u, v) == 1) && ((v % 2) == 0)) {
                                    c2wout=pow(cos2w, ((float)c2wupout / (float)c2wdownout));
#endif
                                    mult=multstart;
                                    for (i=0; i<*nummult; i++) { 
                                      // test multiplier against left coefficient
                                      if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == leftinvexp) || (s2wdownout == (leftinvexp * 2))) {
                                        if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == leftinvexp) || (c2wdownout == (leftinvexp * 2))) {
                                          multiplierout=cleft * leftmassterm * s2wout * c2wout * e2out * mult->mult[abs(leftinvexp)];
                                          if (interesting(range, multiplierout)) {
                                            complexity=mult->multcomplexity\
                                                  + (upout + downout)\
                                                  + (abs(piupout) + pidownout)\
                                                  + (abs(aupout) + adownout)\
                                                  + (abs(e2upout) + e2downout);
/*
                                                  + (abs(s2wupout) + s2wdownout)\
                                                  + (abs(c2wupout) + c2wdownout);
*/
                                            coffhit[0]=leftinvexp;
                                            match->invexp=leftinvexp;
                                            match->massratio=massratio;
                                            match->upin=mult->upin;
                                            match->downin=mult->downin;
                                            match->piupin=mult->piupin;
                                            match->pidownin=mult->pidownin;
                                            match->aupin=mult->aupin;
                                            match->adownin=mult->adownin;
                                            match->e2upin=mult->e2upin;
                                            match->e2downin=mult->e2downin;
                                            match->nbvupin=mult->nbvupin;
                                            match->nbsupin=mult->nbsupin;
                                            match->upout=upout;
                                            match->downout=downout;
                                            match->piupout=piupout;
                                            match->pidownout=pidownout;
                                            match->aupout=aupout;
                                            match->adownout=adownout;
                                            match->e2upout=e2upout;
                                            match->e2downout=e2downout;
                                            match->s2wupout=s2wupout;
                                            match->s2wdownout=s2wdownout;
                                            match->c2wupout=c2wupout;
                                            match->c2wdownout=c2wdownout;
                                            match->matchcomplexity=complexity;
                                            match->match=multiplierout;
                                            initUses(&match->uses);
                                            addUses(&match->uses, &mult->uses);
                                            if (massratio == 0) {
                                              match->uses.G=1;
                                            } else if (massratio == 1) {
                                              match->uses.v=1;
                                            } else if (massratio == 2) {
                                              match->uses.mz=1;
                                            } else if (massratio == 3) {
                                              match->uses.mw=1;
                                            } else if (massratio == 4) {
                                              match->uses.mh0=1;
                                            }
                                            if (aupout != 0) {
                                              match->uses.alpha_em=1;
                                            }
                                            if ((s2wupout != 0) || (c2wupout != 0)) {
                                              match->uses.sin2w=1;
                                            }
                                            match->matchmult=e2out * mult->mult[abs(leftinvexp)];
                                            *nummatches=*nummatches+1;
                                            *matchptr=*matchptr+1;
                                            match=*matchptr;
                                          }  // if interesting left
                                        } // sanity check c2w left
                                      } // sanity check s2w left

                                      // test multiplier against middle coefficient
                                      if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == middleinvexp) || (s2wdownout == (middleinvexp * 2))) {
                                        if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == middleinvexp) || (c2wdownout == (middleinvexp * 2))) {
                                          multiplierout=cmiddle * middlemassterm * s2wout * c2wout * e2out * mult->mult[abs(middleinvexp)];
                                          if (interesting(range, multiplierout)) { 
                                            complexity=mult->multcomplexity\
                                                  + (upout + downout)\
                                                  + (abs(piupout) + pidownout)\
                                                  + (abs(aupout) + adownout)\
                                                  + (abs(e2upout) + e2downout);
/*
                                                  + (abs(s2wupout) + s2wdownout)\
                                                  + (abs(c2wupout) + c2wdownout);
*/
                                            coffhit[1]=middleinvexp;
                                            match->invexp=middleinvexp;
                                            match->massratio=massratio;
                                            match->upin=mult->upin;
                                            match->downin=mult->downin;
                                            match->piupin=mult->piupin;
                                            match->pidownin=mult->pidownin;
                                            match->aupin=mult->aupin;
                                            match->adownin=mult->adownin;
                                            match->e2upin=mult->e2upin;
                                            match->e2downin=mult->e2downin;
                                            match->nbvupin=mult->nbvupin;
                                            match->nbsupin=mult->nbsupin;
                                            match->upout=upout;
                                            match->downout=downout;
                                            match->piupout=piupout;
                                            match->pidownout=pidownout;
                                            match->aupout=aupout;
                                            match->adownout=adownout;
                                            match->e2upout=e2upout;
                                            match->e2downout=e2downout;
                                            match->s2wupout=s2wupout;
                                            match->s2wdownout=s2wdownout;
                                            match->c2wupout=c2wupout;
                                            match->c2wdownout=c2wdownout;
                                            match->matchcomplexity=complexity;
                                            match->match=multiplierout;
                                            initUses(&match->uses);
                                            addUses(&match->uses, &mult->uses);
                                            if (massratio == 0) {
                                              match->uses.G=1;
                                            } else if (massratio == 1) {
                                              match->uses.v=1;
                                            } else if (massratio == 2) {
                                              match->uses.mz=1;
                                            } else if (massratio == 3) {
                                              match->uses.mw=1;
                                            } else if (massratio == 4) {
                                              match->uses.mh0=1;
                                            }
                                            if (aupout != 0) {
                                              match->uses.alpha_em=1;
                                            }
                                            if ((s2wupout != 0) || (c2wupout != 0)) {
                                              match->uses.sin2w=1;
                                            }
                                            match->matchmult=e2out * mult->mult[abs(middleinvexp)];
                                            *nummatches=*nummatches+1;
                                            *matchptr=*matchptr+1;
                                            match=*matchptr;
                                          }  // if interesting middle
                                        } // sanity check c2w middle
                                      } // sanity check s2w middle

                                      // test multiplier against right coefficient
                                      if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == rightinvexp) || (s2wdownout == (rightinvexp * 2))) {
                                        if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == rightinvexp) || (c2wdownout == (rightinvexp * 2))) {
                                          multiplierout=cright * rightmassterm * s2wout * c2wout * e2out * mult->mult[abs(rightinvexp)];
                                          if (interesting(range, multiplierout)) { 
                                            complexity=mult->multcomplexity\
                                                  + (upout + downout)\
                                                  + (abs(piupout) + pidownout)\
                                                  + (abs(aupout) + adownout)\
                                                  + (abs(e2upout) + e2downout);
/*
                                                  + (abs(s2wupout) + s2wdownout)\
                                                  + (abs(c2wupout) + c2wdownout);
*/
                                            coffhit[2]=rightinvexp;
                                            match->invexp=rightinvexp;
                                            match->massratio=massratio;
                                            match->upin=mult->upin;
                                            match->downin=mult->downin;
                                            match->piupin=mult->piupin;
                                            match->pidownin=mult->pidownin;
                                            match->aupin=mult->aupin;
                                            match->adownin=mult->adownin;
                                            match->e2upin=mult->e2upin;
                                            match->e2downin=mult->e2downin;
                                            match->nbvupin=mult->nbvupin;
                                            match->nbsupin=mult->nbsupin;
                                            match->upout=upout;
                                            match->downout=downout;
                                            match->piupout=piupout;
                                            match->pidownout=pidownout;
                                            match->aupout=aupout;
                                            match->adownout=adownout;
                                            match->e2upout=e2upout;
                                            match->e2downout=e2downout;
                                            match->s2wupout=s2wupout;
                                            match->s2wdownout=s2wdownout;
                                            match->c2wupout=c2wupout;
                                            match->c2wdownout=c2wdownout;
                                            match->matchcomplexity=complexity;
                                            match->match=multiplierout;
                                            initUses(&match->uses);
                                            addUses(&match->uses, &mult->uses);
                                            if (massratio == 0) {
                                              match->uses.G=1;
                                            } else if (massratio == 1) {
                                              match->uses.v=1;
                                            } else if (massratio == 2) {
                                              match->uses.mz=1;
                                            } else if (massratio == 3) {
                                              match->uses.mw=1;
                                            } else if (massratio == 4) {
                                              match->uses.mh0=1;
                                            }
                                            if (aupout != 0) {
                                              match->uses.alpha_em=1;
                                            }
                                            if ((s2wupout != 0) || (c2wupout != 0)) {
                                              match->uses.sin2w=1;
                                            }
                                            match->matchmult=e2out * mult->mult[abs(rightinvexp)];
                                            *nummatches=*nummatches+1;
                                            *matchptr=*matchptr+1;
                                            match=*matchptr;
                                          }  // if interesting right
                                        } // sanity check c2w right
                                      } // sanity check s2w right
                                      mult++;
                                    }  // for i
#ifdef SIN2W
                                  } // gcd c2wout
                                } // c2wdownout
                              } // c2wupout
                            } // gcd s2wout
                          } // s2wdownout
                        } // s2wupout
#endif
                      } // e2upout
                    } // gcd aout
                  } // adownout
                } // aupout
/*
              } // gcd piout
            } // pidownout
*/
          } // piupout
        } // gcd updownout
      } // downout
    } // upout
  } // for massratio
#ifdef SHOWSTATUS
  printf(".\n");
#endif
}

int solvePolyforCoefficients(multipliers *multstart, int *nummult, matches **matchptr, int *nummatches, int *coffhit, int random_input_count, random_input random_inputs, char *poly, int leftinvexp, int middleinvexp, int rightinvexp, int range) {
  /***********/
  /* Phase 1 */
  /***********/
  // solve a three term polynomial-like formula for the unknown coefficients given known roots (particle masses)
  int i;
  long int samples=0;
  struct timespec starttime;
  struct timespec endtime;
  double elapsedtime;
  int nummatchstart;
  int nummatchend;
  int goodc;
  int pml, pmr, plr;

  clock_gettime(CLOCK_REALTIME, &starttime);

  //  mc test vars
  double r;
  double e_test_1;
  double e_test_2;
  double e_test_3;
  double e_test=0;
  double u_test_1;
  double u_test_2;
  double u_test_3;
  double u_test=0;
  double t_test_1;
  double t_test_2;
  double t_test_3;
  double t_test=0;
  double worst_test;
  double rangefactor;
  double massleft_e=0;
  double massmiddle_e=0;
  double massright_e=0;
  double massleft_u=0;
  double massmiddle_u=0;
  double massright_u=0;
  double massleft_t=0;
  double massmiddle_t=0;
  double massright_t=0;
  double precision;
  double rangemultiplier[6];

  // these tunings affect speed and reliability, adjust with extreme care
  double testratio=25.0;              // acceptable ratios of e_test/u_test/t_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
  int ratiograceperiod=15;            // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
  int stalledlimit=500000;            // most polyforms can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
  double defaultrangemultiplier=5.0;  // lowest practical range multiplier, fastest for most polyforms
  double stalledrangemultiplier=17.0; // this value works better for slow to solve polyforms and fast polyforms that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve polyforms this will continuously trigger
  int slowcheckpoint=1000000;         // progress point to check on slow processes
  double stuckprecision=1.0E-2;       // if precision is not past this level by slowcheckpoint, try resetting

  // mc outputs
  int progress[6];
  int stalled[6];
  int ordering;
  int best_ordering;
  double best_precision_last;
  double precision_last[6];
  double cleft[6];
  double cleft_center[6];
  double cleft_range[6];
  double cleft_range_new;
  double cmiddle[6];
  double cmiddle_center[6];
  double cmiddle_range[6];
  double cmiddle_range_new;
  double cright[6];
  double cright_center[6];
  double cright_range[6];
  double cright_range_new;

  // starting v3.39 denominator is always electron mass for this step, actual mass ratio factors are now added in cscanner
  massleft_e=  1.0;
  massleft_u=  pow(((double)mu_ref                   / me_ref), (1.0 / (double)leftinvexp));
  massleft_t=  pow(((double)random_inputs.tau_sample / me_ref), (1.0 / (double)leftinvexp));
  massmiddle_e=1.0;
  massmiddle_u=pow(((double)mu_ref                   / me_ref), (1.0 / (double)middleinvexp));
  massmiddle_t=pow(((double)random_inputs.tau_sample / me_ref), (1.0 / (double)middleinvexp));
  massright_e= 1.0;
  massright_u= pow(((double)mu_ref                   / me_ref), (1.0 / (double)rightinvexp));
  massright_t= pow(((double)random_inputs.tau_sample / me_ref), (1.0 / (double)rightinvexp));

#ifdef SHOWSTATUS
  printf("status, Solving phase 1 formula for coefficients, polyform: %s\n", poly);
  fflush(stdout);
#endif

  //  solve formula for coefficients
  best_precision_last=1.0E99;
  while (best_precision_last > 1.0E-11) {
    //  init outputs
    for (ordering=0; ordering<=5; ordering++) {
      cleft[ordering]=1.0E3;
      cleft_center[ordering]=1.0E3;
      cleft_range[ordering]=0.5E3;
      cmiddle[ordering]=1.0E3;
      cmiddle_center[ordering]=1.0E3;
      cmiddle_range[ordering]=0.5E3;
      cright[ordering]=1.0E3;
      cright_center[ordering]=1.0E3;
      cright_range[ordering]=0.5E3;
      precision_last[ordering]=1.0E99;
      progress[ordering]=0;
      stalled[ordering]=0;
      rangemultiplier[ordering]=defaultrangemultiplier;
    }
    best_precision_last=1.0E99;
    best_ordering=-1;
    i=0;
    for (samples=0; best_precision_last > 1.0E-11; samples++) {
      for (ordering=0; ordering<=5; ordering++) {
        if ((samples > 1) && ((samples % slowcheckpoint) == 0)) { // check on slow processes
#ifdef DEBUG10
          if ((samples % 10000000) == 0) { // rate limit periodic debug prints
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);

            printf("debug, polyform: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9e, best_ordering: %d, precision_last: %.9e, rangefactor: %.9e, i: %d, e_test: %.9e, u_test: %.9e, t_test: %.9e, cleft: %.9e, cmiddle: %.9e, cright: %.9e, cleftrange: %.9e, cmiddlerange: %.9e, crightrange: %.9e\n", poly, samples, elapsedtime, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], rangefactor, i, e_test, u_test, t_test, cleft_center[ordering], cmiddle_center[ordering], cright_center[ordering], cleft_range[ordering], cmiddle_range[ordering], cright_range[ordering]);
            fflush(stdout);
          }
#endif
          if ((progress[ordering] == ratiograceperiod) || (precision_last[ordering] > stuckprecision)) { // it's stuck, try resetting
#ifdef DEBUG10
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
            printf("debug, polyform: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, resetting\n", poly, samples, elapsedtime, ordering, progress[ordering]);
            fflush(stdout);
#endif
            cleft[ordering]=1.0E3;
            cleft_center[ordering]=1.0E3;
            cleft_range[ordering]=0.5E3;
            cmiddle[ordering]=1.0E3;
            cmiddle_center[ordering]=1.0E3;
            cmiddle_range[ordering]=0.5E3;
            cright[ordering]=1.0E3;
            cright_center[ordering]=1.0E3;
            cright_range[ordering]=0.5E3;
            precision_last[ordering]=1.0E99;
            progress[ordering]=0;
            stalled[ordering]=0;
            rangemultiplier[ordering]=defaultrangemultiplier;
          }
        }
        if (((best_precision_last > 1.0E-7) || (ordering == best_ordering)) && (cleft_center[ordering] > 1.0E-9) && (cmiddle_center[ordering] > 1.0E-9) && (cright_center[ordering] > 1.0E-9)) { // skip other ordings if one is far enough along and skip oderings where any coefficient has gone out of bounds
        if (ordering == 0) {
          pml=1;
          pmr=1;
          plr=1;
        } else if (ordering == 1) {
          pml=1;
          pmr=1;
          plr=0;
        } else if (ordering == 2) {
          pml=1;
          pmr=0;
          plr=0;
        } else if (ordering == 3) {
          pml=0;
          pmr=1;
          plr=1;
        } else if (ordering == 4) {
          pml=0;
          pmr=0;
          plr=1;
        } else if (ordering == 5) {
          pml=0;
          pmr=0;
          plr=0;
        } // end if ordering
        if (stalled[ordering] == stalledlimit) {
          rangemultiplier[ordering]=stalledrangemultiplier; // may be a slow solution, try bigger multiplier
          //stalled[ordering]=0;
#ifdef DEBUG10
          clock_gettime(CLOCK_REALTIME, &endtime);
          elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
          printf("debug, polyform: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, stalled\n", poly, samples, elapsedtime, ordering, progress[ordering]);
          fflush(stdout);
#endif
        }
        // generate random coefficients within limits
        stalled[ordering]++;
        i=0;
        goodc=0;
        while ((i < 200) && (goodc == 0)) {
          i++;
          cleft[ordering]=0.0;
          cmiddle[ordering]=0.0;
          cright[ordering]=0.0;
          while (cleft[ordering] <= 0.0) {
            r=(double)drand48();
            cleft[ordering]=((cleft_center[ordering] - cleft_range[ordering]) + (r * 2.0 * cleft_range[ordering]));
          }
          while (cmiddle[ordering] <= 0.0) {
            r=(double)drand48();
            cmiddle[ordering]=((cmiddle_center[ordering] - cmiddle_range[ordering]) + (r * 2.0 * cmiddle_range[ordering]));
          }                 
          while (cright[ordering] <= 0.0) {
            r=(double)drand48();
            cright[ordering]=((cright_center[ordering] - cright_range[ordering]) + (r * 2.0 * cright_range[ordering]));
          }
          // verify coefficient relationships
          goodc=1;
          if (pml == 1) {
            if (cmiddle[ordering] < cleft[ordering]) {
              goodc=0;
            }
          } else {
            if (cmiddle[ordering] > cleft[ordering]) {
              goodc=0;
            }
          }
          if (pmr == 1) {
            if (cmiddle[ordering] < cright[ordering]) {
              goodc=0;
            }
          } else {
            if (cmiddle[ordering] > cright[ordering]) {
              goodc=0;
            }
          }
          if (plr == 1) {
            if (cleft[ordering] < cright[ordering]) {
              goodc=0;
            }
          } else {
            if (cleft[ordering] > cright[ordering]) {
              goodc=0;
            }
          } // end if pml
        }
#ifdef DEBUG11
        if (i > 100) { 
          printf("debug, i: %d\n", i);
        }
#endif
        if (goodc == 1) {
          e_test_1=cleft[ordering] * massleft_e;
          e_test_2=cmiddle[ordering] * massmiddle_e;
          e_test_3=cright[ordering] * massright_e;
  
          u_test_1=cleft[ordering] * massleft_u;
          u_test_2=cmiddle[ordering] * massmiddle_u;
          u_test_3=cright[ordering] * massright_u;

          t_test_1=cleft[ordering] * massleft_t;
          t_test_2=cmiddle[ordering] * massmiddle_t;
          t_test_3=cright[ordering] * massright_t;

          e_test=e_test_1 - e_test_2 + e_test_3 - 1.0;
          u_test=u_test_1 - u_test_2 + u_test_3 - 1.0;
          if ((progress[ordering] < ratiograceperiod) || (((fabs(e_test) / fabs(u_test)) < testratio) && ((fabs(u_test) / fabs(e_test)) < testratio))) {
            t_test=t_test_1 - t_test_2 + t_test_3 - 1.0;
              if ((progress[ordering] < ratiograceperiod) || (((fabs(e_test) / fabs(t_test)) < testratio) && ((fabs(t_test) / fabs(e_test)) < testratio) &&\
                                                 ((fabs(u_test) / fabs(t_test)) < testratio) && ((fabs(t_test) / fabs(u_test)) < testratio))) {
#ifdef DEBUG12
              printf("debug, etest: %.3e, utest: %.3e, ttest: %.3e, e_test_1: %.3e, e_test_2: %.3e, e_test_3: %.3e, u_test_1: %.3e, u_test_2: %.3e, u_test_3: %.3e, t_test_1: %.3e, t_test_2: %.3e, t_test_3: %.3e\n", e_test, u_test, t_test, e_test_1, e_test_2, e_test_3, u_test_1, u_test_2, u_test_3, t_test_1, t_test_2, t_test_3);
              fflush(stdout);
#endif

              precision=(fabs(e_test) + fabs(u_test) + fabs(t_test));
              if (precision < precision_last[ordering]) {
                progress[ordering]++;
                stalled[ordering]=0;
                //rangemultiplier[ordering]=defaultrangemultiplier;
                precision_last[ordering]=precision;
                if (precision_last[ordering] < best_precision_last) {
                  best_precision_last=precision;
                  best_ordering=ordering;
                }
                if (fabs(e_test) > fabs(u_test)) {
                  worst_test=fabs(e_test);
                } else {
                  worst_test=fabs(u_test);
                }
                if (fabs(t_test) > worst_test) {
                  worst_test=fabs(t_test);
                }
                rangefactor=worst_test * rangemultiplier[ordering];
                if (rangefactor > 1.0) {
                  cleft_range[ordering]=cleft[ordering] / 2.0;
                  cmiddle_range[ordering]=cmiddle[ordering] / 2.0;
                  cright_range[ordering]=cright[ordering] / 2.0;
                } else {
/*
                  // option 1 - comparable to option 2, maybe slightly slower
                  cleft_range[ordering]=cleft[ordering] * rangefactor;
                  cmiddle_range[ordering]=cmiddle[ordering] * rangefactor;
                  cright_range[ordering]=cright[ordering] * rangefactor;
*/
                  // option 2 - best for now
                  cleft_range_new=cleft[ordering] * rangefactor;
                  cleft_range[ordering]=((cleft_range[ordering] + cleft_range_new + cleft_range_new) / 3.0);
                  cmiddle_range_new=cmiddle[ordering] * rangefactor;
                  cmiddle_range[ordering]=((cmiddle_range[ordering] + cmiddle_range_new + cmiddle_range_new) / 3.0);
                  cright_range_new=cright[ordering] * rangefactor;
                  cright_range[ordering]=((cright_range[ordering] + cright_range_new + cright_range_new) / 3.0);
/*
                  // option 3 - consistently slower than option 2 and sometimes dramatically slows some normally fast polyforms
                  cleft_range_new=cleft[ordering] * rangefactor;
                  cleft_range[ordering]=((cleft_range[ordering] + cleft_range_new) / 2.0);
                  cmiddle_range_new=cmiddle[ordering] * rangefactor;
                  cmiddle_range[ordering]=((cmiddle_range[ordering] + cmiddle_range_new) / 2.0);
                  cright_range_new=cright[ordering] * rangefactor;
                  cright_range[ordering]=((cright_range[ordering] + cright_range_new) / 2.0);
*/
/*
                  // option 4 - comparable to option 2 but sometimes slows normaly fast polyforms
                  cleft_range_new=(fabs(cleft[ordering] - cleft_center[ordering]) * 2.20);
                  cleft_range[ordering]=((cleft_range[ordering] + cleft_range_new + (cleft[ordering] * rangefactor)) / 3.0);
                  cmiddle_range_new=(fabs(cmiddle[ordering] - cmiddle_center[ordering]) * 2.20);
                  cmiddle_range[ordering]=((cmiddle_range[ordering] + cmiddle_range_new + (cmiddle[ordering] * rangefactor)) / 3.0);
                  cright_range_new=(fabs(cright[ordering] - cright_center[ordering]) * 2.20);
                  cright_range[ordering]=((cright_range[ordering] + cright_range_new + (cright[ordering] * rangefactor)) / 3.0);
*/
                }
                cleft_center[ordering]=cleft[ordering];
                cmiddle_center[ordering]=cmiddle[ordering];
                cright_center[ordering]=cright[ordering];

#ifdef DEBUG11
                clock_gettime(CLOCK_REALTIME, &endtime);
                elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                printf("debug, polyform: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9e, best_ordering: %d, precision_last: %.9e, rangefactor: %.9e, i: %d, e_test: %.9e, u_test: %.9e, t_test: %.9e, cleft: %.9e, cmiddle: %.9e, cright: %.9e, cleftrange: %.9e, cmiddlerange: %.9e, crightrange: %.9e\n", poly, samples, elapsedtime, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], rangefactor, i, e_test, u_test, t_test, cleft[ordering], cmiddle[ordering], cright[ordering], cleft_range[ordering], cmiddle_range[ordering], cright_range[ordering]);
                fflush(stdout);
#endif
              } // end if precision
            } // end ttest ratios
          } // end utest ratios
        }  // end if goodc
      } // end if best
      } // end for ordering
    }  // end for samples
  }  // end while precision_last

  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);

#ifdef SHOWSTATUS
  printf("status, Solved  phase 1 formula for coefficients, random input: %i, polyform:  %s, Tau: %.9e, samples: %ld, ordering: %d, precision: %.3e (%6.4fs)\n", random_input_count, poly, random_inputs.tau_sample, samples, best_ordering, precision, elapsedtime);
  fflush(stdout);
#endif
#ifdef SHOWSTATUS
  printf("status, +-----------------+------------+------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
  printf("status, |    tau mass     |    poly    | mass ratio | polynomial term |   coefficient   | c * polyterm(e) | c * polyterm(u) | c * polyterm(t) |\n");
  printf("status, +-----------------+------------+------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
  printf("status, | %.9e | %10s |    M/me    |     +left       | %.9e | %.9e | %.9e | %.9e |\n", random_inputs.tau_sample, poly, cleft_center[best_ordering], e_test_1, u_test_1, t_test_1);
  printf("status, | %.9e | %10s |    M/me    |     -middle     | %.9e | %.9e | %.9e | %.9e |\n", random_inputs.tau_sample, poly, cmiddle_center[best_ordering], e_test_2, u_test_2, t_test_2);
  printf("status, | %.9e | %10s |    M/me    |     +right      | %.9e | %.9e | %.9e | %.9e |\n", random_inputs.tau_sample, poly, cright_center[best_ordering], e_test_3, u_test_3, t_test_3);
  printf("status, +-----------------+------------+------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
  fflush(stdout);
#endif

  // for debugging very long phase 1 solutions
  //exit(0);

  clock_gettime(CLOCK_REALTIME, &starttime);
  nummatchstart=*nummatches;

#ifdef SHOWSTATUS
  printf("status, Scanning for coefficient multipliers that match interesting integer or simple rational numbers and mass ratio: ");
  fflush(stdout);
#endif
  cscanner(multstart, nummult, matchptr, nummatches, coffhit, range, poly, leftinvexp, middleinvexp, rightinvexp, cleft_center[best_ordering], cmiddle_center[best_ordering], cright_center[best_ordering], random_inputs);
  nummatchend=*nummatches;
  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
#ifdef SHOWSTATUS
  printf("status, Found %d interesting coefficient multipliers (%6.4fs)\n", (nummatchend-nummatchstart), elapsedtime);
  fflush(stdout);
#endif

  return(0);
}

int main(int argc, char **argv) {
  struct timespec t;
  long seed;
  double testrand;
  double r;
  int nummatches;
  int i;
  int coffhit[3];
  int numcoff;
  long exseed;
  long seedsec;
  long seedus;
  int nummult;
  multipliers *multstart;
  matches *matchstart;
  matches *matchptr;
  int leftinvexp;
  int middleinvexp;
  int rightinvexp;
  int numdimensions;
  int invexp1;
  int invexp2;
  int invexp3;
  random_input random_inputs;
  char poly[20];
  int maxcomplexity;
  int minsymmetry;

  int random_input_count=0;
  int range;

  if (argc == 6) {
    exseed=atol(argv[1]);
    numdimensions=atoi(argv[2]);
    range=atoi(argv[3]);
    minsymmetry=atoi(argv[4]);
    maxcomplexity=atoi(argv[5]);
  } else {
    printf("\n\
lepton version %s\n", version);
printf("\n\
usage: polylepton <seed> <exponentlimit> <phase1filter> <minsymmetry> <maxcomplexity>\n\
\n\
 example: polylepton 0 12 5 70 65 (no external seed, exponents from -1/12 to +1/12, phase1 filter of 1x10^-5,\n\
                             minimum phase 2 symmetry of 70, maximum phase 2 complexity of 65. (recommended values).\n\
\n\
 example: polylepton 1000000 26 6 30 150 (with an external seed, exponents from -1/26 to +1/26, phase1 fitler of 1x10^-6,\n\
                                       minimum phase 2 symmetry of 30, maximum phase 2 complexity of 150).\n\
\n\
 *** Warning: very rough work in progress ***\n\
\n\
 This program searches for polynomial-like formulas (polyforms) that might generate the charged lepton masses from simple coefficients.\n\
 These three term (plus constant of -1) formulas use positive and negative rational exponents less than 1 and are designed to give three positive real roots\n\
 representing the charged lepton masses. Unlike integer exponents, rational exponents <= +/-1 generate positive-only roots with a spectrum similar to that\n\
 of the observed masses using relatively small coefficients on each term.  To ensure each term is dimensionless, the charged lepton masses are combined with\n\
 a reference mass inside each term as a mass ratio.   Teh reference masses currently supported are the Planck mass \"mp\", Higgs vacuum expectation value \"v\",\n\
 Z boson mass \"mz\", W boson mass \"mw\" and Higgs boson mass \"mh0\".\n\
\n\
 Processing is broken into two phases.   Phase 1 starts with the known charged lepton masses, solves a proposed polyform for the three coefficients and searches\n\
 for interesting multipliers to those coefficients. Phase 2 reverses the process by converting the close multiplier matches from phase 1 to exact formulas\n\
 and solving them for the charged lepton masses and other outputs.  These results are then compared to their experimental values and uncertainties.\n\
\n\
 Phase 1:\n\
   A random combination of three exponents are selected.   The three charged lepton masses are used as inputs for each term with the Tau mass being varied\n\
   randomly within it's experimental uncertainty range each time phase 1 is run.  The electron mass is temporarily used as the reference mass in each term.\n\
\n\
   This formula is then solved for the coefficients for each of the three exponent terms (left, middle, right).  Each coefficient is then multiplied by a\n\
   series of test multipliers (including the actual reference masses) to see if the result is close to a relatively small integer rational number. The phase 1\n\
   limit determines which results are close enough to pass to phase 2.\n\
\n\
 Phase 2:\n\
   formulas are constructed with the same exponents and various interesting multipliers found in phase 1 and then solved for the three charged lepton masses and\n\
   other variables.  In phase 2 each exponent term can have a different mass ratio, as long as the product of the mass ratio and multiplier\n\
   generates the correct individual polynomial term.\n\
\n\
   Before solving each proposed formula, all of the experimentally known variables are ranked by relative standard uncertainty and the three with the highest\n\
   uncertainty are used as outputs (solved for) with the rest used as inputs. The resulting outputs are then checked against their experimental ranges.\n\
\n\
   Symmetry and complexity scores are also assigned to each result based on the numerical values in the multiplier, with higher symmetry and lower complexity\n\
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
    return(1);
  }

  // init pseudorandom number generator from external seed and clock
  clock_gettime(CLOCK_REALTIME, &t);
  seedsec=(t.tv_sec % 1000000000);
  seedus=(t.tv_nsec / 1000); 
  seed=exseed ^ (seedsec + seedus);
  srand48(seed);
  testrand=drand48();
  printf("init, version: %s, external seed: %ld, seedsec: %ld, seedus: %ld, seed: %ld, firstrand: %.9e\n", version, exseed, seedsec, seedus, seed, testrand);

  // init matches array
  matchstart = (matches *)malloc(5000000 * sizeof(matches));

  // init mult array
  nummult=0;
  multstart=(multipliers *)malloc(1000000 * sizeof(multipliers));
  initMultiplierArray(multstart, &nummult);

  while (1) {
    random_input_count++;
    // generate random tau mass
    r=drand48();
    random_inputs.tau_sample=((tau_ref - tau_ref_error) + (r * 2.0 * tau_ref_error));
    // use static reference tau mass
    //random_inputs.tau_sample=tau_ref_center;

    // generate random G
    r=drand48();
    random_inputs.G_sample=((G_ref - G_ref_error) + (r * 2.0 * G_ref_error));
    // use static reference G
    //G_sample=G_ref;
    random_inputs.mp_sample=kg_to_ev * sqrt(hbar_ref * c_ref / random_inputs.G_sample);

    // generate random mz
    r=drand48();
    random_inputs.mz_sample=((mz_ref - mz_ref_error) + (r * 2.0 * mz_ref_error));
    // use static reference mz
    //random_inputs.mz_sample=mz_ref;

    // generate random mw
    r=drand48();
    random_inputs.mw_sample=((mw_ref - mw_ref_error) + (r * 2.0 * mw_ref_error));
    // use static reference mw 
    //random_inputs.mw_sample=mw_ref;

    // generate random mh0
    r=drand48();
    random_inputs.mh0_sample=((mh0_ref - mh0_ref_error) + (r * 2.0 * mh0_ref_error));
    // use static reference mh0
    //random_inputs.mh0_sample=mh0_ref;

    // generate random sin2w
    r=drand48();
    random_inputs.sin2w_sample=((sin2w_ref - sin2w_ref_error) + (r * 2.0 * sin2w_ref_error));
    // use static reference sin2w
    //random_inputs.sin2w_sample=sin2w_ref;

#ifdef NEGATIVEEXP
    // allow negative
    // select random polyform
    r=drand48();
    invexp1=(int)(r * 2 * ((double)numdimensions + 0.5)) - numdimensions;
    while (invexp1 == 0) {
      r=drand48();
      invexp1=(int)(r * 2 * ((double)numdimensions + 0.5)) - numdimensions;
    }
    r=drand48();
    invexp2=(int)(r * 2 * ((double)numdimensions + 0.5)) - numdimensions;
    while ((invexp2 == 0) || (invexp2 == invexp1)) {
      r=drand48();
      invexp2=(int)(r * 2 * ((double)numdimensions + 0.5)) - numdimensions;
    }
    r=drand48();
    invexp3=(int)(r * 2 * ((double)numdimensions + 0.5)) - numdimensions;
    while ((invexp3 == 0) || (invexp3 == invexp1) || (invexp3 == invexp2)) {
      r=drand48();
      invexp3=(int)(r * 2 * ((double)numdimensions + 0.5)) - numdimensions;
    }
#else
    // only non-negative
    // select random polyform
    r=drand48();
    invexp1=(int)(r * ((double)numdimensions - 0.5)) + 1;
    while (invexp1 == 0) {
      r=drand48();
      invexp1=(int)(r * ((double)numdimensions - 0.5)) + 1;
    }
    r=drand48();
    invexp2=(int)(r * ((double)numdimensions - 0.5)) + 1;
    while ((invexp2 == 0) || (invexp2 == invexp1)) {
      r=drand48();
      invexp2=(int)(r * ((double)numdimensions - 0.5)) + 1;
    }
    r=drand48();
    invexp3=(int)(r * ((double)numdimensions - 0.5)) + 1;
    while ((invexp3 == 0) || (invexp3 == invexp1) || (invexp3 == invexp2)) {
      r=drand48();
      invexp3=(int)(r * ((double)numdimensions + 0.5)) + 1;
    }
#endif

    // sort exponents to ensure real roots
    if ((invexp1 < invexp2) && (invexp1 < invexp3)) {
      leftinvexp = invexp1;
    }
    if ((invexp2 < invexp1) && (invexp2 < invexp3)) {
      leftinvexp = invexp2;
    }
    if ((invexp3 < invexp1) && (invexp3 < invexp2)) {
      leftinvexp = invexp3;
    }
    if (((invexp1 < invexp2) && (invexp1 > invexp3)) || ((invexp1 > invexp2) && (invexp1 < invexp3))) {
      middleinvexp = invexp1;
    }
    if (((invexp2 < invexp1) && (invexp2 > invexp3)) || ((invexp2 > invexp1) && (invexp2 < invexp3))) {
      middleinvexp = invexp2;
    }
    if (((invexp3 < invexp1) && (invexp3 > invexp2)) || ((invexp3 > invexp1) && (invexp3 < invexp2))) {
      middleinvexp = invexp3;
    }
    if ((invexp1 > invexp2) && (invexp1 > invexp3)) {
      rightinvexp = invexp1;
    }
    if ((invexp2 > invexp1) && (invexp2 > invexp3)) {
      rightinvexp = invexp2;
    }
    if ((invexp3 > invexp1) && (invexp3 > invexp2)) {
      rightinvexp = invexp3;
    }

/*
    // force polyform for debugging or focused searches
    leftinvexp=+1;
    middleinvexp=+2;
    rightinvexp=+3;
*/
/*
    // force polyform for debugging or focused searches
    leftinvexp=-18;
    middleinvexp=+17;
    rightinvexp=+18;
*/

    // init mass, poly terms and strings
    poly[19]=0;
    sprintf(poly, "P%+d%+d%+d", leftinvexp, middleinvexp, rightinvexp);
    matchptr=matchstart;
    nummatches=0;
    for (i=0; i<=2; i++) {
      coffhit[i]=0;
    }

    // phase 1
    solvePolyforCoefficients(multstart, &nummult, &matchptr, &nummatches, coffhit, random_input_count, random_inputs, poly, leftinvexp, middleinvexp, rightinvexp, range);
    if (nummatches > 0) {
      numcoff=0;
      for (i=0; i<=2; i++) {
        if (coffhit[i] != 0) {
          numcoff++;
        }
      }
      if (numcoff == 3) {
        // phase 2
        verifyMatches(matchstart, &nummatches, poly, leftinvexp, middleinvexp, rightinvexp, random_input_count, minsymmetry, maxcomplexity);
      } else {
#ifdef SHOWSTATUS
        printf("status, No complete three-term phase 2 formulas to solve, coffhit: %d, %d, %d\n", coffhit[0], coffhit[1], coffhit[2]);
        fflush(stdout);
#endif
      }
    } else {
#ifdef SHOWSTATUS
      printf("status, No interesting coefficient multipliers found\n");
      fflush(stdout);
#endif
    } // end nummatches
  } // end while 1
  return(0);
}   
