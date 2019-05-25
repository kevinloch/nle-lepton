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

const float version=3.37;

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define UPLOAD
#define MIXEDMASSMODE
#define SHOWSTATUS
#define SHOWSEARCH
#define NEGATIVEEXP

// pick one G reference
#define CODATA_G
//#define ROSI_G 
//#define WIDE_G // wide=larger uncertainty

//#define DEBUG10
//#define DEBUG11
//#define DEBUG12

//#define DEBUG20
//#define DEBUG21
//#define DEBUG22

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
  double mult[13];
  input_use uses;
} multipliers;

typedef struct {
  double tau_sample;
  double G_sample;
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
  char massstr[20];
  char massstrinv[20];
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
  printf("\nuses\n");
  printf("------------------------------\n");
  printf("alpha_em: %d\n", uses->alpha_em);
  printf("v:        %d\n", uses->v);
  printf("G:        %d\n", uses->G);
  printf("mz:       %d\n", uses->mz);
  printf("mw:       %d\n", uses->mw);
  printf("mh0:      %d\n", uses->mh0);
  printf("sin2w:    %d\n", uses->sin2w);
  printf("------------------------------\n");
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


double solvePolyforMasses(int leftinvexp, int middleinvexp, int rightinvexp, matches *leftmatchptr, matches *middlematchptr, matches *rightmatchptr, input_use *alluses, int maxcomplexity) {
  // solve polynomial like function for particle masses using the supplied coefficients, exponents,
  // the electron mass and fine-structure constant as inputs.  Return the relative difference between
  // the computed and experimental muon mass (as an indication of accuracty)
  long int samples=0;
  int i,j;
  double r;
  struct timespec t;
  long seed;
  long seedsec;
  long seedus;
  double mp;
  double left, middle, right;
  double leftstatic, middlestatic, rightstatic;
  double lefts2w, middles2w, rights2w;
  double leftc2w, middlec2w, rightc2w;
  double leftmassterm, middlemassterm, rightmassterm;
  double leftmeterm, middlemeterm, rightmeterm;
  double leftmuterm, middlemuterm, rightmuterm;
  double leftmtterm, middlemtterm, rightmtterm;
  int arange, merange, murange, vrange, grange, hrange, mzrange, mwrange, mh0range, taurange, sin2wrange;
  struct timespec starttime;
  struct timespec endtime;
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
  char execstr[352];
  char outstr01[320];
  char outstr02[320];
  char outstr03[320];
  char outstr04[320];
  char outstr05[320];
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
  long int samplelimit=100000000;
  int invalid=0;
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
  double e_test=0;
  double e_test_last=0;
  double u_test=0;
  double u_test_last=0;
  double t_test=0;
  double t_test_last=0;
  double precision=0;
  double precision_last=0;
  double avg_rel_range_new;
  double a9;
  double a9mp;

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

//  printf("leftexp: %.3e, middleexp: %.3e, rightexp: %.3e\n", leftexp, middleexp, rightexp);

  // float the three used parameters with the highest experimental uncertainty, use the rest as inputs
  unknowns=0;
  // mH0 relative uncertainty 1.3E-3
  if (alluses->mh0 == 1) {
    floatmh0=1;
    unknowns++;
  } else {
    floatmh0=0;
    mh0_center=mh0_ref;
    mh0_range=mh0_ref_error;
  }

  // sin2w relative uncertainty 1.2E-3, sinw: 5.9E-4, cosw: 1.7E-4
  if (alluses->sin2w == 1) {
    if ((alluses->mw == 1) || (alluses->mz == 1)) {  // check if we are explicitly using mw or mz
      mwmzmode=1; // sin2w will be derived from mw and mz
      floatsin2w=0;
      sin2w_center=sin2w_ref;
      sin2w_range=sin2w_ref_error;
    } else {
      mwmzmode=0;
      floatsin2w=1;
      unknowns++;
    }
  } else {
    mwmzmode=0;
    floatsin2w=0;
    sin2w_center=sin2w_ref;
    sin2w_range=sin2w_ref_error;
  }

#ifdef WIDE_G // higher uncertainty puts G here when this is selected
  // G(wide) relative uncertainty 6.0E-4, mp(wide): 3.0E-4
  if (alluses->G == 1) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=G_ref;
    G_range=G_ref_error;
  }
#endif

  // after this point we must also check for three unknowns

  // mw relative uncertainty 1.5E-4
  if (((alluses->mw == 1) || (mwmzmode == 1)) && (unknowns < 3)) {
    floatmw=1;
    unknowns++;
  } else {
    floatmw=0;
    mw_center=mw_ref;
    mw_range=mw_ref_error;
  }

#ifdef ROSI_G // higher uncertainty puts G here when this is selected
  // G(rosi) relative uncertainty 1.5E-4, mp(rosi): 7.4E-5
  if ((alluses->G == 1) && (unknowns < 3)) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=G_ref;
    G_range=G_ref_error;
  }
#endif

  // tau relative uncertainty 6.8E-5
  if (unknowns < 3) { // tau is always used but only floated if we have less than 3 uknowns so far
    floattau=1;
    unknowns++;
  } else {
    floattau=0;
    tau_center=mu_ref;
    tau_range=mu_ref_error;
  }

#ifdef CODATA_G
  // G(codata) relative uncertainty 4.6E-5, mp(codata): 2.3E-5
  if ((alluses->G == 1) && (unknowns < 3)) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=G_ref;
    G_range=G_ref_error;
  }
#endif

  // mz relative uncertainty 2.3E-5
  if (((alluses->mz == 1) || (mwmzmode == 1)) && (unknowns < 3)) {
    floatmz=1;
    unknowns++;
  } else {
    floatmz=0;
    mz_center=mz_ref;
    mz_range=mz_ref_error;
  }

  // v relative uncertainty 2.6E-7
  if ((alluses->v == 1) && (unknowns < 3)) {
    floatv=1;
    unknowns++;
  } else {
    floatv=0;
    v_center=v_ref;
    v_range=v_ref_error;
  }

  // mu relative uncertainty 2.3E-8
  if (unknowns < 3) {  // mu is always used but only floated if we have less than 3 uknowns so far
    floatmu=1;
    unknowns++;
  } else {
    floatmu=0;
    mu_center=mu_ref;
    mu_range=mu_ref_error;
  }

  // we always have three unknowns at this point so the rest are never floated

  // me relative uncertainty 3.0E-10
  me_center=me_ref;
  me_range=me_ref_error;

  // alpha relative uncertainty 1.5E-10
  alpha_center=alpha_ref;
  alpha_range=alpha_ref_error;

  // systematically try all non-floated input extremes
#ifdef DEBUG20
  printf("debug, Begin phase 2 input loops, unknowns: %d, floatmu: %d, floatv: %d, floatmz: %d, floatg: %d, floattau: %d, floatmw: %d, floatsin2w: %d, mwmzmode: %d, floatmh0: %d\n", unknowns, floatmu, floatv, floatmz, floatg, floattau, floatmw, floatsin2w, mwmzmode, floatmh0);

  printUses(alluses);
#endif
  for (arange=!alluses->alpha_em; ((arange<=1) && (! invalid)); arange++) {
    // alpha is never floated, it is always used as high precision input when needed
    if (arange == 0) {
      alpha=(alpha_center - alpha_range);
    } else {
      alpha=(alpha_center + alpha_range);
    }  
    a9=pow(alpha, 9.00);
    for (merange=0; ((merange<=1) && (! invalid)); merange++) {
      // me is never floated it is always used as high precision input
      if (merange == 0) {
        me=(me_center - me_range);
      } else {
        me=(me_center + me_range);
      }
        for (murange=floatmu; ((murange<=1) && (! invalid)); murange++) {
          // mu is always used but only floated if necessary
          if (floatmu == 0) {
            if (murange == 0) {
              mu=(mu_center - mu_range);
            } else {
              mu=(mu_center + mu_range);
            }
          }
          for (vrange=!(alluses->v || floatv); ((vrange<=1) && (! invalid)); vrange++) {
            if (floatv == 0) {
              if (vrange == 0) {
                v=(v_center - v_range);
              } else {
                v=(v_center + v_range);
              }
            }
            for (mzrange=!((alluses->mz || mwmzmode) || floatmz); ((mzrange<=1) && (! invalid)); mzrange++) {
              if (floatmz == 0) {
                if (mzrange == 0) {
                  mz=(mz_center - mz_range);
                } else {
                  mz=(mz_center + mz_range);
                }
              }
              for (grange=!(alluses->G || floatg); ((grange<=1) && (! invalid)); grange++) {
                if (floatg == 0) { 
                  if (grange == 0) {
                    G=(G_center - G_range);
                  } else {
                    G=(G_center + G_range);
                  } 
                } 
                for (taurange=floattau; ((taurange<=1) && (! invalid)); taurange++) {
                  // tau is always used but only floated if necessary
                  if (floattau == 0) {
                    if (taurange == 0) {
                      tau=(tau_center - tau_range);
                    } else {
                      tau=(tau_center + tau_range);
                    }
                  }
                  for (mwrange=!((alluses->mw || mwmzmode) || floatmw); ((mwrange<=1) && (! invalid)); mwrange++) {
                    if (floatmw == 0) {
                      if (mwrange == 0) {
                        mw=(mw_center - mw_range);
                      } else {
                        mw=(mw_center + mw_range);
                      }
                    }
                    for (sin2wrange=!(alluses->sin2w || floatsin2w); ((sin2wrange<=1) && (! invalid)); sin2wrange++) {
                      if (floatsin2w == 0) {
                        if (sin2wrange == 0) {
                          sin2w=(sin2w_center - sin2w_range);
                        } else {
                          sin2w=(sin2w_center + sin2w_range);
                        }
                      }
                      for (mh0range=!(alluses->mh0 || floatmh0); ((mh0range<=1) && (! invalid)); mh0range++) {
                        if (floatmh0 == 0) {
                          if (mh0range == 0) {
                            mh0=(mh0_center - mh0_range);
                          } else {
                            mh0=(mh0_center + mh0_range);
                          }
                        }
#ifdef DEBUG20
                printf("debug, Begin phase 2 samples loop, arange: %d, merange: %d, hrange: %d, murange: %d, vrange: %d, mzrange: %d, grange: %d, taurange: %d, mwrange: %d, sin2wrange: %d, mh0range: %d\n", arange, merange, hrange, murange, vrange, mzrange, grange, taurange, mwrange, sin2wrange, mh0range);
#endif
                precision_last=1.0E99;
                //  reset mc test vars and outputs
                if (floatmu == 1) {
                  mu_last=0;
                  mu_center=mu_ref;
                  mu_range=mu_ref * 0.001;
                  mu_range_new=0;
                }
                if (floatv == 1) {
                  v_last=0;
                  v_center=v_ref; 
                  v_range=v_ref * 0.001;
                  v_range_new=0;
                }
                if (floatmz == 1) {
                  mz_last=0;
                  mz_center=mz_ref;
                  mz_range=mz_ref * 0.01;
                  mz_range_new=0;
                }
                if (floatg == 1) {
                  G_last=0;
                  G_center=G_ref; 
                  G_range=G_ref * 0.01;
                  G_range_new=0;
                }
                if (floattau == 1) {
                  tau_last=0;
                  tau_center=tau_ref; 
                  tau_range=tau_ref * 0.01;
                  tau_range_new=0;
                }
                if (floatmw == 1) {
                  mw_last=0;
                  mw_center=mw_ref; 
                  mw_range=mw_ref * 0.1;
                  mw_range_new=0;
                }
                if (floatsin2w == 1) {
                  sin2w_last=0;
                  sin2w_center=sin2w_ref;
                  sin2w_range=sin2w_ref * 0.1;
                  sin2w_range_new=0;
                }
                if (floatmh0 == 1) {
                  mh0_last=0;
                  mh0_center=mh0_ref;
                  mh0_range=mh0_ref * 0.1;
                  mh0_range_new=0;
                }

                precision=0;
                precision_last=1E99;
                e_test=0;
                e_test_last=1E99;
                u_test=0;
                u_test_last=1E99;
                t_test=0;
                t_test_last=1E99;
                for (samples = 0; ((samples < samplelimit) && (precision_last > 1.0E-11) && !((samples > 50000000) && (precision_last > 1.0E-6))); samples++) {
                  // guess random values for mc outputs
                  if (floatmu == 1) {
                    r=drand48();
                    mu=((mu_center - mu_range) + (r * 2 * mu_range));
                    i=0;
                    while ((mu < mu_ref * 0.999) || (mu > (mu_ref * 1.001))) { // sanity check on mu mass to help solve correct root and speed convergance
                      if (i > 50) { // safety valve in case search gets out of bounds
                        i=0;
                        mu_last=0;
                        mu_center=mu_ref;
                        mu_range=mu_ref * 0.001;
                        mu_range_new=0;
                      }
                      r=drand48();
                      mu=((mu_center - mu_range) + (r * 2 * mu_range));
                      i++;
                    }
                  }
                  if (floatv == 1) {
                    r=drand48();
                    v=((v_center - v_range) + (r * 2 * v_range));
                    i=0;
                    while ((v < (v_ref * 0.999)) || (v > (v_ref * 1.001))) { // sanity check to help convergance
                      if (i > 50) { // safety valve in case search gets out of bounds
                        i=0;
                        v_last=0;
                        v_center=v_ref;
                        v_range=v_ref * 0.001;
                        v_range_new=0;
                      }
                      r=drand48();
                      v=((v_center - v_range) + (r * 2 * v_range));
                      i++;
                    }
                  }
                  if (floatmz == 1) {
                    r=drand48();
                    mz=((mz_center - mz_range) + (r * 2 * mz_range));
                    i=0;
                    while ((mz < mz_ref * 0.99) || (mz > (mz_ref * 1.01))) { // sanity check on mz mass to help solve correct root and speed convergance
                      if (i > 50) { // safety valve in case search gets out of bounds
                        i=0;
                        mz_last=0;
                        mz_center=mz_ref;
                        mz_range=mz_ref * 0.01;
                        mz_range_new=0;
                      }
                      r=drand48();
                      mz=((mz_center - mz_range) + (r * 2 * mz_range));
                      i++;
                    }
                  }
                  if (floatg == 1) {
                    r=drand48();
                    G=((G_center - G_range) + (r * 2 * G_range));
                    i=0;
                    while ((G < (G_ref * 0.99)) || (G > (G_ref * 1.01))) {  // sanity check to help convergance
                      if (i > 50) {  // safety valve in case search gets out of bounds
                        i=0;
                        G_last=0;
                        G_center=G_ref;
                        G_range=G_ref * 0.01;
                        G_range_new=0;
                      }
                      r=drand48();
                      G=((G_center - G_range) + (r * 2 * G_range));
                      i++;
                    }
                    a9mp=a9 * kg_to_ev * sqrt(hbar_ref * c_ref / G);
                  }
                  if (floattau == 1) {
                    r=drand48();
                    tau=((tau_center - tau_range) + (r * 2 * tau_range));
                    i=0;
                    while ((tau < tau_ref * 0.99) || (tau > (tau_ref * 1.01))) { // sanity check on tau mass to help solve correct root and speed convergance
                      if (i > 50) { // safety valve in case search gets out of bounds
                        i=0;
                        tau_last=0;
                        tau_center=tau_ref;
                        tau_range=tau_ref * 0.01;
                        tau_range_new=0;
                      }
                      r=drand48();
                      tau=((tau_center - tau_range) + (r * 2 * tau_range));
                      i++;
                    }
                  }
                  if (floatmw == 1) {
                    r=drand48();
                    mw=((mw_center - mw_range) + (r * 2 * mw_range));
                    i=0;
                    while ((mw < mw_ref * 0.9) || (mw > (mw_ref * 1.1))) { // sanity check on mw mass to help solve correct root and speed convergance
                      if (i > 50) { // safety valve in case search gets out of bounds
                        i=0;
                        mw_last=0;
                        mw_center=mw_ref;
                        mw_range=mw_ref * 0.1;
                        mw_range_new=0;
                      }
                      r=drand48();
                      mw=((mw_center - mw_range) + (r * 2 * mw_range));
                      i++;
                    }
                  }
                  if (alluses->sin2w == 1) {  
                    if (mwmzmode == 0) {  // float sin2w
                      r=drand48();
                      sin2w=((sin2w_center - sin2w_range) + (r * 2 * sin2w_range));
                      i=0;
                      while ((sin2w < sin2w_ref * 0.9) || (sin2w > (sin2w_ref * 1.1))) { // sanity check on sin2w to help solve correct root and speed convergance
                        if (i > 50) { // safety valve in case search gets out of bounds
                          i=0;
                          sin2w_last=0;
                          sin2w_center=sin2w_ref; 
                          sin2w_range=sin2w_ref * 0.1;
                          sin2w_range_new=0;
                        }
                        r=drand48();
                        sin2w=((sin2w_center - sin2w_range) + (r * 2 * sin2w_range));
                        i++;
                      }
                      cos2w=1.0 - sin2w;
                    } else { // derive sin2w from mw/mz
                      cos2w=pow((mw/mz), 2.0);
                      sin2w=1.0 - cos2w;
                    } // end mwmzmode
                  } // end alluses sin2w
                  if (floatmh0 == 1) {
                    r=drand48();
                    mh0=((mh0_center - mh0_range) + (r * 2 * mh0_range));
                    i=0;
                    while ((mh0 < mh0_ref * 0.9) || (mh0 > (mh0_ref * 1.1))) { // sanity check on mh0 mass to help solve correct root and speed convergance
                      if (i > 50) { // safety valve in case search gets out of bounds
                        i=0;
                        mh0_last=0;
                        mh0_center=mh0_ref;
                        mh0_range=mh0_ref * 0.1;
                        mh0_range_new=0;
                      }
                      r=drand48();
                      mh0=((mh0_center - mh0_range) + (r * 2 * mh0_range));
                      i++;
                    }
                  }

                  if (leftmatchptr->massratio == 0) {
                    leftmassterm=a9mp;
                  } else if (leftmatchptr->massratio == 1) {
                    leftmassterm=v / sqrt(2);
                  } else if (leftmatchptr->massratio == 2) {
                    leftmassterm=mz;
                  } else if (leftmatchptr->massratio == 3) {
                    leftmassterm=mw;
                  } else if (leftmatchptr->massratio == 4) {
                    leftmassterm=mh0;
                  }
                  if (middlematchptr->massratio == 0) {
                    middlemassterm=a9mp;
                  } else if (middlematchptr->massratio == 1) {
                    middlemassterm=v / sqrt(2);
                  } else if (middlematchptr->massratio == 2) {
                    middlemassterm=mz;
                  } else if (middlematchptr->massratio == 3) {
                    middlemassterm=mw;
                  } else if (middlematchptr->massratio == 4) {
                    middlemassterm=mh0;
                  }
                  if (rightmatchptr->massratio == 0) {
                    rightmassterm=a9mp;
                  } else if (rightmatchptr->massratio == 1) {
                    rightmassterm=v / sqrt(2);
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

#ifdef DEBUG22
                  printf("sample: %ld, left: %.6e, leftstatic: %.6e, lefts2w: %.6e, leftc2w: %.6e, lefts2wupout: %d, lefts2wdownout: %d, leftc2wupout: %d, leftc2wdownout: %d\n", samples, left, leftstatic, lefts2w, leftc2w, leftmatchptr->s2wupout, leftmatchptr->s2wdownout, leftmatchptr->c2wupout, leftmatchptr->c2wdownout);
                  printf("sample: %ld, middle: %.6e, middlestatic: %.6e, middles2w: %.6e, middlec2w: %.6e, middles2wupout: %d, middles2wdownout: %d, middlec2wupout: %d, middlec2wdownout: %d\n", samples, middle, middlestatic, middles2w, middlec2w, middlematchptr->s2wupout, middlematchptr->s2wdownout, middlematchptr->c2wupout, middlematchptr->c2wdownout);
                  printf("sample: %ld, right: %.6e, rightstatic: %.6e, rights2w: %.6e, rightc2w: %.6e, rights2wupout: %d, rights2wdownout: %d, rightc2wupout: %d, rightc2wdownout: %d\n", samples, right, rightstatic, rights2w, rightc2w, rightmatchptr->s2wupout, rightmatchptr->s2wdownout, rightmatchptr->c2wupout, rightmatchptr->c2wdownout);
                  printf("sample: %ld, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, left: %.3e, middle: %.3e, right: %.3e\n", samples, e_test, u_test, t_test, left, middle, right);
                  fflush(stdout);
#endif
                  if (e_test > 0.0) {
                    if (u_test > 0.0) {
                      if (((e_test / u_test) < 25.0) && ((u_test / e_test) < 25.0)) {
                        if (t_test > 0.0) {
                          if (((e_test / t_test) < 25.0) && ((t_test / e_test) < 25.0) &&\
                              ((u_test / t_test) < 25.0) && ((t_test / u_test) < 25.0)) {
                            precision=e_test + u_test + t_test;
#ifdef DEBUG22
                            printf ("%ld: precision:  %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, left: %.3e, middle: %.3e, right: %.3e\n", samples, precision, e_test, u_test, t_test, left, middle, right);
                            fflush(stdout);
#endif
                            if ((precision > 0) && (precision < (precision_last * 0.99))) { // 99999
#ifdef DEBUG21
                              printf ("%ld: precision:  %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e\n", samples, precision, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range);
                              fflush(stdout);
#endif
                              precision_last=precision;
                              e_test_last=e_test;
                              u_test_last=u_test;
                              t_test_last=t_test;
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
                              if (floatmu == 1) {
                                mu_range=(mu * precision);
                                mu_center=mu;
                              }
                              if (floatv == 1) {
                                v_range=(v * precision);
                                v_center=v;
                              }
                              if (floatg == 1) {
                                G_range=(G * precision);
                                G_center=G;
                              }
                              if (floatmz == 1) {
                                mz_range=(mz * precision);
                                mz_center=mz;
                              }
                              if (floatmw == 1) {
                                mw_range=(mw * precision);
                                mw_center=mw;
                              }
                              if (floatmh0 == 1) {
                                mh0_range=(mh0 * precision);
                                mh0_center=mh0;
                              }
                              tau_range=(tau * precision);
                              tau_center=tau;
                              sin2w_range=(sin2w * precision);
                              sin2w_center=sin2w;

                            } // end if  precision < precision_last
                          } // end if t_test/u_test/e_test
                        } // end if t_test > 0
                      } // end if e_test / u_test
                    } // end if u_test > 0
                  } // end if e_test > 0
                } // end for samples 

#ifdef DEBUG20
                printf("debug, Finished phase 2 samples loop, samples: %ld, mass mode: %d%d%d, precision: %.6e\n", samples, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, precision_last);
#endif
                if ((samples >= (samplelimit-1)) || (precision_last >= 1.0E-11)) {
                  invalid=1;
#ifdef SHOWSTATUS
                  printf ("status, Failed to solve phase 2 formula for masses within sample limit, samples: %ld, e_test:  %.3e, u_test:  %.3e, t_test: %.3e tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e, precision: %.3e\n", samples, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range, precision_last);
                  fflush(stdout);
#endif
                } else {

#ifdef DEBUG20
                  //printf("debug, phase 2 samples loop valid, updating output vars\n");
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

                } // end if precision < 1E-11
                      } // end mh0range
                    } // end sin2wrange
                  } // end mwrange
                } // end taurange
              } // end grange
            } // end mzrange
          } // end vrange
        }  // end murange
    } // end merange
  }  // end arange

  if (invalid == 0) {
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsedtime=(((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9)) - 1.0;

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
      if (invexp == 1) {
        sprintf(leftformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  matchptr->massstr);
      } else if (invexp == -1) {
        sprintf(leftformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstrinv);
      } else if (invexp > 1) {
        sprintf(leftformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstr, invexp);
      } else {
        sprintf(leftformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstrinv, -invexp);
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
      if (invexp == 1) {
        sprintf(middleformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  matchptr->massstr);
      } else if (invexp == -1) {
        sprintf(middleformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstrinv);
      } else if (invexp > 1) {
        sprintf(middleformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstr, invexp);
      } else {
        sprintf(middleformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstrinv, -invexp);
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
      if (invexp == 1) {
        sprintf(rightformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  matchptr->massstr);
      } else if (invexp == -1) {
        sprintf(rightformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstrinv);
      } else if (invexp > 1) {
        sprintf(rightformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstr, invexp);
      } else {
        sprintf(rightformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, matchptr->massstrinv, -invexp);
      }
/* end create formula output */

      sprintf(outstr01, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
      printf("%s\n", outstr01);
      sprintf(outstr02, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
      printf("%s\n", outstr02);
      sprintf(outstr03, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
      printf("%s\n", outstr03);
      if (alluses->alpha_em == 1) {
        sprintf(usedasinput, "*");
        sprintf(outstr04, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 04, | alpha_em   || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |                |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, alpha_out_c, alpha_out_error, alpha_out_relerror, alpha_ref, alpha_ref_error, alpha_ref_relerror, alpha_out_diff, alpha_out_reldiff, usedasinput);
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
        sprintf(outstr06, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 06, | v          || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, v_out_c, v_out_error, v_out_relerror, v_ref, v_ref_error, v_ref_relerror, v_out_diff, v_out_reldiff, usedasinput, usedasoutput);
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
        sprintf(outstr07, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 08, | mZ         || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mz_out_c, mz_out_error, mz_out_relerror, mz_ref, mz_ref_error, mz_ref_relerror, mz_out_diff, mz_out_reldiff, usedasinput, usedasoutput);
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
        sprintf(outstr08, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 07, | G          || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, G_out_c, G_out_error, G_out_relerror, G_ref, G_ref_error, G_ref_relerror, G_out_diff, G_out_reldiff, usedasinput, usedasoutput);
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
        sprintf(outstr09, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 09, | mW         || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mw_out_c, mw_out_error, mw_out_relerror, mw_ref, mw_ref_error, mw_ref_relerror, mw_out_diff, mw_out_reldiff, usedasinput, usedasoutput);
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
        sprintf(outstr10, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 11, | sin2w      || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |               |       %s        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, sin2w_out_c, sin2w_out_error, sin2w_out_relerror, sin2w_ref, sin2w_ref_error, sin2w_ref_relerror, sin2w_out_diff, sin2w_out_reldiff, usedasoutput);
        printf("%s\n", outstr10);
      } else {
        outstr10[0]=0;
      }
      if (alluses->mh0 == 1) {
        sprintf(usedasoutput, "*");
        sprintf(outstr11, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 10, | mH0        || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |               |       %s        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mh0_out_c, mh0_out_error, mh0_out_relerror, mh0_ref, mh0_ref_error, mh0_ref_relerror, mh0_out_diff, mh0_out_reldiff, usedasoutput);
        printf("%s\n", outstr11);
      } else {
        outstr11[0]=0;
      }
      sprintf(outstr12, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 12, | Electron   || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       *       |                |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, me_out_c, me_out_error, me_out_relerror, me_ref, me_ref_error, me_ref_relerror, me_out_diff, me_out_reldiff);
      printf("%s\n", outstr12);
      if (floatmu == 1) {
        sprintf(usedasinput, " ");
        sprintf(usedasoutput, "*");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr13, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 13, | Muon       || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mu_out_c, mu_out_error, mu_out_relerror, mu_ref, mu_ref_error, mu_ref_relerror, mu_out_diff, mu_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr13);
      sprintf(outstr14, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 14, | Tau        || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |               |       *        |", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, tau_out_c, tau_out_error, tau_out_relerror, tau_ref, tau_ref_error, tau_ref_relerror, tau_out_diff, tau_out_reldiff);
      printf("%s\n", outstr14);
      sprintf(outstr15, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
      printf("%s\n", outstr15);
      sprintf(outstr16, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 16, %5d, %s", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, leftmatchptr->matchcomplexity, leftformulastr);
      printf("%s\n", outstr16);
      sprintf(outstr17, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 17, %5d, %s", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, middlematchptr->matchcomplexity, middleformulastr);
      printf("%s\n", outstr17);
      sprintf(outstr18, "result, %.4f, %3d, %3d, P%+d%+d%+d, M%d%d%d, %12lld, 18, %5d, %s", combinedscore, symmetry, complexity, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, rightmatchptr->matchcomplexity, rightformulastr);
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
        if (outstr05[0] != 0) {
          sprintf(execstr, "curl -s \"http://localhost/lepton/%s\" > /dev/null 2>&1\n", underscore(outstr05, 320));
          system(execstr);
        }
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
  } // end if not invalid
  if (invalid == 1) {
    return(1.0E+99);
  } else {
    return(precision_last);
  }
}

int verifyMatches(matches *matchstart, int *nummatches, int leftinvexp, int middleinvexp, int rightinvexp, int random_input_count, int minsymmetry, int maxcomplexity) {
  //  For polynomials with interesting coefficients on all three exponent terms,
  //  separate the match table into a separate list for each exponent, then test all unique combinations of coefficients
  //  for accuracy by comparing the computed muon mass to it's experimental value.  Print results matching a minimum threshold of interest.
  int i,j;
  matches *match;
  double muon_reldiff;
  matches *leftmatches;
  matches *middlematches;
  matches *rightmatches;
  matches *leftmatchptr;
  matches *middlematchptr;
  matches *rightmatchptr;
  matches *tmpmatchptr;
  int numleftmatches, nummiddlematches, numrightmatches;
  int l,m,r,h;
  char hash[52];
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
            strncpy(tmpmatchptr->massstr, match->massstr, 20);
            strncpy(tmpmatchptr->massstrinv, match->massstrinv, 20);
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
        strncpy(leftmatchptr->massstr, match->massstr, 20);
        strncpy(leftmatchptr->massstrinv, match->massstrinv, 20);
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
            strncpy(tmpmatchptr->massstr, match->massstr, 20);
            strncpy(tmpmatchptr->massstrinv, match->massstrinv, 20);
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
        strncpy(middlematchptr->massstr, match->massstr, 20);
        strncpy(middlematchptr->massstrinv, match->massstrinv, 20);
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
            strncpy(tmpmatchptr->massstr, match->massstr, 20);
            strncpy(tmpmatchptr->massstrinv, match->massstrinv, 20);
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
        strncpy(rightmatchptr->massstr, match->massstr, 20);
        strncpy(rightmatchptr->massstrinv, match->massstrinv, 20);
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
  printf("status, Solving phase 2 formulas for masses, random input: %d, polyform P%+d%+d%+d,                 progress: total (0/%ld) left (0/%d) middle (0/%d) right (0/%d)\n", random_input_count, leftinvexp, middleinvexp, rightinvexp, totalcombos, numleftmatches, nummiddlematches, numrightmatches);
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
          //printUses(&leftuses);
          //printUses(&middleuses);
          //printUses(&rightuses);
          //printUses(&alluses);
          clock_gettime(CLOCK_REALTIME, &starttime);
          precision=solvePolyforMasses(leftinvexp, middleinvexp, rightinvexp, leftmatchptr, middlematchptr, rightmatchptr, &alluses, maxcomplexity);
#ifdef SHOWSTATUS
          clock_gettime(CLOCK_REALTIME, &endtime);
          elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
//          if ((m == 0) && (r ==0)) {
            if (precision > 1.0E-11) {
              printf("status, Failed to solve phase 2 formula for masses within sample limit, random input: %d, polyform P%+d%+d%+d, mass mode: %d%d%d, progress: total (%ld/%ld) left (%d/%d) middle (%d/%d) right (%d/%d), precision: %.3e, (%6.4fs)\n", random_input_count, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, combo, totalcombos, l+1, numleftmatches, m+1, nummiddlematches, r+1, numrightmatches, precision, elapsedtime);
              fflush(stdout);
            } else {
              printf("status, Solved  phase 2 formula  for masses, random input: %d, polyform P%+d%+d%+d, mass mode: %d%d%d, progress: total (%ld/%ld) left (%d/%d) middle (%d/%d) right (%d/%d), precision: %.3e, (%6.4fs)\n", random_input_count, leftinvexp, middleinvexp, rightinvexp, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, combo, totalcombos, l+1, numleftmatches, m+1, nummiddlematches, r+1, numrightmatches, precision, elapsedtime);
              fflush(stdout);
//            }
          }
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
*/
  if (range == 4) {
    rangehigh=1.0001;
    rangelow=0.9999;
  } else if (range == 5) {
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
  int phiupin=0, phidownin=1;
  int piupin=0, pidownin=1;
  int aupin=0, adownin=1;
  int e2upin=0, e2downin=2;
  int e3upin=0, e3downin=2;
  int e5upin=0, e5downin=2;
  int nbvupin=0;
  int nbsupin=0;
  double updownin;
  double phiin;
  double piin;
  double ain;
  double e2in;
  double e3in;
  double e5in;
  double phi;
  double nbv[13];
  double nbs[13];

  unsigned int u, v;
  int i;
  multipliers *m;

  printf("init, Initializing pre-computed multiplier accelerator array (a subset of multiplier terms that easily fits in on-chip cache)\n");
  clock_gettime(CLOCK_REALTIME, &starttime);

  phi=(1.0 + sqrt(5.0)) / 2.0;

  // nball volume (index is term exponent and n-ball dimension)
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

  // nball surface area (index is term exponend and n-ball dimension - 1)
  // 2pi * nbv[n]
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
/*
            for (phiupin = -1; phiupin <=1; phiupin++) {
              for (phidownin = 1; phidownin <= 2; phidownin++) {
                u=abs(phiupin);
                v=phidownin;
                if (gcd(u, v) == 1) {
                  phiin=updownin * pow(phi, ((float)phiupin / (float)phidownin));
*/

            for (piupin = -2; piupin <=2; piupin++) {
//              for (pidownin=1; pidownin <= 1; pidownin++) {
//                u=abs(piupin);
//                v=pidownin;
//                if (gcd(u, v) == 1) {
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
/*
                              for (e3upin=-1; e3upin<=1; e3upin++) {
                                for (e3downin=2; e3downin <= 2; e3downin++) {
                                  u=abs(e3upin);
                                  v=e3downin;
                                  if (gcd2(u, v) == 1) {
                                    e3in=e2in * pow(3.0, ((float)e3upin / (float)e3downin));
                                    for (e5upin=-1; e5upin<=1; e5upin++) {
                                      for (e5downin=2; e5downin <= 2; e5downin++) {
                                        u=abs(e5upin);
                                        v=e5downin;
                                        if (gcd2(u, v) == 1) {
                                          e5in=e3in * pow(5.0, ((float)e5upin / (float)e5downin));
*/
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
/*
                                                 + (abs(phiupin) + phidownin)\
*/
                                                 + (abs(piupin) + pidownin)\
                                                 + (abs(aupin) + adownin)\
                                                 + (abs(e2upin) + e2downin);
/*
                                                 + (abs(e3upin) + e3downin)\
                                                 + (abs(e5upin) + e5downin);
*/
                                          for (i=1; i<=12; i++) {
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
                                            }  // if i
                                          } // for i
                                          initUses(&mult->uses);
                                          if (aupin != 0) {
                                            mult->uses.alpha_em=1;
                                          }
                                          *nummult=*nummult+1;
                                          mult++;
/*
                                        } // e5up gcd
                                      } // e5downin
                                    } // e5upin
                                  } // e3up gcd
                                } // e3downin
                              } // e3upin
*/
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
/*
}  // gcd phiin
} // phidownin
} // phiupin
*/
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

void cscanner(multipliers *multstart, int *nummult, matches **matchptr, int *nummatches, int *coffhit, int range, char *poly, char *massstr, char *massstrinv, int massratio,  int leftinvexp, int middleinvexp, int rightinvexp, double cleft, double cmiddle, double cright, input_use *search_uses, random_input random_inputs) {
  //  Each coefficient is multiplied by various numbers and the resulting value is tested to see if it is
  // close to an interesting integer or simple rational number.  The results are stored in the match table.
  double rangehigh=0;
  double rangelow=0;
  int i,j;
  double multiplierout;
  double multiplierin;
  char outs[21];
  double outld;
  multipliers *mult;
  unsigned int u, v;
  matches *match;
  int complexity;

  int upout=1, downout=1;
  int phiupout=0, phidownout=1;
  int piupout=0, pidownout=1;
  int aupout=1, adownout=1;
  int e2upout=0, e2downout=2;
  int s2wupout=0, s2wdownout=1;
  int c2wupout=0, c2wdownout=1;
  double updownout;
  double phiout;
  double piout;
  double aout;
  double e2out;
  double s2wout=1.0;
  double c2wout=1.0;
  double s2win=1.0;
  double c2win=1.0;

  double phi;
  double sin2w;
  double cos2w;

  phi=(1.0 + sqrt(5.0)) / 2.0;

  sin2w=random_inputs.sin2w_sample;
  cos2w=1.0 - sin2w;

  match=*matchptr;
  for (upout=1; upout <= 3; upout++) {
#ifdef SHOWSTATUS
//    printf("status, Scanning for coefficient multipliers that match interesting integer or simple rational numbers, progress: %2d%%\n", ((upout - 1) * 33));
//    fflush(stdout);
#endif
    for (downout=1; downout <=3; downout++) {
      u=upout;
      v=downout;
      if (gcd(u, v) == 1) {
        updownout=(double)upout / (double)downout;
/*
        for (phiupout=-1; phiupout <= 1; phiupout++) {
          for (phidownout=1; phidownout <= 2; phidownout++) {
            u=abs(phiupout);
            v=phidownout;
            if (gcd(u, v) == 1) {
              phiout=updownout * pow(phi, ((float)phiupout / (float)phidownout));
*/

        for (piupout=-2; piupout <= 2; piupout++) {
//          for (pidownout=1; pidownout <= 1; pidownout++) {
//            u=abs(piupout);
//            v=pidownout;
//            if (gcd(u, v) == 1) {
              piout=updownout * pow(M_PI, ((float)piupout / (float)pidownout));
              for (aupout=-2; aupout <= 2; aupout++) {
                for (adownout=1; adownout <= 2; adownout++) {
                  u=abs(aupout);
                  v=adownout;
                  if (gcd(u, v) == 1) {
                    aout=piout * pow(alpha_ref, ((float)aupout / (float)adownout));
                    for (e2upout=-1; e2upout <= 1; e2upout++) {
//                      for (e2downout=2; e2downout <= 2; e2downout++) {
//                        u=abs(e2upout);
//                        v=e2downout;
//                        if (gcd2(u, v) == 1) {
                          multiplierout=aout * pow(2.0, ((float)e2upout / (float)e2downout));
// s2w and c2w are calculated separately since they contain significant uncertainty and are reconstructed separately in phase 2
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
                                      mult=multstart;
                                      for (i=0; i<*nummult; i++) { 
                                        // test multiplier against left coefficient
//if (!((s2wupout != 0) && (c2wupout != 0))) {
if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == leftinvexp) || (s2wdownout == (leftinvexp * 2))) {
  if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == leftinvexp) || (c2wdownout == (leftinvexp * 2))) {
//    if ((massratio != 1) || (mult->e2upin == 0)) {
                                        outld=cleft * s2wout * c2wout * multiplierout * mult->mult[abs(leftinvexp)];
//printf("lefttest: %e, multiplierout: %e, mult: %e, cleft: %e\n", outld, multiplierout, mult->mult[abs(leftinvexp)], cleft);
                                        if (interesting(range, outld)) {
//printf("*********************lefttest: %e, multiplierout: %e, mult: %e, cleft: %e\n", outld, multiplierout, mult->mult[abs(leftinvexp)], cleft);
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
                                          match->match=outld;
                                          initUses(&match->uses);
                                          addUses(&match->uses, search_uses);
                                          addUses(&match->uses, &mult->uses);
                                          if (aupout != 0) {
                                            match->uses.alpha_em=1;
                                          }
                                          if ((s2wupout != 0) || (c2wupout != 0)) {
                                            match->uses.sin2w=1;
                                          }
                                          match->matchmult=multiplierout * mult->mult[abs(leftinvexp)];
                                          strncpy(match->massstr, massstr, 20);
                                          strncpy(match->massstrinv, massstrinv, 20);
                                          *nummatches=*nummatches+1;
                                          *matchptr=*matchptr+1;
                                          match=*matchptr;
//    } // sanity check sqrt(2) in
  } // sanity check c2w left
} // sanity check s2w left
//} // sanity check not both s2w and c2w
                                        }  // if interesting left

                                        // test multiplier against middle coefficient
//if (!((s2wupout != 0) && (c2wupout != 0))) {
if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == middleinvexp) || (s2wdownout == (middleinvexp * 2))) {
  if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == middleinvexp) || (c2wdownout == (middleinvexp * 2))) {
//    if ((massratio != 1) || (mult->e2upin == 0)) {
                                        outld=cmiddle * s2wout * c2wout * multiplierout * mult->mult[abs(middleinvexp)];
//printf("middletest: %e, multiplierout: %e, mult: %e, cmiddle: %e\n", outld, multiplierout, mult->mult[abs(middleinvexp)], cmiddle);
                                        if (interesting(range, outld)) { 
//printf("*********************middletest: %e, multiplierout: %e, mult: %e, cmiddle: %e\n", outld, multiplierout, mult->mult[abs(middleinvexp)], cmiddle);
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
                                          match->match=outld;
                                          initUses(&match->uses);
                                          addUses(&match->uses, search_uses);
                                          addUses(&match->uses, &mult->uses);
                                          if (aupout != 0) {
                                            match->uses.alpha_em=1;
                                          }
                                          if ((s2wupout != 0) || (c2wupout != 0)) {
                                            match->uses.sin2w=1;
                                          }
                                          match->matchmult=multiplierout * mult->mult[abs(middleinvexp)];
                                          strncpy(match->massstr, massstr, 20);
                                          strncpy(match->massstrinv, massstrinv, 20);
                                          *nummatches=*nummatches+1;
                                          *matchptr=*matchptr+1;
                                          match=*matchptr;
//    } // sanity check sqrt(2) in
  } // sanity check c2w middle
} // sanity check s2w middle
//} // sanity check not both s2w and c2w
                                        }  // if interesting middle

                                        // test multiplier against right coefficient
//if (!((s2wupout != 0) && (c2wupout != 0))) {
if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == rightinvexp) || (s2wdownout == (rightinvexp * 2))) {
  if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == rightinvexp) || (c2wdownout == (rightinvexp * 2))) {
//    if ((massratio != 1) || (mult->e2upin == 0)) {
                                        outld=cright * s2wout * c2wout * multiplierout * mult->mult[abs(rightinvexp)];
//printf("righttest: %e, multiplierout: %e, mult: %e, cright: %e\n", outld, multiplierout, mult->mult[abs(rightinvexp)], cright);
                                        if (interesting(range, outld)) { 
//printf("*********************righttest: %e, multiplierout: %e, mult: %e, cright: %e\n", outld, multiplierout, mult->mult[abs(rightinvexp)], cright);
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
                                          match->match=outld;
                                          initUses(&match->uses);
                                          addUses(&match->uses, search_uses);
                                          addUses(&match->uses, &mult->uses);
                                          if (aupout != 0) {
                                            match->uses.alpha_em=1;
                                          }
                                          if ((s2wupout != 0) || (c2wupout != 0)) {
                                            match->uses.sin2w=1;
                                          }
                                          match->matchmult=multiplierout * mult->mult[abs(rightinvexp)];
                                          strncpy(match->massstr, massstr, 20);
                                          strncpy(match->massstrinv, massstrinv, 20);
                                          *nummatches=*nummatches+1;
                                          *matchptr=*matchptr+1;
                                          match=*matchptr;
//    } // sanity check sqrt(2) in
  } // sanity check c2w right
} // sanity check s2w right
//} // sanity check not both s2w and c2w
                                        }  // if interesting right

                                        mult++;
                                      }  // for i
                                    } // gcd c2wout
                                  } // c2wdownout
                                } // c2wupout
/*
                              } // gcd s2wout
                            } // s2wdownout
*/
                          } // s2wupout
                        } // gcd e2out
                      } // e2downout
                    } // e2upout
                  } // gcd aout
                } // adownout
              } // aupout
/*
            } // gcd piout
          } // pidownout
*/
       } // piupout
/*
} // gcd phiout
} // phidownout
} // phiupout
*/
      } // gcd updownout
    } // downout
  } // upout
}

int solvePolyforCoefficients(multipliers *multstart, int *nummult, matches **matchptr, int *nummatches, int *coffhit, int random_input_count, random_input random_inputs, double massterm, char *massstr, char *massstrinv, char *poly, int leftinvexp, int middleinvexp, int rightinvexp, int massratio, int range, input_use *search_uses) {
  // solve a three term polynomial-like formula for the unknown coefficients given known roots (particle masses)
  int i;
  long int samples=0;
  long int runs=0;
  char s1[21];
  char s2[21];
  char s3[21];
  char s4[21];
  char s5[21];
  char s6[21];
  char s7[21];
  char s12[21];
  struct timespec starttime;
  struct timespec endtime;
  double elapsedtime;
  int nummatchstart;
  int nummatchend;
  int goodc;
  int pml, pmr, plr;
  int progress;
  int stalled;
  int stalledboost;

  clock_gettime(CLOCK_REALTIME, &starttime);

  //  mc test vars
  long double r;
  long double e_test_1;
  long double e_test_2;
  long double e_test_3;
  long double e_test=0;
  long double e_test_last=9999999999;
  long double u_test_1;
  long double u_test_2;
  long double u_test_3;
  long double u_test=0;
  long double u_test_last=9999999999;
  long double t_test_1;
  long double t_test_2;
  long double t_test_3;
  long double t_test=0;
  long double t_test_last=9999999999;
  long double precision=0;
  long double precision_last=9999999999;
  long double massleft_e=0;
  long double massmiddle_e=0;
  long double massright_e=0;
  long double massleft_u=0;
  long double massmiddle_u=0;
  long double massright_u=0;
  long double massleft_t=0;
  long double massmiddle_t=0;
  long double massright_t=0;
  double mp=0;

  // mc outputs
  long double cleft=0;
  long double cleft_last=0;
  long double cleft_center=0;
  long double cleft_range=0;
  long double cleft_range_new=0;
  long double cmiddle=0;
  long double cmiddle_last=0;
  long double cmiddle_center=0;
  long double cmiddle_range=0;
  long double cmiddle_range_new=0;
  long double cright=0;
  long double cright_last=0;
  long double cright_center=0;
  long double cright_range=0;
  long double cright_range_new=0;

  massleft_e=powl(((long double)me_ref / (long double)massterm), (1.0 / (long double)leftinvexp));
  massleft_u=powl(((long double)mu_ref / (long double)massterm), (1.0 / (long double)leftinvexp));
  massleft_t=powl(((long double)random_inputs.tau_sample / (long double)massterm), (1.0 / (long double)leftinvexp));
  massmiddle_e=powl(((long double)me_ref / (long double)massterm), (1.0 / (long double)middleinvexp));
  massmiddle_u=powl(((long double)mu_ref / (long double)massterm), (1.0 / (long double)middleinvexp));
  massmiddle_t=powl(((long double)random_inputs.tau_sample / (long double)massterm), (1.0 / (long double)middleinvexp));
  massright_e=powl(((long double)me_ref / (long double)massterm), (1.0 / (long double)rightinvexp));
  massright_u=powl(((long double)mu_ref / (long double)massterm), (1.0 / (long double)rightinvexp));
  massright_t=powl(((long double)random_inputs.tau_sample / (long double)massterm), (1.0 / (long double)rightinvexp));

#ifdef SHOWSTATUS
  printf("status, Solving phase 1 formula for coefficients, polyform: %s, mass ratio: %s\n", poly, massstr);
  fflush(stdout);
#endif

  //  solve formula for coefficients
  runs=0;
  precision_last=1.0E30;
  while (precision_last > 1.0E-11) {
    runs++;
#ifdef SHOWSTATUS
//printf("status, run: %ld\n", runs);
//fflush(stdout);
#endif
    // scan all 8 possible coefficient relationships
    for (pml=1; ((pml >= 0) && (precision_last > 1.0E-11)); pml--) {
      for (pmr=1; ((pmr >= 0) && (precision_last > 1.0E-11)); pmr--) {
        for (plr=1; ((plr >= 0) && (precision_last > 1.0E-11)); plr--) {
#ifdef DEBUG10
          printf("pml: %d, pmr: %d, plr: %d\n", pml, pmr, plr);
          fflush(stdout);
#endif
          if (leftinvexp < 0) {
            cleft_center=1.0E1;
            cleft_range=1.0E1;
          } else {
            cleft_center=1.0E2;
            cleft_range=1.0E2;
          }
          cleft_range_new=0;
          if (middleinvexp < 0) {
            cmiddle_center=1.0E1;
            cmiddle_range=1.0E1;
          } else {
            cmiddle_center=1.0E2;
            cmiddle_range=1.0E2;
          }
          cmiddle_range_new=0;
          if (rightinvexp < 0) {
            cright_center=1.0E1;
            cright_range=1.0E1;
          } else {
            cright_center=1.0E2;
            cright_range=1.0E2;
          }
          cright_range_new=0;
          precision_last=1.0E4;
          i=0;
          progress=0;
          stalled=0;
          stalledboost=0;
          //for (samples = 0; ((samples < 100000000) && (precision_last > 1.0E-11) && !((samples > 100000) && (progress == 0)) && !(((samples % 5000) == 0) && (i == 200)) && !((stalled > 2000000) && (precision_last > 1.0E-1))); samples++) {
          for (samples = 0; ((samples < 100000000) && (precision_last > 1.0E-11) && !((samples > 100000) && (progress == 0)) && !(((samples % 5000) == 0) && (i == 200)) && !((stalledboost > 5)) && !((stalled > 2000000) && (precision_last > 1.0E-1))); samples++) {
            // generate random coefficients within limits
            stalled++;
            if ((stalled % 2000000) == 0) {
              stalledboost++;
              cleft_range = cleft_range *     0.1;
              cmiddle_range = cmiddle_range * 0.1;
              cright_range = cright_range *   0.1;
#ifdef DEBUG12
              printf("stalled boost\n");
              printf("stalled: %d, precision_last: %.9Le, pre-filter precision: %.9Le, e_test: %.9Le, u_test: %.9Le, t_test: %.9Le, cleft: %.9Le, cmiddle: %.9Le, cright: %.9Le, cleftrange: %.9Le, cmiddlerange: %.9Le, crightrange: %.9Le\n", stalled, precision_last, precision, e_test, u_test, t_test, cleft, cmiddle, cright, cleft_range, cmiddle_range, cright_range);
              fflush(stdout);
#endif

            }
            i=0;
            goodc=0;
            while ((i < 200) && (goodc == 0)) {
              i++;
              cleft=0.0;
              cmiddle=0.0;
              cright=0.0;
              while (cleft <= 0.0) {
                r=(long double)drand48();
                cleft=((cleft_center - cleft_range) + (r * 2.0 * cleft_range));
              }
              while (cmiddle <= 0.0) {
                r=(long double)drand48();
                cmiddle=((cmiddle_center - cmiddle_range) + (r * 2.0 * cmiddle_range));
              }                 
              while (cright <= 0.0) {
                r=(long double)drand48();
                cright=((cright_center - cright_range) + (r * 2.0 * cright_range));
              }
              // verify coefficient relationships
              goodc=1;
              if (pml == 1) {
                if (cmiddle < cleft) {
                  goodc=0;
                }
              } else {
                if (cmiddle > cleft) {
                  goodc=0;
                }
              }
              if (pmr == 1) {
                if (cmiddle < cright) {
                  goodc=0;
                }
              } else {
                if (cmiddle > cright) {
                  goodc=0;
                }
              }
              if (plr == 1) {
                if (cleft < cright) {
                  goodc=0;
                }
              } else {
                if (cleft > cright) {
                  goodc=0;
                }
              }
            } // end while goodc
            if (goodc == 1) {
              e_test_1=cleft * massleft_e;
              e_test_2=cmiddle * massmiddle_e;
              e_test_3=cright * massright_e;

              u_test_1=cleft * massleft_u;
              u_test_2=cmiddle * massmiddle_u;
              u_test_3=cright * massright_u;

              t_test_1=cleft * massleft_t;
              t_test_2=cmiddle * massmiddle_t;
              t_test_3=cright * massright_t;

              e_test=e_test_1 - e_test_2 + e_test_3 - 1.0;
//printf("etest: %e\n", e_test);
              if (e_test > 0.0) {
                u_test=u_test_1 - u_test_2 + u_test_3 - 1.0;
//printf("utest: %e\n", u_test);
                if (u_test > 0.0) {
                  if ((progress < 15) || (((e_test / u_test) < 25.0) && ((u_test / e_test) < 25.0))) {
                    t_test=t_test_1 - t_test_2 + t_test_3 - 1.0;
//printf("ttest: %e\n", t_test);
                    if (t_test > 0.0) {
//printf("etest: %.3e, utest: %.3e, ttest: %.3e, e_test_1: %.3e, e_test_2: %.3e, e_test_3: %.3e, u_test_1: %.3e, u_test_2: %.3e, u_test_3: %.3e, t_test_1: %.3e, t_test_2: %.3e, t_test_3: %.3e\n", e_test, u_test, t_test, e_test_1, e_test_2, e_test_3, u_test_1, u_test_2, u_test_3, t_test_1, t_test_2, t_test_3);
                      if ((progress < 15) || (((e_test / t_test) < 25.0) && ((t_test / e_test) < 25.0) &&\
                                              ((u_test / t_test) < 25.0) && ((t_test / u_test) < 25.0))) {
//printf("etest: %.3e, utest: %.3e, ttest: %.3e, e_test_1: %.3e, e_test_2: %.3e, e_test_3: %.3e, u_test_1: %.3e, u_test_2: %.3e, u_test_3: %.3e, t_test_1: %.3e, t_test_2: %.3e, t_test_3: %.3e\n", e_test, u_test, t_test, e_test_1, e_test_2, e_test_3, u_test_1, u_test_2, u_test_3, t_test_1, t_test_2, t_test_3);
                        precision=(e_test + u_test + t_test);
                        if ((precision > 0) && (precision < (precision_last * 0.99))) {
                        //if ((precision > 0) && (precision < precision_last)) {
                          progress++;
                          stalled=0;
                          stalledboost=0;
                          precision_last=precision;
                          e_test_last=e_test;
                          u_test_last=u_test;
                          t_test_last=t_test;
                          cleft_last=cleft;
                          cmiddle_last=cmiddle;
                          cright_last=cright;
                          cleft_range_new=(fabsl(cleft - cleft_center) * 2.20);
                          cleft_range=((cleft_range + cleft_range_new + cleft_range_new) / 3.0);
                          cmiddle_range_new=(fabsl(cmiddle - cmiddle_center) * 2.20);
                          cmiddle_range=((cmiddle_range + cmiddle_range_new + cmiddle_range_new) / 3.0);
                          cright_range_new=(fabsl(cright - cright_center) * 2.20);
                          cright_range=((cright_range + cright_range_new + cright_range_new) / 3.0);
                          cleft_center=cleft;
                          cmiddle_center=cmiddle;
                          cright_center=cright;
#ifdef DEBUG11
                          printf("progress: %d, samples: %ld, i: %d, post-filter precision: %.9Le, e_test: %.9Le, u_test: %.9Le, t_test: %.9Le, cleft: %.9Le, cmiddle: %.9Le, cright: %.9Le, cleftrange: %.9Le, cmiddlerange: %.9Le, crightrange: %.9Le\n", progress, samples, i, precision, e_test, u_test, t_test, cleft, cmiddle, cright, cleft_range, cmiddle_range, cright_range);
                          fflush(stdout);
#endif
                        } // end if precision
                      } // end ttest ratios
                    } //  end ttest > 0
                  } // end utest ratios
                } // end utest > 0
              }  // end etest > 0
            }  // end if goodc
          }  // end for samples
        }  // end for plr
      }  // end for pmr
    } // end for pml
  }  // end while precision_last

  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);

#ifdef SHOWSTATUS
  printf("status, Solved  phase 1 formula for coefficients, random input: %i, polyform:  %s, mass ratio: %s, Tau: %.9e, sin2w: %.9e, massterm: %.9e, runs: %ld, samples: %ld, precision: %.3Le (%6.4fs)\n", random_input_count, poly, massstr, random_inputs.tau_sample, random_inputs.sin2w_sample, massterm, runs, samples, precision, elapsedtime);
  fflush(stdout);
#endif
#ifdef SHOWSEARCH
  printf("search, +-----------------+-----------------+---------+--------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
  printf("search, |    tau mass     |    mass term    |  poly   |  mass ratio  | polynomial term |   coefficient   | c * polyterm(e) | c * polyterm(u) | c * polyterm(t) |\n");
  printf("search, +-----------------+-----------------+---------+--------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
  printf("search, | %.9e | %.9e | %s | %s |     +left       | %.9Le | %.9Le | %.9Le | %.9Le |\n", random_inputs.tau_sample, massterm, poly, massstr, cleft, e_test_1, u_test_1, t_test_1);
  printf("search, | %.9e | %.9e | %s | %s |     -middle     | %.9Le | %.9Le | %.9Le | %.9Le |\n", random_inputs.tau_sample, massterm, poly, massstr, cmiddle, e_test_2, u_test_2, t_test_2);
  printf("search, | %.9e | %.9e | %s | %s |     +right      | %.9Le | %.9Le | %.9Le | %.9Le |\n", random_inputs.tau_sample, massterm, poly, massstr, cright, e_test_3, u_test_3, t_test_3);
  printf("search, +-----------------+-----------------+---------+--------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
  fflush(stdout);
#endif

  clock_gettime(CLOCK_REALTIME, &starttime);
  nummatchstart=*nummatches;

#ifdef SHOWSTATUS
  printf("status, Scanning for coefficient multipliers that match interesting integer or simple rational numbers\n");
#endif
  cscanner(multstart, nummult, matchptr, nummatches, coffhit, range, poly, massstr, massstrinv, massratio,  leftinvexp, middleinvexp, rightinvexp, cleft, cmiddle, cright, search_uses, random_inputs);
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
  int massratio;
  struct timespec t;
  long seed;
  double testrand;
  double r;
  int nummatches;
  int i;
  char execstr[320];
  int coffhit[3];
  int numcoff;
  long exseed;
  long seedsec;
  long seedus;
  int nummult;
  multipliers *multstart;
  multipliers *mult;
  matches *matchstart;
  matches *matchptr;
  double muon_reldiff;
  int leftinvexp;
  int middleinvexp;
  int rightinvexp;
  int numdimensions;
  int invexp1;
  int invexp2;
  int invexp3;
  random_input random_inputs;
  char massstr[20];
  char massstrinv[20];
  double massterm;
  char poly[20];
  double mp;
  input_use search_uses;
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
lepton version %.2f\n", version);
printf("\n\
usage: lepton <seed> <exponentlimit> <phase1filter> <minsymmetry> <maxcomplexity>\n\
\n\
 example: lepton 0 4 5 50 70 (no external seed, max exponent -1/4 or +1/4, phase1 filter of 1x10^-5,\n\
                             minimum phase 2 symmetry of 50, maximum phase 2 complexity of 70. (recommended values).\n\
\n\
 example: lepton 1000000 3 6 7 30 150 (with an external seed, max exponent -1/3 or +1/3, phase1 fitler of 1x10^-6,\n\
                                       minimum phase 2 symmetry of 30, maximum phase 2 complexity of 150).\n\
\n");
printf("\
 *** Warning: very rough work in progress ***\n\
\n\
 This program searches for polynomial-like formulas (polyforms) that might generate the charged lepton masses from simple coefficients.\n\
 These three term (plus constant of -1) formulas use positive and negative rational exponents less than 1 and are designed to give exactly three positive real\n\
 roots representing the lepton masses. Unlike integer exponents, rational exponents less than +/-1 generate positive-only roots with a spectrum similar to that\n\
 of the observed masses using relatively small coefficients on each term.\n\
\n");
printf("\
 Processing is broken into two phases.   Phase 1 starts with the known charged lepton masses, solves a proposed polyform for the three coefficients and searches\n\
 for interesting multipliers to those coefficients. Phase 2 reverses the process by converting the close multiplier matches from phase 1 to exact formulas\n\
 and solving them for the charged lepton masses and other outputs.  These results are then filtered by comparing them to their experimental values and\n\
 uncertainties.\n\
\n");
printf("\
 Phase 1:\n\
   A random combination of three exponents are selected along with a mass ratio (like M/v or M/mZ etc.).   The three charged lepton masses are used as inputs\n\
   along with the other mass used in the selected mass ratio. Any input with experimental uncertainty greater than 1x10^-7 is varied randomly within it's\n\
   uncertainty range for each run to enable usefully tight filter limits in phase1.\n\
\n");
printf("\
   This formula is then solved for the coefficients for each of the three exponent terms (left, middle, right).  Each coefficient is then multiplied by a\n\
   series of test multipliers to see if the result is a relatively small integer rational number.   These interesting (passing the phase 1 filter)\n\
   coefficient * multiplier matches are then stored in a table for phase 2.  Each polyform (set of three exponents) is run through phase 1 for every mass ratio\n\
   and all of the matches for a given polyform  are saved in the same table for phase 2 to enable mixed mass ratio formulas.\n\
\n");
printf("\
 Phase 2:\n\
   formulas are constructed with the same exponents and various interesting multipliers from phase 1 and then solved for the three charged lepton masses and\n\
   other variables.  Unlike in phase 1, in phase 2 each exponent term can have a different mass ratio, as long as the product of the mass ratio and multiplier\n\
   generates the correct individual polynomial term.\n\
\n");
printf("\
   Before solving each proposed formula, all of the experimentally known variables are ranked by relative standard uncertainty and the three with the highest\n\
   uncertainty are used as outputs (solved for) with the rest used as inputs. The resulting outputs are then checked against their experimental ranges.\n\
\n");
printf("\
   Symmetry and complexity scores are also\n\
   assigned to each result based on the numerical values in the multiplier, with higher symmetry and lower complexity generally meaning a simpler formula.\n\
   Minimum symmetry and maximum complexity\n\
   filters are provided to restrict the formulas allowed to be solved.   This can greatly speed up finding the most interesting\n\
   formulas as cpu time is not wasted solving the much more common higher complexity formulas.   The optimal maxcomplexity setting varies by exponent limit\n\
   and the phase1 limit.\n\
\n");
printf("\
 seed            External integer used as part of a seed for srand48() along with second and microsecond clock data.\n\
                 This can help separate threads started at the same time on the same machine have different seeds.\n\
                 Recommended value is <cpu_index>*1000000 where 1 <= cpu_index <= 999.\n\
\n");
printf("\
 exponent limit  Maximum inverse exponent (3 allows for exponents of -1, +1, -1/2, +1/2, -1/3, +1/3).\n\
\n\
 phase1 limit    Threshold for selecting (multiplier * coefficient) for further processing if the result is\n\
                 this close to an integer or simple rational number. (threshold = 1x10^-inlimit)\n\
\n\
 minsymmetry:    Minimum multiplier symmetry allowed to be processed in phase 2.\n\
\n\
 maxcomplexity:  Maximum multiplier complexity allowed to be processed in phase 2.\n\
\n");
    return(1);
  }

  // init pseudorandom number generator from external seed and clock
  clock_gettime(CLOCK_REALTIME, &t);
  seedsec=(t.tv_sec % 1000000000);
  seedus=(t.tv_nsec / 1000); 
  seed=exseed ^ (seedsec + seedus);
  srand48(seed);
  testrand=drand48();
  printf("init, version: %.2f, external seed: %ld, seedsec: %ld, seedus: %ld, seed: %ld, firstrand: %.9e\n", version, exseed, seedsec, seedus, seed, testrand);

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
    mp=kg_to_ev * sqrt(hbar_ref * c_ref / random_inputs.G_sample);

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

    // sort exponents
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

    // init mass, poly terms and strings
    massstr[19]=0;
    massstrinv[19]=0;
    poly[19]=0;
    sprintf(poly, "P%+d%+d%+d", leftinvexp, middleinvexp, rightinvexp);

#ifdef MIXEDMASSMODE 
    matchptr=matchstart;
    nummatches=0;
    for (i=0; i<=2; i++) {
      coffhit[i]=0;
    }
    for (massratio=0; massratio<=4; massratio++) {
      initUses(&search_uses);
      if (massratio == 0) {
        massterm=pow(alpha_ref, 9.00) * mp;
        sprintf(massstr,    "M/(a^9.00)mp");
        sprintf(massstrinv,    "(a^9.00)mp/M");
        search_uses.G=1;
      } else if (massratio == 1) {
        massterm=v_ref / sqrt(2);
        sprintf(massstr,    " sqrt(2)M/v ");
        sprintf(massstrinv,    " v/sqrt(2)M ");
        search_uses.v=1;
      } else if (massratio == 2) {
        massterm=random_inputs.mz_sample;
        sprintf(massstr,    "    M/mZ    ");
        sprintf(massstrinv,    "    mZ/M    ");
        search_uses.mz=1;
      } else if (massratio == 3) {
        massterm=random_inputs.mw_sample;
        sprintf(massstr,    "    M/mW    ");
        sprintf(massstrinv,    "    mW/M    ");
        search_uses.mw=1;
      } else if (massratio == 4) {
        massterm=random_inputs.mh0_sample;
        sprintf(massstr,    "    M/mH0   ");
        sprintf(massstrinv,    "    mH0/M   ");
        search_uses.mh0=1;
      }
      solvePolyforCoefficients(multstart, &nummult, &matchptr, &nummatches, coffhit, random_input_count, random_inputs, massterm, massstr, massstrinv, poly, leftinvexp, middleinvexp, rightinvexp, massratio, range, &search_uses);
    }
    if (nummatches > 0) {
      numcoff=0;
      for (i=0; i<=2; i++) {
        if (coffhit[i] != 0) {
          numcoff++;
        }
      }
      if (numcoff == 3) {
        verifyMatches(matchstart, &nummatches, leftinvexp, middleinvexp, rightinvexp, random_input_count, minsymmetry, maxcomplexity);
      } else {
#ifdef SHOWSTATUS
        printf("status, No complete three-term phase 2 formulas to solve, coffhit: %d%d%d\n", coffhit[0], coffhit[1], coffhit[2]);
        fflush(stdout);
#endif
      }
    } else {
#ifdef SHOWSTATUS
      printf("status, No interesting coefficient multipliers found\n");
      fflush(stdout);
#endif
    }
#else   // not mixedmassmode
    for (massratio=0; massratio<=4; massratio++) {
      initUses(&search_uses);
      if (massratio == 0) {
        massterm=pow(alpha_ref, 9.00) * mp;
        sprintf(massstr,    "M/(a^9.00)mp");
        sprintf(massstrinv,    "(a^9.00)mp/M");
        search_uses.G=1;
      } else if (massratio == 1) {
        massterm=v_ref / sqrt(2);
        sprintf(massstr,    " sqrt(2)M/v ");
        sprintf(massstrinv,    " v/sqrt(2)M ");
        search_uses.v=1;
      } else if (massratio == 2) {
        massterm=random_inputs.mz_sample;
        sprintf(massstr,    "    M/mZ    ");
        sprintf(massstrinv,    "    mZ/M    ");
        search_uses.mz=1;
      } else if (massratio == 3) {
        massterm=random_inputs.mw_sample;
        sprintf(massstr,    "    M/mW    ");
        sprintf(massstrinv,    "    mW/M    ");
        search_uses.mw=1;
      } else if (massratio == 4) {
        massterm=random_inputs.mh0_sample;
        sprintf(massstr,    "    M/mH0   ");
        sprintf(massstrinv,    "    mH0/M   ");
        search_uses.mh0=1;
      }
      matchptr=matchstart;
      nummatches=0;
      for (i=0; i<=2; i++) {
        coffhit[i]=0;
      }
      solvePolyforCoefficients(multstart, &nummult, &matchptr, &nummatches, coffhit, random_input_count, random_inputs, massterm, massstr, massstrinv, poly, leftinvexp, middleinvexp, rightinvexp, massratio, range, &search_uses);
      if (nummatches > 0) {
        numcoff=0;
        for (i=0; i<=2; i++) {
          if (coffhit[i] != 0) {
            numcoff++;
          }
        }
        if (numcoff == 3) {
          verifyMatches(matchstart, &nummatches, leftinvexp, middleinvexp, rightinvexp, random_input_count, minsymmetry, maxcomplexity);
        } else {
#ifdef SHOWSTATUS
          printf("status, No complete three-term phase 2 formulas to solve.\n");
          fflush(stdout);
#endif
        }
      } else {
#ifdef SHOWSTATUS
        printf("status, No interesting coefficient multipliers found, no phase 2 formulas to solve\n");
        fflush(stdout);
#endif
      } // end nummatches >0
    }  // end for massratio
#endif // end mixedmassmode
  } // end while 1
  return(0);
}   
