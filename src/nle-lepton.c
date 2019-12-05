/*
 * Copyright (c) 2019, Kevin M. Loch
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nle-lepton.h"
#include "reference.h"
#include "usage.h"
#include "initMultiplierArray.h"
#include "phase1.h"
#include "verifyMatches.h"

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
  char exponents[20];
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
    printUsage();
    return(1);
  }

  // init pseudorandom number generator from external seed and clock
  clock_gettime(CLOCK_REALTIME, &t);
  seedsec=(t.tv_sec % 1000000000);
  seedus=(t.tv_nsec / 1000); 
  seed=exseed ^ (seedsec + seedus);
  srand48(seed);
  testrand=drand48();
  printf("init, version: %s, external seed: %ld, seedsec: %ld, seedus: %ld, seed: %ld, firstrand: %.9e\n", NLE_VERSION, exseed, seedsec, seedus, seed, testrand);

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
    // select random exponents
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
    // select random exponents
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
    // force exponents for debugging or focused searches
    leftinvexp=+1;
    middleinvexp=+2;
    rightinvexp=+3;
*/
/*
    // force exponents for debugging or focused searches
    leftinvexp=-18;
    middleinvexp=+17;
    rightinvexp=+18;
*/

    // init strings
    exponents[19]=0;
    sprintf(exponents, "E%+d%+d%+d", leftinvexp, middleinvexp, rightinvexp);
    matchptr=matchstart;
    nummatches=0;
    for (i=0; i<=2; i++) {
      coffhit[i]=0;
    }

    // phase 1
    solveNLEforCoefficients(multstart, &nummult, &matchptr, &nummatches, coffhit, random_input_count, random_inputs, exponents, leftinvexp, middleinvexp, rightinvexp, range);
    if (nummatches > 0) {
      numcoff=0;
      for (i=0; i<=2; i++) {
        if (coffhit[i] != 0) {
          numcoff++;
        }
      }
      if (numcoff == 3) {
        // phase 2
        verifyMatches(matchstart, &nummatches, exponents, leftinvexp, middleinvexp, rightinvexp, random_input_count, minsymmetry, maxcomplexity);
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
