#include <stdio.h>
#include <math.h> // pow, M_PI
#include <stdlib.h> // abs
#include <time.h>
#include "nle-lepton.h"
#include "reference.h"
#include "util.h"

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
