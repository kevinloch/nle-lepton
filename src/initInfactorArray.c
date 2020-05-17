#include <stdio.h>
#include <math.h> // pow, M_PI
#include <stdlib.h> // abs
#include <time.h>
#include "nle-lepton.h"
#include "util.h"

//#define DEBUG_INFACTOR

void initInfactorArray(nle_config_t *nle_config, nle_state_t *nle_state) {
  // store pre-computed static multiplier terms that go inside radical

  struct timespec start_time;
  struct timespec end_time;
  double elapsed_time;

  int upin=1, downin=1;
  int piupin=0, pidownin=1;
  int aupin=0, adownin=1;
  int e2upin=0, e2downin=2;
  int nbvupin=0;
  int nssupin=0;
  int userupin=0, userdownin=1;
  double updownin;
  double piin;
  double ain;
  double e2in;
  double nbv[27];
  double nss[27];
  double userin;
  double infactor;
  int upcomplexity, downcomplexity;

  unsigned int u, v;
  int i;
  nle_infactor_precomputed_t *multiplier;

  clock_gettime(CLOCK_REALTIME, &start_time);

  // nball volume (index is term exponent and n-ball dimension)
  // nbv[n] = nss[n-1]/n (n > 0)
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
  // nss[n] = 2pi * nbv[n-1] (n > 0)
  // nss[n] = 2pi * nss[n-2]/(n-1) (n > 1)
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

  multiplier=nle_state->infactors_precomputed_start;
  nle_state->infactors_precomputed_count=0;
  for (upin=1; upin <= nle_config->infactor_rational_max; upin++) {
    for (downin=1; downin <= nle_config->infactor_rational_max; downin++) {
      if ((!nle_config->infactor_rational_filter) || (upin <= 4) || ((downin == 1) && ((upin == 8) || (upin == 9) || (upin == 12) || (upin == 16) || (upin == 18) || (upin == 24) || (upin == 27) || (upin == 32)))) {
        if ((!nle_config->infactor_rational_filter) || (downin <= 4) || ((upin == 1) && ((downin == 8) || (downin == 9) || (downin == 12) || (downin == 16) || (downin == 18) || (downin == 24) || (downin == 27) || (downin == 32)))) { 
          u=upin;
          v=downin;
          if (gcd(u, v) == 1) {
            updownin=(double)upin / (double)downin;

            for (e2upin=-nle_config->infactor_2_exp_up_max; e2upin <= nle_config->infactor_2_exp_up_max; e2upin++) {
              for (e2downin=1; e2downin <= nle_config->infactor_2_exp_down_max; e2downin++) {
                u=abs(e2upin);
                v=e2downin;
                if (!((u == 1) && (v == 1)) && (gcd(u, v) == 1)) {  // extra checks to prevent 2 or 1/2
                  e2in=pow(2.0, ((double)e2upin / (double)e2downin));

                  for (aupin=-nle_config->infactor_alpha_exp_up_max; aupin <= nle_config->infactor_alpha_exp_up_max; aupin++) {
                    for (adownin=1; adownin <= nle_config->infactor_alpha_exp_down_max; adownin++) {
                      u=abs(aupin);
                      v=adownin;
                      if (gcd(u, v) == 1) {
                        ain=pow(nle_config->ref_alpha_em, ((double)aupin / (double)adownin));

                        for (piupin=-nle_config->infactor_pi_exp_up_max; piupin <= nle_config->infactor_pi_exp_up_max; piupin++) {
                          for (pidownin=1; pidownin <= nle_config->infactor_pi_exp_down_max; pidownin++) {
                            u=abs(piupin);
                            v=pidownin;
                            if (gcd(u, v) == 1) {
                              piin=pow(M_PI, ((double)piupin / (double)pidownin));
 
                              for (nssupin=-nle_config->infactor_nss_enable; nssupin <= nle_config->infactor_nss_enable; nssupin++) {
                                for (nbvupin=-nle_config->infactor_nbv_enable; nbvupin <= nle_config->infactor_nbv_enable; nbvupin++) {
                                  if (((nbvupin == 0) && (nssupin == 0)) || ((nbvupin != 0) && (nssupin == 0)) || ((nbvupin == 0) && (nssupin != 0))) {

                                    for (userupin=-nle_config->infactor_user_exp_up_max; userupin <= nle_config->infactor_user_exp_up_max; userupin++) {
                                      for (userdownin=1; userdownin <= nle_config->infactor_user_exp_down_max; userdownin++) {
                                        u=abs(userupin);
                                        v=userdownin;
                                        if (gcd(u, v) == 1) {
                                          userin=pow(nle_config->infactor_user, ((double)userupin / (double)userdownin));
                                          infactor=updownin * e2in * ain * piin * userin;
#ifdef DEBUG_INFACTOR
                                          printf("debug, infactorInit, up: %d, down: %d, e2up: %d, e2down: %d, aup: %d, adown: %d, piup: %d, pidown: %d, nssup: %d, nbvup: %d, userup: %d, userdown: %d, infactor: %.9e, updown: %.9e, e2: %.9e, a: %.9e, pi: %.9e, user: %.9e\n", upin, downin, e2upin, e2downin, aupin, adownin, piupin, pidownin, nssupin, nbvupin, userupin, userdownin, infactor, updownin, e2in, ain, piin, userin);
                                          fflush(stdout);
#endif
                                          multiplier->infactor_rational_up=upin;
                                          multiplier->infactor_rational_down=downin;
                                          multiplier->infactor_2_exp_up=e2upin;
                                          multiplier->infactor_2_exp_down=e2downin;
                                          multiplier->infactor_alpha_exp_up=aupin;
                                          multiplier->infactor_alpha_exp_down=adownin;
                                          multiplier->infactor_pi_exp_up=piupin;
                                          multiplier->infactor_pi_exp_down=pidownin;
                                          multiplier->infactor_nss=nssupin;
                                          multiplier->infactor_nbv=nbvupin;
                                          multiplier->infactor_user_exp_up=userupin;
                                          multiplier->infactor_user_exp_down=userdownin;
                                          // ignore 1 on rationals
                                          if (upin == 1) {
                                            upcomplexity=0;
                                          } else {
                                            upcomplexity=upin;
                                          }
                                          if (downin == 1) {
                                            downcomplexity=0;
                                          } else {
                                            downcomplexity=downin;
                                          }
                                          multiplier->infactor_complexity=\
                                                 upcomplexity + downcomplexity\
                                                 + abs(e2upin) + (e2downin-1)\
                                                 + abs(aupin) + (adownin-1)\
                                                 + abs(piupin) + (pidownin-1)\
                                                 + abs(nssupin)\
                                                 + abs(nbvupin)\
                                                 + abs(userupin) + (userdownin-1);
                                          for (i=1; i<= nle_config->exp_inv_max; i++) {
                                            if (i == 1) {
                                              multiplier->infactor_multiplier[1]=infactor * pow(nbv[1], (double)nbvupin) * pow(nss[1], (double)nssupin);
                                            } else if (i == 2) {
                                              multiplier->infactor_multiplier[2]=pow(infactor * pow(nbv[2], (double)nbvupin) * pow(nss[2], (double)nssupin), (1.0 / 2.0));
                                            } else if (i == 3) {
                                              multiplier->infactor_multiplier[3]=pow(infactor * pow(nbv[3], (double)nbvupin) * pow(nss[3], (double)nssupin), (1.0 / 3.0));
                                            } else if (i == 4) {
                                              multiplier->infactor_multiplier[4]=pow(infactor * pow(nbv[4], (double)nbvupin) * pow(nss[4], (double)nssupin), (1.0 / 4.0));
                                            } else if (i == 5) {
                                              multiplier->infactor_multiplier[5]=pow(infactor * pow(nbv[5], (double)nbvupin) * pow(nss[5], (double)nssupin), (1.0 / 5.0));
                                            } else if (i == 6) {
                                              multiplier->infactor_multiplier[6]=pow(infactor * pow(nbv[6], (double)nbvupin) * pow(nss[6], (double)nssupin), (1.0 / 6.0));
                                            } else if (i == 7) {
                                              multiplier->infactor_multiplier[7]=pow(infactor * pow(nbv[7], (double)nbvupin) * pow(nss[7], (double)nssupin), (1.0 / 7.0));
                                            } else if (i == 8) {
                                              multiplier->infactor_multiplier[8]=pow(infactor * pow(nbv[8], (double)nbvupin) * pow(nss[8], (double)nssupin), (1.0 / 8.0));
                                            } else if (i == 9) {
                                              multiplier->infactor_multiplier[9]=pow(infactor * pow(nbv[9], (double)nbvupin) * pow(nss[9], (double)nssupin), (1.0 / 9.0));
                                            } else if (i == 10) {
                                              multiplier->infactor_multiplier[10]=pow(infactor * pow(nbv[10], (double)nbvupin) * pow(nss[10], (double)nssupin), (1.0 / 10.0));
                                            } else if (i == 11) {
                                              multiplier->infactor_multiplier[11]=pow(infactor * pow(nbv[11], (double)nbvupin) * pow(nss[11], (double)nssupin), (1.0 / 11.0));
                                            } else if (i == 12) {
                                              multiplier->infactor_multiplier[12]=pow(infactor * pow(nbv[12], (double)nbvupin) * pow(nss[12], (double)nssupin), (1.0 / 12.0));
                                            } else if (i == 13) {
                                              multiplier->infactor_multiplier[13]=pow(infactor * pow(nbv[13], (double)nbvupin) * pow(nss[13], (double)nssupin), (1.0 / 13.0));
                                            } else if (i == 14) {
                                              multiplier->infactor_multiplier[14]=pow(infactor * pow(nbv[14], (double)nbvupin) * pow(nss[14], (double)nssupin), (1.0 / 14.0));
                                            } else if (i == 15) {
                                              multiplier->infactor_multiplier[15]=pow(infactor * pow(nbv[15], (double)nbvupin) * pow(nss[15], (double)nssupin), (1.0 / 15.0));
                                            } else if (i == 16) {
                                              multiplier->infactor_multiplier[16]=pow(infactor * pow(nbv[16], (double)nbvupin) * pow(nss[16], (double)nssupin), (1.0 / 16.0));
                                            } else if (i == 17) {
                                              multiplier->infactor_multiplier[17]=pow(infactor * pow(nbv[17], (double)nbvupin) * pow(nss[17], (double)nssupin), (1.0 / 17.0));
                                            } else if (i == 18) {
                                              multiplier->infactor_multiplier[18]=pow(infactor * pow(nbv[18], (double)nbvupin) * pow(nss[18], (double)nssupin), (1.0 / 18.0));
                                            } else if (i == 19) {
                                              multiplier->infactor_multiplier[19]=pow(infactor * pow(nbv[19], (double)nbvupin) * pow(nss[19], (double)nssupin), (1.0 / 19.0));
                                            } else if (i == 20) {
                                              multiplier->infactor_multiplier[20]=pow(infactor * pow(nbv[20], (double)nbvupin) * pow(nss[20], (double)nssupin), (1.0 / 20.0));
                                            } else if (i == 21) {
                                              multiplier->infactor_multiplier[21]=pow(infactor * pow(nbv[21], (double)nbvupin) * pow(nss[21], (double)nssupin), (1.0 / 21.0));
                                            } else if (i == 22) {
                                              multiplier->infactor_multiplier[22]=pow(infactor * pow(nbv[22], (double)nbvupin) * pow(nss[22], (double)nssupin), (1.0 / 22.0));
                                            } else if (i == 23) {
                                              multiplier->infactor_multiplier[23]=pow(infactor * pow(nbv[23], (double)nbvupin) * pow(nss[23], (double)nssupin), (1.0 / 23.0));
                                            } else if (i == 24) {
                                              multiplier->infactor_multiplier[24]=pow(infactor * pow(nbv[24], (double)nbvupin) * pow(nss[24], (double)nssupin), (1.0 / 24.0));
                                            } else if (i == 25) {
                                              multiplier->infactor_multiplier[25]=pow(infactor * pow(nbv[25], (double)nbvupin) * pow(nss[25], (double)nssupin), (1.0 / 25.0));
                                            } else if (i == 26) {
                                              multiplier->infactor_multiplier[26]=pow(infactor * pow(nbv[26], (double)nbvupin) * pow(nss[26], (double)nssupin), (1.0 / 26.0));
                                            }  // if i
#ifdef DEBUG_INFACTOR
                                            printf("debug, infactorInit, i: %d, multiplier: %.9e\n", i, multiplier->infactor_multiplier[i]);
                                            fflush(stdout);
#endif
                                          } // for i
                                          initUses(&multiplier->infactor_uses);
                                          if (aupin != 0) {
                                            multiplier->infactor_uses.alpha_em=1;
                                          }
                                          nle_state->infactors_precomputed_count++;
#ifdef DEBUG_INFACTOR
                                          printf("debug, infactorInit, count: %d\n", nle_state->infactors_precomputed_count);
                                          fflush(stdout);
#endif
                                          multiplier++;
                                        } // gcd user
                                      } // down user
                                    } // up user
                                  } // not both nbvupin and nssupin
                                } // nbv
                              } // nss
                            } // gcd pi
                          } // pidown
                        } // piup
                      } // gcd a
                    } // adown
                  } // aup
                } // gcd e2
              } // e2down
            } // e2up
          } // gcd updown
        } // downfilter
      } // upfilter
    } // downin
  } // upin

  clock_gettime(CLOCK_REALTIME, &end_time);
  elapsed_time=((double)(end_time.tv_sec - 1500000000) + ((double)end_time.tv_nsec / 1.0E9)) - ((double)(start_time.tv_sec - 1500000000) + ((double)start_time.tv_nsec) / 1.0E9);

  printf("init, Initialized %d pre-computed multipliers for static factors inside radical (%6.4fs)\n", nle_state->infactors_precomputed_count, elapsed_time);
  fflush(stdout);
}
