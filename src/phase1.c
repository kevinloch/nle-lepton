#include <stdio.h>
#include <stdlib.h> // drand
#include <math.h> // pow
#include <time.h>
#include "nle-lepton.h"
#include "reference.h"
#include "cscanner.h"

int solveNLEforCoefficients(multipliers *multstart, int *nummult, matches **matchptr, int *nummatches, int *coffhit, int random_input_count, random_input random_inputs, char *exponents, int leftinvexp, int middleinvexp, int rightinvexp, int range) {
  /***********/
  /* Phase 1 */
  /***********/
  // solve a three term polynomial-like non-linear equation for the unknown coefficients given known roots (particle masses)
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
  int stalledlimit=500000;            // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
  double defaultrangemultiplier=5.0;  // lowest practical range multiplier, fastest for most formulas
  double stalledrangemultiplier=17.0; // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
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
  printf("status, Solving phase 1 formula for coefficients, exponents: %s\n", exponents);
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

            printf("debug, exponents: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9e, best_ordering: %d, precision_last: %.9e, rangefactor: %.9e, i: %d, e_test: %.9e, u_test: %.9e, t_test: %.9e, cleft: %.9e, cmiddle: %.9e, cright: %.9e, cleftrange: %.9e, cmiddlerange: %.9e, crightrange: %.9e\n", exponents, samples, elapsedtime, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], rangefactor, i, e_test, u_test, t_test, cleft_center[ordering], cmiddle_center[ordering], cright_center[ordering], cleft_range[ordering], cmiddle_range[ordering], cright_range[ordering]);
            fflush(stdout);
          }
#endif
          if ((progress[ordering] == ratiograceperiod) || (precision_last[ordering] > stuckprecision)) { // it's stuck, try resetting
#ifdef DEBUG10
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
            printf("debug, exponents: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, resetting\n", exponents, samples, elapsedtime, ordering, progress[ordering]);
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
          printf("debug, exponents: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, stalled\n", exponents, samples, elapsedtime, ordering, progress[ordering]);
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
                  // option 3 - consistently slower than option 2 and sometimes dramatically slows some normally fast formulas
                  cleft_range_new=cleft[ordering] * rangefactor;
                  cleft_range[ordering]=((cleft_range[ordering] + cleft_range_new) / 2.0);
                  cmiddle_range_new=cmiddle[ordering] * rangefactor;
                  cmiddle_range[ordering]=((cmiddle_range[ordering] + cmiddle_range_new) / 2.0);
                  cright_range_new=cright[ordering] * rangefactor;
                  cright_range[ordering]=((cright_range[ordering] + cright_range_new) / 2.0);
*/
/*
                  // option 4 - comparable to option 2 but sometimes slows normaly fast formulas
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
                printf("debug, exponents: %s, samples: %10ld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9e, best_ordering: %d, precision_last: %.9e, rangefactor: %.9e, i: %d, e_test: %.9e, u_test: %.9e, t_test: %.9e, cleft: %.9e, cmiddle: %.9e, cright: %.9e, cleftrange: %.9e, cmiddlerange: %.9e, crightrange: %.9e\n", exponents, samples, elapsedtime, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], rangefactor, i, e_test, u_test, t_test, cleft[ordering], cmiddle[ordering], cright[ordering], cleft_range[ordering], cmiddle_range[ordering], cright_range[ordering]);
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
  printf("status, Solved  phase 1 formula for coefficients, random input: %i, exponents:  %s, Tau: %.9e, samples: %ld, ordering: %d, precision: %.3e (%6.4fs)\n", random_input_count, exponents, random_inputs.tau_sample, samples, best_ordering, precision, elapsedtime);
  fflush(stdout);
#endif
#ifdef SHOWSTATUS
  printf("status, +-----------------+------------+------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
  printf("status, |    Tau mass     | Exponents  | Mass ratio | NLE term |   Coefficient   | C * NLEterm(e)  | C * NLEterm(u)  | C * NLEterm(t)  |\n");
  printf("status, +-----------------+------------+------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
  printf("status, | %.9e | %10s |    M/me    |  +left   | %.9e | %.9e | %.9e | %.9e |\n", random_inputs.tau_sample, exponents, cleft_center[best_ordering], e_test_1, u_test_1, t_test_1);
  printf("status, | %.9e | %10s |    M/me    |  -middle | %.9e | %.9e | %.9e | %.9e |\n", random_inputs.tau_sample, exponents, cmiddle_center[best_ordering], e_test_2, u_test_2, t_test_2);
  printf("status, | %.9e | %10s |    M/me    |  +right  | %.9e | %.9e | %.9e | %.9e |\n", random_inputs.tau_sample, exponents, cright_center[best_ordering], e_test_3, u_test_3, t_test_3);
  printf("status, +-----------------+------------+------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
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
  cscanner(multstart, nummult, matchptr, nummatches, coffhit, range, exponents, leftinvexp, middleinvexp, rightinvexp, cleft_center[best_ordering], cmiddle_center[best_ordering], cright_center[best_ordering], random_inputs);
  nummatchend=*nummatches;
  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
#ifdef SHOWSTATUS
  printf("status, Found %d interesting coefficient multipliers (%6.4fs)\n", (nummatchend-nummatchstart), elapsedtime);
  fflush(stdout);
#endif

  return(0);
}
