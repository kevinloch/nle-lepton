#include <stdio.h>
#include <stdlib.h> // drand
#include <math.h> // pow
#include <time.h>
#include "nle-lepton.h"
#include "cscanner.h"

//#define DEBUG10
//#define DEBUG11
//#define DEBUG12

void solveNLEforCoefficients(nle_config_t *nle_config, nle_state_t *nle_state) {
  // solve a three term polynomial-like non-linear equation for the unknown coefficients given known roots (particle masses)
  int i;
  long long samples=0;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  int matches_count_start;
  int matches_count_end;
  int good_coefficients;
  int c2_gt_c1, c2_gt_c3, c1_gt_c3;

  clock_gettime(CLOCK_REALTIME, &starttime);

  //  mc test vars
  double r;
  double sm1_test_term1=0;
  double sm1_test_term2=0;
  double sm1_test_term3=0;
  double sm1_test=0;
  double sm2_test_term1=0;
  double sm2_test_term2=0;
  double sm2_test_term3=0;
  double sm2_test=0;
  double sm3_test_term1=0;
  double sm3_test_term2=0;
  double sm3_test_term3=0;
  double sm3_test=0;
  double worst_test;
  double range_factor;
  double term1_mass_sm1=0;
  double term2_mass_sm1=0;
  double term3_mass_sm1=0;
  double term1_mass_sm2=0;
  double term2_mass_sm2=0;
  double term3_mass_sm2=0;
  double term1_mass_sm3=0;
  double term2_mass_sm3=0;
  double term3_mass_sm3=0;
  double precision=1.0E99;
  double range_multiplier[6];

  // these tunings affect speed and reliability, adjust with extreme care
  double precision_target=1.0E-11;     // solve NLE to this level of precision
  double test_ratio=25.0;              // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
  int ratio_grace_period=25;           // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
  long long stalled_limit=500000;      // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
  double defaultrange_multiplier=5.0;  // lowest practical range multiplier, fastest for most formulas
  double stalledrange_multiplier=17.0; // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
  int slowcheckpoint=1000000;          // progress point to check on slow processes
  double stuckprecision=1.0E-2;        // if precision is not past this level by slowcheckpoint, try resetting

  // mc outputs
  int progress[6];
  int stalled[6];
  int ordering;
  int best_ordering;
  double best_precision_last;
  double precision_last[6];
  double c1[6];
  double c1_center[6];
  double c1_range[6];
  double c1_range_new;
  double c2[6];
  double c2_center[6];
  double c2_range[6];
  double c2_range_new;
  double c3[6];
  double c3_center[6];
  double c3_range[6];
  double c3_range_new;

  // starting v4.0 denominator is always v for this step, other mass ratio factors are now swapped out in cscanner()
  term1_mass_sm1=  pow(((double)nle_state->random_sample_sm1 / nle_state->random_sample_v), (1.0 / (double)nle_state->term1.exp_inv));
  term1_mass_sm2=  pow(((double)nle_state->random_sample_sm2 / nle_state->random_sample_v), (1.0 / (double)nle_state->term1.exp_inv));
  term1_mass_sm3=  pow(((double)nle_state->random_sample_sm3 / nle_state->random_sample_v), (1.0 / (double)nle_state->term1.exp_inv));
  term2_mass_sm1=pow(((double)nle_state->random_sample_sm1 / nle_state->random_sample_v), (1.0 / (double)nle_state->term2.exp_inv));
  term2_mass_sm2=pow(((double)nle_state->random_sample_sm2 / nle_state->random_sample_v), (1.0 / (double)nle_state->term2.exp_inv));
  term2_mass_sm3=pow(((double)nle_state->random_sample_sm3 / nle_state->random_sample_v), (1.0 / (double)nle_state->term2.exp_inv));
  term3_mass_sm1= pow(((double)nle_state->random_sample_sm1 / nle_state->random_sample_v), (1.0 / (double)nle_state->term3.exp_inv));
  term3_mass_sm2= pow(((double)nle_state->random_sample_sm2 / nle_state->random_sample_v), (1.0 / (double)nle_state->term3.exp_inv));
  term3_mass_sm3= pow(((double)nle_state->random_sample_sm3 / nle_state->random_sample_v), (1.0 / (double)nle_state->term3.exp_inv));

  if (nle_config->status_enable == 1) {
    printf("status, Solving phase 1 formula for coefficients, exponents: %s\n", nle_state->exponents_str);
    fflush(stdout);
  }

  //  solve formula for coefficients
  best_precision_last=1.0E99;
  while (best_precision_last > precision_target) {
    //  init outputs
    for (ordering=0; ordering<=5; ordering++) {
      c1[ordering]=1.0E0;
      c1_center[ordering]=1.0E0;
      c1_range[ordering]=0.5E0;
      c2[ordering]=1.0E0;
      c2_center[ordering]=1.0E0;
      c2_range[ordering]=0.5E0;
      c3[ordering]=1.0E0;
      c3_center[ordering]=1.0E0;
      c3_range[ordering]=0.5E0;
      precision_last[ordering]=1.0E99;
      progress[ordering]=0;
      stalled[ordering]=0;
      range_multiplier[ordering]=defaultrange_multiplier;
    }
    best_precision_last=1.0E99;
    best_ordering=-1;
    i=0;
    for (samples=0; best_precision_last > precision_target; samples++) {
      for (ordering=0; ordering<=5; ordering++) {
        if ((samples > 1) && ((samples % slowcheckpoint) == 0)) { // check on slow processes
#ifdef DEBUG10
          if ((samples % 10000000) == 0) { // rate limit periodic debug prints
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);

            printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9e, best_ordering: %d, precision_last: %.9e, range_factor: %.9e, i: %d, sm1_test: %.9e, sm2_test: %.9e, sm3_test: %.9e, c1: %.9e, c2: %.9e, c3: %.9e, c1range: %.9e, c2range: %.9e, c3range: %.9e\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], range_factor, i, sm1_test, sm2_test, sm3_test, c1_center[ordering], c2_center[ordering], c3_center[ordering], c1_range[ordering], c2_range[ordering], c3_range[ordering]);
            fflush(stdout);
          }
#endif
          if ((progress[ordering] == ratio_grace_period) || (precision_last[ordering] > stuckprecision)) { // it's stuck, try resetting
#ifdef DEBUG11
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
            printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, resetting\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering]);
            fflush(stdout);
#endif
            c1[ordering]=1.0E0;
            c1_center[ordering]=1.0E0;
            c1_range[ordering]=0.5E0;
            c2[ordering]=1.0E0;
            c2_center[ordering]=1.0E0;
            c2_range[ordering]=0.5E0;
            c3[ordering]=1.0E0;
            c3_center[ordering]=1.0E0;
            c3_range[ordering]=0.5E0;
            precision_last[ordering]=1.0E99;
            progress[ordering]=0;
            stalled[ordering]=0;
            range_multiplier[ordering]=defaultrange_multiplier;
            if (ordering == best_ordering) {
              best_precision_last=1.0E99;
            }
          }
        }

        if ((best_precision_last > 1.0E-7) || (ordering == best_ordering)) { // skip other ordings if one is far enough along
          if (ordering == 0) {
            c2_gt_c1=1;
            c2_gt_c3=1;
            c1_gt_c3=1;
          } else if (ordering == 1) {
            c2_gt_c1=1;
            c2_gt_c3=1;
            c1_gt_c3=0;
          } else if (ordering == 2) {
            c2_gt_c1=1;
            c2_gt_c3=0;
            c1_gt_c3=0;
          } else if (ordering == 3) {
            c2_gt_c1=0;
            c2_gt_c3=1;
            c1_gt_c3=1;
          } else if (ordering == 4) {
            c2_gt_c1=0;
            c2_gt_c3=0;
            c1_gt_c3=1;
          } else if (ordering == 5) {
            c2_gt_c1=0;
            c2_gt_c3=0;
            c1_gt_c3=0;
          } // end if ordering
          if (stalled[ordering] == stalled_limit) {
            range_multiplier[ordering]=stalledrange_multiplier; // may be a slow solution, try bigger multiplier
            //stalled[ordering]=0;
#ifdef DEBUG11
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
            printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, stalled\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering]);
            fflush(stdout);
#endif
          }

          // generate random coefficients within limits
          stalled[ordering]++;
          i=0;
          good_coefficients=0;
          while ((i < 200) && (good_coefficients == 0)) {
            i++;
            c1[ordering]=0.0;
            c2[ordering]=0.0;
            c3[ordering]=0.0;
            while (c1[ordering] <= 0.0) {
              r=(double)drand48();
              c1[ordering]=((c1_center[ordering] - c1_range[ordering]) + (r * 2.0 * c1_range[ordering]));
            }
            while (c2[ordering] <= 0.0) {
              r=(double)drand48();
              c2[ordering]=((c2_center[ordering] - c2_range[ordering]) + (r * 2.0 * c2_range[ordering]));
            }                 
            while (c3[ordering] <= 0.0) {
              r=(double)drand48();
              c3[ordering]=((c3_center[ordering] - c3_range[ordering]) + (r * 2.0 * c3_range[ordering]));
            }

            // verify coefficient relationships
            good_coefficients=1;
            if (c2_gt_c1 == 1) {
              if (c2[ordering] < c1[ordering]) {
                good_coefficients=0;
              }
            } else {
              if (c2[ordering] > c1[ordering]) {
                good_coefficients=0;
              }
            }
            if (c2_gt_c3 == 1) {
              if (c2[ordering] < c3[ordering]) {
                good_coefficients=0;
              }
            } else {
              if (c2[ordering] > c3[ordering]) {
                good_coefficients=0;
              }
            }
            if (c1_gt_c3 == 1) {
              if (c1[ordering] < c3[ordering]) {
                good_coefficients=0;
              }
            } else {
              if (c1[ordering] > c3[ordering]) {
                good_coefficients=0;
              }
            } // end if c2_gt_c1
          }
#ifdef DEBUG11
          if (i > 100) { 
            printf("debug, i: %d\n", i);
          }
#endif
          if (good_coefficients == 1) {
            sm1_test_term1=c1[ordering] * term1_mass_sm1;
            sm1_test_term2=c2[ordering] * term2_mass_sm1;
            sm1_test_term3=c3[ordering] * term3_mass_sm1;
  
            sm2_test_term1=c1[ordering] * term1_mass_sm2;
            sm2_test_term2=c2[ordering] * term2_mass_sm2;
            sm2_test_term3=c3[ordering] * term3_mass_sm2;

            sm3_test_term1=c1[ordering] * term1_mass_sm3;
            sm3_test_term2=c2[ordering] * term2_mass_sm3;
            sm3_test_term3=c3[ordering] * term3_mass_sm3;

            sm1_test=sm1_test_term1 - sm1_test_term2 + sm1_test_term3 - 1.0;
            sm2_test=sm2_test_term1 - sm2_test_term2 + sm2_test_term3 - 1.0;
            if ((progress[ordering] < ratio_grace_period) || (((fabs(sm1_test) / fabs(sm2_test)) < test_ratio) && ((fabs(sm2_test) / fabs(sm1_test)) < test_ratio))) {
              sm3_test=sm3_test_term1 - sm3_test_term2 + sm3_test_term3 - 1.0;
                if ((progress[ordering] < ratio_grace_period) || (((fabs(sm1_test) / fabs(sm3_test)) < test_ratio) && ((fabs(sm3_test) / fabs(sm1_test)) < test_ratio) &&\
                                                                ((fabs(sm2_test) / fabs(sm3_test)) < test_ratio) && ((fabs(sm3_test) / fabs(sm2_test)) < test_ratio))) {
#ifdef DEBUG12
                printf("debug, etest: %.3e, utest: %.3e, ttest: %.3e, sm1_test_term1: %.3e, sm1_test_term2: %.3e, sm1_test_term3: %.3e, sm2_test_term1: %.3e, sm2_test_term2: %.3e, sm2_test_term3: %.3e, sm3_test_term1: %.3e, sm3_test_term2: %.3e, sm3_test_term3: %.3e\n", sm1_test, sm2_test, sm3_test, sm1_test_term1, sm1_test_term2, sm1_test_term3, sm2_test_term1, sm2_test_term2, sm2_test_term3, sm3_test_term1, sm3_test_term2, sm3_test_term3);
                fflush(stdout);
#endif

                precision=(fabs(sm1_test) + fabs(sm2_test) + fabs(sm3_test));
                if (precision < precision_last[ordering]) {
                  progress[ordering]++;
                  stalled[ordering]=0;

                  // update best_precision_last and best_ordering
                  precision_last[ordering]=precision;
                  if (precision_last[ordering] < best_precision_last) {
                    best_precision_last=precision;
                    best_ordering=ordering;
                  }

                  // determine new search range for each coefficient
                  if (fabs(sm1_test) > fabs(sm2_test)) {
                    worst_test=fabs(sm1_test);
                  } else {
                    worst_test=fabs(sm2_test);
                  }
                  if (fabs(sm3_test) > worst_test) {
                    worst_test=fabs(sm3_test);
                  }
                  range_factor=worst_test * range_multiplier[ordering];
                  if (range_factor > 1.0) {
                    c1_range[ordering]=c1[ordering] / 2.0;
                    c2_range[ordering]=c2[ordering] / 2.0;
                    c3_range[ordering]=c3[ordering] / 2.0;
                  } else {
                    c1_range_new=c1[ordering] * range_factor;
                    c1_range[ordering]=((c1_range[ordering] + c1_range_new + c1_range_new) / 3.0);
                    c2_range_new=c2[ordering] * range_factor;
                    c2_range[ordering]=((c2_range[ordering] + c2_range_new + c2_range_new) / 3.0);
                    c3_range_new=c3[ordering] * range_factor;
                    c3_range[ordering]=((c3_range[ordering] + c3_range_new + c3_range_new) / 3.0);
                  }
                  c1_center[ordering]=c1[ordering];
                  c2_center[ordering]=c2[ordering];
                  c3_center[ordering]=c3[ordering];

#ifdef DEBUG11
                  clock_gettime(CLOCK_REALTIME, &endtime);
                  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                  printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9e, best_ordering: %d, precision_last: %.9e, range_factor: %.9e, i: %d, sm1_test: %.9e, sm2_test: %.9e, sm3_test: %.9e, c1: %.9e, c2: %.9e, c3: %.9e, c1range: %.9e, c2range: %.9e, c3range: %.9e\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], range_factor, i, sm1_test, sm2_test, sm3_test, c1[ordering], c2[ordering], c3[ordering], c1_range[ordering], c2_range[ordering], c3_range[ordering]);
                  fflush(stdout);
#endif
                } // end if precision
              } // end ttest ratios
            } // end utest ratios
          } // end if good_coefficients
        } // end if best
      } // end for ordering
    }  // end for samples
  }  // end while precision_last

  clock_gettime(CLOCK_REALTIME, &endtime);
  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);

  if (nle_config->status_enable ==1) {
    printf("status, Solved  phase 1 formula for coefficients, random input: %i, exponents:  %s, sm3: %.9e, samples: %lld, ordering: %d, precision: %.3e (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->random_sample_sm3, samples, best_ordering, precision, elapsed_time);
    printf("status, +------------+------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
    printf("status, | Exponents  | Mass ratio | NLE term |   Coefficient   |  C * term(sm1)  |  C * term(sm2)  |  C * term(sm3)  |\n");
    printf("status, +------------+------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
    printf("status, | %10s |    M/v     |  term1   | %.9e | %.9e | %.9e | %.9e |\n", nle_state->exponents_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
    printf("status, | %10s |    M/v     |  term2   | %.9e | %.9e | %.9e | %.9e |\n", nle_state->exponents_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
    printf("status, | %10s |    M/v     |  term3   | %.9e | %.9e | %.9e | %.9e |\n", nle_state->exponents_str, c3_center[best_ordering], sm1_test_term3, sm2_test_term3, sm3_test_term3);
    printf("status, +------------+------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
    fflush(stdout);

    // for debugging very long phase 1 solutions
    //exit(0);

    printf("status, Scanning for coefficient multipliers that match interesting integer or simple rational numbers and mass ratio: ");
    fflush(stdout);
  }
  nle_state->term1.coefficient=c1_center[best_ordering];
  nle_state->term2.coefficient=c2_center[best_ordering];
  nle_state->term3.coefficient=c3_center[best_ordering];
  matches_count_start=nle_state->phase1_matches_count;
  clock_gettime(CLOCK_REALTIME, &starttime);
  cscanner(nle_config, nle_state);
  clock_gettime(CLOCK_REALTIME, &endtime);
  matches_count_end=nle_state->phase1_matches_count;
  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
  if (nle_config->status_enable ==1) {
    printf("status, Found %d interesting coefficient multipliers (%6.4fs)\n", (matches_count_end-matches_count_start), elapsed_time);
    fflush(stdout);
  }
}
