#include <stdio.h>
#include <stdlib.h> // drand
#include <math.h> // powl
#include <time.h>
#include "nle-lepton.h"
#include "util.h"
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
  long double r;
  long double sm1_test_term1=0;
  long double sm1_test_term2=0;
  long double sm1_test_term3=0;
  long double sm1_test=0;
  long double sm2_test_term1=0;
  long double sm2_test_term2=0;
  long double sm2_test_term3=0;
  long double sm2_test=0;
  long double sm3_test_term1=0;
  long double sm3_test_term2=0;
  long double sm3_test_term3=0;
  long double sm3_test=0;
  long double worst_test;
  long double range_factor;
  long double term1_mass_sm1=0;
  long double term2_mass_sm1=0;
  long double term3_mass_sm1=0;
  long double term1_mass_sm2=0;
  long double term2_mass_sm2=0;
  long double term3_mass_sm2=0;
  long double term1_mass_sm3=0;
  long double term2_mass_sm3=0;
  long double term3_mass_sm3=0;
  long double precision=1.0E99;
  long double range_multiplier[6];
  long double dynamicrange_c[6];
  long double dynamicrange_sm1[6];
  long double dynamicrange_sm2[6];
  long double dynamicrange_sm3[6];
  long double dr_high;
  long double dr_low;
  int dr_exception;
  int unsolveable_exception;
  int active_ordering_count;

  // mc outputs
  int progress[6];
  int stalled[6];
  int ordering;
  int best_ordering;
  long double best_precision_last;
  long double precision_last[6];
  long double c1[6];
  long double c1_center[6];
  long double c1_range[6];
  long double c1_range_new;
  long double c2[6];
  long double c2_center[6];
  long double c2_range[6];
  long double c2_range_new;
  long double c3[6];
  long double c3_center[6];
  long double c3_range[6];
  long double c3_range_new;
  long double two_term_test=0.0;

  // tuneables
  long double precision_target;
  long double test_ratio;
  int ratio_grace_period;
  long long stalled_limit;
  long double defaultrange_multiplier;
  long double stalledrange_multiplier;
  int slowcheckpoint;
  long double stuckprecision;
  long double dr_exception_limit=0.0;
  long double dr_grace_period=0;
  int ratio_checkpoint=0;

  // these tunings affect speed and reliability, adjust with extreme care
  if (nle_config->smrfactor_1minus_enable == 1) {
    // 2-term mode with 1-smr
    precision_target=1.0E-15;     // solve NLE to this level of precision
    test_ratio=25.0;              // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=3;         // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;         // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    defaultrange_multiplier=3.0;  // lowest practical range multiplier, fastest for most formulas
    stalledrange_multiplier=3.0;  // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    slowcheckpoint=10000;         // progress point to check on slow processes
    stuckprecision=1.0E+30;       // if precision is not past this level by slowcheckpoint, try resetting
    dr_exception_limit=1.0E+14;   // Mamimum allowed dynamic range of any term (floating point operation limit)
    dr_grace_period=500000;       // don't check dynamic range until this many samples
    ratio_checkpoint=50000;       // don't check for unsolveable exception until this many samples
  } else if (nle_config->nle_mode == 2) {
    // 2-term mode without 1-smr
    precision_target=1.0E-15;     // solve NLE to this level of precision
    test_ratio=25.0;              // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=25;        // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;         // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    defaultrange_multiplier=1.0;  // lowest practical range multiplier, fastest for most formulas
    stalledrange_multiplier=2.0;  // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    slowcheckpoint=1000000;       // progress point to check on slow processes
    stuckprecision=1.0E+30;       // if precision is not past this level by slowcheckpoint, try resetting
  } else {
    // 3-term mode
    precision_target=1.0E-15;     // solve NLE to this level of precision
    test_ratio=25.0;              // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=25;        // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;         // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    defaultrange_multiplier=5.0;  // lowest practical range multiplier, fastest for most formulas
    stalledrange_multiplier=17.0; // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    slowcheckpoint=1000000;       // progress point to check on slow processes
    stuckprecision=1.0E-2;        // if precision is not past this level by slowcheckpoint, try resetting
  }

  // starting v4.0 denominator is always v for this step, other mass ratio factors are now swapped out in cscanner()
  if (nle_state->term1.smrfactor_1minus == 1) {
    term1_mass_sm1=powl((1.0 - ((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm1 / (long double)nle_state->input_sample_v)), (1.0 / (long double)nle_state->term1.exp_inv));
    term1_mass_sm2=powl((1.0 - ((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm2 / (long double)nle_state->input_sample_v)), (1.0 / (long double)nle_state->term1.exp_inv));
    // check if (1-smr) is negative for heaviest mass state
    if ((1.0 - ((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm3)) < 0) {
      term1_mass_sm3=-powl(-(1.0 - ((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm3 / (long double)nle_state->input_sample_v)), (1.0 / (long double)nle_state->term1.exp_inv));
    } else {
      term1_mass_sm3=powl((1.0 - ((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm3 / (long double)nle_state->input_sample_v)), (1.0 / (long double)nle_state->term1.exp_inv));
    }
    // use the same smrfactor for both terms in 1-smr mode
    term2_mass_sm1=powl(((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm1 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term2.exp_inv));
    term2_mass_sm2=powl(((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm2 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term2.exp_inv));
    term2_mass_sm3=powl(((long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm3 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term2.exp_inv));
  } else {
    term1_mass_sm1=powl(((long double)nle_state->input_sample_sm1 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term1.exp_inv));
    term1_mass_sm2=powl(((long double)nle_state->input_sample_sm2 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term1.exp_inv));
    term1_mass_sm3=powl(((long double)nle_state->input_sample_sm3 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term1.exp_inv));
    term2_mass_sm1=powl(((long double)nle_state->input_sample_sm1 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term2.exp_inv));
    term2_mass_sm2=powl(((long double)nle_state->input_sample_sm2 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term2.exp_inv));
    term2_mass_sm3=powl(((long double)nle_state->input_sample_sm3 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term2.exp_inv));
  }
  if (nle_config->nle_mode > 2) {
    term3_mass_sm1= powl(((long double)nle_state->input_sample_sm1 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term3.exp_inv));
    term3_mass_sm2= powl(((long double)nle_state->input_sample_sm2 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term3.exp_inv));
    term3_mass_sm3= powl(((long double)nle_state->input_sample_sm3 / (long double)nle_state->input_sample_v), (1.0 / (long double)nle_state->term3.exp_inv));
  }

  if (nle_config->status_enable == 1) {
    printf("status, Solving phase 1 formula for coefficients, exponents: %s\n", nle_state->exponents_str);
    fflush(stdout);
  }

  //  solve formula for coefficients
  best_precision_last=1.0E99;
  dr_exception=0;
  unsolveable_exception=0;
  while ((best_precision_last > precision_target) && (dr_exception == 0) && (unsolveable_exception == 0)) {
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
    for (samples=0; ((best_precision_last > precision_target) && (dr_exception == 0) && (unsolveable_exception == 0)); samples++) {
      active_ordering_count=0;
      for (ordering=0; ((ordering <= 5) && (dr_exception == 0) && (unsolveable_exception == 0)); ordering++) {
        if ((best_precision_last > 1.0E-7) || (ordering == best_ordering)) { // skip other ordings if one is far enough along
          active_ordering_count++;
          if ((samples > 1) && ((samples % slowcheckpoint) == 0)) { // check on slow processes
            if ((nle_config->smrfactor_1minus_enable == 1) && (samples == ratio_checkpoint)) {
              if (active_ordering_count > 3) {
                unsolveable_exception=1;
#ifdef DEBUG10
                printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, active_ordering_count: %d, unsolveable exception\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], active_ordering_count);
                fflush(stdout);
#endif
              }
            }
#ifdef DEBUG10
            if ((samples % 10000000) == 0) { // rate limit periodic debug prints
              clock_gettime(CLOCK_REALTIME, &endtime);
              elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);

              printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9Le, best_ordering: %d, precision_last: %.9Le, range_factor: %.9Le, i: %d, sm1_test: %.9Le, sm2_test: %.9Le, sm3_test: %.9Le, c1: %.9Le, c2: %.9Le, c3: %.9Le, c1range: %.9Le, c2range: %.9Le, c3range: %.9Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], range_factor, i, sm1_test, sm2_test, sm3_test, c1_center[ordering], c2_center[ordering], c3_center[ordering], c1_range[ordering], c2_range[ordering], c3_range[ordering]);
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
              r=(long double)drand48();
              c1[ordering]=((c1_center[ordering] - c1_range[ordering]) + (r * 2.0 * c1_range[ordering]));
            }
            while (c2[ordering] <= 0.0) {
              r=(long double)drand48();
              c2[ordering]=((c2_center[ordering] - c2_range[ordering]) + (r * 2.0 * c2_range[ordering]));
            }                 
            while (c3[ordering] <= 0.0) {
              r=(long double)drand48();
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
            if (nle_config->nle_mode == 2) {
              // for 2-term mixed mode, we will use term3 as a pseudo term for the mixing of terms 1 and 2
              // This enables us to always solve the NLE regardless of whether the mass spectrum is allowed by the exponents
              // Once solved we will test whether c3 has the correct relationship to c1 and c2
              sm1_test_term1=c1[ordering] * term1_mass_sm1 * term1_mass_sm1;
              sm1_test_term2=c2[ordering] * term2_mass_sm1 * term2_mass_sm1;
              sm1_test_term3=c3[ordering] * term1_mass_sm1 * term2_mass_sm1;
 
              sm2_test_term1=c1[ordering] * term1_mass_sm2 * term1_mass_sm2;
              sm2_test_term2=c2[ordering] * term2_mass_sm2 * term2_mass_sm2;
              sm2_test_term3=c3[ordering] * term1_mass_sm2 * term2_mass_sm2;

              sm3_test_term1=c1[ordering] * term1_mass_sm3 * term1_mass_sm3;
              sm3_test_term2=c2[ordering] * term2_mass_sm3 * term2_mass_sm3;
              sm3_test_term3=c3[ordering] * term1_mass_sm3 * term2_mass_sm3;
              if (nle_config->smrfactor_1minus_enable == 1) {
                sm1_test=sm1_test_term1 + sm1_test_term2 + sm1_test_term3 - 1.0;
                sm2_test=sm2_test_term1 + sm2_test_term2 + sm2_test_term3 - 1.0;
              } else {
                sm1_test=sm1_test_term1 + sm1_test_term2 - sm1_test_term3 - 1.0;
                sm2_test=sm2_test_term1 + sm2_test_term2 - sm2_test_term3 - 1.0;
              }
            } else if (nle_config->nle_mode == 3) {
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
            }

            if ((progress[ordering] < ratio_grace_period) || (((fabsl(sm1_test) / fabsl(sm2_test)) < test_ratio) && ((fabsl(sm2_test) / fabsl(sm1_test)) < test_ratio))) {
              if ((nle_config->nle_mode == 2) && (nle_config->smrfactor_1minus_enable == 1)) {
                sm3_test=sm3_test_term1 + sm3_test_term2 + sm3_test_term3 - 1.0;
              } else if (nle_config->nle_mode == 2) {
                sm3_test=sm3_test_term1 + sm3_test_term2 - sm3_test_term3 - 1.0;
              } else if (nle_config->nle_mode == 3) {
                sm3_test=sm3_test_term1 - sm3_test_term2 + sm3_test_term3 - 1.0;
              }
              if ((progress[ordering] < ratio_grace_period) || (((fabsl(sm1_test) / fabsl(sm3_test)) < test_ratio) && ((fabsl(sm3_test) / fabsl(sm1_test)) < test_ratio) &&\
                                                                ((fabsl(sm2_test) / fabsl(sm3_test)) < test_ratio) && ((fabsl(sm3_test) / fabsl(sm2_test)) < test_ratio))) {
#ifdef DEBUG12
                printf("debug, etest: %.3Le, utest: %.3Le, ttest: %.3Le, sm1_test_term1: %.3Le, sm1_test_term2: %.3Le, sm1_test_term3: %.3Le, sm2_test_term1: %.3Le, sm2_test_term2: %.3Le, sm2_test_term3: %.3Le, sm3_test_term1: %.3Le, sm3_test_term2: %.3Le, sm3_test_term3: %.3Le\n", sm1_test, sm2_test, sm3_test, sm1_test_term1, sm1_test_term2, sm1_test_term3, sm2_test_term1, sm2_test_term2, sm2_test_term3, sm3_test_term1, sm3_test_term2, sm3_test_term3);
                fflush(stdout);
#endif

                precision=(fabsl(sm1_test) + fabsl(sm2_test) + fabsl(sm3_test));
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
                  if (fabsl(sm1_test) > fabsl(sm2_test)) {
                    worst_test=fabsl(sm1_test);
                  } else {
                    worst_test=fabsl(sm2_test);
                  }
                  if (fabsl(sm3_test) > worst_test) {
                    worst_test=fabsl(sm3_test);
                  }
                  range_factor=worst_test * range_multiplier[ordering];
                  if (range_factor > 1.0) {
                    c1_range[ordering]=c1[ordering] / 2.0;
                    c2_range[ordering]=c2[ordering] / 2.0;
                    c3_range[ordering]=c3[ordering] / 2.0;
                  } else {
                    if (nle_config->smrfactor_1minus_enable == 1) {
                      // this is super fast and reliable in (1-smr) mode for some reason
                      c1_range[ordering]=range_factor;
                      c2_range[ordering]=range_factor;
                      c3_range[ordering]=range_factor;
                    } else {
                      c1_range_new=c1[ordering] * range_factor;
                      c1_range[ordering]=((c1_range[ordering] + c1_range_new + c1_range_new) / 3.0);
                      c2_range_new=c2[ordering] * range_factor;
                      c2_range[ordering]=((c2_range[ordering] + c2_range_new + c2_range_new) / 3.0);
                      c3_range_new=c3[ordering] * range_factor;
                      c3_range[ordering]=((c3_range[ordering] + c3_range_new + c3_range_new) / 3.0);
                    }
                  }
                  c1_center[ordering]=c1[ordering];
                  c2_center[ordering]=c2[ordering];
                  c3_center[ordering]=c3[ordering];

                  if (nle_config->smrfactor_1minus_enable == 1) {
                    // calculate dynamic range of coefficients and terms for each particle to see if it exceeds dr_exception_limit
                    dr_high=1.0E-99;
                    dr_low=1.0E+99;
                    if (fabsl(c1[ordering]) > dr_high) {
                      dr_high=fabsl(c1[ordering]);
                    }
                    if (fabsl(c1[ordering]) < dr_low) {
                      dr_low=fabsl(c1[ordering]);
                    }
                    if (fabsl(c2[ordering]) > dr_high) {
                      dr_high=fabsl(c2[ordering]);
                    }
                    if (fabsl(c2[ordering]) < dr_low) {
                      dr_low=fabsl(c2[ordering]);
                    }
                    if (fabsl(c3[ordering]) > dr_high) {
                      dr_high=fabsl(c3[ordering]);
                    }
                    if (fabsl(c3[ordering]) < dr_low) {
                      dr_low=fabsl(c3[ordering]);
                    }
                    dynamicrange_c[ordering]=(dr_high / dr_low);
                    if ((samples > dr_grace_period) && (ordering == best_ordering) && (dynamicrange_c[ordering] > dr_exception_limit)) {
                      dr_exception=1;
#ifdef DEBUG10
                      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, dr_exception (coefficients): %.3Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], dynamicrange_c[ordering]);
                      fflush(stdout);
#endif
                    }
                    dr_high=1.0E-99;
                    dr_low=1.0E+99;
                    if (fabsl(sm1_test_term1) > dr_high) {
                      dr_high=fabsl(sm1_test_term1);
                    }
                    if (fabsl(sm1_test_term1) < dr_low) {
                      dr_low=fabsl(sm1_test_term1);
                    }
                    if (fabsl(sm1_test_term2) > dr_high) {
                      dr_high=fabsl(sm1_test_term2);
                    }
                    if (fabsl(sm1_test_term2) < dr_low) {
                      dr_low=fabsl(sm1_test_term2);
                    }
                    if (fabsl(sm1_test_term3) > dr_high) {
                      dr_high=fabsl(sm1_test_term3);
                    }
                    if (fabsl(sm1_test_term3) < dr_low) {
                      dr_low=fabsl(sm1_test_term3);
                    }
                    dynamicrange_sm1[ordering]=(dr_high / dr_low);
                    if ((samples > dr_grace_period) && (ordering == best_ordering) && (dynamicrange_sm1[ordering] > dr_exception_limit)) {
                      dr_exception=1;
#ifdef DEBUG10
                      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, dr_exception (sm1_test): %.3Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], dynamicrange_sm1[ordering]);
                      fflush(stdout);
#endif
                    }
                    dr_high=1.0E-99;
                    dr_low=1.0E+99;
                    if (fabsl(sm2_test_term1) > dr_high) {
                      dr_high=fabsl(sm2_test_term1);
                    }
                    if (fabsl(sm2_test_term1) < dr_low) {
                      dr_low=fabsl(sm2_test_term1);
                    }
                    if (fabsl(sm2_test_term2) > dr_high) {
                      dr_high=fabsl(sm2_test_term2);
                    }
                    if (fabsl(sm2_test_term2) < dr_low) {
                      dr_low=fabsl(sm2_test_term2);
                    }
                    if (fabsl(sm2_test_term3) > dr_high) {
                      dr_high=fabsl(sm2_test_term3);
                    }
                    if (fabsl(sm2_test_term3) < dr_low) {
                      dr_low=fabsl(sm2_test_term3);
                    }
                    dynamicrange_sm2[ordering]=(dr_high / dr_low);
                    if ((samples > dr_grace_period) && (ordering == best_ordering) && (dynamicrange_sm2[ordering] > dr_exception_limit)) {
                      dr_exception=1;
#ifdef DEBUG10
                      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, dr_exception (sm2_test): %.3Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], dynamicrange_sm2[ordering]);
                      fflush(stdout);
#endif
                    }
                    dr_high=1.0E-99;
                    dr_low=1.0E+99;
                    if (fabsl(sm3_test_term1) > dr_high) {
                      dr_high=fabsl(sm3_test_term1);
                    }
                    if (fabsl(sm3_test_term1) < dr_low) {
                      dr_low=fabsl(sm3_test_term1);
                    }
                    if (fabsl(sm3_test_term2) > dr_high) {
                      dr_high=fabsl(sm3_test_term2);
                    }
                    if (fabsl(sm3_test_term2) < dr_low) {
                      dr_low=fabsl(sm3_test_term2);
                    }
                    if (fabsl(sm3_test_term3) > dr_high) {
                      dr_high=fabsl(sm3_test_term3);
                    }
                    if (fabsl(sm3_test_term3) < dr_low) {
                      dr_low=fabsl(sm3_test_term3);
                    }
                    dynamicrange_sm3[ordering]=(dr_high / dr_low);
                    if ((samples > dr_grace_period) && (ordering == best_ordering) && (dynamicrange_sm3[ordering] > dr_exception_limit)) {
                      dr_exception=1;
#ifdef DEBUG10
                      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, dr_exception (sm3_test): %.3Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], dynamicrange_sm3[ordering]);
                      fflush(stdout);
#endif
                    }
                  } // end if smrfactor_1minus_enable == 1

#ifdef DEBUG11
                  clock_gettime(CLOCK_REALTIME, &endtime);
                  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                  printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9Le, best_ordering: %d, precision_last: %.9Le, range_factor: %.9Le, i: %d, sm1_test: %.9Le, sm2_test: %.9Le, sm3_test: %.9Le, c1: %.9Le, c2: %.9Le, c3: %.9Le, c1range: %.9Le, c2range: %.9Le, c3range: %.9Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], range_factor, i, sm1_test, sm2_test, sm3_test, c1[ordering], c2[ordering], c3[ordering], c1_range[ordering], c2_range[ordering], c3_range[ordering]);
                  fflush(stdout);
#endif
                } // end if precision
              } // end sm3_test ratios
            } // end sm2_test ratios
          } // end if good_coefficients
        } // end if best
      } // end for ordering
    }  // end for samples
  }  // end while precision_last

  if ((best_precision_last < precision_target) && (isnan(best_precision_last) == 0)) {
    two_term_test=c3_center[best_ordering] / (sqrtl(c1_center[best_ordering] * c2_center[best_ordering]));
    if (nle_config->status_enable == 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      if (nle_config->nle_mode == 2) {
        if (nle_config->smrfactor_1minus_enable == 1) {
          printf("status, Solved  phase 1 formula for coefficients, input sample: %i, exponents:  %s, sm3: %.14e, term1.smrfactor: %.14e, two_term_test: %.14Le, sqrt(c1): %.14Le, sqrt(c2): %.14Le, samples: %lld, ordering: %d, precision: %.3Le (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, nle_state->term1.smrfactor, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]), samples, best_ordering, best_precision_last, elapsed_time);
          // verbose for testing
          //printf("status, Solved  phase 1 formula for coefficients, input sample: %i, exponents:  %s, sm3: %.14e, term1.smrfactor: %.14e, alpha_w: %.14e, two_term_test: %.14Le, sqrt(c1): %.14Le, sqrt(c2): %.14Le, sm1_test_1: %.9Le, sm1_test_2: %.9Le, sm1_test_3: %.9Le, sm2_test_1: %.9Le, sm2_test_2: %.9Le, sm2_test_3: %.9Le, sm3_test_1: %.9Le, sm3_test_2: %.9Le, sm3_test_3: %.9Le, dr_c: %.9Le, dr_e: %.9Le, dr_u: %.9Le, dr_t: %.9Le, samples: %lld, ordering: %d, precision: %.3Le (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, nle_state->term1.smrfactor, nle_state->input_sample_alpha_w, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]), sm1_test_term1, sm1_test_term2, sm1_test_term3, sm2_test_term1, sm2_test_term2, sm2_test_term3, sm3_test_term1, sm3_test_term2, sm3_test_term3, dynamicrange_c[best_ordering], dynamicrange_sm1[best_ordering], dynamicrange_sm2[best_ordering], dynamicrange_sm3[best_ordering], samples, best_ordering, best_precision_last, elapsed_time);
        } else {
          printf("status, Solved  phase 1 formula for coefficients, input sample: %i, exponents:  %s, sm3: %.9e, samples: %lld, ordering: %d, two_term_test: %.9Le, precision: %.3Le (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, samples, best_ordering, two_term_test, precision, elapsed_time);
        }
      } else if (nle_config->nle_mode == 3) {
        printf("status, Solved  phase 1 formula for coefficients, input sample: %i, exponents:  %s, sm3: %.9e, samples: %lld, ordering: %d, precision: %.3Le (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, samples, best_ordering, precision, elapsed_time);
      }
      printf("status, +------------+--------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
    printf("status, | Exponents  |  Mass ratio  | NLE term |   Coefficient   |  C * term(sm1)  |  C * term(sm2)  |  C * term(sm3)  |\n");
      printf("status, +------------+--------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
      if (nle_config->nle_mode == 2) {
        if (nle_state->term1.smrfactor_1minus == 1) {
          printf("status, | %10s | 1-(smrf M/v) |  t1^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
          printf("status, | %10s |    smrf M/v  |  t2^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
        } else {
          printf("status, | %10s |     M/v      |  t1^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
          printf("status, | %10s |     M/v      |  t2^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
        }
        printf("status, | %10s |              | t1 * t2  | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c3_center[best_ordering], sm1_test_term3, sm2_test_term3, sm3_test_term3);
      } else if (nle_config->nle_mode == 3) {
        printf("status, | %10s |     M/v      |  term1   | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
        printf("status, | %10s |     M/v      |  term2   | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
        printf("status, | %10s |     M/v      |  term3   | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c3_center[best_ordering], sm1_test_term3, sm2_test_term3, sm3_test_term3);
      }
      printf("status, +------------+--------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
      fflush(stdout);

      // for debugging very long phase 1 solutions
      //exit(0);
    } // end if status_enable

    // only run cscanner if mode=3 or two_term_test is an interesting integer match
    // for this interesting() check  we use a fixed filter of 3.  In cscanner two_term_test is evaluated with phase1_filter on c3 which may be more restrictive
    // This will help identify interesting geometries for further inspection
    if ((nle_config->nle_mode != 2) || ((two_term_test >= 0.98) && interesting(3, nle_config->phase1_int_match_max, nle_config->phase1_int_match_filter, two_term_test))) {
      matches_count_start=nle_state->phase1_matches_count;
      clock_gettime(CLOCK_REALTIME, &starttime);
      if (nle_config->nle_mode == 2) {
        // in two term mode take the square root of c1 and c2 and use two_term_test as c3
        nle_state->term1.coefficient=(double)sqrtl(c1_center[best_ordering]);
        nle_state->term2.coefficient=(double)sqrtl(c2_center[best_ordering]);
        nle_state->term3.coefficient=(double)two_term_test;
        if (nle_config->status_enable == 1) {
          if (nle_config->smrfactor_1minus_enable == 1) {
            printf("status, Found interesting two_term_test, input_sample: %i, exponents: %s, sm3: %.14e, term1_smrfactor: %.14e, two_term_test: %.14Le, sqrt(c1): %.14Le, sqrt(c2): %.14Le\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, nle_state->term1.smrfactor, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]));
            fflush(stdout);
          } else {
            printf("status, Found interesting two_term_test, input_sample: %i, exponents: %s, sm3: %.14e, two_term_test: %.14Le, sqrt(c1): %.14Le, sqrt(c2): %.14Le\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]));
            fflush(stdout);
          }
        }
      } else {
        nle_state->term1.coefficient=(double)c1_center[best_ordering];
        nle_state->term2.coefficient=(double)c2_center[best_ordering];
        nle_state->term3.coefficient=(double)c3_center[best_ordering];
      }

      cscanner(nle_config, nle_state);

      clock_gettime(CLOCK_REALTIME, &endtime);
      matches_count_end=nle_state->phase1_matches_count;
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      if (nle_config->status_enable ==1 ) {
        printf("status, Found %d interesting coefficient multipliers (%6.4fs)\n", (matches_count_end-matches_count_start), elapsed_time);
        fflush(stdout);
      }
    } else {
      printf("status, two_term_test was not close enough to an interesting integer, skipping factoring process.\n");
      fflush(stdout);
    }
  } else { // faled to solve, should only happen in 1-smr mode
    if (nle_config->status_enable == 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      two_term_test=c3_center[best_ordering] / (sqrtl(c1_center[best_ordering] * c2_center[best_ordering]));
      printf("status, Failed to solve  phase 1 formula for coefficients, input sample: %i, exponents:  %s, sm3: %.14e, term1.smrfactor: %.14e (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, nle_state->term1.smrfactor, elapsed_time);
    }
  } // end if best_precision_last
}
