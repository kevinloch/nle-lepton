#include <stdio.h>
#include <stdlib.h> // drand
#include <math.h> // powl
#include <time.h>
#include "nle-lepton.h"
#include "util.h"
#include "cscanner.h"
#include "getFormulaStr.h"

//#define DEBUG10
//#define DEBUG11
//#define DEBUG12

int solveNLEforCoefficients(nle_config_t *nle_config, nle_state_t *nle_state) {
  // solve a three term polynomial-like non-linear equation for the unknown coefficients given known roots (particle masses)
  // returns 0 if equation was solved, 1 if failed
  int i;
  long long samples=0;
  struct timespec starttime;
  struct timespec endtime;
  double elapsed_time;
  int matches_count_start;
  int matches_count_end;
  int good_coefficients;
  int c2_gt_c1, c2_gt_c3, c1_gt_c3;
  char smrf_str[80];
  char smrfactor_mass_str[32];
  char out_str_01[512];
  char signature_str[4];
  char exec_str[512];
  long double term1_exp;
  long double term2_exp;
  long double term3_exp;

  clock_gettime(CLOCK_REALTIME, &starttime);

  //  mc test vars
  long double r;
  long double smrfactor_mass=0;
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
  long double range_factor=0;
  long double smrf_sm1=0;
  long double term1_mass_sm1=0;
  long double term2_mass_sm1=0;
  long double term3_mass_sm1=0;
  long double smrf_sm2=0;
  long double term1_mass_sm2=0;
  long double term2_mass_sm2=0;
  long double term3_mass_sm2=0;
  long double smrf_sm3=0;
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
  int unsolvable_exception;
  char unsolvable_exception_str[128];

  // mc outputs
  int progress[6];
  int last_progress[6];
  int stalled[6];
  int ordering;
  int best_ordering=-1;
  int ordering_enabled[6];
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
  long double best_ordering_threshold;
  long double test_ratio;
  int ratio_grace_period;
  long long stalled_limit;
  long double defaultrange_multiplier;
  long double stalledrange_multiplier;
  int slowcheckpoint;
  long double stuckprecision;
  long double dr_exception_limit=0.0;
  long double dr_grace_period=0;
  long double two_term_precision;

  // these tunings affect speed and reliability, adjust with extreme care
  if (nle_config->smrfactor_1minus_enable == 1) {
    // 2-term mode with 1-smr
    precision_target=1.0E-15;        // solve NLE to this level of precision
    best_ordering_threshold=1.0E-99; // Only process best ordering after it reaches this precision.  This speeds up phase1 processing.  Should be disabled for (1-smr) mode
    test_ratio=25.0;                 // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=3;            // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;            // most formulas can make further progress with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    defaultrange_multiplier=2.0;     // lowest practical range multiplier, fastest for most formulas
    stalledrange_multiplier=3.0;     // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    slowcheckpoint=10000;            // number of samples to check on slow processes for reporting and/or reset
    stuckprecision=2.0E+99;          // if precision is not past this level by slowcheckpoint, try resetting
    dr_exception_limit=1.0E+16;      // Mamimum allowed dynamic range of any term.  This is used in 1-smr mode only
    dr_grace_period=25;              // don't check dynamic range until this much progress
    two_term_precision=1.0E-3;       // check two_term_test after reaching this precision.  Disabling current ordering if not close enough to an integer
  } else if (nle_config->nle_mode == 2) {
    // 2-term mode without 1-smr
    precision_target=1.0E-15;        // solve NLE to this level of precision
    best_ordering_threshold=1.0E-7;  // Only process best ordering after it reaches this precision.  This speeds up phase1 processing
    test_ratio=25.0;                 // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=25;           // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;            // most formulas can make further progress with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    defaultrange_multiplier=1.0;     // lowest practical range multiplier, fastest for most formulas
    stalledrange_multiplier=2.0;     // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    slowcheckpoint=10000;            // number of samples to check on slow processes for reporting and/or reset
    stuckprecision=2.0E+99;          // if precision is not past this level by slowcheckpoint, try resetting
  } else {
    // 3-term mode
    precision_target=1.0E-15;        // solve NLE to this level of precision
    best_ordering_threshold=1.0E-7;  // Only process best ordering after it reaches this precision.  This speeds up phase1 processing
    test_ratio=25.0;                 // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=25;           // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;            // most formulas can make further progress with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    defaultrange_multiplier=5.0;     // lowest practical range multiplier, fastest for most formulas
    stalledrange_multiplier=17.0;    // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    slowcheckpoint=10000;            // number of samples to check on slow processes for reporting and/or reset
    stuckprecision=1.0E-1;           // if precision is not past this level by slowcheckpoint, try resetting
  }

  term1_exp = 1.0 / (long double)nle_state->term1.exp_inv;
  term2_exp = 1.0 / (long double)nle_state->term2.exp_inv;
  term3_exp = 1.0 / (long double)nle_state->term3.exp_inv;

  // load smrfactor string, and set solution mass ratio reference mass if (1-smr) is enabled
  if (nle_config->smrfactor_1minus_enable == 1) {
    getSmrfStr(nle_config, smrf_str, nle_state->term1.current_smrfactors, nle_state->term1.smrfactor);
    if (nle_state->term1.smrfactor_mass_id == 0) {
      smrfactor_mass=(long double)nle_state->input_sample_mp;
      sprintf(smrfactor_mass_str, "mP   ");
    } else if (nle_state->term1.smrfactor_mass_id == 1) {
      smrfactor_mass=(long double)nle_state->input_sample_v;
      sprintf(smrfactor_mass_str, "v    ");
    } else if (nle_state->term1.smrfactor_mass_id == 2) {
      smrfactor_mass=(long double)nle_state->input_sample_mz;
      sprintf(smrfactor_mass_str, "mz   ");
    } else if (nle_state->term1.smrfactor_mass_id == 3) {
      smrfactor_mass=(long double)nle_state->input_sample_mw;
      sprintf(smrfactor_mass_str, "mw   ");
    } else if (nle_state->term1.smrfactor_mass_id == 4) {
      smrfactor_mass=(long double)nle_state->input_sample_mh0;
      sprintf(smrfactor_mass_str, "mh0  ");
    } else if (nle_state->term1.smrfactor_mass_id == 5) {
      smrfactor_mass=(long double)nle_state->input_sample_muser;
      sprintf(smrfactor_mass_str, "muser");
    }
    if (nle_config->phase1_status_enable == 1) { // not needed in (1-smr) mode as p1 should quickly solve or abort
      printf("status, Solving          phase 1 formula for coefficients, input smaple: %lld, exponents:  %s, mixing polarity: %s, mass config: %s, sm3: %.14e, smrfactor mass: %.14Le, smrf: %s\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->nle_mixing_polarity_str, nle_state->smrfactor_mass_configuration_str, nle_state->input_sample_sm3, smrfactor_mass, smrf_str);
      fflush(stdout);
    }
  } else {
    smrf_str[0]=0;
    smrfactor_mass=(long double)nle_state->input_sample_v;
    if (nle_config->phase1_status_enable == 1) { // not needed in (1-smr) mode as p1 should quickly solve or abort
      printf("status, Solving phase 1 formula for coefficients, input sample: %lld, exponents:  %s\n", nle_state->phase1_seq, nle_state->exponents_str);
      fflush(stdout);
    }
  }

  // set mass component of each term
  unsolvable_exception=0;
  sprintf(unsolvable_exception_str, "(null)");
  if (nle_config->smrfactor_1minus_enable == 1) {
    if (nle_state->smrfactor_mass_configuration == 1) {
      smrf_sm1=(long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm1 / smrfactor_mass;
      smrf_sm2=(long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm2 / smrfactor_mass;
      smrf_sm3=(long double)nle_state->term1.smrfactor * (long double)nle_state->input_sample_sm3 / smrfactor_mass;
    } else if (nle_state->smrfactor_mass_configuration == 0) {
      smrf_sm1=(long double)nle_state->term1.smrfactor * smrfactor_mass / (long double)nle_state->input_sample_sm1;
      smrf_sm2=(long double)nle_state->term1.smrfactor * smrfactor_mass / (long double)nle_state->input_sample_sm2;
      smrf_sm3=(long double)nle_state->term1.smrfactor * smrfactor_mass / (long double)nle_state->input_sample_sm3;
    }
    // check if (1-smr) is negative for sm1 and invert inside and outside radical
    if ((1.0 - smrf_sm1) < 0) {
      term1_mass_sm1=-powl(-(1.0 - smrf_sm1), term1_exp);
    } else {
      term1_mass_sm1=powl((1.0 - smrf_sm1), term1_exp);
    }
    // check if (1-smr) is negative for sm2 and invert inside and outside radical
    if ((1.0 - smrf_sm2) < 0) {
      term1_mass_sm2=-powl(-(1.0 - smrf_sm2), term1_exp);
    } else {
      term1_mass_sm2=powl((1.0 - smrf_sm2), term1_exp);
    }
    // check if (1-smr) is negative for sm3 and invert inside and outside radical
    if ((1.0 - smrf_sm3) < 0) {
      term1_mass_sm3=-powl(-(1.0 - smrf_sm3), term1_exp);
    } else {
      term1_mass_sm3=powl((1.0 - smrf_sm3), term1_exp);
    }
    term2_mass_sm1=powl(smrf_sm1, term2_exp);
    term2_mass_sm2=powl(smrf_sm2, term2_exp);
    term2_mass_sm3=powl(smrf_sm3, term2_exp);
    // check if term 1 is even exponent and any (1-smr) is negative
    if (((nle_state->term1.exp_inv % 2) == 0) && (((1.0 - smrf_sm1) <= 0) || ((1.0 - smrf_sm2) <= 0) || ((1.0 - smrf_sm3) <= 0))) {
      unsolvable_exception=1;
      sprintf(unsolvable_exception_str, "term 1 has negative (1-smr) with even exponent                                  ");
#ifdef DEBUG10
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, %s\n", nle_state->exponents_str, samples, elapsed_time, unsolvable_exception_str);
      fflush(stdout);
#endif
    }
    // check if polarity makes formula unsolvable
    if ((nle_state->nle_mixing_polarity == 0) && ((1.0 - smrf_sm1) <= 0) && ((1.0 - smrf_sm2) <= 0) && ((1.0 - smrf_sm3) <= 0)) {
      unsolvable_exception=1;
      sprintf(unsolvable_exception_str, "mixed term polarity is unsolvable                                               ");
#ifdef DEBUG10  
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, %s\n", nle_state->exponents_str, samples, elapsed_time, unsolvable_exception_str);
      fflush(stdout);
#endif
    }
    if ((nle_state->nle_mixing_polarity == 1) && ((1.0 - smrf_sm1) >= 0) && ((1.0 - smrf_sm2) >= 0) && ((1.0 - smrf_sm3) >= 0)) {
      unsolvable_exception=1;
      sprintf(unsolvable_exception_str, "mixed term polarity is unsolvable                                               ");
#ifdef DEBUG10
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, %s\n", nle_state->exponents_str, samples, elapsed_time, unsolvable_exception_str);
      fflush(stdout);
#endif
    }
  } else {
    term1_mass_sm1=powl(((long double)nle_state->input_sample_sm1 / smrfactor_mass), term1_exp);
    term1_mass_sm2=powl(((long double)nle_state->input_sample_sm2 / smrfactor_mass), term1_exp);
    term1_mass_sm3=powl(((long double)nle_state->input_sample_sm3 / smrfactor_mass), term1_exp);
    term2_mass_sm1=powl(((long double)nle_state->input_sample_sm1 / smrfactor_mass), term2_exp);
    term2_mass_sm2=powl(((long double)nle_state->input_sample_sm2 / smrfactor_mass), term2_exp);
    term2_mass_sm3=powl(((long double)nle_state->input_sample_sm3 / smrfactor_mass), term2_exp);
  }
  if (nle_config->nle_mode == 3) {
    term3_mass_sm1=powl(((long double)nle_state->input_sample_sm1 / smrfactor_mass), term3_exp);
    term3_mass_sm2=powl(((long double)nle_state->input_sample_sm2 / smrfactor_mass), term3_exp);
    term3_mass_sm3=powl(((long double)nle_state->input_sample_sm3 / smrfactor_mass), term3_exp);
  }

  //  solve formula for coefficients
  best_precision_last=1.0E99;
  while ((best_precision_last > precision_target) && (unsolvable_exception == 0)) {
    //  init outputs
    for (ordering=0; ordering <= 5; ordering++) {
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
      last_progress[ordering]=0;
      stalled[ordering]=0;
      range_multiplier[ordering]=defaultrange_multiplier;
      ordering_enabled[ordering]=1;
    }
    best_precision_last=1.0E99;
    best_ordering=-1;
    for (samples=0; ((best_precision_last > precision_target) && (unsolvable_exception == 0)); samples++) {

      // samples limit
      if ((nle_config->smrfactor_1minus_enable == 1) && (samples > nle_config->phase1_mc_samples_limit)) {
        unsolvable_exception=1;
        sprintf(unsolvable_exception_str, "exceeded samples limit: %12lld                                             ", nle_config->phase1_mc_samples_limit);
#ifdef DEBUG10  
        clock_gettime(CLOCK_REALTIME, &endtime);
        elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
        printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, %s\n", nle_state->exponents_str, samples, elapsed_time, unsolvable_exception_str);
        fflush(stdout);
#endif
      }

      for (ordering=0; ((ordering <= 5) && (unsolvable_exception == 0)); ordering++) {
        if (ordering_enabled[ordering] == 1) {

          // in (1-smr) mode, periodically check if any progress has been made since last unsolvable_checkpoint
          if ((nle_config->smrfactor_1minus_enable == 1) && (samples > 1) && ((samples % nle_config->phase1_unsolvable_checkpoint) == 0)) {
            // check progress
            if (ordering == best_ordering) {
              // best ordering, check if any progress has been made on best_ordering since last unsolvable_checkpoint
              if ((progress[ordering] - last_progress[ordering]) > 0) {
                last_progress[ordering]=progress[ordering];
              } else {
                // no progress, abort processing
                unsolvable_exception=1;
                sprintf(unsolvable_exception_str, "best ordering did not make any progress since last unsolvable_checkpoint        ");
#ifdef DEBUG10
                clock_gettime(CLOCK_REALTIME, &endtime);
                elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, %s\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], unsolvable_exception_str);
                fflush(stdout);
#endif
              }
            } else {
              // not best ordering, just copy progress to last_progress
              last_progress[ordering]=progress[ordering];
            }
          } // end unsolvable_checkpoint

          // periodically check on slow processes for reporting, coefficient dynamic range, two_term_test, and if stuck
          if ((samples > 1) && ((samples % slowcheckpoint) == 0)) {
#ifdef DEBUG10
            if ((samples % 10000000) == 0) { // rate limit periodic debug prints
              clock_gettime(CLOCK_REALTIME, &endtime);
              elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
              two_term_test=c3_center[ordering] / (sqrtl(c1_center[ordering] * c2_center[ordering]));
              printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9Le, best_ordering: %d, precision_last: %.9Le, range_factor: %.9Le, two-term test: %.9Le, sm1_test: %.9Le, sm2_test: %.9Le, sm3_test: %.9Le, c1: %.9Le, c2: %.9Le, c3: %.9Le, c1range: %.9Le, c2range: %.9Le, c3range: %.9Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], range_factor, two_term_test, sm1_test, sm2_test, sm3_test, c1_center[ordering], c2_center[ordering], c3_center[ordering], c1_range[ordering], c2_range[ordering], c3_range[ordering]);
              fflush(stdout);
            }
#endif

            if (best_precision_last < best_ordering_threshold) {
              // best_ordering has sufficient precision, disable all other orderings
              for (i=0; i <= 5; i++) {
                if ((i != best_ordering) && (ordering_enabled[i] == 1)) {
                  ordering_enabled[i]=0;
#ifdef DEBUG10
                  printf("debug, best_ordering: %d, best_precision_last: %.9Le, disabling ordering: %d\n", best_ordering, best_precision_last, i);
                  fflush(stdout);
#endif
                }
              }
            } // end of best_ordering check

            // In (1-smr) mode, if precision is high enough, abort processing if two_term_test is out of bounds (not even close to interesting)
            if ((nle_config->smrfactor_1minus_enable == 1) && (precision_last[ordering] < two_term_precision) && (ordering == best_ordering)) {
              two_term_test=c3_center[ordering] / (sqrtl(c1_center[ordering] * c2_center[ordering]));
              if ((two_term_test < (long double)nle_config->phase1_two_term_test_min) || (two_term_test > (long double)nle_config->phase1_two_term_test_max)) {
                unsolvable_exception=1;
                sprintf(unsolvable_exception_str, "two-term test is not within a reasonable range: %.14Le            ", two_term_test);
#ifdef DEBUG10
              clock_gettime(CLOCK_REALTIME, &endtime);
              elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
              printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, %s\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], unsolvable_exception_str);
              fflush(stdout);
#endif
              }
            } // end of two_term_test check

            // in (1-smr) mode, calculate dynamic range of coefficients and terms for each particle and abort processing it exceeds dr_exception_limit
            // we test coefficients here at the frequent slowcheckpoint since we have them stored between samples. Individual term components are tested after progress is made since those are not stored between samples
            if ((nle_config->smrfactor_1minus_enable == 1) && (ordering == best_ordering) && (progress[ordering] > dr_grace_period)) {
              dr_high=1.0E-99;
              dr_low=1.0E+99;
              if (fabsl(c1_center[ordering]) > dr_high) {
                dr_high=fabsl(c1_center[ordering]);
              }
              if (fabsl(c1_center[ordering]) < dr_low) {
                dr_low=fabsl(c1_center[ordering]);
              }
              if (fabsl(c2_center[ordering]) > dr_high) {
                dr_high=fabsl(c2_center[ordering]);
              }
              if (fabsl(c2_center[ordering]) < dr_low) {
                dr_low=fabsl(c2_center[ordering]);
              }
              if (fabsl(c3_center[ordering]) > dr_high) {
                dr_high=fabsl(c3_center[ordering]);
              }
              if (fabsl(c3_center[ordering]) < dr_low) {
                dr_low=fabsl(c3_center[ordering]);
              }
              dynamicrange_c[ordering]=(dr_high / dr_low);
              if (dynamicrange_c[ordering] > dr_exception_limit) {
                unsolvable_exception=1;
                sprintf(unsolvable_exception_str, "coefficient dynamic range is too high: %.3Le                                 ", dynamicrange_c[ordering]);
#ifdef DEBUG10  
                clock_gettime(CLOCK_REALTIME, &endtime);
                elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, %s\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], unsolvable_exception_str);
                fflush(stdout);
#endif
              }
            } // end coefficient dr check

            // check if this ordering is stuck and needs to be reset
            if ((progress[ordering] == ratio_grace_period) || (precision_last[ordering] > stuckprecision)) {
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
            } // end if stuck
          } // end if slowcheckpoint

          // check if stalled (too many samples since last progress increment)
          if (stalled[ordering] == stalled_limit) {
            range_multiplier[ordering]=stalledrange_multiplier; // may be a slow solution, try bigger multiplier
            //stalled[ordering]=0;
#ifdef DEBUG11
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
            printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, stalled\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering]);
            fflush(stdout);
#endif
          } // end stalled check
          stalled[ordering]++;

          // generate random coefficients within limits
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
          i=0;
          good_coefficients=0;
          while ((i < 200) && (good_coefficients == 0)) {
            i++;
            c1[ordering]=0.0;
            c2[ordering]=0.0;
            c3[ordering]=0.0;
            while (c1[ordering] <= 0.0) {
              r=pcg_ldrand64(nle_state);
              c1[ordering]=c1_center[ordering] - c1_range[ordering] + (r * c1_range[ordering] * 2.0);
            }
            while (c2[ordering] <= 0.0) {
              r=pcg_ldrand64(nle_state);
              c2[ordering]=c2_center[ordering] - c2_range[ordering] + (r * c2_range[ordering] * 2.0);
            }                 
            while (c3[ordering] <= 0.0) {
              r=pcg_ldrand64(nle_state);
              c3[ordering]=c3_center[ordering] - c3_range[ordering] + (r * c3_range[ordering] * 2.0);
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
              sm1_test_term1=c1[ordering] * term1_mass_sm1 * term1_mass_sm1;
              sm1_test_term2=c2[ordering] * term2_mass_sm1 * term2_mass_sm1;
              sm1_test_term3=c3[ordering] * term1_mass_sm1 * term2_mass_sm1;
              sm2_test_term1=c1[ordering] * term1_mass_sm2 * term1_mass_sm2;
              sm2_test_term2=c2[ordering] * term2_mass_sm2 * term2_mass_sm2;
              sm2_test_term3=c3[ordering] * term1_mass_sm2 * term2_mass_sm2;
              sm3_test_term1=c1[ordering] * term1_mass_sm3 * term1_mass_sm3;
              sm3_test_term2=c2[ordering] * term2_mass_sm3 * term2_mass_sm3;
              sm3_test_term3=c3[ordering] * term1_mass_sm3 * term2_mass_sm3;
              if (nle_config->smrfactor_1minus_enable == 1) { // two different mixing polarity options for (1-smr) mode
                if (nle_state->nle_mixing_polarity == 0) {
                  sm1_test=sm1_test_term1 + sm1_test_term2 - sm1_test_term3 - 1.0;
                  sm2_test=sm2_test_term1 + sm2_test_term2 - sm2_test_term3 - 1.0;
                  sm3_test=sm3_test_term1 + sm3_test_term2 - sm3_test_term3 - 1.0;
                } else if (nle_state->nle_mixing_polarity == 1) {
                  sm1_test=sm1_test_term1 + sm1_test_term2 + sm1_test_term3 - 1.0;
                  sm2_test=sm2_test_term1 + sm2_test_term2 + sm2_test_term3 - 1.0;
                  sm3_test=sm3_test_term1 + sm3_test_term2 + sm3_test_term3 - 1.0;
                }
              } else { // 2-term without (1-smr)
                sm1_test=sm1_test_term1 + sm1_test_term2 - sm1_test_term3 - 1.0;
                sm2_test=sm2_test_term1 + sm2_test_term2 - sm2_test_term3 - 1.0;
                sm3_test=sm3_test_term1 + sm3_test_term2 - sm3_test_term3 - 1.0;
              }
            } else if (nle_config->nle_mode == 3) {
              // three independent terms
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
              sm3_test=sm3_test_term1 - sm3_test_term2 + sm3_test_term3 - 1.0;
            }

            if ((progress[ordering] < ratio_grace_period) || (((fabsl(sm1_test) / fabsl(sm2_test)) < test_ratio) && ((fabsl(sm2_test) / fabsl(sm1_test)) < test_ratio))) {
              if ((progress[ordering] < ratio_grace_period) || (((fabsl(sm1_test) / fabsl(sm3_test)) < test_ratio) && ((fabsl(sm3_test) / fabsl(sm1_test)) < test_ratio) &&\
                                                              ((fabsl(sm2_test) / fabsl(sm3_test)) < test_ratio) && ((fabsl(sm3_test) / fabsl(sm2_test)) < test_ratio))) {
#ifdef DEBUG12
                printf("debug, sm1_test: %21.14Le, sm2_test: %21.14Le, sm3_test: %21.14Le, sm1_test_term1: %21.14Le, sm1_test_term2: %21.14Le, sm1_test_term3: %21.14Le, sm2_test_term1: %21.14Le, sm2_test_term2: %21.14Le, sm2_test_term3: %21.14Le, sm3_test_term1: %21.14Le, sm3_test_term2: %21.14Le, sm3_test_term3: %21.14Le, c1: %21.14Le, c2: %21.14Le, c3: %21.14Le\n", sm1_test, sm2_test, sm3_test, sm1_test_term1, sm1_test_term2, sm1_test_term3, sm2_test_term1, sm2_test_term2, sm2_test_term3, sm3_test_term1, sm3_test_term2, sm3_test_term3, c1[ordering], c2[ordering], c3[ordering]);
                fflush(stdout);
#endif

                precision=(fabsl(sm1_test) + fabsl(sm2_test) + fabsl(sm3_test));
                if (precision < (precision_last[ordering] * 0.9999)) {
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

                  // check dynamic range of individual term components and abort if out of range
                  if ((nle_config->smrfactor_1minus_enable == 1) && (ordering == best_ordering) && (samples > dr_grace_period)) {
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
                    if (dynamicrange_sm1[ordering] > dr_exception_limit) {
                      unsolvable_exception=1;
                      sprintf(unsolvable_exception_str, "sm1_test dynamic range is too high: %.3Le                                   ", dynamicrange_sm1[ordering]);
#ifdef DEBUG10  
                      clock_gettime(CLOCK_REALTIME, &endtime);
                      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, %s\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], unsolvable_exception_str);
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
                    if (dynamicrange_sm2[ordering] > dr_exception_limit) {
                      unsolvable_exception=1;
                      sprintf(unsolvable_exception_str, "sm2_test dynamic range is too high: %.3Le                                   ", dynamicrange_sm2[ordering]);
#ifdef DEBUG10        
                      clock_gettime(CLOCK_REALTIME, &endtime);
                      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, %s\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], unsolvable_exception_str);
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
                    if (dynamicrange_sm3[ordering] > dr_exception_limit) {
                      unsolvable_exception=1;
                      sprintf(unsolvable_exception_str, "sm3_test dynamic range is too high: %.3Le                                   ", dynamicrange_sm3[ordering]);
#ifdef DEBUG10        
                      clock_gettime(CLOCK_REALTIME, &endtime);
                      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, %s\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], unsolvable_exception_str);
                      fflush(stdout);
#endif
                    }
                  } // end if smrfactor_1minus_enable == 1

#ifdef DEBUG11
                  clock_gettime(CLOCK_REALTIME, &endtime);
                  elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
                  two_term_test=c3_center[ordering] / (sqrtl(c1_center[ordering] * c2_center[ordering]));
                  printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, ordering: %d, progress: %6d, best_precision_last: %.9Le, best_ordering: %d, precision_last: %.9Le, range_factor: %.9Le, two-term test: %.9Le, sm1_test: %.9Le, sm2_test: %.9Le, sm3_test: %.9Le, c1: %.9Le, c2: %.9Le, c3: %.9Le, c1range: %.9Le, c2range: %.9Le, c3range: %.9Le\n", nle_state->exponents_str, samples, elapsed_time, ordering, progress[ordering], best_precision_last, best_ordering, precision_last[ordering], range_factor, two_term_test, sm1_test, sm2_test, sm3_test, c1_center[ordering], c2_center[ordering], c3_center[ordering], c1_range[ordering], c2_range[ordering], c3_range[ordering]);
                  fflush(stdout);
#endif
                } // end if precision
              } // end sm3_test ratios
            } // end sm2_test ratios
          } // end if good_coefficients
        } // end if ordering_enabled
      } // end for ordering
    }  // end for samples
  }  // end while precision_last

  if ((best_precision_last < precision_target) && (isnan(best_precision_last) == 0)) {
    two_term_test=c3_center[best_ordering] / (sqrtl(c1_center[best_ordering] * c2_center[best_ordering]));
    clock_gettime(CLOCK_REALTIME, &endtime);
    elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
    if (nle_config->nle_mode == 2) {
      if (sm1_test_term3 > 0) {
        signature_str[0]='+';
      } else {
        signature_str[0]='-';
      }
      if (sm2_test_term3 > 0) {
        signature_str[1]='+';
      } else {
        signature_str[1]='-';
      }
      if (sm3_test_term3 > 0) {
        signature_str[2]='+';
      } else {
        signature_str[2]='-';
      }
      signature_str[3]=0;
      if (nle_config->smrfactor_1minus_enable == 1) {
        printf("status, Solved           phase 1 formula for coefficients, input sample: %lld, exponents:  %s, mixing polarity: %s, mass config: %s, sm3: %.14e, smrfactor mass: %.14Le, smrf: %s, term1 signature: %s, two-term test: %.14Le, c1: %.14Le, c2: %.14Le, c3: %.14Le, samples: %lld, ordering: %d, precision: %.3Le (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->nle_mixing_polarity_str, nle_state->smrfactor_mass_configuration_str, nle_state->input_sample_sm3, smrfactor_mass, smrf_str, signature_str, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]), c3_center[best_ordering], samples, best_ordering, best_precision_last, elapsed_time);
      } else {
        printf("status, Solved  phase 1 formula for coefficients, input sample: %lld, exponents:  %s, mixing polarity: -, sm3: %.14e, two-term test: %.14Le, c1: %.14Le, c2: %.14Le, c3: %.14Le, samples: %lld, ordering: %d, precision: %.3Le (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]), c3_center[best_ordering], samples, best_ordering, best_precision_last, elapsed_time);
      }
    } else if (nle_config->nle_mode == 3) {
      printf("status, Solved  phase 1 formula for coefficients, input sample: %lld, exponents:  %s, sm3: %.14e, samples: %lld, ordering: %d, precision: %.3Le (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, samples, best_ordering, precision, elapsed_time);
    }
    fflush(stdout); 

    if (nle_config->phase1_solution_detail == 1) {
      printf("status, +------------+------------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
    printf("status, | Exponents  |    Mass ratio    | NLE term |   Coefficient   |  C * term(sm1)  |  C * term(sm2)  |  C * term(sm3)  |\n");
      printf("status, +------------+------------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
      if (nle_config->nle_mode == 2) {
        if (nle_config->smrfactor_1minus_enable == 1) {
          if (nle_state->smrfactor_mass_configuration == 1) {
            printf("status, | %10s | 1-(smrf*M/%s) |  t1^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, smrfactor_mass_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
            printf("status, | %10s |    smrf*M/%s  |  t2^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, smrfactor_mass_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
          } else if (nle_state->smrfactor_mass_configuration == 0) {
            printf("status, | %10s | 1-(smrf*%s/M) |  t1^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, smrfactor_mass_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
            printf("status, | %10s |    smrf*%s/M) |  t2^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, smrfactor_mass_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
          }
        } else {
          printf("status, | %10s |       M/v        |  t1^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
          printf("status, | %10s |       M/v        |  t2^2    | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
        }
        printf("status, | %10s |                  | t1 * t2  | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c3_center[best_ordering], sm1_test_term3, sm2_test_term3, sm3_test_term3);
      } else if (nle_config->nle_mode == 3) {
        printf("status, | %10s |       M/v        |  term1   | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c1_center[best_ordering], sm1_test_term1, sm2_test_term1, sm3_test_term1);
        printf("status, | %10s |       M/v        |  term2   | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c2_center[best_ordering], sm1_test_term2, sm2_test_term2, sm3_test_term2);
        printf("status, | %10s |       M/v        |  term3   | %.9Le | %.9Le | %.9Le | %.9Le |\n", nle_state->exponents_str, c3_center[best_ordering], sm1_test_term3, sm2_test_term3, sm3_test_term3);
      }
      printf("status, +------------+------------------+----------+-----------------+-----------------+-----------------+-----------------+\n");
      fflush(stdout);

      // for debugging very long phase 1 solutions
      //exit(0);
    } // end if phase1_solution_detail

    // only run cscanner if mode=3 or two_term_test is an interesting integer match
    // for this interesting() check  we use a fixed filter of 3.  In cscanner two_term_test is evaluated with phase1_filter on c3 which may be more restrictive
    // This will help identify interesting geometries for further inspection
    if ((nle_config->nle_mode != 2) || ((two_term_test >= nle_config->phase1_two_term_test_min) && (two_term_test <= nle_config->phase1_two_term_test_max) && interesting(nle_config->phase1_filter, nle_config->phase1_int_match_max, nle_config->phase1_int_match_filter, two_term_test))) {
      matches_count_start=nle_state->phase1_matches_count;
      clock_gettime(CLOCK_REALTIME, &starttime);
      if (nle_config->nle_mode == 2) {
        // in two term mode take the square root of c1 and c2 and use two_term_test as c3
        nle_state->term1.coefficient=(double)sqrtl(c1_center[best_ordering]);
        nle_state->term2.coefficient=(double)sqrtl(c2_center[best_ordering]);
        nle_state->term3.coefficient=(double)two_term_test;
        if (nle_config->smrfactor_1minus_enable == 1) {
          sprintf(out_str_01, "status, Found interesting two-term test, exponents:  %s, mixing polarity: %s, mass config: %s, sm3: %.14e, smrfactor_mass: %.14Le, term1 signature: %s, two-term test: %.14Le, c1: %.14Le, c2: %.14Le, c3: %.14Le, smrf: %s", nle_state->exponents_str, nle_state->nle_mixing_polarity_str, nle_state->smrfactor_mass_configuration_str, nle_state->input_sample_sm3, smrfactor_mass, signature_str, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]), c3_center[best_ordering], smrf_str);
          printf("%s\n", out_str_01);
          fflush(stdout);
        } else {
          sprintf(out_str_01, "status, Found interesting two-term test, exponents:  %s, mixing polarity: -, sm3: %.14e, smrfactor_mass: %.14Le, two-term test: %.14Le, c1: %.14Le, c2: %.14Le", nle_state->exponents_str, nle_state->input_sample_sm3, smrfactor_mass, two_term_test, sqrtl(c1_center[best_ordering]), sqrtl(c2_center[best_ordering]));
          printf("%s\n", out_str_01);
          fflush(stdout);
        } // end if 1-minus
        if (nle_config->upload_results_enable == 1) {
          // upload interesting two_term_test as these are rare and significant
          sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_01, 320));
          system(exec_str);
        } // end if upload enable
      } else {
        // 3-term mode
        nle_state->term1.coefficient=(double)c1_center[best_ordering];
        nle_state->term2.coefficient=(double)c2_center[best_ordering];
        nle_state->term3.coefficient=(double)c3_center[best_ordering];
      } // end if mode==2

      if (nle_config->phase2_enable == 1) {
        // send coefficients (and two_term_test in 2-term mode) to factoring engine
        cscanner(nle_config, nle_state);

        clock_gettime(CLOCK_REALTIME, &endtime);
        matches_count_end=nle_state->phase1_matches_count;
        elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
        if (nle_config->phase1_status_enable ==1 ) {
          printf("status, Found %d interesting coefficient multipliers (%6.4fs)\n", (matches_count_end-matches_count_start), elapsed_time);
          fflush(stdout);
        }
      } else {
        // we solved p1 but phase 2 is disabled so we return failed code because cscanner was not run
        return(1);
      }
    } else {
      if ((nle_config->phase1_status_enable == 1) && (nle_config->phase2_enable == 1)) {
        printf("status, two-term test was out of range or not close enough to an interesting integer, skipping factoring process\n");
        fflush(stdout);
      }
      // we solved p1 but two_term_test was not interesting so return failed code because cscanner was not run
      return(1);
    }
  } else { // faled to solve, should only happen in 1-smr mode
#ifdef DEBUG10
    if (isnan(best_precision_last) == 1) {
      printf("debug, exponents: %s, samples: %10lld, time: %6.4fs, best_ordering: %d, progress: %6d, NaN exception\n", nle_state->exponents_str, samples, elapsed_time, best_ordering, progress[best_ordering]);
      fflush(stdout);
    }
#endif
    if (nle_config->phase1_status_enable == 1) {
      clock_gettime(CLOCK_REALTIME, &endtime);
      elapsed_time=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
      if (nle_config->smrfactor_1minus_enable == 1) {
        printf("status, Aborted solving  phase 1 formula for coefficients, input sample: %lld, exponents:  %s, mixing polarity: %s, mass config: %s, sm3: %.14e, smrfactor mass: %.14Le, smrf: %s, %s (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->nle_mixing_polarity_str, nle_state->smrfactor_mass_configuration_str, nle_state->input_sample_sm3, smrfactor_mass, smrf_str, unsolvable_exception_str, elapsed_time);
      } else {
        printf("status, Failed to solve  phase 1 formula for coefficients, input sample: %lld, exponents:  %s, sm3: %.14e (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, nle_state->input_sample_sm3, elapsed_time);
      }
      fflush(stdout);
    }
    // we failed, return 1
    return(1);
  } // end if best_precision_last
  // success, return 0
  return(0);
}
