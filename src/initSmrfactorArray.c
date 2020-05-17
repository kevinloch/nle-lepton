#include <stdio.h>
#include <math.h> // pow, M_PI
#include <stdlib.h> // abs
#include <time.h>
#include "nle-lepton.h"
#include "util.h"

//#define DEBUG_SMRFACTOR

void initSmrfactorArray(nle_config_t *nle_config, nle_state_t *nle_state) {
  // store pre-computed static multiplier terms that go inside solution mass ratio on (1-smr) terms

  struct timespec start_time;
  struct timespec end_time;
  double elapsed_time;

  int upsmr=1, downsmr=1;
  int piupsmr=0, pidownsmr=1;
  int aupsmr=0, adownsmr=1;
  int e2upsmr=0, e2downsmr=2;
  int userupsmr=0, userdownsmr=1;
  double updownsmr;
  double pismr;
  double asmr;
  double e2smr;
  double usersmr;
  double smrfactor;
  int upcomplexity, downcomplexity;

  unsigned int u, v;
  nle_smrfactor_precomputed_t *multiplier;

  clock_gettime(CLOCK_REALTIME, &start_time);

  multiplier=nle_state->smrfactors_precomputed_start;
  nle_state->smrfactors_precomputed_count=0;
  for (upsmr=1; upsmr <= nle_config->smrfactor_rational_max; upsmr++) {
    for (downsmr=1; downsmr <= nle_config->smrfactor_rational_max; downsmr++) {
      if ((!nle_config->smrfactor_rational_filter) || (upsmr <= 4) || ((downsmr == 1) && ((upsmr == 8) || (upsmr == 9) || (upsmr == 12) || (upsmr == 16) || (upsmr == 18) || (upsmr == 24) || (upsmr == 27) || (upsmr == 32)))) {
        if ((!nle_config->smrfactor_rational_filter) || (downsmr <= 4) || ((upsmr == 1) && ((downsmr == 8) || (downsmr == 9) || (downsmr == 12) || (downsmr == 16) || (downsmr == 18) || (downsmr == 24) || (downsmr == 27) || (downsmr == 32)))) { 
          u=upsmr;
          v=downsmr;
          if (gcd(u, v) == 1) {
            updownsmr=(double)upsmr / (double)downsmr;

            for (e2upsmr=-nle_config->smrfactor_2_exp_up_max; e2upsmr <= nle_config->smrfactor_2_exp_up_max; e2upsmr++) {
              for (e2downsmr=1; e2downsmr <= nle_config->smrfactor_2_exp_down_max; e2downsmr++) {
                u=abs(e2upsmr);
                v=e2downsmr;
                if (!((u == 1) && (v == 1)) && (gcd(u, v) == 1)) {  // extra checks to prevent 2 or 1/2
                  e2smr=pow(2.0, ((double)e2upsmr / (double)e2downsmr));

                  for (piupsmr=-nle_config->smrfactor_pi_exp_up_max; piupsmr <= nle_config->smrfactor_pi_exp_up_max; piupsmr++) {
                    for (pidownsmr=1; pidownsmr <= nle_config->smrfactor_pi_exp_down_max; pidownsmr++) {
                      u=abs(piupsmr);
                      v=pidownsmr;
                      if (gcd(u, v) == 1) {
                        pismr=pow(M_PI, ((double)piupsmr / (double)pidownsmr));

                        for (aupsmr=-nle_config->smrfactor_alpha_exp_up_max; aupsmr <= nle_config->smrfactor_alpha_exp_up_max; aupsmr++) {
                          for (adownsmr=1; adownsmr <= nle_config->smrfactor_alpha_exp_down_max; adownsmr++) {
                            u=abs(aupsmr);
                            v=adownsmr;
                            if (gcd(u, v) == 1) {
                              asmr=pow(nle_config->ref_alpha_em, ((double)aupsmr / (double)adownsmr));

                              for (userupsmr=-nle_config->smrfactor_user_exp_up_max; userupsmr <= nle_config->smrfactor_user_exp_up_max; userupsmr++) {
                                for (userdownsmr=1; userdownsmr <= nle_config->smrfactor_user_exp_down_max; userdownsmr++) {
                                  u=abs(userupsmr);
                                  v=userdownsmr;
                                  if (gcd(u, v) == 1) {
                                    usersmr=pow(nle_config->smrfactor_user, ((double)userupsmr / (double)userdownsmr));
                                    smrfactor=updownsmr * e2smr * asmr * pismr * usersmr;
#ifdef DEBUG_SMRFACTOR
                                    printf("debug, up: %d, down: %d, e2up: %d, e2down: %d, aup: %d, adown: %d, piup: %d, pidown: %d, userup: %d, userdown: %d, smrf: %.9e, updown: %.9e, e2: %.9e, a: %.9e, pi: %.9e, user: %.9e\n", upsmr, downsmr, e2upsmr, e2downsmr, aupsmr, adownsmr, piupsmr, pidownsmr, userupsmr, userdownsmr, smrfactor, updownsmr, e2smr, asmr, pismr, usersmr);
                                    fflush(stdout);
#endif
                                    multiplier->smrfactor_rational_up=upsmr;
                                    multiplier->smrfactor_rational_down=downsmr;
                                    multiplier->smrfactor_2_exp_up=e2upsmr;
                                    multiplier->smrfactor_2_exp_down=e2downsmr;
                                    multiplier->smrfactor_alpha_exp_up=aupsmr;
                                    multiplier->smrfactor_alpha_exp_down=adownsmr;
                                    multiplier->smrfactor_pi_exp_up=piupsmr;
                                    multiplier->smrfactor_pi_exp_down=pidownsmr;
                                    multiplier->smrfactor_user_exp_up=userupsmr;
                                    multiplier->smrfactor_user_exp_down=userdownsmr;
                                    // ignore 1 on rationals
                                    if (upsmr == 1) {
                                      upcomplexity=0;
                                    } else {
                                      upcomplexity=upsmr;
                                    }
                                    if (downsmr == 1) {
                                      downcomplexity=0;
                                    } else {
                                      downcomplexity=downsmr;
                                    }
                                    multiplier->smrfactor_complexity=\
                                           upcomplexity + downcomplexity\
                                           + abs(e2upsmr) + (e2downsmr-1)\
                                           + abs(aupsmr) + (adownsmr-1)\
                                           + abs(piupsmr) + (pidownsmr-1)\
                                           + abs(userupsmr) + (userdownsmr-1);
                                    multiplier->smrfactor_multiplier=smrfactor; 
                                    initUses(&multiplier->smrfactor_uses);
                                    if (aupsmr != 0) {
                                      multiplier->smrfactor_uses.alpha_em=1;
                                    }
                                    nle_state->smrfactors_precomputed_count++;
#ifdef DEBUG_SMRFACTOR
                                    printf("debug, count: %d\n", nle_state->smrfactors_precomputed_count);
                                    fflush(stdout);
#endif
                                    multiplier++;
                                  } // gcd user
                                } // down user
                              } // up user
                            } // gcd a
                          } // adown
                        } // aup
                      } // gcd pi
                    } // pidown
                  } // piup
                } // gcd e2
              } // e2down
            } // e2up
          } // gcd updown
        } // downfilter
      } // upfilter
    } // downsmr
  } // upsmr

  clock_gettime(CLOCK_REALTIME, &end_time);
  elapsed_time=((double)(end_time.tv_sec - 1500000000) + ((double)end_time.tv_nsec / 1.0E9)) - ((double)(start_time.tv_sec - 1500000000) + ((double)start_time.tv_nsec) / 1.0E9);

  printf("init, Initialized %d pre-computed multipliers for static factors inside solution mass ratio on (1-smr) terms (%6.4fs)\n", nle_state->smrfactors_precomputed_count, elapsed_time);
  fflush(stdout);
}
