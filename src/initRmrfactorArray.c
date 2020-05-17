#include <stdio.h>
#include <math.h> // pow, M_PI
#include <stdlib.h> // abs
#include <time.h>
#include "nle-lepton.h"
#include "util.h"

//#define DEBUG_RMRFACTOR

void initRmrfactorArray(nle_config_t *nle_config, nle_state_t *nle_state) {
  // store pre-computed static multiplier terms that go inside reference mass ratio on (1-rmr-smr) terms

  struct timespec start_time;
  struct timespec end_time;
  double elapsed_time;

  int uprmr=1, downrmr=1;
  int piuprmr=0, pidownrmr=1;
  int auprmr=0, adownrmr=1;
  int e2uprmr=0, e2downrmr=2;
  int useruprmr=0, userdownrmr=1;
  double updownrmr;
  double pirmr;
  double armr;
  double e2rmr;
  double userrmr;
  double rmrfactor;
  int upcomplexity, downcomplexity;

  unsigned int u, v;
  nle_rmrfactor_precomputed_t *multiplier;

  clock_gettime(CLOCK_REALTIME, &start_time);

  multiplier=nle_state->rmrfactors_precomputed_start;
  nle_state->rmrfactors_precomputed_count=0;
  for (uprmr=1; uprmr <= nle_config->rmrfactor_rational_max; uprmr++) {
    for (downrmr=1; downrmr <= nle_config->rmrfactor_rational_max; downrmr++) {
      if ((!nle_config->rmrfactor_rational_filter) || (uprmr <= 4) || ((downrmr == 1) && ((uprmr == 8) || (uprmr == 9) || (uprmr == 12) || (uprmr == 16) || (uprmr == 18) || (uprmr == 24) || (uprmr == 27) || (uprmr == 32)))) {
        if ((!nle_config->rmrfactor_rational_filter) || (downrmr <= 4) || ((uprmr == 1) && ((downrmr == 8) || (downrmr == 9) || (downrmr == 12) || (downrmr == 16) || (downrmr == 18) || (downrmr == 24) || (downrmr == 27) || (downrmr == 32)))) { 
          u=uprmr;
          v=downrmr;
          if (gcd(u, v) == 1) {
            updownrmr=(double)uprmr / (double)downrmr;

            for (e2uprmr=-nle_config->rmrfactor_2_exp_up_max; e2uprmr <= nle_config->rmrfactor_2_exp_up_max; e2uprmr++) {
              for (e2downrmr=1; e2downrmr <= nle_config->rmrfactor_2_exp_down_max; e2downrmr++) {
                u=abs(e2uprmr);
                v=e2downrmr;
                if (!((u == 1) && (v == 1)) && (gcd(u, v) == 1)) {  // extra checks to prevent 2 or 1/2
                  e2rmr=pow(2.0, ((double)e2uprmr / (double)e2downrmr));

                  for (piuprmr=-nle_config->rmrfactor_pi_exp_up_max; piuprmr <= nle_config->rmrfactor_pi_exp_up_max; piuprmr++) {
                    for (pidownrmr=1; pidownrmr <= nle_config->rmrfactor_pi_exp_down_max; pidownrmr++) {
                      u=abs(piuprmr);
                      v=pidownrmr;
                      if (gcd(u, v) == 1) {
                        pirmr=pow(M_PI, ((double)piuprmr / (double)pidownrmr));

                        for (auprmr=-nle_config->rmrfactor_alpha_exp_up_max; auprmr <= nle_config->rmrfactor_alpha_exp_up_max; auprmr++) {
                          for (adownrmr=1; adownrmr <= nle_config->rmrfactor_alpha_exp_down_max; adownrmr++) {
                            u=abs(auprmr);
                            v=adownrmr;
                            if (gcd(u, v) == 1) {
                              armr=pow(nle_config->ref_alpha_em, ((double)auprmr / (double)adownrmr));

                              for (useruprmr=-nle_config->rmrfactor_user_exp_up_max; useruprmr <= nle_config->rmrfactor_user_exp_up_max; useruprmr++) {
                                for (userdownrmr=1; userdownrmr <= nle_config->rmrfactor_user_exp_down_max; userdownrmr++) {
                                  u=abs(useruprmr);
                                  v=userdownrmr;
                                  if (gcd(u, v) == 1) {
                                    userrmr=pow(nle_config->rmrfactor_user, ((double)useruprmr / (double)userdownrmr));
                                    rmrfactor=updownrmr * e2rmr * armr * pirmr * userrmr;
#ifdef DEBUG_RMRFACTOR
                                    printf("debug, up: %d, down: %d, e2up: %d, e2down: %d, aup: %d, adown: %d, piup: %d, pidown: %d, userup: %d, userdown: %d, rmrf: %.9e, updown: %.9e, e2: %.9e, a: %.9e, pi: %.9e, user: %.9e\n", uprmr, downrmr, e2uprmr, e2downrmr, auprmr, adownrmr, piuprmr, pidownrmr, useruprmr, userdownrmr, rmrfactor, updownrmr, e2rmr, armr, pirmr, userrmr);
                                    fflush(stdout);
#endif
                                    multiplier->rmrfactor_rational_up=uprmr;
                                    multiplier->rmrfactor_rational_down=downrmr;
                                    multiplier->rmrfactor_2_exp_up=e2uprmr;
                                    multiplier->rmrfactor_2_exp_down=e2downrmr;
                                    multiplier->rmrfactor_alpha_exp_up=auprmr;
                                    multiplier->rmrfactor_alpha_exp_down=adownrmr;
                                    multiplier->rmrfactor_pi_exp_up=piuprmr;
                                    multiplier->rmrfactor_pi_exp_down=pidownrmr;
                                    multiplier->rmrfactor_user_exp_up=useruprmr;
                                    multiplier->rmrfactor_user_exp_down=userdownrmr;
                                    // ignore 1 on rationals
                                    if (uprmr == 1) {
                                      upcomplexity=0;
                                    } else {
                                      upcomplexity=uprmr;
                                    }
                                    if (downrmr == 1) {
                                      downcomplexity=0;
                                    } else {
                                      downcomplexity=downrmr;
                                    }
                                    multiplier->rmrfactor_complexity=\
                                           upcomplexity + downcomplexity\
                                           + abs(e2uprmr) + (e2downrmr-1)\
                                           + abs(auprmr) + (adownrmr-1)\
                                           + abs(piuprmr) + (pidownrmr-1)\
                                           + abs(useruprmr) + (userdownrmr-1);
                                    multiplier->rmrfactor_multiplier=rmrfactor; 
                                    initUses(&multiplier->rmrfactor_uses);
                                    if (auprmr != 0) {
                                      multiplier->rmrfactor_uses.alpha_em=1;
                                    }
                                    nle_state->rmrfactors_precomputed_count++;
#ifdef DEBUG_RMRFACTOR
                                    printf("debug, count: %d\n", nle_state->rmrfactors_precomputed_count);
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
    } // downrmr
  } // uprmr

  clock_gettime(CLOCK_REALTIME, &end_time);
  elapsed_time=((double)(end_time.tv_sec - 1500000000) + ((double)end_time.tv_nsec / 1.0E9)) - ((double)(start_time.tv_sec - 1500000000) + ((double)start_time.tv_nsec) / 1.0E9);

  printf("init, Initialized %d pre-computed multipliers for static factors inside reference mass ratio on (1-rmr-smr) terms (%6.4fs)\n", nle_state->rmrfactors_precomputed_count, elapsed_time);
  fflush(stdout);
}
