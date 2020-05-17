#include <stdio.h>
#include <math.h> // pow, M_PI
#include <stdlib.h> // abs
#include <time.h>
#include "nle-lepton.h"
#include "util.h"

//#define DEBUG_OUTFACTOR

void initOutfactorArray(nle_config_t *nle_config, nle_state_t *nle_state) {
  // store pre-computed static multiplier terms that go outside radical

  struct timespec start_time;
  struct timespec end_time;
  double elapsed_time;

  int upout=1, downout=1;
  int piupout=0, pidownout=1;
  int aupout=0, adownout=1;
  int e2upout=0, e2downout=2;
  int user1upout=0, user1downout=1;
  int user2upout=0, user2downout=1;
  int user3upout=0, user3downout=1;
  double updownout;
  double piout;
  double aout;
  double e2out;
  double user1out;
  double user2out;
  double user3out;
  double outfactor;
  int upcomplexity, downcomplexity;

  unsigned int u, v;
  nle_outfactor_precomputed_t *multiplier;

  clock_gettime(CLOCK_REALTIME, &start_time);

  multiplier=nle_state->outfactors_precomputed_start;
  nle_state->outfactors_precomputed_count=0;
  for (upout=1; upout <= nle_config->outfactor_rational_max; upout++) {
    for (downout=1; downout <= nle_config->outfactor_rational_max; downout++) {
      if ((!nle_config->outfactor_rational_filter) || (upout <= 4) || ((downout == 1) && ((upout == 8) || (upout == 9) || (upout == 12) || (upout == 16) || (upout == 18) || (upout == 24) || (upout == 27) || (upout == 32)))) {
        if ((!nle_config->outfactor_rational_filter) || (downout <= 4) || ((upout == 1) && ((downout == 8) || (downout == 9) || (downout == 12) || (downout == 16) || (downout == 18) || (downout == 24) || (downout == 27) || (downout == 32)))) { 
          u=upout;
          v=downout;
          if (gcd(u, v) == 1) {
            updownout=(double)upout / (double)downout;

            for (e2upout=-nle_config->outfactor_2_exp_up_max; e2upout <= nle_config->outfactor_2_exp_up_max; e2upout++) {
              for (e2downout=1; e2downout <= nle_config->outfactor_2_exp_down_max; e2downout++) {
                u=abs(e2upout);
                v=e2downout;
                if (!((u == 1) && (v == 1)) && (gcd(u, v) == 1)) {  // extra checks to prevent 2 or 1/2
                  e2out=pow(2.0, ((double)e2upout / (double)e2downout));

                  for (aupout=-nle_config->outfactor_alpha_exp_up_max; aupout <= nle_config->outfactor_alpha_exp_up_max; aupout++) {
                    for (adownout=1; adownout <= nle_config->outfactor_alpha_exp_down_max; adownout++) {
                      u=abs(aupout);
                      v=adownout;
                      if (gcd(u, v) == 1) {
                        aout=pow(nle_config->ref_alpha_em, ((double)aupout / (double)adownout));

                        for (piupout=-nle_config->outfactor_pi_exp_up_max; piupout <= nle_config->outfactor_pi_exp_up_max; piupout++) {
                          for (pidownout=1; pidownout <= nle_config->outfactor_pi_exp_down_max; pidownout++) {
                            u=abs(piupout);
                            v=pidownout;
                            if (gcd(u, v) == 1) {
                              piout=pow(M_PI, ((double)piupout / (double)pidownout));
 
                              for (user1upout=-nle_config->outfactor_user1_exp_up_max; user1upout <= nle_config->outfactor_user1_exp_up_max; user1upout++) {
                                for (user1downout=1; user1downout <= nle_config->outfactor_user1_exp_down_max; user1downout++) {
                                  u=abs(user1upout);
                                  v=user1downout;
                                  if (gcd(u, v) == 1) {
                                    user1out=pow(nle_config->outfactor_user1, ((double)user1upout / (double)user1downout));

                                    for (user2upout=-nle_config->outfactor_user2_exp_up_max; user2upout <= nle_config->outfactor_user2_exp_up_max; user2upout++) {
                                      for (user2downout=1; user2downout <= nle_config->outfactor_user2_exp_down_max; user2downout++) {
                                        u=abs(user2upout);
                                        v=user2downout;
                                        if (gcd(u, v) == 1) {
                                          user2out=pow(nle_config->outfactor_user2, ((double)user2upout / (double)user2downout));

                                          for (user3upout=-nle_config->outfactor_user3_exp_up_max; user3upout <= nle_config->outfactor_user3_exp_up_max; user3upout++) {
                                            for (user3downout=1; user3downout <= nle_config->outfactor_user3_exp_down_max; user3downout++) {
                                              u=abs(user3upout);
                                              v=user3downout;
                                              if (gcd(u, v) == 1) {
                                                user3out=pow(nle_config->outfactor_user3, ((double)user3upout / (double)user3downout));

                                                outfactor=updownout * e2out * aout * piout * user1out * user2out * user3out;
#ifdef DEBUG_OUTFACTOR
                                                printf("debug, outfactorInit, up: %d, down: %d, e2up: %d, e2down: %d, aup: %d, adown: %d, piup: %d, pidown: %d, use1rup: %d, user1down: %d, use2rup: %d, user2down: %d, use3rup: %d, user3down: %d, outfactor: %.9e, updown: %.9e, e2: %.9e, a: %.9e, pi: %.9e, user1: %.9e, user2: %.9e, user3: %.9e\n", upout, downout, e2upout, e2downout, aupout, adownout, piupout, pidownout, user1upout, user1downout, user2upout, user2downout, user3upout, user3downout, outfactor, updownout, e2out, aout, piout, user1out, user2out, user3out);
                                                fflush(stdout);
#endif
                                                multiplier->outfactor_rational_up=upout;
                                                multiplier->outfactor_rational_down=downout;
                                                multiplier->outfactor_2_exp_up=e2upout;
                                                multiplier->outfactor_2_exp_down=e2downout;
                                                multiplier->outfactor_alpha_exp_up=aupout;
                                                multiplier->outfactor_alpha_exp_down=adownout;
                                                multiplier->outfactor_pi_exp_up=piupout;
                                                multiplier->outfactor_pi_exp_down=pidownout;
                                                multiplier->outfactor_user1_exp_up=user1upout;
                                                multiplier->outfactor_user1_exp_down=user1downout;
                                                multiplier->outfactor_user2_exp_up=user2upout;
                                                multiplier->outfactor_user2_exp_down=user2downout;
                                                multiplier->outfactor_user3_exp_up=user3upout;
                                                multiplier->outfactor_user3_exp_down=user3downout;
                                                // ignore 1 on rationals
                                                if (upout == 1) {
                                                  upcomplexity=0;
                                                } else {
                                                  upcomplexity=upout;
                                                }
                                                if (downout == 1) {
                                                  downcomplexity=0;
                                                } else {
                                                  downcomplexity=downout;
                                                }
                                                multiplier->outfactor_complexity=\
                                                       upcomplexity + downcomplexity\
                                                       + abs(e2upout) + (e2downout-1)\
                                                       + abs(aupout) + (adownout-1)\
                                                       + abs(piupout) + (pidownout-1)\
                                                       + abs(user1upout) + (user1downout-1)\
                                                       + abs(user2upout) + (user2downout-1)\
                                                       + abs(user3upout) + (user3downout-1);
                                                multiplier->outfactor_multiplier=outfactor; 
                                                initUses(&multiplier->outfactor_uses);
                                                if (aupout != 0) {
                                                  multiplier->outfactor_uses.alpha_em=1;
                                                }
                                                nle_state->outfactors_precomputed_count++;
#ifdef DEBUG_OUTFACTOR
                                                printf("debug, outfactorInit, count: %d\n", nle_state->outfactors_precomputed_count);
                                                fflush(stdout);
#endif
                                                multiplier++;
                                              } // gcd user3
                                            } // down user3
                                          } // up user3
                                        } // gcd user2
                                      } // down user2
                                    } // up user2
                                  } // gcd user1
                                } // down user1
                              } // up user1
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
    } // downout
  } // upout

  clock_gettime(CLOCK_REALTIME, &end_time);
  elapsed_time=((double)(end_time.tv_sec - 1500000000) + ((double)end_time.tv_nsec / 1.0E9)) - ((double)(start_time.tv_sec - 1500000000) + ((double)start_time.tv_nsec) / 1.0E9);

  printf("init, Initialized %d pre-computed multipliers for static factors outside radical (%6.4fs)\n", nle_state->outfactors_precomputed_count, elapsed_time);
  fflush(stdout);
}
