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
  int userupout=0, userdownout=1;
  double updownout;
  double piout;
  double aout;
  double e2out;
  double userout;
  double outfactor;

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
                  e2out=pow(2.0, ((float)e2upout / (float)e2downout));

                  for (aupout=-nle_config->outfactor_alpha_exp_up_max; aupout <= nle_config->outfactor_alpha_exp_up_max; aupout++) {
                    for (adownout=1; adownout <= nle_config->outfactor_alpha_exp_down_max; adownout++) {
                      u=abs(aupout);
                      v=adownout;
                      if (gcd(u, v) == 1) {
                        aout=pow(nle_config->ref_alpha, ((float)aupout / (float)adownout));

                        for (piupout=-nle_config->outfactor_pi_exp_up_max; piupout <= nle_config->outfactor_pi_exp_up_max; piupout++) {
                          for (pidownout=1; pidownout <= nle_config->outfactor_pi_exp_down_max; pidownout++) {
                            u=abs(piupout);
                            v=pidownout;
                            if (gcd(u, v) == 1) {
                              piout=pow(M_PI, ((float)piupout / (float)pidownout));
 
                              for (userupout=-nle_config->outfactor_user_exp_up_max; userupout <= nle_config->outfactor_user_exp_up_max; userupout++) {
                                for (userdownout=1; userdownout <= nle_config->outfactor_user_exp_down_max; userdownout++) {
                                  u=abs(userupout);
                                  v=userdownout;
                                  if (gcd(u, v) == 1) {
                                    userout=pow(nle_config->outfactor_user, ((float)userupout / (float)userdownout));
                                    outfactor=updownout * e2out * aout * piout * userout;
#ifdef DEBUG_OUTFACTOR
                                    printf("debug, up: %d, down: %d, e2up: %d, e2down: %d, aup: %d, adown: %d, piup: %d, pidown: %d, userup: %d, userdown: %d, outfactor: %.9e, updown: %.9e, e2: %.9e, a: %.9e, pi: %.9e, user: %.9e\n", upout, downout, e2upout, e2downout, aupout, adownout, piupout, pidownout, userupout, userdownout, outfactor, updownout, e2out, aout, piout, userout);
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
                                    multiplier->outfactor_user_exp_up=userupout;
                                    multiplier->outfactor_user_exp_down=userdownout;
                                    multiplier->outfactor_complexity=\
                                           upout + downout\
                                           + abs(e2upout) + (e2downout-1)\
                                           + abs(aupout) + (adownout-1)\
                                           + abs(piupout) + (pidownout-1)\
                                           + abs(userupout) + (userdownout-1);
                                    multiplier->outfactor_multiplier=outfactor; 
                                    initUses(&multiplier->outfactor_uses);
                                    if (aupout != 0) {
                                      multiplier->outfactor_uses.alpha_em=1;
                                    }
                                    nle_state->outfactors_precomputed_count++;
#ifdef DEBUG_OUTFACTOR
                                    printf("debug, count: %d\n", nle_state->outfactors_precomputed_count);
                                    fflush(stdout);
#endif
                                    multiplier++;
                                  } // gcd user
                                } // down user
                              } // up user
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
