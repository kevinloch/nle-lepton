#include <stdio.h>
#include <math.h> // pow, M_PI
#include <stdlib.h> // abs
#include <time.h>
#include "nle-lepton.h"
#include "util.h"

//#define DEBUG_DYNAMICFACTOR

void initDynamicfactorArray(nle_config_t *nle_config, nle_state_t *nle_state, nle_dynamicfactor_precomputed_t *dynamicfactors_precomputed_start, int *dynamicfactors_precomputed_count) {
  // store pre-computed dynamic multiplier terms that go outside the radical

  struct timespec start_time;
  struct timespec end_time;

  double elapsed_time;
  unsigned int u, v;
  int rmr_mass_enabled[6];
  int sin2w_exp_up=0, sin2w_exp_down=1;
  int cos2w_exp_up=0, cos2w_exp_down=1;
  int rmr_exp_up=0, rmr_exp_down=1;
  int rmr_mass_id_up=0, rmr_mass_id_down=0;
  double rmr_mass_up=1.0;
  double rmr_mass_down=1.0;
  double outfactor_sin2w=1.0;
  double outfactor_cos2w=1.0;
  double outfactor_rmr=1.0;
  double sin2w;
  double cos2w;
  double dynamicfactor;
  nle_dynamicfactor_precomputed_t *multiplier;
  int skip;

  clock_gettime(CLOCK_REALTIME, &start_time);

  sin2w=nle_state->input_sample_sin2w;
  cos2w=1.0 - sin2w;

  /*
    Index of mass id's and potential output variables
    0:  G
    1:  v
    2:  mz
    3:  mw
    4:  mH0
    5:  m_user
    6:  sm1
    7:  sm2
    8:  sm3
    9:  sin2w
    10: alpha_em
    11: alpha_w
  */

  // determine which rmr masses are enabled
  if (nle_config->outfactor_rmr_mp_enable == 1) {
    rmr_mass_enabled[0]=1;
  } else {
    rmr_mass_enabled[0]=0;
  }
  if (nle_config->outfactor_rmr_v_enable == 1) {
    rmr_mass_enabled[1]=1;
  } else {
    rmr_mass_enabled[1]=0;
  }
  if (nle_config->outfactor_rmr_mz_enable == 1) {
    rmr_mass_enabled[2]=1;
  } else {
    rmr_mass_enabled[2]=0;
  }
  if (nle_config->outfactor_rmr_mw_enable == 1) {
    rmr_mass_enabled[3]=1;
  } else {
    rmr_mass_enabled[3]=0;
  }
  if (nle_config->outfactor_rmr_mh0_enable == 1) {
    rmr_mass_enabled[4]=1;
  } else {
    rmr_mass_enabled[4]=0;
  }
  if (nle_config->outfactor_rmr_user_enable == 1) {
    rmr_mass_enabled[5]=1;
  } else {
    rmr_mass_enabled[5]=0;
  }

  multiplier=dynamicfactors_precomputed_start;
  *dynamicfactors_precomputed_count=0;

  for (sin2w_exp_up=-nle_config->outfactor_weak_exp_up_max; sin2w_exp_up <= nle_config->outfactor_weak_exp_up_max; sin2w_exp_up++) {
    for (sin2w_exp_down=1; sin2w_exp_down <= nle_config->outfactor_weak_exp_down_max; sin2w_exp_down++) {
      u=abs(sin2w_exp_up);
      v=sin2w_exp_down;
      if (gcd(u, v) == 1) {
        outfactor_sin2w=pow(sin2w, ((float)sin2w_exp_up / (float)sin2w_exp_down));

        for (cos2w_exp_up=-nle_config->outfactor_weak_exp_up_max; cos2w_exp_up <= nle_config->outfactor_weak_exp_up_max; cos2w_exp_up++) {
          for (cos2w_exp_down=1; cos2w_exp_down <= nle_config->outfactor_weak_exp_down_max; cos2w_exp_down++) {
            u=abs(cos2w_exp_up);
            v=cos2w_exp_down;
            if (gcd(u, v) == 1) {
              outfactor_cos2w=pow(cos2w, ((float)cos2w_exp_up / (float)cos2w_exp_down));

              for (rmr_exp_up=0; rmr_exp_up <= nle_config->outfactor_rmr_exp_up_max; rmr_exp_up++) {
                for (rmr_exp_down=1; rmr_exp_down <= nle_config->outfactor_rmr_exp_down_max; rmr_exp_down++) {
                  u=rmr_exp_up;
                  v=rmr_exp_down;
                  if (gcd(u, v) == 1) {
                    skip=0;
                    for (rmr_mass_id_up=0; rmr_mass_id_up <= 5; rmr_mass_id_up++) {
                      for (rmr_mass_id_down=0; rmr_mass_id_down <= 5; rmr_mass_id_down++) {
                        if ((rmr_mass_enabled[rmr_mass_id_up] == 1) && (rmr_mass_enabled[rmr_mass_id_down] == 1) && (rmr_mass_id_up != rmr_mass_id_down)\
                           && ((rmr_exp_up != 0) || (skip == 0))) {
                          if (rmr_exp_up == 0) {
                            skip=1;
                          }
                          if (rmr_mass_id_up == 0) {
                            rmr_mass_up=nle_state->input_sample_mp;
                          } else if (rmr_mass_id_up == 1) {
                            rmr_mass_up=nle_state->input_sample_v;
                          } else if (rmr_mass_id_up == 2) {
                            rmr_mass_up=nle_state->input_sample_mz;
                          } else if (rmr_mass_id_up == 3) {
                            rmr_mass_up=nle_state->input_sample_mw;
                          } else if (rmr_mass_id_up == 4) {
                            rmr_mass_up=nle_state->input_sample_mh0;
                          } else if (rmr_mass_id_up == 5) {
                            rmr_mass_up=nle_state->input_sample_muser;
                          }
                          if (rmr_mass_id_down == 0) {
                            rmr_mass_down=nle_state->input_sample_mp;
                          } else if (rmr_mass_id_down == 1) {
                            rmr_mass_down=nle_state->input_sample_v;
                          } else if (rmr_mass_id_down == 2) {
                            rmr_mass_down=nle_state->input_sample_mz;
                          } else if (rmr_mass_id_down == 3) {
                            rmr_mass_down=nle_state->input_sample_mw;
                          } else if (rmr_mass_id_down == 4) {
                            rmr_mass_down=nle_state->input_sample_mh0;
                          } else if (rmr_mass_id_down == 5) {
                            rmr_mass_down=nle_state->input_sample_muser;
                          }
                          outfactor_rmr=pow((rmr_mass_up / rmr_mass_down), ((float)rmr_exp_up / (float)rmr_exp_down));

                          dynamicfactor=outfactor_sin2w * outfactor_cos2w * outfactor_rmr;

                          multiplier->outfactor_sin2w_exp_up=sin2w_exp_up;
                          multiplier->outfactor_sin2w_exp_down=sin2w_exp_down;
                          multiplier->outfactor_cos2w_exp_up=cos2w_exp_up;
                          multiplier->outfactor_cos2w_exp_down=cos2w_exp_down;
                          multiplier->outfactor_rmr_exp_up=rmr_exp_up;
                          multiplier->outfactor_rmr_exp_down=rmr_exp_down;
                          multiplier->outfactor_rmr_mass_id_up=rmr_mass_id_up;
                          multiplier->outfactor_rmr_mass_id_down=rmr_mass_id_down;
                          multiplier->dynamicfactor_complexity=\
                                              + (abs(sin2w_exp_up) + (sin2w_exp_down-1))\
                                              + (abs(cos2w_exp_up) + (cos2w_exp_down-1))\
                                              + (abs(rmr_exp_up) + (rmr_exp_down-1));
                          multiplier->dynamicfactor_multiplier=dynamicfactor;
                          initUses(&multiplier->dynamicfactor_uses);

                          if ((sin2w_exp_up != 0) || (cos2w_exp_up != 0)) {
                            multiplier->dynamicfactor_uses.sin2w=1;
                          }
                          if (rmr_exp_up != 0) {
                            if ((rmr_mass_id_up == 0) || (rmr_mass_id_down == 0)) {
                              multiplier->dynamicfactor_uses.G=1;
                            }
                            if ((rmr_mass_id_up == 1) || (rmr_mass_id_down == 1)) {
                              multiplier->dynamicfactor_uses.v=1;
                            }
                            if ((rmr_mass_id_up == 2) || (rmr_mass_id_down == 2)) { 
                              multiplier->dynamicfactor_uses.mz=1;
                            }
                            if ((rmr_mass_id_up == 3) || (rmr_mass_id_down == 3)) {
                              multiplier->dynamicfactor_uses.mw=1;
                            }
                            if ((rmr_mass_id_up == 4) || (rmr_mass_id_down == 4)) {
                              multiplier->dynamicfactor_uses.mh0=1;
                            }
                            if ((rmr_mass_id_up == 5) || (rmr_mass_id_down == 5)) {
                              multiplier->dynamicfactor_uses.m_user=1;
                            }
                          }
                          *dynamicfactors_precomputed_count=*dynamicfactors_precomputed_count + 1;
#ifdef DEBUG_DYNAMICFACTOR
                          printf("debug, count: %d, sin2w_exp_up: %d, sin2w_exp_down: %d, cos2w_exp_up: %d, cos2w_exp_down: %d, rmr_exp_up: %d, rmr_exp_down: %d, rmr_mass_id_up: %d, rmr_mass_id_down: %d, dynamicfactor_multiplier: %.9e\n", *dynamicfactors_precomputed_count, sin2w_exp_up, sin2w_exp_down, cos2w_exp_up, cos2w_exp_down, rmr_exp_up, rmr_exp_down, rmr_mass_id_up, rmr_mass_id_down, multiplier->dynamicfactor_multiplier);
                          fflush(stdout);
#endif
                          multiplier++;
                        } // rmr_mass_id !=
                      } // rmr_mass_id_down
                    } // rmr_mass_id_up
                  } // gcd rmr
                } // rmr down
              } // rmr_up
            } // gcd outfactor_cos2w
          } // cos2w_exp_down
        } // cos2w_exp_up
      } // gcd outfactor_sin2w
    } // sin2w_exp_down
  } // sin2w_exp_up

  clock_gettime(CLOCK_REALTIME, &end_time);
  elapsed_time=((double)(end_time.tv_sec - 1500000000) + ((double)end_time.tv_nsec / 1.0E9)) - ((double)(start_time.tv_sec - 1500000000) + ((double)start_time.tv_nsec) / 1.0E9);

  printf("init, Initialized %d pre-computed multipliers for dynamic factors outside radical (%6.4fs)\n", *dynamicfactors_precomputed_count, elapsed_time);
  fflush(stdout);
}
