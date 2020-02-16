#include <stdio.h>
#include <stdlib.h> // abs
#include <math.h>   // pow
#include "nle-lepton.h"
#include "util.h"

//#define DEBUG_CSCANNER

void cscanner(nle_config_t *nle_config, nle_state_t *nle_state) {
  //  Each coefficient is multiplied by various numbers and the resulting value is tested to see if it is
  // close to an interesting integer or simple rational number.  The results are stored in the match table.
  int infactor, outfactor;
  double multiplier;
  nle_infactor_precomputed_t *infactors;
  nle_outfactor_precomputed_t *outfactors;
  unsigned int u, v;
  nle_phase1_match_t *match;
  int phase1_filter;
  int max_int;
  int filter_int;
  int mass_ratio_enabled;

  int mass_ratio_id;
  double mass_ratio;
  double term1_mass_ratio_factor;
  double term2_mass_ratio_factor;
  double term3_mass_ratio_factor;
  int sin2w_exp_up=0, sin2w_exp_down=1;
  int cos2w_exp_up=0, cos2w_exp_down=1;
  double outfactor_sin2w=1.0;
  double outfactor_cos2w=1.0;
  double sin2w;
  double cos2w;

  sin2w=nle_state->random_sample_sin2w;
  cos2w=1.0 - sin2w;

  phase1_filter=nle_config->phase1_filter;
  max_int=nle_config->phase1_int_match_max;
  filter_int=nle_config->phase1_int_match_filter;

  match=nle_state->phase1_matches_start;
  // here we substitute the M/v mass ratio used in phase 1 with the actual test mass ratio
  for (mass_ratio_id=0; mass_ratio_id<=5; mass_ratio_id++) {
    mass_ratio_enabled=1;
    if (mass_ratio_id == 0) {
      if (nle_config->smrfactor_mass_mp_enable == 1) {
        mass_ratio=nle_state->random_sample_mp/nle_config->ref_v;
        if (nle_config->status_enable ==1) {
          printf(" M/mp");
          fflush(stdout);
        }
      } else {
        mass_ratio_enabled=0;
      }
    } else if (mass_ratio_id == 1) {
      if (nle_config->smrfactor_mass_v_enable == 1) {
        mass_ratio=1.0;
        if (nle_config->status_enable ==1) {
          printf(" M/v");              
          fflush(stdout);               
        }
      } else {
        mass_ratio_enabled=0;
      }
    } else if (mass_ratio_id == 2) {
      if (nle_config->smrfactor_mass_mz_enable == 1) {
        mass_ratio=nle_state->random_sample_mz/nle_config->ref_v;
        if (nle_config->status_enable ==1) {
          printf(" M/mz");             
          fflush(stdout);               
        }
      } else {
        mass_ratio_enabled=0;
      } 
    } else if (mass_ratio_id == 3) {
      if (nle_config->smrfactor_mass_mw_enable == 1) {
        mass_ratio=nle_state->random_sample_mw/nle_config->ref_v;
        if (nle_config->status_enable ==1) {
          printf(" M/mw");             
          fflush(stdout);               
        }
      } else {
        mass_ratio_enabled=0;
      }
    } else if (mass_ratio_id == 4) {
      if (nle_config->smrfactor_mass_mh0_enable == 1) {
        mass_ratio=nle_state->random_sample_mh0/nle_config->ref_v;
        if (nle_config->status_enable ==1) {
          printf(" M/mh0");             
          fflush(stdout);               
        }
      } else {
        mass_ratio_enabled=0;
      }
    } else if (mass_ratio_id == 5) {
      if ((nle_config->smrfactor_mass_user_enable == 1) && (mass_ratio_id == 5)) {
        mass_ratio=nle_state->random_sample_muser/nle_config->ref_v;
        if (nle_config->status_enable ==1) {
          printf(" M/m_user");
          fflush(stdout);
        }
      } else {
        mass_ratio_enabled=0;
      }
    }                             

    if (mass_ratio_enabled == 1) {
      term1_mass_ratio_factor=pow(mass_ratio, (1.0 / (double)nle_state->term1.exp_inv));
      term2_mass_ratio_factor=pow(mass_ratio, (1.0 / (double)nle_state->term2.exp_inv));
      term3_mass_ratio_factor=pow(mass_ratio, (1.0 / (double)nle_state->term3.exp_inv));

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
#ifdef DEBUG_CSCANNER
                  printf("sin2w_exp_up: %d, sin2w_exp_down: %d, cos2w_exp_up: %d, cos2w_exp_down: %d\n", sin2w_exp_up, sin2w_exp_down, cos2w_exp_up, cos2w_exp_down);
                  fflush(stdout);
#endif
                  outfactors=nle_state->outfactors_precomputed_start;
                  for (outfactor=0; outfactor < nle_state->outfactors_precomputed_count; outfactor++) {

                    infactors=nle_state->infactors_precomputed_start;
                    for (infactor=0; infactor < nle_state->infactors_precomputed_count; infactor++) { 

                      // test multiplier against term1 coefficient
                      if ((sin2w_exp_down == 1) || (sin2w_exp_down == 2) || (sin2w_exp_down == nle_state->term1.exp_inv) || (sin2w_exp_down == (nle_state->term1.exp_inv * 2))) {
                        if ((cos2w_exp_down == 1) || (cos2w_exp_down == 2) || (cos2w_exp_down == nle_state->term1.exp_inv) || (cos2w_exp_down == (nle_state->term1.exp_inv * 2))) {
                        multiplier=nle_state->term1.coefficient * term1_mass_ratio_factor * outfactor_sin2w * outfactor_cos2w * outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term1.exp_inv)];
                          if (interesting(phase1_filter, max_int, filter_int, multiplier)) {
                            nle_state->terms_matched[0]=nle_state->term1.exp_inv;
                            match->exp_inv=nle_state->term1.exp_inv;
                            match->smrfactor_mass=mass_ratio_id;
                            match->infactor_rational_up=infactors->infactor_rational_up;
                            match->infactor_rational_down=infactors->infactor_rational_down;
                            match->infactor_2_exp_up=infactors->infactor_2_exp_up;
                            match->infactor_2_exp_down=infactors->infactor_2_exp_down;
                            match->infactor_alpha_exp_up=infactors->infactor_alpha_exp_up;
                            match->infactor_alpha_exp_down=infactors->infactor_alpha_exp_down;
                            match->infactor_pi_exp_up=infactors->infactor_pi_exp_up;
                            match->infactor_pi_exp_down=infactors->infactor_pi_exp_down;
                            match->infactor_nss=infactors->infactor_nss;
                            match->infactor_nbv=infactors->infactor_nbv;
                            match->infactor_user_exp_up=infactors->infactor_user_exp_up;
                            match->infactor_user_exp_down=infactors->infactor_user_exp_down;
                            match->outfactor_rational_up=outfactors->outfactor_rational_up;
                            match->outfactor_rational_down=outfactors->outfactor_rational_down;
                            match->outfactor_2_exp_up=outfactors->outfactor_2_exp_up;
                            match->outfactor_2_exp_down=outfactors->outfactor_2_exp_down;
                            match->outfactor_alpha_exp_up=outfactors->outfactor_alpha_exp_up;
                            match->outfactor_alpha_exp_down=outfactors->outfactor_alpha_exp_down;
                            match->outfactor_pi_exp_up=outfactors->outfactor_pi_exp_up;
                            match->outfactor_pi_exp_down=outfactors->outfactor_pi_exp_down;
                            match->outfactor_sin2w_exp_up=sin2w_exp_up;
                            match->outfactor_sin2w_exp_down=sin2w_exp_down;
                            match->outfactor_cos2w_exp_up=cos2w_exp_up;
                            match->outfactor_cos2w_exp_down=cos2w_exp_down;
                            match->outfactor_user_exp_up=outfactors->outfactor_user_exp_up;
                            match->outfactor_user_exp_down=outfactors->outfactor_user_exp_down;
                            match->static_multiplier=outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term1.exp_inv)]; // excludes sin2w and cos2w that have uncertainty
                            match->match=multiplier;
                            match->match_complexity=outfactors->outfactor_complexity\
                                  + infactors->infactor_complexity\
                                  + (abs(sin2w_exp_up) + (sin2w_exp_down-1))\
                                  + (abs(cos2w_exp_up) + (cos2w_exp_down-1));
                            initUses(&match->match_uses);
                            addUses(&match->match_uses, &outfactors->outfactor_uses);
                            addUses(&match->match_uses, &infactors->infactor_uses);
                            if (mass_ratio_id == 0) {
                              match->match_uses.G=1;
                            } else if (mass_ratio_id == 1) {
                              match->match_uses.v=1;
                            } else if (mass_ratio_id == 2) {
                              match->match_uses.mz=1;
                            } else if (mass_ratio_id == 3) {
                              match->match_uses.mw=1;
                            } else if (mass_ratio_id == 4) {
                              match->match_uses.mh0=1;
                            } else if (mass_ratio_id == 5) {
                              match->match_uses.m_user=1;
                            }
                            if ((sin2w_exp_up != 0) || (cos2w_exp_up != 0)) {
                              match->match_uses.sin2w=1;
                            }
                            nle_state->phase1_matches_count++;
                            match++;
                          }  // if interesting term1
                        } // sanity check c2w term1
                      } // sanity check s2w term1

                      // test multiplier against term2 coefficient
                      if ((sin2w_exp_down == 1) || (sin2w_exp_down == 2) || (sin2w_exp_down == nle_state->term2.exp_inv) || (sin2w_exp_down == (nle_state->term2.exp_inv * 2))) {
                        if ((cos2w_exp_down == 1) || (cos2w_exp_down == 2) || (cos2w_exp_down == nle_state->term2.exp_inv) || (cos2w_exp_down == (nle_state->term2.exp_inv * 2))) {
                        multiplier=nle_state->term2.coefficient * term2_mass_ratio_factor * outfactor_sin2w * outfactor_cos2w * outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term2.exp_inv)];
                          if (interesting(phase1_filter, max_int, filter_int, multiplier)) {
                            nle_state->terms_matched[1]=nle_state->term2.exp_inv;
                            match->exp_inv=nle_state->term2.exp_inv;
                            match->smrfactor_mass=mass_ratio_id;
                            match->infactor_rational_up=infactors->infactor_rational_up;
                            match->infactor_rational_down=infactors->infactor_rational_down;
                            match->infactor_2_exp_up=infactors->infactor_2_exp_up;
                            match->infactor_2_exp_down=infactors->infactor_2_exp_down;
                            match->infactor_alpha_exp_up=infactors->infactor_alpha_exp_up;
                            match->infactor_alpha_exp_down=infactors->infactor_alpha_exp_down;
                            match->infactor_pi_exp_up=infactors->infactor_pi_exp_up;
                            match->infactor_pi_exp_down=infactors->infactor_pi_exp_down;
                            match->infactor_nss=infactors->infactor_nss;
                            match->infactor_nbv=infactors->infactor_nbv;
                            match->infactor_user_exp_up=infactors->infactor_user_exp_up;
                            match->infactor_user_exp_down=infactors->infactor_user_exp_down;
                            match->outfactor_rational_up=outfactors->outfactor_rational_up;
                            match->outfactor_rational_down=outfactors->outfactor_rational_down;
                            match->outfactor_2_exp_up=outfactors->outfactor_2_exp_up;
                            match->outfactor_2_exp_down=outfactors->outfactor_2_exp_down;
                            match->outfactor_alpha_exp_up=outfactors->outfactor_alpha_exp_up;
                            match->outfactor_alpha_exp_down=outfactors->outfactor_alpha_exp_down;
                            match->outfactor_pi_exp_up=outfactors->outfactor_pi_exp_up;
                            match->outfactor_pi_exp_down=outfactors->outfactor_pi_exp_down;
                            match->outfactor_sin2w_exp_up=sin2w_exp_up;
                            match->outfactor_sin2w_exp_down=sin2w_exp_down;
                            match->outfactor_cos2w_exp_up=cos2w_exp_up;
                            match->outfactor_cos2w_exp_down=cos2w_exp_down;
                            match->outfactor_user_exp_up=outfactors->outfactor_user_exp_up;
                            match->outfactor_user_exp_down=outfactors->outfactor_user_exp_down;
                            match->static_multiplier=outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term2.exp_inv)]; // excludes sin2w and cos2w that have uncertainty
                            match->match=multiplier;
                            match->match_complexity=outfactors->outfactor_complexity\
                                  + infactors->infactor_complexity\
                                  + (abs(sin2w_exp_up) + (sin2w_exp_down-1))\
                                  + (abs(cos2w_exp_up) + (cos2w_exp_down-1));
                            initUses(&match->match_uses);
                            addUses(&match->match_uses, &outfactors->outfactor_uses);
                            addUses(&match->match_uses, &infactors->infactor_uses);
                            if (mass_ratio_id == 0) {
                              match->match_uses.G=1;
                            } else if (mass_ratio_id == 1) {
                              match->match_uses.v=1;
                            } else if (mass_ratio_id == 2) {
                              match->match_uses.mz=1;
                            } else if (mass_ratio_id == 3) {
                              match->match_uses.mw=1;
                            } else if (mass_ratio_id == 4) {
                              match->match_uses.mh0=1;
                            } else if (mass_ratio_id == 5) {
                              match->match_uses.m_user=1;
                            }
                            if ((sin2w_exp_up != 0) || (cos2w_exp_up != 0)) {
                              match->match_uses.sin2w=1;
                            }
                            nle_state->phase1_matches_count++;
                            match++;
                          }  // if interesting term2
                        } // sanity check c2w term2
                      } // sanity check s2w term2

                      // test multiplier against term3 coefficient
                      if ((sin2w_exp_down == 1) || (sin2w_exp_down == 2) || (sin2w_exp_down == nle_state->term3.exp_inv) || (sin2w_exp_down == (nle_state->term3.exp_inv * 2))) {
                        if ((cos2w_exp_down == 1) || (cos2w_exp_down == 2) || (cos2w_exp_down == nle_state->term3.exp_inv) || (cos2w_exp_down == (nle_state->term3.exp_inv * 2))) {
                        multiplier=nle_state->term3.coefficient * term3_mass_ratio_factor * outfactor_sin2w * outfactor_cos2w * outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term3.exp_inv)];
                          if (interesting(phase1_filter, max_int, filter_int, multiplier)) {
                            nle_state->terms_matched[2]=nle_state->term3.exp_inv;
                            match->exp_inv=nle_state->term3.exp_inv;
                            match->smrfactor_mass=mass_ratio_id;
                            match->infactor_rational_up=infactors->infactor_rational_up;
                            match->infactor_rational_down=infactors->infactor_rational_down;
                            match->infactor_2_exp_up=infactors->infactor_2_exp_up;
                            match->infactor_2_exp_down=infactors->infactor_2_exp_down;
                            match->infactor_alpha_exp_up=infactors->infactor_alpha_exp_up;
                            match->infactor_alpha_exp_down=infactors->infactor_alpha_exp_down;
                            match->infactor_pi_exp_up=infactors->infactor_pi_exp_up;
                            match->infactor_pi_exp_down=infactors->infactor_pi_exp_down;
                            match->infactor_nss=infactors->infactor_nss;
                            match->infactor_nbv=infactors->infactor_nbv;
                            match->infactor_user_exp_up=infactors->infactor_user_exp_up;
                            match->infactor_user_exp_down=infactors->infactor_user_exp_down;
                            match->outfactor_rational_up=outfactors->outfactor_rational_up;
                            match->outfactor_rational_down=outfactors->outfactor_rational_down;
                            match->outfactor_2_exp_up=outfactors->outfactor_2_exp_up;
                            match->outfactor_2_exp_down=outfactors->outfactor_2_exp_down;
                            match->outfactor_alpha_exp_up=outfactors->outfactor_alpha_exp_up;
                            match->outfactor_alpha_exp_down=outfactors->outfactor_alpha_exp_down;
                            match->outfactor_pi_exp_up=outfactors->outfactor_pi_exp_up;
                            match->outfactor_pi_exp_down=outfactors->outfactor_pi_exp_down;
                            match->outfactor_sin2w_exp_up=sin2w_exp_up;
                            match->outfactor_sin2w_exp_down=sin2w_exp_down;
                            match->outfactor_cos2w_exp_up=cos2w_exp_up;
                            match->outfactor_cos2w_exp_down=cos2w_exp_down;
                            match->outfactor_user_exp_up=outfactors->outfactor_user_exp_up;
                            match->outfactor_user_exp_down=outfactors->outfactor_user_exp_down;
                            match->static_multiplier=outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term3.exp_inv)]; // excludes sin2w and cos2w that have uncertainty
                            match->match=multiplier;
                            match->match_complexity=outfactors->outfactor_complexity\
                                  + infactors->infactor_complexity\
                                  + (abs(sin2w_exp_up) + (sin2w_exp_down-1))\
                                  + (abs(cos2w_exp_up) + (cos2w_exp_down-1));
                            initUses(&match->match_uses);
                            addUses(&match->match_uses, &outfactors->outfactor_uses);
                            addUses(&match->match_uses, &infactors->infactor_uses);
                            if (mass_ratio_id == 0) {
                              match->match_uses.G=1;
                            } else if (mass_ratio_id == 1) {
                              match->match_uses.v=1;
                            } else if (mass_ratio_id == 2) {
                              match->match_uses.mz=1;
                            } else if (mass_ratio_id == 3) {
                              match->match_uses.mw=1;
                            } else if (mass_ratio_id == 4) {
                              match->match_uses.mh0=1;
                            } else if (mass_ratio_id == 5) {
                              match->match_uses.m_user=1;
                            }
                            if ((sin2w_exp_up != 0) || (cos2w_exp_up != 0)) {
                              match->match_uses.sin2w=1;
                            }
                            nle_state->phase1_matches_count++;
                            match++;
                          }  // if interesting term3
                        } // sanity check c2w term3
                      } // sanity check s2w term3

                      infactors++;
                    } // for infactor
                    outfactors++;
                  } // for outfactor
                } // gcd outfactor_cos2w
              } // cos2w_exp_down
            } // cos2w_exp_up
          } // gcd outfactor_sin2w
        } // sin2w_exp_down
      } // sin2w_exp_up
    } // end if mass_ratio_enabled
  } // for mass_ratio_id
  if (nle_config->status_enable ==1) {
    printf(".\n");
  }
}
