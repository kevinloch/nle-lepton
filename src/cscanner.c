#include <stdio.h>
#include <stdlib.h> // abs
#include <math.h>   // pow
#include "nle-lepton.h"
#include "util.h"
#include "initDynamicfactorArray.h"

void cscanner(nle_config_t *nle_config, nle_state_t *nle_state) {
  //  Each coefficient is multiplied by various numbers and the resulting value is tested to see if it is
  // close to an interesting integer or simple rational number.  The results are stored in the match table.
  int infactor, outfactor, dynamicfactor;
  double multiplier;
  nle_infactor_precomputed_t *infactors;
  nle_outfactor_precomputed_t *outfactors;
  nle_dynamicfactor_precomputed_t *dynamicfactors;
  nle_dynamicfactor_precomputed_t *dynamicfactors_precomputed_start;
  int dynamicfactors_precomputed_count;
  nle_phase1_match_t *match;
  int phase1_filter;
  int max_int;
  int filter_int;
  int mass_ratio_id;
  int mass_ratio_id_max;
  int mass_ratio_enabled;
  char mass_str[32];

  double mass_ratio=1.0;
  double term1_mass_ratio_factor;
  double term2_mass_ratio_factor;
  double term3_mass_ratio_factor;

  phase1_filter=nle_config->phase1_filter;
  max_int=nle_config->phase1_int_match_max;
  filter_int=nle_config->phase1_int_match_filter;

  // allocate memory for and initialize precomputed dynamicfactor array
  dynamicfactors_precomputed_start=(nle_dynamicfactor_precomputed_t *)malloc(1000000 * sizeof(nle_dynamicfactor_precomputed_t));
  dynamicfactors_precomputed_count=0;
  initDynamicfactorArray(nle_config, nle_state, dynamicfactors_precomputed_start, &dynamicfactors_precomputed_count);

  if (nle_config->phase1_status_enable == 1) {
    printf("status, Scanning for coefficient multipliers that match interesting integer or simple rational numbers and reference mass: ");
    fflush(stdout);
  }

  // initialize match pointer to beginning and reset count
  match=nle_state->phase1_matches_start;
  nle_state->phase1_matches_count=0;

  // sequence through mass ratio reference masses, or run once if (1-smr) enabled
  if (nle_config->smrfactor_1minus_enable == 1) {
    mass_ratio_id_max=0;
  } else {
    mass_ratio_id_max=5;
  }
  // here we substitute the M/v mass ratio used in phase 1 with the actual test mass ratio
  for (mass_ratio_id=0; mass_ratio_id <= mass_ratio_id_max; mass_ratio_id++) {
    mass_ratio_enabled=1;
    if (nle_config->smrfactor_1minus_enable == 0) {
      if (mass_ratio_id == 0) {
        if (nle_config->smrfactor_mass_mp_enable == 1) {
          mass_ratio=nle_state->input_sample_mp/nle_config->ref_v;
          if (nle_config->phase1_status_enable == 1) {
            printf(" M/mp");
            fflush(stdout);
          }
        } else {
          mass_ratio_enabled=0;
        }
      } else if (mass_ratio_id == 1) {
        if (nle_config->smrfactor_mass_v_enable == 1) {
          mass_ratio=1.0;
          if (nle_config->phase1_status_enable == 1) {
            printf(" M/v");              
            fflush(stdout);               
          }
        } else {
          mass_ratio_enabled=0;
        }
      } else if (mass_ratio_id == 2) {
        if (nle_config->smrfactor_mass_mz_enable == 1) {
          mass_ratio=nle_state->input_sample_mz/nle_config->ref_v;
          if (nle_config->phase1_status_enable == 1) {
            printf(" M/mz");             
            fflush(stdout);               
          }
        } else {
          mass_ratio_enabled=0;
        } 
      } else if (mass_ratio_id == 3) {
        if (nle_config->smrfactor_mass_mw_enable == 1) {
          mass_ratio=nle_state->input_sample_mw/nle_config->ref_v;
          if (nle_config->phase1_status_enable == 1) {
            printf(" M/mw");             
            fflush(stdout);               
          }
        } else {
          mass_ratio_enabled=0;
        }
      } else if (mass_ratio_id == 4) {
        if (nle_config->smrfactor_mass_mh0_enable == 1) {
          mass_ratio=nle_state->input_sample_mh0/nle_config->ref_v;
          if (nle_config->phase1_status_enable == 1) {
            printf(" M/mh0");             
            fflush(stdout);               
          }
        } else {
          mass_ratio_enabled=0;
        }
      } else if (mass_ratio_id == 5) {
        if (nle_config->smrfactor_mass_user_enable == 1) {
          mass_ratio=nle_state->input_sample_muser/nle_config->ref_v;
          if (nle_config->phase1_status_enable == 1) {
            printf(" M/m_user");
            fflush(stdout);
          }
        } else {
          mass_ratio_enabled=0;
        }
      }                             
    } else {
      // overwrite mass_ratio_id with smrfactor mass ratio
      mass_ratio_id=nle_state->term1.smrfactor_mass;
      if (nle_state->term1.smrfactor_mass == 0) {
        sprintf(mass_str, "mP   ");
      } else if (nle_state->term1.smrfactor_mass == 1) {
        sprintf(mass_str, "v    ");
      } else if (nle_state->term1.smrfactor_mass == 2) {
        sprintf(mass_str, "mz   ");
      } else if (nle_state->term1.smrfactor_mass == 3) {
        sprintf(mass_str, "mw   ");
      } else if (nle_state->term1.smrfactor_mass == 4) {
        sprintf(mass_str, "mh0  ");
      } else if (nle_state->term1.smrfactor_mass == 5) {
        sprintf(mass_str, "muser");
      }
      if (nle_config->phase1_status_enable == 1) {
        printf("%s", mass_str);
      }
    } // end if 1minus_enable

    if (mass_ratio_enabled == 1) {
      if (nle_state->term1.smrfactor_1minus) {
        term1_mass_ratio_factor=1.0;
        term2_mass_ratio_factor=1.0;
        term3_mass_ratio_factor=1.0;
      } else {
        term1_mass_ratio_factor=pow(mass_ratio, (1.0 / (double)nle_state->term1.exp_inv));
        term2_mass_ratio_factor=pow(mass_ratio, (1.0 / (double)nle_state->term2.exp_inv));
        term3_mass_ratio_factor=pow(mass_ratio, (1.0 / (double)nle_state->term3.exp_inv));
      }

      // pre-computed static outfactors
      dynamicfactors=dynamicfactors_precomputed_start;
      for (dynamicfactor=0; dynamicfactor < dynamicfactors_precomputed_count; dynamicfactor++) {

        // pre-computed static outfactors
        outfactors=nle_state->outfactors_precomputed_start;
        for (outfactor=0; outfactor < nle_state->outfactors_precomputed_count; outfactor++) {

          // pre-computed static infactors
          infactors=nle_state->infactors_precomputed_start;
          for (infactor=0; infactor < nle_state->infactors_precomputed_count; infactor++) { 

            // test multiplier against term1 coefficient
            if (((dynamicfactors->outfactor_sin2w_exp_down == 1) || (dynamicfactors->outfactor_sin2w_exp_down == 2) || (dynamicfactors->outfactor_sin2w_exp_down == nle_state->term1.exp_inv) || (dynamicfactors->outfactor_sin2w_exp_down == (nle_state->term1.exp_inv * 2)))\
             && ((dynamicfactors->outfactor_cos2w_exp_down == 1) || (dynamicfactors->outfactor_cos2w_exp_down == 2) || (dynamicfactors->outfactor_cos2w_exp_down == nle_state->term1.exp_inv) || (dynamicfactors->outfactor_cos2w_exp_down == (nle_state->term1.exp_inv * 2)))) {
              multiplier=nle_state->term1.coefficient * term1_mass_ratio_factor * dynamicfactors->dynamicfactor_multiplier * outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term1.exp_inv)];
              if (interesting(phase1_filter, max_int, filter_int, multiplier)) {
                nle_state->terms_matched[0]=nle_state->term1.exp_inv;
                match->term_id=1;
                match->exp_inv=nle_state->term1.exp_inv;
                if (nle_config->smrfactor_1minus_enable == 1) {
                  match->smrfactor_mass=nle_state->term1.smrfactor_mass;
                } else {
                  match->smrfactor_mass=mass_ratio_id;
                }
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
                match->outfactor_sin2w_exp_up=dynamicfactors->outfactor_sin2w_exp_up;
                match->outfactor_sin2w_exp_down=dynamicfactors->outfactor_sin2w_exp_down;
                match->outfactor_cos2w_exp_up=dynamicfactors->outfactor_cos2w_exp_up;
                match->outfactor_cos2w_exp_down=dynamicfactors->outfactor_cos2w_exp_down;
                match->outfactor_rmr_exp_up=dynamicfactors->outfactor_rmr_exp_up;
                match->outfactor_rmr_exp_down=dynamicfactors->outfactor_rmr_exp_down;
                match->outfactor_rmr_mass_id_up=dynamicfactors->outfactor_rmr_mass_id_up;
                match->outfactor_rmr_mass_id_down=dynamicfactors->outfactor_rmr_mass_id_down;
                match->outfactor_user1_exp_up=outfactors->outfactor_user1_exp_up;
                match->outfactor_user1_exp_down=outfactors->outfactor_user1_exp_down;
                match->outfactor_user2_exp_up=outfactors->outfactor_user2_exp_up;
                match->outfactor_user2_exp_down=outfactors->outfactor_user2_exp_down;
                match->outfactor_user3_exp_up=outfactors->outfactor_user3_exp_up;
                match->outfactor_user3_exp_down=outfactors->outfactor_user3_exp_down;
                match->static_multiplier=outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term1.exp_inv)]; // excludes sin2w and cos2w and rmr that have uncertainty
                match->match=multiplier;
                match->match_complexity=outfactors->outfactor_complexity\
                      + infactors->infactor_complexity\
                      + dynamicfactors->dynamicfactor_complexity;
                initUses(&match->match_uses);
                if (nle_config->smrfactor_1minus_enable == 1) {
                  addUses(&match->match_uses, &nle_state->term1.current_smrfactors->smrfactor_uses);
                }
                addUses(&match->match_uses, &outfactors->outfactor_uses);
                addUses(&match->match_uses, &infactors->infactor_uses);
                addUses(&match->match_uses, &dynamicfactors->dynamicfactor_uses);
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
                nle_state->phase1_matches_count++;
                match++;
              }  // if interesting term1
            } // sanity check s2w term1

            // test multiplier against term2 coefficient
            if (((dynamicfactors->outfactor_sin2w_exp_down == 1) || (dynamicfactors->outfactor_sin2w_exp_down == 2) || (dynamicfactors->outfactor_sin2w_exp_down == nle_state->term2.exp_inv) || (dynamicfactors->outfactor_sin2w_exp_down == (nle_state->term2.exp_inv * 2)))\
             && ((dynamicfactors->outfactor_cos2w_exp_down == 1) || (dynamicfactors->outfactor_cos2w_exp_down == 2) || (dynamicfactors->outfactor_cos2w_exp_down == nle_state->term2.exp_inv) || (dynamicfactors->outfactor_cos2w_exp_down == (nle_state->term2.exp_inv * 2)))) {
              multiplier=nle_state->term2.coefficient * term2_mass_ratio_factor * dynamicfactors->dynamicfactor_multiplier * outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term2.exp_inv)];
              if (interesting(phase1_filter, max_int, filter_int, multiplier)) {
                nle_state->terms_matched[1]=nle_state->term2.exp_inv;
                match->term_id=2;
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
                match->outfactor_sin2w_exp_up=dynamicfactors->outfactor_sin2w_exp_up;
                match->outfactor_sin2w_exp_down=dynamicfactors->outfactor_sin2w_exp_down;
                match->outfactor_cos2w_exp_up=dynamicfactors->outfactor_cos2w_exp_up;
                match->outfactor_cos2w_exp_down=dynamicfactors->outfactor_cos2w_exp_down;
                match->outfactor_rmr_exp_up=dynamicfactors->outfactor_rmr_exp_up;
                match->outfactor_rmr_exp_down=dynamicfactors->outfactor_rmr_exp_down;
                match->outfactor_rmr_mass_id_up=dynamicfactors->outfactor_rmr_mass_id_up;
                match->outfactor_rmr_mass_id_down=dynamicfactors->outfactor_rmr_mass_id_down;
                match->outfactor_user1_exp_up=outfactors->outfactor_user1_exp_up;
                match->outfactor_user1_exp_down=outfactors->outfactor_user1_exp_down;
                match->outfactor_user2_exp_up=outfactors->outfactor_user2_exp_up;
                match->outfactor_user2_exp_down=outfactors->outfactor_user2_exp_down;
                match->outfactor_user3_exp_up=outfactors->outfactor_user3_exp_up;
                match->outfactor_user3_exp_down=outfactors->outfactor_user3_exp_down;
                match->static_multiplier=outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term2.exp_inv)]; // excludes sin2w and cos2w and rmr that have uncertainty
                match->match=multiplier;
                match->match_complexity=outfactors->outfactor_complexity\
                      + infactors->infactor_complexity\
                      + dynamicfactors->dynamicfactor_complexity;
                initUses(&match->match_uses);
                if (nle_config->smrfactor_1minus_enable == 1) {
                  addUses(&match->match_uses, &nle_state->term2.current_smrfactors->smrfactor_uses);
                }
                addUses(&match->match_uses, &outfactors->outfactor_uses);
                addUses(&match->match_uses, &infactors->infactor_uses);
                addUses(&match->match_uses, &dynamicfactors->dynamicfactor_uses);
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
                nle_state->phase1_matches_count++;
                match++;
              }  // if interesting term2
            } // sanity check s2w term2

            // test multiplier against term3 coefficient
            if ((nle_config->nle_mode == 3) || ((nle_config->nle_mode == 2) && (infactors->infactor_nss == 0) && (infactors->infactor_nbv == 0)\
               && (infactors->infactor_2_exp_up == 0) && (infactors->infactor_alpha_exp_up == 0) && (infactors->infactor_pi_exp_up == 0)\
               && (infactors->infactor_user_exp_up == 0) && (outfactors->outfactor_2_exp_up == 0) && (outfactors->outfactor_alpha_exp_up == 0)\
               && (outfactors->outfactor_pi_exp_up == 0) && (dynamicfactors->outfactor_sin2w_exp_up == 0) && (dynamicfactors->outfactor_cos2w_exp_up == 0)\
               && (outfactors->outfactor_user1_exp_up == 0) && (outfactors->outfactor_user2_exp_up == 0) && (outfactors->outfactor_user3_exp_up == 0)\
               && (dynamicfactors->outfactor_rmr_exp_up == 0))) { // if nle-mode == 2 only test rationals, no other factors
              if (((dynamicfactors->outfactor_sin2w_exp_down == 1) || (dynamicfactors->outfactor_sin2w_exp_down == 2) || (dynamicfactors->outfactor_sin2w_exp_down == nle_state->term3.exp_inv) || (dynamicfactors->outfactor_sin2w_exp_down == (nle_state->term3.exp_inv * 2)))\
               && ((dynamicfactors->outfactor_cos2w_exp_down == 1) || (dynamicfactors->outfactor_cos2w_exp_down == 2) || (dynamicfactors->outfactor_cos2w_exp_down == nle_state->term3.exp_inv) || (dynamicfactors->outfactor_cos2w_exp_down == (nle_state->term3.exp_inv * 2)))) {
                multiplier=nle_state->term3.coefficient * term3_mass_ratio_factor * dynamicfactors->dynamicfactor_multiplier * outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term3.exp_inv)];
                if (interesting(phase1_filter, max_int, filter_int, multiplier)) {
                  nle_state->terms_matched[2]=nle_state->term3.exp_inv;
                  match->term_id=3;
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
                  match->outfactor_sin2w_exp_up=dynamicfactors->outfactor_sin2w_exp_up;
                  match->outfactor_sin2w_exp_down=dynamicfactors->outfactor_sin2w_exp_down;
                  match->outfactor_cos2w_exp_up=dynamicfactors->outfactor_cos2w_exp_up;
                  match->outfactor_cos2w_exp_down=dynamicfactors->outfactor_cos2w_exp_down;
                  match->outfactor_rmr_exp_up=dynamicfactors->outfactor_rmr_exp_up;
                  match->outfactor_rmr_exp_down=dynamicfactors->outfactor_rmr_exp_down;
                  match->outfactor_rmr_mass_id_up=dynamicfactors->outfactor_rmr_mass_id_up;
                  match->outfactor_rmr_mass_id_down=dynamicfactors->outfactor_rmr_mass_id_down;
                  match->outfactor_user1_exp_up=outfactors->outfactor_user1_exp_up;
                  match->outfactor_user1_exp_down=outfactors->outfactor_user1_exp_down;
                  match->outfactor_user2_exp_up=outfactors->outfactor_user2_exp_up;
                  match->outfactor_user2_exp_down=outfactors->outfactor_user2_exp_down;
                  match->outfactor_user3_exp_up=outfactors->outfactor_user3_exp_up;
                  match->outfactor_user3_exp_down=outfactors->outfactor_user3_exp_down;
                  match->static_multiplier=outfactors->outfactor_multiplier * infactors->infactor_multiplier[abs(nle_state->term3.exp_inv)]; // excludes sin2w and cos2w and rmr that have uncertainty
                  match->match=multiplier;
                  match->match_complexity=outfactors->outfactor_complexity\
                        + infactors->infactor_complexity\
                        + dynamicfactors->dynamicfactor_complexity;
                  initUses(&match->match_uses);
                  addUses(&match->match_uses, &outfactors->outfactor_uses);
                  addUses(&match->match_uses, &infactors->infactor_uses);
                  addUses(&match->match_uses, &dynamicfactors->dynamicfactor_uses);
                  if (mass_ratio_id == 0) {
                    match->match_uses.G=1;
                  } else if ((mass_ratio_id == 1) && (nle_config->nle_mode != 2)) {
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
                  nle_state->phase1_matches_count++;
                  match++;
                }  // if interesting term3
              } // sanity check s2w term3
            } // end term 3 mode 2 check

            infactors++;
          } // end for infactor
          outfactors++;
        } // end for outfactor
        dynamicfactors++;
      } // end for dynamicfactor
    } // end if mass_ratio_enabled
  } // for mass_ratio_id
  if (nle_config->phase1_status_enable == 1) {
    printf(".\n");
  }
  free(dynamicfactors_precomputed_start);
}
