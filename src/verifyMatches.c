#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "nle-lepton.h"
#include "util.h"
#include "phase2.h"

//#define DEBUG_VERIFY

#ifdef DEBUG_VERIFY
#include "getFormulaStr.h"
#endif

void verifyMatches(nle_config_t *nle_config, nle_state_t *nle_state) {
  //  separate the match table into a separate list for each exponent, then test all unique combinations of coefficients
  int i,j;
  nle_phase1_match_t *phase1_match;
  nle_phase1_match_t *term1_match;
  nle_phase1_match_t *term2_match;
  nle_phase1_match_t *term3_match;
  nle_phase1_match_t *temp_match;
  int t1,t2,t3;
  int dupe;
  struct timespec start_time;
  struct timespec end_time;
  double elapsed_time;
  double precision=1.0E99;
  int tmpmatchup;
  int tmpmatchdown;
  int tmpmatchcomplexity;
  int upcomplexity, downcomplexity;
  long long tmphash;
  long combo_count;
  long combo;
  nle_input_use_t term1_uses;
  nle_input_use_t term2_uses;
  nle_input_use_t term3_uses;
  int complexity;
  int symmetry;
#ifdef DEBUG_VERIFY
  char term1_formula_str[288];
  char term2_formula_str[288];
  char term3_formula_str[288];
#endif

  // extract all term1 coefficients
  phase1_match=nle_state->phase1_matches_start;
  term1_match=nle_state->term1.matches_start;
  nle_state->term1.matches_count=0;
  for (i=0; i < nle_state->phase1_matches_count; i++) {
    if (phase1_match->term_id == 1) {
     // determine integer/rational match value
      if (phase1_match->match > 1.0) {
       tmpmatchup=(int)(phase1_match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / phase1_match->match) + 0.5);
      }
      tmphash=(long long)phase1_match->smrfactor_mass ^ ((long long)((((double)tmpmatchup / (double)tmpmatchdown) * (1.0 / phase1_match->static_multiplier) * 1.0E9) + 0.5));
      // ignore 1 on rationals
      if (tmpmatchup == 1) {
        upcomplexity=0;
      } else {
        upcomplexity=tmpmatchup;
      }
      if (tmpmatchdown == 1) {
        downcomplexity=0;
      } else {
        downcomplexity=tmpmatchdown;
      }
      tmpmatchcomplexity=(phase1_match->match_complexity + upcomplexity + downcomplexity);
#ifdef DEBUG_VERIFY
      printf("debug, term1, match_seq: %d, tmpmatchup: %d, tmpmatchdown: %d, tmpmatchcomplexity: %d, tmphash: %lld\n", i, tmpmatchup, tmpmatchdown, tmpmatchcomplexity, tmphash);
      getFormulaStr(nle_config, term1_formula_str, phase1_match);
      printf("debug, term1, phase1_match=%s\n", term1_formula_str);
      fflush(stdout);
#endif
      // search existing match table for dupes and see if we have lower complexity
      temp_match=nle_state->term1.matches_start;
      dupe=0;
      for (j=0; j< nle_state->term1.matches_count; j++) {
        if (tmphash == temp_match->match_hash) {
          if (tmpmatchcomplexity <= temp_match->match_complexity) {
#ifdef DEBUG_VERIFY
            printf("debug, term1, replace, tmpmatchcomplexity: %d, temp_match->match_complexity: %d, tmphash: %lld\n", tmpmatchcomplexity, temp_match->match_complexity, tmphash);
            getFormulaStr(nle_config, term1_formula_str, temp_match);
            printf("debug, term1, replace, old:=%s\n", term1_formula_str);
#endif
            // replace
            temp_match->term_id=1;
            temp_match->exp_inv=phase1_match->exp_inv;
            temp_match->smrfactor_mass=phase1_match->smrfactor_mass;
            temp_match->infactor_rational_up=phase1_match->infactor_rational_up;
            temp_match->infactor_rational_down=phase1_match->infactor_rational_down;
            temp_match->infactor_2_exp_up=phase1_match->infactor_2_exp_up;
            temp_match->infactor_2_exp_down=phase1_match->infactor_2_exp_down;
            temp_match->infactor_alpha_exp_up=phase1_match->infactor_alpha_exp_up;
            temp_match->infactor_alpha_exp_down=phase1_match->infactor_alpha_exp_down;
            temp_match->infactor_pi_exp_up=phase1_match->infactor_pi_exp_up;
            temp_match->infactor_pi_exp_down=phase1_match->infactor_pi_exp_down;
            temp_match->infactor_nss=phase1_match->infactor_nss;
            temp_match->infactor_nbv=phase1_match->infactor_nbv;
            temp_match->infactor_user_exp_up=phase1_match->infactor_user_exp_up;
            temp_match->infactor_user_exp_down=phase1_match->infactor_user_exp_down;
            temp_match->outfactor_rational_up=phase1_match->outfactor_rational_up;
            temp_match->outfactor_rational_down=phase1_match->outfactor_rational_down;
            temp_match->outfactor_2_exp_up=phase1_match->outfactor_2_exp_up;
            temp_match->outfactor_2_exp_down=phase1_match->outfactor_2_exp_down;
            temp_match->outfactor_alpha_exp_up=phase1_match->outfactor_alpha_exp_up;
            temp_match->outfactor_alpha_exp_down=phase1_match->outfactor_alpha_exp_down;
            temp_match->outfactor_pi_exp_up=phase1_match->outfactor_pi_exp_up;
            temp_match->outfactor_pi_exp_down=phase1_match->outfactor_pi_exp_down;
            temp_match->outfactor_sin2w_exp_up=phase1_match->outfactor_sin2w_exp_up;
            temp_match->outfactor_sin2w_exp_down=phase1_match->outfactor_sin2w_exp_down;
            temp_match->outfactor_cos2w_exp_up=phase1_match->outfactor_cos2w_exp_up;
            temp_match->outfactor_cos2w_exp_down=phase1_match->outfactor_cos2w_exp_down;
            temp_match->outfactor_user1_exp_up=phase1_match->outfactor_user1_exp_up;
            temp_match->outfactor_user1_exp_down=phase1_match->outfactor_user1_exp_down;
            temp_match->outfactor_user2_exp_up=phase1_match->outfactor_user2_exp_up;
            temp_match->outfactor_user2_exp_down=phase1_match->outfactor_user2_exp_down;
            temp_match->outfactor_user3_exp_up=phase1_match->outfactor_user3_exp_up;
            temp_match->outfactor_user3_exp_down=phase1_match->outfactor_user3_exp_down;
            temp_match->static_multiplier=phase1_match->static_multiplier;
            temp_match->match=phase1_match->match;
            temp_match->match_up=tmpmatchup;
            temp_match->match_down=tmpmatchdown;
            temp_match->match_complexity=tmpmatchcomplexity;
            temp_match->match_hash=tmphash;
            initUses(&temp_match->match_uses);
            addUses(&temp_match->match_uses, &phase1_match->match_uses);
#ifdef DEBUG_VERIFY
            getFormulaStr(nle_config, term1_formula_str, temp_match);
            printf("debug, term1, replace, new:=%s\n", term1_formula_str);
            fflush(stdout);
#endif
          }
          dupe=1;
          break;
        }
        temp_match++;
      } // end for j
      if (dupe == 0) {
        // insert
        term1_match->term_id=1;
        term1_match->exp_inv=phase1_match->exp_inv;
        term1_match->smrfactor_mass=phase1_match->smrfactor_mass;
        term1_match->infactor_rational_up=phase1_match->infactor_rational_up;
        term1_match->infactor_rational_down=phase1_match->infactor_rational_down;
        term1_match->infactor_2_exp_up=phase1_match->infactor_2_exp_up;
        term1_match->infactor_2_exp_down=phase1_match->infactor_2_exp_down;
        term1_match->infactor_alpha_exp_up=phase1_match->infactor_alpha_exp_up;
        term1_match->infactor_alpha_exp_down=phase1_match->infactor_alpha_exp_down;
        term1_match->infactor_pi_exp_up=phase1_match->infactor_pi_exp_up;
        term1_match->infactor_pi_exp_down=phase1_match->infactor_pi_exp_down;
        term1_match->infactor_nss=phase1_match->infactor_nss;
        term1_match->infactor_nbv=phase1_match->infactor_nbv;
        term1_match->infactor_user_exp_up=phase1_match->infactor_user_exp_up;
        term1_match->infactor_user_exp_down=phase1_match->infactor_user_exp_down;
        term1_match->outfactor_rational_up=phase1_match->outfactor_rational_up;
        term1_match->outfactor_rational_down=phase1_match->outfactor_rational_down;
        term1_match->outfactor_2_exp_up=phase1_match->outfactor_2_exp_up;
        term1_match->outfactor_2_exp_down=phase1_match->outfactor_2_exp_down;
        term1_match->outfactor_alpha_exp_up=phase1_match->outfactor_alpha_exp_up;
        term1_match->outfactor_alpha_exp_down=phase1_match->outfactor_alpha_exp_down;
        term1_match->outfactor_pi_exp_up=phase1_match->outfactor_pi_exp_up;
        term1_match->outfactor_pi_exp_down=phase1_match->outfactor_pi_exp_down;
        term1_match->outfactor_sin2w_exp_up=phase1_match->outfactor_sin2w_exp_up;
        term1_match->outfactor_sin2w_exp_down=phase1_match->outfactor_sin2w_exp_down;
        term1_match->outfactor_cos2w_exp_up=phase1_match->outfactor_cos2w_exp_up;
        term1_match->outfactor_cos2w_exp_down=phase1_match->outfactor_cos2w_exp_down;
        term1_match->outfactor_user1_exp_up=phase1_match->outfactor_user1_exp_up;
        term1_match->outfactor_user1_exp_down=phase1_match->outfactor_user1_exp_down;
        term1_match->outfactor_user2_exp_up=phase1_match->outfactor_user2_exp_up;
        term1_match->outfactor_user2_exp_down=phase1_match->outfactor_user2_exp_down;
        term1_match->outfactor_user3_exp_up=phase1_match->outfactor_user3_exp_up;
        term1_match->outfactor_user3_exp_down=phase1_match->outfactor_user3_exp_down;
        term1_match->static_multiplier=phase1_match->static_multiplier;
        term1_match->match=phase1_match->match;
        term1_match->match_up=tmpmatchup;
        term1_match->match_down=tmpmatchdown;
        term1_match->match_complexity=tmpmatchcomplexity;
        term1_match->match_hash=tmphash;
        initUses(&term1_match->match_uses);
        addUses(&term1_match->match_uses, &phase1_match->match_uses);
#ifdef DEBUG_VERIFY
        getFormulaStr(nle_config, term1_formula_str, term1_match);
        printf("debug, term1, addnew term1=%s\n", term1_formula_str);
       fflush(stdout);
#endif
        nle_state->term1.matches_count++;
        term1_match++;
      } // end if not dupe
    }  // end if invexp
    phase1_match++;
  } // end for i

  // extract all term2 term coefficients
  phase1_match=nle_state->phase1_matches_start;
  term2_match=nle_state->term2.matches_start;
  nle_state->term2.matches_count=0;
  for (i=0; i < nle_state->phase1_matches_count; i++) {
    if (phase1_match->term_id == 2) {
     // determine integer/rational match value
      if (phase1_match->match > 1.0) {
       tmpmatchup=(int)(phase1_match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / phase1_match->match) + 0.5);
      }
      tmphash=(long long)phase1_match->smrfactor_mass ^ ((long long)((((double)tmpmatchup / (double)tmpmatchdown) * (1.0 / phase1_match->static_multiplier) * 1.0E9) + 0.5));
      // ignore 1 on rationals
      if (tmpmatchup == 1) {
        upcomplexity=0;                           
      } else {                                  
        upcomplexity=tmpmatchup;
      }                                         
      if (tmpmatchdown == 1) {
        downcomplexity=0;                         
      } else {                                  
        downcomplexity=tmpmatchdown;
      }                                         
      tmpmatchcomplexity=(phase1_match->match_complexity + upcomplexity + downcomplexity);
#ifdef DEBUG_VERIFY
      printf("debug, term2, match_seq: %d, tmpmatchup: %d, tmpmatchdown: %d, tmpmatchcomplexity: %d, tmphash: %lld\n", i, tmpmatchup, tmpmatchdown, tmpmatchcomplexity, tmphash);
      getFormulaStr(nle_config, term2_formula_str, phase1_match);
      printf("debug, term2, phase1_match=%s\n", term2_formula_str);
      fflush(stdout);
#endif
      // search existing match table for dupes and see if we have lower complexity
      temp_match=nle_state->term2.matches_start;
      dupe=0;
      for (j=0; j< nle_state->term2.matches_count; j++) {
        if (tmphash == temp_match->match_hash) {
          if (tmpmatchcomplexity <= temp_match->match_complexity) {
#ifdef DEBUG_VERIFY
            printf("debug, term2, replace, tmpmatchcomplexity: %d, temp_match->match_complexity: %d, tmpmash: %lld\n", tmpmatchcomplexity, temp_match->match_complexity, tmphash);
            getFormulaStr(nle_config, term2_formula_str, temp_match);
            printf("debug, term2, replace, old:=%s\n", term2_formula_str);
#endif
            // replace
            temp_match->term_id=2;
            temp_match->exp_inv=phase1_match->exp_inv;
            temp_match->smrfactor_mass=phase1_match->smrfactor_mass;
            temp_match->infactor_rational_up=phase1_match->infactor_rational_up;
            temp_match->infactor_rational_down=phase1_match->infactor_rational_down;
            temp_match->infactor_2_exp_up=phase1_match->infactor_2_exp_up;
            temp_match->infactor_2_exp_down=phase1_match->infactor_2_exp_down;
            temp_match->infactor_alpha_exp_up=phase1_match->infactor_alpha_exp_up;
            temp_match->infactor_alpha_exp_down=phase1_match->infactor_alpha_exp_down;
            temp_match->infactor_pi_exp_up=phase1_match->infactor_pi_exp_up;
            temp_match->infactor_pi_exp_down=phase1_match->infactor_pi_exp_down;
            temp_match->infactor_nss=phase1_match->infactor_nss;
            temp_match->infactor_nbv=phase1_match->infactor_nbv;
            temp_match->infactor_user_exp_up=phase1_match->infactor_user_exp_up;
            temp_match->infactor_user_exp_down=phase1_match->infactor_user_exp_down;
            temp_match->outfactor_rational_up=phase1_match->outfactor_rational_up;
            temp_match->outfactor_rational_down=phase1_match->outfactor_rational_down;
            temp_match->outfactor_2_exp_up=phase1_match->outfactor_2_exp_up;
            temp_match->outfactor_2_exp_down=phase1_match->outfactor_2_exp_down;
            temp_match->outfactor_alpha_exp_up=phase1_match->outfactor_alpha_exp_up;
            temp_match->outfactor_alpha_exp_down=phase1_match->outfactor_alpha_exp_down;
            temp_match->outfactor_pi_exp_up=phase1_match->outfactor_pi_exp_up;
            temp_match->outfactor_pi_exp_down=phase1_match->outfactor_pi_exp_down;
            temp_match->outfactor_sin2w_exp_up=phase1_match->outfactor_sin2w_exp_up;
            temp_match->outfactor_sin2w_exp_down=phase1_match->outfactor_sin2w_exp_down;
            temp_match->outfactor_cos2w_exp_up=phase1_match->outfactor_cos2w_exp_up;
            temp_match->outfactor_cos2w_exp_down=phase1_match->outfactor_cos2w_exp_down;
            temp_match->outfactor_user1_exp_up=phase1_match->outfactor_user1_exp_up;
            temp_match->outfactor_user1_exp_down=phase1_match->outfactor_user1_exp_down;
            temp_match->outfactor_user2_exp_up=phase1_match->outfactor_user2_exp_up;
            temp_match->outfactor_user2_exp_down=phase1_match->outfactor_user2_exp_down;
            temp_match->outfactor_user3_exp_up=phase1_match->outfactor_user3_exp_up;
            temp_match->outfactor_user3_exp_down=phase1_match->outfactor_user3_exp_down;
            temp_match->static_multiplier=phase1_match->static_multiplier;
            temp_match->match=phase1_match->match;
            temp_match->match_up=tmpmatchup;
            temp_match->match_down=tmpmatchdown;
            temp_match->match_complexity=tmpmatchcomplexity;
            temp_match->match_hash=tmphash;
            initUses(&temp_match->match_uses);
            addUses(&temp_match->match_uses, &phase1_match->match_uses);
#ifdef DEBUG_VERIFY
            getFormulaStr(nle_config, term2_formula_str, temp_match);
            printf("debug, term2, replace, new:=%s\n", term2_formula_str);
           fflush(stdout);
#endif
          }
          dupe=1;
          break;
        }
        temp_match++;
      } // end for j
      if (dupe ==0) {
        // insert
        term2_match->term_id=2;
        term2_match->exp_inv=phase1_match->exp_inv;
        term2_match->smrfactor_mass=phase1_match->smrfactor_mass;
        term2_match->infactor_rational_up=phase1_match->infactor_rational_up;
        term2_match->infactor_rational_down=phase1_match->infactor_rational_down;
        term2_match->infactor_2_exp_up=phase1_match->infactor_2_exp_up;
        term2_match->infactor_2_exp_down=phase1_match->infactor_2_exp_down;
        term2_match->infactor_alpha_exp_up=phase1_match->infactor_alpha_exp_up;
        term2_match->infactor_alpha_exp_down=phase1_match->infactor_alpha_exp_down;
        term2_match->infactor_pi_exp_up=phase1_match->infactor_pi_exp_up;
        term2_match->infactor_pi_exp_down=phase1_match->infactor_pi_exp_down;
        term2_match->infactor_nss=phase1_match->infactor_nss;
        term2_match->infactor_nbv=phase1_match->infactor_nbv;
        term2_match->infactor_user_exp_up=phase1_match->infactor_user_exp_up;
        term2_match->infactor_user_exp_down=phase1_match->infactor_user_exp_down;
        term2_match->outfactor_rational_up=phase1_match->outfactor_rational_up;
        term2_match->outfactor_rational_down=phase1_match->outfactor_rational_down;
        term2_match->outfactor_2_exp_up=phase1_match->outfactor_2_exp_up;
        term2_match->outfactor_2_exp_down=phase1_match->outfactor_2_exp_down;
        term2_match->outfactor_alpha_exp_up=phase1_match->outfactor_alpha_exp_up;
        term2_match->outfactor_alpha_exp_down=phase1_match->outfactor_alpha_exp_down;
        term2_match->outfactor_pi_exp_up=phase1_match->outfactor_pi_exp_up;
        term2_match->outfactor_pi_exp_down=phase1_match->outfactor_pi_exp_down;
        term2_match->outfactor_sin2w_exp_up=phase1_match->outfactor_sin2w_exp_up;
        term2_match->outfactor_sin2w_exp_down=phase1_match->outfactor_sin2w_exp_down;
        term2_match->outfactor_cos2w_exp_up=phase1_match->outfactor_cos2w_exp_up;
        term2_match->outfactor_cos2w_exp_down=phase1_match->outfactor_cos2w_exp_down;
        term2_match->outfactor_user1_exp_up=phase1_match->outfactor_user1_exp_up;
        term2_match->outfactor_user1_exp_down=phase1_match->outfactor_user1_exp_down;
        term2_match->outfactor_user2_exp_up=phase1_match->outfactor_user2_exp_up;
        term2_match->outfactor_user2_exp_down=phase1_match->outfactor_user2_exp_down;
        term2_match->outfactor_user3_exp_up=phase1_match->outfactor_user3_exp_up;
        term2_match->outfactor_user3_exp_down=phase1_match->outfactor_user3_exp_down;
        term2_match->static_multiplier=phase1_match->static_multiplier;
        term2_match->match=phase1_match->match;
        term2_match->match_up=tmpmatchup;
        term2_match->match_down=tmpmatchdown;
        term2_match->match_complexity=tmpmatchcomplexity;
        term2_match->match_hash=tmphash;
        initUses(&term2_match->match_uses);
        addUses(&term2_match->match_uses, &phase1_match->match_uses);
#ifdef DEBUG_VERIFY
        getFormulaStr(nle_config, term2_formula_str, term2_match);
        printf("debug, term2, addnew term2=%s\n", term2_formula_str);
        fflush(stdout);
#endif
        nle_state->term2.matches_count++;
        term2_match++;
      } // end if not dupe
    }  // end if invexp
    phase1_match++;
  } // end for i

  // extract all term3 coefficients
  phase1_match=nle_state->phase1_matches_start;
  term3_match=nle_state->term3.matches_start;
  nle_state->term3.matches_count=0;
  for (i=0; i < nle_state->phase1_matches_count; i++) {
    if (phase1_match->term_id == 3) {
     // determine integer/rational match value
      if (phase1_match->match > 1.0) {
       tmpmatchup=(int)(phase1_match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / phase1_match->match) + 0.5);
      }
      tmphash=(long long)phase1_match->smrfactor_mass ^ ((long long)((((double)tmpmatchup / (double)tmpmatchdown) * (1.0 / phase1_match->static_multiplier) * 1.0E9) + 0.5));
      // ignore 1 on rationals
      if (tmpmatchup == 1) {
        upcomplexity=0;                           
      } else {                                  
        upcomplexity=tmpmatchup;
      }                                         
      if (tmpmatchdown == 1) {
        downcomplexity=0;                         
      } else {                                  
        downcomplexity=tmpmatchdown;
      }                                         
      tmpmatchcomplexity=(phase1_match->match_complexity + upcomplexity + downcomplexity);
#ifdef DEBUG_VERIFY
      printf("debug, term3, match_seq: %d, tmpmatchup: %d, tmpmatchdown: %d, tmpmatchcomplexity: %d, tmphash: %lld\n", i, tmpmatchup, tmpmatchdown, tmpmatchcomplexity, tmphash);
      getFormulaStr(nle_config, term3_formula_str, phase1_match);
      printf("debug, term3, phase1_match=%s\n", term3_formula_str);
      fflush(stdout);
#endif
      // search existing match table for dupes and see if we have lower complexity
      temp_match=nle_state->term3.matches_start;
      dupe=0;
      for (j=0; j< nle_state->term3.matches_count; j++) {
        if (tmphash == temp_match->match_hash) {
          if (tmpmatchcomplexity <= temp_match->match_complexity) {
#ifdef DEBUG_VERIFY
            printf("debug, term3, replace, tmpmatchcomplexity: %d, temp_match->match_complexity: %d, tmphash: %lld\n", tmpmatchcomplexity, temp_match->match_complexity, tmphash);
            getFormulaStr(nle_config, term3_formula_str, temp_match);
            printf("debug, term3, replace, old:=%s\n", term3_formula_str);
#endif
            // replace
            temp_match->term_id=3;
            temp_match->exp_inv=phase1_match->exp_inv;
            temp_match->smrfactor_mass=phase1_match->smrfactor_mass;
            temp_match->infactor_rational_up=phase1_match->infactor_rational_up;
            temp_match->infactor_rational_down=phase1_match->infactor_rational_down;
            temp_match->infactor_2_exp_up=phase1_match->infactor_2_exp_up;
            temp_match->infactor_2_exp_down=phase1_match->infactor_2_exp_down;
            temp_match->infactor_alpha_exp_up=phase1_match->infactor_alpha_exp_up;
            temp_match->infactor_alpha_exp_down=phase1_match->infactor_alpha_exp_down;
            temp_match->infactor_pi_exp_up=phase1_match->infactor_pi_exp_up;
            temp_match->infactor_pi_exp_down=phase1_match->infactor_pi_exp_down;
            temp_match->infactor_nss=phase1_match->infactor_nss;
            temp_match->infactor_nbv=phase1_match->infactor_nbv;
            temp_match->infactor_user_exp_up=phase1_match->infactor_user_exp_up;
            temp_match->infactor_user_exp_down=phase1_match->infactor_user_exp_down;
            temp_match->outfactor_rational_up=phase1_match->outfactor_rational_up;
            temp_match->outfactor_rational_down=phase1_match->outfactor_rational_down;
            temp_match->outfactor_2_exp_up=phase1_match->outfactor_2_exp_up;
            temp_match->outfactor_2_exp_down=phase1_match->outfactor_2_exp_down;
            temp_match->outfactor_alpha_exp_up=phase1_match->outfactor_alpha_exp_up;
            temp_match->outfactor_alpha_exp_down=phase1_match->outfactor_alpha_exp_down;
            temp_match->outfactor_pi_exp_up=phase1_match->outfactor_pi_exp_up;
            temp_match->outfactor_pi_exp_down=phase1_match->outfactor_pi_exp_down;
            temp_match->outfactor_sin2w_exp_up=phase1_match->outfactor_sin2w_exp_up;
            temp_match->outfactor_sin2w_exp_down=phase1_match->outfactor_sin2w_exp_down;
            temp_match->outfactor_cos2w_exp_up=phase1_match->outfactor_cos2w_exp_up;
            temp_match->outfactor_cos2w_exp_down=phase1_match->outfactor_cos2w_exp_down;
            temp_match->outfactor_user1_exp_up=phase1_match->outfactor_user1_exp_up;
            temp_match->outfactor_user1_exp_down=phase1_match->outfactor_user1_exp_down;
            temp_match->outfactor_user2_exp_up=phase1_match->outfactor_user2_exp_up;
            temp_match->outfactor_user2_exp_down=phase1_match->outfactor_user2_exp_down;
            temp_match->outfactor_user3_exp_up=phase1_match->outfactor_user3_exp_up;
            temp_match->outfactor_user3_exp_down=phase1_match->outfactor_user3_exp_down;
            temp_match->static_multiplier=phase1_match->static_multiplier;
            temp_match->match=phase1_match->match;
            temp_match->match_up=tmpmatchup;
            temp_match->match_down=tmpmatchdown;
            temp_match->match_complexity=tmpmatchcomplexity;
            temp_match->match_hash=tmphash;
            initUses(&temp_match->match_uses);
            addUses(&temp_match->match_uses, &phase1_match->match_uses);
#ifdef DEBUG_VERIFY
            getFormulaStr(nle_config, term3_formula_str, temp_match);
            printf("debug, term3, replace, new:=%s\n", term3_formula_str);
            fflush(stdout);
#endif
          }
          dupe=1;
          break;
        }
        temp_match++;
      } // end for j
      if (dupe ==0) {
        // insert
        term3_match->term_id=3;
        term3_match->exp_inv=phase1_match->exp_inv;
        term3_match->smrfactor_mass=phase1_match->smrfactor_mass;
        term3_match->infactor_rational_up=phase1_match->infactor_rational_up;
        term3_match->infactor_rational_down=phase1_match->infactor_rational_down;
        term3_match->infactor_2_exp_up=phase1_match->infactor_2_exp_up;
        term3_match->infactor_2_exp_down=phase1_match->infactor_2_exp_down;
        term3_match->infactor_alpha_exp_up=phase1_match->infactor_alpha_exp_up;
        term3_match->infactor_alpha_exp_down=phase1_match->infactor_alpha_exp_down;
        term3_match->infactor_pi_exp_up=phase1_match->infactor_pi_exp_up;
        term3_match->infactor_pi_exp_down=phase1_match->infactor_pi_exp_down;
        term3_match->infactor_nss=phase1_match->infactor_nss;
        term3_match->infactor_nbv=phase1_match->infactor_nbv;
        term3_match->infactor_user_exp_up=phase1_match->infactor_user_exp_up;
        term3_match->infactor_user_exp_down=phase1_match->infactor_user_exp_down;
        term3_match->outfactor_rational_up=phase1_match->outfactor_rational_up;
        term3_match->outfactor_rational_down=phase1_match->outfactor_rational_down;
        term3_match->outfactor_2_exp_up=phase1_match->outfactor_2_exp_up;
        term3_match->outfactor_2_exp_down=phase1_match->outfactor_2_exp_down;
        term3_match->outfactor_alpha_exp_up=phase1_match->outfactor_alpha_exp_up;
        term3_match->outfactor_alpha_exp_down=phase1_match->outfactor_alpha_exp_down;
        term3_match->outfactor_pi_exp_up=phase1_match->outfactor_pi_exp_up;
        term3_match->outfactor_pi_exp_down=phase1_match->outfactor_pi_exp_down;
        term3_match->outfactor_sin2w_exp_up=phase1_match->outfactor_sin2w_exp_up;
        term3_match->outfactor_sin2w_exp_down=phase1_match->outfactor_sin2w_exp_down;
        term3_match->outfactor_cos2w_exp_up=phase1_match->outfactor_cos2w_exp_up;
        term3_match->outfactor_cos2w_exp_down=phase1_match->outfactor_cos2w_exp_down;
        term3_match->outfactor_user1_exp_up=phase1_match->outfactor_user1_exp_up;
        term3_match->outfactor_user1_exp_down=phase1_match->outfactor_user1_exp_down;
        term3_match->outfactor_user2_exp_up=phase1_match->outfactor_user2_exp_up;
        term3_match->outfactor_user2_exp_down=phase1_match->outfactor_user2_exp_down;
        term3_match->outfactor_user3_exp_up=phase1_match->outfactor_user3_exp_up;
        term3_match->outfactor_user3_exp_down=phase1_match->outfactor_user3_exp_down;
        term3_match->static_multiplier=phase1_match->static_multiplier;
        term3_match->match=phase1_match->match;
        term3_match->match_up=tmpmatchup;
        term3_match->match_down=tmpmatchdown;
        term3_match->match_complexity=tmpmatchcomplexity;
        term3_match->match_hash=tmphash;
        initUses(&term3_match->match_uses);
        addUses(&term3_match->match_uses, &phase1_match->match_uses);
#ifdef DEBUG_VERIFY
        getFormulaStr(nle_config, term3_formula_str, term3_match);
        printf("debug, term3, addnew term3=%s\n", term3_formula_str);
        fflush(stdout);
#endif
        nle_state->term3.matches_count++;
        term3_match++;
      } // end if not dupe
    }  // end if invexp
    phase1_match++;
  } // end for i

  // generate all combinations and send to phase2
  term1_match=nle_state->term1.matches_start;
  combo=0;
  combo_count=(nle_state->term1.matches_count * nle_state->term2.matches_count * nle_state->term3.matches_count);
  if (nle_config->status_enable == 1) {
    printf("status, Solving phase 2 formulas for masses, input sample: %d, exponents: %s,                 progress: total (0/%ld) term1 (0/%d) term2 (0/%d) term3 (0/%d)\n", nle_state->phase1_seq, nle_state->exponents_str, combo_count, nle_state->term1.matches_count, nle_state->term2.matches_count, nle_state->term3.matches_count);
    fflush(stdout);
  }
  for (t1=0; t1 < nle_state->term1.matches_count; t1++) {
#ifdef DEBUG_VERIFY
    nle_state->term1.current_match=term1_match;
    getFormulaStr(nle_config, term1_formula_str, nle_state->term1.current_match);
    printf("term1=%s\n", term1_formula_str);
    fflush(stdout);
#endif
    initUses(&term1_uses);
    addUses(&term1_uses, &term1_match->match_uses);
    term2_match=nle_state->term2.matches_start;
    for (t2=0; t2 < nle_state->term2.matches_count; t2++) {
#ifdef DEBUG_VERIFY
      nle_state->term2.current_match=term2_match;
      getFormulaStr(nle_config, term2_formula_str, nle_state->term2.current_match);
      printf("term2=%s\n", term2_formula_str);
      fflush(stdout);
#endif
      initUses(&term2_uses);
      addUses(&term2_uses, &term2_match->match_uses);
      term3_match=nle_state->term3.matches_start;
      for (t3=0; t3 < nle_state->term3.matches_count; t3++) {
        combo++;
#ifdef DEBUG_VERIFY
        nle_state->term3.current_match=term3_match;
        getFormulaStr(nle_config, term3_formula_str, nle_state->term3.current_match);
        printf("term3=%s\n", term3_formula_str);
        fflush(stdout);
#endif

        // calculate complexity score
        complexity=term1_match->match_complexity + term2_match->match_complexity + term3_match->match_complexity;

        // calculate symmetry score.   This measures how many factors are identical or inverse identical between the dfferent terms
        symmetry=0;
        if (nle_config->nle_mode == 2) {
          checkSymmetry2(&symmetry, (term1_match->match_up * term1_match->outfactor_rational_down), (term2_match->match_up * term2_match->outfactor_rational_down));
          checkSymmetry2(&symmetry, (term1_match->match_down * term1_match->outfactor_rational_up), (term2_match->match_down * term2_match->outfactor_rational_up));
          checkSymmetry2(&symmetry, (term1_match->outfactor_2_exp_up * term1_match->outfactor_2_exp_down), (term2_match->outfactor_2_exp_up * term2_match->outfactor_2_exp_down));
          checkSymmetry2(&symmetry, (term1_match->outfactor_pi_exp_up * term1_match->outfactor_pi_exp_down), (term2_match->outfactor_pi_exp_up * term2_match->outfactor_pi_exp_down));
          checkSymmetry2(&symmetry, (term1_match->outfactor_alpha_exp_up * term1_match->outfactor_alpha_exp_down), (term2_match->outfactor_alpha_exp_up * term2_match->outfactor_alpha_exp_down));
          checkSymmetry2(&symmetry, (term1_match->outfactor_sin2w_exp_up * term1_match->outfactor_sin2w_exp_down), (term2_match->outfactor_sin2w_exp_up * term2_match->outfactor_sin2w_exp_down));
          checkSymmetry2(&symmetry, (term1_match->outfactor_cos2w_exp_up * term1_match->outfactor_cos2w_exp_down), (term2_match->outfactor_cos2w_exp_up * term2_match->outfactor_cos2w_exp_down));
          checkSymmetry2(&symmetry, (term1_match->outfactor_user1_exp_up * term1_match->outfactor_user1_exp_down), (term2_match->outfactor_user1_exp_up * term2_match->outfactor_user1_exp_down));
          checkSymmetry2(&symmetry, (term1_match->outfactor_user2_exp_up * term1_match->outfactor_user2_exp_down), (term2_match->outfactor_user2_exp_up * term2_match->outfactor_user2_exp_down));
          checkSymmetry2(&symmetry, (term1_match->outfactor_user3_exp_up * term1_match->outfactor_user3_exp_down), (term2_match->outfactor_user3_exp_up * term2_match->outfactor_user3_exp_down));
          checkSymmetry2(&symmetry, term1_match->infactor_rational_up, term2_match->infactor_rational_up);
          checkSymmetry2(&symmetry, term1_match->infactor_rational_down, term2_match->infactor_rational_down);
          checkSymmetry2(&symmetry, term1_match->infactor_nbv, term2_match->infactor_nbv);
          checkSymmetry2(&symmetry, term1_match->infactor_nss, term2_match->infactor_nss);
          checkSymmetry2(&symmetry, (term1_match->infactor_2_exp_up * term1_match->infactor_2_exp_down), (term2_match->infactor_2_exp_up * term2_match->infactor_2_exp_down));
          checkSymmetry2(&symmetry, (term1_match->infactor_pi_exp_up * term1_match->infactor_pi_exp_down), (term2_match->infactor_pi_exp_up * term2_match->infactor_pi_exp_down));
          checkSymmetry2(&symmetry, (term1_match->infactor_alpha_exp_up * term1_match->infactor_alpha_exp_down), (term2_match->infactor_alpha_exp_up * term2_match->infactor_alpha_exp_down));
          checkSymmetry2(&symmetry, (term1_match->infactor_user_exp_up * term1_match->infactor_user_exp_down), (term2_match->infactor_user_exp_up * term2_match->infactor_user_exp_down));
        } else if (nle_config->nle_mode == 3) {
          checkSymmetry3(&symmetry, (term1_match->match_up * term1_match->outfactor_rational_down), (term2_match->match_up * term2_match->outfactor_rational_down), (term3_match->match_up * term3_match->outfactor_rational_down));
          checkSymmetry3(&symmetry, (term1_match->match_down * term1_match->outfactor_rational_up), (term2_match->match_down * term2_match->outfactor_rational_up), (term3_match->match_down * term3_match->outfactor_rational_up));
          checkSymmetry3(&symmetry, (term1_match->outfactor_2_exp_up * term1_match->outfactor_2_exp_down), (term2_match->outfactor_2_exp_up * term2_match->outfactor_2_exp_down), (term3_match->outfactor_2_exp_up * term3_match->outfactor_2_exp_down));
          checkSymmetry3(&symmetry, (term1_match->outfactor_pi_exp_up * term1_match->outfactor_pi_exp_down), (term2_match->outfactor_pi_exp_up * term2_match->outfactor_pi_exp_down), (term3_match->outfactor_pi_exp_up * term3_match->outfactor_pi_exp_down));
          checkSymmetry3(&symmetry, (term1_match->outfactor_alpha_exp_up * term1_match->outfactor_alpha_exp_down), (term2_match->outfactor_alpha_exp_up * term2_match->outfactor_alpha_exp_down), (term3_match->outfactor_alpha_exp_up * term3_match->outfactor_alpha_exp_down));
          checkSymmetry3(&symmetry, (term1_match->outfactor_sin2w_exp_up * term1_match->outfactor_sin2w_exp_down), (term2_match->outfactor_sin2w_exp_up * term2_match->outfactor_sin2w_exp_down), (term3_match->outfactor_sin2w_exp_up * term3_match->outfactor_sin2w_exp_down));
          checkSymmetry3(&symmetry, (term1_match->outfactor_cos2w_exp_up * term1_match->outfactor_cos2w_exp_down), (term2_match->outfactor_cos2w_exp_up * term2_match->outfactor_cos2w_exp_down), (term3_match->outfactor_cos2w_exp_up * term3_match->outfactor_cos2w_exp_down));
          checkSymmetry3(&symmetry, (term1_match->outfactor_user1_exp_up * term1_match->outfactor_user1_exp_down), (term2_match->outfactor_user1_exp_up * term2_match->outfactor_user1_exp_down), (term3_match->outfactor_user1_exp_up * term3_match->outfactor_user1_exp_down));
          checkSymmetry3(&symmetry, (term1_match->outfactor_user2_exp_up * term1_match->outfactor_user2_exp_down), (term2_match->outfactor_user2_exp_up * term2_match->outfactor_user2_exp_down), (term3_match->outfactor_user2_exp_up * term3_match->outfactor_user2_exp_down));
          checkSymmetry3(&symmetry, (term1_match->outfactor_user3_exp_up * term1_match->outfactor_user3_exp_down), (term2_match->outfactor_user3_exp_up * term2_match->outfactor_user3_exp_down), (term3_match->outfactor_user3_exp_up * term3_match->outfactor_user3_exp_down));
          checkSymmetry3(&symmetry, term1_match->infactor_rational_up, term2_match->infactor_rational_up, term3_match->infactor_rational_up);
          checkSymmetry3(&symmetry, term1_match->infactor_rational_down, term2_match->infactor_rational_down, term3_match->infactor_rational_down);
          checkSymmetry3(&symmetry, term1_match->infactor_nbv, term2_match->infactor_nbv, term3_match->infactor_nbv);
          checkSymmetry3(&symmetry, term1_match->infactor_nss, term2_match->infactor_nss, term3_match->infactor_nss);
          checkSymmetry3(&symmetry, (term1_match->infactor_2_exp_up * term1_match->infactor_2_exp_down), (term2_match->infactor_2_exp_up * term2_match->infactor_2_exp_down), (term3_match->infactor_2_exp_up * term3_match->infactor_2_exp_down));
          checkSymmetry3(&symmetry, (term1_match->infactor_pi_exp_up * term1_match->infactor_pi_exp_down), (term2_match->infactor_pi_exp_up * term2_match->infactor_pi_exp_down), (term3_match->infactor_pi_exp_up * term3_match->infactor_pi_exp_down));
          checkSymmetry3(&symmetry, (term1_match->infactor_alpha_exp_up * term1_match->infactor_alpha_exp_down), (term2_match->infactor_alpha_exp_up * term2_match->infactor_alpha_exp_down), (term3_match->infactor_alpha_exp_up * term3_match->infactor_alpha_exp_down));
          checkSymmetry3(&symmetry, (term1_match->infactor_user_exp_up * term1_match->infactor_user_exp_down), (term2_match->infactor_user_exp_up * term2_match->infactor_user_exp_down), (term3_match->infactor_user_exp_up * term3_match->infactor_user_exp_down));
        }
        if ((symmetry >= nle_config->phase2_symmetry_min) && (complexity <= nle_config->phase2_complexity_max)) {
         if ((nle_config->phase2_check_nbv_nss == 0) || ((term1_match->infactor_nbv == term2_match->infactor_nbv) && ((nle_config->nle_mode == 2) || (term1_match->infactor_nbv == term3_match->infactor_nbv))\
                                                      && (term1_match->infactor_nss == term2_match->infactor_nss) && ((nle_config->nle_mode == 2) || (term1_match->infactor_nss == term3_match->infactor_nss)))) { // consistency check
           if ((nle_config->phase2_check_weak == 0) || ((term1_match->outfactor_sin2w_exp_up == term2_match->outfactor_sin2w_exp_up) && ((nle_config->nle_mode == 2) || (term1_match->outfactor_sin2w_exp_up == term3_match->outfactor_sin2w_exp_up))\
                                                     && (term1_match->outfactor_sin2w_exp_down == term2_match->outfactor_sin2w_exp_down) && ((nle_config->nle_mode == 2) || (term1_match->outfactor_sin2w_exp_down == term3_match->outfactor_sin2w_exp_down))\
                                                     && (term1_match->outfactor_cos2w_exp_up == term2_match->outfactor_cos2w_exp_up) && ((nle_config->nle_mode == 2) || (term1_match->outfactor_cos2w_exp_up == term3_match->outfactor_cos2w_exp_up))\
                                                     && (term1_match->outfactor_cos2w_exp_down == term2_match->outfactor_cos2w_exp_down) && ((nle_config->nle_mode == 2) || (term1_match->outfactor_cos2w_exp_down == term3_match->outfactor_cos2w_exp_down)))) {
              initUses(&term3_uses);
              addUses(&term3_uses, &term3_match->match_uses);
              initUses(&nle_state->all_uses);
              addUses(&nle_state->all_uses, &term1_uses);
              addUses(&nle_state->all_uses, &term2_uses);
              addUses(&nle_state->all_uses, &term3_uses);
              nle_state->term1.current_match=term1_match;
              nle_state->term2.current_match=term2_match;
              nle_state->term3.current_match=term3_match;
              nle_state->current_symmetry=symmetry;
              clock_gettime(CLOCK_REALTIME, &start_time);

              // send to phase2 to verify formula
              precision=solveNLEforMasses(nle_config, nle_state);

              if (nle_config->status_enable ==1) {
                clock_gettime(CLOCK_REALTIME, &end_time);
                elapsed_time=((double)(end_time.tv_sec - 1500000000) + ((double)end_time.tv_nsec / 1.0E9)) - ((double)(start_time.tv_sec - 1500000000) + ((double)start_time.tv_nsec) / 1.0E9);
                if (precision < 1.0E30) {
                  printf("status, Solved  phase 2 formula  for masses, input sample: %d, exponents: %s, mass mode: %d%d%d, progress: total (%ld/%ld) term1 (%d/%d) term2 (%d/%d) term3 (%d/%d), precision: %.3e, (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, term1_match->smrfactor_mass, term2_match->smrfactor_mass, term3_match->smrfactor_mass, combo, combo_count, t1+1, nle_state->term1.matches_count, t2+1, nle_state->term2.matches_count, t3+1, nle_state->term3.matches_count, precision, elapsed_time);
                  fflush(stdout);
                } else {
                  printf("status, Failed to solve phase 2 formula  for masses, input sample: %d, exponents: %s, mass mode: %d%d%d, progress: total (%ld/%ld) term1 (%d/%d) term2 (%d/%d) term3 (%d/%d), precision: %.3e, (%6.4fs)\n", nle_state->phase1_seq, nle_state->exponents_str, term1_match->smrfactor_mass, term2_match->smrfactor_mass, term3_match->smrfactor_mass, combo, combo_count, t1+1, nle_state->term1.matches_count, t2+1, nle_state->term2.matches_count, t3+1, nle_state->term3.matches_count, precision, elapsed_time);
                  fflush(stdout);
                } // end precision
              } // end status_enable
            } // end weak consistency check
          } // end nbv/nss consistency check
        } // end if symmetry and complexity check
        term3_match++;
      } // for t3
      term2_match++;
    } // for t2
    term1_match++;
  } // for t1
}
