#include "nle-lepton.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "util.h"
#include "getFormulaStr.h"
#include "selectOutputs.h"

//#define DEBUG20
//#define DEBUG21
//#define DEBUG22
//#define DEBUG23

long double solveNLEforMasses(nle_config_t *nle_config, nle_state_t *nle_state) {
  // solve polynomial-like non-linear equation for particle masses using the supplied coefficients, exponent and factors
  long long samples=0;
  int i;
  int alpha_seq, sm1_seq, sm2_seq, v_seq, G_seq, mz_seq, mw_seq, mh0_seq, sm3_seq, sin2w_seq, muser_seq;
  struct timespec start_time, start_time2;
#ifdef DEBUG20
  struct timespec end_time;
  double elapsed_time;
  int unknowns=0;
#endif
  int valid_result;
  long double results_window;
  char exec_str[512];
  char mass_str[32];
  char out_str_01[512];
  char out_str_02[512];
  char out_str_03[512];
  char out_str_04[512];
  char out_str_05[512];
  char out_str_06[512];
  char out_str_07[512];
  char out_str_08[512];
  char out_str_09[512];
  char out_str_10[512];
  char out_str_11[512];
  char out_str_12[512];
  char out_str_13[512];
  char out_str_14[512];
  char out_str_15[512];
  char out_str_16[512];
  char out_str_17[512];
  char out_str_18[512];
  char out_str_19[512];
  char out_str_20[512];
  char used_as_input[5];
  char used_as_output[5];
  char user1_out_str[32];
  char user2_out_str[32];
  char user3_out_str[32];
  char user_in_str[32];
  long long result_hash;
  int complexity;
  long double term1_exp;
  long double term2_exp;
  long double term3_exp;
  char term1_formula_str[288];
  char term2_formula_str[288];
  char term3_formula_str[288];
  char smrf_str[80];
  int symmetry;
  float combined_score;
 
  //  mc test vars
  long double r;
  long double sm1_test_term1=0.0l;
  long double sm1_test_term2=0.0l;
  long double sm1_test_term3=0.0l;
  long double sm1_test=0.0l;
  long double sm2_test_term1=0.0l;
  long double sm2_test_term2=0.0l;
  long double sm2_test_term3=0.0l;
  long double sm2_test=0.0l;
  long double sm3_test_term1=0.0l;
  long double sm3_test_term2=0.0l;
  long double sm3_test_term3=0.0l;
  long double sm3_test=0.0l;
  long double precision=0.0l;
  long double precision_last=0.0l;
  long double term1_coefficient, term2_coefficient, term3_coefficient;
  long double term1_static, term2_static, term3_static;
  long double outfactor_rmr_mass_up=1.0l;
  long double outfactor_rmr_mass_down=1.0l;
  long double term1_rmr, term2_rmr, term3_rmr;
  long double term1_sin2w, term2_sin2w, term3_sin2w;
  long double term1_cos2w, term2_cos2w, term3_cos2w;
  long double term1_smrfactor_mass=0.0l;
  long double term2_smrfactor_mass=0.0l;
  long double term3_smrfactor_mass=0.0l;
  long double smrf_sm1=0.0l;
  long double smrf_sm2=0.0l;
  long double smrf_sm3=0.0l;
  long double term1_mass_sm1=0.0l;
  long double term2_mass_sm1=0.0l;
  long double term3_mass_sm1=0.0l;
  long double term1_mass_sm2=0.0l;
  long double term2_mass_sm2=0.0l;
  long double term3_mass_sm2=0.0l;
  long double term1_mass_sm3=0.0l;
  long double term2_mass_sm3=0.0l;
  long double term3_mass_sm3=0.0l;
  long double mp=0.0l;
  long double worst_test;
  long double rangefactor;
  long double range_multiplier;
  long long stalled;
  int progress;

  // mc outputs
  long double alpha=0.0l;
  long double alpha_last=0.0l;
  long double alpha_center=0.0l;
  long double alpha_range=0.0l;
  long double alpha_range_new=0.0l;
  long double sm1=0.0l;
  long double sm1_last=0.0l;
  long double sm1_center=0.0l;
  long double sm1_range=0.0l;
  long double sm1_range_new=0.0l;
  long double sm2=0.0l;
  long double sm2_last=0.0l;
  long double sm2_center=0.0l;
  long double sm2_range=0.0l;
  long double sm2_range_new=0.0l;
  long double v=0.0l;
  long double v_last=0.0l;
  long double v_center=0.0l;
  long double v_range=0.0l;
  long double v_range_new=0.0l;
  long double sm3=0.0l;
  long double sm3_last=0.0l;
  long double sm3_center=0.0l;
  long double sm3_range=0.0l;
  long double sm3_range_new=0.0l;
  long double G=0.0l;
  long double G_last=0.0l;
  long double G_center=0.0l;
  long double G_range=0.0l;
  long double G_range_new=0.0l;
  long double mz=0.0l;
  long double mz_last=0.0l;
  long double mz_center=0.0l;
  long double mz_range=0.0l;
  long double mz_range_new=0.0l;
  long double mw=0.0l;
  long double mw_last=0.0l;
  long double mw_center=0.0l;
  long double mw_range=0.0l;
  long double mw_range_new=0.0l;
  long double mh0=0.0l;
  long double mh0_last=0.0l;
  long double mh0_center=0.0l;
  long double mh0_range=0.0l;
  long double mh0_range_new=0.0l;
  long double sin2w=0.0l;
  long double sin2w_last=0.0l;
  long double sin2w_center=0.0l;
  long double sin2w_range=0.0l;
  long double sin2w_range_new=0.0l;
  long double cos2w=0.0l;
  long double muser=0.0l;
  long double muser_last=0.0l;
  long double muser_center=0.0l;
  long double muser_range=0.0l;
  long double muser_range_new=0.0l;

  // for reporting
  long double alpha_out=0.0l;
  long double alpha_out_low=1.0E30l;
  long double alpha_out_high=-1.0E30l;
  long double alpha_out_c=0.0l;
  long double alpha_out_error=0.0l;
  long double alpha_out_relerror=0.0l;
  long double alpha_out_diff=0.0l;
  long double alpha_out_reldiff=0.0l;
  long double sm1_out=0.0l;
  long double sm1_out_low=1.0E30l;
  long double sm1_out_high=-1.0E30l;
  long double sm1_out_c=0.0l;
  long double sm1_out_error=0.0l;
  long double sm1_out_relerror=0.0l;
  long double sm1_out_diff=0.0l;
  long double sm1_out_reldiff=0.0l;
  long double sm2_out=0.0l;
  long double sm2_out_low=1.0E30l;
  long double sm2_out_high=-1.0E30l;
  long double sm2_out_c=0.0l;
  long double sm2_out_error=0.0l;
  long double sm2_out_relerror=0.0l;
  long double sm2_out_diff=0.0l;
  long double sm2_out_reldiff=0.0l;
  long double sm3_out=0.0l;
  long double sm3_out_low=1.0E30l;
  long double sm3_out_high=-1.0E30l;
  long double sm3_out_c=0.0l;
  long double sm3_out_error=0.0l;
  long double sm3_out_relerror=0.0l;
  long double sm3_out_diff=0.0l;
  long double sm3_out_reldiff=0.0l;
  long double v_out=0.0l;
  long double v_out_low=1.0E30l;
  long double v_out_high=-1.0E30l;
  long double v_out_c=0.0l;
  long double v_out_error=0.0l;
  long double v_out_relerror=0.0l;
  long double v_out_diff=0.0l;
  long double v_out_reldiff=0.0l;
  long double G_out=0.0l;
  long double G_out_low=1.0E30l;
  long double G_out_high=-1.0E30l;
  long double G_out_c=0.0l;
  long double G_out_error=0.0l;
  long double G_out_relerror=0.0l;
  long double G_out_diff=0.0l;
  long double G_out_reldiff=0.0l;
  long double mz_out=0.0l;
  long double mz_out_low=1.0E30l;
  long double mz_out_high=-1.0E30l;
  long double mz_out_c=0.0l;
  long double mz_out_error=0.0l;
  long double mz_out_relerror=0.0l;
  long double mz_out_diff=0.0l;
  long double mz_out_reldiff=0.0l;
  long double mw_out=0.0l;
  long double mw_out_low=1.0E30l;
  long double mw_out_high=-1.0E30l;
  long double mw_out_c=0.0l;
  long double mw_out_error=0.0l;
  long double mw_out_relerror=0.0l;
  long double mw_out_diff=0.0l;
  long double mw_out_reldiff=0.0l;
  long double mh0_out=0.0l;
  long double mh0_out_low=1.0E30l;
  long double mh0_out_high=-1.0E30l;
  long double mh0_out_c=0.0l;
  long double mh0_out_error=0.0l;
  long double mh0_out_relerror=0.0l;
  long double mh0_out_diff=0.0l;
  long double mh0_out_reldiff=0.0l;
  long double sin2w_out=0.0l;
  long double sin2w_out_low=1.0E30l;
  long double sin2w_out_high=-1.0E30l;
  long double sin2w_out_c=0.0l;
  long double sin2w_out_error=0.0l;
  long double sin2w_out_relerror=0.0l;
  long double sin2w_out_diff=0.0l;
  long double sin2w_out_reldiff=0.0l;
  long double muser_out=0.0l;
  long double muser_out_low=1.0E30l;
  long double muser_out_high=-1.0E30l;
  long double muser_out_c=0.0l;
  long double muser_out_error=0.0l;
  long double muser_out_relerror=0.0l;
  long double muser_out_diff=0.0l;
  long double muser_out_reldiff=0.0l;

  // tuneables
  long double precision_target;
  long double test_ratio;
  int ratio_grace_period;
  int stalled_limit;
  long double default_range_multiplier;
  long double stalled_range_multiplier;
  int process_checkpoint;
  long double stuck_precision;

  // these tunings affect speed and reliability, adjust with extreme care
  if (nle_config->smrfactor_1minus_enable == 1) {
    // 2-term mode with 1-smr
    precision_target=1.0E-15l;      // solve NLE to this level of precision
    test_ratio=25.0l;               // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=3;           // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;           // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    default_range_multiplier=3.0l;  // lowest practical range multiplier, fastest for most formulas
    stalled_range_multiplier=3.0l;  // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    process_checkpoint=1000000;     // check on processing every this many samples
    stuck_precision=1.0E+30l;       // if precision is not past this level by process_checkpoint, try resetting
  } else if (nle_config->nle_mode == 2) {
    // 2-term mode without 1-smr
    precision_target=1.0E-15l;      // solve NLE to this level of precision
    test_ratio=25.0l;               // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=3;           // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;           // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    default_range_multiplier=1.0l;  // lowest practical range multiplier, fastest for most formulas
    stalled_range_multiplier=2.0l;  // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    process_checkpoint=1000000;     // check on processing every this many samples
    stuck_precision=1.0E+30l;       // if precision is not past this level by process_checkpoint, try resetting
  } else {
    // 3-term mode
    precision_target=1.0E-15l;      // solve NLE to this level of precision
    test_ratio=25.0l;               // acceptable ratios of sm1_test/sm2_test/sm3_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
    ratio_grace_period=3;           // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
    stalled_limit=500000;           // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
    default_range_multiplier=5.0l;  // lowest practical range multiplier, fastest for most formulas
    stalled_range_multiplier=10.0l; // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
    process_checkpoint=1000000;     // check on processing every this many samples
    stuck_precision=1.0E-2l;        // if precision is not past this level by process_checkpoint, try resetting
  }

  clock_gettime(CLOCK_REALTIME, &start_time);

  term1_exp = 1.0l / (long double)nle_state->term1.exp_inv;
  term2_exp = 1.0l / (long double)nle_state->term2.exp_inv;
  term3_exp = 1.0l / (long double)nle_state->term3.exp_inv;

  // generate formula strings for each term
  getFormulaStr(nle_config, nle_state, term1_formula_str, nle_state->term1.current_match);
  getFormulaStr(nle_config, nle_state, term2_formula_str, nle_state->term2.current_match);
  getFormulaStr(nle_config, nle_state, term3_formula_str, nle_state->term3.current_match);

  // generate smrfactor string if (1-smr) is enabled
  if (nle_config->smrfactor_1minus_enable == 1) {
    getSmrfStr(nle_config, smrf_str, nle_state->term1.current_smrfactors, nle_state->term1.smrfactor);
  } else {
    smrf_str[0]=0;
  }

  term1_static=((long double)nle_state->term1.current_match->match_up / (long double)nle_state->term1.current_match->match_down) / (long double)nle_state->term1.current_match->static_multiplier;
  term2_static=((long double)nle_state->term2.current_match->match_up / (long double)nle_state->term2.current_match->match_down) / (long double)nle_state->term2.current_match->static_multiplier;
  term3_static=((long double)nle_state->term3.current_match->match_up / (long double)nle_state->term3.current_match->match_down) / (long double)nle_state->term3.current_match->static_multiplier;

  // determine which three variables have the highest uncertainty and will be used as outputs (floated)
#ifdef DEBUG20
  unknowns=selectOutputs(nle_config, nle_state);
#else
  selectOutputs(nle_config, nle_state);
#endif

  // set center/range for non-floated variables
  if (nle_state->all_uses.float_muser == 0) {
    muser_center=(long double)nle_config->smrfactor_mass_user;
    muser_range=(long double)nle_config->smrfactor_mass_user_error;
  }
  if (nle_state->all_uses.float_mh0 == 0) {
    mh0_center=(long double)nle_config->ref_mh0;
    mh0_range=(long double)nle_config->ref_mh0_error;
  }
  if (nle_state->all_uses.float_sin2w == 0) {
    sin2w_center=(long double)nle_config->ref_sin2w;
    sin2w_range=(long double)nle_config->ref_sin2w_error;
  }
  if (nle_state->all_uses.float_mw == 0) {
    mw_center=(long double)nle_config->ref_mw;
    mw_range=(long double)nle_config->ref_mw_error;
  }
  if (nle_state->all_uses.float_sm3 == 0) {
    sm3_center=(long double)nle_config->ref_sm3;
    sm3_range=(long double)nle_config->ref_sm3_error;
  }
  if (nle_state->all_uses.float_G == 0) {
    G_center=(long double)nle_config->ref_G;
    G_range=(long double)nle_config->ref_G_error;
  }
  if (nle_state->all_uses.float_mz == 0) {
    mz_center=(long double)nle_config->ref_mz;
    mz_range=(long double)nle_config->ref_mz_error;
  }
  if (nle_state->all_uses.float_v == 0) {
    v_center=(long double)nle_config->ref_v;
    v_range=(long double)nle_config->ref_v_error;
  }
  if (nle_state->all_uses.float_sm2 == 0) {
    sm2_center=(long double)nle_config->ref_sm2;
    sm2_range=(long double)nle_config->ref_sm2_error;
  }
  if (nle_state->all_uses.float_sm1 == 0) {
    sm1_center=(long double)nle_config->ref_sm1;
    sm1_range=(long double)nle_config->ref_sm1_error;
  }
  if (nle_state->all_uses.float_alpha_em == 0) {
    alpha_center=(long double)nle_config->ref_alpha_em;
    alpha_range=(long double)nle_config->ref_alpha_em_error;
  }

  // systematically try all non-floated input extremes
#ifdef DEBUG20
  printf("debug, Begin phase 2 input loops, exponents: %s, unknowns: %d, float_G: %d, float_v: %d, float_mz: %d, float_mw: %d, float_mh0: %d, float_muser: %d, float_sm1: %d, float_sm2: %d, float_sm3: %d, float_sin2w: %d, float_alpha_em: %d,  float_alpha_w: %d, mw_mz_mode: %d\n", nle_state->exponents_str, unknowns, nle_state->all_uses.float_G, nle_state->all_uses.float_v, nle_state->all_uses.float_mz, nle_state->all_uses.float_mw, nle_state->all_uses.float_mh0, nle_state->all_uses.float_muser, nle_state->all_uses.float_sm1, nle_state->all_uses.float_sm2, nle_state->all_uses.float_sm3, nle_state->all_uses.float_sin2w, nle_state->all_uses.float_alpha_em, nle_state->all_uses.float_alpha_em, nle_state->all_uses.mw_mz_mode);
  printf("debug, term1=%s\n", term1_formula_str);
  printf("debug, term2=%s\n", term2_formula_str);
  printf("debug, term3=%s\n", term3_formula_str);
  if (nle_config->smrfactor_1minus_enable == 1) {
    printf("debug, smrf= %s\n", smrf_str);
  }
  printInputSamples(nle_state);
  printUses(&nle_state->all_uses);
  fflush(stdout);
#endif

  for (alpha_seq=(!nle_state->all_uses.alpha_em || nle_state->all_uses.float_alpha_em || nle_config->phase2_ignore_small_rel_unc); alpha_seq <= 1; alpha_seq++) {
    if (nle_state->all_uses.float_alpha_em == 0) {
      if (nle_config->phase2_ignore_small_rel_unc == 1) {
        alpha=(long double)nle_config->ref_alpha_em;
      } else {
        if (alpha_seq == 0) {
          alpha=(alpha_center - alpha_range);
        } else {
          alpha=(alpha_center + alpha_range);
        }  
      }
    }
    for (sm1_seq=(nle_state->all_uses.float_sm1 || nle_config->phase2_ignore_small_rel_unc); sm1_seq <= 1; sm1_seq++) {
      if (nle_state->all_uses.float_sm1 == 0) {
        if (nle_config->phase2_ignore_small_rel_unc == 1) {
          sm1=(long double)nle_config->ref_sm1;
        } else {
          if (sm1_seq == 0) {
            sm1=(sm1_center - sm1_range);
          } else {
            sm1=(sm1_center + sm1_range);
          }
        }
      }
      for (sm2_seq=nle_state->all_uses.float_sm2; sm2_seq <= 1; sm2_seq++) {
        if (nle_state->all_uses.float_sm2 == 0) {
          if (sm2_seq == 0) {
            sm2=(sm2_center - sm2_range);
          } else {
            sm2=(sm2_center + sm2_range);
          }
        }
        for (v_seq=(!nle_state->all_uses.v || nle_state->all_uses.float_v); v_seq <= 1; v_seq++) {
          if (nle_state->all_uses.float_v == 0) {
            if (v_seq == 0) {
              v=(v_center - v_range);
            } else {
              v=(v_center + v_range);
            }
          }
          for (mz_seq=(!nle_state->all_uses.mz || nle_state->all_uses.float_mz); mz_seq <= 1; mz_seq++) {
            if (nle_state->all_uses.float_mz == 0) {
              if (mz_seq == 0) {
                mz=(mz_center - mz_range);
              } else {
                mz=(mz_center + mz_range);
              }
            }
            for (G_seq=(!nle_state->all_uses.G || nle_state->all_uses.float_G); G_seq <= 1; G_seq++) {
              if (nle_state->all_uses.float_G == 0) { 
                if (G_seq == 0) {
                  G=(G_center - G_range);
                } else {
                  G=(G_center + G_range);
                } 
                mp=(long double)nle_config->ref_kg_to_ev * (long double)sqrtl(nle_config->ref_hbar * nle_config->ref_c / G);
              } 
              for (sm3_seq=nle_state->all_uses.float_sm3; sm3_seq <= 1; sm3_seq++) {
                // sm3 is always used but only floated if necessary
                if (nle_state->all_uses.float_sm3 == 0) {
                  if (sm3_seq == 0) {
                    sm3=(sm3_center - sm3_range);
                  } else {
                    sm3=(sm3_center + sm3_range);
                  }
                }
                for (mw_seq=((!nle_state->all_uses.mw) || nle_state->all_uses.float_mw); mw_seq <= 1; mw_seq++) {
                  if (nle_state->all_uses.float_mw == 0) {
                    if (mw_seq == 0) {
                      mw=(mw_center - mw_range);
                    } else {
                      mw=(mw_center + mw_range);
                    }
                  }
                  for (sin2w_seq=(!nle_state->all_uses.sin2w || nle_state->all_uses.float_sin2w || !nle_state->all_uses.mw_mz_mode); sin2w_seq <= 1; sin2w_seq++) {
                    if ((nle_state->all_uses.float_sin2w == 0) && (nle_state->all_uses.mw_mz_mode == 0)) {
                      if (sin2w_seq == 0) {
                        sin2w=(sin2w_center - sin2w_range);
                      } else {
                        sin2w=(sin2w_center + sin2w_range);
                      }
                    }
                    for (mh0_seq=(!nle_state->all_uses.mh0 || nle_state->all_uses.float_mh0); mh0_seq <= 1; mh0_seq++) {
                      if (nle_state->all_uses.float_mh0 == 0) {
                        if (mh0_seq == 0) {
                          mh0=(mh0_center - mh0_range);
                        } else {
                          mh0=(mh0_center + mh0_range);
                        }
                      }
                      for (muser_seq=(!nle_state->all_uses.m_user || nle_state->all_uses.float_muser || (nle_config->smrfactor_mass_user_error == 0)); muser_seq <= 1; muser_seq++) {
                        if (nle_state->all_uses.float_muser == 0) {
                          if (nle_config->smrfactor_mass_user_error == 0) {
                            muser=(long double)nle_config->smrfactor_mass_user;
                          } else {
                            if (muser_seq == 0) {
                              muser=(muser_center - muser_range);
                            } else {
                              muser=(muser_center + muser_range);
                            }
                          }
                        }

#ifdef DEBUG20
                        printf("debug, Begin    phase 2 samples loop, exponents: %s, alpha_seq: %d, sm1_seq: %d, sm2_seq: %d, v_seq: %d, mz_seq: %d, G_seq: %d, sm3_seq: %d, mw_seq: %d, sin2w_seq: %d, mh0_seq: %d, muser_seq: %d\n", nle_state->exponents_str, alpha_seq, sm1_seq, sm2_seq, v_seq, mz_seq, G_seq, sm3_seq, mw_seq, sin2w_seq, mh0_seq, muser_seq);
                        fflush(stdout);
#endif
                        clock_gettime(CLOCK_REALTIME, &start_time2);
                        precision_last=1.0E99l;
                        //  reset mc test vars and outputs
                        if (nle_state->all_uses.float_alpha_em == 1) {
                          alpha_last=(long double)nle_config->ref_alpha_em;
                          alpha_center=(long double)nle_config->ref_alpha_em;
                          alpha=alpha_center;
                          alpha_range=(long double)nle_config->ref_alpha_em * 0.1l;
                        }
                        if (nle_state->all_uses.float_sm1 == 1) {
                          sm1_last=(long double)nle_config->ref_sm1;
                          sm1_center=(long double)nle_config->ref_sm1;
                          sm1=sm1_center;
                          sm1_range=(long double)nle_config->ref_sm1 * 0.1l;
                        }
                        if (nle_state->all_uses.float_sm2 == 1) {
                          sm2_last=(long double)nle_config->ref_sm2;
                          sm2_center=(long double)nle_config->ref_sm2;
                          sm2=sm2_center;
                          sm2_range=(long double)nle_config->ref_sm2 * 0.1l;
                        }
                        if (nle_state->all_uses.float_v == 1) {
                          v_last=(long double)nle_config->ref_v;
                          v_center=(long double)nle_config->ref_v; 
                          v=v_center;
                          v_range=(long double)nle_config->ref_v * 0.1l;
                        }
                        if (nle_state->all_uses.float_mz == 1) {
                          mz_last=(long double)nle_config->ref_mz;
                          mz_center=(long double)nle_config->ref_mz;
                          mz=mz_center;
                          mz_range=(long double)nle_config->ref_mz * 0.1l;
                        }
                        if (nle_state->all_uses.float_G == 1) {
                          G_last=(long double)nle_config->ref_G;
                          G_center=(long double)nle_config->ref_G; 
                          G=G_center;
                          G_range=(long double)nle_config->ref_G * 0.1l;
                        }
                        if (nle_state->all_uses.float_sm3 == 1) {
                          sm3_last=(long double)nle_config->ref_sm3;
                          sm3_center=(long double)nle_config->ref_sm3;
                          sm3=sm3_center;
                          sm3_range=(long double)nle_config->ref_sm3 * 0.1l;
                        }
                        if (nle_state->all_uses.float_mw == 1) {
                          mw_last=(long double)nle_config->ref_mw;
                          mw_center=(long double)nle_config->ref_mw; 
                          mw=mw_center;
                          mw_range=(long double)nle_config->ref_mw * 0.1l;
                        }
                        if (nle_state->all_uses.float_sin2w == 1) {
                          sin2w_last=(long double)nle_config->ref_sin2w;
                          sin2w_center=(long double)nle_config->ref_sin2w;
                          sin2w=sin2w_center;
                          sin2w_range=(long double)nle_config->ref_sin2w * 0.1l;
                        }
                        if (nle_state->all_uses.float_mh0 == 1) {
                          mh0_last=(long double)nle_config->ref_mh0;
                          mh0_center=(long double)nle_config->ref_mh0;
                          mh0=mh0_center;
                          mh0_range=(long double)nle_config->ref_mh0 * 0.1l;
                        }
                        if (nle_state->all_uses.float_muser == 1) {
                          muser_last=(long double)nle_config->smrfactor_mass_user;
                          muser_center=(long double)nle_config->smrfactor_mass_user;
                          muser=muser_center;
                          muser_range=(long double)nle_config->smrfactor_mass_user * 0.1l;
                        }
                        precision_last=1.0E99l;
                        stalled=0;
                        progress=0;
                        range_multiplier=default_range_multiplier;
                        rangefactor=1.0;
                        for (samples=0; (precision_last > precision_target); samples++) {
                          // periodically check on processing for various conditions and reporting
                          if ((samples > 1) && (progress > 0) && ((samples % process_checkpoint) == 0)) {
#ifdef DEBUG20
                            clock_gettime(CLOCK_REALTIME, &end_time);
                            elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                            if ((samples % 10000000) == 0) { // rate limit periodic debug prints
                              printf ("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last:  %.3Le, precision: %.3Le, sm1_test:  %.3Le, sm2_test:  %.3Le, sm3_test: %.3Le, sm3: %.9Le, sm3_range: %.4Le, G: %.9Le, G_range: %.4Le, v: %.9Le, v_range: %.4Le, sm2: %.9Le, sm2_range: %.4Le, mz: %.9Le, mz_range: %.4Le, mw: %.9Le, mw_range: %.4Le, sin2w: %.9Le, sin2w_range: %.4Le, mh0: %.9Le, mh0_range: %.4Le, muser: %.9Le, muser_range: %.9Le\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last, precision, sm1_test, sm2_test, sm3_test, sm3, sm3_range, G, G_range, v, v_range, sm2, sm2_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range, muser, muser_range);
                              fflush(stdout);
                            }
#endif

                            // check if solution is stuck and needs to be reset
                            if (precision_last > stuck_precision) {
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &end_time);
                              elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, resetting\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              //  reset mc test vars and outputs
                              if (nle_state->all_uses.float_alpha_em == 1) {
                                alpha_last=(long double)nle_config->ref_alpha_em;
                                alpha_center=(long double)nle_config->ref_alpha_em;
                                alpha=alpha_center;
                                alpha_range=(long double)nle_config->ref_alpha_em * 0.1l;
                              }
                              if (nle_state->all_uses.float_sm1 == 1) {
                                sm1_last=(long double)nle_config->ref_sm1;
                                sm1_center=(long double)nle_config->ref_sm1;
                                sm1=sm1_center;
                                sm1_range=(long double)nle_config->ref_sm1 * 0.1l;
                              }
                              if (nle_state->all_uses.float_sm2 == 1) {
                                sm2_last=(long double)nle_config->ref_sm2;
                                sm2_center=(long double)nle_config->ref_sm2;
                                sm2=sm2_center;
                                sm2_range=(long double)nle_config->ref_sm2 * 0.1l;
                              }
                              if (nle_state->all_uses.float_v == 1) {
                                v_last=(long double)nle_config->ref_v;
                                v_center=(long double)nle_config->ref_v;
                                v=v_center;
                                v_range=(long double)nle_config->ref_v * 0.1l;
                              }
                              if (nle_state->all_uses.float_mz == 1) {
                                mz_last=(long double)nle_config->ref_mz;
                                mz_center=(long double)nle_config->ref_mz;
                                mz=mz_center;
                                mz_range=(long double)nle_config->ref_mz * 0.1l;
                              }
                              if (nle_state->all_uses.float_G == 1) {
                                G_last=(long double)nle_config->ref_G;
                                G_center=(long double)nle_config->ref_G; 
                                G=G_center;
                                G_range=(long double)nle_config->ref_G * 0.1l;
                              }
                              if (nle_state->all_uses.float_sm3 == 1) {
                                sm3_last=(long double)nle_config->ref_sm3;
                                sm3_center=(long double)nle_config->ref_sm3;
                                sm3=sm3_center;
                                sm3_range=(long double)nle_config->ref_sm3 * 0.1l;
                              }
                              if (nle_state->all_uses.float_mw == 1) {
                                mw_last=(long double)nle_config->ref_mw;
                                mw_center=(long double)nle_config->ref_mw; 
                                mw=mw_center;
                                mw_range=(long double)nle_config->ref_mw * 0.1l;
                              }
                              if (nle_state->all_uses.float_sin2w == 1) {
                                sin2w_last=(long double)nle_config->ref_sin2w;
                                sin2w_center=(long double)nle_config->ref_sin2w;
                                sin2w=sin2w_center;
                                sin2w_range=(long double)nle_config->ref_sin2w * 0.1l;
                              }
                              if (nle_state->all_uses.float_mh0 == 1) {
                                mh0_last=(long double)nle_config->ref_mh0;
                                mh0_center=(long double)nle_config->ref_mh0;
                                mh0=mh0_center;
                                mh0_range=(long double)nle_config->ref_mh0 * 0.1l;
                              }
                              if (nle_state->all_uses.float_muser == 1) {
                                muser_last=(long double)nle_config->smrfactor_mass_user;
                                muser_center=(long double)nle_config->smrfactor_mass_user;
                                muser=muser_center;
                                muser_range=(long double)nle_config->smrfactor_mass_user * 0.1l;
                              }
                              precision_last=1.0E99l;
                              stalled=0;
                              progress=0;
                              range_multiplier=default_range_multiplier;
                              rangefactor=1.0l;
                            } // end check on stuck/reset
                          } // end if process_checkpoint

                          // check if stalled (too many samples since last progress increment)
                          if (stalled == stalled_limit) {
                            range_multiplier=stalled_range_multiplier; // may be a slow solution, try bigger multiplier
                            rangefactor=rangefactor * stalled_range_multiplier / default_range_multiplier;
#ifdef DEBUG20
                            clock_gettime(CLOCK_REALTIME, &end_time);
                            elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                            printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, stalled\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
#endif
                            if ((rangefactor > 0.1l) || (progress <= 3)) {
                              //  use default ranges
                              if (nle_state->all_uses.float_alpha_em == 1) {
                                alpha_range=(long double)nle_config->ref_alpha_em * 0.1l;
                              }
                              if (nle_state->all_uses.float_sm1 == 1) {
                                sm1_range=(long double)nle_config->ref_sm1 * 0.1l;
                              }
                              if (nle_state->all_uses.float_sm2 == 1) {
                                sm2_range=(long double)nle_config->ref_sm2 * 0.1l;
                              }
                              if (nle_state->all_uses.float_v == 1) {
                                v_range=(long double)nle_config->ref_v * 0.1l;
                              }
                              if (nle_state->all_uses.float_mz == 1) {
                                mz_range=(long double)nle_config->ref_mz * 0.1l;
                              }
                              if (nle_state->all_uses.float_G == 1) {
                                G_range=(long double)nle_config->ref_G * 0.1l;
                              }
                              if (nle_state->all_uses.float_sm3 == 1) {
                                sm3_range=(long double)nle_config->ref_sm3 * 0.1l;
                              }
                              if (nle_state->all_uses.float_mw == 1) {
                                mw_range=(long double)nle_config->ref_mw * 0.1l;
                              }
                              if (nle_state->all_uses.float_sin2w == 1) {
                                sin2w_range=(long double)nle_config->ref_sin2w * 0.1l;
                              }
                              if (nle_state->all_uses.float_mh0 == 1) {
                                mh0_range=(long double)nle_config->ref_mh0 * 0.1l;
                              }
                              if (nle_state->all_uses.float_muser == 1) {
                                muser_range=(long double)nle_config->smrfactor_mass_user * 0.1l;
                              }
                            } else {
                              if (nle_state->all_uses.float_alpha_em == 1) {
                                alpha_range=alpha_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_sm1 == 1) {
                                sm1_range=sm1_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_sm3 == 1) {
                                sm3_range=sm3_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_sm2 == 1) {
                                sm2_range=sm2_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_v == 1) {
                                v_range=v_last * rangefactor;
                              }
                                if (nle_state->all_uses.float_G == 1) {
                                G_range=G_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_mz == 1) {
                                mz_range=mz_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_mw == 1) {
                                mw_range=mw_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_sin2w == 1) {
                                sin2w_range=sin2w_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_mh0 == 1) {
                                mh0_range=mh0_last * rangefactor;
                              }
                              if (nle_state->all_uses.float_muser == 1) {
                                muser_range=muser_last * rangefactor;
                              }
                            }
                          } // end check if stalled
                          stalled++;

                          // guess random values for mc outputs
                          if (nle_state->all_uses.float_alpha_em == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            alpha=((alpha_center - alpha_range) + (r * 2.0l * alpha_range));
                            i=0;
                            while ((alpha < ((long double)nle_config->ref_alpha_em * 0.5l)) || (alpha > ((long double)nle_config->ref_alpha_em * 1.5l))) { // sanity check to help convergence
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, alpha range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                alpha_last=(long double)nle_config->ref_alpha_em;
                                alpha_center=(long double)nle_config->ref_alpha_em;
                                alpha_range=(long double)nle_config->ref_alpha_em * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              alpha=((alpha_center - alpha_range) + (r * 2.0l * alpha_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_sm1 == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            sm1=((sm1_center - sm1_range) + (r * 2.0l * sm1_range));
                            i=0;
                            while ((sm1 < ((long double)nle_config->ref_sm1 * 0.5l)) || (sm1 > ((long double)nle_config->ref_sm1 * 1.5l))) { // sanity check to help convergence
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, sm1 range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                sm1_last=(long double)nle_config->ref_sm1;
                                sm1_center=(long double)nle_config->ref_sm1;
                                sm1_range=(long double)nle_config->ref_sm1 * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              sm1=((sm1_center - sm1_range) + (r * 2.0l * sm1_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_sm2 == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            sm2=((sm2_center - sm2_range) + (r * 2.0l * sm2_range));
                            i=0;
                            while ((sm2 < ((long double)nle_config->ref_sm2 * 0.5l)) || (sm2 > ((long double)nle_config->ref_sm2 * 1.5l))) { // sanity check to help convergence
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, sm2 range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                sm2_last=(long double)nle_config->ref_sm2;
                                sm2_center=(long double)nle_config->ref_sm2;
                                sm2_range=(long double)nle_config->ref_sm2 * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              sm2=((sm2_center - sm2_range) + (r * 2.0l * sm2_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_v == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            v=((v_center - v_range) + (r * 2.0l * v_range));
                            i=0;
                            while ((v < ((long double)nle_config->ref_v * 0.5l)) || (v > (nle_config->ref_v * 1.5l))) { // sanity check to help convergence
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, v range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                v_last=(long double)nle_config->ref_v;
                                v_center=(long double)nle_config->ref_v;
                                v_range=(long double)nle_config->ref_v * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              v=((v_center - v_range) + (r * 2.0l * v_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_mz == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            mz=((mz_center - mz_range) + (r * 2.0l * mz_range));
                            i=0;
                            while ((mz < ((long double)nle_config->ref_mz * 0.5l)) || (mz > ((long double)nle_config->ref_mz * 1.5l)) || (mw >= mz) || (mz >= mh0)) { // sanity check to help convergence and prevent mw >= mz or mz >= mh0
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, mz range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                mz_last=(long double)nle_config->ref_mz;
                                mz_center=(long double)nle_config->ref_mz;
                                mz_range=(long double)nle_config->ref_mz * 0.1l;
                                // if mw >= mz and mw is floated then we must reset mw also to clear error
                                if ((mw >= mz) && (nle_state->all_uses.float_mw == 1)) {
                                  mw_last=(long double)nle_config->ref_mw;
                                  mw_center=(long double)nle_config->ref_mw; 
                                  mw=mw_center;
                                  mw_range=(long double)nle_config->ref_mw * 0.1l;
                                }
                                if ((mz >= mh0) && (nle_state->all_uses.float_mh0 == 1)) {
                                  mh0_last=(long double)nle_config->ref_mh0;
                                  mh0_center=(long double)nle_config->ref_mh0;
                                  mh0=mh0_center;
                                  mh0_range=(long double)nle_config->ref_mh0 * 0.1l;
                                }
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              mz=((mz_center - mz_range) + (r * 2.0l * mz_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_G == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            G=((G_center - G_range) + (r * 2.0l * G_range));
                            i=0;
                            while ((G < ((long double)nle_config->ref_G * 0.5l)) || (G > ((long double)nle_config->ref_G * 1.5l))) { // sanity check to help convergence 
                              if (i > 50) {  // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, G range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                G_last=(long double)nle_config->ref_G;
                                G_center=(long double)nle_config->ref_G;
                                G_range=(long double)nle_config->ref_G * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              G=((G_center - G_range) + (r * 2.0l * G_range));
                              i++;
                            }
                            mp=(long double)nle_config->ref_kg_to_ev * (long double)sqrtl(nle_config->ref_hbar * nle_config->ref_c / G);
                          }

                          if (nle_state->all_uses.float_sm3 == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            sm3=((sm3_center - sm3_range) + (r * 2.0l * sm3_range));
                            i=0;
                            while ((sm3 < ((long double)nle_config->ref_sm3 * 0.5l)) || (sm3 > ((long double)nle_config->ref_sm3 * 1.5l))) { // sanity check to help convergence
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, sm3 range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                sm3_last=(long double)nle_config->ref_sm3;
                                sm3_center=(long double)nle_config->ref_sm3;
                                sm3_range=(long double)nle_config->ref_sm3 * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              sm3=((sm3_center - sm3_range) + (r * 2.0l * sm3_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_mw == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            mw=((mw_center - mw_range) + (r * 2.0l * mw_range));
                            i=0;
                            while ((mw < ((long double)nle_config->ref_mw * 0.5l)) || (mw > ((long double)nle_config->ref_mw * 1.5l)) || (mw >= mz) || (mw >= mh0)) { // sanity check to help convergence and prevent mw >= mz or mw >= mh0
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, mw range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                mw_last=(long double)nle_config->ref_mw;
                                mw_center=(long double)nle_config->ref_mw;
                                mw_range=(long double)nle_config->ref_mw * 0.1l;
                                // if mw >= mz and mz is floated then we must reset mz also to clear error
                                if ((mw >= mz) && (nle_state->all_uses.float_mz == 1)) {
                                  mz_last=(long double)nle_config->ref_mz;
                                  mz_center=(long double)nle_config->ref_mz;
                                  mz=mz_center;
                                  mz_range=(long double)nle_config->ref_mz * 0.1l;
                                }
                                if ((mw >= mh0) && (nle_state->all_uses.float_mh0 == 1)) {
                                  mh0_last=(long double)nle_config->ref_mh0;
                                  mh0_center=(long double)nle_config->ref_mh0;
                                  mh0=mh0_center;
                                  mh0_range=(long double)nle_config->ref_mh0 * 0.1l;
                                }
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              mw=((mw_center - mw_range) + (r * 2.0l * mw_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_sin2w == 1) {  
                            r=(long double)pcg_ldrand64(nle_state);
                            sin2w=((sin2w_center - sin2w_range) + (r * 2.0l * sin2w_range));
                            i=0;
                            while ((sin2w < ((long double)nle_config->ref_sin2w * 0.5l)) || (sin2w > ((long double)nle_config->ref_sin2w * 1.5l))) { // sanity check to help convergence
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, sin2w range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                sin2w_last=(long double)nle_config->ref_sin2w;
                                sin2w_center=(long double)nle_config->ref_sin2w; 
                                sin2w_range=(long double)nle_config->ref_sin2w * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              sin2w=((sin2w_center - sin2w_range) + (r * 2.0l * sin2w_range));
                              i++;
                            }
                          } // end nle_state->all_uses sin2w

                          if (nle_state->all_uses.float_mh0 == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            mh0=((mh0_center - mh0_range) + (r * 2.0l * mh0_range));
                            i=0;
                            while ((mh0 < ((long double)nle_config->ref_mh0 * 0.5l)) || (mh0 > ((long double)nle_config->ref_mh0 * 1.5l)) || (mz >= mh0) || (mw >= mh0)) { // sanity check to help convergence and prevent mz >= mh0 or mw >= mh0
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, mh0 range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                mh0_last=(long double)nle_config->ref_mh0;
                                mh0_center=(long double)nle_config->ref_mh0;
                                mh0_range=(long double)nle_config->ref_mh0 * 0.1l;
                                if ((mz >= mh0) && (nle_state->all_uses.float_mz == 1)) {
                                  mz_last=(long double)nle_config->ref_mz;
                                  mz_center=(long double)nle_config->ref_mz;
                                  mz=mz_center;
                                  mz_range=(long double)nle_config->ref_mz * 0.1l;
                                }
                                if ((mw >= mh0) && (nle_state->all_uses.float_mw == 1)) {
                                  mw_last=(long double)nle_config->ref_mw;
                                  mw_center=(long double)nle_config->ref_mw;
                                  mw=mw_center;
                                  mw_range=(long double)nle_config->ref_mw * 0.1l;
                                }
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              mh0=((mh0_center - mh0_range) + (r * 2.0l * mh0_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.float_muser == 1) {
                            r=(long double)pcg_ldrand64(nle_state);
                            muser=((muser_center - muser_range) + (r * 2.0l * muser_range));
                            i=0;
                            while ((muser < ((long double)nle_config->smrfactor_mass_user * 0.5l)) || (muser > ((long double)nle_config->smrfactor_mass_user * 1.5l))) { // sanity check to help convergence
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last: %.3Le, muser range error\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                muser_last=(long double)nle_config->smrfactor_mass_user;
                                muser_center=(long double)nle_config->smrfactor_mass_user;
                                muser_range=(long double)nle_config->smrfactor_mass_user * 0.1l;
                              }
                              r=(long double)pcg_ldrand64(nle_state);
                              muser=((muser_center - muser_range) + (r * 2.0l * muser_range));
                              i++;
                            }
                          }

                          if (nle_state->all_uses.mw_mz_mode == 1) {
                            // cos2w derived from mw/mz, sin2w derived from cos2w
                            cos2w=powl((mw / mz), 2.0l);
                            sin2w=1.0l - cos2w;
                          }  else {
                            // cos2w derived from sin2w , sin2w is either sequenced or floated above
                            cos2w=1.0l - sin2w;
                          }

                          // set solution mass ratio reference mass and factors
                          if (nle_state->term1.current_match->smrfactor_mass_id == 0) {
                            term1_smrfactor_mass=mp;
                          } else if (nle_state->term1.current_match->smrfactor_mass_id == 1) {
                            term1_smrfactor_mass=v;
                          } else if (nle_state->term1.current_match->smrfactor_mass_id == 2) {
                            term1_smrfactor_mass=mz;
                          } else if (nle_state->term1.current_match->smrfactor_mass_id == 3) {
                            term1_smrfactor_mass=mw;
                          } else if (nle_state->term1.current_match->smrfactor_mass_id == 4) {
                            term1_smrfactor_mass=mh0;
                          } else if (nle_state->term1.current_match->smrfactor_mass_id == 5) {
                            term1_smrfactor_mass=muser;
                          }
                          if (nle_state->term2.current_match->smrfactor_mass_id == 0) {
                            term2_smrfactor_mass=mp;
                          } else if (nle_state->term2.current_match->smrfactor_mass_id == 1) {
                            term2_smrfactor_mass=v;
                          } else if (nle_state->term2.current_match->smrfactor_mass_id == 2) {
                            term2_smrfactor_mass=mz;
                          } else if (nle_state->term2.current_match->smrfactor_mass_id == 3) {
                            term2_smrfactor_mass=mw;
                          } else if (nle_state->term2.current_match->smrfactor_mass_id == 4) {
                            term2_smrfactor_mass=mh0;
                          } else if (nle_state->term2.current_match->smrfactor_mass_id == 5) {
                            term2_smrfactor_mass=muser;
                          }
                          if (nle_state->term3.current_match->smrfactor_mass_id == 0) {
                            term3_smrfactor_mass=mp;
                          } else if (nle_state->term3.current_match->smrfactor_mass_id == 1) {
                            term3_smrfactor_mass=v;
                          } else if (nle_state->term3.current_match->smrfactor_mass_id == 2) {
                            term3_smrfactor_mass=mz;
                          } else if (nle_state->term3.current_match->smrfactor_mass_id == 3) {
                            term3_smrfactor_mass=mw;
                          } else if (nle_state->term3.current_match->smrfactor_mass_id == 4) {
                            term3_smrfactor_mass=mh0;
                          } else if (nle_state->term3.current_match->smrfactor_mass_id == 5) {
                            term3_smrfactor_mass=muser;
                          }

                          // set mass component of each term
                          if (nle_config->smrfactor_1minus_enable == 1) {
                            if (nle_state->smrfactor_mass_configuration == 1) {
                              smrf_sm1=nle_state->term1.smrfactor * sm1 / term1_smrfactor_mass;
                              smrf_sm2=nle_state->term1.smrfactor * sm2 / term1_smrfactor_mass;
                              smrf_sm3=nle_state->term1.smrfactor * sm3 / term1_smrfactor_mass;
                            } else if (nle_state->smrfactor_mass_configuration == 0) {
                              smrf_sm1=nle_state->term1.smrfactor * term1_smrfactor_mass / sm1;
                              smrf_sm2=nle_state->term1.smrfactor * term1_smrfactor_mass / sm2;
                              smrf_sm3=nle_state->term1.smrfactor * term1_smrfactor_mass / sm3;
                            }
                            // check if (1-smr) is negative for sm1 and invert inside and outside radical
                            if ((1.0l - smrf_sm1) < 0.0l) {
                              term1_mass_sm1=-powl(-(1.0l - smrf_sm1), term1_exp);
                            } else {
                              term1_mass_sm1=powl((1.0l - smrf_sm1), term1_exp);
                            }
                            // check if (1-smr) is negative for sm2 and invert inside and outside radical
                            if ((1.0l - smrf_sm2) < 0.0l) {
                              term1_mass_sm2=-powl(-(1.0l - smrf_sm2), term1_exp);
                            } else {
                              term1_mass_sm2=powl((1.0l - smrf_sm2), term1_exp);
                            }
                            // check if (1-smr) is negative for sm3 and invert inside and outside radical
                            if ((1.0l - smrf_sm3) < 0.0l) {
                              term1_mass_sm3=-powl(-(1.0l - smrf_sm3), term1_exp);
                            } else {
                              term1_mass_sm3=powl((1.0l - smrf_sm3), term1_exp);
                            }
                            term2_mass_sm1=powl(smrf_sm1, term2_exp);
                            term2_mass_sm2=powl(smrf_sm2, term2_exp);
                            term2_mass_sm3=powl(smrf_sm3, term2_exp);
                          } else {
                            term1_mass_sm1=powl((sm1 / term1_smrfactor_mass), term1_exp);
                            term1_mass_sm2=powl((sm2 / term1_smrfactor_mass), term1_exp);
                            term1_mass_sm3=powl((sm3 / term1_smrfactor_mass), term1_exp);
                            term2_mass_sm1=powl((sm1 / term2_smrfactor_mass), term2_exp);
                            term2_mass_sm2=powl((sm2 / term2_smrfactor_mass), term2_exp);
                            term2_mass_sm3=powl((sm3 / term2_smrfactor_mass), term2_exp);
                          }
                          if (nle_config->nle_mode == 3) {
                            term3_mass_sm1=powl((sm1 / term3_smrfactor_mass), term3_exp);
                            term3_mass_sm2=powl((sm2 / term3_smrfactor_mass), term3_exp);
                            term3_mass_sm3=powl((sm3 / term3_smrfactor_mass), term3_exp);
                          }

                          // set outfactor_rmr masses for each term if configured (independent from and not to be confused with rmrfactor masses in (1-rmr-smr) mode)
                          if (nle_state->term1.current_match->outfactor_rmr_exp_up != 0) {
                            if (nle_state->term1.current_match->outfactor_rmr_mass_id_up == 0) {
                              outfactor_rmr_mass_up=mp;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_up == 1) {
                              outfactor_rmr_mass_up=v;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_up == 2) {
                              outfactor_rmr_mass_up=mz;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_up == 3) {
                              outfactor_rmr_mass_up=mw;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_up == 4) {
                              outfactor_rmr_mass_up=mh0;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_up == 5) {
                              outfactor_rmr_mass_up=muser;
                            }
                            if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 0) {
                              outfactor_rmr_mass_down=mp;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 1) {
                              outfactor_rmr_mass_down=v;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 2) {
                              outfactor_rmr_mass_down=mz;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 3) {
                              outfactor_rmr_mass_down=mw;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 4) {
                              outfactor_rmr_mass_down=mh0;
                            } else if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 5) {
                              outfactor_rmr_mass_down=muser;
                            }
                            term1_rmr=powl((outfactor_rmr_mass_down / outfactor_rmr_mass_up), ((long double)nle_state->term1.current_match->outfactor_rmr_exp_up / (long double)nle_state->term1.current_match->outfactor_rmr_exp_down));
                          } else {
                            term1_rmr=1.0l;
                          }
                          if (nle_state->term1.current_match->outfactor_sin2w_exp_up != 0) {
                            term1_sin2w=powl(sin2w, ((long double)nle_state->term1.current_match->outfactor_sin2w_exp_up / (long double)nle_state->term1.current_match->outfactor_sin2w_exp_down));
                          } else {
                            term1_sin2w=1.0l;
                          }
                          if (nle_state->term1.current_match->outfactor_cos2w_exp_up != 0) {
                            term1_cos2w=powl(cos2w, ((long double)nle_state->term1.current_match->outfactor_cos2w_exp_up / (long double)nle_state->term1.current_match->outfactor_cos2w_exp_down));
                          } else {
                            term1_cos2w=1.0l;
                          }
                          if (nle_state->term2.current_match->outfactor_rmr_exp_up != 0) {
                            if (nle_state->term2.current_match->outfactor_rmr_mass_id_up == 0) {
                              outfactor_rmr_mass_up=mp;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_up == 1) {
                              outfactor_rmr_mass_up=v;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_up == 2) {
                              outfactor_rmr_mass_up=mz;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_up == 3) {
                              outfactor_rmr_mass_up=mw;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_up == 4) {
                              outfactor_rmr_mass_up=mh0;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_up == 5) {
                              outfactor_rmr_mass_up=muser;
                            } 
                            if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 0) { 
                              outfactor_rmr_mass_down=mp;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_down == 1) {
                              outfactor_rmr_mass_down=v;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_down == 2) {
                              outfactor_rmr_mass_down=mz;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_down == 3) {
                              outfactor_rmr_mass_down=mw;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_down == 4) {
                              outfactor_rmr_mass_down=mh0;
                            } else if (nle_state->term2.current_match->outfactor_rmr_mass_id_down == 5) {
                              outfactor_rmr_mass_down=muser;
                            } 
                            term2_rmr=powl((outfactor_rmr_mass_down / outfactor_rmr_mass_up), ((long double)nle_state->term2.current_match->outfactor_rmr_exp_up / (long double)nle_state->term2.current_match->outfactor_rmr_exp_down));
                          } else {
                            term2_rmr=1.0l;
                          }
                          if (nle_state->term2.current_match->outfactor_sin2w_exp_up != 0) {
                            term2_sin2w=powl(sin2w, ((long double)nle_state->term2.current_match->outfactor_sin2w_exp_up / (long double)nle_state->term2.current_match->outfactor_sin2w_exp_down));
                          } else {
                            term2_sin2w=1.0l;
                          }
                          if (nle_state->term2.current_match->outfactor_cos2w_exp_up != 0) {
                            term2_cos2w=powl(cos2w, ((long double)nle_state->term2.current_match->outfactor_cos2w_exp_up / (long double)nle_state->term2.current_match->outfactor_cos2w_exp_down));
                          } else {
                            term2_cos2w=1.0l;
                          }
                          if (nle_state->term3.current_match->outfactor_rmr_exp_up != 0) {
                            if (nle_state->term3.current_match->outfactor_rmr_mass_id_up == 0) {
                              outfactor_rmr_mass_up=mp;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_up == 1) {
                              outfactor_rmr_mass_up=v;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_up == 2) {
                              outfactor_rmr_mass_up=mz;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_up == 3) {
                              outfactor_rmr_mass_up=mw;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_up == 4) {
                              outfactor_rmr_mass_up=mh0;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_up == 5) {
                              outfactor_rmr_mass_up=muser;
                            } 
                            if (nle_state->term1.current_match->outfactor_rmr_mass_id_down == 0) { 
                              outfactor_rmr_mass_down=mp;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_down == 1) {
                              outfactor_rmr_mass_down=v;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_down == 2) {
                              outfactor_rmr_mass_down=mz;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_down == 3) {
                              outfactor_rmr_mass_down=mw;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_down == 4) {
                              outfactor_rmr_mass_down=mh0;
                            } else if (nle_state->term3.current_match->outfactor_rmr_mass_id_down == 5) {
                              outfactor_rmr_mass_down=muser;
                            } 
                            term3_rmr=powl((outfactor_rmr_mass_down / outfactor_rmr_mass_up), ((long double)nle_state->term3.current_match->outfactor_rmr_exp_up / (long double)nle_state->term3.current_match->outfactor_rmr_exp_down));
                          } else {
                            term3_rmr=1.0l;
                          }
                          if (nle_state->term3.current_match->outfactor_sin2w_exp_up != 0) {
                            term3_sin2w=powl(sin2w, ((long double)nle_state->term3.current_match->outfactor_sin2w_exp_up / (long double)nle_state->term3.current_match->outfactor_sin2w_exp_down));
                          } else {
                            term3_sin2w=1.0l;
                          }
                          if (nle_state->term3.current_match->outfactor_cos2w_exp_up != 0) {
                            term3_cos2w=powl(cos2w, ((long double)nle_state->term3.current_match->outfactor_cos2w_exp_up / (long double)nle_state->term3.current_match->outfactor_cos2w_exp_down));
                          } else {
                            term3_cos2w=1.0l;
                          }
                          term1_coefficient=(term1_static * term1_rmr / term1_sin2w) / term1_cos2w;
                          term2_coefficient=(term2_static * term2_rmr / term2_sin2w) / term2_cos2w;
                          term3_coefficient=(term3_static * term3_rmr / term3_sin2w) / term3_cos2w;

                          // Combine factors for each term and generate test value for each solution mass
                          if (nle_config->nle_mode == 2) {
                            // for 2-term mixed mode, we will use term3 as a pseudo term for the mixing of terms 1 and 2
                            sm1_test_term1=term1_coefficient * term1_mass_sm1 * term1_coefficient * term1_mass_sm1;
                            sm1_test_term2=term2_coefficient * term2_mass_sm1 * term2_coefficient * term2_mass_sm1;
                            sm1_test_term3=term3_coefficient * term1_coefficient * term1_mass_sm1 * term2_coefficient * term2_mass_sm1;
                            sm2_test_term1=term1_coefficient * term1_mass_sm2 * term1_coefficient * term1_mass_sm2;
                            sm2_test_term2=term2_coefficient * term2_mass_sm2 * term2_coefficient * term2_mass_sm2;
                            sm2_test_term3=term3_coefficient * term1_coefficient * term1_mass_sm2 * term2_coefficient * term2_mass_sm2;
                            sm3_test_term1=term1_coefficient * term1_mass_sm3 * term1_coefficient * term1_mass_sm3;
                            sm3_test_term2=term2_coefficient * term2_mass_sm3 * term2_coefficient * term2_mass_sm3;
                            sm3_test_term3=term3_coefficient * term1_coefficient * term1_mass_sm3 * term2_coefficient * term2_mass_sm3;
                            if (nle_config->smrfactor_1minus_enable == 1) { // two different mixing polarity options for (1-smr) mode
                              if (nle_state->nle_mixing_polarity == 0) {
                                sm1_test=sm1_test_term1 + sm1_test_term2 - sm1_test_term3 - 1.0l;
                                sm2_test=sm2_test_term1 + sm2_test_term2 - sm2_test_term3 - 1.0l;
                                sm3_test=sm3_test_term1 + sm3_test_term2 - sm3_test_term3 - 1.0l;
                              } else if (nle_state->nle_mixing_polarity == 1) {
                                sm1_test=sm1_test_term1 + sm1_test_term2 + sm1_test_term3 - 1.0l;
                                sm2_test=sm2_test_term1 + sm2_test_term2 + sm2_test_term3 - 1.0l;
                                sm3_test=sm3_test_term1 + sm3_test_term2 + sm3_test_term3 - 1.0l;
                              }
                            } else { // 2-term without (1-smr)
                              sm1_test=sm1_test_term1 + sm1_test_term2 - sm1_test_term3 - 1.0l;
                              sm2_test=sm2_test_term1 + sm2_test_term2 - sm2_test_term3 - 1.0l;
                              sm3_test=sm3_test_term1 + sm3_test_term2 - sm3_test_term3 - 1.0l;
                            }
                          } else if (nle_config->nle_mode == 3) {
                            sm1_test_term1=term1_coefficient * term1_mass_sm1;
                            sm1_test_term2=term2_coefficient * term2_mass_sm1;
                            sm1_test_term3=term3_coefficient * term3_mass_sm1;
                            sm2_test_term1=term1_coefficient * term1_mass_sm2;
                            sm2_test_term2=term2_coefficient * term2_mass_sm2;
                            sm2_test_term3=term3_coefficient * term3_mass_sm2;
                            sm3_test_term1=term1_coefficient * term1_mass_sm3;
                            sm3_test_term2=term2_coefficient * term2_mass_sm3;
                            sm3_test_term3=term3_coefficient * term3_mass_sm3;
                            sm1_test=sm1_test_term1 - sm1_test_term2 + sm1_test_term3 - 1.0l;
                            sm2_test=sm2_test_term1 - sm2_test_term2 + sm2_test_term3 - 1.0l;
                            sm3_test=sm3_test_term1 - sm3_test_term2 + sm3_test_term3 - 1.0l;
                          }

#ifdef DEBUG23
                          printf("debug, exponents: %s, samples: %lld, sm1_test_term1: %.6Le, sm1_test_term2: %.6Le, sm1_test_term3: %.6Le, sm2_test_term1: %.6Le, sm2_test_term2: %.6Le, sm2_test_term3: %.6Le, sm3_test_term1: %.6Le, sm3_test_term2: %.6Le, sm3_test_term3: %.6Le\n", nle_state->exponents_str, samples, sm1_test_term1, sm1_test_term2, sm1_test_term3, sm2_test_term1, sm2_test_term2, sm2_test_term3, sm3_test_term1, sm3_test_term2, sm3_test_term3); 
                          //printf("debug, exponents: %s, samples: %lld, term1_c:   %.6Le, term1static:   %.6Le, term1s2w:   %.6Le, term1c2w:   %.6Le, term1s2wupout:   %d, term1s2wdownout:   %d, term1c2wupout:   %d, term1c2wdownout:   %d\n", nle_state->exponents_str, samples, term1_coefficient, term1_static, term1_sin2w, term1_sin2w, nle_state->term1.current_match->outfactor_sin2w_exp_up, nle_state->term1.current_match->outfactor_sin2w_exp_down, nle_state->term1.current_match->outfactor_cos2w_exp_up, nle_state->term1.current_match->outfactor_cos2w_exp_down);
                          //printf("debug, exponents: %s, samples: %lld, term2_c: %.6Le, term2static: %.6Le, term2s2w: %.6Le, term2c2w: %.6Le, term2s2wupout: %d, term2s2wdownout: %d, term2c2wupout: %d, term2c2wdownout: %d\n", nle_state->exponents_str, samples, term2_coefficient, term2_static, term2_sin2w, term2_sin2w, nle_state->term2.current_match->outfactor_sin2w_exp_up, nle_state->term2.current_match->outfactor_sin2w_exp_down, nle_state->term2.current_match->outfactor_cos2w_exp_up, nle_state->term2.current_match->outfactor_cos2w_exp_down);
                          //printf("debug, exponents: %s, samples: %lld, term3_c:  %.6Le, term3static:  %.6Le, term3s2w:  %.6Le, term3c2w:  %.6Le, term3s2wupout:  %d, term3s2wdownout:  %d, term3c2wupout:  %d, term3c2wdownout:  %d\n", nle_state->exponents_str, samples, term3_coefficient, term3_static, term3_sin2w, term3_sin2w, nle_state->term3.current_match->outfactor_sin2w_exp_up, nle_state->term3.current_match->outfactor_sin2w_exp_down, nle_state->term3.current_match->outfactor_cos2w_exp_up, nle_state->term3.current_match->outfactor_cos2w_exp_down);
                          //printf("debug, exponents: %s, samples: %lld, term1_mass_sm1: %.6Le, term1_mass_sm2: %.6Le, term1_mass_sm3: %.6Le, term2_mass_sm1: %.6Le, term2_mass_sm2: %.6Le, term2_mass_sm3: %.6Le, term3_mass_sm1: %.6Le, term3_mass_sm2: %.6Le, term3_mass_sm3: %.6Le\n", nle_state->exponents_str, samples, term1_mass_sm1, term1_mass_sm2, term1_mass_sm3, term2_mass_sm1, term2_mass_sm2, term2_mass_sm3, term3_mass_sm1, term3_mass_sm2, term3_mass_sm3);
                          //printf("debug, exponents: %s, samples: %lld, sm1_test:  %.3Le, sm2_test:  %.3Le, sm3_test: %.3Le, term1_c: %.3Le, term2_c: %.3Le, term3_c: %.3Le\n", nle_state->exponents_str, samples, sm1_test, sm2_test, sm3_test, term1_coefficient, term2_coefficient, term3_coefficient);
                          //printf("debug, exponents: %s, samples: %lld, term1 up: %d, term1 down: %d, term2 up: %d, term2 down: %d, term3 up: %d, term3 down: %d,  term1 static_multiplier: %.9e, , term2 static_multiplier: %.9e, term3 static_multiplier: %.9e\n", nle_state->exponents_str, samples, nle_state->term1.current_match->match_up, nle_state->term1.current_match->match_down, nle_state->term2.current_match->match_up, nle_state->term2.current_match->match_down, nle_state->term3.current_match->match_up, nle_state->term3.current_match->match_down, nle_state->term1.current_match->static_multiplier, nle_state->term2.current_match->static_multiplier, nle_state->term3.current_match->static_multiplier);

                          fflush(stdout);
#endif
                          if ((progress < ratio_grace_period) || (((fabsl(sm1_test) / fabsl(sm2_test)) < test_ratio) && ((fabsl(sm2_test) / fabsl(sm1_test)) < test_ratio))) {
                            if ((progress < ratio_grace_period) || (((fabsl(sm1_test) / fabsl(sm3_test)) < test_ratio) && ((fabsl(sm3_test) / fabsl(sm1_test)) < test_ratio) &&\
                               ((fabsl(sm2_test) / fabsl(sm3_test)) < test_ratio) && ((fabsl(sm3_test) / fabsl(sm2_test)) < test_ratio))) {
                              precision=fabsl(sm1_test) + fabsl(sm2_test) + fabsl(sm3_test);
#ifdef DEBUG22
                              clock_gettime(CLOCK_REALTIME, &end_time);
                              elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                              printf ("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: %d, rangefactor: %.9Le, precision_last:  %.3Le, precision: %.3Le, sm1_test:  %.3Le, sm2_test:  %.3Le, sm3_test: %.3Le, alpha: %.9Le, alpha_range: %.9Le, sm1: %.9Le, sm1_range: %.9Le, sm3: %.9Le, sm3_range: %.4Le, G: %.9Le, G_range: %.4Le, v: %.9Le, v_range: %.4Le, sm2: %.9Le, sm2_range: %.4Le, mz: %.9Le, mz_range: %.4Le, mw: %.9Le, mw_range: %.4Le, sin2w: %.9Le, sin2w_range: %.4Le, mh0: %.9Le, mh0_range: %.4Le, muser: %.9Le, muser_range: %.9Le\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last, precision, sm1_test, sm2_test, sm3_test, alpha, alpha_range, sm1, sm1_range, sm3, sm3_range, G, G_range, v, v_range, sm2, sm2_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range, muser, muser_range);
                              fflush(stdout);
#endif
                              if (precision < precision_last) {
                                progress++;
                                stalled=0;

                                // update precision and output variables state
                                precision_last=precision;
                                alpha_last=alpha;
                                sm1_last=sm1;
                                sm2_last=sm2;
                                sm3_last=sm3;
                                v_last=v;
                                G_last=G;
                                mz_last=mz;
                                mw_last=mw;
                                sin2w_last=sin2w;
                                mh0_last=mh0;
                                muser_last=muser;

                                // determine new search range for each output variable
                                if (fabsl(sm1_test) > fabsl(sm2_test)) {
                                  worst_test=fabsl(sm1_test);
                                } else {
                                  worst_test=fabsl(sm2_test);
                                }
                                if (fabsl(sm3_test) > worst_test) {
                                  worst_test=fabsl(sm3_test);
                                }
                                range_multiplier=default_range_multiplier;
                                rangefactor=worst_test * range_multiplier;
                                if ((rangefactor > 0.1l) || (progress <= 3)) {
                                  //  use default ranges
                                  if (nle_state->all_uses.float_alpha_em == 1) {
                                    alpha_center=alpha_last;
                                    alpha_range=(long double)nle_config->ref_alpha_em * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_sm1 == 1) {
                                    sm1_center=sm1_last;
                                    sm1_range=(long double)nle_config->ref_sm1 * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_sm2 == 1) {
                                  sm2_center=sm2_last;
                                    sm2_range=(long double)nle_config->ref_sm2 * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_v == 1) {
                                    v_center=v_last;
                                    v_range=(long double)nle_config->ref_v * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_mz == 1) {
                                    mz_center=mz_last;
                                    mz_range=(long double)nle_config->ref_mz * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_G == 1) {
                                    G_center=G_last;
                                    G_range=(long double)nle_config->ref_G * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_sm3 == 1) {
                                    sm3_center=sm3_last;
                                    sm3_range=(long double)nle_config->ref_sm3 * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_mw == 1) {
                                    mw_center=mw_last;
                                    mw_range=(long double)nle_config->ref_mw * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_sin2w == 1) {
                                    sin2w_center=sin2w_last;
                                    sin2w_range=(long double)nle_config->ref_sin2w * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_mh0 == 1) {
                                    mh0_center=mh0_last;
                                    mh0_range=(long double)nle_config->ref_mh0 * 0.1l;
                                  }
                                  if (nle_state->all_uses.float_muser == 1) {
                                    muser_center=muser_last;
                                    muser_range=(long double)nle_config->smrfactor_mass_user * 0.1l;
                                  }
                                } else { 
                                  if (nle_state->all_uses.float_alpha_em == 1) {
                                    alpha_center=alpha_last;
                                    alpha_range_new=alpha_last * rangefactor;
                                    alpha_range=((alpha_range + alpha_range_new + alpha_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_sm1 == 1) {
                                    sm1_center=sm1_last;
                                    sm1_range_new=sm1_last * rangefactor;
                                    sm1_range=((sm1_range + sm1_range_new + sm1_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_sm3 == 1) {
                                    sm3_center=sm3_last;
                                    sm3_range_new=sm3_last * rangefactor;
                                    sm3_range=((sm3_range + sm3_range_new + sm3_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_sm2 == 1) {
                                    sm2_center=sm2_last;
                                    sm2_range_new=sm2_last * rangefactor;
                                    sm2_range=((sm2_range + sm2_range_new + sm2_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_v == 1) {
                                    v_center=v_last;
                                    v_range_new=v_last * rangefactor;
                                    v_range=((v_range + v_range_new + v_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_G == 1) {
                                    G_center=G_last;
                                    G_range_new=G_last * rangefactor;
                                    G_range=((G_range + G_range_new + G_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_mz == 1) {
                                    mz_center=mz_last;
                                    mz_range_new=mz_last * rangefactor;
                                    mz_range=((mz_range + mz_range_new + mz_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_mw == 1) {
                                    mw_center=mw_last;
                                    mw_range_new=mw_last * rangefactor;
                                    mw_range=((mw_range + mw_range_new + mw_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_sin2w == 1) {
                                    sin2w_center=sin2w_last;
                                    sin2w_range_new=sin2w_last * rangefactor;
                                    sin2w_range=((sin2w_range + sin2w_range_new + sin2w_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_mh0 == 1) {
                                    mh0_center=mh0_last;
                                    mh0_range_new=mh0_last * rangefactor;
                                    mh0_range=((mh0_range + mh0_range_new + mh0_range_new) / 3.0l);
                                  }
                                  if (nle_state->all_uses.float_muser == 1) {
                                    muser_center=muser_last;
                                    muser_range_new=muser_last * rangefactor;
                                    muser_range=((muser_range + muser_range_new + muser_range_new) / 3.0l);
                                  }
                                }
#ifdef DEBUG21
                                clock_gettime(CLOCK_REALTIME, &end_time);
                                elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                                printf ("debug, exponents: %s, samples: %lld, time: %6.4fs, progress: +%d, rangefactor: %.9Le, precision_last:  %.3Le, precision: %.3Le, sm1_test:  %.3Le, sm2_test:  %.3Le, sm3_test: %.3Le, alpha: %.9Le, alpha_range: %.9Le, sm1: %.9Le, sm1_range: %.9Le, sm3: %.9Le, sm3_range: %.4Le, G: %.9Le, G_range: %.4Le, v: %.9Le, v_range: %.4Le, sm2: %.9Le, sm2_range: %.4Le, mz: %.9Le, mz_range: %.4Le, mw: %.9Le, mw_range: %.4Le, sin2w: %.9Le, sin2w_range: %.4Le, mh0: %.9Le, mh0_range: %.4Le, muser: %.9Le, muser_range: %.9Le\n", nle_state->exponents_str, samples, elapsed_time, progress, rangefactor, precision_last, precision, sm1_test, sm2_test, sm3_test, alpha, alpha_range, sm1, sm1_range, sm3, sm3_range, G, G_range, v, v_range, sm2, sm2_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range, muser, muser_range);
                                fflush(stdout);
#endif
                              } // end if  precision < precision_last
                            } // end if sm3_test/sm2_test/sm1_test
                          } // end if sm1_test / sm2_test
                        } // end for samples 

#ifdef DEBUG21
                        clock_gettime(CLOCK_REALTIME, &end_time);
                        elapsed_time=((long double)(end_time.tv_sec - 1500000000) + ((long double)end_time.tv_nsec / 1.0E9)) - ((long double)(start_time2.tv_sec - 1500000000) + ((long double)start_time2.tv_nsec) / 1.0E9);
                        printf("debug, Finished phase 2 samples loop, exponents: %s, samples: %lld, mass mode: %d%d%d, precision: %.6Le (%6.4fs)\n", nle_state->exponents_str, samples, nle_state->term1.current_match->smrfactor_mass_id, nle_state->term2.current_match->smrfactor_mass_id, nle_state->term3.current_match->smrfactor_mass_id, precision_last, elapsed_time);
                        fflush(stdout);
#endif
                        // determine output values/ranges
                        alpha_out=alpha_last;
                        if (alpha_out < alpha_out_low) {
                        alpha_out_low=alpha_out;
                        }
                        if (alpha_out > alpha_out_high) {
                          alpha_out_high=alpha_out;
                        }
                        alpha_out_c=((alpha_out_high + alpha_out_low) / 2.0l);
                        alpha_out_error=(alpha_out_high - alpha_out_c);
                        alpha_out_relerror=alpha_out_error / alpha_out_c;
                        alpha_out_diff=alpha_out_c - (long double)nle_config->ref_alpha_em;
                        alpha_out_reldiff=alpha_out_diff / (long double)nle_config->ref_alpha_em;

                        sm1_out=sm1_last;
                        if (sm1_out < sm1_out_low) {
                          sm1_out_low=sm1_out;
                        }
                        if (sm1_out > sm1_out_high) {
                          sm1_out_high=sm1_out;
                        }
                        sm1_out_c=((sm1_out_high + sm1_out_low) / 2.0l);
                        sm1_out_error=(sm1_out_high - sm1_out_c);
                        sm1_out_relerror=sm1_out_error / sm1_out_c;
                        sm1_out_diff=sm1_out_c - (long double)nle_config->ref_sm1;
                        sm1_out_reldiff=sm1_out_diff / (long double)nle_config->ref_sm1;

                        sm2_out=sm2_last;
                        if (sm2_out < sm2_out_low) {
                          sm2_out_low=sm2_out;
                        }
                        if (sm2_out > sm2_out_high) {
                          sm2_out_high=sm2_out;
                        }
                        sm2_out_c=((sm2_out_high + sm2_out_low) / 2.0l);
                        sm2_out_error=(sm2_out_high - sm2_out_c);
                        sm2_out_relerror=sm2_out_error / sm2_out_c;
                        sm2_out_diff=sm2_out_c - (long double)nle_config->ref_sm2;
                        sm2_out_reldiff=sm2_out_diff / (long double)nle_config->ref_sm2;

                        sm3_out=sm3_last;
                        if (sm3_out < sm3_out_low) {
                          sm3_out_low=sm3_out;
                        }
                        if (sm3_out > sm3_out_high) {
                          sm3_out_high=sm3_out;
                        }
                        sm3_out_c=((sm3_out_high + sm3_out_low) / 2.0l);
                        sm3_out_error=(sm3_out_high - sm3_out_c);
                        sm3_out_relerror=sm3_out_error / sm3_out_c;
                        sm3_out_diff=sm3_out_c - (long double)nle_config->ref_sm3;
                        sm3_out_reldiff=sm3_out_diff / (long double)nle_config->ref_sm3;

                        v_out=v_last;
                        if (v_out < v_out_low) {
                          v_out_low=v_out;
                        }
                        if (v_out > v_out_high) {
                          v_out_high=v_out;
                        }
                        v_out_c=((v_out_high + v_out_low) / 2.0l);
                        v_out_error=(v_out_high - v_out_c);
                        v_out_relerror=v_out_error / v_out_c;
                        v_out_diff=v_out_c - (long double)nle_config->ref_v;
                        v_out_reldiff=v_out_diff / (long double)nle_config->ref_v;

                        G_out=G_last;
                        if (G_out < G_out_low) {
                          G_out_low=G_out;
                        }
                        if (G_out > G_out_high) {
                          G_out_high=G_out;
                        }
                        G_out_c=((G_out_high + G_out_low) / 2.0l);
                        G_out_error=(G_out_high - G_out_c);
                        G_out_relerror=G_out_error / G_out_c;
                        G_out_diff=G_out_c - (long double)nle_config->ref_G;
                        G_out_reldiff=G_out_diff / (long double)nle_config->ref_G;

                        mz_out=mz_last;
                        if (mz_out < mz_out_low) {
                          mz_out_low=mz_out;
                        }
                        if (mz_out > mz_out_high) {
                          mz_out_high=mz_out;
                        }
                        mz_out_c=((mz_out_high + mz_out_low) / 2.0l);
                        mz_out_error=(mz_out_high - mz_out_c);
                        mz_out_relerror=mz_out_error / mz_out_c;
                        mz_out_diff=mz_out_c - (long double)nle_config->ref_mz;
                        mz_out_reldiff=mz_out_diff / (long double)nle_config->ref_mz;

                        mw_out=mw_last;
                        if (mw_out < mw_out_low) {
                          mw_out_low=mw_out;
                        }
                        if (mw_out > mw_out_high) {
                          mw_out_high=mw_out;
                        }
                        mw_out_c=((mw_out_high + mw_out_low) / 2.0l);
                        mw_out_error=(mw_out_high - mw_out_c);
                        mw_out_relerror=mw_out_error / mw_out_c;
                        mw_out_diff=mw_out_c - (long double)nle_config->ref_mw;
                        mw_out_reldiff=mw_out_diff / (long double)nle_config->ref_mw;

                        sin2w_out=sin2w_last;
                        if (sin2w_out < sin2w_out_low) {
                          sin2w_out_low=sin2w_out;
                        }
                        if (sin2w_out > sin2w_out_high) {
                          sin2w_out_high=sin2w_out;
                        }
                        sin2w_out_c=((sin2w_out_high + sin2w_out_low) / 2.0l);
                        sin2w_out_error=(sin2w_out_high - sin2w_out_c);
                        sin2w_out_relerror=sin2w_out_error / sin2w_out_c;
                        sin2w_out_diff=sin2w_out_c - (long double)nle_config->ref_sin2w;
                        sin2w_out_reldiff=sin2w_out_diff / (long double)nle_config->ref_sin2w;

                        mh0_out=mh0_last;
                        if (mh0_out < mh0_out_low) {
                          mh0_out_low=mh0_out;
                        }
                        if (mh0_out > mh0_out_high) {
                          mh0_out_high=mh0_out;
                        }
                        mh0_out_c=((mh0_out_high + mh0_out_low) / 2.0l);
                        mh0_out_error=(mh0_out_high - mh0_out_c);
                        mh0_out_relerror=mh0_out_error / mh0_out_c;
                        mh0_out_diff=mh0_out_c - (long double)nle_config->ref_mh0;
                        mh0_out_reldiff=mh0_out_diff / (long double)nle_config->ref_mh0;

                        muser_out=muser_last;
                        if (muser_out < muser_out_low) {
                          muser_out_low=muser_out;
                        }
                        if (muser_out > muser_out_high) {
                          muser_out_high=muser_out;
                        }
                        muser_out_c=((muser_out_high + muser_out_low) / 2.0l);
                        muser_out_error=(muser_out_high - muser_out_c);
                        muser_out_relerror=muser_out_error / muser_out_c;
                        muser_out_diff=muser_out_c - (long double)nle_config->smrfactor_mass_user;
                        muser_out_reldiff=muser_out_diff / (long double)nle_config->smrfactor_mass_user;
                      } // end muser_seq
                    } // end mh0_seq
                  } // end sin2w_seq
                } // end mw_seq
              } // end sm3_seq
            } // end G_seq
          } // end mz_seq
        } // end v_seq
      } // end sm2_seq
    } // end sm1_seq
  } // end alpha_seq

  // verify outputs against experimental uncertainties * phase2_results_window
  results_window=(long double)nle_config->phase2_results_window;
  valid_result=1;
  if ((nle_state->all_uses.float_G == 1) && ((fabsl(G_out_reldiff) / fmaxl(G_out_relerror, (results_window * (long double)nle_config->ref_G_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_v == 1) && ((fabsl(v_out_reldiff) / fmaxl(v_out_relerror, (results_window * (long double)nle_config->ref_v_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_mz == 1) && ((fabsl(mz_out_reldiff) / fmaxl(mz_out_relerror, (results_window * (long double)nle_config->ref_mz_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_mw == 1) && ((fabsl(mw_out_reldiff) / fmaxl(mw_out_relerror, (results_window * (long double)nle_config->ref_mw_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_mh0 == 1) && ((fabsl(mh0_out_reldiff) / fmaxl(mh0_out_relerror, (results_window * (long double)nle_config->ref_mh0_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_muser == 1) && ((fabsl(muser_out_reldiff) / fmaxl(muser_out_relerror, (results_window * (long double)nle_config->smrfactor_mass_user_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_sm1 == 1) && ((fabsl(sm1_out_reldiff) / fmaxl(sm1_out_relerror, (results_window * (long double)nle_config->ref_sm1_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_sm2 == 1) && ((fabsl(sm2_out_reldiff) / fmaxl(sm2_out_relerror, (results_window * (long double)nle_config->ref_sm2_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_sm3 == 1) && ((fabsl(sm3_out_reldiff) / fmaxl(sm3_out_relerror, (results_window * (long double)nle_config->ref_sm3_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_sin2w == 1) && ((fabsl(sin2w_out_reldiff) / fmaxl(sin2w_out_relerror, (results_window * (long double)nle_config->ref_sin2w_relerror))) > 1.0l)) {
    valid_result=0;
  }
  if ((nle_state->all_uses.float_alpha_em == 1) && ((fabsl(alpha_out_reldiff) / fmaxl(alpha_out_relerror, (results_window * (long double)nle_config->ref_alpha_em_relerror))) > 1.0l)) {
    valid_result=0;
  }

  if ((valid_result == 1) || (nle_config->phase2_results_always == 1)) {
    complexity=nle_state->term1.current_match->match_complexity + nle_state->term2.current_match->match_complexity + nle_state->term3.current_match->match_complexity;
    symmetry=nle_state->current_symmetry;
    //result_hash=lrand48();
    result_hash=(nle_state->term1.current_match->match_hash ^ nle_state->term2.current_match->match_hash) ^ nle_state->term3.current_match->match_hash;

    combined_score = (float)complexity / (float)symmetry;

    if (nle_config->nle_mode == 2) {
      sprintf(mass_str, "M%d%d", nle_state->term1.current_match->smrfactor_mass_id, nle_state->term2.current_match->smrfactor_mass_id);
    } else if (nle_config->nle_mode == 3) {
      sprintf(mass_str, "M%d%d%d", nle_state->term1.current_match->smrfactor_mass_id, nle_state->term2.current_match->smrfactor_mass_id, nle_state->term3.current_match->smrfactor_mass_id);
    } else {
      mass_str[0]=0;
    }

    sprintf(out_str_01, "result, %.4f, %3d, %3d, %s, %s, %12lld, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash);
    printf("%s\n", out_str_01);
    sprintf(out_str_02, "result, %.4f, %3d, %3d, %s, %s, %12lld, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash);
    printf("%s\n", out_str_02);
    sprintf(out_str_03, "result, %.4f, %3d, %3d, %s, %s, %12lld, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash);
    printf("%s\n", out_str_03);
    if (nle_state->all_uses.alpha_em == 1) {
      if (nle_state->all_uses.float_alpha_em == 1 ) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_04, "result, %.4f, %3d, %3d, %s, %s, %12lld, 04, | alpha_em   || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, alpha_out_c, alpha_out_error, alpha_out_relerror, nle_config->ref_alpha_em, nle_config->ref_alpha_em_error, nle_config->ref_alpha_em_relerror, alpha_out_diff, alpha_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_04);
    } else {
      out_str_04[0]=0;
    }
    if (nle_state->all_uses.v == 1) {
      if (nle_state->all_uses.float_v == 1 ) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_05, "result, %.4f, %3d, %3d, %s, %s, %12lld, 05, | v          || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, v_out_c, v_out_error, v_out_relerror, nle_config->ref_v, nle_config->ref_v_error, nle_config->ref_v_relerror, v_out_diff, v_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_05);
    } else {
      out_str_05[0]=0;
    }
    if ((nle_state->all_uses.mz == 1) || (nle_state->all_uses.mw_mz_mode == 1)) {
      if (nle_state->all_uses.float_mz == 1 ) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_06, "result, %.4f, %3d, %3d, %s, %s, %12lld, 06, | mZ         || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, mz_out_c, mz_out_error, mz_out_relerror, nle_config->ref_mz, nle_config->ref_mz_error, nle_config->ref_mz_relerror, mz_out_diff, mz_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_06);
    } else {
      out_str_06[0]=0;
    }
    if (nle_state->all_uses.G == 1) {
      if (nle_state->all_uses.float_G == 1 ) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_07, "result, %.4f, %3d, %3d, %s, %s, %12lld, 07, | G          || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, G_out_c, G_out_error, G_out_relerror, nle_config->ref_G, nle_config->ref_G_error, nle_config->ref_G_relerror, G_out_diff, G_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_07);
    } else {
      out_str_07[0]=0;
    }
    if ((nle_state->all_uses.mw == 1) || (nle_state->all_uses.mw_mz_mode)) {
      if (nle_state->all_uses.float_mw == 1 ) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_08, "result, %.4f, %3d, %3d, %s, %s, %12lld, 08, | mW         || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, mw_out_c, mw_out_error, mw_out_relerror, nle_config->ref_mw, nle_config->ref_mw_error, nle_config->ref_mw_relerror, mw_out_diff, mw_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_08);
    } else {
      out_str_08[0]=0;
    }
    if (nle_state->all_uses.sin2w == 1) {
      if (nle_state->all_uses.float_sin2w == 1) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else if (nle_state->all_uses.mw_mz_mode == 1) {
        sprintf(used_as_output, "D");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_09, "result, %.4f, %3d, %3d, %s, %s, %12lld, 09, | sin2w      || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, sin2w_out_c, sin2w_out_error, sin2w_out_relerror, nle_config->ref_sin2w, nle_config->ref_sin2w_error, nle_config->ref_sin2w_relerror, sin2w_out_diff, sin2w_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_09);
    } else {
      out_str_09[0]=0;
    }
    if (nle_state->all_uses.mh0 == 1) {
      if (nle_state->all_uses.float_mh0 == 1 ) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_10, "result, %.4f, %3d, %3d, %s, %s, %12lld, 10, | mH0        || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, mh0_out_c, mh0_out_error, mh0_out_relerror, nle_config->ref_mh0, nle_config->ref_mh0_error, nle_config->ref_mh0_relerror, mh0_out_diff, mh0_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_10);
    } else {
      out_str_10[0]=0;
    }
    if (nle_state->all_uses.m_user == 1) {
      if (nle_state->all_uses.float_muser == 1 ) {
        sprintf(used_as_output, "*");
        sprintf(used_as_input, " ");
      } else {
        sprintf(used_as_input, "*");
        sprintf(used_as_output, " ");
      }
      sprintf(out_str_11, "result, %.4f, %3d, %3d, %s, %s, %12lld, 11, | m_user     || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, muser_out_c, muser_out_error, muser_out_relerror, nle_config->smrfactor_mass_user, nle_config->smrfactor_mass_user_error, nle_config->smrfactor_mass_user_relerror, muser_out_diff, muser_out_reldiff, used_as_input, used_as_output);
      printf("%s\n", out_str_11);
    } else {
      out_str_11[0]=0;
    }
    if (nle_state->all_uses.float_sm1 == 1 ) {
      sprintf(used_as_output, "*");
      sprintf(used_as_input, " ");
    } else {
      sprintf(used_as_input, "*");
      sprintf(used_as_output, " ");
    }
    sprintf(out_str_12, "result, %.4f, %3d, %3d, %s, %s, %12lld, 12, | sm1        || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, sm1_out_c, sm1_out_error, sm1_out_relerror, nle_config->ref_sm1, nle_config->ref_sm1_error, nle_config->ref_sm1_relerror, sm1_out_diff, sm1_out_reldiff, used_as_input, used_as_output);
    printf("%s\n", out_str_12);
    if (nle_state->all_uses.float_sm2 == 1) {
      sprintf(used_as_input, " ");
      sprintf(used_as_output, "*");
    } else {
      sprintf(used_as_input, "*");
      sprintf(used_as_output, " ");
    }
    sprintf(out_str_13, "result, %.4f, %3d, %3d, %s, %s, %12lld, 13, | sm2        || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, sm2_out_c, sm2_out_error, sm2_out_relerror, nle_config->ref_sm2, nle_config->ref_sm2_error, nle_config->ref_sm2_relerror, sm2_out_diff, sm2_out_reldiff, used_as_input, used_as_output);
    printf("%s\n", out_str_13);
    if (nle_state->all_uses.float_sm3 == 1) {
      sprintf(used_as_input, " ");
      sprintf(used_as_output, "*");
    } else {
      sprintf(used_as_input, "*");
      sprintf(used_as_output, " ");
    }
    sprintf(out_str_14, "result, %.4f, %3d, %3d, %s, %s, %12lld, 14, | sm3        || %.15Le | %.3Le | %.3Le || %.15e | %.3e | %.3e || %11.4Le | %11.4Le |       %s       |       %s        |", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, sm3_out_c, sm3_out_error, sm3_out_relerror, nle_config->ref_sm3, nle_config->ref_sm3_error, nle_config->ref_sm3_relerror, sm3_out_diff, sm3_out_reldiff, used_as_input, used_as_output);
    printf("%s\n", out_str_14);
    sprintf(out_str_15, "result, %.4f, %3d, %3d, %s, %s, %12lld, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash);
    printf("%s\n", out_str_15);

    // if outfactor_user* appears in any term display it's value on next line
    if ((nle_state->term1.current_match->outfactor_user1_exp_up != 0) || (nle_state->term2.current_match->outfactor_user1_exp_up != 0) || (nle_state->term3.current_match->outfactor_user1_exp_up != 0)) {
      sprintf(user1_out_str, "uout1=%.9e,", nle_config->outfactor_user1);
    } else {
      sprintf(user1_out_str, "                     ");
    }
    if ((nle_state->term1.current_match->outfactor_user2_exp_up != 0) || (nle_state->term2.current_match->outfactor_user2_exp_up != 0) || (nle_state->term3.current_match->outfactor_user2_exp_up != 0)) {
      sprintf(user2_out_str, "uout2=%.9e,", nle_config->outfactor_user2);
    } else {
      sprintf(user2_out_str, "                     ");
    }
    if ((nle_state->term1.current_match->outfactor_user3_exp_up != 0) || (nle_state->term2.current_match->outfactor_user3_exp_up != 0) || (nle_state->term3.current_match->outfactor_user3_exp_up != 0)) {
      sprintf(user3_out_str, "uout3=%.9e,", nle_config->outfactor_user3);
    } else {
      sprintf(user3_out_str, "                     ");
    }

    // if infactor_user appears in any term display it's value on next line
    if ((nle_state->term1.current_match->infactor_user_exp_up != 0) || (nle_state->term2.current_match->infactor_user_exp_up != 0) || (nle_state->term3.current_match->infactor_user_exp_up != 0)) {
      sprintf(user_in_str, "uin=%.9e,", nle_config->infactor_user);
    } else {
      sprintf(user_in_str, "                    ");
    }

    if (nle_config->nle_mode == 2) {
      if (nle_state->nle_mixing_polarity == 0) {
        sprintf(out_str_16, "result, %.4f, %3d, %3d, %s, %s, %12lld, 16, formula: term1^2 - (term3 * term1 * term2) + term2^2 - 1 = 0, %s %s %s %s", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, user1_out_str, user2_out_str, user3_out_str, user_in_str);
      } else {
        sprintf(out_str_16, "result, %.4f, %3d, %3d, %s, %s, %12lld, 16, formula: term1^2 + (term3 * term1 * term2) + term2^2 - 1 = 0, %s %s %s %s", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, user1_out_str, user2_out_str, user3_out_str, user_in_str);
      }
    } else if (nle_config->nle_mode == 3) {
      sprintf(out_str_16, "result, %.4f, %3d, %3d, %s, %s, %12lld, 16, formula: term1 - term2 + term3 - 1 = 0, %s %s %s %s", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, user1_out_str, user2_out_str, user3_out_str, user_in_str);
    }
    printf("%s\n", out_str_16);
    sprintf(out_str_17, "result, %.4f, %3d, %3d, %s, %s, %12lld, 17, term1=%s", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, term1_formula_str);
    printf("%s\n", out_str_17);
    sprintf(out_str_18, "result, %.4f, %3d, %3d, %s, %s, %12lld, 18, term2=%s", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, term2_formula_str);
    printf("%s\n", out_str_18);
    sprintf(out_str_19, "result, %.4f, %3d, %3d, %s, %s, %12lld, 19, term3=%s", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, term3_formula_str);
    printf("%s\n", out_str_19);
    if (nle_config->smrfactor_1minus_enable == 1) {
      sprintf(out_str_20, "result, %.4f, %3d, %3d, %s, %s, %12lld, 20, smrf= %s", combined_score, symmetry, complexity, nle_state->exponents_str, mass_str, result_hash, smrf_str);
      printf("%s\n", out_str_20);
    }
    fflush(stdout);
    if (nle_config->upload_results_enable == 1) {
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_01, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_02, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_03, 511));
      system(exec_str);
      if (out_str_04[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_04, 511));
        system(exec_str);
      }
      if (out_str_05[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_05, 511));
        system(exec_str);
      }
      if (out_str_06[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_06, 511));
        system(exec_str);
      }
      if (out_str_07[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_07, 511));
        system(exec_str);
      }
      if (out_str_08[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_08, 511));
        system(exec_str);
      }
      if (out_str_09[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_09, 511));
        system(exec_str);
      }
      if (out_str_10[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_10, 511));
        system(exec_str);
      }
      if (out_str_11[0] != 0) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_11, 511));
        system(exec_str);
      }
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_12, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_13, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_14, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_15, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_16, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_17, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_18, 511));
      system(exec_str);
      sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_19, 511));
      system(exec_str);
      if (nle_config->smrfactor_1minus_enable == 1) {
        sprintf(exec_str, "curl -s \"%s/%s\" > /dev/null 2>&1\n", nle_config->upload_url, underscore(out_str_20, 511));
        system(exec_str);
      }
    } // end if upload_results_enable
  } // end if score
  return(precision_last);
}
