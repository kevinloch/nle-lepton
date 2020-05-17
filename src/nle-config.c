#define _GNU_SOURCE // needed for strcasestr in string.h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nle-lepton.h"

void initConfig(nle_config_t *nle_config) {
  // sets default values for configuration options.
  // these can be overridden by nle-lepton.cfg
  nle_config->external_seed=0;
  strcpy(nle_config->config_file_name, "./nle-lepton.cfg");
  nle_config->upload_results_enable=0;
  strcpy(nle_config->upload_url, "http://localhost/lepton");
  nle_config->phase1_status_enable=1;
  nle_config->phase1_solution_detail=1;
  nle_config->phase2_status_enable=1;
  nle_config->ref_sm1=0.51099895000E6;
  nle_config->ref_sm1_error=0.00000000015E6;
  nle_config->ref_sm2=105.6583755E6;
  nle_config->ref_sm2_error=0.0000023E6;
  nle_config->ref_sm3=1776.86E6;
  nle_config->ref_sm3_error=0.12E6;
  nle_config->nle_mode=3;
  nle_config->nle_mixing_polarity=0;
  nle_config->exp_inv_max=9;
  nle_config->exp_inv_term1_force=0;
  nle_config->exp_inv_term2_force=0;
  nle_config->exp_inv_term3_force=0;
  nle_config->exp_inv_term4_force=0;
  nle_config->exp_inv_2seq_limit=12;
  nle_config->exp_inv_3seq_limit=9;
  nle_config->exp_inv_4seq_limit=6;
  nle_config->exp_neg_enable=1;
  nle_config->phase1_random_samples_enable=1;
  nle_config->phase1_filter=5;
  nle_config->phase1_int_match_max=16;
  nle_config->phase1_int_match_filter=1;
  nle_config->phase1_two_term_test_min=1.75;
  nle_config->phase1_two_term_test_max=4.25;
  nle_config->phase1_unsolvable_checkpoint=100000;
  nle_config->phase1_mc_samples_limit=100000000;
  nle_config->phase2_enable=1;
  nle_config->phase2_symmetry_min=65;
  nle_config->phase2_complexity_max=75;
  nle_config->phase2_check_nbv_nss=1;
  nle_config->phase2_check_weak=1;
  nle_config->phase2_results_window=1.1;
  nle_config->phase2_results_always=0;
  nle_config->phase2_ignore_small_rel_unc=0;
  nle_config->phase2_check_rmr=1;
  nle_config->smrfactor_mass_mp_enable=1;
  nle_config->smrfactor_mass_v_enable=1;
  nle_config->smrfactor_mass_mz_enable=1;
  nle_config->smrfactor_mass_mw_enable=1;
  nle_config->smrfactor_mass_mh0_enable=1;
  nle_config->smrfactor_mass_user_enable=0;
  nle_config->smrfactor_mass_user=1.0E9;
  nle_config->smrfactor_mass_user_error=0.0;
  nle_config->smrfactor_1minus_enable=0;
  nle_config->smrfactor_1minus_random=0;
  nle_config->smrfactor_rational_max=32;
  nle_config->smrfactor_rational_filter=1;
  nle_config->smrfactor_2_exp_up_max=1;
  nle_config->smrfactor_2_exp_down_max=2;
  nle_config->smrfactor_alpha_exp_up_max=2;
  nle_config->smrfactor_alpha_exp_down_max=2;
  nle_config->smrfactor_pi_exp_up_max=2;
  nle_config->smrfactor_pi_exp_down_max=2;
  nle_config->smrfactor_user=1.0;
  nle_config->smrfactor_user_exp_up_max=0;
  nle_config->smrfactor_user_exp_down_max=1;
  nle_config->rmrfactor_smrfactor_sync=1;
  nle_config->rmrfactor_mass_mp_enable=1;
  nle_config->rmrfactor_mass_v_enable=1;
  nle_config->rmrfactor_mass_mz_enable=1;
  nle_config->rmrfactor_mass_mw_enable=1;
  nle_config->rmrfactor_mass_mh0_enable=1;
  nle_config->rmrfactor_mass_user_enable=0;
  nle_config->rmrfactor_mass_user=1.0E9;
  nle_config->rmrfactor_mass_user_error=0.0;
  nle_config->rmrfactor_1minus_enable=0;
  nle_config->rmrfactor_1minus_random=0;
  nle_config->rmrfactor_rational_max=32;
  nle_config->rmrfactor_rational_filter=1;
  nle_config->rmrfactor_2_exp_up_max=1;
  nle_config->rmrfactor_2_exp_down_max=2;
  nle_config->rmrfactor_alpha_exp_up_max=2;
  nle_config->rmrfactor_alpha_exp_down_max=2;
  nle_config->rmrfactor_pi_exp_up_max=2;
  nle_config->rmrfactor_pi_exp_down_max=2;
  nle_config->rmrfactor_user=1.0;
  nle_config->rmrfactor_user_exp_up_max=0;
  nle_config->rmrfactor_user_exp_down_max=1;
  nle_config->infactor_rational_max=16;
  nle_config->infactor_rational_filter=1;
  nle_config->infactor_2_exp_up_max=1;
  nle_config->infactor_2_exp_down_max=2;
  nle_config->infactor_alpha_exp_up_max=1;
  nle_config->infactor_alpha_exp_down_max=1;
  nle_config->infactor_pi_exp_up_max=3;
  nle_config->infactor_pi_exp_down_max=1;
  nle_config->infactor_nss_enable=1;
  nle_config->infactor_nbv_enable=1;
  nle_config->infactor_user=1.0;
  nle_config->infactor_user_exp_up_max=0;
  nle_config->infactor_user_exp_down_max=1;
  nle_config->outfactor_rational_max=3;
  nle_config->outfactor_rational_filter=1;
  nle_config->outfactor_2_exp_up_max=1;
  nle_config->outfactor_2_exp_down_max=2;
  nle_config->outfactor_alpha_exp_up_max=2;
  nle_config->outfactor_alpha_exp_down_max=2;
  nle_config->outfactor_pi_exp_up_max=2;
  nle_config->outfactor_pi_exp_down_max=2;
  nle_config->outfactor_weak_exp_up_max=0;
  nle_config->outfactor_weak_exp_down_max=1;
  nle_config->outfactor_user1=1.0;
  nle_config->outfactor_user1_exp_up_max=0;
  nle_config->outfactor_user1_exp_down_max=1;
  nle_config->outfactor_user2=1.0;
  nle_config->outfactor_user2_exp_up_max=0;
  nle_config->outfactor_user2_exp_down_max=1;
  nle_config->outfactor_user3=1.0;
  nle_config->outfactor_user3_exp_up_max=0;
  nle_config->outfactor_user3_exp_down_max=1;
  nle_config->outfactor_rmr_exp_up_max=0;
  nle_config->outfactor_rmr_exp_down_max=1;
  nle_config->outfactor_rmr_mp_enable=0;
  nle_config->outfactor_rmr_v_enable=0;
  nle_config->outfactor_rmr_mz_enable=0;
  nle_config->outfactor_rmr_mw_enable=0;
  nle_config->outfactor_rmr_mh0_enable=0;
  nle_config->outfactor_rmr_user_enable=0;
  nle_config->ref_c=2.997924580000E+08;
  nle_config->ref_h=6.62607015E-34;
  nle_config->ref_hbar=1.05457181764616E-34;
  nle_config->ref_e=1.602176634E-19;
  nle_config->ref_ev_to_kg=1.78266192162790E-36;
  nle_config->ref_kg_to_ev=5.60958860380445E35;
  nle_config->ref_alpha_em=7.2973525693E-3;
  nle_config->ref_alpha_em_error=0.0000000011E-3;
  nle_config->ref_alpha_w=3.0E-7;
  nle_config->ref_alpha_w_error=2.0E-7;
  nle_config->ref_v=246.219651E9;
  nle_config->ref_v_error=0.000063E9;
  nle_config->ref_mz=91.1876E9;
  nle_config->ref_mz_error=0.0021E9;
  nle_config->ref_mw=80.379E9;
  nle_config->ref_mw_error=0.012E9;
  nle_config->ref_mh0=125.35E9;
  nle_config->ref_mh0_error=0.15E9;
  nle_config->ref_sin2w=0.22311;
  nle_config->ref_sin2w_error=0.00027;
  nle_config->ref_G=6.67430E-11;
  nle_config->ref_G_error=0.00015E-11;
}

void checkOptionBool(int *config_int, char *option, char *value, char *matchstr) {
  size_t matchstr_length;

  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) != NULL) && (option[matchstr_length] != '_')) {
    //printf("Boolean option: %s, value: %s\n", option, value);
    if (strcasestr(value, "yes") != NULL) {
      *config_int=1;
    } else {
      *config_int=0;
    }
    //printf("set value: %d\n", *config_int);
    //fflush(stdout);
  }
}

void checkOptionInt(int *config_int, char *option, char *value, char *matchstr) {
  size_t matchstr_length;
  
  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) != NULL) && (option[matchstr_length] != '_')) {
    //printf("Integer option: %s, value: %s\n", option, value);
    *config_int=strtol(value, NULL, 10);
    //printf("set value: %d\n", *config_int);
    //fflush(stdout);
  }
}

void checkOptionDouble(double *config_double, char *option, char *value, char *matchstr) {
  size_t matchstr_length;
  
  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) != NULL) && (option[matchstr_length] != '_')) {
    //printf("Double option: %s, value: %s\n", option, value);
    *config_double=strtod(value, NULL);
    //printf("set value: %.16e\n", *config_double);
    //fflush(stdout);
  }
}

void checkOptionStr(char *config_str,  char *option, char *value, char *matchstr) {
  size_t matchstr_length;
  
  matchstr_length=strlen(matchstr);
  if ((strstr(option, matchstr) != NULL) && (option[matchstr_length] != '_')) {
    //printf("String option: %s, value: %s\n", option, value);
    strcpy(config_str, value);
    //printf("set value: %s\n", config_str);
    //fflush(stdout);
  }
}

void setOptionValue(nle_config_t *nle_config, char *option, char *value) {
  checkOptionBool(&nle_config->upload_results_enable, option, value, "upload_results_enable");
  checkOptionStr(nle_config->upload_url, option, value, "upload_url");
  checkOptionBool(&nle_config->phase1_status_enable, option, value, "phase1_status_enable");
  checkOptionBool(&nle_config->phase1_solution_detail, option, value, "phase1_solution_detail");
  checkOptionBool(&nle_config->phase2_status_enable, option, value, "phase2_status_enable");
  checkOptionDouble(&nle_config->ref_sm1, option, value, "ref_sm1");
  checkOptionDouble(&nle_config->ref_sm1_error, option, value, "ref_sm1_error");
  checkOptionDouble(&nle_config->ref_sm2, option, value, "ref_sm2");
  checkOptionDouble(&nle_config->ref_sm2_error, option, value, "ref_sm2_error");
  checkOptionDouble(&nle_config->ref_sm3, option, value, "ref_sm3");
  checkOptionDouble(&nle_config->ref_sm3_error, option, value, "ref_sm3_error");
  checkOptionInt(&nle_config->nle_mode, option, value, "nle_mode");
  checkOptionInt(&nle_config->nle_mixing_polarity, option, value, "nle_mixing_polarity");
  checkOptionInt(&nle_config->exp_inv_max, option, value, "exp_inv_max");
  checkOptionInt(&nle_config->exp_inv_term1_force, option, value, "exp_inv_term1_force");
  checkOptionInt(&nle_config->exp_inv_term2_force, option, value, "exp_inv_term2_force");
  checkOptionInt(&nle_config->exp_inv_term3_force, option, value, "exp_inv_term3_force");
  checkOptionInt(&nle_config->exp_inv_term4_force, option, value, "exp_inv_term4_force");
  checkOptionInt(&nle_config->exp_inv_2seq_limit, option, value, "exp_inv_2seq_limit");
  checkOptionInt(&nle_config->exp_inv_3seq_limit, option, value, "exp_inv_3seq_limit");
  checkOptionInt(&nle_config->exp_inv_4seq_limit, option, value, "exp_inv_4seq_limit");
  checkOptionBool(&nle_config->exp_neg_enable, option, value, "exp_neg_enable");
  checkOptionBool(&nle_config->phase1_random_samples_enable, option, value, "phase1_random_samples_enable");
  checkOptionInt(&nle_config->phase1_filter, option, value, "phase1_filter");
  checkOptionInt(&nle_config->phase1_int_match_max, option, value, "phase1_int_match_max");
  checkOptionBool(&nle_config->phase1_int_match_filter, option, value, "phase1_int_match_filter");
  checkOptionDouble(&nle_config->phase1_two_term_test_min, option, value, "phase1_two_term_test_min");
  checkOptionDouble(&nle_config->phase1_two_term_test_max, option, value, "phase1_two_term_test_max");
  checkOptionInt(&nle_config->phase1_unsolvable_checkpoint, option, value, "phase1_unsolvable_checkpoint");
  checkOptionInt(&nle_config->phase1_mc_samples_limit, option, value, "phase1_mc_samples_limit");
  checkOptionBool(&nle_config->phase2_enable, option, value, "phase2_enable");
  checkOptionInt(&nle_config->phase2_symmetry_min, option, value, "phase2_symmetry_min");
  checkOptionInt(&nle_config->phase2_complexity_max, option, value, "phase2_complexity_max");
  checkOptionBool(&nle_config->phase2_check_nbv_nss, option, value, "phase2_check_nbv_nss");
  checkOptionBool(&nle_config->phase2_check_weak, option, value, "phase2_check_weak");
  checkOptionDouble(&nle_config->phase2_results_window, option, value, "phase2_results_window");
  checkOptionBool(&nle_config->phase2_results_always, option, value, "phase2_results_always");
  checkOptionBool(&nle_config->phase2_ignore_small_rel_unc, option, value, "phase2_ignore_small_rel_unc");
  checkOptionBool(&nle_config->phase2_check_rmr, option, value, "phase2_check_rmr");
  checkOptionBool(&nle_config->smrfactor_mass_mp_enable, option, value, "smrfactor_mass_mp_enable");
  checkOptionBool(&nle_config->smrfactor_mass_v_enable, option, value, "smrfactor_mass_v_enable");
  checkOptionBool(&nle_config->smrfactor_mass_mz_enable, option, value, "smrfactor_mass_mz_enable");
  checkOptionBool(&nle_config->smrfactor_mass_mw_enable, option, value, "smrfactor_mass_mw_enable");
  checkOptionBool(&nle_config->smrfactor_mass_mh0_enable, option, value, "smrfactor_mass_mh0_enable");
  checkOptionBool(&nle_config->smrfactor_mass_user_enable, option, value, "smrfactor_mass_user_enable");
  checkOptionDouble(&nle_config->smrfactor_mass_user, option, value, "smrfactor_mass_user");
  checkOptionDouble(&nle_config->smrfactor_mass_user_error, option, value, "smrfactor_mass_user_error");
  checkOptionBool(&nle_config->smrfactor_1minus_enable, option, value, "smrfactor_1minus_enable");
  checkOptionBool(&nle_config->smrfactor_1minus_random, option, value, "smrfactor_1minus_random");
  checkOptionInt(&nle_config->smrfactor_rational_max, option, value, "smrfactor_rational_max");
  checkOptionBool(&nle_config->smrfactor_rational_filter, option, value, "smrfactor_rational_filter");
  checkOptionInt(&nle_config->smrfactor_2_exp_up_max, option, value, "smrfactor_2_exp_up_max");
  checkOptionInt(&nle_config->smrfactor_2_exp_down_max, option, value, "smrfactor_2_exp_down_max");
  checkOptionInt(&nle_config->smrfactor_alpha_exp_up_max, option, value, "smrfactor_alpha_exp_up_max");
  checkOptionInt(&nle_config->smrfactor_alpha_exp_down_max, option, value, "smrfactor_alpha_exp_down_max");
  checkOptionInt(&nle_config->smrfactor_pi_exp_up_max, option, value, "smrfactor_pi_exp_up_max");
  checkOptionInt(&nle_config->smrfactor_pi_exp_down_max, option, value, "smrfactor_pi_exp_down_max");
  checkOptionDouble(&nle_config->smrfactor_user, option, value, "smrfactor_user");
  checkOptionInt(&nle_config->smrfactor_user_exp_up_max, option, value, "smrfactor_user_exp_up_max");
  checkOptionInt(&nle_config->smrfactor_user_exp_down_max, option, value, "smrfactor_user_exp_down_max");
  checkOptionBool(&nle_config->rmrfactor_mass_mp_enable, option, value, "rmrfactor_mass_mp_enable");
  checkOptionBool(&nle_config->rmrfactor_mass_v_enable, option, value, "rmrfactor_mass_v_enable");
  checkOptionBool(&nle_config->rmrfactor_mass_mz_enable, option, value, "rmrfactor_mass_mz_enable");
  checkOptionBool(&nle_config->rmrfactor_mass_mw_enable, option, value, "rmrfactor_mass_mw_enable");
  checkOptionBool(&nle_config->rmrfactor_mass_mh0_enable, option, value, "rmrfactor_mass_mh0_enable");
  checkOptionBool(&nle_config->rmrfactor_mass_user_enable, option, value, "rmrfactor_mass_user_enable");
  checkOptionDouble(&nle_config->rmrfactor_mass_user, option, value, "rmrfactor_mass_user");
  checkOptionDouble(&nle_config->rmrfactor_mass_user_error, option, value, "rmrfactor_mass_user_error");
  checkOptionBool(&nle_config->rmrfactor_1minus_enable, option, value, "rmrfactor_1minus_enable");
  checkOptionBool(&nle_config->rmrfactor_1minus_random, option, value, "rmrfactor_1minus_random");
  checkOptionBool(&nle_config->rmrfactor_smrfactor_sync, option, value, "rmrfactor_smrfactor_sync");
  checkOptionInt(&nle_config->rmrfactor_rational_max, option, value, "rmrfactor_rational_max");
  checkOptionBool(&nle_config->rmrfactor_rational_filter, option, value, "rmrfactor_rational_filter");
  checkOptionInt(&nle_config->rmrfactor_2_exp_up_max, option, value, "rmrfactor_2_exp_up_max");
  checkOptionInt(&nle_config->rmrfactor_2_exp_down_max, option, value, "rmrfactor_2_exp_down_max");
  checkOptionInt(&nle_config->rmrfactor_alpha_exp_up_max, option, value, "rmrfactor_alpha_exp_up_max");
  checkOptionInt(&nle_config->rmrfactor_alpha_exp_down_max, option, value, "rmrfactor_alpha_exp_down_max");
  checkOptionInt(&nle_config->rmrfactor_pi_exp_up_max, option, value, "rmrfactor_pi_exp_up_max");
  checkOptionInt(&nle_config->rmrfactor_pi_exp_down_max, option, value, "rmrfactor_pi_exp_down_max");
  checkOptionDouble(&nle_config->rmrfactor_user, option, value, "rmrfactor_user");
  checkOptionInt(&nle_config->rmrfactor_user_exp_up_max, option, value, "rmrfactor_user_exp_up_max");
  checkOptionInt(&nle_config->rmrfactor_user_exp_down_max, option, value, "rmrfactor_user_exp_down_max");
  checkOptionInt(&nle_config->infactor_rational_max, option, value, "infactor_rational_max");
  checkOptionBool(&nle_config->infactor_rational_filter, option, value, "infactor_rational_filter");
  checkOptionInt(&nle_config->infactor_2_exp_up_max, option, value, "infactor_2_exp_up_max");
  checkOptionInt(&nle_config->infactor_2_exp_down_max, option, value, "infactor_2_exp_down_max");
  checkOptionInt(&nle_config->infactor_alpha_exp_up_max, option, value, "infactor_alpha_exp_up_max");
  checkOptionInt(&nle_config->infactor_alpha_exp_down_max, option, value, "infactor_alpha_exp_down_max");
  checkOptionInt(&nle_config->infactor_pi_exp_up_max, option, value, "nfactor_pi_exp_up_max");
  checkOptionInt(&nle_config->infactor_pi_exp_down_max, option, value, "infactor_pi_exp_down_max");
  checkOptionBool(&nle_config->infactor_nss_enable, option, value, "infactor_nss_enable");
  checkOptionBool(&nle_config->infactor_nbv_enable, option, value, "infactor_nbv_enable");
  checkOptionDouble(&nle_config->infactor_user, option, value, "infactor_user");
  checkOptionInt(&nle_config->infactor_user_exp_up_max, option, value, "infactor_user_exp_up_max");
  checkOptionInt(&nle_config->infactor_user_exp_down_max, option, value, "infactor_user_exp_down_max");
  checkOptionInt(&nle_config->outfactor_rational_max, option, value, "outfactor_rational_max");
  checkOptionBool(&nle_config->outfactor_rational_filter, option, value, "outfactor_rational_filter");
  checkOptionInt(&nle_config->outfactor_2_exp_up_max, option, value, "outfactor_2_exp_up_max");
  checkOptionInt(&nle_config->outfactor_2_exp_down_max, option, value, "outfactor_2_exp_down_max");
  checkOptionInt(&nle_config->outfactor_alpha_exp_up_max, option, value, "outfactor_alpha_exp_up_max");
  checkOptionInt(&nle_config->outfactor_alpha_exp_down_max, option, value, "outfactor_alpha_exp_down_max");
  checkOptionInt(&nle_config->outfactor_pi_exp_up_max, option, value, "outfactor_pi_exp_up_max");
  checkOptionInt(&nle_config->outfactor_pi_exp_down_max, option, value, "outfactor_pi_exp_down_max");
  checkOptionInt(&nle_config->outfactor_weak_exp_up_max, option, value, "outfactor_weak_exp_up_max");
  checkOptionInt(&nle_config->outfactor_weak_exp_down_max, option, value, "outfactor_weak_exp_down_max");
  checkOptionDouble(&nle_config->outfactor_user1, option, value, "outfactor_user1");
  checkOptionInt(&nle_config->outfactor_user1_exp_up_max, option, value, "outfactor_user1_exp_up_max");
  checkOptionInt(&nle_config->outfactor_user1_exp_down_max, option, value, "outfactor_user1_exp_down_max");
  checkOptionDouble(&nle_config->outfactor_user2, option, value, "outfactor_user2");
  checkOptionInt(&nle_config->outfactor_user2_exp_up_max, option, value, "outfactor_user2_exp_up_max");
  checkOptionInt(&nle_config->outfactor_user2_exp_down_max, option, value, "outfactor_user2_exp_down_max");
  checkOptionDouble(&nle_config->outfactor_user3, option, value, "outfactor_user3");
  checkOptionInt(&nle_config->outfactor_user3_exp_up_max, option, value, "outfactor_user3_exp_up_max");
  checkOptionInt(&nle_config->outfactor_user3_exp_down_max, option, value, "outfactor_user3_exp_down_max");
  checkOptionInt(&nle_config->outfactor_rmr_exp_up_max, option, value, "outfactor_rmr_exp_up_max");
  checkOptionInt(&nle_config->outfactor_rmr_exp_down_max, option, value, "outfactor_rmr_exp_down_max");
  checkOptionBool(&nle_config->outfactor_rmr_mp_enable, option, value, "outfactor_rmr_mp_enable");
  checkOptionBool(&nle_config->outfactor_rmr_v_enable, option, value, "outfactor_rmr_v_enable");
  checkOptionBool(&nle_config->outfactor_rmr_mz_enable, option, value, "outfactor_rmr_mz_enable");
  checkOptionBool(&nle_config->outfactor_rmr_mw_enable, option, value, "outfactor_rmr_mw_enable");
  checkOptionBool(&nle_config->outfactor_rmr_mh0_enable, option, value, "outfactor_rmr_mh0_enable");
  checkOptionBool(&nle_config->outfactor_rmr_user_enable, option, value, "outfactor_rmr_user_enable");
  checkOptionDouble(&nle_config->ref_c, option, value, "ref_c");
  checkOptionDouble(&nle_config->ref_h, option, value, "ref_h");
  checkOptionDouble(&nle_config->ref_hbar, option, value, "ref_hbar");
  checkOptionDouble(&nle_config->ref_e, option, value, "ref_e");
  checkOptionDouble(&nle_config->ref_ev_to_kg, option, value, "ref_ev_to_kg");
  checkOptionDouble(&nle_config->ref_kg_to_ev, option, value, "ref_kg_to_ev");
  checkOptionDouble(&nle_config->ref_alpha_em, option, value, "ref_alpha_em");
  checkOptionDouble(&nle_config->ref_alpha_em_error, option, value, "ref_alpha_em_error");
  checkOptionDouble(&nle_config->ref_alpha_w, option, value, "ref_alpha_w");
  checkOptionDouble(&nle_config->ref_alpha_w_error, option, value, "ref_alpha_w_error");
  checkOptionDouble(&nle_config->ref_v, option, value, "ref_v");
  checkOptionDouble(&nle_config->ref_v_error, option, value, "ref_v_error");
  checkOptionDouble(&nle_config->ref_mz, option, value, "ref_mz");
  checkOptionDouble(&nle_config->ref_mz_error, option, value, "ref_mz_error");
  checkOptionDouble(&nle_config->ref_mw, option, value, "ref_mw");
  checkOptionDouble(&nle_config->ref_mw_error, option, value, "ref_mw_error");
  checkOptionDouble(&nle_config->ref_mh0, option, value, "ref_mh0");
  checkOptionDouble(&nle_config->ref_mh0_error, option, value, "ref_mh0_error");
  checkOptionDouble(&nle_config->ref_sin2w, option, value, "ref_sin2w");
  checkOptionDouble(&nle_config->ref_sin2w_error, option, value, "ref_sin2w_error");
  checkOptionDouble(&nle_config->ref_G, option, value, "ref_G");
  checkOptionDouble(&nle_config->ref_G_error, option, value, "ref_G_error");
}

int loadConfig(nle_config_t *nle_config) {
  FILE *config_file;
  char *input_line_p;
  char *symbol_p;
  size_t input_line_length;
  char input_line[256];
  char input_line_trimmed[256];
  size_t input_line_trimmed_length;
  char option[256];
  size_t option_length;
  char value[256];
  size_t value_length;
  
  // attempt to open config file
  printf("init, Loading configuration from %s\n", nle_config->config_file_name);
  fflush(stdout);
  config_file=fopen(nle_config->config_file_name, "r");
  if (config_file == NULL) {
    printf("init, Error: could not open %s\n", nle_config->config_file_name);
    fflush(stdout);
    return(1);
  }

  // read and process each line of config file
  input_line_p=fgets(input_line, 256, config_file);
  while(input_line_p != NULL) {
    input_line_length=strlen(input_line);

    // search for comment symbol and remove any comments, otherwise just remove newline
    symbol_p=strchr(input_line, '#');
    if (symbol_p != NULL) {
      input_line_trimmed_length=(symbol_p - input_line);
    } else {
      input_line_trimmed_length=input_line_length-1;
    } 
    strncpy(input_line_trimmed, input_line, input_line_trimmed_length);
    input_line_trimmed[input_line_trimmed_length]=0;

    // search for option/value delimiter and split option and value strings
    symbol_p=strchr(input_line_trimmed, '=');
    if (symbol_p != NULL) {
      option_length=(symbol_p - input_line_trimmed);
      value_length=(input_line_trimmed_length - option_length);
      if ((option_length <= 254) && (value_length <= 254)) {
        strncpy(option, input_line_trimmed, option_length);
        option[option_length]=0;
        strncpy(value, (symbol_p+=1), value_length);
        value[value_length]=0;

        // send to option value processing fucntion
        setOptionValue(nle_config, option, value);

      } // end option_length and value_length checks
    } // end symbol_p check
    input_line_p=fgets(input_line, 256, config_file);
  } // end while input_line_raw

  // set relative uncertainties
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
  nle_config->ref_G_relerror=nle_config->ref_G_error / nle_config->ref_G;
  nle_config->relerror[0]=nle_config->ref_G_relerror;
  nle_config->ref_v_relerror=nle_config->ref_v_error / nle_config->ref_v;
  nle_config->relerror[1]=nle_config->ref_v_relerror;
  nle_config->ref_mz_relerror=nle_config->ref_mz_error / nle_config->ref_mz;
  nle_config->relerror[2]=nle_config->ref_mz_relerror;
  nle_config->ref_mw_relerror=nle_config->ref_mw_error / nle_config->ref_mw;
  nle_config->relerror[3]=nle_config->ref_mw_relerror;
  nle_config->ref_mh0_relerror=nle_config->ref_mh0_error / nle_config->ref_mh0;
  nle_config->relerror[4]=nle_config->ref_mh0_relerror;
  nle_config->smrfactor_mass_user_relerror=nle_config->smrfactor_mass_user_error / nle_config->smrfactor_mass_user;
  nle_config->relerror[5]=nle_config->smrfactor_mass_user_relerror;
  nle_config->ref_sm1_relerror=nle_config->ref_sm1_error / nle_config->ref_sm1;
  nle_config->relerror[6]=nle_config->ref_sm1_relerror;
  nle_config->ref_sm2_relerror=nle_config->ref_sm2_error / nle_config->ref_sm2;
  nle_config->relerror[7]=nle_config->ref_sm2_relerror;
  nle_config->ref_sm3_relerror=nle_config->ref_sm3_error / nle_config->ref_sm3;
  nle_config->relerror[8]=nle_config->ref_sm3_relerror;
  nle_config->ref_sin2w_relerror=nle_config->ref_sin2w_error / nle_config->ref_sin2w;
  nle_config->relerror[9]=nle_config->ref_sin2w_relerror;
  nle_config->ref_alpha_em_relerror=nle_config->ref_alpha_em_error / nle_config->ref_alpha_em;
  nle_config->relerror[10]=nle_config->ref_alpha_em_relerror;
  nle_config->ref_alpha_w_relerror=nle_config->ref_alpha_w_error / nle_config->ref_alpha_w;
  nle_config->relerror[11]=nle_config->ref_alpha_w_relerror;

  return(0);
}
