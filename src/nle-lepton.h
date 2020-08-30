#ifndef NLE_LEPTON_H
#define NLE_LEPTON_H

#define NLE_VERSION "4.3.7.1"

typedef struct {
  int G;
  int float_G;
  int v;
  int float_v;
  int mz;
  int float_mz;
  int mw;
  int float_mw;
  int mh0;
  int float_mh0;
  int m_user;
  int float_muser;
  int float_sm1;
  int float_sm2;
  int float_sm3;
  int sin2w;
  int float_sin2w;
  int alpha_em;
  int float_alpha_em;
  int alpha_w;
  int float_alpha_w;
  int mw_mz_mode;
} nle_input_use_t;

typedef struct {
  int smrfactor_rational_up;
  int smrfactor_rational_down;
  int smrfactor_2_exp_up;
  int smrfactor_2_exp_down;
  int smrfactor_alpha_exp_up;
  int smrfactor_alpha_exp_down;
  int smrfactor_pi_exp_up;
  int smrfactor_pi_exp_down;
  int smrfactor_user_exp_up;
  int smrfactor_user_exp_down;
  int smrfactor_complexity;
  nle_input_use_t smrfactor_uses;
  double smrfactor_multiplier;
} nle_smrfactor_precomputed_t;

typedef struct {
  int infactor_rational_up;
  int infactor_rational_down;
  int infactor_2_exp_up;
  int infactor_2_exp_down;
  int infactor_alpha_exp_up;
  int infactor_alpha_exp_down;
  int infactor_pi_exp_up;
  int infactor_pi_exp_down;
  int infactor_nss;
  int infactor_nbv;
  int infactor_user_exp_up;
  int infactor_user_exp_down;
  int infactor_complexity;
  nle_input_use_t infactor_uses;
  double infactor_multiplier[27];
} nle_infactor_precomputed_t;

typedef struct {
  int outfactor_rational_up;
  int outfactor_rational_down;
  int outfactor_2_exp_up;
  int outfactor_2_exp_down;
  int outfactor_alpha_exp_up;
  int outfactor_alpha_exp_down;
  int outfactor_pi_exp_up;
  int outfactor_pi_exp_down;
  int outfactor_user1_exp_up;
  int outfactor_user1_exp_down;
  int outfactor_user2_exp_up;
  int outfactor_user2_exp_down;
  int outfactor_user3_exp_up;
  int outfactor_user3_exp_down;
  int outfactor_complexity;
  nle_input_use_t outfactor_uses;
  double outfactor_multiplier;
} nle_outfactor_precomputed_t;

typedef struct {
  int outfactor_sin2w_exp_up;
  int outfactor_sin2w_exp_down;
  int outfactor_cos2w_exp_up;
  int outfactor_cos2w_exp_down;
  int outfactor_rmr_exp_up;
  int outfactor_rmr_exp_down;
  int outfactor_rmr_mass_id_up;
  int outfactor_rmr_mass_id_down;
  int dynamicfactor_complexity;
  nle_input_use_t dynamicfactor_uses;
  double dynamicfactor_multiplier;
} nle_dynamicfactor_precomputed_t;

typedef struct {
  int term_id;
  int exp_inv;
  int smrfactor_mass_id;
  int infactor_rational_up;
  int infactor_rational_down;
  int infactor_2_exp_up;
  int infactor_2_exp_down;
  int infactor_alpha_exp_up;
  int infactor_alpha_exp_down;
  int infactor_pi_exp_up;
  int infactor_pi_exp_down;
  int infactor_nss;
  int infactor_nbv;
  int infactor_user_exp_up;
  int infactor_user_exp_down;
  int outfactor_rational_up;
  int outfactor_rational_down;
  int outfactor_2_exp_up;
  int outfactor_2_exp_down;
  int outfactor_alpha_exp_up;
  int outfactor_alpha_exp_down;
  int outfactor_pi_exp_up;
  int outfactor_pi_exp_down;
  int outfactor_sin2w_exp_up;
  int outfactor_sin2w_exp_down;
  int outfactor_cos2w_exp_up;
  int outfactor_cos2w_exp_down;
  int outfactor_rmr_exp_up;
  int outfactor_rmr_exp_down;
  int outfactor_rmr_mass_id_up;
  int outfactor_rmr_mass_id_down;
  int outfactor_user1_exp_up;
  int outfactor_user1_exp_down;
  int outfactor_user2_exp_up;
  int outfactor_user2_exp_down;
  int outfactor_user3_exp_up;
  int outfactor_user3_exp_down;
  double static_multiplier;
  double match;
  int match_up;
  int match_down;
  int match_complexity;
  long long match_hash;
  nle_input_use_t match_uses;
} nle_phase1_match_t;

typedef struct {
  int exp_inv;
  int smrfactor_mass_id;
  double smrfactor;
  double coefficient;
  nle_phase1_match_t *matches_start;
  int matches_count;
  nle_phase1_match_t *current_match;
  nle_smrfactor_precomputed_t *current_smrfactors;
} nle_term_state_t;

typedef struct {
  __uint128_t pcg_state;
  long long phase1_seq;
  double input_sample_sm1;
  double input_sample_sm2;
  double input_sample_sm3;
  double input_sample_alpha_em;
  double input_sample_alpha_w;
  double input_sample_v;
  double input_sample_mz;
  double input_sample_mw;
  double input_sample_sin2w;
  double input_sample_G;
  double input_sample_mp;
  double input_sample_mh0;
  double input_sample_muser;
  char exponents_str[20];
  int nle_mixing_polarity;
  char nle_mixing_polarity_str[2];
  int smrfactor_mass_configuration;
  char smrfactor_mass_configuration_str[5];
  nle_term_state_t term1;
  nle_term_state_t term2;
  nle_term_state_t term3;
  nle_term_state_t term4;
  nle_smrfactor_precomputed_t *smrfactors_precomputed_start;
  int smrfactors_precomputed_count;
  nle_infactor_precomputed_t *infactors_precomputed_start;
  int infactors_precomputed_count;
  nle_outfactor_precomputed_t *outfactors_precomputed_start;
  int outfactors_precomputed_count;
  nle_phase1_match_t *phase1_matches_start;
  int phase1_matches_count;
  int terms_matched[5];
  nle_input_use_t all_uses;
  int current_symmetry;
} nle_state_t;

typedef struct {
  long external_seed;
  char config_file_name[256];
  int upload_results_enable;
  char upload_url[256];
  int phase1_status_enable;
  int phase1_solution_detail;
  int phase2_status_enable;
  int nle_mode;
  int smrfactor_1minus_enable;
  int nle_mixing_polarity;
  int smrfactor_mass_configuration;
  int exp_inv_min;
  int exp_inv_max;
  int exp_inv_include;
  int exp_inv_term1_force;
  int exp_inv_term2_force;
  int exp_inv_term3_force;
  int exp_inv_term4_force;
  int exp_inv_2seq_limit;
  int exp_inv_3seq_limit;
  int exp_inv_4seq_limit;
  int exp_pos_enable;
  int exp_neg_enable;
  int phase1_run_continuous;
  int phase1_random_samples_enable;
  int phase1_filter;
  int phase1_int_match_max;
  int phase1_int_match_filter;
  double phase1_two_term_test_min;
  double phase1_two_term_test_max;
  int phase1_unsolvable_checkpoint;
  long long phase1_mc_samples_limit;
  int phase2_enable;
  int phase2_symmetry_min;
  int phase2_complexity_max;
  int phase2_check_nbv_nss;
  int phase2_check_weak;
  int phase2_check_rmr;
  double phase2_results_window;
  int phase2_results_always; 
  int phase2_ignore_small_rel_unc;
  int smrfactor_mass_mp_enable;
  int smrfactor_mass_v_enable;
  int smrfactor_mass_mz_enable;
  int smrfactor_mass_mw_enable;
  int smrfactor_mass_mh0_enable;
  int smrfactor_mass_user_enable;
  int smrfactor_mass_user_random;
  int smrfactor_mass_user_scan;
  double smrfactor_mass_user_min;
  double smrfactor_mass_user_max;
  double smrfactor_mass_user_step;
  double smrfactor_mass_user;
  double smrfactor_mass_user_error;
  double smrfactor_mass_user_relerror;
  int smrfactor_gt_sm3;
  int smrfactor_lt_sm1;
  int smrfactor_rational_max;
  int smrfactor_rational_filter;
  int smrfactor_2_exp_up_max;
  int smrfactor_2_exp_down_max;
  int smrfactor_alpha_exp_up_max;
  int smrfactor_alpha_exp_down_max;
  int smrfactor_pi_exp_up_max;
  int smrfactor_pi_exp_down_max;
  double smrfactor_user;
  int smrfactor_user_exp_up_max;
  int smrfactor_user_exp_down_max;
  int infactor_rational_max;
  int infactor_rational_filter;
  int infactor_2_exp_up_max;
  int infactor_2_exp_down_max;
  int infactor_alpha_exp_up_max;
  int infactor_alpha_exp_down_max;
  int infactor_pi_exp_up_max;
  int infactor_pi_exp_down_max;
  int infactor_nss_enable;
  int infactor_nbv_enable;
  double infactor_user;
  int infactor_user_exp_up_max;
  int infactor_user_exp_down_max;
  int outfactor_rational_max;
  int outfactor_rational_filter;
  int outfactor_2_exp_up_max;
  int outfactor_2_exp_down_max;
  int outfactor_alpha_exp_up_max;
  int outfactor_alpha_exp_down_max;
  int outfactor_pi_exp_up_max;
  int outfactor_pi_exp_down_max;
  int outfactor_weak_exp_up_max;
  int outfactor_weak_exp_down_max;
  double outfactor_user1;
  int outfactor_user1_exp_up_max;
  int outfactor_user1_exp_down_max;
  double outfactor_user2;
  int outfactor_user2_exp_up_max;
  int outfactor_user2_exp_down_max;
  double outfactor_user3;
  int outfactor_user3_exp_up_max;
  int outfactor_user3_exp_down_max;
  int outfactor_rmr_exp_up_max;
  int outfactor_rmr_exp_down_max;
  int outfactor_rmr_mp_enable;
  int outfactor_rmr_v_enable;
  int outfactor_rmr_mz_enable;
  int outfactor_rmr_mw_enable;
  int outfactor_rmr_mh0_enable;
  int outfactor_rmr_user_enable;
  double ref_sm1;
  double ref_sm1_error;
  double ref_sm1_relerror;
  double ref_sm2;
  double ref_sm2_error;
  double ref_sm2_relerror;
  double ref_sm3;
  double ref_sm3_error;
  double ref_sm3_relerror;
  double ref_c;
  double ref_h;
  double ref_hbar;
  double ref_e;
  double ref_ev_to_kg;
  double ref_kg_to_ev;
  double ref_alpha_em;
  double ref_alpha_em_error;
  double ref_alpha_em_relerror;
  double ref_alpha_w;
  double ref_alpha_w_error;
  double ref_alpha_w_relerror;
  double ref_v;
  double ref_v_error;
  double ref_v_relerror;
  double ref_mz;
  double ref_mz_error;
  double ref_mz_relerror;
  double ref_mw;
  double ref_mw_error;
  double ref_mw_relerror;
  double ref_mh0;
  double ref_mh0_error;
  double ref_mh0_relerror;
  double ref_sin2w;
  double ref_sin2w_error;
  double ref_sin2w_relerror;
  double ref_G;
  double ref_G_error;
  double ref_G_relerror;
  double relerror[12];
} nle_config_t;

#endif // NLE_LEPTON_H
