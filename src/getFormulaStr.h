#ifndef NLE_FORMULA_H
#define NLE_FORMULA_H

void getSmrfStr(nle_config_t *nle_config, char *smrf_str, nle_smrfactor_precomputed_t *current_smrfactors, double smrfactor);
void getRmrfStr(nle_config_t *nle_config, char *rmrf_str, nle_rmrfactor_precomputed_t *current_rmrfactors, double rmrfactor);
void getFormulaStr(nle_config_t *nle_config, nle_state_t *nle_state, char *formula_str, nle_phase1_match_t *current_match);

#endif // NLE_FORMULA_H
