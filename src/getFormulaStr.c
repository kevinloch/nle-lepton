#include <stdio.h>
#include "nle-lepton.h"

void getSmrfStr(nle_config_t *nle_config, char *smrf_str, nle_smrfactor_precomputed_t *current_smrfactors, double smrfactor) {
  // creates a text string representation of the current solution mass ratio factors

  int upsmr;
  int downsmr;
  char e2smr[32];
  char updownsmr[32];
  char pismr[32];
  char asmr[32];
  char usmr[32];

  upsmr=current_smrfactors->smrfactor_rational_up;
  downsmr=current_smrfactors->smrfactor_rational_down;
  if ((upsmr == 1) && (downsmr == 1)) {
    sprintf(updownsmr, "       ");
  } else {
    sprintf(updownsmr, "(%2d/%2d)", upsmr, downsmr);
  }

  if (current_smrfactors->smrfactor_2_exp_up == 0) {
    sprintf(e2smr, "        ");
  } else {
    sprintf(e2smr, "2^(%2d/%d)", current_smrfactors->smrfactor_2_exp_up, current_smrfactors->smrfactor_2_exp_down);
  }

  if (current_smrfactors->smrfactor_pi_exp_up == 0) {
    sprintf(pismr, "         ");
  } else {
    if ((current_smrfactors->smrfactor_pi_exp_up == 1) && (current_smrfactors->smrfactor_pi_exp_down == 1)) {
      sprintf(pismr, "pi       ");
    } else {
      sprintf(pismr, "pi^(%2d/%d)", current_smrfactors->smrfactor_pi_exp_up, current_smrfactors->smrfactor_pi_exp_down);
    }
  }

  if (current_smrfactors->smrfactor_alpha_exp_up == 0) {
    sprintf(asmr, "        ");
  } else {
    if ((current_smrfactors->smrfactor_alpha_exp_up == 1) && (current_smrfactors->smrfactor_alpha_exp_down == 1)) {
      sprintf(asmr, "a       ");
    } else {
      sprintf(asmr, "a^(%2d/%d)", current_smrfactors->smrfactor_alpha_exp_up, current_smrfactors->smrfactor_alpha_exp_down);
    }
  }

  if (current_smrfactors->smrfactor_user_exp_up == 0) {
    sprintf(usmr, "           ");
  } else {
    if ((current_smrfactors->smrfactor_user_exp_up == 1) && (current_smrfactors->smrfactor_user_exp_down == 1)) {
      sprintf(usmr, "usmr       ");
    } else {
      sprintf(usmr, "usmr^(%2d/%d)", current_smrfactors->smrfactor_user_exp_up, current_smrfactors->smrfactor_user_exp_down);
    }
  }

  sprintf(smrf_str, "'%s %s %s %s %s'", updownsmr, e2smr, pismr, asmr, usmr);
}

void getRmrfStr(nle_config_t *nle_config, char *rmrf_str, nle_rmrfactor_precomputed_t *current_rmrfactors, double rmrfactor) {
  // creates a text string representation of the current solution mass ratio factors

  int uprmr;
  int downrmr;
  char e2rmr[32];
  char updownrmr[32];
  char pirmr[32];
  char armr[32];
  char urmr[32];

  uprmr=current_rmrfactors->rmrfactor_rational_up;
  downrmr=current_rmrfactors->rmrfactor_rational_down;
  if ((uprmr == 1) && (downrmr == 1)) {
    sprintf(updownrmr, "       ");
  } else {
    sprintf(updownrmr, "(%2d/%2d)", uprmr, downrmr);
  }

  if (current_rmrfactors->rmrfactor_2_exp_up == 0) {
    sprintf(e2rmr, "        ");
  } else {
    sprintf(e2rmr, "2^(%2d/%d)", current_rmrfactors->rmrfactor_2_exp_up, current_rmrfactors->rmrfactor_2_exp_down);
  }

  if (current_rmrfactors->rmrfactor_pi_exp_up == 0) {
    sprintf(pirmr, "         ");
  } else {
    if ((current_rmrfactors->rmrfactor_pi_exp_up == 1) && (current_rmrfactors->rmrfactor_pi_exp_down == 1)) {
      sprintf(pirmr, "pi       ");
    } else {
      sprintf(pirmr, "pi^(%2d/%d)", current_rmrfactors->rmrfactor_pi_exp_up, current_rmrfactors->rmrfactor_pi_exp_down);
    }
  }

  if (current_rmrfactors->rmrfactor_alpha_exp_up == 0) {
    sprintf(armr, "        ");
  } else {
    if ((current_rmrfactors->rmrfactor_alpha_exp_up == 1) && (current_rmrfactors->rmrfactor_alpha_exp_down == 1)) {
      sprintf(armr, "a       ");
    } else {
      sprintf(armr, "a^(%2d/%d)", current_rmrfactors->rmrfactor_alpha_exp_up, current_rmrfactors->rmrfactor_alpha_exp_down);
    }
  }

  if (current_rmrfactors->rmrfactor_user_exp_up == 0) {
    sprintf(urmr, "           ");
  } else {
    if ((current_rmrfactors->rmrfactor_user_exp_up == 1) && (current_rmrfactors->rmrfactor_user_exp_down == 1)) {
      sprintf(urmr, "urmr       ");
    } else {
      sprintf(urmr, "urmr^(%2d/%d)", current_rmrfactors->rmrfactor_user_exp_up, current_rmrfactors->rmrfactor_user_exp_down);
    }
  }

  sprintf(rmrf_str, "'%s %s %s %s %s'", updownrmr, e2rmr, pirmr, armr, urmr);
}

void getFormulaStr(nle_config_t *nle_config, nle_state_t *nle_state, char *formula_str, nle_phase1_match_t *current_match) {
  // creates a complete text string representation of a single formula term
  int upout;
  int downout;
  char e2out[32];
  char updownout[32];
  char piout[32];
  char aout[32];
  char uout1[32];
  char uout2[32];
  char uout3[32];
  char s2wout[32];
  char c2wout[32];
  char outfactor_rmr_mass_up[16];
  char outfactor_rmr_mass_down[16];
  char rmrout[64];
  char updownin[32];
  char nbin[32];
  char e2in[32];
  char piin[32];
  char ain[32];
  char uin[32];
  char massstr[48];
  char massstrinv[48];
  char rmrfactor_mass_str_up[32];
  char rmrfactor_mass_str_down[32];
  char smrfactor_mass_str[32];

  // note: all terms except match_up and match_down are inverted here as they represent offsets to the real coefficient
  upout=current_match->match_up * current_match->outfactor_rational_down;
  downout=current_match->match_down * current_match->outfactor_rational_up;
  if ((upout == 1) && (downout == 1)) {
    sprintf(updownout, "       ");
  } else {
    sprintf(updownout, "(%2d/%2d)", upout, downout);
  }

  if (current_match->outfactor_2_exp_up == 0) {
    sprintf(e2out, "        ");
  } else {
    sprintf(e2out, "2^(%2d/%d)", -(current_match->outfactor_2_exp_up), current_match->outfactor_2_exp_down);
  }

  if (current_match->outfactor_pi_exp_up == 0) {
    sprintf(piout, "         ");
  } else {
    if ((current_match->outfactor_pi_exp_up == -1) && (current_match->outfactor_pi_exp_down == 1)) {
      sprintf(piout, "pi       ");
    } else {
      sprintf(piout, "pi^(%2d/%d)", -(current_match->outfactor_pi_exp_up), current_match->outfactor_pi_exp_down);
    }
  }

  if (current_match->outfactor_alpha_exp_up == 0) {
    sprintf(aout, "        ");
  } else {
    if ((current_match->outfactor_alpha_exp_up == -1) && (current_match->outfactor_alpha_exp_down == 1)) {
      sprintf(aout, "a       ");
    } else {
      sprintf(aout, "a^(%2d/%d)", -(current_match->outfactor_alpha_exp_up), current_match->outfactor_alpha_exp_down);
    }
  }

  if (current_match->outfactor_user1_exp_up == 0) {
    sprintf(uout1, "            ");
  } else {
    if ((current_match->outfactor_user1_exp_up == -1) && (current_match->outfactor_user1_exp_down == 1)) {
      sprintf(uout1, "uout1       ");
    } else {
      sprintf(uout1, "uout1^(%2d/%d)", -(current_match->outfactor_user1_exp_up), current_match->outfactor_user1_exp_down);
    }
  }

  if (current_match->outfactor_user2_exp_up == 0) {
    sprintf(uout2, "            ");
  } else {
    if ((current_match->outfactor_user2_exp_up == -1) && (current_match->outfactor_user2_exp_down == 1)) {
      sprintf(uout2, "uout2       ");
    } else {
      sprintf(uout2, "uout2^(%2d/%d)", -(current_match->outfactor_user2_exp_up), current_match->outfactor_user2_exp_down);
    }
  }

  if (current_match->outfactor_user3_exp_up == 0) {
    sprintf(uout3, "            ");
  } else {
    if ((current_match->outfactor_user3_exp_up == -1) && (current_match->outfactor_user3_exp_down == 1)) {
      sprintf(uout3, "uout3       ");
    } else {
      sprintf(uout3, "uout3^(%2d/%d)", -(current_match->outfactor_user3_exp_up), current_match->outfactor_user3_exp_down);
    }
  }

  if (current_match->outfactor_rmr_exp_up == 0) {
    sprintf(rmrout, "                     ");
  } else {
    if (current_match->outfactor_rmr_mass_id_up == 0) {
      sprintf(outfactor_rmr_mass_up, "mP    ");
    } else if (current_match->outfactor_rmr_mass_id_up == 1) {
      sprintf(outfactor_rmr_mass_up, "v     ");
    } else if (current_match->outfactor_rmr_mass_id_up == 2) {
      sprintf(outfactor_rmr_mass_up, "mz    ");
    } else if (current_match->outfactor_rmr_mass_id_up == 3) {
      sprintf(outfactor_rmr_mass_up, "mw    ");
    } else if (current_match->outfactor_rmr_mass_id_up == 4) {
      sprintf(outfactor_rmr_mass_up, "mh0   ");
    } else if (current_match->outfactor_rmr_mass_id_up == 5) {
      sprintf(outfactor_rmr_mass_up, "m_user");
    }
    if (current_match->outfactor_rmr_mass_id_down == 0) {
      sprintf(outfactor_rmr_mass_down, "    mP");
    } else if (current_match->outfactor_rmr_mass_id_down == 1) {
      sprintf(outfactor_rmr_mass_down, "     v");
    } else if (current_match->outfactor_rmr_mass_id_down == 2) {
      sprintf(outfactor_rmr_mass_down, "    mz");
    } else if (current_match->outfactor_rmr_mass_id_down == 3) {
      sprintf(outfactor_rmr_mass_down, "    mw");
    } else if (current_match->outfactor_rmr_mass_id_down == 4) {
      sprintf(outfactor_rmr_mass_down, "   mh0");
    } else if (current_match->outfactor_rmr_mass_id_down == 5) {
      sprintf(outfactor_rmr_mass_down, "m_user");
    }
    if ((current_match->outfactor_rmr_exp_up == 1) && (current_match->outfactor_rmr_exp_down == 1)) {
      sprintf(rmrout, "(%s/%s)      ", outfactor_rmr_mass_down, outfactor_rmr_mass_up); // down/up is reversed since we are translating coefficient factors back into the original formula and we don't allow negative exponents on rmr
    } else {
      sprintf(rmrout, "(%s/%s)^(%d/%d)", outfactor_rmr_mass_down, outfactor_rmr_mass_up, current_match->outfactor_rmr_exp_up, current_match->outfactor_rmr_exp_down); // down/up is reversed since we are translating coefficient factors back into the original formula and we don't allow negative exponents on rmr
    }
  }

  if (current_match->outfactor_sin2w_exp_up == 0) {
    sprintf(s2wout, "            ");
  } else {
    if ((current_match->outfactor_sin2w_exp_up == -1) && (current_match->outfactor_sin2w_exp_down == 1)) {
      sprintf(s2wout, "sin2w       ");
    } else {
      sprintf(s2wout, "sin2w^(%2d/%d)", -(current_match->outfactor_sin2w_exp_up), current_match->outfactor_sin2w_exp_down);
    }
  }

  if (current_match->outfactor_cos2w_exp_up == 0) {
    sprintf(c2wout, "            ");
  } else {
    if ((current_match->outfactor_cos2w_exp_up == -1) && (current_match->outfactor_cos2w_exp_down == 1)) {
      sprintf(c2wout, "cos2w       ");
    } else {
      sprintf(c2wout, "cos2w^(%2d/%d)", -(current_match->outfactor_cos2w_exp_up), current_match->outfactor_cos2w_exp_down);
    }
  }

  if ((current_match->infactor_rational_up == 1) && (current_match->infactor_rational_down == 1)) {
    sprintf(updownin, "       ");
  } else {
    sprintf(updownin, "(%2d/%2d)", current_match->infactor_rational_down, current_match->infactor_rational_up);
  }

  if ((current_match->infactor_nbv == 0) && (current_match->infactor_nss == 0)) {
    sprintf(nbin, "      ");
  } else {
    if (current_match->infactor_nbv == -1) {
      sprintf(nbin, "nbv   ");
    } else if (current_match->infactor_nbv == 1) {
      sprintf(nbin, "nbv^-1");
    } else if (current_match->infactor_nss == -1) {
      sprintf(nbin, "nss   ");
    } else {
      sprintf(nbin, "nss^-1");
    }
  }

  if (current_match->infactor_2_exp_up == 0) {
    sprintf(e2in, "        ");
  } else {
    sprintf(e2in, "2^(%2d/%d)", -(current_match->infactor_2_exp_up), current_match->infactor_2_exp_down);
  }

  if (current_match->infactor_pi_exp_up == 0) {
    sprintf(piin, "         ");
  } else if ((current_match->infactor_pi_exp_up == -1) && (current_match->infactor_pi_exp_down == 1)) {
    sprintf(piin, "pi       ");
  } else {
    sprintf(piin, "pi^(%2d/%d)", -(current_match->infactor_pi_exp_up), current_match->infactor_pi_exp_down);
  }

  if (current_match->infactor_alpha_exp_up == 0) {
    sprintf(ain, "        ");
  } else if ((current_match->infactor_alpha_exp_up == -1) && (current_match->infactor_alpha_exp_down == 1)) {
    sprintf(ain, "a       ");
  } else {
    sprintf(ain, "a^(%2d/%d)", -(current_match->infactor_alpha_exp_up), current_match->infactor_alpha_exp_down);
  }

  if (current_match->infactor_user_exp_up == 0) {
    sprintf(uin, "          ");
  } else if ((current_match->infactor_user_exp_up == -1) && (current_match->infactor_user_exp_down == 1)) {
    sprintf(uin, "uin       ");
  } else {
    sprintf(uin, "uin^(%2d/%d)", -(current_match->infactor_user_exp_up), current_match->infactor_user_exp_down);
  }

  // load rmrfactor strings if (1-rmr-smr) is enabled
  if (nle_config->rmrfactor_1minus_enable == 1) {
    if (nle_state->term1.rmrfactor_mass_id_up == 0) {
      sprintf(rmrfactor_mass_str_up, "mP");
    } else if (nle_state->term1.rmrfactor_mass_id_up == 1) {
      sprintf(rmrfactor_mass_str_up, "v");
    } else if (nle_state->term1.rmrfactor_mass_id_up == 2) {
      sprintf(rmrfactor_mass_str_up, "mz");
    } else if (nle_state->term1.rmrfactor_mass_id_up == 3) {
      sprintf(rmrfactor_mass_str_up, "mw");
    } else if (nle_state->term1.rmrfactor_mass_id_up == 4) {
      sprintf(rmrfactor_mass_str_up, "mh0");
    } else if (nle_state->term1.rmrfactor_mass_id_up == 5) {
      sprintf(rmrfactor_mass_str_up, "muser");
    }
    if (nle_state->term1.rmrfactor_mass_id_down == 0) {
      sprintf(rmrfactor_mass_str_down, "mP");
    } else if (nle_state->term1.rmrfactor_mass_id_down == 1) {
      sprintf(rmrfactor_mass_str_down, "v");
    } else if (nle_state->term1.rmrfactor_mass_id_down == 2) {
      sprintf(rmrfactor_mass_str_down, "mz");
    } else if (nle_state->term1.rmrfactor_mass_id_down == 3) {
      sprintf(rmrfactor_mass_str_down, "mw");
    } else if (nle_state->term1.rmrfactor_mass_id_down == 4) {
      sprintf(rmrfactor_mass_str_down, "mh0");
    } else if (nle_state->term1.rmrfactor_mass_id_down == 5) {
      sprintf(rmrfactor_mass_str_down, "muser");
    }
  }

  // load smrfactor mass string
  if (current_match->smrfactor_mass_id == 0) {
    sprintf(smrfactor_mass_str, "mP");
  } else if (current_match->smrfactor_mass_id == 1) {
    sprintf(smrfactor_mass_str, "v");
  } else if (current_match->smrfactor_mass_id == 2) {
    sprintf(smrfactor_mass_str, "mz");
  } else if (current_match->smrfactor_mass_id == 3) {
    sprintf(smrfactor_mass_str, "mw");
  } else if (current_match->smrfactor_mass_id == 4) {
    sprintf(smrfactor_mass_str, "mh0");
  } else if (current_match->smrfactor_mass_id == 5) {
    sprintf(smrfactor_mass_str, "muser");
  }

  if ((nle_config->rmrfactor_1minus_enable == 1) && (current_match->term_id == 1)) { 
    sprintf(massstr,       "(1-(rmrf*%s/%s)-(smrf*M/%s))", rmrfactor_mass_str_up, rmrfactor_mass_str_down, smrfactor_mass_str);
    sprintf(massstrinv,    "1/(1-(rmrf*%s/%s)-(smrf*M/%s))", rmrfactor_mass_str_up, rmrfactor_mass_str_down, smrfactor_mass_str);
  } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 1)) {
    sprintf(massstr,       "(1-(smrf*M/%s))", smrfactor_mass_str);
    sprintf(massstrinv,    "1/(1-(smrf*M/%s))", smrfactor_mass_str);
  } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 2)) {
    sprintf(massstr,       "       smrf*M/%s     ", smrfactor_mass_str);
    sprintf(massstrinv,    "      %s/(smrf*M)    ", smrfactor_mass_str);
  } else if ((nle_config->nle_mode == 2) && (current_match->term_id == 3)) {
    sprintf(massstr,       "                     ");
    sprintf(massstrinv,    "                     ");
  } else {
    sprintf(massstr,       "         M/%s        ", smrfactor_mass_str);
    sprintf(massstrinv,    "         %s/M        ", smrfactor_mass_str);
  }

  if (current_match->exp_inv == 1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s'", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstr);
  } else if (current_match->exp_inv == -1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s'", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstrinv);
  } else if (current_match->exp_inv > 1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s (%s %s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstr, current_match->exp_inv);
  } else {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s (%s %s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstrinv, -current_match->exp_inv);
  }
}
