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

void getFormulaStr(nle_config_t *nle_config, char *formula_str, nle_phase1_match_t *current_match) {
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
  char rmr_mass_up[16];
  char rmr_mass_down[16];
  char rmrout[64];
  char updownin[32];
  char nbin[32];
  char e2in[32];
  char piin[32];
  char ain[32];
  char uin[32];
  char massstr[32];
  char massstrinv[32];

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
      sprintf(rmr_mass_up, "mP    ");
    } else if (current_match->outfactor_rmr_mass_id_up == 1) {
      sprintf(rmr_mass_up, "v     ");
    } else if (current_match->outfactor_rmr_mass_id_up == 2) {
      sprintf(rmr_mass_up, "mz    ");
    } else if (current_match->outfactor_rmr_mass_id_up == 3) {
      sprintf(rmr_mass_up, "mw    ");
    } else if (current_match->outfactor_rmr_mass_id_up == 4) {
      sprintf(rmr_mass_up, "mh0   ");
    } else if (current_match->outfactor_rmr_mass_id_up == 5) {
      sprintf(rmr_mass_up, "m_user");
    }
    if (current_match->outfactor_rmr_mass_id_down == 0) {
      sprintf(rmr_mass_down, "    mP");
    } else if (current_match->outfactor_rmr_mass_id_down == 1) {
      sprintf(rmr_mass_down, "     v");
    } else if (current_match->outfactor_rmr_mass_id_down == 2) {
      sprintf(rmr_mass_down, "    mz");
    } else if (current_match->outfactor_rmr_mass_id_down == 3) {
      sprintf(rmr_mass_down, "    mw");
    } else if (current_match->outfactor_rmr_mass_id_down == 4) {
      sprintf(rmr_mass_down, "   mh0");
    } else if (current_match->outfactor_rmr_mass_id_down == 5) {
      sprintf(rmr_mass_down, "m_user");
    }
    if ((current_match->outfactor_rmr_exp_up == 1) && (current_match->outfactor_rmr_exp_down == 1)) {
      sprintf(rmrout, "(%s/%s)      ", rmr_mass_down, rmr_mass_up); // down/up is reversed since we are translating coefficient factors back into the original formula and we don't allow negative exponents on rmr
    } else {
      sprintf(rmrout, "(%s/%s)^(%d/%d)", rmr_mass_down, rmr_mass_up, current_match->outfactor_rmr_exp_up, current_match->outfactor_rmr_exp_down); // down/up is reversed since we are translating coefficient factors back into the original formula and we don't allow negative exponents on rmr
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

  if (current_match->smrfactor_mass == 0) {
    if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 1)) {
      sprintf(massstr,       "   (1-(smrf*M/mP))   ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 2)) {
      sprintf(massstr,       "       smrf*M/mP     ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->nle_mode == 2) && (current_match->term_id == 3)) {
      sprintf(massstr,       "                     ");
      sprintf(massstrinv,    "                     ");
    } else {
      sprintf(massstr,       "         M/mP        ");
      sprintf(massstrinv,    "         mP/M        ");
    }
  } else if (current_match->smrfactor_mass == 1) {
    if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 1)) {
      sprintf(massstr,       "    (1-(smrf*M/v))   ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 2)) {
      sprintf(massstr,       "        smrf*M/v     ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->nle_mode == 2) && (current_match->term_id == 3)) {
      sprintf(massstr,       "                     ");
      sprintf(massstrinv,    "                     ");
    } else {
      sprintf(massstr,       "         M/v         ");
      sprintf(massstrinv,    "         v/M         ");
    }
  } else if (current_match->smrfactor_mass == 2) {
    if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 1)) {
      sprintf(massstr,       "   (1-(smrf*M/mZ))   ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 2)) {
      sprintf(massstr,       "       smrf*M/mZ     ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->nle_mode == 2) && (current_match->term_id == 3)) {
      sprintf(massstr,       "                     ");
      sprintf(massstrinv,    "                     ");
    } else {
      sprintf(massstr,       "         M/mZ        ");
      sprintf(massstrinv,    "         mZ/M        ");
    }
  } else if (current_match->smrfactor_mass == 3) {
    if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 1)) {
      sprintf(massstr,       "   (1-(smrf*M/mW))   ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 2)) {
      sprintf(massstr,       "       smrf*M/mW     ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->nle_mode == 2) && (current_match->term_id == 3)) {
      sprintf(massstr,       "                     ");
      sprintf(massstrinv,    "                     ");
    } else {
      sprintf(massstr,       "         M/mW        ");
      sprintf(massstrinv,    "         mW/M        ");
    }
  } else if (current_match->smrfactor_mass == 4) {
    if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 1)) {
      sprintf(massstr,       "   (1-(smrf*M/mH0))  ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 2)) {
      sprintf(massstr,       "      smrf*M/mH0     ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->nle_mode == 2) && (current_match->term_id == 3)) {
      sprintf(massstr,       "                     ");
      sprintf(massstrinv,    "                     ");
    } else {
      sprintf(massstr,       "        M/mH0        ");
      sprintf(massstrinv,    "        mH0/M        ");
    }
  } else if (current_match->smrfactor_mass == 5) {
    if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 1)) {
      sprintf(massstr,       " (1-(smrf*M/m_user)) ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->smrfactor_1minus_enable == 1) && (current_match->term_id == 2)) {
      sprintf(massstr,       "     smrf*M/m_user   ");
      sprintf(massstrinv,    "                     ");
    } else if ((nle_config->nle_mode == 2) && (current_match->term_id == 3)) {
      sprintf(massstr,       "                     ");
      sprintf(massstrinv,    "                     ");
    } else {
      sprintf(massstr,       "       M/m_user      ");
      sprintf(massstrinv,    "       m_user/M      ");
    }
  }

  if (current_match->exp_inv == 1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s        '", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstr);
  } else if (current_match->exp_inv == -1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s        '", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstrinv);
  } else if (current_match->exp_inv > 1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s (%s %s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstr, current_match->exp_inv);
  } else {
    sprintf(formula_str, "'%s %s %s %s %s %s %s %s %s %s (%s %s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, uout1, uout2, uout3, rmrout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstrinv, -current_match->exp_inv);
  }
}
