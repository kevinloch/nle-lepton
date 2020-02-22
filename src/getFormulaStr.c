#include <stdio.h>
#include "nle-lepton.h"

void getFormulaStr(char *formula_str, nle_phase1_match_t *current_match, int nle_mode) {
  // creates a complete text string representation of a single formula term
  int upout;
  int downout;
  char e2out[32];
  char updownout[32];
  char piout[32];
  char aout[32];
  char uout[32];
  char s2wout[32];
  char c2wout[32];
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

  if (current_match->outfactor_user_exp_up == 0) {
    sprintf(uout, "           ");
  } else {
    if ((current_match->outfactor_user_exp_up == -1) && (current_match->outfactor_user_exp_down == 1)) {
      sprintf(uout, "uout       ");
    } else {
      sprintf(uout, "uout^(%2d/%d)", -(current_match->outfactor_user_exp_up), current_match->outfactor_user_exp_down);
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
    sprintf(massstr,       "    M/mp    ");
    sprintf(massstrinv,    "    mp/M    ");
  } else if (current_match->smrfactor_mass == 1) {
    if ((nle_mode == 2) && (current_match->term_id == 3)) {
      // for 2-term mixed mode, we use smrfactor_mass == 1 as placeholder for third term so we don't want it to show in equation str
      sprintf(massstr,       "            ");
      sprintf(massstrinv,    "            ");
    } else {
      sprintf(massstr,       "    M/v     ");
      sprintf(massstrinv,    "    v/M     ");
    }
  } else if (current_match->smrfactor_mass == 2) {
    sprintf(massstr,       "    M/mZ    ");
    sprintf(massstrinv,    "    mZ/M    ");
  } else if (current_match->smrfactor_mass == 3) {
    sprintf(massstr,       "    M/mW    ");
    sprintf(massstrinv,    "    mW/M    ");
  } else if (current_match->smrfactor_mass == 4) {
    sprintf(massstr,       "    M/mH0   ");
    sprintf(massstrinv,    "    mH0/M   ");
  } else if (current_match->smrfactor_mass == 5) {
    sprintf(massstr,       "  M/m_user  ");
    sprintf(massstrinv,    "  m_user/M  ");
  }

  if (current_match->exp_inv == 1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s  %s %s %s %s %s %s %s       '", updownout, e2out, piout, aout, uout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstr);
  } else if (current_match->exp_inv == -1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s  %s %s %s %s %s %s %s       '", updownout, e2out, piout, aout, uout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstrinv);
  } else if (current_match->exp_inv > 1) {
    sprintf(formula_str, "'%s %s %s %s %s %s %s (%s %s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, uout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstr, current_match->exp_inv);
  } else {
    sprintf(formula_str, "'%s %s %s %s %s %s %s (%s %s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, uout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, uin, massstrinv, -current_match->exp_inv);
  }
}
