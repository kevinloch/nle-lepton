#include "nle-lepton.h"
#include <stdio.h>
#include <stdlib.h>
#include "util.h"

void generateExponents(nle_config_t *nle_config, nle_state_t *nle_state) {
  int valid;
  double r;
  int exp_inv_1=1;
  int exp_inv_2=1;
  int exp_inv_3=1;

  valid=0;
  while (!valid) {
    valid=1; // assumme exponents are valid unless checks fail

    // generate 2 or 3 random exponents
    r=pcg_ldrand64(nle_state);
    exp_inv_1=(int)(r * 2 * ((double)nle_config->exp_inv_max + 0.5)) - nle_config->exp_inv_max;
    while (exp_inv_1 == 0) {
      r=pcg_ldrand64(nle_state);
      exp_inv_1=(int)(r * 2 * ((double)nle_config->exp_inv_max + 0.5)) - nle_config->exp_inv_max;
    }
    r=pcg_ldrand64(nle_state);
    exp_inv_2=(int)(r * 2 * ((double)nle_config->exp_inv_max + 0.5)) - nle_config->exp_inv_max;
    while ((exp_inv_2 == 0) || ((nle_config->smrfactor_1minus_enable == 0) && (exp_inv_2 == exp_inv_1))) { // same exponents are ok in (1-smr) mode
      r=pcg_ldrand64(nle_state);
      exp_inv_2=(int)(r * 2 * ((double)nle_config->exp_inv_max + 0.5)) - nle_config->exp_inv_max;
    }
    if (nle_config->nle_mode > 2) {
      r=pcg_ldrand64(nle_state);
      exp_inv_3=(int)(r * 2 * ((double)nle_config->exp_inv_max + 0.5)) - nle_config->exp_inv_max;
      while ((exp_inv_3 == 0) || (exp_inv_3 == exp_inv_1) || (exp_inv_3 == exp_inv_2)) {
        r=pcg_ldrand64(nle_state);
        exp_inv_3=(int)(r * 2 * ((double)nle_config->exp_inv_max + 0.5)) - nle_config->exp_inv_max;
      }
    }

    // set nle_state exponent variables, sort if nle_mode=3
    if (nle_config->nle_mode == 2) {
      nle_state->term1.exp_inv=exp_inv_1;
      nle_state->term2.exp_inv=exp_inv_2;
      nle_state->term3.exp_inv=1; // force term3 to 1 as it is used for synthetic middle coefficient in 2-term mixed mode
    } else if (nle_config->nle_mode == 3) {
      // sort exponents to ensure real roots
      if ((exp_inv_1 < exp_inv_2) && (exp_inv_1 < exp_inv_3)) {
        nle_state->term1.exp_inv = exp_inv_1;
      }
      if ((exp_inv_2 < exp_inv_1) && (exp_inv_2 < exp_inv_3)) {
        nle_state->term1.exp_inv = exp_inv_2;
      }
      if ((exp_inv_3 < exp_inv_1) && (exp_inv_3 < exp_inv_2)) {
        nle_state->term1.exp_inv = exp_inv_3;
      }
      if (((exp_inv_1 < exp_inv_2) && (exp_inv_1 > exp_inv_3)) || ((exp_inv_1 > exp_inv_2) && (exp_inv_1 < exp_inv_3))) {
        nle_state->term2.exp_inv = exp_inv_1;
      }
      if (((exp_inv_2 < exp_inv_1) && (exp_inv_2 > exp_inv_3)) || ((exp_inv_2 > exp_inv_1) && (exp_inv_2 < exp_inv_3))) {
        nle_state->term2.exp_inv = exp_inv_2;
      }
      if (((exp_inv_3 < exp_inv_1) && (exp_inv_3 > exp_inv_2)) || ((exp_inv_3 > exp_inv_1) && (exp_inv_3 < exp_inv_2))) {
        nle_state->term2.exp_inv = exp_inv_3;
      }
      if ((exp_inv_1 > exp_inv_2) && (exp_inv_1 > exp_inv_3)) {
        nle_state->term3.exp_inv = exp_inv_1;
      }
      if ((exp_inv_2 > exp_inv_1) && (exp_inv_2 > exp_inv_3)) {
        nle_state->term3.exp_inv = exp_inv_2;
      }
      if ((exp_inv_3 > exp_inv_1) && (exp_inv_3 > exp_inv_2)) {
        nle_state->term3.exp_inv = exp_inv_3;
      }
    }

    // enforce exponent sign restrictions
    if (nle_config->nle_mode == 2) {
      if ((nle_config->exp_pos_enable == 0) && ((nle_state->term1.exp_inv > 0) || (nle_state->term2.exp_inv > 0))) {
        valid=0;
      }
      if ((nle_config->exp_neg_enable == 0) && ((nle_state->term1.exp_inv < 0) || (nle_state->term2.exp_inv < 0))) {
        valid=0;
      }
    } else {
      if ((nle_config->exp_pos_enable == 0) && ((nle_state->term1.exp_inv > 0) || (nle_state->term2.exp_inv > 0) || (nle_state->term3.exp_inv > 0))) {
        valid=0;
      }
      if ((nle_config->exp_neg_enable == 0) && ((nle_state->term1.exp_inv < 0) || (nle_state->term2.exp_inv < 0) || (nle_state->term3.exp_inv < 0))) {
        valid=0;
      }
    }

    // check if exp_inv_include is set and enforce
    if (nle_config->exp_inv_include != 0) {
      if ((nle_state->term1.exp_inv != nle_config->exp_inv_include) && (nle_state->term1.exp_inv != -nle_config->exp_inv_include) && (nle_state->term2.exp_inv != nle_config->exp_inv_include)&& (nle_state->term2.exp_inv != -nle_config->exp_inv_include)) {
        valid=0;
      }
    }

    // enforce exp_inv_min
    if ((nle_state->term1.exp_inv < nle_config->exp_inv_min) && (nle_state->term1.exp_inv > -nle_config->exp_inv_min) && (nle_state->term2.exp_inv < nle_config->exp_inv_min)&& (nle_state->term2.exp_inv > -nle_config->exp_inv_min)) {
      valid=0;
    }

    // special exponent checks for 2-term mode
    if (nle_config->nle_mode == 2) {
      if (nle_config->smrfactor_1minus_enable == 0) {
        if ((nle_state->term1.exp_inv * nle_state->term2.exp_inv) < 0) {
          // if terms have opposite sign, this is not supported in 2-term mode without (1-smr)
          valid=0;
        }
        if ((abs(nle_state->term1.exp_inv) == 1) || (abs(nle_state->term2.exp_inv) == 1)) {
          // if either term is +/- 1, this is not supported in 2-term mode without 1-smr mode
          valid=0;
        }
      } else if (nle_config->smrfactor_1minus_enable == 1) {
        // special exponent checks for 1-smr mode
        // suppress trivial geometries
        if ((abs(nle_state->term1.exp_inv) == 1) && (abs(nle_state->term2.exp_inv) == 1)) {
          valid=0;
        }
        if ((abs(nle_state->term1.exp_inv) == 2) && (abs(nle_state->term2.exp_inv) == 2)) {
          valid=0;
        }
      }
    }

    // check if exponents violate sequential limits
    if ((nle_state->term2.exp_inv == nle_state->term1.exp_inv+1) && ((abs(nle_state->term1.exp_inv) >= nle_config->exp_inv_2seq_limit) || ((abs(nle_state->term2.exp_inv) >= nle_config->exp_inv_2seq_limit)))) {
      valid=0;
    }
    if ((nle_config->nle_mode > 2) && (nle_state->term3.exp_inv == nle_state->term2.exp_inv+1) && ((abs(nle_state->term2.exp_inv) >= nle_config->exp_inv_2seq_limit) || ((abs(nle_state->term3.exp_inv) >= nle_config->exp_inv_2seq_limit)))) {
      valid=0;
    }
    if ((nle_config->nle_mode > 2) && (nle_state->term3.exp_inv == nle_state->term2.exp_inv+1) && (nle_state->term2.exp_inv == nle_state->term1.exp_inv+1)\
        && ((abs(nle_state->term1.exp_inv) >= nle_config->exp_inv_2seq_limit) || ((abs(nle_state->term2.exp_inv) >= nle_config->exp_inv_3seq_limit)) || ((abs(nle_state->term3.exp_inv) >= nle_config->exp_inv_3seq_limit)) )) {
      valid=0;
    }

    // optionally force exponents for debugging or focused searches
    if (nle_config->exp_inv_term1_force != 0) {
      nle_state->term1.exp_inv=nle_config->exp_inv_term1_force;
      valid=1; // override failed checks if forced
    }
    if (nle_config->exp_inv_term2_force != 0) {
      nle_state->term2.exp_inv=nle_config->exp_inv_term2_force;
      valid=1; // override failed checks if forced
    }
    if (nle_config->exp_inv_term3_force != 0) {
      nle_state->term3.exp_inv=nle_config->exp_inv_term3_force;
      valid=1; // override failed checks if forced
    }
  } // end while !valid

  // set exponents string
  nle_state->exponents_str[19]=0;
  if (nle_config->nle_mode == 2) {
    sprintf(nle_state->exponents_str, "E%+d%+d", nle_state->term1.exp_inv, nle_state->term2.exp_inv);
  } else if (nle_config->nle_mode == 3) {
    sprintf(nle_state->exponents_str, "E%+d%+d%+d", nle_state->term1.exp_inv, nle_state->term2.exp_inv, nle_state->term3.exp_inv);
  }
  nle_state->phase1_matches_count=0;
}
