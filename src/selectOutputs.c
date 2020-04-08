#include <stdio.h>
#include "nle-lepton.h"
#include "util.h"

//#define DEBUG_SELECT_OUTPUTS

int selectOutputs(nle_config_t *nle_config, nle_state_t *nle_state) {
  // select the three used variables with the highest experimental uncertainty
  // return the number of unknowns (outputs to solve for in phase 2)

  int unknowns;
  int i;
  int i_1;
  int i_2;
  int i_3;
  double u_1;
  double u_2;
  double u_3;
  double rel_error;
  int reference_mass;

  /*
    Index of variables
    0: alpha_em
    1:  v
    2:  G
    3:  mz
    4:  mw
    5:  mH0
    6:  sin2w
    7:  m_user
    8:  sm1
    9:  sm2
    10: sm3
  */

#ifdef DEBUG_SELECT_OUTPUTS
  printUses(&nle_state->all_uses);
#endif

  // determine if sin2w should be solved directly or derived from mw and mz
  if (nle_state->all_uses.sin2w == 1) {
    if ((nle_state->all_uses.mw == 1) && (nle_state->all_uses.mz == 1)) {
      // sin2 and both mw or mz are used, derive sin2 from mw and mz
      nle_state->all_uses.mw_mz_mode=1;
    } else {
      // sin2 is used but not both mw and mz, sin2 can be used directly
      nle_state->all_uses.mw_mz_mode=0;
    }
  } else {
    // sin2w is not used at all
    nle_state->all_uses.mw_mz_mode=0;
  }

  // find highest uncertainty variable
  u_1=0;
  i_1=-1;
  for (i=0; i<= 10; i++) {
    if ((i == 0) && (nle_state->all_uses.alpha_em == 1)) {
      rel_error=nle_config->ref_alpha_em_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i; 
      }
    } else if ((i == 1) && (nle_state->all_uses.v == 1)) {
      rel_error=nle_config->ref_v_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if ((i == 2) && (nle_state->all_uses.G == 1)) {
      rel_error=nle_config->ref_G_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if ((i == 3) && (nle_state->all_uses.mz == 1)) {
      rel_error=nle_config->ref_mz_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if ((i == 4) && (nle_state->all_uses.mw == 1)) {
      rel_error=nle_config->ref_mw_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if ((i == 5) && (nle_state->all_uses.mh0 == 1)) {
      rel_error=nle_config->ref_mh0_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if ((i == 6) && (nle_state->all_uses.sin2w == 1) && (nle_state->all_uses.mw_mz_mode == 0)) {
      rel_error=nle_config->ref_sin2w_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if ((i == 7) && (nle_state->all_uses.m_user == 1)) {
      rel_error=nle_config->smrfactor_mass_user_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if (i == 8) {
      rel_error=nle_config->ref_sm1_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if (i == 9) {
      rel_error=nle_config->ref_sm2_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } else if (i == 10) {
      rel_error=nle_config->ref_sm3_relerror;
      if (rel_error > u_1) {
        u_1=rel_error;
        i_1=i;
      }
    } // end for i
  }

  // find second highest uncertainty variable
  u_2=0;
  i_2=-1;
  for (i=0; i<= 10; i++) {
    if ((i == 0) && (nle_state->all_uses.alpha_em == 1)) {
      rel_error=nle_config->ref_alpha_em_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if ((i == 1) && (nle_state->all_uses.v == 1)) {
      rel_error=nle_config->ref_v_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if ((i == 2) && (nle_state->all_uses.G == 1)) {
      rel_error=nle_config->ref_G_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if ((i == 3) && (nle_state->all_uses.mz == 1)) {
      rel_error=nle_config->ref_mz_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if ((i == 4) && (nle_state->all_uses.mw == 1)) {
      rel_error=nle_config->ref_mw_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if ((i == 5) && (nle_state->all_uses.mh0 == 1)) {
      rel_error=nle_config->ref_mh0_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if ((i == 6) && (nle_state->all_uses.sin2w == 1) && (nle_state->all_uses.mw_mz_mode == 0)) {
      rel_error=nle_config->ref_sin2w_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if ((i == 7) && (nle_state->all_uses.m_user == 1)) {
      rel_error=nle_config->smrfactor_mass_user_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if (i == 8) {
      rel_error=nle_config->ref_sm1_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if (i == 9) {
      rel_error=nle_config->ref_sm2_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    } else if (i == 10) {
      rel_error=nle_config->ref_sm3_relerror;
      if ((rel_error > u_2) && (rel_error < u_1)) {
        u_2=rel_error;
        i_2=i;
      }
    }
  } // end for i

  // find third highest uncertainty variable
  u_3=0;
  i_3=-1;
  for (i=0; i<= 10; i++) {
    if ((i == 0) && (nle_state->all_uses.alpha_em == 1)) {
      rel_error=nle_config->ref_alpha_em_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if ((i == 1) && (nle_state->all_uses.v == 1)) {
      rel_error=nle_config->ref_v_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if ((i == 2) && (nle_state->all_uses.G == 1)) {
      rel_error=nle_config->ref_G_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if ((i == 3) && (nle_state->all_uses.mz == 1)) {
      rel_error=nle_config->ref_mz_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if ((i == 4) && (nle_state->all_uses.mw == 1)) {
      rel_error=nle_config->ref_mw_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if ((i == 5) && (nle_state->all_uses.mh0 == 1)) {
      rel_error=nle_config->ref_mh0_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if ((i == 6) && (nle_state->all_uses.sin2w == 1) && (nle_state->all_uses.mw_mz_mode == 0)) {
      rel_error=nle_config->ref_sin2w_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if ((i == 7) && (nle_state->all_uses.m_user == 1)) {
      rel_error=nle_config->smrfactor_mass_user_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if (i == 8) {
      rel_error=nle_config->ref_sm1_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if (i == 9) {
      rel_error=nle_config->ref_sm2_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    } else if (i == 10) {
      rel_error=nle_config->ref_sm3_relerror;
      if ((rel_error > u_3) && (rel_error < u_2)) {
        u_3=rel_error;
        i_3=i;
      }
    }
  } // end for i

  // set float flag on top three highest uncertainty variables
  // include sin2w only if mw_mz_mode == 0
  unknowns=0;
  if ((i_1 == 0) || (i_2 == 0) || (i_3 == 0)) {
    nle_state->all_uses.float_alpha_em=1;
    unknowns++;
  }
  if ((i_1 == 1) || (i_2 == 1) || (i_3 == 1)) {
    nle_state->all_uses.float_v=1;
    unknowns++;
  }
  if ((i_1 == 2) || (i_2 == 2) || (i_3 == 2)) {
    nle_state->all_uses.float_G=1;
    unknowns++;
  }
  if ((i_1 == 3) || (i_2 == 3) || (i_3 == 3)) {
    nle_state->all_uses.float_mz=1;
    unknowns++;
  }
  if ((i_1 == 4) || (i_2 == 4) || (i_3 == 4)) {
    nle_state->all_uses.float_mw=1;
    unknowns++;
  }
  if ((i_1 == 5) || (i_2 == 5) || (i_3 == 5)) {
    nle_state->all_uses.float_mh0=1;
    unknowns++;
  }
  if (((i_1 == 6) || (i_2 == 6) || (i_3 == 6)) && (nle_state->all_uses.mw_mz_mode == 0)) {
    nle_state->all_uses.float_sin2w=1;
    unknowns++;
  }
  if ((i_1 == 7) || (i_2 == 7) || (i_3 == 7)) {
    nle_state->all_uses.float_muser=1;
    unknowns++;
  }
  if ((i_1 == 8) || (i_2 == 8) || (i_3 == 8)) {
    nle_state->all_uses.float_sm1=1;
    unknowns++;
  }
  if ((i_1 == 9) || (i_2 == 9) || (i_3 == 9)) {
    nle_state->all_uses.float_sm2=1;
    unknowns++;
  }
  if ((i_1 == 10) || (i_2 == 10) || (i_3 == 10)) {
    nle_state->all_uses.float_sm3=1;
    unknowns++;
  }

  if (nle_config->nle_mode == 2) {
    // special checks for 2-term mode

    reference_mass=nle_state->term1.smrfactor_mass;
    if ((i_1 != reference_mass) && (i_2 != reference_mass) && (i_3 != reference_mass)) {
      // if mode == 2 and reference_mass is not floated, unset third unknown and set reference_mass in it's place
      // 2-term mixed will not converge unless the reference mass is floated
      if (i_3 == 0) {
        nle_state->all_uses.float_alpha_em=0;
        unknowns--;
      } else if (i_3 == 1) {
        nle_state->all_uses.float_v=0;
        unknowns--;
      } else if (i_3 == 2) {
        nle_state->all_uses.float_G=0;
        unknowns--;
      } else if (i_3 == 3) {
        nle_state->all_uses.float_mz=0;
        unknowns--;
      } else if (i_3 == 4) {
        nle_state->all_uses.float_mw=0;
        unknowns--;
      } else if (i_3 == 5) {
        nle_state->all_uses.float_mh0=0;
        unknowns--;
      } else if (i_3 == 6) {
        nle_state->all_uses.float_sin2w=0;
        unknowns--;
      } else if (i_3 == 7) {
        nle_state->all_uses.float_muser=0;
        unknowns--;
      } else if (i_3 == 8) {
        nle_state->all_uses.float_sm1=0;
        unknowns--;
      } else if (i_3 == 9) {
        nle_state->all_uses.float_sm2=0;
        unknowns--;
      } else if (i_3 == 10) {
        nle_state->all_uses.float_sm3=0;
        unknowns--;
      }
      i_3=reference_mass;
      if (reference_mass == 0) {
        nle_state->all_uses.float_alpha_em=1;
      } else if (reference_mass == 1) {
        nle_state->all_uses.float_v=1;
      } else if (reference_mass == 2) {
        nle_state->all_uses.float_G=1;
      } else if (reference_mass == 3) {
        nle_state->all_uses.float_mz=1;
      } else if (reference_mass == 4) {
        nle_state->all_uses.float_mw=1;
      } else if (reference_mass == 5) {
        nle_state->all_uses.float_mh0=1;
      } else if (reference_mass == 6) {
        nle_state->all_uses.float_sin2w=1;
      } else if (reference_mass == 7) {
        nle_state->all_uses.float_muser=1;
      } else if (reference_mass == 8) {
        nle_state->all_uses.float_sm1=1;
      } else if (reference_mass == 9) {
        nle_state->all_uses.float_sm2=1;
      } else if (reference_mass == 10) {
        nle_state->all_uses.float_sm3=1;
      }
      unknowns++;
    }

    if (nle_state->all_uses.float_sm3 != 1) {
      // if mode == 2 and sm3 is not floated, unset second unknown and set sm3 in it's place
      // 2-term mixed will not converge unless one of the solution masses is floated
      if (i_2 == 0) {
        nle_state->all_uses.float_alpha_em=0;
        unknowns--;
      } else if (i_2 == 1) {
        nle_state->all_uses.float_v=0;
        unknowns--;
      } else if (i_2 == 2) {
        nle_state->all_uses.float_G=0;
        unknowns--;
      } else if (i_2 == 3) {
        nle_state->all_uses.float_mz=0;
        unknowns--;
      } else if (i_2 == 4) {
        nle_state->all_uses.float_mw=0;
        unknowns--;
      } else if (i_2 == 5) {
        nle_state->all_uses.float_mh0=0;
        unknowns--;
      } else if (i_2 == 6) {
        nle_state->all_uses.float_sin2w=0;
        unknowns--;
      } else if (i_2 == 7) {
        nle_state->all_uses.float_muser=0;
        unknowns--;
      } else if (i_2 == 8) {
        nle_state->all_uses.float_sm1=0;
        unknowns--;
      } else if (i_2 == 9) {
        nle_state->all_uses.float_sm2=0;
        unknowns--;
      }
      i_2=10;    
      nle_state->all_uses.float_sm3=1;
      unknowns++;
    }
  }

#ifdef DEBUG_SELECT_OUTPUTS
  printf("debug, u_1: %.3e, i_1: %d, u_2: %.3e, i_2: %d, u_3: %.3e, i_3: %3d, mw_mz_mode: %d\n", u_1, i_1, u_2, i_2, u_3, i_3, nle_state->all_uses.mw_mz_mode);
  fflush(stdout);
#endif
  return(unknowns);
}
