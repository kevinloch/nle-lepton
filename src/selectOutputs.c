#include <stdio.h>
#include "nle-lepton.h"
#include "util.h"
#include "math.h"

void selectHighestUncertaintyOutput(nle_config_t *nle_config, nle_state_t *nle_state, int *unchosen_outputs, int *unknowns) {

  int i;
  double u_1;
  int i_1;

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

  // find highest uncertainty variable from unchosen outputs
  u_1=0;
  i_1=-1;
  for (i=0; i < 12; i++) {
    if (unchosen_outputs[i] == 1) {
      if (i == 0) {
        // we are really evaluating mp so we check sqrt(G_relerror)
        if (sqrt(nle_config->relerror[i]) > u_1) {
          u_1=sqrt(nle_config->relerror[i]);
          i_1=i;
        }
      } else if (i == 9) {
        // only float sin2w if mw_mz_mode == 0
        if ((nle_state->all_uses.mw_mz_mode == 0) && (nle_config->relerror[i] > u_1)) {
          u_1=nle_config->relerror[i];
          i_1=i;
        }
      } else  if (nle_config->relerror[i] > u_1) {
        u_1=nle_config->relerror[i];
        i_1=i;
      }
    } // end if unchosen_outputs
  } // end for i

  // set selected output as chosen and floated
  if (i_1 == 0) {
    nle_state->all_uses.float_G=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 1) {
    nle_state->all_uses.float_v=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 2) {
    nle_state->all_uses.float_mz=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 3) {
    nle_state->all_uses.float_mw=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 4) {
    nle_state->all_uses.float_mh0=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 5) {
    nle_state->all_uses.float_muser=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 6) {
    nle_state->all_uses.float_sm1=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 7) {
    nle_state->all_uses.float_sm2=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 8) {
    nle_state->all_uses.float_sm3=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 9) {
    nle_state->all_uses.float_sin2w=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 10) {
    nle_state->all_uses.float_alpha_em=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
  if (i_1 == 11) {
    nle_state->all_uses.float_alpha_w=1;
    unchosen_outputs[i_1]=0;
    *unknowns=*unknowns+1;
  }
}

int selectOutputs(nle_config_t *nle_config, nle_state_t *nle_state) {
  // select the three used variables with the highest experimental uncertainty
  // return the number of unknowns (outputs to solve for in phase 2)

  int unknowns;
  int unchosen_outputs[12];
  int i;
  int reference_mass;

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

  // start with zero unknowns
  unknowns=0;

  // initialize list of unchosen outputs
  for (i=0; i < 12; i++) {
    if (i == 0) {
      if (nle_state->all_uses.G == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      }
    } else if (i == 1) {
      if (nle_state->all_uses.v == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      } 
    } else if (i == 2) {
      if (nle_state->all_uses.mz == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      } 
    } else if (i == 3) {
      if (nle_state->all_uses.mw == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      } 
    } else if (i == 4) {
      if (nle_state->all_uses.mh0 == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      } 
    } else if (i == 5) {
      if (nle_state->all_uses.m_user == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      } 
    } else if (i == 6) {
      unchosen_outputs[i]=1;
    } else if (i == 7) {
      unchosen_outputs[i]=1;
    } else if (i == 8) {
      unchosen_outputs[i]=1;
    } else if (i == 9) {
      if (nle_state->all_uses.sin2w == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      }
    } else if (i == 10) {
      if (nle_state->all_uses.alpha_em == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      }
    } else if (i == 11) {
      if (nle_state->all_uses.alpha_w == 1) {
        unchosen_outputs[i]=1;
      } else {
        unchosen_outputs[i]=0;
      }
    } // end if i
  } // end for i

  // 2-term mixed will not converge unless at least one reference mass is floated
  if (nle_config->nle_mode == 2) {
    if (nle_config->smrfactor_1minus_enable == 1) {
      reference_mass=nle_state->term1.smrfactor_mass;
    } else {
      if (nle_state->term1.current_match->smrfactor_mass == nle_state->term2.current_match->smrfactor_mass) {
        reference_mass=nle_state->term1.current_match->smrfactor_mass;
      } else {
        // find which term has the highest uncertainty reference mass
        if (nle_config->relerror[nle_state->term1.current_match->smrfactor_mass] > nle_config->relerror[nle_state->term2.current_match->smrfactor_mass]) {
          reference_mass=nle_state->term1.current_match->smrfactor_mass;
        } else {
          reference_mass=nle_state->term2.current_match->smrfactor_mass;
        }
      }
    }
    if ((reference_mass == 0) && (nle_state->all_uses.float_G == 0)) {
      nle_state->all_uses.float_G=1;
      unchosen_outputs[0]=0;
      unknowns++;
    } else if ((reference_mass == 1) && (nle_state->all_uses.float_v == 0)) {
      nle_state->all_uses.float_v=1;
      unchosen_outputs[1]=0;
      unknowns++;
    } else if ((reference_mass == 2) && (nle_state->all_uses.float_mz == 0)) {
      nle_state->all_uses.float_mz=1;
      unchosen_outputs[2]=0;
      unknowns++;
    } else if ((reference_mass == 3) && (nle_state->all_uses.float_mw == 0)) {
      nle_state->all_uses.float_mw=1;
      unchosen_outputs[3]=0;
      unknowns++;
    } else if ((reference_mass == 4) && (nle_state->all_uses.float_mh0 == 0)) {
      nle_state->all_uses.float_mh0=1;
      unchosen_outputs[4]=0;
      unknowns++;
    } else if ((reference_mass == 5) && (nle_state->all_uses.float_muser == 0)) {
      nle_state->all_uses.float_muser=1;
      unchosen_outputs[5]=0;
      unknowns++;
    } else if ((reference_mass == 6) && (nle_state->all_uses.float_sm1 == 0)) {
      nle_state->all_uses.float_sm1=1;
      unchosen_outputs[6]=0;
      unknowns++;
    } else if ((reference_mass == 7) && (nle_state->all_uses.float_sm2 == 0)) {
      nle_state->all_uses.float_sm2=1;
      unchosen_outputs[7]=0;
      unknowns++;
    } else if ((reference_mass == 7) && (nle_state->all_uses.float_sm3 == 0)) {
      nle_state->all_uses.float_sm3=1;
      unchosen_outputs[8]=0;
      unknowns++;
    }
  }

  // 2-term mixed will not converge unless at least one solution mass is floated
  if (nle_config->nle_mode == 2) {
    // find which solution mass has the highest relative uncertainty
    // this assumes all three have different relative uncertainties, or all three have the same relative uncertainty
    if ((nle_config->relerror[6] > nle_config->relerror[7]) && (nle_config->relerror[6] > nle_config->relerror[8])) {
      nle_state->all_uses.float_sm1=1;
      unchosen_outputs[6]=0;
      unknowns++;
    } else if ((nle_config->relerror[7] > nle_config->relerror[6]) && (nle_config->relerror[7] > nle_config->relerror[8])) {
      nle_state->all_uses.float_sm2=1;
      unchosen_outputs[7]=0;
      unknowns++;
    } else {
      nle_state->all_uses.float_sm3=1;
      unchosen_outputs[8]=0;
      unknowns++;
    }
  }

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

  // select unknowns until we have three total
  while (unknowns < 3) {
    selectHighestUncertaintyOutput(nle_config, nle_state, unchosen_outputs, &unknowns);
  }

  return(unknowns);
}
