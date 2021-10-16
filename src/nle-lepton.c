/*
 * Copyright (c) 2019, Kevin M. Loch
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "nle-lepton.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "util.h"
#include "nle-config.h"
#include "usage.h"
#include "initInfactorArray.h"
#include "initOutfactorArray.h"
#include "initSmrfactorArray.h"
#include "phase1.h"
#include "verifyMatches.h"
#include "generateExponents.h"

//#define DEBUG_SMRFACTOR // prints smrfactor_seq and optionally select specific smrfactor_seq
//#define DEBUG_SM3       // to scan for predicted sm3 once an interesting two-term test has been found

int processCmdArgs(nle_config_t *nle_config, int argc, char **argv) {
  int i;
  char *option_start;

  if (argc == 1) {
    return(0);
  } else {
    for (i=1; i <= (argc - 1); i++) {
      if (argv[i][1] == 'c') {
        // configuration file name
        if (argv[i][2] != 0) {
          // option concatenated onto switch
          option_start=argv[i];
          strcpy(nle_config->config_file_name, (option_start + (size_t)2)); 
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          strcpy(nle_config->config_file_name, option_start);
        } // end if no space
      } else if (argv[i][1] == 'h') {
        // print help
        printUsage();
        exit(0);
      } else if (argv[i][1] == 's') {
        if (argv[i][2] != 0) {
          // option concatenated onto switch
          option_start=argv[i];
          nle_config->external_seed=strtol((option_start + (size_t)2), NULL, 10);
        } else if ((argc >= (i + 1)) && (argv[i + 1][0] != '-')) {
          // option is probably next argv
          option_start=argv[i + 1];
          nle_config->external_seed=strtol(option_start, NULL, 10);
        } // end if no space
      } // end which option
    } // end for argc
  } // end if any options
  return(0);
}

void checkConfig(nle_config_t *nle_config) {
  // various checks on configuration file that will cause immediate exit if failed
  int i;

  // Check operating mode
  if ((nle_config->nle_mode != 3) && (nle_config->nle_mode != 2)) {
    printf("init, Error: in nle-lepton.cfg, nle_mode=%d is unsupported.  2 and 3 are the only supported modes in this version.\n", nle_config->nle_mode);
    fflush(stdout);
    exit(1);
  }

  // Check 2-term and 1-minus are used together
  if ((nle_config->smrfactor_1minus_enable == 1) && (nle_config->nle_mode != 2)) {
    printf("init, Error: in nle-lepton.cfg, nle_mode must be set to 2 when smrfactor_1minus_enable=yes\n");
    fflush(stdout);
    exit(1);
  }

  // Check if forced exponents exceed exp_inv_max
  if ((abs(nle_config->exp_inv_term1_force) > nle_config->exp_inv_max) || (abs(nle_config->exp_inv_term2_force) > nle_config->exp_inv_max) || (abs(nle_config->exp_inv_term2_force) > nle_config->exp_inv_max)) {
    printf("init, Error: in nle-lepton.cfg, forced exponents exceed exp_inv_max\n");
    fflush(stdout);
    exit(1);
  }

  // Check if all forced exponents or none
  i=0;
  if (nle_config->exp_inv_term1_force != 0) {
    i++;
  }
  if (nle_config->exp_inv_term2_force != 0) {
    i++;
  }
  if ((i > 0) && (i != 2) && (nle_config->nle_mode == 2)) {
    printf("init, Error: in nle-lepton.cfg, when nle_mode=2 term1 and term2 exponents must be forced or all set to random (exp_inv_term1_force...)\n");
    fflush(stdout);
    exit(1);
  }
  if (nle_config->exp_inv_term3_force != 0) {
    i++;
  }
  if ((i > 0) && (i != 3) && (nle_config->nle_mode == 3)) {
    printf("init, Error: in nle-lepton.cfg, when nle_mode=3 term1, term2, and term3 exponents must be forced or all set to random (exp_inv_term1_force...)\n");
    fflush(stdout);
    exit(1);
  }
  if (nle_config->exp_inv_term4_force != 0) {
    i++;
  }
  if ((i > 0) && (i != 4) && (nle_config->nle_mode == 4)) {
    printf("init, Error: in nle-lepton.cfg, when nle_mode=4 all exponents must be forced or all set to random (exp_inv_term1_force...)\n");
    fflush(stdout);
    exit(1);
  }

  // Check if phase 2 and smrfactor_mass_random mode are both enabled
  if ((nle_config->phase2_enable == 1) && (nle_config->smrfactor_mass_user_random == 1)) {
    printf("init, Error: in nle-lepton.cfg, setting phase2_enable=yes while smrfactor_mass_user_random=yes is not supported.\n");
    fflush(stdout);
    exit(1);
  }

  // Check if phase 2 and smrfactor_mass_scan mode are both enabled
  if ((nle_config->phase2_enable == 1) && (nle_config->smrfactor_mass_user_scan == 1)) {
    printf("init, Error: in nle-lepton.cfg, setting phase2_enable=yes while smrfactor_mass_user_scan=yes is not supported.\n");
    fflush(stdout);
    exit(1);
  }

  // Check if smrfactor_1minus is diabled and smrfactor_mass_random mode is enabled
  if ((nle_config->smrfactor_1minus_enable == 0) && (nle_config->smrfactor_mass_user_random == 1)) {
    printf("init, Error: in nle-lepton.cfg, setting smrfactor_1minus_enable=no while smrfactor_mass_user_random=yes is not supported.\n");
    fflush(stdout);
    exit(1);
  }

  // Check if smrfactor_1minus is diabled and smrfactor_mass_scan mode is enabled
  if ((nle_config->smrfactor_1minus_enable == 0) && (nle_config->smrfactor_mass_user_scan == 1)) {
    printf("init, Error: in nle-lepton.cfg, setting smrfactor_1minus_enable=no while smrfactor_mass_user_scan=yes is not supported.\n");
    fflush(stdout);
    exit(1);
  }

  //Check if both smrfactor_gt_sm3 and smrfactor_lt_sm1 are enabled
  if ((nle_config->smrfactor_gt_sm3 == 1) && (nle_config->smrfactor_lt_sm1 == 1)) {
    printf("init, Error: in nle-lepton.cfg, smrfactor_gt_sm3 and smrfactor_lt_sm1 cannot both be enabled at the same time\n");
    fflush(stdout);
    exit(1);
  }

  // Warn if phase 2 is disabled
  if (nle_config->phase2_enable == 0) {
    printf("init, Warning, phase 2 processing is disabled!  Set phase2_enabled=yes in nle-lepton.cfg if you want to re-enable\n");
  }

}

int main(int argc, char **argv) {
  nle_config_t nle_config;
  nle_state_t nle_state;
  struct timespec t;
  long seed;
  double testrand;
  double r;
  int i;
  int run_once;
  int smrfactor_seq;
  nle_smrfactor_precomputed_t *smrfactors;
  int coefficients_matched;
  long seedsec;
  long seedus;
  int failed;
  int smr_reference_mass_id;
  int smr_reference_mass_id_max;
  int smr_mass_ratio_enabled;
  int run;
  int polarity_seq;
  int valid_polarity;
  int mass_config_seq;
  int valid_mass_config;
  double smrfactor_mass_user;
  int mass_user_exp;
  int mass_user_exp_min;
  int mass_user_exp_max;
  int mass_user_exp_range;
  int valid_random_user_mass;
#ifdef DEBUG_SM3
  double sm3;
#endif

  // initialize nle_config to default values
  initConfig(&nle_config);

  // process command line arguments and load command line options into nle_conifg (including external seed and config file name)
  processCmdArgs(&nle_config, argc, argv);

  // initialize pseudorandom number generator from external seed and clock
  clock_gettime(CLOCK_REALTIME, &t);
  seedsec=(t.tv_sec % 1000000000);
  seedus=(t.tv_nsec / 1000);
  seed=nle_config.external_seed ^ (seedsec + seedus);
  nle_state.pcg_state=((__uint128_t)seed) << 64;
  testrand=pcg_ldrand64(&nle_state);

  // print version, config file nameand random seed data
  printf("init, version: %s, external seed: %ld, seconds seed: %ld, microseconds seed: %ld, composite seed: %ld, first random number: %.9e\n", NLE_VERSION, nle_config.external_seed, seedsec, seedus, seed, testrand);
  fflush(stdout);

  // load config file
  if (loadConfig(&nle_config) == 1) {
    exit(1);
  }

  // run sanity checks on config
  checkConfig(&nle_config);

  // allocate memory for  matches arrays
  nle_state.phase1_matches_start = (nle_phase1_match_t *)malloc(5000000 * sizeof(nle_phase1_match_t));
  nle_state.term1.matches_start = (nle_phase1_match_t *)malloc(100000 * sizeof(nle_phase1_match_t));
  nle_state.term2.matches_start = (nle_phase1_match_t *)malloc(100000 * sizeof(nle_phase1_match_t));
  nle_state.term3.matches_start = (nle_phase1_match_t *)malloc(100000 * sizeof(nle_phase1_match_t));

  // allocate memory for and initialize precomupted smrfactor array if smrfactor_1minus_enable == 1
  // otherwise set smrfactor count to zero
  if (nle_config.smrfactor_1minus_enable == 1) {
    nle_state.smrfactors_precomputed_start=(nle_smrfactor_precomputed_t *)malloc(1000000 * sizeof(nle_smrfactor_precomputed_t));
    initSmrfactorArray(&nle_config, &nle_state);
  } else {
    nle_state.smrfactors_precomputed_count=0;
  }

  // allocate memory for and initialize precomupted infactor array
  nle_state.infactors_precomputed_start=(nle_infactor_precomputed_t *)malloc(1000000 * sizeof(nle_infactor_precomputed_t));
  initInfactorArray(&nle_config, &nle_state);

  // allocate memory for and initialize precomputed outfactor array
  nle_state.outfactors_precomputed_start=(nle_outfactor_precomputed_t *)malloc(1000000 * sizeof(nle_outfactor_precomputed_t));
  initOutfactorArray(&nle_config, &nle_state);

  // main operating loop
  nle_state.phase1_seq=0;
  run=1;
  while (run) {
    if (nle_config.phase1_run_continuous == 0) {
      run=0; // only run once 
    }

    nle_state.phase1_seq++;

    // generate valid exponents
    generateExponents(&nle_config, &nle_state);

    if (nle_config.phase1_random_samples_enable == 1) {
      // generate random samples of experimental values each time phase 1 is run
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_sm1=((nle_config.ref_sm1 - nle_config.ref_sm1_error) + (r * 2.0 * nle_config.ref_sm1_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_sm2=((nle_config.ref_sm2 - nle_config.ref_sm2_error) + (r * 2.0 * nle_config.ref_sm2_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_sm3=((nle_config.ref_sm3 - nle_config.ref_sm3_error) + (r * 2.0 * nle_config.ref_sm3_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_v=((nle_config.ref_v - nle_config.ref_v_error) + (r * 2.0 * nle_config.ref_v_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_alpha_em=((nle_config.ref_alpha_em - nle_config.ref_alpha_em_error) + (r * 2.0 * nle_config.ref_alpha_em_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_alpha_w=((nle_config.ref_alpha_w - nle_config.ref_alpha_w_error) + (r * 2.0 * nle_config.ref_alpha_w_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_G=((nle_config.ref_G - nle_config.ref_G_error) + (r * 2.0 * nle_config.ref_G_error));
      nle_state.input_sample_mp=nle_config.ref_kg_to_ev * sqrt(nle_config.ref_hbar * nle_config.ref_c / nle_state.input_sample_G);
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_mz=((nle_config.ref_mz - nle_config.ref_mz_error) + (r * 2.0 * nle_config.ref_mz_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_mw=((nle_config.ref_mw - nle_config.ref_mw_error) + (r * 2.0 * nle_config.ref_mw_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_sin2w=((nle_config.ref_sin2w - nle_config.ref_sin2w_error) + (r * 2.0 * nle_config.ref_sin2w_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_mh0=((nle_config.ref_mh0 - nle_config.ref_mh0_error) + (r * 2.0 * nle_config.ref_mh0_error));
      r=pcg_ldrand64(&nle_state);
      nle_state.input_sample_muser=((nle_config.smrfactor_mass_user - nle_config.smrfactor_mass_user_error) + (r * 2.0 * nle_config.smrfactor_mass_user_error));
    } else {
      // set input_sample to center reference values
      nle_state.input_sample_sm1=nle_config.ref_sm1;
      nle_state.input_sample_sm2=nle_config.ref_sm2;
      nle_state.input_sample_sm3=nle_config.ref_sm3;
      nle_state.input_sample_v=nle_config.ref_v;
      nle_state.input_sample_alpha_em=nle_config.ref_alpha_em;
      nle_state.input_sample_alpha_w=nle_config.ref_alpha_w;
      nle_state.input_sample_G=nle_config.ref_G;
      nle_state.input_sample_mp=nle_config.ref_kg_to_ev * sqrt(nle_config.ref_hbar * nle_config.ref_c / nle_state.input_sample_G);
      nle_state.input_sample_mz=nle_config.ref_mz;
      nle_state.input_sample_mw=nle_config.ref_mw;
      nle_state.input_sample_sin2w=nle_config.ref_sin2w;
      nle_state.input_sample_mh0=nle_config.ref_mh0;
      nle_state.input_sample_muser=nle_config.smrfactor_mass_user;
    }

    // sequence through user mass range if smrfactor_mass_user and smrfactor_mass_user_scan are enabled and smrfactor_mass_user_random is disabled. Otherwise run once
    run_once=1;
    for (smrfactor_mass_user=nle_config.smrfactor_mass_user_min; (((nle_config.smrfactor_mass_user_enable == 1) && (nle_config.smrfactor_mass_user_scan == 1) && (nle_config.smrfactor_mass_user_random == 0) && (smrfactor_mass_user <= nle_config.smrfactor_mass_user_max)) || (run_once == 1)); smrfactor_mass_user+=(smrfactor_mass_user * nle_config.smrfactor_mass_user_step)) {
      run_once=0;

      if (nle_config.smrfactor_mass_user_enable == 1) {
        // set random or scan user mass if configured
        if (nle_config.smrfactor_mass_user_random == 1) {
          valid_random_user_mass=0;
          mass_user_exp_min=(int)log10(nle_config.smrfactor_mass_user_min);
          mass_user_exp_max=(int)log10(nle_config.smrfactor_mass_user_max) + 1;
          mass_user_exp_range=mass_user_exp_max - mass_user_exp_min;
          while(valid_random_user_mass == 0) {
            r=pcg_ldrand64(&nle_state);
            mass_user_exp=(int)(mass_user_exp_min + (r * mass_user_exp_range));
            r=pcg_ldrand64(&nle_state);
            smrfactor_mass_user=pow(10.0, (r + (double)mass_user_exp));
            if ((smrfactor_mass_user >= nle_config.smrfactor_mass_user_min) && (smrfactor_mass_user <= nle_config.smrfactor_mass_user_max)) {
              valid_random_user_mass=1;
            }
          }
          nle_state.input_sample_muser=smrfactor_mass_user;
        } else if (nle_config.smrfactor_mass_user_scan == 1) {
          nle_state.input_sample_muser=smrfactor_mass_user;
        }
      }

      // sequence through solution mass ratio reference masses if (1-smr) enabled, run once otherwise
      if (nle_config.smrfactor_1minus_enable == 1) {
        smr_reference_mass_id_max=5;
      } else {
        // run once
        smr_reference_mass_id_max=0;
      }
      for (smr_reference_mass_id=0; smr_reference_mass_id <= smr_reference_mass_id_max; smr_reference_mass_id++) {
        smr_mass_ratio_enabled=1;
        if (smr_reference_mass_id == 0) {
          if (nle_config.smrfactor_mass_mp_enable == 0) {
            smr_mass_ratio_enabled=0;
          }
        } else if (smr_reference_mass_id == 1) {
          if (nle_config.smrfactor_mass_v_enable == 0) {
            smr_mass_ratio_enabled=0;
          }
        } else if (smr_reference_mass_id == 2) {
          if (nle_config.smrfactor_mass_mz_enable == 0) {
            smr_mass_ratio_enabled=0;
          } 
        } else if (smr_reference_mass_id == 3) {
          if (nle_config.smrfactor_mass_mw_enable == 0) {
            smr_mass_ratio_enabled=0;
          }
        } else if (smr_reference_mass_id == 4) {
          if (nle_config.smrfactor_mass_mh0_enable == 0) {
            smr_mass_ratio_enabled=0;
          }
        } else if (smr_reference_mass_id == 5) {
          if (nle_config.smrfactor_mass_user_enable == 0) {
            smr_mass_ratio_enabled=0;
          }
        }                             
        if ((nle_config.smrfactor_1minus_enable == 0) || (smr_mass_ratio_enabled == 1)) {
          if (nle_config.smrfactor_1minus_enable == 1) {
              nle_state.term1.smrfactor_mass_id=smr_reference_mass_id;
              nle_state.term2.smrfactor_mass_id=smr_reference_mass_id;
          }
          // sequence through smrfactors if smrfactor_1minus is enabled, otherwise run once when smrfactor_seq == 0
          smrfactors=nle_state.smrfactors_precomputed_start;
          for (smrfactor_seq=0; smrfactor_seq <= nle_state.smrfactors_precomputed_count; smrfactor_seq++) {
            if ((nle_state.smrfactors_precomputed_count > 0) && (smrfactor_seq < nle_state.smrfactors_precomputed_count)) {
              nle_state.term1.current_smrfactors=smrfactors;
              nle_state.term2.current_smrfactors=smrfactors;
              nle_state.term1.smrfactor=smrfactors->smrfactor_multiplier;
              nle_state.term2.smrfactor=smrfactors->smrfactor_multiplier;
            }
            if ((nle_state.smrfactors_precomputed_count == 0) || (smrfactor_seq < nle_state.smrfactors_precomputed_count)) {
#ifdef DEBUG_SMRFACTOR
            //if (smrfactor_seq == 956) { // use to select specific smrfactors for debugging
              printf("debug, smrfactor_seq: %d\n", smrfactor_seq);
              fflush(stdout);
#endif
              // sequence through enabled mixing polarities, or run once if not (1-smr) mode   
              for (polarity_seq=0; polarity_seq <= 1; polarity_seq++) {
                // check polarity 
                valid_polarity=0;
                if ((nle_config.smrfactor_1minus_enable == 0) && (polarity_seq == 0)) { // run once f not (1-smr) mode
                  nle_state.nle_mixing_polarity=0; // not used if not (1-smr) but set anyway
                  valid_polarity=1;
                } else if (nle_config.smrfactor_1minus_enable == 1) {
                  if ((polarity_seq == 0) && ((nle_config.nle_mixing_polarity == -1) || (nle_config.nle_mixing_polarity == 0))) { 
                    nle_state.nle_mixing_polarity=0;
                    valid_polarity=1;
                  } else if ((polarity_seq == 1) && ((nle_config.nle_mixing_polarity == -1) || (nle_config.nle_mixing_polarity == 1))) {
                    nle_state.nle_mixing_polarity=1;
                    valid_polarity=1;
                  }
                } // end if 1-minus
                if (valid_polarity == 1) {

                  // sequence through enabled mass configurations, or run once if not (1-smr) mode
                  for (mass_config_seq=0; mass_config_seq <=1; mass_config_seq++) {
                    // check mass configuration 
                    valid_mass_config=0;
                    if ((nle_config.smrfactor_1minus_enable == 0) && (mass_config_seq == 0)) { // run once f not (1-smr) mode
                      nle_state.smrfactor_mass_configuration=1; // not used if not (1-smr) but set anyway
                      valid_mass_config=1;
                    } else if (nle_config.smrfactor_1minus_enable == 1) {
                      if ((mass_config_seq == 0) && ((nle_config.smrfactor_mass_configuration == -1) || (nle_config.smrfactor_mass_configuration == 0))) {
                        nle_state.smrfactor_mass_configuration=0;
                        valid_mass_config=1;
                      } else if ((mass_config_seq == 1) && ((nle_config.smrfactor_mass_configuration == -1) || (nle_config.smrfactor_mass_configuration == 1))) {
                        nle_state.smrfactor_mass_configuration=1;
                        valid_mass_config=1;
                      }
                    } // end if 1-minus
                    if (valid_mass_config == 1) {
#ifdef DEBUG_SM3
                    // this is for finding predicted sm3 once an interesting two-term test has been found.
                    for (sm3=1783.4089539E6; sm3 < 1785.862E6; sm3+=1.0E-3) {
                      nle_state.input_sample_sm3=sm3;
#endif

                      // set polarity and mass config strings used in status and other outputs
                      if (nle_state.nle_mixing_polarity == 0) {
                        nle_state.nle_mixing_polarity_str[0]='-';
                      } else if (nle_state.nle_mixing_polarity == 1) {
                        nle_state.nle_mixing_polarity_str[0]='+';
                      }
                      nle_state.nle_mixing_polarity_str[1]=0;

                      if (nle_state.smrfactor_mass_configuration == 0) {
                        sprintf(nle_state.smrfactor_mass_configuration_str, "mr/M");
                      } else if (nle_state.smrfactor_mass_configuration == 1) {
                        sprintf(nle_state.smrfactor_mass_configuration_str, "M/mr");
                      }
                      nle_state.smrfactor_mass_configuration_str[4]=0;

                      // clear terms_matched
                      for (i=0; i<=2; i++) {
                        nle_state.terms_matched[i]=0;
                      }
                      failed=solveNLEforCoefficients(&nle_config, &nle_state);
                      if (failed == 0) {
                        if (nle_state.phase1_matches_count > 0) {
                          coefficients_matched=0;
                          for (i=0; i <= 2; i++) {
                            if (nle_state.terms_matched[i] != 0) {
                              coefficients_matched++;
                            }
                          }
                          if (coefficients_matched == 3) {
                            // phase 2
                            verifyMatches(&nle_config, &nle_state);
                          } else {
                            if (nle_config.phase1_status_enable == 1) {
                              printf("status, No complete three-term phase 2 formulas to solve, terms with matches: %d, %d, %d\n", nle_state.terms_matched[0], nle_state.terms_matched[1], nle_state.terms_matched[2]);
                              fflush(stdout);
                            }
                          }
                        } else {
                          if (nle_config.phase1_status_enable == 1) {
                            printf("status, No interesting coefficient multipliers found\n");
                            fflush(stdout);
                          }
                        } // end nummatches
                      } // end if failed
#ifdef DEBUG_SM3
                    }
#endif
                    } // end if valid mass_config
                  } // end for mass_config_seq
                } // end if valid polarity
              } // end for polarity_seq
#ifdef DEBUG_SMRFACTOR
            //} // end if smrfactor_seq
#endif
              smrfactors++;
            } // end if smrfactor_seq
          } // end for smrfactor_seq
        } // end if smr_mass_ratio_enabled
      } // end for smr_reference_mass_id
    } // end for smrfactor_mass_user
  } // end while run
  exit(0);
}   
