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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nle-lepton.h"
#include "util.h"
#include "nle-config.h"
#include "usage.h"
#include "initInfactorArray.h"
#include "initOutfactorArray.h"
#include "initSmrfactorArray.h"
#include "initRmrfactorArray.h"
#include "phase1.h"
#include "verifyMatches.h"
#include "generateExponents.h"

//#define DEBUG_SMRFACTOR // prints smrfactor_seq and optionally select specific smrfactor_seq
//#define DEBUG_RMRFACTOR // prints rmrfactor_seq and optionally select specific rmrfactor_seq

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

  // check operating mode
  if ((nle_config->nle_mode != 3) && (nle_config->nle_mode != 2)) {
    printf("init, Error: in nle-lepton.cfg, nle_mode=%d is unsupported.  2 and 3 are the only supported modes in this version.\n", nle_config->nle_mode);
    fflush(stdout);
    exit(1);
  }

  // check 2-term and 1-minus are used together
  if ((nle_config->smrfactor_1minus_enable == 1) && (nle_config->nle_mode != 2)) {
    printf("init, Error: in nle-lepton.cfg, nle_mode must be set to 2 when smrfactor_1minus_enable=yes\n");
    fflush(stdout);
    exit(1);
  }

  // check if forced exponents exceed exp_inv_max
  if ((abs(nle_config->exp_inv_term1_force) > nle_config->exp_inv_max) || (abs(nle_config->exp_inv_term2_force) > nle_config->exp_inv_max) || (abs(nle_config->exp_inv_term2_force) > nle_config->exp_inv_max)) {
    printf("init, Error: in nle-lepton.cfg, forced exponents exceed exp_inv_max\n");
    fflush(stdout);
    exit(1);
  }

  // check if all forced exponents or none
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

  // warn if phase 2 is disabled
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
  int smrfactor_seq;
  int smrfactor_random_id=0;
  int smrfactor_seq_max;
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
  int rmr_mass_id_max;
  int rmr_mass_id_up;
  int rmr_mass_id_down;
  int rmr_mass_ratio_enabled;
  int rmrfactor_seq;
  int rmrfactor_random_id=0;
  int rmrfactor_seq_max;
  nle_rmrfactor_precomputed_t *rmrfactors;

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

  // allocate memory for and initialize precomupted rmrfactor array if rmrfactor_1minus_enable == 1
  // otherwise set rmrfactor count to zero
  if (nle_config.rmrfactor_1minus_enable == 1) {
    nle_state.rmrfactors_precomputed_start=(nle_rmrfactor_precomputed_t *)malloc(1000000 * sizeof(nle_rmrfactor_precomputed_t));
    initRmrfactorArray(&nle_config, &nle_state);
  } else {
    nle_state.rmrfactors_precomputed_count=0;
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
    if (nle_config.phase1_random_samples_enable == 0) {
      run=0; // only run once if random samples is disabled
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

    // sequence through solution mass ratio reference masses if (1-smr) enabled, run once otherwise
    if (nle_config.smrfactor_1minus_enable == 1) {
      smr_reference_mass_id_max=5;
    } else {
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
        // select random smrfactors or sequence through smrfactors if smrfactor_1minus is enabled, otherwise run once when smrfactor_seq == 0
        if (nle_config.smrfactor_1minus_random == 1) {
          r=pcg_ldrand64(&nle_state);
          smrfactor_random_id=(int)((r * (double)nle_state.smrfactors_precomputed_count) + 0.5);
          smrfactors=nle_state.smrfactors_precomputed_start + smrfactor_random_id;
          smrfactor_seq_max=0;
        } else {
          smrfactors=nle_state.smrfactors_precomputed_start;
          smrfactor_seq_max=nle_state.smrfactors_precomputed_count;
        }
        for (smrfactor_seq=0; smrfactor_seq <= smrfactor_seq_max; smrfactor_seq++) {
          if ((nle_state.smrfactors_precomputed_count > 0) && (smrfactor_seq < nle_state.smrfactors_precomputed_count)) {
            nle_state.term1.current_smrfactors=smrfactors;
            nle_state.term2.current_smrfactors=smrfactors;
            nle_state.term1.smrfactor=smrfactors->smrfactor_multiplier;
            nle_state.term2.smrfactor=smrfactors->smrfactor_multiplier;
          }
          if ((nle_state.smrfactors_precomputed_count == 0) || (smrfactor_seq < nle_state.smrfactors_precomputed_count)) {
#ifdef DEBUG_SMRFACTOR
            //if (smrfactor_seq == 294) { // use to select specific smrfactors for debugging
            printf("debug, smrfactor_seq: %d, smr_reference_mass_id: %d, smrfactor_random_id: %d\n", smrfactor_seq, smr_reference_mass_id, smrfactor_random_id);
            fflush(stdout);
#endif

            // sequence through reference mass ratio masses if (1-rmr-smr) enabled, run once otherwise
            if (nle_config.rmrfactor_1minus_enable == 1) {
              rmr_mass_id_max=5;
            } else {
              rmr_mass_id_max=0;
            }
            for (rmr_mass_id_up=0; rmr_mass_id_up <= rmr_mass_id_max; rmr_mass_id_up++) {
              for (rmr_mass_id_down=0; rmr_mass_id_down <= rmr_mass_id_max; rmr_mass_id_down++) {
                rmr_mass_ratio_enabled=1;
                if ((rmr_mass_id_up == 0) || (rmr_mass_id_down == 0)) {
                  if (nle_config.rmrfactor_mass_mp_enable == 0) {
                    rmr_mass_ratio_enabled=0;
                  }
                }
                if ((rmr_mass_id_up == 1) || (rmr_mass_id_down == 1)) {
                  if (nle_config.rmrfactor_mass_v_enable == 0) {
                    rmr_mass_ratio_enabled=0;
                  }
                }
                if ((rmr_mass_id_up == 2) || (rmr_mass_id_down == 2)) {
                  if (nle_config.rmrfactor_mass_mz_enable == 0) {
                    rmr_mass_ratio_enabled=0;
                  }
                }
                if ((rmr_mass_id_up == 3) || (rmr_mass_id_down == 3)) {
                  if (nle_config.rmrfactor_mass_mw_enable == 0) {
                    rmr_mass_ratio_enabled=0;
                  }
                }
                if ((rmr_mass_id_up == 4) || (rmr_mass_id_down == 4)) {
                  if (nle_config.rmrfactor_mass_mh0_enable == 0) {
                    rmr_mass_ratio_enabled=0;
                  }
                }
                if ((rmr_mass_id_up == 5) || (rmr_mass_id_down == 5)) {
                  if (nle_config.rmrfactor_mass_user_enable == 0) {
                    rmr_mass_ratio_enabled=0;
                  }
                }
                // skip if same mass up/down
                if (rmr_mass_id_up == rmr_mass_id_down) {
                  rmr_mass_ratio_enabled=0;
                }
                // skip if rmr/smr sync is enabled and rmr_mass_id_down != smr_reference_mass_id
                if ((nle_config.rmrfactor_smrfactor_sync == 1) && (rmr_mass_id_down != smr_reference_mass_id)) {
                  rmr_mass_ratio_enabled=0;
                }
                if ((nle_config.rmrfactor_1minus_enable == 0) || (rmr_mass_ratio_enabled == 1)) {
                  if (nle_config.rmrfactor_1minus_enable == 1) {
                      nle_state.term1.rmrfactor_mass_id_up=rmr_mass_id_up;
                      nle_state.term1.rmrfactor_mass_id_down=rmr_mass_id_down;
                  }
                  // select random rmrfactors or sequence through rmrfactors if rmrfactor_1minus is enabled, otherwise run once when rmrfactor_seq == 0
                  if (nle_config.rmrfactor_1minus_random == 1) {
                    r=pcg_ldrand64(&nle_state);
                    rmrfactor_random_id=(int)((r * (double)nle_state.rmrfactors_precomputed_count) + 0.5);
                    rmrfactors=nle_state.rmrfactors_precomputed_start + rmrfactor_random_id;
                    rmrfactor_seq_max=0;
                  } else {
                    rmrfactors=nle_state.rmrfactors_precomputed_start;
                    rmrfactor_seq_max=nle_state.rmrfactors_precomputed_count;
                  }
                  for (rmrfactor_seq=0; rmrfactor_seq <= rmrfactor_seq_max; rmrfactor_seq++) {
                    if ((nle_state.rmrfactors_precomputed_count > 0) && (rmrfactor_seq < nle_state.rmrfactors_precomputed_count)) {
                      nle_state.term1.current_rmrfactors=rmrfactors;
                      nle_state.term1.rmrfactor=rmrfactors->rmrfactor_multiplier;
                    }
                    if ((nle_state.rmrfactors_precomputed_count == 0) || (rmrfactor_seq < nle_state.rmrfactors_precomputed_count)) {
#ifdef DEBUG_RMRFACTOR
                      //if (rmrfactor_seq == 103) { // use to select specific rmrfactors for debugging
                      printf("debug, rmrfactor_seq: %d, rmr_mass_id_up: %d, rmr_mass_id_down: %d, smrfactmr_seq: %d, smr_reference_mass_id: %d, rmrfactor_random_id: %d, smrfactor_random_id: %d\n", rmrfactor_seq, rmr_mass_id_up, rmr_mass_id_down, smrfactor_seq, smr_reference_mass_id, rmrfactor_random_id, smrfactor_random_id);
                      fflush(stdout);
#endif

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
                          // phase 1
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
                        } // end if valid polarity
                      } // end for polarity_seq

#ifdef DEBUG_RMRFACTOR
                      //} // end if rmrfactor_seq ==
#endif
                      if (nle_config.rmrfactor_1minus_random == 0) {
                        rmrfactors++;
                      }
                    } // end if rmrfactor_seq <
                  } // end for rmrfactor_seq
                } // end if rmr_mass_ratio_enabled
              } // end for rmr_mass_id_down
            } // end for rmr_mass_id_up

#ifdef DEBUG_SMRFACTOR
            //} // end if smrfactor_seq ==
#endif
            if (nle_config.smrfactor_1minus_random == 0) {
              smrfactors++;
            }
          } // end if smrfactor_seq <
        } // end for smrfactor_seq
      } // end if smr_mass_ratio_enabled
    } // end for smr_reference_mass_id

  } // end while run
  exit(0);
}   
