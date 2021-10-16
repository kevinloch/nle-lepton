#include "nle-lepton.h"
#include <stdio.h>
#include <math.h>

// fast 64-bit random number generator to make full use of x86 80-bit extended prescision long doubles
// Adapted from https://www.pcg-random.org/posts/does-it-beat-the-minimal-standard.html by M. E. Oneil
long double pcg_ldrand64(nle_state_t *nle_state) {
    const __uint128_t MULTIPLIER = (((__uint128_t)0x2d99787926d46932) << 64) | ((__uint128_t)0xa4c1f32680f70c55);
    nle_state->pcg_state *= MULTIPLIER;
    nle_state->pcg_state += MULTIPLIER;
    return (long double)(nle_state->pcg_state >> 64) / (long double)((__uint128_t)1 << 64);
}

void initUses(nle_input_use_t *uses) {
  uses->G=0;
  uses->float_G=0;
  uses->v=0;
  uses->float_v=0;
  uses->mz=0;
  uses->float_mz=0;
  uses->mw=0;
  uses->float_mw=0;
  uses->mh0=0;
  uses->float_mh0=0;
  uses->m_user=0;
  uses->float_muser=0;
  uses->sin2w=0;
  uses->float_sin2w=0;
  uses->float_sm1=0;
  uses->float_sm2=0;
  uses->float_sm3=0;
  uses->alpha_em=0;
  uses->float_alpha_em=0;
  uses->alpha_w=0;
  uses->float_alpha_w=0;
  uses->mw_mz_mode=0;
}

void addUses(nle_input_use_t *dest, nle_input_use_t *src) {
    dest->G=               dest->G || src->G;
    dest->v=               dest->v || src->v;
    dest->mz=             dest->mz || src->mz;
    dest->mw=             dest->mw || src->mw;
    dest->mh0=           dest->mh0 || src->mh0;
    dest->m_user=     dest->m_user || src->m_user;
    dest->sin2w=       dest->sin2w || src->sin2w;
    dest->alpha_em= dest->alpha_em || src->alpha_em;
    dest->alpha_w=   dest->alpha_w || src->alpha_w;
}

void printUses(nle_input_use_t *uses) {
  printf("debug, uses:\n");
  printf("debug, ------------------------------\n");
  printf("debug, G:        %d\n", uses->G);
  printf("debug, v:        %d\n", uses->v);
  printf("debug, mz:       %d\n", uses->mz);
  printf("debug, mw:       %d\n", uses->mw);
  printf("debug, mh0:      %d\n", uses->mh0);
  printf("debug, m_user:   %d\n", uses->m_user);
  printf("debug, sin2w:    %d\n", uses->sin2w);
  printf("debug, alpha_em: %d\n", uses->alpha_em);
  printf("debug, alpha_w:  %d\n", uses->alpha_w);
  printf("debug, ------------------------------\n");
}

void printInputSamples(nle_state_t *nle_state) {
  printf("debug, input samples (labeled as ref_* to copy/paste into config for debugging)\n");
  printf("debug, ------------------------------\n");
  printf("debug, ref_sm1=%.14e\n", nle_state->input_sample_sm1);
  printf("debug, ref_sm2=%.14e\n", nle_state->input_sample_sm2);
  printf("debug, ref_sm3=%.14e\n", nle_state->input_sample_sm3);
  printf("debug, ref_v=%.14e\n", nle_state->input_sample_v);
  printf("debug, ref_alpha_em=%.14e\n", nle_state->input_sample_alpha_em);
  printf("debug, ref_alpha_w=%.14e\n", nle_state->input_sample_alpha_w);
  printf("debug, ref_G=%.14e\n", nle_state->input_sample_G);
  printf("debug, ref_mp=%.14e\n", nle_state->input_sample_mp);
  printf("debug, ref_mz=%.14e\n", nle_state->input_sample_mz);
  printf("debug, ref_mw=%.14e\n", nle_state->input_sample_mw);
  printf("debug, ref_sin2w=%.14e\n", nle_state->input_sample_sin2w);
  printf("debug, ref_mh0=%.14e\n", nle_state->input_sample_mh0);
  printf("debug, ref_mass_user=%.14e\n", nle_state->input_sample_muser);
  printf("debug, ------------------------------\n");
}

char *underscore(char *str, int len) {
  int i;

  for (i=0; i<len; i++) {
    if (str[i] == 32) {
      str[i]=95;
    }
  }
  return(str);
}

// from https://en.wikipedia.org/wiki/Binary_GCD_algorithm
unsigned int gcd(unsigned int u, unsigned int v)
{
    // simple cases (termination)
    if (u == v)
        return u;

    if (u == 0)
        return v;

    if (v == 0)
        return u;

    // look for factors of 2
    if (~u & 1) // u is even
    {
        if (v & 1) // v is odd
            return gcd(u >> 1, v);
        else // both u and v are even
            return gcd(u >> 1, v >> 1) << 1;
    }

    if (~v & 1) // u is odd, v is even
        return gcd(u, v >> 1);

    // reduce larger argument
    if (u > v)
        return gcd((u - v) >> 1, v);

    return gcd((v - u) >> 1, u);
}

void checkSymmetry2(int *symmetry, int term1, int term2) {

  if (term1 == term2) {
    *symmetry+=2;
  }
  if (term1 == -term2) {
    *symmetry+=1;
  }
}

void checkSymmetry3(int *symmetry, int term1, int term2, int term3) {

  if (term1 == term2) {
    *symmetry+=2;
  }
  if (term2 == term3) {
    *symmetry+=2;
  }
  if (term1 == term3) {
    *symmetry+=2;
  }
  if (term1 == -term2) {
    *symmetry+=1;
  }
  if (term2 == -term3) {
    *symmetry+=1;
  }
  if (term1 == -term3) {
    *symmetry+=1;
  }
}

int interesting(int range, int max_int, int filter_int, double inld) {
  double offset=0.01;
  double max_remainder=0.02;
  double testld;
  int testint;

  if (range == 3) {
    offset=       0.001; 
    max_remainder=0.002;
  } else if (range == 4) {
    offset=       0.0001;
    max_remainder=0.0002;
  } else if (range == 5) {
    offset=       0.00001;
    max_remainder=0.00002;
  } else if (range == 6) {
    offset=       0.000001;
    max_remainder=0.000002;
/*
  } else if (range == 7) {
    offset=       0.0000001;
    max_remainder=0.0000002;
  } else if (range == 8) {
    offset=       0.00000001
    max_remainder=0.00000002;
*/
  }

  // test if too big or small
  if ((inld < (1.0 / (max_int + 0.1))) || (inld > (max_int + 0.1))) {
    return(0);
  }

  if (inld > 1.0) {
    if (filter_int == 1) {
      // test if int > 4 are not divisible by 2 or 3
      testint=(int)(inld + 0.498);
      if ((testint > 4) && ((testint % 2) != 0) && ((testint % 3) != 0)) {
        return(0);
      }
    }

    // test if close to int
    testld=fmod((inld + offset), 1.0); 
    if (testld <= max_remainder) {
      return(1);
    }
  } else {
    if (filter_int == 1) {
      //  test if 1/int > 4 are divisible by 2 or 3
      testint=(int)((1.0 / inld) + 0.498);
      if ((testint > 4) && ((testint % 2) != 0) && ((testint % 3) != 0)) {
        return(0);
      }
    }

    // test if small number matches 1/int
    testld=fmod(((1.0 / inld) + offset), 1.0);
    if (testld <= max_remainder) {
      return(1);
    }
  }

  return(0);
}
