#include <stdio.h>
#include <math.h>
#include "nle-lepton.h"

void initUses(nle_input_use_t *uses) {
  uses->alpha_em=0;
  uses->float_alpha_em=0;
  uses->v=0;
  uses->float_v=0;
  uses->G=0;
  uses->float_G=0;
  uses->mz=0;
  uses->float_mz=0;
  uses->mw=0;
  uses->float_mw=0;
  uses->mh0=0;
  uses->float_mh0=0;
  uses->sin2w=0;
  uses->float_sin2w=0;
  uses->m_user=0;
  uses->float_muser=0;
  uses->mw_mz_mode=0;
  uses->float_sm1=0;
  uses->float_sm2=0;
  uses->float_sm3=0;
}

void addUses(nle_input_use_t *dest, nle_input_use_t *src) {
    dest->alpha_em= dest->alpha_em || src->alpha_em;
    dest->v=               dest->v || src->v;
    dest->G=               dest->G || src->G;
    dest->mz=             dest->mz || src->mz;
    dest->mw=             dest->mw || src->mw;
    dest->mh0=           dest->mh0 || src->mh0;
    dest->sin2w=       dest->sin2w || src->sin2w;
    dest->m_user=     dest->m_user || src->m_user;
}

void printUses(nle_input_use_t *uses) {
  printf("debug, uses:\n");
  printf("debug, ------------------------------\n");
  printf("debug, alpha_em: %d\n", uses->alpha_em);
  printf("debug, v:        %d\n", uses->v);
  printf("debug, G:        %d\n", uses->G);
  printf("debug, mz:       %d\n", uses->mz);
  printf("debug, mw:       %d\n", uses->mw);
  printf("debug, mh0:      %d\n", uses->mh0);
  printf("debug, sin2w:    %d\n", uses->sin2w);
  printf("debug, m_user:   %d\n", uses->m_user);
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
  double rangelow=0.9;
  double rangehigh=1.1;
  double testld;
  int testint;

/*
  if (range == 2) {
    rangehigh=1.01;
    rangelow=0.99;
*/
  if (range == 3) {
    rangehigh=1.001;
    rangelow=0.999;
  } else if (range == 4) {
    rangehigh=1.0001;
    rangelow=0.9999;
  } else if (range == 5) {
    rangehigh=1.00001;
    rangelow=0.99999;
  } else if (range == 6) {
    rangehigh=1.000001;
    rangelow=0.999999;
/*
  } else if (range == 7) {
    rangehigh=1.0000001;
    rangelow=0.9999999;
  } else if (range == 8) {
    rangehigh=1.00000001;
    rangelow=0.99999999;
*/
  }

  // test if too big or small
  if ((inld < (1.0 / (max_int + 0.1))) || (inld > (max_int + 0.1))) {
    return(0);
  }

  if (filter_int == 1) {
    // test if int > 4 are not divisible by 2 or 3
    testint=(int)(inld + 0.498);
    if ((testint > 4) && ((testint % 2) != 0) && ((testint %3) != 0)) {
      return(0);
    }
  }

  // test if close to int
  testld=fmodl(inld, 1.0); 
  if ((testld >= rangelow) && (testld <= rangehigh)) {
    return(1);
  }

  if  ((filter_int == 1) && (inld < 1.0)) {
    //  test if 1/int > 4 are divisible by 2 or 3
    testint=(int)((1.0 / inld) + 0.498);
    if ((testint > 4) && ((testint % 2) != 0) && ((testint %3) != 0)) {
      return(0);
    }

    // test if small number matches 1/int
    testld=fmodl((1.0 / inld), 1.0);
    if ((testld >= rangelow) && (testld <= rangehigh)) {
      return(1);
    }
  }

  return(0);
}
