#include <stdio.h>
#include <math.h>
#include "nle-lepton.h"

void initUses(input_use *uses) {
  uses->alpha_em=0;
  uses->v=0;
  uses->G=0;
  uses->mz=0;
  uses->mw=0;
  uses->mh0=0;
  uses->sin2w=0;
}

void addUses(input_use *dest, input_use *src) {
    dest->alpha_em= dest->alpha_em || src->alpha_em;
    dest->v=               dest->v || src->v;
    dest->G=               dest->G || src->G;
    dest->mz=             dest->mz || src->mz;
    dest->mw=             dest->mw || src->mw;
    dest->mh0=           dest->mh0 || src->mh0;
    dest->sin2w=       dest->sin2w || src->sin2w;
}

void printUses(input_use *uses) {
  printf("debug, uses:\n");
  printf("debug, ------------------------------\n");
  printf("debug, alpha_em: %d\n", uses->alpha_em);
  printf("debug, v:        %d\n", uses->v);
  printf("debug, G:        %d\n", uses->G);
  printf("debug, mz:       %d\n", uses->mz);
  printf("debug, mw:       %d\n", uses->mw);
  printf("debug, mh0:      %d\n", uses->mh0);
  printf("debug, sin2w:    %d\n", uses->sin2w);
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

/*
    if ((v == 5) || (v == 7) || (v == 9) || (v == 10) || (v == 11)) {
        return -1;
    }
*/

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

// from https://en.wikipedia.org/wiki/Binary_GCD_algorithm
unsigned int gcd2(unsigned int u, unsigned int v)
{
    // simple cases (termination)
    if (u == v)
        return u;

    if (u == 0)
        return v-1;

    if (v == 0)
        return u;

/*
    if ((v == 5) || (v == 7) || (v == 9) || (v == 10) || (v == 11)) {
        return -1;
    }
*/

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

void checkSymmetry(int *symmetry, int left, int middle, int right) {

  if (left == middle) {
    *symmetry+=2;
  }
  if (middle == right) {
    *symmetry+=2;
  }
  if (left == right) {
    *symmetry+=2;
  }
  if (left == -middle) {
    *symmetry+=1;
  }
  if (middle == -right) {
    *symmetry+=1;
  }
  if (left == -right) {
    *symmetry+=1;
  }
}
int interesting(int range, double inld) {
  double rangelow;
  double rangehigh;
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
  if ((inld < 0.620) || (inld > 16.1)) {  // 27=0.037 / 64=0.0156 // 9=0.1 // 16=.0620 // 32=0.0311
    return(0);
  }

  // test if int > 4 are not divisible by 2 or 3
  testint=(int)(inld + 0.5);
  if ((testint > 4) && ((testint % 2) != 0) && ((testint %3) != 0)) {
    return(0);
  }

  // test if close to int
  testld=fmodl(inld, 1.0); 
  if ((testld >= rangelow) && (testld <=rangehigh)) {
    return(1);
  }

  if  (inld < 1.0) {
    //  test if 1/int > 4 are divisible by 2 or 3
    testint=(int)((1.0 / inld) + 0.5);
    if ((testint > 4) && ((testint % 2) != 0) && ((testint %3) != 0)) {
      return(0);
    }

    // test if small number matches 1/int
    testld=fmodl((1.0 / inld), 1.0);
    if ((testld >= rangelow) && (testld <=rangehigh)) {
      return(1);
    }
  }

  return(0);
}
