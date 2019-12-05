#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "nle-lepton.h"
#include "util.h"
#include "phase2.h"

int verifyMatches(matches *matchstart, int *nummatches, char *exponents, int leftinvexp, int middleinvexp, int rightinvexp, int random_input_count, int minsymmetry, int maxcomplexity) {
  //  For formulas with interesting coefficients on all three exponent terms,
  //  separate the match table into a separate list for each exponent, then test all unique combinations of coefficients
  //  for accuracy by comparing the computed muon mass to it's experimental value.  Print results matching a minimum threshold of interest.
  int i,j;
  matches *match;
  matches *leftmatches;
  matches *middlematches;
  matches *rightmatches;
  matches *leftmatchptr;
  matches *middlematchptr;
  matches *rightmatchptr;
  matches *tmpmatchptr;
  int numleftmatches, nummiddlematches, numrightmatches;
  int l,m,r;
  int dupe;
  struct timespec starttime;
  struct timespec endtime;
  double elapsedtime;
  double precision;
  int tmpmatchup;
  int tmpmatchdown;
  int tmpmatchcomplexity;
  long long tmphash;
  long totalcombos;
  long combo;
  input_use leftuses;
  input_use middleuses;
  input_use rightuses;
  input_use alluses;
  int complexity;
  int symmetry;

  leftmatches = (matches *)malloc(100000 * sizeof(matches));
  middlematches = (matches *)malloc(100000 * sizeof(matches));
  rightmatches = (matches *)malloc(100000 * sizeof(matches));

  // extract all left term coefficients
  match=matchstart;
  leftmatchptr=leftmatches;
  numleftmatches=0;
  for (i=0; i<*nummatches; i++) {
    if (match->invexp == leftinvexp) {
     // determine integer/rational match value
      if (match->match > 1.0) {
       tmpmatchup=(int)(match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / match->match) + 0.5);
      }
      tmphash=(long long)match->massratio ^ ((long long)((((tmpmatchup / tmpmatchdown) * (1.0 / match->matchmult)) * (long long)1.0E9) + 0.5));
      tmpmatchcomplexity=(match->matchcomplexity + (tmpmatchup * match->downout) + (tmpmatchdown * match->upout));
      // search existing match table for dupes and see if we have lower complexity
      tmpmatchptr=leftmatches;
      dupe=0;
      for (j=0; j< numleftmatches; j++) {
        if (tmphash == tmpmatchptr->matchhash) {
          if (tmpmatchcomplexity < tmpmatchptr->matchcomplexity) {
            // replace
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->matchup=tmpmatchup;
            tmpmatchptr->matchdown=tmpmatchdown;
            tmpmatchptr->invexp=match->invexp;
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->upin=match->upin;
            tmpmatchptr->downin=match->downin;
            tmpmatchptr->piupin=match->piupin;
            tmpmatchptr->pidownin=match->pidownin;
            tmpmatchptr->aupin=match->aupin;
            tmpmatchptr->adownin=match->adownin;
            tmpmatchptr->e2upin=match->e2upin;
            tmpmatchptr->e2downin=match->e2downin;
            tmpmatchptr->nbvupin=match->nbvupin;
            tmpmatchptr->nbsupin=match->nbsupin;
            tmpmatchptr->upout=match->upout;
            tmpmatchptr->downout=match->downout;
            tmpmatchptr->piupout=match->piupout;
            tmpmatchptr->pidownout=match->pidownout;
            tmpmatchptr->aupout=match->aupout;
            tmpmatchptr->adownout=match->adownout;
            tmpmatchptr->e2upout=match->e2upout;
            tmpmatchptr->e2downout=match->e2downout;
            tmpmatchptr->s2wupout=match->s2wupout;
            tmpmatchptr->s2wdownout=match->s2wdownout;
            tmpmatchptr->c2wupout=match->c2wupout;
            tmpmatchptr->c2wdownout=match->c2wdownout;
            tmpmatchptr->matchcomplexity=tmpmatchcomplexity;
            tmpmatchptr->matchmult=match->matchmult;
            tmpmatchptr->match=match->match;
            tmpmatchptr->matchhash=tmphash;
            initUses(&tmpmatchptr->uses);
            addUses(&tmpmatchptr->uses, &match->uses);
          }
          dupe=1;
          break;
        }
        tmpmatchptr++;
      } // end for j
      if (dupe ==0) {
        leftmatchptr->massratio=match->massratio;
        leftmatchptr->matchup=tmpmatchup;
        leftmatchptr->matchdown=tmpmatchdown;
        leftmatchptr->invexp=match->invexp;
        leftmatchptr->massratio=match->massratio;
        leftmatchptr->upin=match->upin;
        leftmatchptr->downin=match->downin;
        leftmatchptr->piupin=match->piupin;
        leftmatchptr->pidownin=match->pidownin;
        leftmatchptr->aupin=match->aupin;
        leftmatchptr->adownin=match->adownin;
        leftmatchptr->e2upin=match->e2upin;
        leftmatchptr->e2downin=match->e2downin;
        leftmatchptr->nbvupin=match->nbvupin;
        leftmatchptr->nbsupin=match->nbsupin;
        leftmatchptr->upout=match->upout;
        leftmatchptr->downout=match->downout;
        leftmatchptr->piupout=match->piupout;
        leftmatchptr->pidownout=match->pidownout;
        leftmatchptr->aupout=match->aupout;
        leftmatchptr->adownout=match->adownout;
        leftmatchptr->e2upout=match->e2upout;
        leftmatchptr->e2downout=match->e2downout;
        leftmatchptr->s2wupout=match->s2wupout;
        leftmatchptr->s2wdownout=match->s2wdownout;
        leftmatchptr->c2wupout=match->c2wupout;
        leftmatchptr->c2wdownout=match->c2wdownout;
        leftmatchptr->matchcomplexity=tmpmatchcomplexity;
        leftmatchptr->matchmult=match->matchmult;
        leftmatchptr->match=match->match;
        leftmatchptr->matchhash=tmphash;
        initUses(&tmpmatchptr->uses);
        addUses(&tmpmatchptr->uses, &match->uses);
        numleftmatches++;
        leftmatchptr++;
      } // end if not dupe
    }  // end if invexp
    match++;
  } // end for i

  // extract all middle term coefficients
  match=matchstart;
  middlematchptr=middlematches;
  nummiddlematches=0;
  for (i=0; i<*nummatches; i++) {
    if (match->invexp == middleinvexp) {
     // determine integer/rational match value
      if (match->match > 1.0) {
       tmpmatchup=(int)(match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / match->match) + 0.5);
      }
      tmphash=(long long)match->massratio ^ ((long long)((((tmpmatchup / tmpmatchdown) * (1.0 / match->matchmult)) * (long long)1.0E9) + 0.5));
      tmpmatchcomplexity=(match->matchcomplexity + (tmpmatchup * match->downout) + (tmpmatchdown * match->upout));
      // search existing match table for dupes and see if we have lower complexity
      tmpmatchptr=middlematches;
      dupe=0;
      for (j=0; j< nummiddlematches; j++) {
        if (tmphash == tmpmatchptr->matchhash) {
          if (tmpmatchcomplexity < tmpmatchptr->matchcomplexity) {
            // replace
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->matchup=tmpmatchup;
            tmpmatchptr->matchdown=tmpmatchdown;
            tmpmatchptr->invexp=match->invexp;
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->upin=match->upin;
            tmpmatchptr->downin=match->downin;
            tmpmatchptr->piupin=match->piupin;
            tmpmatchptr->pidownin=match->pidownin;
            tmpmatchptr->aupin=match->aupin;
            tmpmatchptr->adownin=match->adownin;
            tmpmatchptr->e2upin=match->e2upin;
            tmpmatchptr->e2downin=match->e2downin;
            tmpmatchptr->nbvupin=match->nbvupin;
            tmpmatchptr->nbsupin=match->nbsupin;
            tmpmatchptr->upout=match->upout;
            tmpmatchptr->downout=match->downout;
            tmpmatchptr->piupout=match->piupout;
            tmpmatchptr->pidownout=match->pidownout;
            tmpmatchptr->aupout=match->aupout;
            tmpmatchptr->adownout=match->adownout;
            tmpmatchptr->e2upout=match->e2upout;
            tmpmatchptr->e2downout=match->e2downout;
            tmpmatchptr->s2wupout=match->s2wupout;
            tmpmatchptr->s2wdownout=match->s2wdownout;
            tmpmatchptr->c2wupout=match->c2wupout;
            tmpmatchptr->c2wdownout=match->c2wdownout;
            tmpmatchptr->matchcomplexity=tmpmatchcomplexity;
            tmpmatchptr->matchmult=match->matchmult;
            tmpmatchptr->match=match->match;
            tmpmatchptr->matchhash=tmphash;
            initUses(&tmpmatchptr->uses);
            addUses(&tmpmatchptr->uses, &match->uses);
          }
          dupe=1;
          break;
        }
        tmpmatchptr++;
      } // end for j
      if (dupe ==0) {
        middlematchptr->massratio=match->massratio;
        middlematchptr->matchup=tmpmatchup;
        middlematchptr->matchdown=tmpmatchdown;
        middlematchptr->invexp=match->invexp;
        middlematchptr->massratio=match->massratio;
        middlematchptr->upin=match->upin;
        middlematchptr->downin=match->downin;
        middlematchptr->piupin=match->piupin;
        middlematchptr->pidownin=match->pidownin;
        middlematchptr->aupin=match->aupin;
        middlematchptr->adownin=match->adownin;
        middlematchptr->e2upin=match->e2upin;
        middlematchptr->e2downin=match->e2downin;
        middlematchptr->nbvupin=match->nbvupin;
        middlematchptr->nbsupin=match->nbsupin;
        middlematchptr->upout=match->upout;
        middlematchptr->downout=match->downout;
        middlematchptr->piupout=match->piupout;
        middlematchptr->pidownout=match->pidownout;
        middlematchptr->aupout=match->aupout;
        middlematchptr->adownout=match->adownout;
        middlematchptr->e2upout=match->e2upout;
        middlematchptr->e2downout=match->e2downout;
        middlematchptr->s2wupout=match->s2wupout;
        middlematchptr->s2wdownout=match->s2wdownout;
        middlematchptr->c2wupout=match->c2wupout;
        middlematchptr->c2wdownout=match->c2wdownout;
        middlematchptr->matchcomplexity=tmpmatchcomplexity;
        middlematchptr->matchmult=match->matchmult;
        middlematchptr->match=match->match;
        middlematchptr->matchhash=tmphash;
        initUses(&tmpmatchptr->uses);
        addUses(&tmpmatchptr->uses, &match->uses);
        nummiddlematches++;
        middlematchptr++;
      } // end if not dupe
    }  // end if invexp
    match++;
  } // end for i

  // extract all right term coefficients
  match=matchstart;
  rightmatchptr=rightmatches;
  numrightmatches=0;
  for (i=0; i<*nummatches; i++) {
    if (match->invexp == rightinvexp) {
     // determine integer/rational match value
      if (match->match > 1.0) {
       tmpmatchup=(int)(match->match + 0.5);
       tmpmatchdown=1;
      } else {
        tmpmatchup=1;
        tmpmatchdown=(int)((1.0 / match->match) + 0.5);
      }
      tmphash=(long long)match->massratio ^ ((long long)((((tmpmatchup / tmpmatchdown) * (1.0 / match->matchmult)) * (long long)1.0E9) + 0.5));
      tmpmatchcomplexity=(match->matchcomplexity + (tmpmatchup * match->downout) + (tmpmatchdown * match->upout));
      // search existing match table for dupes and see if we have lower complexity
      tmpmatchptr=rightmatches;
      dupe=0;
      for (j=0; j< numrightmatches; j++) {
        if (tmphash == tmpmatchptr->matchhash) {
          if (tmpmatchcomplexity < tmpmatchptr->matchcomplexity) {
            // replace
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->matchup=tmpmatchup;
            tmpmatchptr->matchdown=tmpmatchdown;
            tmpmatchptr->invexp=match->invexp;
            tmpmatchptr->massratio=match->massratio;
            tmpmatchptr->upin=match->upin;
            tmpmatchptr->downin=match->downin;
            tmpmatchptr->piupin=match->piupin;
            tmpmatchptr->pidownin=match->pidownin;
            tmpmatchptr->aupin=match->aupin;
            tmpmatchptr->adownin=match->adownin;
            tmpmatchptr->e2upin=match->e2upin;
            tmpmatchptr->e2downin=match->e2downin;
            tmpmatchptr->nbvupin=match->nbvupin;
            tmpmatchptr->nbsupin=match->nbsupin;
            tmpmatchptr->upout=match->upout;
            tmpmatchptr->downout=match->downout;
            tmpmatchptr->piupout=match->piupout;
            tmpmatchptr->pidownout=match->pidownout;
            tmpmatchptr->aupout=match->aupout;
            tmpmatchptr->adownout=match->adownout;
            tmpmatchptr->e2upout=match->e2upout;
            tmpmatchptr->e2downout=match->e2downout;
            tmpmatchptr->s2wupout=match->s2wupout;
            tmpmatchptr->s2wdownout=match->s2wdownout;
            tmpmatchptr->c2wupout=match->c2wupout;
            tmpmatchptr->c2wdownout=match->c2wdownout;
            tmpmatchptr->matchcomplexity=tmpmatchcomplexity;
            tmpmatchptr->matchmult=match->matchmult;
            tmpmatchptr->match=match->match;
            tmpmatchptr->matchhash=tmphash;
            initUses(&tmpmatchptr->uses);
            addUses(&tmpmatchptr->uses, &match->uses);
          }
          dupe=1;
          break;
        }
        tmpmatchptr++;
      } // end for j
      if (dupe ==0) {
        rightmatchptr->massratio=match->massratio;
        rightmatchptr->matchup=tmpmatchup;
        rightmatchptr->matchdown=tmpmatchdown;
        rightmatchptr->invexp=match->invexp;
        rightmatchptr->massratio=match->massratio;
        rightmatchptr->upin=match->upin;
        rightmatchptr->downin=match->downin;
        rightmatchptr->piupin=match->piupin;
        rightmatchptr->pidownin=match->pidownin;
        rightmatchptr->aupin=match->aupin;
        rightmatchptr->adownin=match->adownin;
        rightmatchptr->e2upin=match->e2upin;
        rightmatchptr->e2downin=match->e2downin;
        rightmatchptr->nbvupin=match->nbvupin;
        rightmatchptr->nbsupin=match->nbsupin;
        rightmatchptr->upout=match->upout;
        rightmatchptr->downout=match->downout;
        rightmatchptr->piupout=match->piupout;
        rightmatchptr->pidownout=match->pidownout;
        rightmatchptr->aupout=match->aupout;
        rightmatchptr->adownout=match->adownout;
        rightmatchptr->e2upout=match->e2upout;
        rightmatchptr->e2downout=match->e2downout;
        rightmatchptr->s2wupout=match->s2wupout;
        rightmatchptr->s2wdownout=match->s2wdownout;
        rightmatchptr->c2wupout=match->c2wupout;
        rightmatchptr->c2wdownout=match->c2wdownout;
        rightmatchptr->matchcomplexity=tmpmatchcomplexity;
        rightmatchptr->matchmult=match->matchmult;
        rightmatchptr->match=match->match;
        rightmatchptr->matchhash=tmphash;
        initUses(&tmpmatchptr->uses);
        addUses(&tmpmatchptr->uses, &match->uses);
        numrightmatches++;
        rightmatchptr++;
      } // end if not dupe
    }  // end if invexp
    match++;
  } // end for i

  // send all combinations of left, middle, right to solution function
  leftmatchptr=leftmatches;
  combo=0;
#ifdef SHOWSTATUS
  totalcombos=(numleftmatches * nummiddlematches * numrightmatches);
#endif
#ifdef SHOWSTATUS
  printf("status, Solving phase 2 formulas for masses, random input: %d, exponents: %s,                 progress: total (0/%ld) left (0/%d) middle (0/%d) right (0/%d)\n", random_input_count, exponents, totalcombos, numleftmatches, nummiddlematches, numrightmatches);
  fflush(stdout);
#endif
  for (l=0; l<numleftmatches; l++) {
    initUses(&leftuses);
    addUses(&leftuses, &leftmatchptr->uses);
    middlematchptr=middlematches;
    for (m=0; m<nummiddlematches; m++) {
      initUses(&middleuses);
      addUses(&middleuses, &middlematchptr->uses);
      rightmatchptr=rightmatches;
      for (r=0; r<numrightmatches; r++) {
        combo++;
        complexity=leftmatchptr->matchcomplexity + middlematchptr->matchcomplexity + rightmatchptr->matchcomplexity;
        // calculate symmetry score.   This measures how many factors are identical or inverse identical between the left, middle and right terms
        symmetry=0;
        checkSymmetry(&symmetry, (leftmatchptr->matchup * leftmatchptr->downout), (middlematchptr->matchup * middlematchptr->downout), (rightmatchptr->matchup * rightmatchptr->downout));
        checkSymmetry(&symmetry, (leftmatchptr->matchdown * leftmatchptr->upout), (middlematchptr->matchdown * middlematchptr->upout), (rightmatchptr->matchdown * rightmatchptr->upout));
        checkSymmetry(&symmetry, leftmatchptr->e2upout, middlematchptr->e2upout, rightmatchptr->e2upout);
        checkSymmetry(&symmetry, leftmatchptr->piupout, middlematchptr->piupout, rightmatchptr->piupout);
        checkSymmetry(&symmetry, (leftmatchptr->aupout * leftmatchptr->adownout), (middlematchptr->aupout * middlematchptr->adownout), (rightmatchptr->aupout * rightmatchptr->adownout));
        checkSymmetry(&symmetry, (leftmatchptr->s2wupout * leftmatchptr->s2wdownout), (middlematchptr->s2wupout * middlematchptr->s2wdownout), (rightmatchptr->s2wupout * rightmatchptr->s2wdownout));
        checkSymmetry(&symmetry, (leftmatchptr->c2wupout * leftmatchptr->c2wdownout), (middlematchptr->c2wupout * middlematchptr->c2wdownout), (rightmatchptr->c2wupout * rightmatchptr->c2wdownout));
        checkSymmetry(&symmetry, leftmatchptr->upin, middlematchptr->upin, rightmatchptr->upin);
        checkSymmetry(&symmetry, leftmatchptr->downin, middlematchptr->downin, rightmatchptr->downin);
        checkSymmetry(&symmetry, leftmatchptr->nbvupin, middlematchptr->nbvupin, rightmatchptr->nbvupin);
        checkSymmetry(&symmetry, leftmatchptr->nbsupin, middlematchptr->nbsupin, rightmatchptr->nbsupin);
        checkSymmetry(&symmetry, leftmatchptr->e2upin, middlematchptr->e2upin, rightmatchptr->e2upin);
        checkSymmetry(&symmetry, leftmatchptr->piupin, middlematchptr->piupin, rightmatchptr->piupin);
        checkSymmetry(&symmetry, (leftmatchptr->aupin * leftmatchptr->adownin), (middlematchptr->aupin * middlematchptr->adownin), (rightmatchptr->aupin * rightmatchptr->adownin));
        if ((symmetry >= minsymmetry) && (complexity <= maxcomplexity)) {
         if ((leftmatchptr->nbvupin == middlematchptr->nbvupin) && (leftmatchptr->nbvupin == rightmatchptr-> nbvupin) && (leftmatchptr->nbsupin == middlematchptr->nbsupin) && (leftmatchptr->nbsupin == rightmatchptr-> nbsupin)) { // consistency check
            initUses(&rightuses);
            addUses(&rightuses, &rightmatchptr->uses);
            initUses(&alluses);
            addUses(&alluses, &leftuses);
            addUses(&alluses, &middleuses);
            addUses(&alluses, &rightuses);
            clock_gettime(CLOCK_REALTIME, &starttime);
            precision=solveNLEforMasses(exponents, leftinvexp, middleinvexp, rightinvexp, leftmatchptr, middlematchptr, rightmatchptr, &alluses, maxcomplexity);
#ifdef SHOWSTATUS
            clock_gettime(CLOCK_REALTIME, &endtime);
            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime.tv_sec - 1500000000) + ((double)starttime.tv_nsec) / 1.0E9);
            printf("status, Solved  phase 2 formula  for masses, random input: %d, exponents: %s, mass mode: %d%d%d, progress: total (%ld/%ld) left (%d/%d) middle (%d/%d) right (%d/%d), precision: %.3e, (%6.4fs)\n", random_input_count, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, combo, totalcombos, l+1, numleftmatches, m+1, nummiddlematches, r+1, numrightmatches, precision, elapsedtime);
            fflush(stdout);
#endif
          } // end consistency check
        } // end if symmetry and complexity
        rightmatchptr++;
      } // for r
      middlematchptr++;
    } // for m
    leftmatchptr++;
  } // for l
  free(leftmatches);
  free(middlematches);
  free(rightmatches);
  return(0);
}
