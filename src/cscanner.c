#include <stdio.h>
#include <stdlib.h> // abs
#include <math.h>   // pow
#include "nle-lepton.h"
#include "reference.h"
#include "util.h"

void cscanner(multipliers *multstart, int *nummult, matches **matchptr, int *nummatches, int *coffhit, int range, char *exponents, int leftinvexp, int middleinvexp, int rightinvexp, double cleft, double cmiddle, double cright, random_input random_inputs) {
  /*************/
  /* Phase 1.5 */
  /*************/
  //  Each coefficient is multiplied by various numbers and the resulting value is tested to see if it is
  // close to an interesting integer or simple rational number.  The results are stored in the match table.
  int i;
  double multiplierout;
  multipliers *mult;
  unsigned int u, v;
  matches *match;
  int complexity;

  int massratio;
  double massterm;
  double leftmassterm;
  double middlemassterm;
  double rightmassterm;
  int upout=1, downout=1;
  int piupout=0, pidownout=1;
  int aupout=1, adownout=1;
  int e2upout=0, e2downout=2;
  int s2wupout=0, s2wdownout=1;
  int c2wupout=0, c2wdownout=1;
  double updownout;
  double piout;
  double aout;
  double e2out;
  double s2wout=1.0;
  double c2wout=1.0;

  double sin2w;
  double cos2w;

  sin2w=random_inputs.sin2w_sample;
  cos2w=1.0 - sin2w;

  match=*matchptr;
  // here we substitute the M/me mass ratio used in phase 1 with the actual test mass ratio
  for (massratio=0; massratio<=4; massratio++) {
    if (massratio == 0) {         
      massterm=random_inputs.mp_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mp, ");             
      fflush(stdout);               
#endif
    } else if (massratio == 1) {  
      massterm=v_ref/me_ref;
#ifdef SHOWSTATUS
      printf("M/v, ");              
      fflush(stdout);               
#endif
    } else if (massratio == 2) {  
      massterm=random_inputs.mz_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mz, ");             
      fflush(stdout);               
#endif
    } else if (massratio == 3) {  
      massterm=random_inputs.mw_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mw, ");             
      fflush(stdout);               
#endif
    } else if (massratio == 4) {  
      massterm=random_inputs.mh0_sample/me_ref;
#ifdef SHOWSTATUS
      printf("M/mh0");             
      fflush(stdout);               
#endif
    }                             
    leftmassterm=pow(massterm, (1.0 / (double)leftinvexp));
    middlemassterm=pow(massterm, (1.0 / (double)middleinvexp));
    rightmassterm=pow(massterm, (1.0 / (double)rightinvexp));

    for (upout=1; upout <= 3; upout++) {
      for (downout=1; downout <=3; downout++) {
        u=upout;
        v=downout;
        if (gcd(u, v) == 1) {
          updownout=(double)upout / (double)downout;
          for (piupout=-2; piupout <= 2; piupout++) {
/*
            for (pidownout=1; pidownout <= 1; pidownout++) {
              u=abs(piupout);
              v=pidownout;
              if (gcd(u, v) == 1) {
*/
                piout=updownout * pow(M_PI, ((float)piupout / (float)pidownout));
                for (aupout=-2; aupout <= 2; aupout++) {
                  for (adownout=1; adownout <= 2; adownout++) {
                    u=abs(aupout);
                    v=adownout;
                    if (gcd(u, v) == 1) {
                      aout=piout * pow(alpha_ref, ((float)aupout / (float)adownout));
                      for (e2upout=-1; e2upout <= 1; e2upout++) {
                        e2out=aout * pow(2.0, ((float)e2upout / (float)e2downout));
// s2w and c2w are calculated separately since they contain significant uncertainty and are reconstructed separately in phase 2
#ifdef SIN2W
                        for (s2wupout=-1; s2wupout <= 1; s2wupout++) {
                          for (s2wdownout=1; s2wdownout <= 2; s2wdownout++) {
                            u=abs(s2wupout);
                            v=s2wdownout;
                            if ((gcd2(u, v) == 1) && ((v % 2) == 0)) {
                              s2wout=pow(sin2w, ((float)s2wupout / (float)s2wdownout));
                              for (c2wupout=-1; c2wupout <= 1; c2wupout++) {
                                for (c2wdownout=1; c2wdownout <= 2; c2wdownout++) {
                                  u=abs(c2wupout);
                                  v=c2wdownout;
                                  if ((gcd2(u, v) == 1) && ((v % 2) == 0)) {
                                    c2wout=pow(cos2w, ((float)c2wupout / (float)c2wdownout));
#endif
                                    mult=multstart;
                                    for (i=0; i<*nummult; i++) { 
                                      // test multiplier against left coefficient
                                      if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == leftinvexp) || (s2wdownout == (leftinvexp * 2))) {
                                        if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == leftinvexp) || (c2wdownout == (leftinvexp * 2))) {
                                          multiplierout=cleft * leftmassterm * s2wout * c2wout * e2out * mult->mult[abs(leftinvexp)];
                                          if (interesting(range, multiplierout)) {
                                            complexity=mult->multcomplexity\
                                                  + (upout + downout)\
                                                  + (abs(piupout) + pidownout)\
                                                  + (abs(aupout) + adownout)\
                                                  + (abs(e2upout) + e2downout);
/*
                                                  + (abs(s2wupout) + s2wdownout)\
                                                  + (abs(c2wupout) + c2wdownout);
*/
                                            coffhit[0]=leftinvexp;
                                            match->invexp=leftinvexp;
                                            match->massratio=massratio;
                                            match->upin=mult->upin;
                                            match->downin=mult->downin;
                                            match->piupin=mult->piupin;
                                            match->pidownin=mult->pidownin;
                                            match->aupin=mult->aupin;
                                            match->adownin=mult->adownin;
                                            match->e2upin=mult->e2upin;
                                            match->e2downin=mult->e2downin;
                                            match->nbvupin=mult->nbvupin;
                                            match->nbsupin=mult->nbsupin;
                                            match->upout=upout;
                                            match->downout=downout;
                                            match->piupout=piupout;
                                            match->pidownout=pidownout;
                                            match->aupout=aupout;
                                            match->adownout=adownout;
                                            match->e2upout=e2upout;
                                            match->e2downout=e2downout;
                                            match->s2wupout=s2wupout;
                                            match->s2wdownout=s2wdownout;
                                            match->c2wupout=c2wupout;
                                            match->c2wdownout=c2wdownout;
                                            match->matchcomplexity=complexity;
                                            match->match=multiplierout;
                                            initUses(&match->uses);
                                            addUses(&match->uses, &mult->uses);
                                            if (massratio == 0) {
                                              match->uses.G=1;
                                            } else if (massratio == 1) {
                                              match->uses.v=1;
                                            } else if (massratio == 2) {
                                              match->uses.mz=1;
                                            } else if (massratio == 3) {
                                              match->uses.mw=1;
                                            } else if (massratio == 4) {
                                              match->uses.mh0=1;
                                            }
                                            if (aupout != 0) {
                                              match->uses.alpha_em=1;
                                            }
                                            if ((s2wupout != 0) || (c2wupout != 0)) {
                                              match->uses.sin2w=1;
                                            }
                                            match->matchmult=e2out * mult->mult[abs(leftinvexp)];
                                            *nummatches=*nummatches+1;
                                            *matchptr=*matchptr+1;
                                            match=*matchptr;
                                          }  // if interesting left
                                        } // sanity check c2w left
                                      } // sanity check s2w left

                                      // test multiplier against middle coefficient
                                      if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == middleinvexp) || (s2wdownout == (middleinvexp * 2))) {
                                        if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == middleinvexp) || (c2wdownout == (middleinvexp * 2))) {
                                          multiplierout=cmiddle * middlemassterm * s2wout * c2wout * e2out * mult->mult[abs(middleinvexp)];
                                          if (interesting(range, multiplierout)) { 
                                            complexity=mult->multcomplexity\
                                                  + (upout + downout)\
                                                  + (abs(piupout) + pidownout)\
                                                  + (abs(aupout) + adownout)\
                                                  + (abs(e2upout) + e2downout);
/*
                                                  + (abs(s2wupout) + s2wdownout)\
                                                  + (abs(c2wupout) + c2wdownout);
*/
                                            coffhit[1]=middleinvexp;
                                            match->invexp=middleinvexp;
                                            match->massratio=massratio;
                                            match->upin=mult->upin;
                                            match->downin=mult->downin;
                                            match->piupin=mult->piupin;
                                            match->pidownin=mult->pidownin;
                                            match->aupin=mult->aupin;
                                            match->adownin=mult->adownin;
                                            match->e2upin=mult->e2upin;
                                            match->e2downin=mult->e2downin;
                                            match->nbvupin=mult->nbvupin;
                                            match->nbsupin=mult->nbsupin;
                                            match->upout=upout;
                                            match->downout=downout;
                                            match->piupout=piupout;
                                            match->pidownout=pidownout;
                                            match->aupout=aupout;
                                            match->adownout=adownout;
                                            match->e2upout=e2upout;
                                            match->e2downout=e2downout;
                                            match->s2wupout=s2wupout;
                                            match->s2wdownout=s2wdownout;
                                            match->c2wupout=c2wupout;
                                            match->c2wdownout=c2wdownout;
                                            match->matchcomplexity=complexity;
                                            match->match=multiplierout;
                                            initUses(&match->uses);
                                            addUses(&match->uses, &mult->uses);
                                            if (massratio == 0) {
                                              match->uses.G=1;
                                            } else if (massratio == 1) {
                                              match->uses.v=1;
                                            } else if (massratio == 2) {
                                              match->uses.mz=1;
                                            } else if (massratio == 3) {
                                              match->uses.mw=1;
                                            } else if (massratio == 4) {
                                              match->uses.mh0=1;
                                            }
                                            if (aupout != 0) {
                                              match->uses.alpha_em=1;
                                            }
                                            if ((s2wupout != 0) || (c2wupout != 0)) {
                                              match->uses.sin2w=1;
                                            }
                                            match->matchmult=e2out * mult->mult[abs(middleinvexp)];
                                            *nummatches=*nummatches+1;
                                            *matchptr=*matchptr+1;
                                            match=*matchptr;
                                          }  // if interesting middle
                                        } // sanity check c2w middle
                                      } // sanity check s2w middle

                                      // test multiplier against right coefficient
                                      if ((s2wdownout == 1) || (s2wdownout == 2) || (s2wdownout == rightinvexp) || (s2wdownout == (rightinvexp * 2))) {
                                        if ((c2wdownout == 1) || (c2wdownout == 2) || (c2wdownout == rightinvexp) || (c2wdownout == (rightinvexp * 2))) {
                                          multiplierout=cright * rightmassterm * s2wout * c2wout * e2out * mult->mult[abs(rightinvexp)];
                                          if (interesting(range, multiplierout)) { 
                                            complexity=mult->multcomplexity\
                                                  + (upout + downout)\
                                                  + (abs(piupout) + pidownout)\
                                                  + (abs(aupout) + adownout)\
                                                  + (abs(e2upout) + e2downout);
/*
                                                  + (abs(s2wupout) + s2wdownout)\
                                                  + (abs(c2wupout) + c2wdownout);
*/
                                            coffhit[2]=rightinvexp;
                                            match->invexp=rightinvexp;
                                            match->massratio=massratio;
                                            match->upin=mult->upin;
                                            match->downin=mult->downin;
                                            match->piupin=mult->piupin;
                                            match->pidownin=mult->pidownin;
                                            match->aupin=mult->aupin;
                                            match->adownin=mult->adownin;
                                            match->e2upin=mult->e2upin;
                                            match->e2downin=mult->e2downin;
                                            match->nbvupin=mult->nbvupin;
                                            match->nbsupin=mult->nbsupin;
                                            match->upout=upout;
                                            match->downout=downout;
                                            match->piupout=piupout;
                                            match->pidownout=pidownout;
                                            match->aupout=aupout;
                                            match->adownout=adownout;
                                            match->e2upout=e2upout;
                                            match->e2downout=e2downout;
                                            match->s2wupout=s2wupout;
                                            match->s2wdownout=s2wdownout;
                                            match->c2wupout=c2wupout;
                                            match->c2wdownout=c2wdownout;
                                            match->matchcomplexity=complexity;
                                            match->match=multiplierout;
                                            initUses(&match->uses);
                                            addUses(&match->uses, &mult->uses);
                                            if (massratio == 0) {
                                              match->uses.G=1;
                                            } else if (massratio == 1) {
                                              match->uses.v=1;
                                            } else if (massratio == 2) {
                                              match->uses.mz=1;
                                            } else if (massratio == 3) {
                                              match->uses.mw=1;
                                            } else if (massratio == 4) {
                                              match->uses.mh0=1;
                                            }
                                            if (aupout != 0) {
                                              match->uses.alpha_em=1;
                                            }
                                            if ((s2wupout != 0) || (c2wupout != 0)) {
                                              match->uses.sin2w=1;
                                            }
                                            match->matchmult=e2out * mult->mult[abs(rightinvexp)];
                                            *nummatches=*nummatches+1;
                                            *matchptr=*matchptr+1;
                                            match=*matchptr;
                                          }  // if interesting right
                                        } // sanity check c2w right
                                      } // sanity check s2w right
                                      mult++;
                                    }  // for i
#ifdef SIN2W
                                  } // gcd c2wout
                                } // c2wdownout
                              } // c2wupout
                            } // gcd s2wout
                          } // s2wdownout
                        } // s2wupout
#endif
                      } // e2upout
                    } // gcd aout
                  } // adownout
                } // aupout
/*
              } // gcd piout
            } // pidownout
*/
          } // piupout
        } // gcd updownout
      } // downout
    } // upout
  } // for massratio
#ifdef SHOWSTATUS
  printf(".\n");
#endif
}
