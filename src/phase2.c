#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "nle-lepton.h"
#include "reference.h"
#include "util.h"

double solveNLEforMasses(char *exponents, int leftinvexp, int middleinvexp, int rightinvexp, matches *leftmatchptr, matches *middlematchptr, matches *rightmatchptr, input_use *alluses, int maxcomplexity) {
  /***********/
  /* Phase 2 */
  /***********/
  // solve polynomial-like non-linear equation for particle masses using the supplied coefficients, exponents,
  // the electron mass and fine-structure constant as inputs.  Return the relative difference between
  // the computed and experimental muon mass (as an indication of accuracty)
  long int samples=0;
  int i;
  int arange, merange, murange, vrange, grange, mzrange, mwrange, mh0range, taurange, sin2wrange;
  struct timespec starttime, starttime2, endtime;
  double elapsedtime;
  int floata;
  int floatme;
  int floatmu;
  int floatv;
  int floattau;
  int floatg;
  int floatmz;
  int floatmw;
  int floatmh0;
  int floatsin2w;
  double score;
  char massstr[20];
  char massstrinv[20];
  char execstr[352];
  char outstr01[320];
  char outstr02[320];
  char outstr03[320];
  char outstr04[320];
  char outstr05[320];
  char outstr06[320];
  char outstr07[320];
  char outstr08[320];
  char outstr09[320];
  char outstr10[320];
  char outstr11[320];
  char outstr12[320];
  char outstr13[320];
  char outstr14[320];
  char outstr15[320];
  char outstr16[320];
  char outstr17[320];
  char outstr18[320];
  char usedasinput[5];
  char usedasoutput[5];
  long long resulthash;
  int complexity;
  double leftexp;
  double middleexp;
  double rightexp;
  matches *matchptr;
  char leftformulastr [288];
  char middleformulastr [288];
  char rightformulastr [288];
  int upout;
  int downout;
  char e2out[16];
  char updownout[16];
  char piout[16];
  char aout[16];
  char s2wout[16];
  char c2wout[16];
  int invexp;
  char updownin[16];
  char nbin[16];
  char e2in[16];
  char piin[16];
  char ain[16];
  int unknowns;
  int mwmzmode;
  int symmetry;
  float combinedscore;
 
  double alpha_ref_relerror=alpha_ref_error / alpha_ref;
  double me_ref_relerror=me_ref_error / me_ref;
  double mu_ref_relerror=mu_ref_error / mu_ref;  
  double tau_ref_relerror=tau_ref_error / tau_ref;
  double v_ref_relerror=v_ref_error / v_ref;
  double G_ref_relerror=G_ref_error / G_ref;
  double mz_ref_relerror=mz_ref_error / mz_ref;
  double mw_ref_relerror=mw_ref_error / mw_ref;
  double mh0_ref_relerror=mh0_ref_error / mh0_ref;
  double sin2w_ref_relerror=sin2w_ref_error / sin2w_ref;

  //  mc test vars
  double r;
  double e_test=0;
  double u_test=0;
  double t_test=0;
  double precision=0;
  double precision_last=0;
  double left, middle, right;
  double leftstatic, middlestatic, rightstatic;
  double lefts2w, middles2w, rights2w;
  double leftc2w, middlec2w, rightc2w;
  double leftmassterm, middlemassterm, rightmassterm;
  double leftmeterm, middlemeterm, rightmeterm;
  double leftmuterm, middlemuterm, rightmuterm;
  double leftmtterm, middlemtterm, rightmtterm;
  double mp;
  double worst_test;
  double rangefactor;
  double rangemultiplier=5.0;
  int stalled;
  int progress;

  // these tunings affect speed and reliability, adjust with extreme care
  double testratio=25.0;              // acceptable ratios of e_test/u_test/t_test, coefficient search ranges are guided by the least precise term so keeping test term ratios relatively close together optimizes search ranges for all coefficients
  int ratiograceperiod=15;            // ignore test ratio until this much progress has been achieved.   Ratios are typically way off at the beginning.   Search ranges need to be able to find solutions within the ratio limits before this trigger
  int stalledlimit=500000;            // most formulas can be solved with less than 500,000 samples, if not then it is probably hard to solve (like P+12+13+14, P+24+25+26, etc.)
  double defaultrangemultiplier=5.0;  // lowest practical range multiplier, fastest for most formulas
  double stalledrangemultiplier=10.0; // this value works better for slow to solve formulas and fast formulas that get stuck.  Will automatically revert to default if just temporarily stuck.  For slow to solve formulas this will continuously trigger
  int slowcheckpoint=1000000;         // progress point to check on slow processes
  double stuckprecision=1.0E-2;       // if precision is not past this level by slowcheckpoint, try resetting

  // mc outputs
  double alpha=0;
  double alpha_last=0;
  double alpha_center=0;
  double alpha_range=0;
  double alpha_range_new=0;
  double me=0;
  double me_last=0;
  double me_center=0;
  double me_range=0;
  double me_range_new=0;
  double mu=0;
  double mu_last=0;
  double mu_center=0;
  double mu_range=0;
  double mu_range_new=0;
  double v=0;
  double v_last=0;
  double v_center=0;
  double v_range=0;
  double v_range_new=0;
  double tau=0;
  double tau_last=0;
  double tau_center=0;
  double tau_range=0;
  double tau_range_new=0;
  double G=0;
  double G_last=0;
  double G_center=0;
  double G_range=0;
  double G_range_new=0;
  double mz=0;
  double mz_last=0;
  double mz_center=0;
  double mz_range=0;
  double mz_range_new=0;
  double mw=0;
  double mw_last=0;
  double mw_center=0;
  double mw_range=0;
  double mw_range_new=0;
  double mh0=0;
  double mh0_last=0;
  double mh0_center=0;
  double mh0_range=0;
  double mh0_range_new=0;
  double sin2w=0;
  double sin2w_last=0;
  double sin2w_center=0;
  double sin2w_range=0;
  double sin2w_range_new=0;
  double cos2w=0;

  // for reporting
  double alpha_out=0;
  double alpha_out_low=1.0E30;
  double alpha_out_high=-1.0E30;
  double alpha_out_c=0;
  double alpha_out_error=0;
  double alpha_out_relerror=0;
  double alpha_out_diff=0;
  double alpha_out_reldiff=0;
  double me_out=0;
  double me_out_low=1.0E30;
  double me_out_high=-1.0E30;
  double me_out_c=0;
  double me_out_error=0;
  double me_out_relerror=0;
  double me_out_diff=0;
  double me_out_reldiff=0;
  double mu_out=0;
  double mu_out_low=1.0E30;
  double mu_out_high=-1.0E30;
  double mu_out_c=0;
  double mu_out_error=0;
  double mu_out_relerror=0;
  double mu_out_diff=0;
  double mu_out_reldiff=0;
  double tau_out=0;
  double tau_out_low=1.0E30;
  double tau_out_high=-1.0E30;
  double tau_out_c=0;
  double tau_out_error=0;
  double tau_out_relerror=0;
  double tau_out_diff=0;
  double tau_out_reldiff=0;
  double v_out=0;
  double v_out_low=1.0E30;
  double v_out_high=-1.0E30;
  double v_out_c=0;
  double v_out_error=0;
  double v_out_relerror=0;
  double v_out_diff=0;
  double v_out_reldiff=0;
  double G_out=0;
  double G_out_low=1.0E30;
  double G_out_high=-1.0E30;
  double G_out_c=0;
  double G_out_error=0;
  double G_out_relerror=0;
  double G_out_diff=0;
  double G_out_reldiff=0;
  double mz_out=0;
  double mz_out_low=1.0E30;
  double mz_out_high=-1.0E30;
  double mz_out_c=0;
  double mz_out_error=0;
  double mz_out_relerror=0;
  double mz_out_diff=0;
  double mz_out_reldiff=0;
  double mw_out=0;
  double mw_out_low=1.0E30;
  double mw_out_high=-1.0E30;
  double mw_out_c=0;
  double mw_out_error=0;
  double mw_out_relerror=0;
  double mw_out_diff=0;
  double mw_out_reldiff=0;
  double mh0_out=0;
  double mh0_out_low=1.0E30;
  double mh0_out_high=-1.0E30;
  double mh0_out_c=0;
  double mh0_out_error=0;
  double mh0_out_relerror=0;
  double mh0_out_diff=0;
  double mh0_out_reldiff=0;
  double sin2w_out=0;
  double sin2w_out_low=1.0E30;
  double sin2w_out_high=-1.0E30;
  double sin2w_out_c=0;
  double sin2w_out_error=0;
  double sin2w_out_relerror=0;
  double sin2w_out_diff=0;
  double sin2w_out_reldiff=0;

  clock_gettime(CLOCK_REALTIME, &starttime);

  leftstatic=((double)leftmatchptr->matchup / (double)leftmatchptr->matchdown) / leftmatchptr->matchmult;
  middlestatic=((double)middlematchptr->matchup / (double)middlematchptr->matchdown) / middlematchptr->matchmult;
  rightstatic=((double)rightmatchptr->matchup / (double)rightmatchptr->matchdown) / rightmatchptr->matchmult;

  leftexp = 1.0 / (double)leftinvexp;
  middleexp = 1.0 / (double)middleinvexp;
  rightexp = 1.0 / (double)rightinvexp;

  // float the three used parameters with the highest experimental uncertainty, use the rest as inputs
  unknowns=0;
  // mH0 relative uncertainty 1.3E-3
  if (alluses->mh0 == 1) {
    floatmh0=1;
    unknowns++;
  } else {
    floatmh0=0;
    mh0_center=(double)mh0_ref;
    mh0_range=(double)mh0_ref_error;
  }

  // sin2w relative uncertainty 1.2E-3, sinw: 5.9E-4, cosw: 1.7E-4
  if (alluses->sin2w == 1) {
    if ((alluses->mw == 1) || (alluses->mz == 1)) {  // check if we are explicitly using mw or mz
      mwmzmode=1; // sin2w will be derived from mw and mz
      floatsin2w=0;
      sin2w_center=(double)sin2w_ref;
      sin2w_range=(double)sin2w_ref_error;
    } else {
      mwmzmode=0;
      floatsin2w=1;
      unknowns++;
    }
  } else {
    mwmzmode=0;
    floatsin2w=0;
    sin2w_center=(double)sin2w_ref;
    sin2w_range=(double)sin2w_ref_error;
  }

#ifdef WIDE_G // higher uncertainty puts G here when this is selected
  // G(wide) relative uncertainty 6.0E-4, mp(wide): 3.0E-4
  if (alluses->G == 1) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=(double)G_ref;
    G_range=(double)G_ref_error;
  }
#endif

  // after this point we must also check for three unknowns

  // mw relative uncertainty 1.5E-4
  if (((alluses->mw == 1) || (mwmzmode == 1)) && (unknowns < 3)) {
    floatmw=1;
    unknowns++;
  } else {
    floatmw=0;
    mw_center=(double)mw_ref;
    mw_range=(double)mw_ref_error;
  }

#ifdef ROSI_G // higher uncertainty puts G here when this is selected
  // G(rosi) relative uncertainty 1.5E-4, mp(rosi): 7.4E-5
  if ((alluses->G == 1) && (unknowns < 3)) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=(double)G_ref;
    G_range=(double)G_ref_error;
  }
#endif

  // tau relative uncertainty 6.8E-5
  if (unknowns < 3) { // tau is always used but only floated if we have less than 3 uknowns so far
    floattau=1;
    unknowns++;
  } else {
    floattau=0;
    tau_center=(double)mu_ref;
    tau_range=(double)mu_ref_error;
  }

#ifdef CODATA_G
  // G(codata) relative uncertainty 4.6E-5, mp(codata): 2.3E-5
  if ((alluses->G == 1) && (unknowns < 3)) {
    floatg=1;
    unknowns++;
  } else {
    floatg=0;
    G_center=(double)G_ref;
    G_range=(double)G_ref_error;
  }
#endif

  // mz relative uncertainty 2.3E-5
  if (((alluses->mz == 1) || (mwmzmode == 1)) && (unknowns < 3)) {
    floatmz=1;
    unknowns++;
  } else {
    floatmz=0;
    mz_center=(double)mz_ref;
    mz_range=(double)mz_ref_error;
  }

  // v relative uncertainty 2.6E-7
  if ((alluses->v == 1) && (unknowns < 3)) {
    floatv=1;
    unknowns++;
  } else {
    floatv=0;
    v_center=(double)v_ref;
    v_range=(double)v_ref_error;
  }

  // mu relative uncertainty 2.3E-8
  if (unknowns < 3) {  // mu is always used but only floated if we have less than 3 uknowns so far
    floatmu=1;
    unknowns++;
  } else {
    floatmu=0;
    mu_center=(double)mu_ref;
    mu_range=(double)mu_ref_error;
  }

  // me relative uncertainty 3.0E-10
  if (unknowns < 3) {  // me is always used but only floated if we have less than 3 uknowns so far
    floatme=1;
    unknowns++;
  } else {
    floatme=0;
    me_center=(double)me_ref;
    me_range=(double)me_ref_error;
  }

  // alpha relative uncertainty 1.5E-10
  if (unknowns < 3) {  
    floata=1;
    unknowns++;
  } else {
    floata=0;
    alpha_center=(double)alpha_ref;
    alpha_range=(double)alpha_ref_error;
  }

  // systematically try all non-floated input extremes
#ifdef DEBUG21
  printf("debug, Begin phase 2 input loops, exponents: %s, unknowns: %d, floatmu: %d, floatv: %d, floatmz: %d, floatg: %d, floattau: %d, floatmw: %d, floatsin2w: %d, mwmzmode: %d, floatmh0: %d\n", exponents, unknowns, floatmu, floatv, floatmz, floatg, floattau, floatmw, floatsin2w, mwmzmode, floatmh0);
  printUses(alluses);
  fflush(stdout);
#endif

#ifdef IGNORE_SMALL_UNCERTAINTIES
  // this speeds up phase 2 by up to 4x, don't use if me or a would be floated
  alpha=(double)alpha_ref;
  alpha_range=alpha_ref_error;
  arange=1;
  me=(double)me_ref;
  me_range=me_ref_error;
  merange=1;
#else
  for (arange=!(alluses->alpha_em || floata); arange <= 1; arange++) {
    if (arange == 0) {
      alpha=(alpha_center - alpha_range);
    } else {
      alpha=(alpha_center + alpha_range);
    }  
    for (merange=floatme; merange <= 1; merange++) {
      if (merange == 0) {
        me=(me_center - me_range);
      } else {
        me=(me_center + me_range);
      }
#endif
      for (murange=floatmu; murange <= 1; murange++) {
        if (floatmu == 0) {
          if (murange == 0) {
            mu=(mu_center - mu_range);
          } else {
            mu=(mu_center + mu_range);
          }
        }
        for (vrange=!(alluses->v || floatv); vrange <= 1; vrange++) {
          if (floatv == 0) {
            if (vrange == 0) {
              v=(v_center - v_range);
            } else {
              v=(v_center + v_range);
            }
          }
          for (mzrange=!((alluses->mz || mwmzmode) || floatmz); mzrange <= 1; mzrange++) {
            if (floatmz == 0) {
              if (mzrange == 0) {
                mz=(mz_center - mz_range);
              } else {
                mz=(mz_center + mz_range);
              }
            }
            for (grange=!(alluses->G || floatg); grange <= 1; grange++) {
              if (floatg == 0) { 
                if (grange == 0) {
                  G=(G_center - G_range);
                } else {
                  G=(G_center + G_range);
                } 
                mp=(double)kg_to_ev * (double)sqrt(hbar_ref * c_ref / G);
              } 
              for (taurange=floattau; taurange <= 1; taurange++) {
                // tau is always used but only floated if necessary
                if (floattau == 0) {
                  if (taurange == 0) {
                    tau=(tau_center - tau_range);
                  } else {
                    tau=(tau_center + tau_range);
                  }
                }
                for (mwrange=!((alluses->mw || mwmzmode) || floatmw); mwrange <= 1; mwrange++) {
                  if (floatmw == 0) {
                    if (mwrange == 0) {
                      mw=(mw_center - mw_range);
                    } else {
                      mw=(mw_center + mw_range);
                    }
                  }
                  for (sin2wrange=!(alluses->sin2w || floatsin2w); sin2wrange <= 1; sin2wrange++) {
                    if (floatsin2w == 0) {
                      if (sin2wrange == 0) {
                        sin2w=(sin2w_center - sin2w_range);
                      } else {
                        sin2w=(sin2w_center + sin2w_range);
                      }
                    }
                    for (mh0range=!(alluses->mh0 || floatmh0); mh0range <= 1; mh0range++) {
                      if (floatmh0 == 0) {
                        if (mh0range == 0) {
                          mh0=(mh0_center - mh0_range);
                        } else {
                          mh0=(mh0_center + mh0_range);
                        }
                      }
#ifdef DEBUG21
                      printf("debug, Begin    phase 2 samples loop, exponents: %s, arange: %d, merange: %d, murange: %d, vrange: %d, mzrange: %d, grange: %d, taurange: %d, mwrange: %d, sin2wrange: %d, mh0range: %d\n", exponents, arange, merange, murange, vrange, mzrange, grange, taurange, mwrange, sin2wrange, mh0range);
                      fflush(stdout);
#endif
                      clock_gettime(CLOCK_REALTIME, &starttime2);
                      precision_last=1.0E99;
                      //  reset mc test vars and outputs
                      if (floata == 1) {
                        alpha_last=(double)alpha_ref;
                        alpha_center=(double)alpha_ref;
                        alpha_range=(double)alpha_ref * 0.1;
                      }
                      if (floatme == 1) {
                        me_last=(double)me_ref;
                        me_center=(double)me_ref;
                        me_range=(double)me_ref * 0.1;
                      }
                      if (floatmu == 1) {
                        mu_last=(double)mu_ref;
                        mu_center=(double)mu_ref;
                        mu_range=(double)mu_ref * 0.1;
                      }
                      if (floatv == 1) {
                        v_last=(double)v_ref;
                        v_center=(double)v_ref; 
                        v_range=(double)v_ref * 0.1;
                      }
                      if (floatmz == 1) {
                        mz_last=(double)mz_ref;
                        mz_center=(double)mz_ref;
                        mz_range=(double)mz_ref * 0.1;
                      }
                      if (floatg == 1) {
                        G_last=(double)G_ref;
                        G_center=(double)G_ref; 
                        G_range=(double)G_ref * 0.1;
                      }
                      if (floattau == 1) {
                        tau_last=(double)tau_ref;
                        tau_center=(double)tau_ref; 
                        tau_range=(double)tau_ref * 0.1;
                      }
                      if (floatmw == 1) {
                        mw_last=(double)mw_ref;
                        mw_center=(double)mw_ref; 
                        mw_range=(double)mw_ref * 0.1;
                      }
                      if (floatsin2w == 1) {
                        sin2w_last=(double)sin2w_ref;
                        sin2w_center=(double)sin2w_ref;
                        sin2w_range=(double)sin2w_ref * 0.1;
                      }
                      if (floatmh0 == 1) {
                        mh0_last=(double)mh0_ref;
                        mh0_center=(double)mh0_ref;
                        mh0_range=(double)mh0_ref * 0.1;
                      }
                      precision_last=1E99;
                      stalled=0;
                      progress=0;
                      rangemultiplier=defaultrangemultiplier;
                      for (samples=0; precision_last > 1.0E-11; samples++) {
                        if ((samples > 1) && ((samples % slowcheckpoint) == 0)) { // check on slow processes
#ifdef DEBUG20
                          clock_gettime(CLOCK_REALTIME, &endtime);
                          elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                          if ((samples % 10000000) == 0) { // rate limit periodic debug prints
                            printf ("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last:  %.3e, precision: %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e, mz: %.9e, mz_range: %.4e, mw: %.9e, mw_range: %.4e, sin2w: %.9e, sin2w_range: %.4e, mh0: %.9e, mh0_range: %.4e\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last, precision, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range);
                            fflush(stdout);
                          }
#endif

                          if ((progress == ratiograceperiod) || (precision_last > stuckprecision)) { // it's stuck, try resetting
#ifdef DEBUG20
                            clock_gettime(CLOCK_REALTIME, &endtime);
                            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                            printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, resetting\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                            fflush(stdout);
#endif
                            //  reset mc test vars and outputs
                            if (floata == 1) {
                              alpha_last=(double)alpha_ref;
                              alpha_center=(double)alpha_ref;
                              alpha_range=(double)alpha_ref * 0.1;
                            }
                            if (floatme == 1) {
                              me_last=(double)me_ref;
                              me_center=(double)me_ref;
                              me_range=(double)me_ref * 0.1;
                            }
                            if (floatmu == 1) {
                              mu_last=(double)mu_ref;
                              mu_center=(double)mu_ref;
                              mu_range=(double)mu_ref * 0.1;
                            }
                            if (floatv == 1) {
                              v_last=(double)v_ref;
                              v_center=(double)v_ref;
                              v_range=(double)v_ref * 0.1;
                            }
                            if (floatmz == 1) {
                              mz_last=(double)mz_ref;
                              mz_center=(double)mz_ref;
                              mz_range=(double)mz_ref * 0.1;
                            }
                            if (floatg == 1) {
                              G_last=(double)G_ref;
                              G_center=(double)G_ref;
                              G_range=(double)G_ref * 0.1;
                            }
                            if (floattau == 1) {
                              tau_last=(double)tau_ref;
                              tau_center=(double)tau_ref;
                              tau_range=(double)tau_ref * 0.1;
                            }
                            if (floatmw == 1) {
                              mw_last=(double)mw_ref;
                              mw_center=(double)mw_ref;
                              mw_range=(double)mw_ref * 0.1;
                            }
                            if (floatsin2w == 1) {
                              sin2w_last=(double)sin2w_ref;
                              sin2w_center=(double)sin2w_ref;
                              sin2w_range=(double)sin2w_ref * 0.1;
                            }
                            if (floatmh0 == 1) {
                              mh0_last=(double)mh0_ref;
                              mh0_center=(double)mh0_ref;
                              mh0_range=(double)mh0_ref * 0.1;
                            }
                            precision_last=1E99;
                            stalled=0;
                            progress=0;
                            rangemultiplier=defaultrangemultiplier;
                          }
                        } // end check on slow processes
                        if (stalled == stalledlimit) {
                          rangemultiplier=stalledrangemultiplier; // may be a slow solution, try bigger multiplier
                          rangefactor=worst_test * rangemultiplier;
#ifdef DEBUG20
                          clock_gettime(CLOCK_REALTIME, &endtime);
                          elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                          printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, stalled\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
#endif
                          if (floata == 1) {
                            alpha_range=alpha_last * rangefactor;
                          }
                          if (floatme == 1) {
                            me_range=me_last * rangefactor;
                          }
                          if (floattau == 1) {
                            tau_range=tau_last * rangefactor;
                          }
                          if (floatmu == 1) {
                            mu_range=mu_last * rangefactor;
                          }
                          if (floatv == 1) {
                            v_range=v_last * rangefactor;
                          }
                          if (floatg == 1) {
                            G_range=G_last * rangefactor;
                          }
                          if (floatmz == 1) {
                            mz_range=mz_last * rangefactor;
                          }
                          if (floatmw == 1) {
                            mw_range=mw_last * rangefactor;
                          }
                          if (floatmh0 == 1) {
                            mh0_range=mh0_last * rangefactor;
                          }
                          if (floatsin2w == 1) {
                            sin2w_range=sin2w_last * rangefactor;
                          }
                        }
                        stalled++;
                        // guess random values for mc outputs
                        if (floata == 1) {
                          r=(double)drand48();
                          alpha=((alpha_center - alpha_range) + (r * 2.0 * alpha_range));
                          i=0;
                          while ((alpha < alpha_ref * 0.9) || (alpha > (alpha_ref * 1.1))) { // sanity check on alpha to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, alpha range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              alpha_last=0;
                              alpha_center=(double)alpha_ref;
                              alpha_range=(double)alpha_ref * 0.1;
                            }
                            r=(double)drand48();
                            alpha=((alpha_center - alpha_range) + (r * 2.0 * alpha_range));
                            i++;
                          }
                        }

                        if (floatme == 1) {
                          r=(double)drand48();
                          me=((me_center - me_range) + (r * 2.0 * me_range));
                          i=0;
                          while ((me < me_ref * 0.9) || (me > (me_ref * 1.1))) { // sanity check on me mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, me range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              me_last=0;
                              me_center=(double)me_ref;
                              me_range=(double)me_ref * 0.1;
                            }
                            r=(double)drand48();
                            me=((me_center - me_range) + (r * 2.0 * me_range));
                            i++;
                          }
                        }
                        if (floatmu == 1) {
                          r=(double)drand48();
                          mu=((mu_center - mu_range) + (r * 2.0 * mu_range));
                          i=0;
                          while ((mu < mu_ref * 0.9) || (mu > (mu_ref * 1.1))) { // sanity check on mu mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mu range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mu_last=0;
                              mu_center=(double)mu_ref;
                              mu_range=(double)mu_ref * 0.1;
                            }
                            r=(double)drand48();
                            mu=((mu_center - mu_range) + (r * 2.0 * mu_range));
                            i++;
                          }
                        }
                        if (floatv == 1) {
                          r=(double)drand48();
                          v=((v_center - v_range) + (r * 2.0 * v_range));
                          i=0;
                          while ((v < (v_ref * 0.9)) || (v > (v_ref * 1.1))) { // sanity check to help convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, v range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              v_last=0;
                              v_center=(double)v_ref;
                              v_range=(double)v_ref * 0.1;
                            }
                            r=(double)drand48();
                            v=((v_center - v_range) + (r * 2.0 * v_range));
                            i++;
                          }
                        }
                        if (floatmz == 1) {
                          r=(double)drand48();
                          mz=((mz_center - mz_range) + (r * 2.0 * mz_range));
                          i=0;
                          while ((mz < mz_ref * 0.9) || (mz > (mz_ref * 1.1))) { // sanity check on mz mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mz range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mz_last=0;
                              mz_center=(double)mz_ref;
                              mz_range=(double)mz_ref * 0.1;
                            }
                            r=(double)drand48();
                            mz=((mz_center - mz_range) + (r * 2.0 * mz_range));
                            i++;
                          }
                        }
                        if (floatg == 1) {
                          r=(double)drand48();
                          G=((G_center - G_range) + (r * 2.0 * G_range));
                          i=0;
                          while ((G < (G_ref * 0.9)) || (G > (G_ref * 1.1))) {  // sanity check to help convergance
                            if (i > 50) {  // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, G range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              G_last=0;
                              G_center=(double)G_ref;
                              G_range=(double)G_ref * 0.1;
                            }
                            r=(double)drand48();
                            G=((G_center - G_range) + (r * 2.0 * G_range));
                            i++;
                          }
                          mp=(double)kg_to_ev * (double)sqrt(hbar_ref * c_ref / G);
                        }
                        if (floattau == 1) {
                          r=(double)drand48();
                          tau=((tau_center - tau_range) + (r * 2.0 * tau_range));
                          i=0;
                          while ((tau < tau_ref * 0.9) || (tau > (tau_ref * 1.1))) { // sanity check on tau mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, tau range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              tau_last=0;
                              tau_center=(double)tau_ref;
                              tau_range=(double)tau_ref * 0.1;
                            }
                            r=(double)drand48();
                            tau=((tau_center - tau_range) + (r * 2.0 * tau_range));
                            i++;
                          }
                        }
                        if (floatmw == 1) {
                          r=(double)drand48();
                          mw=((mw_center - mw_range) + (r * 2.0 * mw_range));
                          i=0;
                          while ((mw < mw_ref * 0.9) || (mw > (mw_ref * 1.1)) || (mw >= mz)) { // sanity check on mw mass to help solve correct root and speed convergance and prevent mw >= mz
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mw range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mw_last=0;
                              mw_center=(double)mw_ref;
                              mw_range=(double)mw_ref * 0.1;
                            }
                            r=(double)drand48();
                            mw=((mw_center - mw_range) + (r * 2.0 * mw_range));
                            i++;
                          }
                        }
                        if (alluses->sin2w == 1) {  
                          if (mwmzmode == 0) {  // float sin2w
                            r=(double)drand48();
                            sin2w=((sin2w_center - sin2w_range) + (r * 2.0 * sin2w_range));
                            i=0;
                            while ((sin2w < sin2w_ref * 0.9) || (sin2w > (sin2w_ref * 1.1))) { // sanity check on sin2w to help solve correct root and speed convergance
                              if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                                clock_gettime(CLOCK_REALTIME, &endtime);
                                elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                                printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, sin2w range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                                fflush(stdout);
#endif
                                i=0;
                                sin2w_last=0;
                                sin2w_center=(double)sin2w_ref; 
                                sin2w_range=(double)sin2w_ref * 0.1;
                              }
                              r=(double)drand48();
                              sin2w=((sin2w_center - sin2w_range) + (r * 2.0 * sin2w_range));
                              i++;
                            }
                            cos2w=1.0 - sin2w;
                          } else { // derive sin2w from mw/mz
                            cos2w=pow((mw/mz), 2.0);
                            sin2w=1.0 - cos2w;
                          } // end mwmzmode
                        } // end alluses sin2w
                        if (floatmh0 == 1) {
                          r=(double)drand48();
                          mh0=((mh0_center - mh0_range) + (r * 2.0 * mh0_range));
                          i=0;
                          while ((mh0 < mh0_ref * 0.9) || (mh0 > (mh0_ref * 1.1))) { // sanity check on mh0 mass to help solve correct root and speed convergance
                            if (i > 50) { // safety valve in case search gets out of bounds
#ifdef DEBUG20
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last: %.3e, mh0 range error\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last);
                              fflush(stdout);
#endif
                              i=0;
                              mh0_last=0;
                              mh0_center=(double)mh0_ref;
                              mh0_range=(double)mh0_ref * 0.1;
                            }
                            r=(double)drand48();
                            mh0=((mh0_center - mh0_range) + (r * 2.0 * mh0_range));
                            i++;
                          }
                        }

                        if (leftmatchptr->massratio == 0) {
                          leftmassterm=mp;
                        } else if (leftmatchptr->massratio == 1) {
                          leftmassterm=v;
                        } else if (leftmatchptr->massratio == 2) {
                          leftmassterm=mz;
                        } else if (leftmatchptr->massratio == 3) {
                          leftmassterm=mw;
                        } else if (leftmatchptr->massratio == 4) {
                          leftmassterm=mh0;
                        }
                        if (middlematchptr->massratio == 0) {
                          middlemassterm=mp;
                        } else if (middlematchptr->massratio == 1) {
                          middlemassterm=v;
                        } else if (middlematchptr->massratio == 2) {
                          middlemassterm=mz;
                        } else if (middlematchptr->massratio == 3) {
                          middlemassterm=mw;
                        } else if (middlematchptr->massratio == 4) {
                          middlemassterm=mh0;
                        }
                        if (rightmatchptr->massratio == 0) {
                          rightmassterm=mp;
                        } else if (rightmatchptr->massratio == 1) {
                          rightmassterm=v;
                        } else if (rightmatchptr->massratio == 2) {
                          rightmassterm=mz;
                        } else if (rightmatchptr->massratio == 3) {
                          rightmassterm=mw;
                        } else if (rightmatchptr->massratio == 4) {
                          rightmassterm=mh0;
                        }

                        leftmeterm=me / leftmassterm;
                        leftmuterm=mu / leftmassterm;
                        leftmtterm=tau / leftmassterm;
                        middlemeterm=me / middlemassterm;
                        middlemuterm=mu / middlemassterm;
                        middlemtterm=tau / middlemassterm;
                        rightmeterm=me / rightmassterm;
                        rightmuterm=mu / rightmassterm;
                        rightmtterm=tau / rightmassterm;

                        if (leftmatchptr->s2wupout != 0) {
                          lefts2w=pow(sin2w, ((double)leftmatchptr->s2wupout / (double)leftmatchptr->s2wdownout));
                        } else {
                          lefts2w=1.0;
                        }
                        if (leftmatchptr->c2wupout != 0) {
                          leftc2w=pow(cos2w, ((double)leftmatchptr->c2wupout / (double)leftmatchptr->c2wdownout));
                        } else {
                          leftc2w=1.0;
                        }

                        if (middlematchptr->s2wupout != 0) {
                          middles2w=pow(sin2w, ((double)middlematchptr->s2wupout / (double)middlematchptr->s2wdownout));
                        } else {
                          middles2w=1.0;
                        }
                        if (middlematchptr->c2wupout != 0) {
                          middlec2w=pow(cos2w, ((double)middlematchptr->c2wupout / (double)middlematchptr->c2wdownout));
                        } else {
                          middlec2w=1.0;
                        }

                        if (rightmatchptr->s2wupout != 0) {
                          rights2w=pow(sin2w, ((double)rightmatchptr->s2wupout / (double)rightmatchptr->s2wdownout));
                        } else {
                          rights2w=1.0;
                        }
                        if (rightmatchptr->c2wupout != 0) {
                          rightc2w=pow(cos2w, ((double)rightmatchptr->c2wupout / (double)rightmatchptr->c2wdownout));
                        } else {
                          rightc2w=1.0;
                        }

                        left=  (leftstatic   / lefts2w  ) / leftc2w;
                        middle=(middlestatic / middles2w) / middlec2w;
                        right= (rightstatic  / rights2w ) / rightc2w;

                        e_test=(left * pow(leftmeterm, leftexp)) - (middle * pow(middlemeterm, middleexp)) + (right * pow(rightmeterm, rightexp)) - 1.0;
                        u_test=(left * pow(leftmuterm, leftexp)) - (middle * pow(middlemuterm, middleexp)) + (right * pow(rightmuterm, rightexp)) - 1.0;
                        t_test=(left * pow(leftmtterm, leftexp)) - (middle * pow(middlemtterm, middleexp)) + (right * pow(rightmtterm, rightexp)) - 1.0;

#ifdef DEBUG23
                        printf("debug, exponents: %s, samples: %ld, left:   %.6e, leftstatic:   %.6e, lefts2w:   %.6e, leftc2w:   %.6e, lefts2wupout:   %d, lefts2wdownout:   %d, leftc2wupout:   %d, leftc2wdownout:   %d\n", exponents, samples, left, leftstatic, lefts2w, leftc2w, leftmatchptr->s2wupout, leftmatchptr->s2wdownout, leftmatchptr->c2wupout, leftmatchptr->c2wdownout);
                        printf("debug, exponents: %s, samples: %ld, middle: %.6e, middlestatic: %.6e, middles2w: %.6e, middlec2w: %.6e, middles2wupout: %d, middles2wdownout: %d, middlec2wupout: %d, middlec2wdownout: %d\n", exponents, samples, middle, middlestatic, middles2w, middlec2w, middlematchptr->s2wupout, middlematchptr->s2wdownout, middlematchptr->c2wupout, middlematchptr->c2wdownout);
                        printf("debug, exponents: %s, samples: %ld, right:  %.6e, rightstatic:  %.6e, rights2w:  %.6e, rightc2w:  %.6e, rights2wupout:  %d, rights2wdownout:  %d, rightc2wupout:  %d, rightc2wdownout:  %d\n", exponents, samples, right, rightstatic, rights2w, rightc2w, rightmatchptr->s2wupout, rightmatchptr->s2wdownout, rightmatchptr->c2wupout, rightmatchptr->c2wdownout);
                        printf("debug, exponents: %s, samples: %ld, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, left: %.3e, middle: %.3e, right: %.3e\n", exponents, samples, e_test, u_test, t_test, left, middle, right);
                        fflush(stdout);
#endif
                        if ((progress < ratiograceperiod) || (((fabs(e_test) / fabs(u_test)) < testratio) && ((fabs(u_test) / fabs(e_test)) < testratio))) {
                          if ((progress < ratiograceperiod) || (((fabs(e_test) / fabs(t_test)) < testratio) && ((fabs(t_test) / fabs(e_test)) < testratio) &&\
                             ((fabs(u_test) / fabs(t_test)) < testratio) && ((fabs(t_test) / fabs(u_test)) < testratio))) {
                            precision=fabs(e_test) + fabs(u_test) + fabs(t_test);
#ifdef DEBUG22
                            clock_gettime(CLOCK_REALTIME, &endtime);
                            elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                            printf ("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last:  %.3e, precision: %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e, mz: %.9e, mz_range: %.4e, mw: %.9e, mw_range: %.4e, sin2w: %.9e, sin2w_range: %.4e, mh0: %.9e, mh0_range: %.4e\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last, precision, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range);
                            fflush(stdout);
#endif
                            if (precision < precision_last) {
                              progress++;
                              stalled=0;
                              precision_last=precision;
                              alpha_last=alpha;
                              me_last=me;
                              mu_last=mu;
                              tau_last=tau;
                              v_last=v;
                              G_last=G;
                              mz_last=mz;
                              mw_last=mw;
                              mh0_last=mh0;
                              sin2w_last=sin2w;
                              if (fabs(e_test) > fabs(u_test)) {
                                worst_test=fabs(e_test);
                              } else {
                                worst_test=fabs(u_test);
                              }
                              if (fabs(t_test) > worst_test) {
                                worst_test=fabs(t_test);
                              }
                              rangemultiplier=defaultrangemultiplier;
                              rangefactor=worst_test * rangemultiplier;
                              if (rangefactor > 0.1) {
                                //  use default ranges
                                if (floata == 1) {
                                  alpha_center=alpha_last;
                                  alpha_range=(double)alpha_ref * 0.1;
                                }
                                if (floatme == 1) {
                                  me_center=me_last;
                                  me_range=(double)me_ref * 0.1;
                                }
                                if (floatmu == 1) {
                                  mu_center=mu_last;
                                  mu_range=(double)mu_ref * 0.1;
                                }
                                if (floatv == 1) {
                                  v_center=v_last;
                                  v_range=(double)v_ref * 0.1;
                                }
                                if (floatmz == 1) {
                                  mz_center=mz_last;
                                  mz_range=(double)mz_ref * 0.1;
                                }
                                if (floatg == 1) {
                                  G_center=G_last;
                                  G_range=(double)G_ref * 0.1;
                                }
                                if (floattau == 1) {
                                  tau_center=tau_last;
                                  tau_range=(double)tau_ref * 0.1;
                                }
                                if (floatmw == 1) {
                                  mw_center=mw_last;
                                  mw_range=(double)mw_ref * 0.1;
                                }
                                if (floatsin2w == 1) {
                                  sin2w_center=sin2w_last;
                                  sin2w_range=(double)sin2w_ref * 0.1;
                                }
                                if (floatmh0 == 1) {
                                  mh0_center=mh0_last;
                                  mh0_range=(double)mh0_ref * 0.1;
                                }
                              } else { 
                                if (floata == 1) {
                                  alpha_center=alpha_last;
                                  alpha_range_new=alpha_last * rangefactor;
                                  alpha_range=((alpha_range + alpha_range_new + alpha_range_new) / 3.0);
                                }
                                if (floatme == 1) {
                                  me_center=me_last;
                                  me_range_new=me_last * rangefactor;
                                  me_range=((me_range + me_range_new + me_range_new) / 3.0);
                                }
                                if (floattau == 1) {
                                  tau_center=tau_last;
                                  tau_range_new=tau_last * rangefactor;
                                  tau_range=((tau_range + tau_range_new + tau_range_new) / 3.0);
                                }
                                if (floatmu == 1) {
                                  mu_center=mu_last;
                                  mu_range_new=mu_last * rangefactor;
                                  mu_range=((mu_range + mu_range_new + mu_range_new) / 3.0);
                                }
                                if (floatv == 1) {
                                  v_center=v_last;
                                  v_range_new=v_last * rangefactor;
                                  v_range=((v_range + v_range_new + v_range_new) / 3.0);
                                }
                                if (floatg == 1) {
                                  G_center=G_last;
                                  G_range_new=G_last * rangefactor;
                                  G_range=((G_range + G_range_new + G_range_new) / 3.0);
                                }
                                if (floatmz == 1) {
                                  mz_center=mz_last;
                                  mz_range_new=mz_last * rangefactor;
                                  mz_range=((mz_range + mz_range_new + mz_range_new) / 3.0);
                                }
                                if (floatmw == 1) {
                                  mw_center=mw_last;
                                  mw_range_new=mw_last * rangefactor;
                                  mw_range=((mw_range + mw_range_new + mw_range_new) / 3.0);
                                }
                                if (floatmh0 == 1) {
                                  mh0_center=mh0_last;
                                  mh0_range_new=mh0_last * rangefactor;
                                  mh0_range=((mh0_range + mh0_range_new + mh0_range_new) / 3.0);
                                }
                                if (floatsin2w == 1) {
                                  sin2w_center=sin2w_last;
                                  sin2w_range_new=sin2w_last * rangefactor;
                                  sin2w_range=((sin2w_range + sin2w_range_new + sin2w_range_new) / 3.0);
                                }
                              }
#ifdef DEBUG21
                              clock_gettime(CLOCK_REALTIME, &endtime);
                              elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                              printf ("debug, exponents: %s, samples: %ld, time: %6.4fs, progress: %d, rangefactor: %.9e, precision_last:  %.3e, precision: %.3e, e_test:  %.3e, u_test:  %.3e, t_test: %.3e, tau: %.9e, tau_range: %.4e, G: %.9e, G_range: %.4e, v: %.9e, v_range: %.4e, mu: %.9e, mu_range: %.4e, mz: %.9e, mz_range: %.4e, mw: %.9e, mw_range: %.4e, sin2w: %.9e, sin2w_range: %.4e, mh0: %.9e, mh0_range: %.4e\n", exponents, samples, elapsedtime, progress, rangefactor, precision_last, precision, e_test, u_test, t_test, tau, tau_range, G, G_range, v, v_range, mu, mu_range, mz, mz_range, mw, mw_range, sin2w, sin2w_range, mh0, mh0_range);
                              fflush(stdout);
#endif
                            } // end if  precision < precision_last
                          } // end if t_test/u_test/e_test
                        } // end if e_test / u_test
                      } // end for samples 

#ifdef DEBUG21
                      clock_gettime(CLOCK_REALTIME, &endtime);
                      elapsedtime=((double)(endtime.tv_sec - 1500000000) + ((double)endtime.tv_nsec / 1.0E9)) - ((double)(starttime2.tv_sec - 1500000000) + ((double)starttime2.tv_nsec) / 1.0E9);
                      printf("debug, Finished phase 2 samples loop, exponents: %s, samples: %ld, mass mode: %d%d%d, precision: %.6e (%6.4fs)\n", exponents, samples, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, precision_last, elapsedtime);
                      fflush(stdout);
#endif
                      alpha_out=alpha_last;
                      if (alpha_out < alpha_out_low) {
                        alpha_out_low=alpha_out;
                      }
                      if (alpha_out > alpha_out_high) {
                        alpha_out_high=alpha_out;
                      }
                      alpha_out_c=((alpha_out_high + alpha_out_low) / 2.0);
                      alpha_out_error=(alpha_out_high - alpha_out_c);
                      alpha_out_relerror=alpha_out_error / alpha_out_c;
                      alpha_out_diff=alpha_out_c - alpha_ref;
                      alpha_out_reldiff=alpha_out_diff / alpha_ref;

                      me_out=me_last;
                      if (me_out < me_out_low) {
                        me_out_low=me_out;
                      }
                      if (me_out > me_out_high) {
                        me_out_high=me_out;
                      }
                      me_out_c=((me_out_high + me_out_low) / 2.0);
                      me_out_error=(me_out_high - me_out_c);
                      me_out_relerror=me_out_error / me_out_c;
                      me_out_diff=me_out_c - me_ref;
                      me_out_reldiff=me_out_diff / me_ref;

                      mu_out=mu_last;
                      if (mu_out < mu_out_low) {
                        mu_out_low=mu_out;
                      }
                      if (mu_out > mu_out_high) {
                        mu_out_high=mu_out;
                      }
                      mu_out_c=((mu_out_high + mu_out_low) / 2.0);
                      mu_out_error=(mu_out_high - mu_out_c);
                      mu_out_relerror=mu_out_error / mu_out_c;
                      mu_out_diff=mu_out_c - mu_ref;
                      mu_out_reldiff=mu_out_diff / mu_ref;

                      tau_out=tau_last;
                      if (tau_out < tau_out_low) {
                        tau_out_low=tau_out;
                      }
                      if (tau_out > tau_out_high) {
                        tau_out_high=tau_out;
                      }
                      tau_out_c=((tau_out_high + tau_out_low) / 2.0);
                      tau_out_error=(tau_out_high - tau_out_c);
                      tau_out_relerror=tau_out_error / tau_out_c;
                      tau_out_diff=tau_out_c - tau_ref;
                      tau_out_reldiff=tau_out_diff / tau_ref;

                      v_out=v_last;
                      if (v_out < v_out_low) {
                        v_out_low=v_out;
                      }
                      if (v_out > v_out_high) {
                        v_out_high=v_out;
                      }
                      v_out_c=((v_out_high + v_out_low) / 2.0);
                      v_out_error=(v_out_high - v_out_c);
                      v_out_relerror=v_out_error / v_out_c;
                      v_out_diff=v_out_c - v_ref;
                      v_out_reldiff=v_out_diff / v_ref;

                      G_out=G_last;
                      if (G_out < G_out_low) {
                        G_out_low=G_out;
                      }
                      if (G_out > G_out_high) {
                        G_out_high=G_out;
                      }
                      G_out_c=((G_out_high + G_out_low) / 2.0);
                      G_out_error=(G_out_high - G_out_c);
                      G_out_relerror=G_out_error / G_out_c;
                      G_out_diff=G_out_c - G_ref;
                      G_out_reldiff=G_out_diff / G_ref;

                      mz_out=mz_last;
                      if (mz_out < mz_out_low) {
                        mz_out_low=mz_out;
                      }
                      if (mz_out > mz_out_high) {
                        mz_out_high=mz_out;
                      }
                      mz_out_c=((mz_out_high + mz_out_low) / 2.0);
                      mz_out_error=(mz_out_high - mz_out_c);
                      mz_out_relerror=mz_out_error / mz_out_c;
                      mz_out_diff=mz_out_c - mz_ref;
                      mz_out_reldiff=mz_out_diff / mz_ref;

                      mw_out=mw_last;
                      if (mw_out < mw_out_low) {
                        mw_out_low=mw_out;
                      }
                      if (mw_out > mw_out_high) {
                        mw_out_high=mw_out;
                      }
                      mw_out_c=((mw_out_high + mw_out_low) / 2.0);
                      mw_out_error=(mw_out_high - mw_out_c);
                      mw_out_relerror=mw_out_error / mw_out_c;
                      mw_out_diff=mw_out_c - mw_ref;
                      mw_out_reldiff=mw_out_diff / mw_ref;

                      mh0_out=mh0_last;
                      if (mh0_out < mh0_out_low) {
                        mh0_out_low=mh0_out;
                      }
                      if (mh0_out > mh0_out_high) {
                        mh0_out_high=mh0_out;
                      }
                      mh0_out_c=((mh0_out_high + mh0_out_low) / 2.0);
                      mh0_out_error=(mh0_out_high - mh0_out_c);
                      mh0_out_relerror=mh0_out_error / mh0_out_c;
                      mh0_out_diff=mh0_out_c - mh0_ref;
                      mh0_out_reldiff=mh0_out_diff / mh0_ref;

                      sin2w_out=sin2w_last;
                      if (sin2w_out < sin2w_out_low) {
                        sin2w_out_low=sin2w_out;
                      }
                      if (sin2w_out > sin2w_out_high) {
                        sin2w_out_high=sin2w_out;
                      }
                      sin2w_out_c=((sin2w_out_high + sin2w_out_low) / 2.0);
                      sin2w_out_error=(sin2w_out_high - sin2w_out_c);
                      sin2w_out_relerror=sin2w_out_error / sin2w_out_c;
                      sin2w_out_diff=sin2w_out_c - sin2w_ref;
                      sin2w_out_reldiff=sin2w_out_diff / sin2w_ref;
                    } // end mh0range
                  } // end sin2wrange
                } // end mwrange
              } // end taurange
            } // end grange
          } // end mzrange
        } // end vrange
      }  // end murange
#ifndef IGNORE_SMALL_UNCERTAINTIES
    } // end merange
  }  // end arange
#endif

  score=((fmax((fabs(mu_out_reldiff) / mu_ref_relerror), 1.0) + fmax((fabs(tau_out_reldiff) / tau_ref_relerror), 1.0) + (alluses->G * fmax((fabs(G_out_reldiff) / G_ref_relerror), 1.0)) + (alluses->v * fmax((fabs(v_out_reldiff) / v_ref_relerror), 1.0)) + ((alluses->mz || mwmzmode) * fmax((fabs(mz_out_reldiff) / mz_ref_relerror), 1.0))\
       + ((alluses->mw || mwmzmode) * fmax((fabs(mw_out_reldiff) / mw_ref_relerror), 1.0)) + (alluses->mh0 * fmax((fabs(mh0_out_reldiff) / mh0_ref_relerror), 1.0)) + (alluses->sin2w * fmax((fabs(sin2w_out_reldiff) / sin2w_ref_relerror), 1.0))) / (2.0 + (double)alluses->G + (double)alluses->v + (double)(alluses->mz || mwmzmode) + (double)(alluses->mw || mwmzmode) + (double)alluses->mh0 + (double)alluses->sin2w)) - 1.0; 
#ifndef ALWAYS_SHOW_RESULTS
  if (score == 0.0) {
#endif
    complexity=leftmatchptr->matchcomplexity + middlematchptr->matchcomplexity + rightmatchptr->matchcomplexity;
    //resulthash=lrand48();
    resulthash=(leftmatchptr->matchhash ^ middlematchptr->matchhash) ^ rightmatchptr->matchhash;

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

    combinedscore = (float)complexity / (float)symmetry;

/* create formula output */
    matchptr=leftmatchptr;
    invexp=leftinvexp;
    // note: all terms except matchup and matchdown are inverted here as they represent offsets to the real coefficient
    upout=matchptr->matchup * matchptr->downout;
    downout=matchptr->matchdown * matchptr->upout;
    if ((upout == 1) && (downout == 1)) {
      sprintf(updownout, "       ");
    } else {
      sprintf(updownout, "(%2d/%2d)", upout, downout);
    }
    if (matchptr->e2upout == 0) {
      sprintf(e2out, "        ");
    } else {
      sprintf(e2out, "2^(%2d/%d)", -(matchptr->e2upout), matchptr->e2downout);
    }
    if (matchptr->piupout == 0) {
      sprintf(piout, "         ");
    } else {
      if ((matchptr->piupout == -1) && (matchptr->pidownout == 1)) {
        sprintf(piout, "pi       ");
      } else {
        sprintf(piout, "pi^(%2d/%d)", -(matchptr->piupout), matchptr->pidownout);
      }
    }
    if (matchptr->aupout == 0) {
      sprintf(aout, "        ");
    } else {
      if ((matchptr->aupout == -1) && (matchptr->adownout == 1)) {
        sprintf(aout, "a       ");
      } else {
        sprintf(aout, "a^(%2d/%d)", -(matchptr->aupout), matchptr->adownout);
      }
    }
    if (matchptr->s2wupout == 0) {
      sprintf(s2wout, "            ");
    } else {
      if ((matchptr->s2wupout == -1) && (matchptr->s2wdownout == 1)) {
        sprintf(s2wout, "sin2w       ");
      } else {
        sprintf(s2wout, "sin2w^(%2d/%d)", -(matchptr->s2wupout), matchptr->s2wdownout);
      }
    }
    if (matchptr->c2wupout == 0) { 
      sprintf(c2wout, "            ");
    } else {
      if ((matchptr->c2wupout == -1) && (matchptr->c2wdownout == 1)) { 
        sprintf(c2wout, "cos2w       ");
      } else {
        sprintf(c2wout, "cos2w^(%2d/%d)", -(matchptr->c2wupout), matchptr->c2wdownout);
      }
    }
    if ((matchptr->upin == 1) && (matchptr->downin == 1)) {
      sprintf(updownin, "       ");
    } else {
      sprintf(updownin, "(%2d/%2d)", matchptr->downin, matchptr->upin);
    }
    if ((matchptr->nbvupin == 0) && (matchptr->nbsupin == 0)) {
      sprintf(nbin, "      ");
    } else {
      if (matchptr->nbvupin == -1) {
        sprintf(nbin, "nbv   ");
      } else if (matchptr->nbvupin == 1) {
        sprintf(nbin, "nbv^-1");
      } else if (matchptr->nbsupin == -1) {
        sprintf(nbin, "nbs   ");
      } else {
        sprintf(nbin, "nbs^-1");
      }
    }
    if (matchptr->e2upin == 0) {
      sprintf(e2in, "        ");
    } else {
      sprintf(e2in, "2^(%2d/%d)", -(matchptr->e2upin), matchptr->e2downin);
    }
    if (matchptr->piupin == 0) {
      sprintf(piin, "         ");
    } else if ((matchptr->piupin == -1) && (matchptr->pidownin == 1)) {
      sprintf(piin, "pi       ");
    } else {
      sprintf(piin, "pi^(%2d/%d)", -(matchptr->piupin), matchptr->pidownin);
    }
    if (matchptr->aupin == 0) {
      sprintf(ain, "        ");
    } else if ((matchptr->aupin == -1) && (matchptr->adownin == 1)) {
      sprintf(ain, "a       ");
    } else {
      sprintf(ain, "a^(%2d/%d)", -(matchptr->aupin), matchptr->adownin);
    }
    if (matchptr->massratio == 0) {
      sprintf(massstr,       "    M/mp    ");
      sprintf(massstrinv,    "    mp/M    ");
    } else if (matchptr->massratio == 1) {
      sprintf(massstr,       "    M/v     ");
      sprintf(massstrinv,    "    v/M     ");
    } else if (matchptr->massratio == 2) {
      sprintf(massstr,       "    M/mZ    ");
      sprintf(massstrinv,    "    mZ/M    ");
    } else if (matchptr->massratio == 3) {
      sprintf(massstr,       "    M/mW    ");
      sprintf(massstrinv,    "    mW/M    ");
    } else if (matchptr->massratio == 4) {
      sprintf(massstr,       "    M/mH0   ");
      sprintf(massstrinv,    "    mH0/M   ");
    }             
    if (invexp == 1) {
      sprintf(leftformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  massstr);
    } else if (invexp == -1) {
      sprintf(leftformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv);
    } else if (invexp > 1) {
      sprintf(leftformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstr, invexp);
    } else {
      sprintf(leftformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv, -invexp);
    }

    matchptr=middlematchptr;
    invexp=middleinvexp;
    // note: all terms except matchup and matchdown are inverted here as they represent offsets to the real coefficient
    upout=matchptr->matchup * matchptr->downout;
    downout=matchptr->matchdown * matchptr->upout;
    if ((upout == 1) && (downout == 1)) {
      sprintf(updownout, "       ");
    } else {
      sprintf(updownout, "(%2d/%2d)", upout, downout);
    }
    if (matchptr->e2upout == 0) {
      sprintf(e2out, "        ");
    } else {
      sprintf(e2out, "2^(%2d/%d)", -(matchptr->e2upout), matchptr->e2downout);
    }
    if (matchptr->piupout == 0) {
      sprintf(piout, "         ");
    } else {
      if ((matchptr->piupout == -1) && (matchptr->pidownout == 1)) {
        sprintf(piout, "pi       ");
      } else {
        sprintf(piout, "pi^(%2d/%d)", -(matchptr->piupout), matchptr->pidownout);
      }
    }
    if (matchptr->aupout == 0) {
      sprintf(aout, "        ");
    } else {
      if ((matchptr->aupout == -1) && (matchptr->adownout == 1)) {
        sprintf(aout, "a       ");
      } else {
        sprintf(aout, "a^(%2d/%d)", -(matchptr->aupout), matchptr->adownout);
      }
    }
    if (matchptr->s2wupout == 0) {
      sprintf(s2wout, "            ");
    } else {
      if ((matchptr->s2wupout == -1) && (matchptr->s2wdownout == 1)) {
        sprintf(s2wout, "sin2w       ");
      } else {
        sprintf(s2wout, "sin2w^(%2d/%d)", -(matchptr->s2wupout), matchptr->s2wdownout);
      }
    }
    if (matchptr->c2wupout == 0) { 
      sprintf(c2wout, "            ");
    } else {
      if ((matchptr->c2wupout == -1) && (matchptr->c2wdownout == 1)) { 
        sprintf(c2wout, "cos2w       ");
      } else {
        sprintf(c2wout, "cos2w^(%2d/%d)", -(matchptr->c2wupout), matchptr->c2wdownout);
      }
    }
    if ((matchptr->upin == 1) && (matchptr->downin == 1)) {
      sprintf(updownin, "       ");
    } else {
      sprintf(updownin, "(%2d/%2d)", matchptr->downin, matchptr->upin);
    }
    if ((matchptr->nbvupin == 0) && (matchptr->nbsupin == 0)) {
      sprintf(nbin, "      ");
    } else {
      if (matchptr->nbvupin == -1) {
        sprintf(nbin, "nbv   ");
      } else if (matchptr->nbvupin == 1) {
        sprintf(nbin, "nbv^-1");
      } else if (matchptr->nbsupin == -1) {
        sprintf(nbin, "nbs   ");
      } else {
        sprintf(nbin, "nbs^-1");
      }
    }
    if (matchptr->e2upin == 0) {
      sprintf(e2in, "        ");
    } else {
      sprintf(e2in, "2^(%2d/%d)", -(matchptr->e2upin), matchptr->e2downin);
    }
    if (matchptr->piupin == 0) {
      sprintf(piin, "         ");
    } else if ((matchptr->piupin == -1) && (matchptr->pidownin == 1)) {
      sprintf(piin, "pi       ");
    } else {
      sprintf(piin, "pi^(%2d/%d)", -(matchptr->piupin), matchptr->pidownin);
    }
    if (matchptr->aupin == 0) {
      sprintf(ain, "        ");
    } else if ((matchptr->aupin == -1) && (matchptr->adownin == 1)) {
      sprintf(ain, "a       ");
    } else {
      sprintf(ain, "a^(%2d/%d)", -(matchptr->aupin), matchptr->adownin);
    }
    if (matchptr->massratio == 0) {
      sprintf(massstr,       "    M/mp    ");
      sprintf(massstrinv,    "    mp/M    ");
    } else if (matchptr->massratio == 1) {
      sprintf(massstr,       "    M/v     ");
      sprintf(massstrinv,    "    v/M     ");
    } else if (matchptr->massratio == 2) {
      sprintf(massstr,       "    M/mZ    ");
      sprintf(massstrinv,    "    mZ/M    ");
    } else if (matchptr->massratio == 3) {
      sprintf(massstr,       "    M/mW    ");
      sprintf(massstrinv,    "    mW/M    ");
    } else if (matchptr->massratio == 4) {
      sprintf(massstr,       "    M/mH0   ");
      sprintf(massstrinv,    "    mH0/M   ");
    }
    if (invexp == 1) {
      sprintf(middleformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  massstr);
    } else if (invexp == -1) {
      sprintf(middleformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv);
    } else if (invexp > 1) {
      sprintf(middleformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstr, invexp);
    } else {
      sprintf(middleformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv, -invexp);
    }

    matchptr=rightmatchptr;
    invexp=rightinvexp;
    // note: all terms except matchup and matchdown are inverted here as they represent offsets to the real coefficient
    upout=matchptr->matchup * matchptr->downout;
    downout=matchptr->matchdown * matchptr->upout;
    if ((upout == 1) && (downout == 1)) {
      sprintf(updownout, "       ");
    } else {
      sprintf(updownout, "(%2d/%2d)", upout, downout);
    }
    if (matchptr->e2upout == 0) {
      sprintf(e2out, "        ");
    } else {
      sprintf(e2out, "2^(%2d/%d)", -(matchptr->e2upout), matchptr->e2downout);
    }
    if (matchptr->piupout == 0) {
      sprintf(piout, "         ");
    } else {
      if ((matchptr->piupout == -1) && (matchptr->pidownout == 1)) {
        sprintf(piout, "pi       ");
      } else {
        sprintf(piout, "pi^(%2d/%d)", -(matchptr->piupout), matchptr->pidownout);
      }
    }
    if (matchptr->aupout == 0) {
      sprintf(aout, "        ");
    } else {
      if ((matchptr->aupout == -1) && (matchptr->adownout == 1)) {
        sprintf(aout, "a       ");
      } else {
        sprintf(aout, "a^(%2d/%d)", -(matchptr->aupout), matchptr->adownout);
      }
    }
    if (matchptr->s2wupout == 0) {
      sprintf(s2wout, "            ");
    } else {
      if ((matchptr->s2wupout == -1) && (matchptr->s2wdownout == 1)) {
        sprintf(s2wout, "sin2w       ");
      } else {
        sprintf(s2wout, "sin2w^(%2d/%d)", -(matchptr->s2wupout), matchptr->s2wdownout);
      }
    }
    if (matchptr->c2wupout == 0) { 
      sprintf(c2wout, "            ");
    } else {
      if ((matchptr->c2wupout == -1) && (matchptr->c2wdownout == 1)) { 
        sprintf(c2wout, "cos2w       ");
      } else {
        sprintf(c2wout, "cos2w^(%2d/%d)", -(matchptr->c2wupout), matchptr->c2wdownout);
      }
    }
    if ((matchptr->upin == 1) && (matchptr->downin == 1)) {
      sprintf(updownin, "       ");
    } else {
      sprintf(updownin, "(%2d/%2d)", matchptr->downin, matchptr->upin);
    }
    if ((matchptr->nbvupin == 0) && (matchptr->nbsupin == 0)) {
      sprintf(nbin, "      ");
    } else {
      if (matchptr->nbvupin == -1) {
        sprintf(nbin, "nbv   ");
      } else if (matchptr->nbvupin == 1) {
        sprintf(nbin, "nbv^-1");
      } else if (matchptr->nbsupin == -1) {
        sprintf(nbin, "nbs   ");
      } else {
        sprintf(nbin, "nbs^-1");
      }
    }
    if (matchptr->e2upin == 0) {
      sprintf(e2in, "        ");
    } else {
      sprintf(e2in, "2^(%2d/%d)", -(matchptr->e2upin), matchptr->e2downin);
    }
    if (matchptr->piupin == 0) {
      sprintf(piin, "         ");
    } else if ((matchptr->piupin == -1) && (matchptr->pidownin == 1)) {
      sprintf(piin, "pi       ");
    } else {
      sprintf(piin, "pi^(%2d/%d)", -(matchptr->piupin), matchptr->pidownin);
    }
    if (matchptr->aupin == 0) {
      sprintf(ain, "        ");
    } else if ((matchptr->aupin == -1) && (matchptr->adownin == 1)) {
      sprintf(ain, "a       ");
    } else {
      sprintf(ain, "a^(%2d/%d)", -(matchptr->aupin), matchptr->adownin);
    }
    if (matchptr->massratio == 0) {
      sprintf(massstr,       "    M/mp    ");
      sprintf(massstrinv,    "    mp/M    ");
    } else if (matchptr->massratio == 1) {
      sprintf(massstr,       "    M/v     ");
      sprintf(massstrinv,    "    v/M     ");
    } else if (matchptr->massratio == 2) {
      sprintf(massstr,       "    M/mZ    ");
      sprintf(massstrinv,    "    mZ/M    ");
    } else if (matchptr->massratio == 3) {
      sprintf(massstr,       "    M/mW    ");
      sprintf(massstrinv,    "    mW/M    ");
    } else if (matchptr->massratio == 4) {
      sprintf(massstr,       "    M/mH0   ");
      sprintf(massstrinv,    "    mH0/M   ");
    }             
    if (invexp == 1) {
      sprintf(rightformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain,  massstr);
    } else if (invexp == -1) {
      sprintf(rightformulastr, "'%s %s %s %s %s %s  %s %s %s %s %s %s       '", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv);
    } else if (invexp > 1) {
      sprintf(rightformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstr, invexp);
    } else {
      sprintf(rightformulastr, "'%s %s %s %s %s %s (%s %s %s %s %s %s)^(1/%2d)'", updownout, e2out, piout, aout, s2wout, c2wout, updownin, nbin, e2in, piin, ain, massstrinv, -invexp);
    }
/* end create formula output */

    sprintf(outstr01, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 01, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr01);
    sprintf(outstr02, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 02, |Parameter   ||         Value         | Std. Err. | Rel. Err. ||       Reference       | Std. Err. | Rel. Err. ||    Diff.    | Rel. Diff.  | Used as input | Used as output |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr02);
    sprintf(outstr03, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 03, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr03);
    if (alluses->alpha_em == 1) {
      if (floata == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr04, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 04, | alpha_em   || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, alpha_out_c, alpha_out_error, alpha_out_relerror, alpha_ref, alpha_ref_error, alpha_ref_relerror, alpha_out_diff, alpha_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr04);
    } else {
      outstr04[0]=0;
    }
    if (alluses->v == 1) {
      if (floatv == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr05, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 06, | v          || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, v_out_c, v_out_error, v_out_relerror, v_ref, v_ref_error, v_ref_relerror, v_out_diff, v_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr05);
    } else {
      outstr05[0]=0;
    }
    if ((alluses->mz == 1) || (mwmzmode == 1)) {
      if (floatmz == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr06, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 08, | mZ         || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mz_out_c, mz_out_error, mz_out_relerror, mz_ref, mz_ref_error, mz_ref_relerror, mz_out_diff, mz_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr06);
    } else {
      outstr06[0]=0;
    }
    if (alluses->G == 1) {
      if (floatg == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr07, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 07, | G          || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, G_out_c, G_out_error, G_out_relerror, G_ref, G_ref_error, G_ref_relerror, G_out_diff, G_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr07);
    } else {
      outstr07[0]=0;
    }
    if ((alluses->mw == 1) || (mwmzmode)) {
      if (floatmw == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr08, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 09, | mW         || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mw_out_c, mw_out_error, mw_out_relerror, mw_ref, mw_ref_error, mw_ref_relerror, mw_out_diff, mw_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr08);
    } else {
      outstr08[0]=0;
    }
    if (alluses->sin2w == 1) {
      if (floatsin2w == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr09, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 11, | sin2w      || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, sin2w_out_c, sin2w_out_error, sin2w_out_relerror, sin2w_ref, sin2w_ref_error, sin2w_ref_relerror, sin2w_out_diff, sin2w_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr09);
    } else {
      outstr09[0]=0;
    }
    if (alluses->mh0 == 1) {
      if (floatmh0 == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr10, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 10, | mH0        || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mh0_out_c, mh0_out_error, mh0_out_relerror, mh0_ref, mh0_ref_error, mh0_ref_relerror, mh0_out_diff, mh0_out_reldiff, usedasinput, usedasoutput);
      printf("%s\n", outstr10);
    } else {
      outstr10[0]=0;
    }
      if (floatme == 1 ) {
        sprintf(usedasoutput, "*");
        sprintf(usedasinput, " ");
      } else {
        sprintf(usedasinput, "*");
        sprintf(usedasoutput, " ");
      }
      sprintf(outstr11, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 12, | Electron   || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, me_out_c, me_out_error, me_out_relerror, me_ref, me_ref_error, me_ref_relerror, me_out_diff, me_out_reldiff, usedasinput, usedasoutput);
    printf("%s\n", outstr11);
    if (floatmu == 1) {
      sprintf(usedasinput, " ");
      sprintf(usedasoutput, "*");
    } else {
      sprintf(usedasinput, "*");
      sprintf(usedasoutput, " ");
    }
      sprintf(outstr12, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 13, | Muon       || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, mu_out_c, mu_out_error, mu_out_relerror, mu_ref, mu_ref_error, mu_ref_relerror, mu_out_diff, mu_out_reldiff, usedasinput, usedasoutput);
    printf("%s\n", outstr12);
    if (floattau == 1) {
      sprintf(usedasinput, " ");
      sprintf(usedasoutput, "*");
    } else {
      sprintf(usedasinput, "*");
      sprintf(usedasoutput, " ");
    }
      sprintf(outstr13, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 14, | Tau        || %.15e | %.3e | %.3e || %.15e | %.3e | %.3e || %11.4e | %11.4e |       %s       |       %s        |", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, tau_out_c, tau_out_error, tau_out_relerror, tau_ref, tau_ref_error, tau_ref_relerror, tau_out_diff, tau_out_reldiff, usedasinput, usedasoutput);
    printf("%s\n", outstr13);
    sprintf(outstr14, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 15, +------------++-----------------------+-----------------------++-----------------------+-----------+-----------++-------------+-------------+---------------+----------------+", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr14);
    sprintf(outstr15, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 16, L-M+R-1=0", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash);
    printf("%s\n", outstr15);
    sprintf(outstr16, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 17, L=%s", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, leftformulastr);
    printf("%s\n", outstr16);
    sprintf(outstr17, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 18, M=%s", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, middleformulastr);
    printf("%s\n", outstr17);
    sprintf(outstr18, "result, %.4f, %3d, %3d, %s, M%d%d%d, %12lld, 19, R=%s", combinedscore, symmetry, complexity, exponents, leftmatchptr->massratio, middlematchptr->massratio, rightmatchptr->massratio, resulthash, rightformulastr);
    printf("%s\n", outstr18);
    fflush(stdout);
#ifdef UPLOAD
    if (complexity <= maxcomplexity) {
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr01, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr02, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr03, 320));
      system(execstr);
      if (outstr04[0] != 0) {
        sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr04, 320));
        system(execstr);
      }
      if (outstr05[0] != 0) {
        sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr05, 320));
        system(execstr);
      }
      if (outstr06[0] != 0) {
        sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr06, 320));
        system(execstr);
      }
      if (outstr07[0] != 0) {
        sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr07, 320));
        system(execstr);
      }
      if (outstr08[0] != 0) {
        sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr08, 320));
        system(execstr);
      }
      if (outstr09[0] != 0) {
        sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr09, 320));
        system(execstr);
      }
      if (outstr10[0] != 0) {
        sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr10, 320));
        system(execstr);
      }
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr11, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr12, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr13, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr14, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr15, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr16, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr17, 320));
      system(execstr);
      sprintf(execstr, "curl -s \"%s/%s\" > /dev/null 2>&1\n", UPLOAD_URL, underscore(outstr18, 320));
      system(execstr);
    } // end if complexity
#endif
#ifndef ALWAYS_SHOW_RESULTS
  } // end if score
#endif


  return(precision_last);
}
