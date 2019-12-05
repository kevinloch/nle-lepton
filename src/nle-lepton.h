#ifndef NLE_LEPTON_H
#define NLE_LEPTON_H

#define NLE_VERSION "3.41"

#define UPLOAD                                // upload results to web server
#define UPLOAD_URL "http://localhost/lepton"  // customize this for your results server.  See scripts/grep-logs.sh for sample script to extract results from web logs
#define SHOWSTATUS                            // print progress updates to stdout
#define NEGATIVEEXP                           // allow negative exponents
#define IGNORE_SMALL_UNCERTAINTIES            // use fixed reference values for electron mass and alpha_em in phase 2
//#define SIN2W                               // enable weak mixing angle terms (much slower)
//#define ALWAYS_SHOW_RESULTS                 // shows phase 2 results even if they don't match experimental uncertainties

// pick one G reference
#define CODATA_G
//#define ROSI_G
//#define WIDE_G // wide=large uncertainty

//#define DEBUG10   // show periodic status and warnings if process is slow
//#define DEBUG11   // show status on each progress step
//#define DEBUG12   // show data on each sample tried.  Caution, extremely noisy

//#define DEBUG20   // show periodic status and warnings if process is slow
//#define DEBUG21   // show status on each progress step
//#define DEBUG22   // show data on each sample tried. Caution, extremely noisy
//#define DEBUG23   // show even more data on each sample tried. Caution, extremely noisy

typedef struct {
  int alpha_em;
  int v;
  int G;
  int mz;
  int mw;
  int mh0;
  int sin2w;
} input_use;

typedef struct {
  int upin;
  int downin;
  int piupin;
  int pidownin;
  int aupin;
  int adownin;
  int e2upin;
  int e2downin;
  int nbvupin;
  int nbsupin;
  int multcomplexity;
  double mult[27];
  input_use uses;
} multipliers;

typedef struct {
  double tau_sample;
  double G_sample;
  double mp_sample;
  double mw_sample;
  double mz_sample;
  double mh0_sample;
  double sin2w_sample;
} random_input;

typedef struct {
  int invexp;
  int massratio;
  int matchup;
  int matchdown;
  int upin;
  int downin;
  int piupin;
  int pidownin;
  int aupin;
  int adownin;
  int e2upin;
  int e2downin;
  int nbvupin;
  int nbsupin;
  int upout;
  int downout;
  int piupout;
  int pidownout;
  int aupout;
  int adownout;
  int e2upout;
  int e2downout;
  int s2wupout;
  int s2wdownout;
  int c2wupout;
  int c2wdownout;
  int matchcomplexity;
  double matchmult;
  double match;
  long long matchhash;
  input_use uses;
} matches;

#endif // NLE_LEPTON_H
