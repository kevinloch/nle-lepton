#ifndef GLOBAL_H
#define GLOBAL_H

#include "nle-lepton.h" // for _G selection

// reference values (CODATA 2018/derived values unless otherwise stated)
// exact values
extern const double c_ref;
extern const double h_ref;
extern const double hbar_ref;
extern const double e_ref;
extern const double ev_to_kg;
extern const double kg_to_ev;

// experimental values
extern const double alpha_ref;
extern const double alpha_ref_error;
extern const double me_ref;
extern const double me_ref_error;
extern const double mu_ref;
extern const double mu_ref_error;
extern const double tau_ref;
extern const double tau_ref_error;
extern const double v_ref;
extern const double v_ref_error;
extern const double mz_ref;
extern const double mz_ref_error;
extern const double mw_ref;
extern const double mw_ref_error;
extern const double mh0_ref;
extern const double mh0_ref_error;
extern const double sin2w_ref;
extern const double sin2w_ref_error;

#ifdef CODATA_G
extern const double G_ref;
extern const double G_ref_error;
#elif defined ROSI_G
extern const double G_ref;
extern const double G_ref_error;
#elif defined WIDE_G
extern const double G_ref;
extern const double G_ref_error;
#endif

#endif // GLOBAL_H
