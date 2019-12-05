#include "nle-lepton.h"

// reference values (CODATA 2018/derived values unless otherwise stated)
// exact values
const double c_ref=             2.997924580000E+08;
const double h_ref=             6.62607015E-34;
const double hbar_ref=          1.05457181764616E-34;
//const double hbar_ref=h_ref/(2.0 * M_PI);
const double e_ref=             1.602176634E-19;
const double ev_to_kg=          1.78266192162790E-36;
const double kg_to_ev=          5.60958860380445E35;

// experimental values
const double alpha_ref=         7.2973525693E-3;
const double alpha_ref_error=   0.0000000011E-3;
const double me_ref=            0.51099895000E6;
const double me_ref_error=      0.00000000015E6;
const double mu_ref=          105.6583755E6;
const double mu_ref_error=      0.0000023E6;
const double tau_ref=        1776.86E6;
const double tau_ref_error=     0.12E6;
const double v_ref=           246.219651E9;      // derived from CODATA 2018 G0F
const double v_ref_error=       0.000063E9;      // derived from CODATA 2018 G0F
const double mz_ref=           91.1876E9;        // pdg
const double mz_ref_error=      0.0021E9;        // pdg
const double mw_ref=           80.379E9;         // pdg
const double mw_ref_error=      0.012E9;         // pdg
const double mh0_ref=         125.18E9;          // pdg
const double mh0_ref_error=     0.16E9;          // pdg
const double sin2w_ref=         0.22311;         // derived from pdg mw/mz
const double sin2w_ref_error=   0.00027;         // derived from pdg mw/mz

#ifdef CODATA_G
const double G_ref=             6.67430E-11;
const double G_ref_error=       0.00015E-11;
#elif defined ROSI_G
const double G_ref=             6.67191E-11;     // rosi
const double G_ref_error=       0.00099E-11;     // rosi
#elif defined WIDE_G
const double G_ref=             6.67400E-11;     // wide
const double G_ref_error=       0.00300E-11;     // wide
#endif
