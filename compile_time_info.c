#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        EOS_GAMMA=(5.0/3.0)\n"
"        TURB_DIFF_ENERGY\n"
"        TURB_DIFF_VELOCITY\n"
"        ADAPTIVE_GRAVSOFT_FORALL=1+2\n"
"        GALSF\n"
"        GALSF_SFR_CRITERION=(0+1+2+4)\n"
"        COOLING\n"
"        OUTPUT_TEMPERATURE\n"
"        MERGESPLIT_HARDCODE_MAX_MASS=(2.0e-10)\n"
"        MERGESPLIT_HARDCODE_MIN_MASS=(2.0e-11)\n"
"\n");
}
