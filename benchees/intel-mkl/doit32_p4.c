#define NAME "intel-mkl32-p4"
#define RIGHT_PROCESSOR (((cpuid_edx(1) & (1 << 26)) && have_sse2()))

#include "doit.c"
