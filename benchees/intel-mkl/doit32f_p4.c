#define NAME "intel-mkl32-f-p4"
#define RIGHT_PROCESSOR (DOUBLE_PRECISION || ((cpuid_edx(1) & (1 << 26)) && have_sse2()))

#include "doitf.c"
