#define NAME "intel-mkl32-f-p3"
#define RIGHT_PROCESSOR (cpuid_edx(1) & (1 << 25))

#include "doitf.c"
