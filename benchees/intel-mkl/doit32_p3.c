#define NAME "intel-mkl32-p3"
#define RIGHT_PROCESSOR (cpuid_edx(1) & (1 << 25))

#include "doit.c"
