#define NAME "ffte-v-sse3"
#define FFTE_VECTOR 1
#define DOC_ITEMS BENCH_DOC("notes", "vector-machine version") BENCH_DOC("notes", "using SSE3 kernels") BENCH_DOC("language", "C")

#include "doit-ffte.c"
