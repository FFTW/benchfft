/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "fftw.h"

/* horrible hack for now */
#define problem fftw_problem
#include "problem.h"
#include "plan.h"
#include "planner.h"
#include "planner-memo.h"
#include "planner-naive.h"
#include "planner-score.h"
#include "planner-estimate.h"
#include "problem-ridft.h"
#undef problem

extern const char *fftw_version;
extern const char *fftw_cc;

static const char *mkvers(void)
{
     return fftw_version;
}

static const char *mkcc(void)
{
     return fftw_cc;
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "fftw3")
BENCH_DOCF("version", mkvers)
BENCH_DOCF("compiled-by", mkcc)
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (sizeof(fftw_real) == sizeof(bench_real) && 
	     p->kind == PROBLEM_COMPLEX);
}

static struct planner *planner;
static struct fftw_problem *problem;
static struct plan *plan;

void setup(struct problem *p)
{
     bench_real *ri, *ii, *ro, *io;
     BENCH_ASSERT(can_do(p));

     planner = fftw_planner_naive_make();
     planner = fftw_planner_memo_make(planner);
     fftw_configuration_dft_standard(planner);

     if (p->sign == -1) {
	  ri = p->in; ii = ri + 1; ro = p->out; io = ro + 1;
     } else {
	  ii = p->in; ri = ii + 1; io = p->out; ro = io + 1;
     }

     problem = fftw_problem_ridft_make_d(
	  fftw_io_tensor_create_rowmajor(p->rank, p->n, p->n, 2, 2),
	  fftw_io_tensor_create_1d(1, 0, 0),
	  ri, ii, ro, io);
     plan = planner->plan(planner, problem);
     BENCH_ASSERT(plan);
     plan->awake(plan, PLAN_AWAKE);
}

void doit(int iter, struct problem *p)
{
     int i;
     struct plan *PLAN = plan;
     struct fftw_problem *PROBLEM = problem;

     for (i = 0; i < iter; ++i) {
	  PLAN->solve(PLAN, PROBLEM);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     fftw_plan_destroy(plan);
     fftw_problem_destroy(problem);
     fftw_planner_destroy(planner);
}
