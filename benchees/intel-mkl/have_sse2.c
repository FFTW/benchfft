/* silly routine to test whether the kernel supports SSE2 */

#include <signal.h>
#include <setjmp.h>

static jmp_buf jb;

static void sighandler(int x)
{
     longjmp(jb, 1);
}

int have_sse(void)
{
     signal(SIGILL, sighandler);
     if (setjmp(jb)) {
	  return 0;
     } else {
	  __asm__ __volatile__ ("xorps %xmm0, %xmm0");
	  return 1;
     }
}
