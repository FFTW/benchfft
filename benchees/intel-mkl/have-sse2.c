/* silly routine to test whether the kernel supports SSE2 */

#include <signal.h>
#include <setjmp.h>

static jmp_buf jb;

static void sighandler(int x)
{
     longjmp(jb, 1);
}

int have_sse2(void)
{
     signal(SIGILL, sighandler);
     if (setjmp(jb)) {
	  return 0;
     } else {
/*	  __asm__ __volatile__ ("xorpd %xmm0, %xmm0"); */
	  __asm__ __volatile__ (".byte 0x66, 0x0f, 0x57, 0xc0");
	  return 1;
     }
}
