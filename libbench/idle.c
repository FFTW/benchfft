#include <utmpx.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <ctype.h>
#include <errno.h>

static time_t idletime(char *tty)
{
     struct stat sbuf;
     if (stat(tty, &sbuf) != 0)
	  return 0;
     return time(NULL) - sbuf.st_atime;
}

static int idlep(int idlemin)
{
     struct utmpx *u;
     int res = 1;
     size_t i;
     time_t idle;
     char tty[5 + sizeof u->ut_line + 1] = "/dev/";

     setutxent();
     for (;;) {
	  u = getutxent();
	  if (!u)
	       break;
	  if (u->ut_type != USER_PROCESS)
	       continue;

	  if ((!u->ut_pid) || (kill(u->ut_pid, 0) == ESRCH)) /* no login process */
	       continue;

	  for (i = 0; i < sizeof(u->ut_line); i++) /* clean up tty if garbled */
	       if (isalnum(u->ut_line[i]) || (u->ut_line[i] == '/'))
		    tty[i + 5] = u->ut_line[i];
	       else
		    tty[i + 5] = '\0';

	  idle = idletime(tty);
	  if (idle < idlemin) {
	       res = 0;
	       break;
	  }
     }
     endutxent();
     return res;
}

void wait_for_idle(int idlemin)
{
     if (idlemin > 0)
	  while (!idlep(idlemin))
	       sleep(1 + idlemin / 2);
}
