/* Simple nrerror() function, needed by rlft3.c.  This was implemented
   without looking the Numerical Recipes nrerror() source code. */

#include <stdlib.h>
#include <stdio.h>

void nrerror(char error_text[])
{
     fprintf(stderr, error_text);
     exit(EXIT_FAILURE);
}
