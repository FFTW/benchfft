#ifndef	__QFTSUBS_H

/*	System include files	*/

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>

/*	QFT function prototypes	*/

void	rfqft(float *R, float *I, float *x, long n, long k);
void	riqft(float *R, float *I, float *x, long n, long k);
void	cfqft(float *R, float *I, float *x, float *y, long n, long k);
void	rdct(float *X, float *x);
void	rdst(float *X, float *x);
int	qftinit(long n);
void	qftclos();

#define	N2	    8		/*	1/2 minimum QFT length	*/
#define	NN	    (N2*2)	/*	Minimum QFT length	*/
#define	MM	    (N2+1)	/*	Minimum pruned length	*/

#define	__QFTSUBS_H

#endif



