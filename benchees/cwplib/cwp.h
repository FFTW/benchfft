/* Copyright (c) Colorado School of Mines, 1995.*/
/* All rights reserved.                       */

/* cwp.h - include file for general purpose CWP stuff */

#ifndef CWP_H
#define CWP_H


/* INCLUDES */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>


/* TYPEDEFS */
 
#ifndef __cplusplus /* if not C++, define the C struct complex */
typedef struct _complexStruct { /* complex number */
	float r,i;
} complex;
#else /* if C++, define the C++ class complex */
#include "Complex.h"
#endif /* C++ */


/* DEFINES */
/* uncomment the next block if you are installing */
/* under ultrix, but not using the GCC compiler */

/*
#ifdef ultrix
#define const
#endef
*/

#ifndef NULL
#define NULL	((void *)0)
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#endif
#ifndef SEEK_SET
#define SEEK_SET (0)
#endif
#ifndef SEEK_CUR
#define SEEK_CUR (1)
#endif
#ifndef SEEK_END
#define SEEK_END (2)
#endif
#ifndef PI
#define PI (3.141592653589793)
#endif
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#ifndef YES
#define YES (1)
#endif
#ifndef NO
#define NO (0)
#endif
#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define CLOSETO(x, y) (ABS((x) - (y)) <= FLT_EPSILON*ABS(y))
#define ISODD(n) ((n) & 01)
#define ISIZE sizeof(int)
#define FSIZE sizeof(float)
#define DSIZE sizeof(double)
#define	STREQ(s,t) (strcmp(s,t) == 0)
#define	STRLT(s,t) (strcmp(s,t) < 0)
#define	STRGT(s,t) (strcmp(s,t) > 0)
#define	DIM(a) (sizeof(a)/sizeof(a[0]))


/* FUNCTION PROTOTYPES */

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

/* allocate and free multi-dimensional arrays */
void *alloc1 (size_t n1, size_t size);
void *realloc1 (void *v, size_t n1, size_t size);
void **alloc2 (size_t n1, size_t n2, size_t size);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);
void free4 (void ****p);
int *alloc1int (size_t n1);
int *realloc1int (int *v, size_t n1);
int **alloc2int (size_t n1, size_t n2);
int ***alloc3int (size_t n1, size_t n2, size_t n3);
float *alloc1float (size_t n1);
float *realloc1float (float *v, size_t n1);
float **alloc2float (size_t n1, size_t n2);
float ***alloc3float (size_t n1, size_t n2, size_t n3);
double *alloc1double (size_t n1);
double *realloc1double (double *v, size_t n1);
double **alloc2double (size_t n1, size_t n2);
double ***alloc3double (size_t n1, size_t n2, size_t n3);
complex *alloc1complex (size_t n1);
complex *realloc1complex (complex *v, size_t n1);
complex **alloc2complex (size_t n1, size_t n2);
complex ***alloc3complex (size_t n1, size_t n2, size_t n3);
void free1int (int *p);
void free2int (int **p);
void free3int (int ***p);
void free1float (float *p);
void free2float (float **p);
void free3float (float ***p);
void free1double (double *p);
void free2double (double **p);
void free3double (double ***p);
void free1complex (complex *p);
void free2complex (complex **p);
void free3complex (complex ***p);

#ifndef __cplusplus /* if not C++, declare C complex functions */
/* complex number manipulation */
complex cadd (complex a, complex b);
complex csub (complex a, complex b);
complex cmul (complex a, complex b);
complex cdiv (complex a, complex b);
float rcabs (complex z);
complex cmplx (float re, float im);
complex conjg (complex z);
complex cneg (complex z);
complex cinv (complex z);
complex csqrt (complex z);
complex cexp (complex z);
complex crmul (complex a, float x);
#endif /* C++ */

/* big matrix handler */
void *bmalloc (int nbpe, int n1, int n2);
void bmfree (void *bm);
void bmread (void *bm, int dir, int k1, int k2, int n, void *v);
void bmwrite (void *bm, int dir, int k1, int k2, int n, void *v);

/* interpolation */
float fsinc (float x);
double dsinc (double x);
void mksinc (float d, int lsinc, float sinc[]);
void ints8r (int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[]);
void ints8c (int nxin, float dxin, float fxin, complex yin[], 
	complex yinl, complex yinr, int nxout, float xout[], complex yout[]);
void intt8r (int ntable, float table[][8],
	int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[]);
void intt8c (int ntable, float table[][8],
	int nxin, float dxin, float fxin, complex yin[], 
	complex yinl, complex yinr, int nxout, float xout[], complex yout[]);
void ress8r (int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, 
	int nxout, float dxout, float fxout, float yout[]);
void ress8c (int nxin, float dxin, float fxin, complex yin[], 
	complex yinl, complex yinr, 
	int nxout, float dxout, float fxout, complex yout[]);
void shfs8r (float dx, int nxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float fxout, float yout[]);
void xindex (int nx, float ax[], float x, int *index);
void intl2b (int nxin, float dxin, float fxin,
	int nyin, float dyin, float fyin, unsigned char *zin,
	int nxout, float dxout, float fxout,
	int nyout, float dyout, float fyout, unsigned char *zout);
void intlin (int nin, float xin[], float yin[], float yinl, float yinr,
	int nout, float xout[], float yout[]);
void intcub (int ideriv, int nin, float xin[], float ydin[][4],
	int nout, float xout[], float yout[]);
void cakima (int n, float x[], float y[], float yd[][4]);
void cmonot (int n, float x[], float y[], float yd[][4]);
void csplin (int n, float x[], float y[], float yd[][4]);
void yxtoxy (int nx, float dx, float fx, float y[], 
	int ny, float dy, float fy, float xylo, float xyhi, float x[]);

/* Butterworth filters */
void bfhighpass (int npoles, float f3db, int n, float p[], float q[]);
void bflowpass (int npoles, float f3db, int n, float p[], float q[]);
void bfdesign (float fpass, float apass, float fstop, float astop,
	int *npoles, float *f3db);

/* differentiator approximations */
void mkdiff (int n, float a, float h, int l, int m, float d[]);
void mkhdiff (float h, int l, float d[]);
void holbergd1 (float e, int n, float d[]);

/* general signal processing */
void conv (int lx, int ifx, float *x, int ly, int ify, float *y,
	int lz, int ifz, float *z);
void xcor (int lx, int ifx, float *x, int ly, int ify, float *y ,
	int lz, int ifz, float *z);
void hilbert (int n, float x[], float y[]);
void antialias (float frac, int phase, int n, float p[], float q[]);

/* Abel transformer */
void *abelalloc (int n);
void abelfree (void *at);
void abel (void *at, float f[], float g[]);

/* Hankel transformer */
void *hankelalloc (int nfft);
void hankelfree (void *ht);
void hankel0 (void *ht, float f[], float h[]);
void hankel1 (void *ht, float f[], float h[]);

/* sorting and searching */
void hpsort (int n, float a[]);
void qksort (int n, float a[]);
void qkfind (int m, int n, float a[]);
void qkisort (int n, float a[], int i[]);
void qkifind (int m, int n, float a[], int i[]);

/* statistics */
float quest (float p, int n, float x[]);
void *questalloc (float p, int n, float x[]);
float questupdate (void *q, int n, float x[]);
void questfree (void *q);

/* PC byte swapping */
void swap_short_2(short *tni2);
void swap_u_short_2(unsigned short *tni2);
void swap_int_4(int *tni4);
void swap_u_int_4(unsigned int *tni4);
void swap_long_4(long *tni4);
void swap_u_long_4(unsigned long *tni4);
void swap_float_4(float *tnf4);
void swap_double_8(double *tndd8);

/* Prime Factor FFTs */
int npfa (int nmin);
int npfao (int nmin, int nmax);
int npfar (int nmin);
int npfaro (int nmin, int nmax);
void pfacc (int isign, int n, complex z[]);
void pfarc (int isign, int n, float rz[], complex cz[]);
void pfacr (int isign, int n, complex cz[], float rz[]);
void pfa2cc (int isign, int idim, int n1, int n2, complex z[]);
void pfa2rc (int isign, int idim, int n1, int n2, float rz[], complex cz[]);
void pfa2cr (int isign, int idim, int n1, int n2, complex cz[], float rz[]);
void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);

/* BLAS (Basic Linear Algebra Subroutines adapted from LINPACK FORTRAN) */
int isamax (int n, float *sx, int incx);
float sasum (int n, float *sx, int incx);
void saxpy (int n, float sa, float *sx, int incx, float *sy, int incy);
void scopy (int n, float *sx, int incx, float *sy, int incy);
float sdot (int n, float *sx, int incx, float *sy, int incy);
float snrm2 (int n, float *sx, int incx);
void sscal (int n, float sa, float *sx, int incx);
void sswap (int n, float *sx, int incx, float *sy, int incy);
int idamax (int n, double *sx, int incx);
double dasum (int n, double *sx, int incx);
void daxpy (int n, double sa, double *sx, int incx, double *sy, int incy);
void dcopy (int n, double *sx, int incx, double *sy, int incy);
double ddot (int n, double *sx, int incx, double *sy, int incy);
double dnrm2 (int n, double *sx, int incx);
void dscal (int n, double sa, double *sx, int incx);
void dswap (int n, double *sx, int incx, double *sy, int incy);

/* LINPACK functions (adapted from LINPACK FORTRAN) */
void sgeco (float **a, int n, int *ipvt, float *rcond, float *z);
void sgefa (float **a, int n, int *ipvt, int *info);
void sgesl (float **a, int n, int *ipvt, float *b, int job);
void sqrdc (float **x, int n, int p, float *qraux, int *jpvt,
	float *work, int job);
void sqrsl (float **x, int n, int k, float *qraux,
	float *y, float *qy, float *qty,
	float *b, float *rsd, float *xb, int job, int *info);
void sqrst (float **x, int n, int p, float *y, float tol,
	float *b, float *rsd, int *k,
	int *jpvt, float *qraux, float *work);
void dgeco (double **a, int n, int *ipvt, double *rcond, double *z);
void dgefa (double **a, int n, int *ipvt, int *info);
void dgesl (double **a, int n, int *ipvt, double *b, int job);

/* other linear system solvers */
void stoepd (int n, double r[], double g[], double f[], double a[]);
void stoepf (int n, float r[], float g[], float f[], float a[]);
void vanded (int n, double v[], double b[], double x[]);
void vandef (int n, float v[], float b[], float x[]);
void tridif (int n, float a[], float b[], float c[], float r[], float u[]);
void tridid (int n, double a[], double b[], double c[], double r[], double u[]);
void tripd(float *d, float *e, float *b, int n);

/* root finding */
int mnewt (int maxiter, float ftol, float dxtol, int n, float *x, void *aux,
	void (*fdfdx)(int n, float *x, float *f, float **dfdx, void *aux));

/* transform rectangular => polar and polar => to rectangular coordinates */
void recttopolar ( int nx, float dx, float fx, int ny, float dy,
	float fy, float **p, int na, float da, float fa, int nr, float dr,
	float fr, float **q);
void polartorect ( int na, float da, float fa, int nr, float dr,
	float fr, float **q, int nx, float dx, float fx, int ny, float dy,
	float fy, float **p);

/* graphics utilities */
void rfwtva (int n, float z[], float zmin, float zmax, float zbase,
	int yzmin, int yzmax, int xfirst, int xlast,
	int wiggle, int nbpr, unsigned char *bits);
void scaxis (float x1, float x2, int *nxnum, float *dxnum, float *fxnum);
int yclip (int nx, float dx, float fx, float y[], float ymin, float ymax,
	float xc[], float yc[]);

/* special functions */
float airya (float x);
float airyb (float x);
float airyap (float x);
float airybp (float x);

/* timers */
float cpusec (void);
float cputime (void);
float wallsec (void);
float walltime (void);

/* pseudo-random numbers */
float franuni (void);
void sranuni (int seed);
float frannor (void);
void srannor (int seed);

/* miscellaneous */
void pp1d (FILE *fp, char *title, int lx, int ifx, float x[]);
void pplot1 (FILE *fp, char *title, int nx, float ax[]);

#ifdef __cplusplus /* if C++, end external linkage specification */
}
#endif

#endif /* CWP_H */
