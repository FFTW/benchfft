#include	"qftsubs.h"

/*	QFT constants	*/

#define	sec0	0.5000000000000000f
#define	sec1	0.5411961001461970f
#define	sec2	0.7071067811865475f
#define	sec3	1.3065629648763763f

#define	sqrt2	sec2

#define	pi	((double) 3.141592653589793)

/*	QFT externals	*/

static	float	**ws;	/*	Workspace pointers	*/
static	float	*sc;	/*	Secant table addr.	*/
static	long	*nc;	/*	Current QFT length	*/
static	long	*ic;	/*	Current secant inc.	*/
static	long	*mc;	/*	Current output pruning	*/
static	long	nn;		/*	Pruned QFT input length	*/
static	long	mm;		/*	Pruned QFT output length	*/
static	long	ii;		/*	Current QFT recursion level	*/

/*	Program: qftinit, Quick Fourier Transform	*/
/*	secant table and workspace initializer.		*/

int	qftinit(long n)

{

	long	 n2, i, m;
	double q;

/*	Compute some constants	*/

	m = 0;
	i = n;
	while(i >>= 1)
		m++;

	n2 = n >> 1;

/*	Allocate the QFT memory	*/

	if ((ws = (float **) malloc(sizeof(float *)*m + 
		 sizeof(long)*3*m + sizeof(float)*(3*(n + m)))) == NULL)
	{
		printf("\n\a Out of memory in qftinit().\n");
		return 0;
	}

	nc = (long *)&ws[m];
	ic = &nc[m];
	mc = &ic[m];
	sc = (float *)&mc[m];
	ws[0] = &sc[n2-1];

/*	Initialize tables	*/
	
	nc[0] = n;
	ic[0] = 1;
	mc[0] = m;

	for (i = 1; i < m; i++)
	{
		nc[i] = nc[i-1] >> 1;
		ic[i] = ic[i-1] << 1;
		
		if (i > 1)
			ws[i] = ws[i-1] + 2*(nc[i-1] + 2);
		else
			ws[i] = ws[i-1] + nc[i] + 1;	
	}	

/*	Compute scaled secants	*/
	
	q = pi / n;

	for (i = 0; i < n2-1; i++)
		sc[i] = 0.5 / cos(q*i);	

	return 1;

}

/*	Program: qftclos, Close a QFT session.		*/

void	qftclos()

{

/*	Deallocate memory space	*/

	free(ws);

}


/*	Program: rfqft, Real Forward Quick Fourier Transform	*/
/*	(uses secant table), length(x) = n = 2^m >= 16. Input	*/
/*	is only computed for x[0],...,x[k-1], M <= k even <= n.	*/
/*	If I == NULL then put spectrum in "packed" format in R.	*/

void	rfqft(float *R, float *I, float *x, long n, long k)

{

	float	*xx, *y;
	long	i, n2;

/*	Set up constants	*/

	if (k < MM)
		nn = MM;
	else
		nn = k;

	mm = n;
	n2 = n >> 1;
	xx = ws[ii = 0];

	for (i = 1; i < mc[0]; i++)
		mc[i] = n;

/*	Form even function	*/

	if (nn <= n2)
		y = x;
	else
	{			
		xx[0] = x[0];

		for (i = 1; i <= n-nn; i += 2)
		{
			xx[i] = x[i];
			xx[i+1] = x[i+1];
		}

		xx[i] = x[i] + x[n-i], i++;	

		for (;i < n2; i += 2)
		{
			xx[i] = x[i] + x[n-i];
			xx[i+1] = x[i+1] + x[n-(i+1)];
		}

		xx[i] = x[i], i++;

		for (;i < MM; i++)
			xx[i] = 0.0;
		
		y = xx;
	}

/*	Do a real DCT	*/

	rdct(R, y);

/*	Form negative odd function	*/

	if (nn <= n2)
	{
		xx[1] = -x[1];
 
		for (i = 2; i < nn; i += 2)
		{
			xx[i] = -x[i];
			xx[i+1] = -x[i+1];
		}
	}
	else
	{			
		for (i = 1; i <= n-nn; i += 2)
		{
			xx[i] = -x[i];
			xx[i+1] = -x[i+1];
		}

		xx[i] = x[n-i] - x[i], i++;	

		for (;i < n2; i += 2)
		{
			xx[i] = x[n-i] - x[i];
			xx[i+1] = x[n-(i+1)] - x[i+1];
		}
	}

	for (; i < MM; i++)	
		xx[i] = 0.0;

/*	Do a real DST	*/

	if (I == NULL)
		rdst(&R[n2+1], &xx[1]);
	else
	{
		rdst(&I[1], &xx[1]);
		I[0] = I[n2] = 0.0;
	}
		
}

/*	Program: cfqft, Complex Forward Quick Fourier Transform	*/
/*	(uses secant table), length(x) = n = 2^m >= 16. Input	*/
/*	is only computed for x[0],...,x[k-1], M <= k even <= n.	*/

void	cfqft(float *X, float *Y, float *x, float *y, long n, long k)

{

	float	rx, ix, ry, iy;
	long	i, n2, n4;

	n2 = n >> 1;
	n4 = n2 >> 1;

/*	Compute packed QFT of the real and imaginary parts	*/

	rfqft(X, NULL, x, n, k);
	rfqft(Y, NULL, y, n, k);

/*	Reverse the complex spectrum of the imaginary part	*/

	for (i = 1; i < n4; i++)
	{
		ix = X[i+n2];
		X[i+n2] = X[n-i];
		X[n-i] = ix;
		iy = Y[i+n2];
		Y[i+n2] = Y[n-i];
		Y[n-i] = iy;
	}		 	

/*	Combine spectra of real and imaginary parts	*/

	for (i = 1; i < n2; i++)
	{
		X[i] = (rx = X[i]) - (iy = Y[n-i]);
		Y[i] = (ix = X[n-i]) + (ry = Y[i]);
		X[n-i] = rx + iy;
		Y[n-i] = ry - ix;
	}

} 

/*	Program: riqft, Real Inverse Quick Fourier Transform	*/
/*	(uses secant table), length(x) = n = 2^m >= 16. Output	*/
/*	is only computed for x[0],...,x[k-1], M <= k even <= n.	*/
/*	If I == NULL then get spectrum in "packed" format in R.	*/

void	riqft(float *R, float *I, float *x, long n, long k)

{

	float	*xx, nr;
	long	i, j, n2;

/*	Set up constants	*/

	if (k < MM)
		mm = MM;
	else
		mm = k;

	nn = n;
	n2 = n >> 1;
	xx = ws[ii = 0];

	for (i = 1; i < mc[0]; i++)
	{
		if ((j = k >> (i-1)) < 1)
			j = 1;

		mc[i] = j;
	}

/*	Do a real DCT and DST	*/

	R[0] *= 0.5;
	R[n2] *= 0.5;
	rdct(x, R);
	R[0] *= 2.0;
	R[n2] *= 2.0;

	if (I == NULL)
		rdst(&xx[1], &R[n2+1]);
	else
		rdst(&xx[1], &I[1]);

/*	Complete the transform	*/

	nr = 2.0 / n;
	x[0] *= nr;
	
	if (mm < n2)
		for (i = 1; i < mm; i++)
			x[i] = (x[i] - xx[i])*nr;	
	else
	{
		for (i = 1; i <= n-mm; i++)
			x[i] = (x[i] - xx[i])*nr;	

		for (; i < n2; i++)
		{
			x[n-i] = (x[i] + xx[i])*nr;
			x[i] = (x[i] - xx[i])*nr;
		}

		x[i++] *= nr;
	}

	for (; i < MM; i++)	
		x[i] = 0.0;

} 

/*	Program: rdct, Real Discrete Cosine Transform	*/
/*	(secant table version), length(x) = n + 1.  On	*/
/*	entry: ii+1 = current recursion level, ws[ii+1]	*/
/*	= workspace pointer, nc[ii+1] = DFT length, and	*/
/*	ic[ii+1] = secant table index increment. All of	*/
/*	the above variables are externally declared.	*/

void	rdct(float *X, float *x)

{

	float	*XX, *xx, t0, t1, t2, t3, t4;
	long	i, j, k, n, nq, n2, n2p1;

/*	Increment level	& start	*/

	n = nc[++ii];
	xx = ws[ii];

/*	Last recursion level ?	*/

	if (n == N2)
	{
	/*	Form even function	*/

	 	xx[0] = x[0] + x[8];
		xx[1] = x[1] + x[7];
		xx[2] = x[2] + x[6];
		xx[3] = x[3] + x[5];

	/*	Save even DCT coefficients	*/

		t0 = xx[0] - x[4];
		t3 = xx[0] + x[4];
		t2 = xx[1] + xx[3];
		t1 = (xx[1] - xx[3])*sqrt2;
		t4 = t3 + xx[2];

		X[0] = t4 + t2;
		X[2] = t0 + t1;
		X[4] = t3 - xx[2];
		X[6] = t0 - t1;
		X[8] = t4 - t2;
	
	/*	Form twiddled odd function	*/

	 	xx[0] = (x[0] - x[8])*sec0;
		xx[1] = (x[1] - x[7])*sec1;
		xx[2] = (x[2] - x[6])*sec2;
		xx[3] = (x[3] - x[5])*sec3;

	/*	Save odd DCT coefficients	*/

		t1 = (xx[1] - xx[3])*sqrt2;	
		t2 = xx[1] + xx[3];
		t4 = xx[0] + xx[2];

		X[1] = t4 + t2 + (t0 = xx[0] + t1);
		X[3] = t0 + (t3 = xx[0] - xx[2]);
		X[5] = t3 + (t1 = xx[0] - t1);
		X[7] = t1 + t4 - t2; 
	}
	else
	{
		n2 = n >> 1;
		n2p1 = n2 + 1;
		nq = ic[ii];  

		if ((k = mc[ii]) > n2)
			k = n2;

	/*	Form even-odd functions	*/

		if (nn <= n2p1)
			for (i = j = 0; i < nn; i++, j += nq)
				xx[i+n2p1] = (xx[i] = x[i])*sc[j];
		else
		{			
			for (i = j = 0; i <= n-nn; i++, j += nq)
				xx[i+n2p1] = (xx[i] = x[i])*sc[j];

			for (; i < n2; i++, j+= nq)
			{
				xx[i] = x[i] + x[n-i];
				xx[i+n2p1] = (x[i] - x[n-i])*sc[j];
			}

			xx[i] = x[i];
			xx[i+n2p1] = 0.0;
		}

	/*	Do recursive DCTs	*/

		rdct((XX = &xx[n+2]), xx);
		rdct(&XX[n2p1], &xx[n2p1]);

	/*	Save DCT coefficients	*/

		t0 = XX[n2p1];

		for (i = 0; i < k; i += 2)
		{	
			X[j = i+i] = XX[i];
			X[j+1] = t0 + (t1 = XX[i+n2p1+1]);
			X[j+2] = XX[i+1];
			X[j+3] = t1 + (t0 = XX[i+n2p1+2]);
		}

		X[i+i] = XX[i];
	}

/*	Decrement level & exit	*/

	ii--;

}

/*	Program: rdst, Real Discrete Sine Transform		*/
/*	(secant table version), length(x) = n-1. On		*/
/*	entry: ii+1 = current recursion level, ws[ii+1]	*/
/*	= workspace pointer, nc[ii+1] = DFT length, and	*/
/*	ic[ii+1] = secant table index increment. All of	*/
/*	the above variables are externally declared.	*/

void	rdst(float *X, float *x)

{

	float	*XX, *xx, t0, t1, t2, t3;
	long	i, j, k, n, nq, n2, n2m1;

/*	Increment level & start	*/

	n = nc[++ii];
	xx = ws[ii];
	
/*	Last recursion level ?	*/

	if (n == N2)
	{
	/*	Form odd function	*/
				
		xx[0] = x[0] - x[6];
		xx[1] = x[1] - x[5];
		xx[2] = x[2] - x[4];

	/*	Save odd DST coefficients	*/

		t0 = (xx[0] + xx[2])*sqrt2;

		X[1] = t0 + xx[1];
		X[3] = xx[0] - xx[2];
		X[5] = t0 - xx[1];
	
	/*	Form twiddled even function	*/

		xx[0] = (x[0] + x[6])*sec1;
		xx[1] = (x[1] + x[5])*sec2;
		xx[2] = (x[2] + x[4])*sec3;

	/*	Save even DST coefficients	*/

		t0 = (xx[0] + xx[2])*sqrt2;
	
  		X[0] = (t1 = t0 + xx[1]) + x[3];
		X[2] = (t2 = xx[0] - xx[2]) + t1 - x[3];
		X[4] = (t0 = t0 - xx[1]) + t2 + x[3];
		X[6] = t0 - x[3];
	}
	else
	{
		n2 = n >> 1;
		n2m1 = n2 - 1;
		nq = ic[ii];

		if ((k = mc[ii]) > n2m1)
			k = n2m1;

	/*	Form even-odd functions	*/

		if (nn <= n2)
			for (i = 0, j = nq; i < nn-1; i++, j += nq)
				xx[i+n2m1] = (xx[i] = x[i])*sc[j];
		else
		{				
			for (i = 0, j = nq; i < n-nn; i++, j += nq)
				xx[i+n2m1] = (xx[i] = x[i])*sc[j];
				
			for (;i < n2m1; i++, j += nq)
			{
				xx[i] = x[i] - x[n-2-i];
				xx[i+n2m1] = (x[i] + x[n-2-i])*sc[j];
			}
		}
	
	/*	Do recursive DSTs	*/

		rdst((XX = &xx[n-2]), xx);
		rdst(&XX[n2m1], &xx[n2m1]);

	/*	Save DST coefficients	*/
		
		X[1] = XX[0];

		if (n2 >= nn)
		{
			t0 = X[0] = XX[n2m1];
 
			for (i = 1; i < k; i += 2)
			{
	 			X[j = i+i] = (t1 = XX[i+n2m1]) + t0;
				X[j+1] = XX[i];
	 			X[j+2] = (t0 = XX[i+n2m1+1]) + t1;
				X[j+3] = XX[i+1];
			}

			X[i+i] = t0;
		}
		else
		{
			X[0] = (t1 = XX[n2m1]) + (t0 = x[n2m1]);

			for (i = 1; i < k; i += 2)
			{
	 			X[j = i+i] = (t2 = XX[i+n2m1]) + t1 - t0;
				X[j+1] = XX[i];
	 			X[j+2] = (t1 = XX[i+n2m1+1]) + t2 + t0;
				X[j+3] = XX[i+1];
			}

			X[i+i] = t1 - t0;
		}
	}

/*	Decrement level	& exit	*/

	ii--;	

}




















