/*
 * Copyright (c) 2002 Matteo Frigo
 * Copyright (c) 2002 Steven G. Johnson
 * Copyright (c) 2023 Paul Caprioli
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "config.h"
#include "bench.h"
#include <math.h>

#define DG unsigned short
#define ACC unsigned long
#define REAL double
#define BITS_IN_REAL 53 /* mantissa */

#define SHFT 16
#define RADIX 65536L
#define IRADIX (1.0 / RADIX)
#define LO(x) ((x) & (RADIX - 1))
#define HI(x) ((x) >> SHFT)
#define HI_SIGNED(x) \
   ((((x) + (ACC)(RADIX >> 1) * RADIX) >> SHFT) - (RADIX >> 1))
#define ZEROEXP (-32768)

#define LEN 10

typedef struct {
     short sign;
     short expt;
     DG d[LEN];
} N[1], COEFF[16];

#define EXA a->expt
#define EXB b->expt
#define EXC c->expt

#define AD a->d
#define BD b->d

#define SGNA a->sign
#define SGNB b->sign

static const N zero = {{ 1, ZEROEXP, {0} }};

static void cpy(const N a, N b)
{
     *b = *a;
}

static void fromreal(REAL x, N a)
{
     int i, e;

     cpy(zero, a);
     if (x == 0.0) return;
     
     if (x > 0) { SGNA = 1; }
     else       { SGNA = -1; x = -x; }

     e = 0;
     while (x >= 1.0) { x *= IRADIX; ++e; }
     while (x < IRADIX) { x *= RADIX; --e; }
     EXA = e;
     
     for (i = LEN - 1; i >= 0 && x != 0.0; --i) {
	  REAL y;

	  x *= RADIX;
	  y = floor(x);
	  AD[i] = (DG)y;
	  x -= y;
     }
}

static void fromshort(int x, N a)
{
     cpy(zero, a);

     if (x < 0) { x = -x; SGNA = -1; } 
     else { SGNA = 1; }
     EXA = 1;
     AD[LEN - 1] = x;
}

static void pack(DG *d, int e, int s, int l, N a)
{
     int i, j;

     for (i = l - 1; i >= 0; --i, --e) 
	  if (d[i] != 0) 
	       break;

     if (i < 0) {
	  /* number is zero */
	  cpy(zero, a);
     } else {
	  EXA = e;
	  SGNA = s;

	  if (i >= LEN - 1) {
	       for (j = LEN - 1; j >= 0; --i, --j)
		    AD[j] = d[i];
	  } else {
	       for (j = LEN - 1; i >= 0; --i, --j)
		    AD[j] = d[i];
	       for ( ; j >= 0; --j)
		    AD[j] = 0;
	  }
     }
}


/* compare absolute values */
static int abscmp(const N a, const N b)
{
     int i;
     if (EXA > EXB) return 1;
     if (EXA < EXB) return -1;
     for (i = LEN - 1; i >= 0; --i) {
	  if (AD[i] > BD[i])
	       return 1;
	  if (AD[i] < BD[i])
	       return -1;
     }
     return 0;
}

static int eq(const N a, const N b)
{
     return (SGNA == SGNB) && (abscmp(a, b) == 0);
}

/* add magnitudes, for |a| >= |b| */
static void addmag0(int s, const N a, const N b, N c)
{
     int ia, ib;
     ACC r = 0;
     DG d[LEN + 1];

     for (ia = 0, ib = EXA - EXB; ib < LEN; ++ia, ++ib) {
	  r += (ACC)AD[ia] + (ACC)BD[ib];
	  d[ia] = LO(r);
	  r = HI(r);
     }     
     for (; ia < LEN; ++ia) {
	  r += (ACC)AD[ia];
	  d[ia] = LO(r);
	  r = HI(r);
     }
     d[ia] = LO(r);
     pack(d, EXA + 1, s * SGNA, LEN + 1, c);
}

static void addmag(int s, const N a, const N b, N c)
{
     if (abscmp(a, b) > 0) addmag0(1, a, b, c); else addmag0(s, b, a, c);
}

/* subtract magnitudes, for |a| >= |b| */
static void submag0(int s, const N a, const N b, N c)
{
     int ia, ib;
     ACC r = 0;
     DG d[LEN];

     for (ia = 0, ib = EXA - EXB; ib < LEN; ++ia, ++ib) {
	  r += (ACC)AD[ia] - (ACC)BD[ib];
	  d[ia] = LO(r);
	  r = HI_SIGNED(r);
     }     
     for (; ia < LEN; ++ia) {
	  r += (ACC)AD[ia];
	  d[ia] = LO(r);
	  r = HI_SIGNED(r);
     }

     pack(d, EXA, s * SGNA, LEN, c);
}

static void submag(int s, const N a, const N b, N c)
{
     if (abscmp(a, b) > 0) submag0(1, a, b, c); else submag0(s, b, a, c);
}

/* c = a + b */
static void add(const N a, const N b, N c)
{
     if (SGNA == SGNB) addmag(1, a, b, c); else submag(1, a, b, c);
}

static void sub(const N a, const N b, N c)
{
     if (SGNA == SGNB) submag(-1, a, b, c); else addmag(-1, a, b, c);
}

static void mul(const N a, const N b, N c)
{
     DG d[2 * LEN];
     int i, j, k;
     ACC r;

     for (i = 0; i < LEN; ++i)
	  d[2 * i] = d[2 * i + 1] = 0;

     for (i = 0; i < LEN; ++i) {
	  ACC ai = AD[i];
	  if (ai) {
	       r = 0;
	       for (j = 0, k = i; j < LEN; ++j, ++k) {
		    r += ai * (ACC)BD[j] + (ACC)d[k];
		    d[k] = LO(r);
		    r = HI(r);
	       }
	       d[k] = LO(r);
	  }
     }

     pack(d, EXA + EXB, SGNA * SGNB, 2 * LEN, c);
}

static REAL toreal(const N a)
{
     REAL h, l, f;
     int i, bits;
     ACC r;
     DG sticky;

     if (EXA != ZEROEXP) {
	  f = IRADIX;
	  i = LEN;

	  bits = 0;
	  h = (r = AD[--i]) * f; f *= IRADIX;
	  for (bits = 0; r > 0; ++bits)
	       r >>= 1;

	  /* first digit */
	  while (bits + SHFT <= BITS_IN_REAL) {
	       h += AD[--i] * f;  f *= IRADIX; bits += SHFT;
	  }

	  /* guard digit (leave one bit for sticky bit, hence `<' instead
	     of `<=') */
	  bits = 0; l = 0.0;
	  while (bits + SHFT < BITS_IN_REAL) {
	       l += AD[--i] * f;  f *= IRADIX; bits += SHFT;
	  }
	  
	  /* sticky bit */
	  sticky = 0;
	  while (i > 0) 
	       sticky |= AD[--i];

	  if (sticky)
	       l += (RADIX / 2) * f;

	  h += l;

	  for (i = 0; i < EXA; ++i) h *= (REAL)RADIX;
	  for (i = 0; i > EXA; --i) h *= IRADIX;
	  if (SGNA == -1) h = -h;
	  return h;
     } else {
	  return 0.0;
     }
}

static void neg(N a)
{
     SGNA = -SGNA;
}

static void inv(const N a, N x)
{
     N w, z, one, two;

     fromreal(1.0 / toreal(a), x); /* initial guess */
     fromshort(1, one);
     fromshort(2, two);

     for (;;) {
	  /* Newton */
	  mul(a, x, w);
	  sub(two, w, z);
	  if (eq(one, z)) break;
	  mul(x, z, x);
     }
}

static const COEFF sinCoeff = {
     {+1, 1,
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
     {-1, 0,
       {43670, 43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690, 10922}},
     {+1, 0,
       {5348, 8738, 8738, 8738, 8738, 8738, 8738, 8738, 8738, 546}},
     {-1, 0,
       {53291, 3325, 208, 53261, 3328, 208, 53261, 3328, 208, 13}},
     {+1, -1,
       {47000, 63602, 58211, 40962, 37148, 21891, 26399, 44430, 51018, 11835}},
     {-1, -1,
       {18840, 1735, 27570, 57567, 8082, 16285, 5006, 40917, 39189, 107}},
     {+1, -2,
       {31480, 49626, 7441, 37348, 36497, 7193, 19429, 17256, 12445, 45202}},
     {-1, -2,
       {23808, 49206, 62861, 20482, 11094, 36547, 49400, 14749, 16287, 215}},
     {+1, -3,
       {32768, 43479, 33807, 9268, 45296, 12316, 21337, 34154, 15233, 51862}},
     {-1, -3,
       {38272, 22718, 30931, 9108, 35035, 9759, 2745, 13322, 42202, 151}},
     {+1, -4,
       {32768, 48630, 341, 31983, 25832, 49254, 50587, 29653, 15323, 23662}},
     {-1, -4,
       {52320, 13135, 19628, 41005, 37680, 3955, 28399, 9772, 50024, 46}},
     {+1, -5,
       {16384, 10745, 42516, 30975, 23226, 4984, 5160, 5212, 52445, 5107}},
     {-1, -5,
       {34508, 16922, 59733, 2551, 12116, 51487, 30682, 28278, 18092, 7}},
     {+1, -6,
       {46592, 62767, 39438, 40391, 47793, 29729, 64929, 38710, 15746, 587}},
     {-1, -7,
       {0, 58473, 50085, 27955, 51332, 49178, 54566, 54370, 45626, 41189}}
};

static void msin(const N a, N b) {
     N a2, r, s;
     int i;

     mul(a, a, a2);

     cpy(&sinCoeff[15], r);
     for (i = 14; i > 1; --i) {
          mul(a2, r, r);
          add(&sinCoeff[i], r, r);
     }
     mul(a2, r, r);
     mul(a2, r, r);
     mul(&sinCoeff[1], a2, a2);
     add(&sinCoeff[0], a2, s);
     sub(s, &sinCoeff[0], b);
     sub(a2, b, b);
     add(b, r, b);
     add(b, s, b);
     mul(a, b, b);
}

static const COEFF cosCoeff = {
     {+1, 1,
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
     {-1, 0,
       {0, 0, 0, 0, 0, 0, 0, 0, 0, 32768}},
     {+1, 0,
       {3160, 43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690, 2730}},
     {-1, 0,
       {27406, 23246, 1456, 45147, 23301, 1456, 45147, 23301, 1456, 91}},
     {+1, 0,
       {20393, 63968, 40985, 6657, 416, 40986, 6657, 416, 40986, 1}},
     {-1, -1,
       {55296, 12180, 19178, 43417, 23375, 48064, 55068, 56871, 37869, 1183}},
     {+1, -1,
       {41992, 41980, 39024, 15712, 27980, 55970, 27723, 30716, 63340, 8}},
     {-1, -2,
       {63488, 19170, 64630, 58133, 35313, 23919, 20112, 24638, 47700, 3228}},
     {+1, -2,
       {57504, 36964, 35902, 16690, 12625, 35052, 56335, 62361, 29689, 13}},
     {-1, -3,
       {40960, 36790, 33204, 39943, 51866, 53802, 52157, 56510, 15409, 2881}},
     {+1, -3,
       {18876, 47799, 48739, 29224, 14511, 12300, 32905, 59648, 38154, 7}},
     {-1, -4,
       {9216, 35552, 65153, 2117, 3423, 4047, 60615, 28157, 36443, 1075}},
     {+1, -4,
       {10087, 7554, 40990, 25715, 64520, 58334, 61223, 406, 62159, 1}},
     {-1, -5,
       {54400, 58044, 40758, 7255, 18031, 22864, 33827, 39607, 29743, 196}},
     {+1, -6,
       {49152, 7848, 60251, 50609, 44223, 6151, 23760, 24871, 61915, 17029}},
     {-1, -6,
       {16224, 5731, 28430, 64684, 54368, 57035, 15183, 64415, 31352, 19}}
};

static void mcos(const N a, N b) {
     N a2, r, s;
     int i;

     mul(a, a, a2);

     cpy(&cosCoeff[15], r);
     for (i = 14; i > 1; --i) {
          mul(a2, r, r);
          add(&cosCoeff[i], r, r);
     }
     mul(a2, r, r);
     mul(a2, r, r);
     mul(&cosCoeff[1], a2, a2);
     add(&cosCoeff[0], a2, s);
     sub(s, &cosCoeff[0], b);
     sub(a2, b, b);
     add(b, r, b);
     add(b, s, b);
}

/* Pi/4.  In radix 64K, avoids using 16 bits to store the leading 6 in 2*Pi. */
static const N quarterPi = {
    {+1, 0,
      {19977, 10498, 7377, 32988, 25227, 50374, 49716, 8552, 55970, 51471}}
};

static void by2pi(REAL m, REAL n, N a)
{
     N b;

     fromreal(n, b);
     inv(b, a);
     fromreal(8*m, b);  /* Note that 8*m is exact in binary floating point. */
     mul(a, b, a);
     mul(quarterPi, a, a);
}

static void sin2pi(REAL m, REAL n, N a);
static void cos2pi(REAL m, REAL n, N a)
{
     N b;
     if (m < 0) cos2pi(-m, n, a);
     else if (m > n * 0.5) cos2pi(n - m, n, a);
     else if (m > n * 0.25) {sin2pi(m - n * 0.25, n, a); neg(a);}
     else if (m > n * 0.125) sin2pi(n * 0.25 - m, n, a);
     else { by2pi(m, n, b); mcos(b, a); }
}

static void sin2pi(REAL m, REAL n, N a)
{
     N b;
     if (m < 0)  {sin2pi(-m, n, a); neg(a);}
     else if (m > n * 0.5) {sin2pi(n - m, n, a); neg(a);}
     else if (m > n * 0.25) {cos2pi(m - n * 0.25, n, a);}
     else if (m > n * 0.125) {cos2pi(n * 0.25 - m, n, a);}
     else {by2pi(m, n, b); msin(b, a);}
}

/*----------------------------------------------------------------------*/
/* FFT stuff */

/* (r0 + i i0)(r1 + i i1) */
static void cmul(N r0, N i0, N r1, N i1, N r2, N i2)
{
     N s, t, q;
     mul(r0, r1, s);
     mul(i0, i1, t);
     sub(s, t, q);
     mul(r0, i1, s);
     mul(i0, r1, t);
     add(s, t, i2);
     cpy(q, r2);
}

/* (r0 - i i0)(r1 + i i1) */
static void cmulj(N r0, N i0, N r1, N i1, N r2, N i2)
{
     N s, t, q;
     mul(r0, r1, s);
     mul(i0, i1, t);
     add(s, t, q);
     mul(r0, i1, s);
     mul(i0, r1, t);
     sub(s, t, i2);
     cpy(q, r2);
}

static void mcexp(int m, int n, N r, N i)
{
     static int cached_n = -1;
     static N w[64][2];
     int k, j;
     if (n != cached_n) {
	  for (j = 1, k = 0; j < n; j += j, ++k) {
	       cos2pi(j, n, w[k][0]);
	       sin2pi(j, n, w[k][1]);
	  }
	  cached_n = n;
     }

     fromshort(1, r);
     fromshort(0, i);
     if (m > 0) {
	  for (k = 0; m; ++k, m >>= 1) 
	       if (m & 1)
		    cmul(w[k][0], w[k][1], r, i, r, i);
     } else {
	  m = -m;
	  for (k = 0; m; ++k, m >>= 1) 
	       if (m & 1)
		    cmulj(w[k][0], w[k][1], r, i, r, i);
     }
}

static void bitrev(int n, N *a)
{
     int i, j, m;
     for (i = j = 0; i < n - 1; ++i) {
	  if (i < j) {
	       N t;
	       cpy(a[2*i], t); cpy(a[2*j], a[2*i]); cpy(t, a[2*j]);
	       cpy(a[2*i+1], t); cpy(a[2*j+1], a[2*i+1]); cpy(t, a[2*j+1]);
	  }

	  /* bit reversed counter */
	  m = n; do { m >>= 1; j ^= m; } while (!(j & m));
     }
}

static void fft0(int n, N *a, int sign)
{
     int i, j, k;

     bitrev(n, a);
     for (i = 1; i < n; i = 2 * i) {
	  for (j = 0; j < i; ++j) {
	       N wr, wi;
	       mcexp(sign * (int)j, 2 * i, wr, wi);
	       for (k = j; k < n; k += 2 * i) {
		    N *a0 = a + 2 * k;
		    N *a1 = a0 + 2 * i;
		    N r0, i0, r1, i1, t0, t1, xr, xi;
		    cpy(a0[0], r0); cpy(a0[1], i0);
		    cpy(a1[0], r1); cpy(a1[1], i1);
		    mul(r1, wr, t0); mul(i1, wi, t1); sub(t0, t1, xr);
		    mul(r1, wi, t0); mul(i1, wr, t1); add(t0, t1, xi);
		    add(r0, xr, a0[0]);  add(i0, xi, a0[1]);
		    sub(r0, xr, a1[0]);  sub(i0, xi, a1[1]);
	       }
	  }
     }
}

/* a[2*k]+i*a[2*k+1] = exp(2*pi*i*k^2/(2*n)) */
static void bluestein_sequence(int n, N *a)
{
     int k, ksq, n2 = 2 * n;

     ksq = 1; /* (-1)^2 */
     for (k = 0; k < n; ++k) {
	  /* careful with overflow */
	  ksq = ksq + 2*k - 1; while (ksq > n2) ksq -= n2;
	  mcexp(ksq, n2, a[2*k], a[2*k+1]);
     }
}

static int pow2_atleast(int x)
{
     int h;
     for (h = 1; h < x; h = 2 * h)
	  ;
     return h;
}

static N *cached_bluestein_w = 0;
static N *cached_bluestein_y = 0;
static int cached_bluestein_n = -1;

static void bluestein(int n, N *a)
{
     int nb = pow2_atleast(2 * n);
     N *b = (N *)bench_malloc(2 * nb * sizeof(N));
     N *w = cached_bluestein_w;
     N *y = cached_bluestein_y;
     N nbinv;
     int i;

     fromreal(1.0 / nb, nbinv); /* exact because nb = 2^k */

     if (cached_bluestein_n != n) {
	  if (w) bench_free(w);
	  if (y) bench_free(y);
	  w = (N *)bench_malloc(2 * n * sizeof(N));
	  y = (N *)bench_malloc(2 * nb * sizeof(N));
	  cached_bluestein_n = n;
	  cached_bluestein_w = w;
	  cached_bluestein_y = y;

	  bluestein_sequence(n, w);
	  for (i = 0; i < 2*nb; ++i)  cpy(zero, y[i]);

	  for (i = 0; i < n; ++i) {
	       cpy(w[2*i], y[2*i]);
	       cpy(w[2*i+1], y[2*i+1]);
	  }
	  for (i = 1; i < n; ++i) {
	       cpy(w[2*i], y[2*(nb-i)]);
	       cpy(w[2*i+1], y[2*(nb-i)+1]);
	  }

	  fft0(nb, y, -1);
     }

     for (i = 0; i < 2*nb; ++i)  cpy(zero, b[i]);
     
     for (i = 0; i < n; ++i) 
	  cmulj(w[2*i], w[2*i+1], a[2*i], a[2*i+1], b[2*i], b[2*i+1]);

     /* scaled convolution b * y */
     fft0(nb, b, -1);

     for (i = 0; i < nb; ++i) 
	  cmul(b[2*i], b[2*i+1], y[2*i], y[2*i+1], b[2*i], b[2*i+1]);
     fft0(nb, b, 1);

     for (i = 0; i < n; ++i) {
	  cmulj(w[2*i], w[2*i+1], b[2*i], b[2*i+1], a[2*i], a[2*i+1]);
	  mul(nbinv, a[2*i], a[2*i]);
	  mul(nbinv, a[2*i+1], a[2*i+1]);
     }

     bench_free(b);
}

static void swapri(int n, N *a)
{
     int i;
     for (i = 0; i < n; ++i) {
	  N t;
	  cpy(a[2 * i], t);
	  cpy(a[2 * i + 1], a[2 * i]);
	  cpy(t, a[2 * i + 1]);
     }
}

static void fft1(int n, N *a, int sign)
{
     if (power_of_two(n)) {
	  fft0(n, a, sign);
     } else {
	  if (sign == 1) swapri(n, a);
	  bluestein(n, a);
	  if (sign == 1) swapri(n, a);
     }
}

static void fromrealv(int n, bench_complex *a, N *b)
{
     int i;

     for (i = 0; i < n; ++i) {
         fromreal((REAL)c_re(a[i]), b[2 * i]);
         fromreal((REAL)c_im(a[i]), b[2 * i + 1]);
     }
}

static void compare(int n, N *a, N *b, double *err)
{
     int i;
     double e1, e2, einf;
     double n1, n2, ninf;

     e1 = e2 = einf = 0.0;
     n1 = n2 = ninf = 0.0;

#    define DO(x1, x2, xinf, var) { 			\
     double d = var;					\
     if (d < 0) d = -d;					\
     x1 += d; x2 += d * d; if (d > xinf) xinf = d;	\
}
	  
     for (i = 0; i < 2 * n; ++i) {
	  N dd;
	  sub(a[i], b[i], dd);
	  DO(n1, n2, ninf, toreal(a[i]));
	  DO(e1, e2, einf, toreal(dd));
     }

#    undef DO
     err[0] = e1 / n1;
     err[1] = sqrt(e2 / n2);
     err[2] = einf / ninf;
}

void fftaccuracy(int n, bench_complex *a, bench_complex *ffta,
		 int sign, double err[6])
{
     N *b = (N *)bench_malloc(2 * n * sizeof(N));
     N *fftb = (N *)bench_malloc(2 * n * sizeof(N));
     N mn, ninv;
     int i;

     fromreal(n, mn); inv(mn, ninv);

     /* forward error */
     fromrealv(n, a, b); fromrealv(n, ffta, fftb);
     fft1(n, b, sign);
     compare(n, b, fftb, err);

     /* backward error */
     fromrealv(n, a, b); fromrealv(n, ffta, fftb);
     for (i = 0; i < 2 * n; ++i) mul(fftb[i], ninv, fftb[i]);
     fft1(n, fftb, -sign);
     compare(n, b, fftb, err + 3);

     bench_free(fftb);
     bench_free(b);
}

void fftaccuracy_done(void)
{
     if (cached_bluestein_w) bench_free(cached_bluestein_w);
     if (cached_bluestein_y) bench_free(cached_bluestein_y);
     cached_bluestein_w = 0;
     cached_bluestein_y = 0;
     cached_bluestein_n = -1;
}
