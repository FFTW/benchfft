/* Modified 8/18/98 by Steven G. Johnson (stevenj@alum.mit.edu) to
   work with benchFFT (http://www.fftw.org/benchfft). */

#include "bench-user.h"
   
/* Received in personal communication with the author:
   
   To: benchfft@fftw.org
   Subject: some FFT/FHT code
   From: cbloodworth@juno.com (Carey E Bloodworth)
   Date: Thu, 13 Aug 1998 18:01:06 EDT
   
   A few days ago, I checked your benchmark results looking for any timings
   for virtual memory FFTs.  I didn't find any, of course, but I did notice
   your timing results for the Mayer FHT code.  I was rather surprised at
   how well it did, because on my system, it's not a high performer.  Both
   my own FFT and FHT code are faster.
   
   So, I thought you might be interested in looking at my code.
   
   The code is original to me.  Although I have copyrighted it, you can
   freely use it for any reason provided my name remains.  The only reason
   it's not PD is because I put quite a bit of work into it and I just like
   the idea of my name being in code that somebody else is using.
   
   [...]
*/

/*
** This Fast Hartley Transform code is Copyrighted by its author, Carey
** Bloodworth, on August 13, 1998.  You may freely use it for any
** purpose provided this notice remains.
*/

/*
** Additional legal note:
**
** It appeares that Bracewell has somehow managed to claim patent on
** the entire FHT concept.  That's quite a stretch considering the Hartley
** transform was developed 50 years ago, and the tehcniques to split that
** type of transform (ie: the Fourier Transform) into the 'Fast' variety
** were developed 30 years ago (Cooley, etc.) and possibly even far before
** that (Gauss).
**
** For the record, I developed my FHT code directly from a recursive Fast
** Fourier program and the Hartley mathematical formula.  It was almost a
** straight, simple conversion to the FHT program simply by replacing the
** complex function with Harltey's 'cas' function.
**
** However, as I am in absolutely no financial position to legally challenge
** Bracewell (note no 'Mr.' in front of his name as a sign of respect....!
** I have a strong opinion of people who file software patents), this code
** may or may not fall under his legal patents.  Use at your own risk because
** I can't afford to even check, much less fight.
**
** Additional note:  As FFT's, even with the real<->complex wrappers are very
** often faster than FHTs, you shouldn't have a compeling need to even use
** a FHT, mine or Bracewell's.  Therefore this legal warning is fairly moot
** since you probably wont be using any FHT code anyway.
**
** Carey Bloodworth, August 21, 1998
*/

#include <math.h>

/*
** In the FwdRealFFT() and RevRealFFT(), you have a choice of doing either
** the DiT style or DiF style FHT.
**
** Output will be in the standard Mayer FHT/FFT format, which is _not_
** normal 'Numerical Recipe FFT output' format.
**
** If you are doing convolutions, see note at end of file.
*/


/*
** Precalculated trig in long bench_real (64 bits of mantissa).
** Perhaps not accurate to the last bit, but it's good enough.
*/
static bench_real SineTable[]={
/* sin(MY_PI     ) */      0.000000000000000000000000000000000000000e+00,
/* sin(MY_PI/1   ) */      0.000000000000000000000000000000000000000e+00,
/* sin(MY_PI/2   ) */      1.000000000000000000000000000000000000000e+00,
/* sin(MY_PI/4   ) */      7.071067811865475027412186737052479656995e-01,
/* sin(MY_PI/8   ) */      3.826834323650897575798401906155277174548e-01,
/* sin(MY_PI/16  ) */      1.950903220161282603309360617060974618653e-01,
/* sin(MY_PI/32  ) */      9.801714032956059818174621156572356994729e-02,
/* sin(MY_PI/64  ) */      4.906767432741801234115375240918410781887e-02,
/* sin(MY_PI/128 ) */      2.454122852291228707409531661909340982675e-02,
/* sin(MY_PI/256 ) */      1.227153828571992560101701352781589093865e-02,
/* sin(MY_PI/512 ) */      6.135884649154475120776119911880641666357e-03,
/* sin(MY_PI/1k  ) */      3.067956762965976150181468540267815114930e-03,
/* sin(MY_PI/2k  ) */      1.533980186284765552615777517431183696317e-03,
/* sin(MY_PI/4k  ) */      7.669903187427044969921158257264437452250e-04,
/* sin(MY_PI/8k  ) */      3.834951875713955740925670268026692610874e-04,
/* sin(MY_PI/16k ) */      1.917475973107032999678822626776764082024e-04,
/* sin(MY_PI/32k ) */      9.587379909597734213219655252657958044438e-05,
/* sin(MY_PI/64k ) */      4.793689960306688268003999509048185245774e-05,
/* sin(MY_PI/128k) */      2.396844980841821779338901565736819065933e-05,
/* sin(MY_PI/256k) */      1.198422490506970595472088780830688392598e-05,
/* sin(MY_PI/512k) */      5.992112452642427609071640315363538320526e-06,
/* sin(MY_PI/1m  ) */      2.996056226334660633689455089267994480906e-06,
/* sin(MY_PI/2m  ) */      1.498028113169011170304617541759739651752e-06,
/* sin(MY_PI/4m  ) */      7.490140565847156918673210856951527603087e-07,
/* sin(MY_PI/8m  ) */      3.745070282923841092784580930619142691285e-07,
/* sin(MY_PI/16m ) */      1.872535141461953375578708413939921229030e-07,
/* sin(MY_PI/32m ) */      9.362675707309807915379451515036635100842e-08,
/* sin(MY_PI/64m ) */      4.681337853654909086833363351942693952878e-08,
/* sin(MY_PI/128m) */      2.340668926827455184830686918395770135248e-08
};


#define MY_PI      3.1415926535897932384626434
#define MY_SQRT_2  0.7071067811865475244008443621
#define MY_SQRT2  1.41421356237309504880

static int Log2(int Num)
{int x=-1;
if (Num==0) return 0;
while (Num) {x++;Num/=2;}
return x;
}

#define NextTrigPow(Pr,Pi,Nr,Ni)    \
 {/*long*/ bench_real temp;             \
   temp = Pr;                       \
   Pr = Pr * Nr - Pi   * Ni + Pr;   \
   Pi = Pi * Nr + temp * Ni + Pi;   \
 }

/* Macro for the DiF FHT */
#define FHT_F2Butterfly(N1,N2,C,S)           \
 {bench_real D1,D2;                              \
  int i1=N1, i2=N2;                          \
  D1=Left[i1];D2=Left[i2];                   \
  {bench_real temp;                              \
   Left[i1] =D1+(temp=Right[i1]);D1=D1-temp; \
   Left[i2] =D2+(temp=Right[i2]);D2=D2-temp; \
  }                                          \
  Right[i1]=D1*(C)+D2*(S);                   \
  Right[i2]=D1*(S)-D2*(C);                   \
 }

void RFHT_F(bench_real *Data,int Len)
/*
** Recursive Decimation in Frequency style Fast Hartley Transform
*/
{int x,Len2,Len4;
 bench_real Sin0,Cos0;
 bench_real Sin,Cos;
 bench_real *Left,*Right;

Len/=2;Left=&Data[0];Right=&Data[Len];
if (Len==2)
  {bench_real d0=Data[0]; bench_real d1=Data[1];
   bench_real d2=Data[2]; bench_real d3=Data[3];
   {bench_real d02=d0+d2; bench_real d13=d1+d3;
    Data[0]=d02+d13; Data[1]=d02-d13;
   }
   {bench_real d02=d0-d2; bench_real d13=d1-d3;
    Data[2]=d02+d13; Data[3]=d02-d13;
   }
   return;
  }

{bench_real t1,t2;
 t1=Left[0];t2=Right[0];
 Left[0]=t1+t2;Right[0]=t1-t2;
 t1=Left[Len/2];t2=Right[Len/2];
 Left[Len/2]=t1+t2;Right[Len/2]=t1-t2;
}

x=Log2(Len)+1;
Sin0=SineTable[x];
Cos0=SineTable[x+1];Cos0=-2.0*Cos0*Cos0;
Sin=Sin0;Cos=1.0+Cos0;

Len2=Len/2;
Len4=Len/4;
for (x=1;x<Len4;x++)
  {
   FHT_F2Butterfly(x,Len-x,Cos,Sin);
   FHT_F2Butterfly(Len2-x,Len2+x,Sin,Cos);
   NextTrigPow(Cos,Sin,Cos0,Sin0);
  }

/* Now do the two Len/4 points the loop missed */
if (Len4)
  {bench_real sq=MY_SQRT_2; /* variable allows optimizations */
   FHT_F2Butterfly(Len4,Len-Len4,sq,sq);
  }

if (Len>=2) RFHT_F(Left, Len);
if (Len>=2) RFHT_F(Right,Len);
}


/* Macro for the DiT FHT */
#define FHT_T2Butterfly(N1,N2,C,S) \
 {bench_real Rx,Ri;                  \
  int i1=N1,i2=N2;               \
  Rx=Right[i1];Ri=Right[i2];     \
  {bench_real cas1,Lx;               \
   cas1=Rx*(C)+Ri*(S);           \
   Lx=Left[i1];                  \
   Left[i1]  = Lx+cas1;          \
   Right[i1] = Lx-cas1;          \
  }                              \
  {bench_real cas2,Li;               \
   cas2=Rx*(S)-Ri*(C);           \
   Li=Left[i2];                  \
   Left[i2]  = Li+cas2;          \
   Right[i2] = Li-cas2;          \
  }                              \
 }

/* Macro for the DiT FHT */
#define FHT_T1Butterfly(N1,N2,C,S)         \
 {int i1=N1,i2=N2;                         \
  bench_real cas1=Right[i1]*(C)+Right[i2]*(S); \
  bench_real temp=Left[i1];                    \
  Left[i1] = temp + cas1;                  \
  Right[i2]= temp - cas1;                  \
 }

void RFHT_T(bench_real *Data,int Len)
/*
** recursive Decimation in Time style Fast Hartley Transform
*/
{int x,Len2,Len4;
 bench_real Sin0,Cos0;
 bench_real Sin,Cos;
 bench_real *Left,*Right;

Len/=2;Right=&Data[Len];Left=&Data[0];
if (Len==4)
  {bench_real d45,d67,sd0123,dd0123;
   {bench_real ss0123,ds0123,ss4567,ds4567;
    {bench_real s01,s23,d01,d23;
     d01 = Data[0] - Data[1];
     s01 = Data[0] + Data[1];
     d23 = Data[2] - Data[3];
     s23 = Data[2] + Data[3];
     ds0123 = (s01 - s23);
     ss0123 = (s01 + s23);
     dd0123 = (d01 - d23);
     sd0123 = (d01 + d23);
    }
    {bench_real s45,s67;
     s45 = Data[4] + Data[5];
     s67 = Data[6] + Data[7];
     d45 = Data[4] - Data[5];
     d67 = Data[6] - Data[7];
     ds4567 = (s45 - s67);
     ss4567 = (s45 + s67);
    }
    Data[4] = ss0123 - ss4567;
    Data[0] = ss0123 + ss4567;
    Data[6] = ds0123 - ds4567;
    Data[2] = ds0123 + ds4567;
   }
   d45 *= MY_SQRT2;
   d67 *= MY_SQRT2;
   Data[5] = sd0123 - d45;
   Data[1] = sd0123 + d45;
   Data[7] = dd0123 - d67;
   Data[3] = dd0123 + d67;
   return;
  }

RFHT_T(&Left[0], Len);
RFHT_T(&Right[0],Len);

/* Do the special x=0 loop below. */
FHT_T1Butterfly(0,0,1.0,0.0);

x=Log2(Len)+1;
Sin0=SineTable[x];
Cos0=SineTable[x+1];Cos0=-2.0*Cos0*Cos0;
Sin=Sin0;Cos=1.0+Cos0;

Len2=Len/2;
Len4=Len/4;
for (x=1;x<Len4;x++)
  {
   FHT_T2Butterfly(x,Len-x,Cos,Sin);
   FHT_T2Butterfly(Len2-x,Len2+x,Sin,Cos);
   NextTrigPow(Cos,Sin,Cos0,Sin0);
  }

/* Now do the two Len/4 points the loop missed */
if (Len4)
  {bench_real sq=MY_SQRT_2; /* variable allows optimizations */
   FHT_T2Butterfly(Len4,Len-Len4,sq,sq);
  }
/* Now do the Len/2 point the loop couldn't do. */
if (Len2)
  FHT_T1Butterfly(Len2,Len2,0.0,1.0);
}

void
FHT_Scramble(bench_real *Data,int Len)
{int k1,k2,k;
 for (k1=1,k2=0;k1<Len;k1++)
    {
     for (k=Len>>1; (!((k2^=k)&k)); k>>=1) ;
     if (k1>k2)
       {bench_real Temp;
         Temp=Data[k1];
         Data[k1]=Data[k2];
         Data[k2]=Temp;
       }
    }
}

void
FwdRealFFT(bench_real *Data,int Len)
{bench_real a,b,c,d;
 int i,j,k;
 FHT_Scramble(Data,Len);RFHT_T(Data,Len);
/* RFHT_F(Data,Len);FHT_Scramble(Data,Len); */
 for (i=1,j=Len-1,k=Len/2;i<k;i++,j--)
   {
    a = Data[i];
    b = Data[j];
    Data[j] = (a-b)*0.5;
    Data[i] = (a+b)*0.5;
   }
}

void
RevRealFFT(bench_real *Data,int Len)
{bench_real a,b,c,d;
 int i,j,k;
 for (i=1,j=Len-1,k=Len/2;i<k;i++,j--)
   {
    a = Data[i];
    b = Data[j];
    Data[j] = (a-b);
    Data[i] = (a+b);
   }
 FHT_Scramble(Data,Len);RFHT_T(Data,Len);
/* RFHT_F(Data,Len);FHT_Scramble(Data,Len); */
}


void
InitFHT(unsigned int Len)
/*
** If you need to use your own trig tables, you can calculate
** them here.
*/
{int x;unsigned int P;
#if (LDBL_DIG > 18) || (LDBL_MANT_DIG > 64)
SineTable[0]=-9.0;
x=1;P=1;
while (P<=Len*4)
  {
   SineTable[x]=sin(MY_PI/P);
   P*=2;
   x++;
  }
#endif
}


/*
** Note: if your only goal is to do a FHT convolution, you can
** save the scramble and the FHT<->FFT wrapper step and do the
** convolution directly with the scrambled data.  This also drastically
** improve virtual memory performance.
**
** This is taken directly from a FHT based multiplication routine
** in a pi program.  Len2 is the total length of the FFT/FHT that
** was performed.  Len was the length of the not-yet zero padded
** numbers that I was multiplying.


  RFHT_F(FFTNum1,Len2);

{int Step=2;
 int Step2=Step*2;
 dn=0.5 / Len2;

 while (Step < Len2)
  { 
  for (x=Step,y=Step2-1;x<Step2;x+=2,y-=2)
    {bench_real h1p,h1m,h2p,h2m;
     bench_real s1,d1;
     h1p=FFTNum1[x];
     h1m=FFTNum1[y];
     s1=h1p+h1m;
     d1=h1p-h1m;
     h2p=FFTNum2[x];
     h2m=FFTNum2[y];
     FFTNum1[x]=(h2p*s1+h2m*d1)*dn;
     FFTNum1[y]=(h2m*s1-h2p*d1)*dn;
    }
  Step*=2;
  Step2*=2;
  }
  FFTNum1[0]   = FFTNum1[0]   * 2.0 * dn * FFTNum1[0];
  FFTNum1[1]   = FFTNum1[1]   * 2.0 * dn * FFTNum1[1];
}

  RFHT_T(FFTNum1,Len2);

*/

