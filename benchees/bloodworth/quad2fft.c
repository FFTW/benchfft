/* Received in personal communication with the author (see fht.c). */

/* Modified 8/18/98 by Steven G. Johnson (stevenj@alum.mit.edu) to
   work with benchFFT (http://www.fftw.org/benchfft). */

#include "bench-user.h"
   
/*
** These "Quad-2" FFTs are Copyrighted by their author, Carey Bloodworth,
** on August 13, 1998.  They may be freely used by anyone for any reason
** provided this notice remains.
*/

#include <math.h>

/*
** A 'Quad-2' FFT is my own name for a Radix-2 FFT that does all four
** trig quadrants at once.  This cuts the trig calculations to only one
** fourth it would otherwise be, plus it reduces trig error since the
** trig recurrance is used less.  And, unlike a Radix-4 or Radix-8, no
** additional registers are required, which makes it nice for the x86.
**
** On systems with very slow memory compared to the CPU, you may get
** better performance by only doing a Radix-2 recursive FFT, instead
** of a Quad-2 recursive FFT.
** **** See below for the two #define's you can do for this ****
**
** Additionally, the FFTs (DiT & DiF) are done in recursive / iterative
** pairs.  By themselves, the recursive formula is slightly slower than
** an iterative FFT, but it has a much better L2 cache locality.  Put the
** two together, and you can see a measurable speed increase.  Excluding
** the scrambling, the recursive/iterative combination also works fairly
** well on virtual memory systems.
**
** The FFT formulas themselves do not do the scrambling or the normalization.
** This is taken care of by the four wrapper functions the user calls:

void FwdRealFFT(bench_real *RealData,int Len)
void RevRealFFT(bench_real *RealData,int Len)
void FwdFFT(Cmplx *Data,int Len)
void RevFFT(Cmplx *Data,int Len)

** If desired (or required), you can call InitFFT(int Len) to initialize
** the trig values.
**
** The four functions use the DiT FFT by defualt, since that seems to run
** slightly faster than a DiF.  Additionally, I have provided hard-wired
** forward and reverse DiT FFTs, since the directional hardwiring provides
** about 5% improvement.
*/

/*
** Whether to do a 'quad' style FFT.  IE: whether to do four radix-2
** butterflies at once.
**
** RQUAD is for the recursive FFT and IQUAD is for the iterative FFT
**
** For 'fast' memory systems, it will likely be faster to do the 'Quad-2'
** style FFT.  For slower memory systems, it'll probably be better to do
** a normal Radix-2 style recursive FFT by not defining RQUAD.  It'll
** almost certainly be faster to leave the IQUAD defined.
*/
#define RQUAD 1
#define IQUAD 1

/*
** The size of the L2 cache, in bytes.  This is used to determine
** the recursive / iterative cross-over.  128k works well on 486 systems,
** even when it has 256k L2 cache, because the 486 TLB only does 128k.
*/
#define L2_CACHE 131072

/*
** Precalculated trig in long bench_real (64 bits of mantissa).
** Perhaps not accurate to the last bit, but it's good enough.
*/
 /* long */  static bench_real SineTable[]={
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

typedef struct {bench_real r,i;} Cmplx;


#define MY_PI      3.1415926535897932384626434
#define MY_SQRT2  0.7071067811865475244008443621

#define POS (+1.0)
#define NEG (-1.0)

static int Log2(int Num)
{int x=-1;
if (Num==0) return 0;
while (Num) {x++;Num/=2;}
return x;
}

#define NextTrigPow(Pr,Pi,Nr,Ni)    \
 { /* long */  bench_real temp;             \
   temp = Pr;                       \
   Pr = Pr * Nr - Pi   * Ni + Pr;   \
   Pi = Pi * Nr + temp * Ni + Pi;   \
 }

/*
** Decimation in Frequency transforms
*/
/*
** On the register poor x86, it seems to be faster to do four
** seperate FFT2_?ButterflyLoop()'s than it is to bring the
** four butterflies into one loop.  Your milage may vary.
*/
#define FFT2_FButterflyLoop(Ndx,Sr,Pr,Si,Pi)             \
   {Cmplx *Left,*Right;                                  \
    Left=&Data[(Ndx)];                                   \
    Right=Left+Step2;                                    \
    while (Left < LastPoint)                             \
      {                                                  \
       {bench_real temp_r,temp_i,t;                          \
        Left->r  = (temp_r = Left->r) + (t = Right->r);  \
        temp_r  -= t;                                    \
        Left->i  = (temp_i = Left->i) + (t = Right->i);  \
        temp_i  -= t;                                    \
        Right->r = temp_r*(Pr)*(Sr) - temp_i*(Pi)*(Si);  \
        Right->i = temp_i*(Pr)*(Sr) + temp_r*(Pi)*(Si);  \
       }                                                 \
       Left += Step; Right += Step;                      \
      }                                                  \
   }

#define RFFT2_FButterfly(Ndx,Sr,Pr,Si,Pi)                   \
    {bench_real temp_r,temp_i,t;int z=(Ndx);                    \
     Left[z].r  = (temp_r = Left[z].r) + (t = Right[z].r);  \
     temp_r  -= t;                                          \
     Left[z].i  = (temp_i = Left[z].i) + (t = Right[z].i);  \
     temp_i  -= t;                                          \
     Right[z].r = temp_r*(Pr)*(Sr) - temp_i*(Pi)*(Si);      \
     Right[z].i = temp_i*(Pr)*(Sr) + temp_r*(Pi)*(Si);      \
    }

static void
FFT2_F(Cmplx *Data, int Len, int Dir)
/*
** Iterative Radix-2 Decimation in Frequency FFT.
*/
{
  int Step, Step2, Step4, Step8;
 /* long */ bench_real Nth_r,Nth_i;
 /* long */ bench_real Pow_r=1.0,Pow_i=0.0;
  int b;
  int L2Len;
  Cmplx *LastPoint=&Data[Len];

  L2Len=Log2(Len);

  Step = Len;
  while (Step > 1)
    {
      Step2 = Step/2;
      Step4 = Step/4;
      Step8 = Step/8;

      Nth_i = SineTable[L2Len];
      Nth_r = SineTable[L2Len+1];Nth_r = -2.0 * Nth_r * Nth_r;
      L2Len--;
      Pow_r = 1.0;Pow_i=0.0;
      FFT2_FButterflyLoop(0,POS,1.0,Dir,0.0);
#ifdef IQUAD
      for (b = 1; b < Step8; b++)
#else
      for (b = 1; b < Step2; b++)
#endif
        {
          NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
          FFT2_FButterflyLoop(b,      POS,Pow_r,Dir,Pow_i);
#ifdef IQUAD
          FFT2_FButterflyLoop(Step2-b,NEG,Pow_r,Dir,Pow_i);
          FFT2_FButterflyLoop(Step4-b,POS,Pow_i,Dir,Pow_r);
          FFT2_FButterflyLoop(Step4+b,NEG,Pow_i,Dir,Pow_r);
#endif
        }
#ifdef IQUAD
      if (Step8)
        { /* long */  bench_real sq=MY_SQRT2;
         FFT2_FButterflyLoop(Step8,      POS,sq,Dir,sq);
         FFT2_FButterflyLoop(Step2-Step8,NEG,sq,Dir,sq);
        }
      if (Step4)
         FFT2_FButterflyLoop(Step4,POS,0.0,Dir,1.0);
#endif
      Step/=2;
    }
}

static void
RFFT2_F(Cmplx *Data,int Len,int Dir)
/*
** Recursive Radix-2 Decimation in Frequency FFT.
*/
{int x;
 Cmplx *Left,*Right;
  /* long */ bench_real Pow_r,Pow_i;
  /* long */ bench_real Nth_r,Nth_i;
 int Len2,Len4,Len8;

if (Len<=(L2_CACHE/sizeof(Cmplx))) {FFT2_F(Data,Len,Dir);return;}

Len2=Len/2;Len4=Len/4;Len8=Len/8;
Left=&Data[0];Right=&Data[Len2];
x=Log2(Len);
Nth_i=SineTable[x];
Nth_r=SineTable[x+1];Nth_r=-2.0*Nth_r*Nth_r;
Pow_r = 1.0;Pow_i = 0.0;

RFFT2_FButterfly(0,POS,1.0,Dir,0.0);
#ifdef RQUAD
for (x=1;x<Len8;x++)
#else
for (x=1;x<Len2;x++)
#endif
  {
    NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
    RFFT2_FButterfly(x,     POS,Pow_r,Dir,Pow_i);
#ifdef RQUAD
    RFFT2_FButterfly(Len2-x,NEG,Pow_r,Dir,Pow_i);
    RFFT2_FButterfly(Len4-x,POS,Pow_i,Dir,Pow_r);
    RFFT2_FButterfly(Len4+x,NEG,Pow_i,Dir,Pow_r);
#endif
  }
#ifdef RQUAD
if (Len8)
  { /* long */ bench_real sq=MY_SQRT2;
   RFFT2_FButterfly(Len8,     POS,sq,Dir,sq);
   RFFT2_FButterfly(Len2-Len8,NEG,sq,Dir,sq);
  }

if (Len4)
   RFFT2_FButterfly(Len4,POS,0.0,Dir,1.0);
#endif

if (Len >= 4)
  {
   RFFT2_F(Data,     Len2,Dir);
   RFFT2_F(Data+Len2,Len2,Dir);
  }
}

/*
** Decimation in time transforms.
*/
/* The iterative butterfly loop. */
/*
** On the register poor x86, it seems to be faster to do four
** seperate FFT2_?ButterflyLoop()'s than it is to bring the
** four butterflies into one loop.  Your milage may vary.
*/
#define FFT2_TButterflyLoop(Ndx,Sr,Pr,Si,Pi)             \
   {Cmplx *Left,*Right;                                  \
    Left=&Data[(Ndx)];                                   \
    Right=Left+Step2;                                    \
    while (Left < LastPoint)                             \
      {                                                  \
       {bench_real tr, ti, t;                                \
        tr = (ti = Right->r) * (Pr) * (Sr) -             \
             Right->i * (Pi) * (Si);                     \
        Right->r = (t = Left->r) - tr;                   \
        Left->r  = t + tr;                               \
        ti = ti * (Pi) * (Si) +                          \
             Right->i * (Pr) *(Sr);                      \
        Right->i = (t = Left->i) - ti;                   \
        Left->i  = t + ti;                               \
       }                                                 \
       Left += Step; Right += Step;                      \
      }                                                  \
   }

/* The butterfly for the recursive FFTs */
#define RFFT2_TButterfly(Ndx,Sr,Pr,Si,Pi)             \
   {int z=(Ndx);                                      \
    {bench_real tr, ti, t;                                \
     tr = (ti = Right[z].r) * (Pr) * (Sr) -           \
          Right[z].i * (Pi) * (Si);                   \
     Right[z].r = (t = Left[z].r) - tr;               \
     Left[z].r  = t + tr;                             \
     ti = ti * (Pi) * (Si) +                          \
          Right[z].i * (Pr) *(Sr);                    \
     Right[z].i = (t = Left[z].i) - ti;               \
     Left[z].i  = t + ti;                             \
    }                                                 \
   }

static void
FFT2_T(Cmplx *Data, int Len,int Dir)
/*
** Iterative Radix-2 Decimation in Time FFT.
*/
{
  int Step, Step2, Step4, Step8;
  /* long */bench_real Nth_r,Nth_i;
  /* long */bench_real Pow_r,Pow_i;
  int b;
  int TrigIndex=1;
  Cmplx *LastPoint=&Data[Len];

  Step = 1;
  while (Step < Len)
    {
      Step *= 2;
      Step2 = Step/2;
      Step4 = Step/4;
      Step8 = Step/8;

      FFT2_TButterflyLoop(0,POS,1.0,Dir,0.0);

      Pow_r = 1.0;Pow_i = 0.0;
      Nth_i = SineTable[TrigIndex++];
      Nth_r = -2.0*SineTable[TrigIndex]*SineTable[TrigIndex];
#ifdef IQUAD
      for (b = 1; b < Step8; b++)
#else
      for (b = 1; b < Step2; b++)
#endif
        {Cmplx *Left,*Right;
          NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
          FFT2_TButterflyLoop(b,      POS,Pow_r,Dir,Pow_i);
#ifdef IQUAD
          FFT2_TButterflyLoop(Step2-b,NEG,Pow_r,Dir,Pow_i);
          FFT2_TButterflyLoop(Step4-b,POS,Pow_i,Dir,Pow_r);
          FFT2_TButterflyLoop(Step4+b,NEG,Pow_i,Dir,Pow_r);
#endif
        }
#ifdef IQUAD
      if (Step8)
        { /* long */  bench_real sq=MY_SQRT2;
         FFT2_TButterflyLoop(Step8,      POS,sq,Dir,sq);
         FFT2_TButterflyLoop(Step2-Step8,NEG,sq,Dir,sq);
        }
      if (Step4)
         FFT2_TButterflyLoop(Step4,POS,0.0,Dir,1.0);
#endif
    }
}

static void
RFFT2_T(Cmplx *Data,int Len,int Dir)
/*
** Recursive Radix-2 Decimation in Time FFT.
*/
{int x,Len2,Len4,Len8;
 Cmplx *Left,*Right;
 /* long */ bench_real Pow_r=1.0,Pow_i=0.0;
 /* long */ bench_real Nth_r,Nth_i;

Len2 = Len/2;
Len4 = Len/4;
Len8 = Len/8;

if (Len<=(L2_CACHE/sizeof(Cmplx))) {FFT2_T(Data,Len,Dir);return;}
if (Len2 >= 2) {RFFT2_T(Data,Len2,Dir);RFFT2_T(Data+Len2,Len2,Dir);}

x=Log2(Len);
Nth_r=-2.0*SineTable[x+1]*SineTable[x+1];
Nth_i=SineTable[x];

Left=&Data[0];
Right=&Data[Len2];
RFFT2_TButterfly(0,POS,1.0,Dir,0.0);
#ifdef RQUAD
for (x=1;x<Len8;x++)
#else
for (x=1;x<Len2;x++)
#endif
  {
    NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
    RFFT2_TButterfly(x,     POS,Pow_r,Dir,Pow_i);
#ifdef RQUAD
    RFFT2_TButterfly(Len2-x,NEG,Pow_r,Dir,Pow_i);
    RFFT2_TButterfly(Len4-x,POS,Pow_i,Dir,Pow_r);
    RFFT2_TButterfly(Len4+x,NEG,Pow_i,Dir,Pow_r);
#endif
  }
#ifdef RQUAD
if (Len8)
   { /* long */ bench_real sq=MY_SQRT2;
    RFFT2_TButterfly(Len8,     POS,sq,Dir,sq);
    RFFT2_TButterfly(Len2-Len8,NEG,sq,Dir,sq);
   }
if (Len4)
    RFFT2_TButterfly(Len4,POS,0,Dir,1.0);
#endif
}


/*
** The normal DiT FFTs, except I've hard-wired them for a specific
** direction.  This is about 5% faster.
*/
static void
RevFFT2_T(Cmplx *Data, int Len)
/*
** Iterative Radix-2 Decimation in Time FFT.
*/
{
  int Step, Step2, Step4, Step8;
  /* long */bench_real Nth_r,Nth_i;
  /* long */bench_real Pow_r,Pow_i;
  int b;
  int TrigIndex=1;
  Cmplx *LastPoint=&Data[Len];

  Step = 1;
  while (Step < Len)
    {
      Step *= 2;
      Step2 = Step/2;
      Step4 = Step/4;
      Step8 = Step/8;

      FFT2_TButterflyLoop(0,POS,1.0,NEG,0.0);

      Pow_r = 1.0;Pow_i = 0.0;
      Nth_i = SineTable[TrigIndex++];
      Nth_r = -2.0*SineTable[TrigIndex]*SineTable[TrigIndex];
#ifdef IQUAD
      for (b = 1; b < Step8; b++)
#else
      for (b = 1; b < Step2; b++)
#endif
        {Cmplx *Left,*Right;
          NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
          FFT2_TButterflyLoop(b,      POS,Pow_r,NEG,Pow_i);
#ifdef IQUAD
          FFT2_TButterflyLoop(Step2-b,NEG,Pow_r,NEG,Pow_i);
          FFT2_TButterflyLoop(Step4-b,POS,Pow_i,NEG,Pow_r);
          FFT2_TButterflyLoop(Step4+b,NEG,Pow_i,NEG,Pow_r);
#endif
        }
#ifdef IQUAD
      if (Step8)
        { /* long */  bench_real sq=MY_SQRT2;
         FFT2_TButterflyLoop(Step8,      POS,sq,NEG,sq);
         FFT2_TButterflyLoop(Step2-Step8,NEG,sq,NEG,sq);
        }
      if (Step4)
         FFT2_TButterflyLoop(Step4,POS,0.0,NEG,1.0);
#endif
    }
}

static void
RevRFFT2_T(Cmplx *Data,int Len)
/*
** Recursive Radix-2 Decimation in Time FFT.
*/
{int x,Len2,Len4,Len8;
 Cmplx *Left,*Right;
 /* long */ bench_real Pow_r=1.0,Pow_i=0.0;
 /* long */ bench_real Nth_r,Nth_i;

Len2 = Len/2;
Len4 = Len/4;
Len8 = Len/8;

if (Len<=(L2_CACHE/sizeof(Cmplx))) {RevFFT2_T(Data,Len);return;}
if (Len2 >= 2) {RevRFFT2_T(Data,Len2);RevRFFT2_T(Data+Len2,Len2);}

x=Log2(Len);
Nth_r=-2.0*SineTable[x+1]*SineTable[x+1];
Nth_i=SineTable[x];

Left=&Data[0];
Right=&Data[Len2];
RFFT2_TButterfly(0,POS,1.0,NEG,0.0);
#ifdef RQUAD
for (x=1;x<Len8;x++)
#else
for (x=1;x<Len2;x++)
#endif
  {
    NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
    RFFT2_TButterfly(x,     POS,Pow_r,NEG,Pow_i);
#ifdef RQUAD
    RFFT2_TButterfly(Len2-x,NEG,Pow_r,NEG,Pow_i);
    RFFT2_TButterfly(Len4-x,POS,Pow_i,NEG,Pow_r);
    RFFT2_TButterfly(Len4+x,NEG,Pow_i,NEG,Pow_r);
#endif
  }
#ifdef RQUAD
if (Len8)
   { /* long */ bench_real sq=MY_SQRT2;
    RFFT2_TButterfly(Len8,     POS,sq,NEG,sq);
    RFFT2_TButterfly(Len2-Len8,NEG,sq,NEG,sq);
   }
if (Len4)
    RFFT2_TButterfly(Len4,POS,0,NEG,1.0);
#endif
}


static void
FwdFFT2_T(Cmplx *Data, int Len)
/*
** Iterative Radix-2 Decimation in Time FFT.
*/
{
  int Step, Step2, Step4, Step8;
  /* long */bench_real Nth_r,Nth_i;
  /* long */bench_real Pow_r,Pow_i;
  int b;
  int TrigIndex=1;
  Cmplx *LastPoint=&Data[Len];

  Step = 1;
  while (Step < Len)
    {
      Step *= 2;
      Step2 = Step/2;
      Step4 = Step/4;
      Step8 = Step/8;

      FFT2_TButterflyLoop(0,POS,1.0,POS,0.0);

      Pow_r = 1.0;Pow_i = 0.0;
      Nth_i = SineTable[TrigIndex++];
      Nth_r = -2.0*SineTable[TrigIndex]*SineTable[TrigIndex];
#ifdef IQUAD
      for (b = 1; b < Step8; b++)
#else
      for (b = 1; b < Step2; b++)
#endif
        {Cmplx *Left,*Right;
          NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
          FFT2_TButterflyLoop(b,      POS,Pow_r,POS,Pow_i);
#ifdef IQUAD
          FFT2_TButterflyLoop(Step2-b,NEG,Pow_r,POS,Pow_i);
          FFT2_TButterflyLoop(Step4-b,POS,Pow_i,POS,Pow_r);
          FFT2_TButterflyLoop(Step4+b,NEG,Pow_i,POS,Pow_r);
#endif
        }
#ifdef IQUAD
      if (Step8)
        { /* long */  bench_real sq=MY_SQRT2;
         FFT2_TButterflyLoop(Step8,      POS,sq,POS,sq);
         FFT2_TButterflyLoop(Step2-Step8,NEG,sq,POS,sq);
        }
      if (Step4)
         FFT2_TButterflyLoop(Step4,POS,0.0,POS,1.0);
#endif
    }
}

static void
FwdRFFT2_T(Cmplx *Data,int Len)
/*
** Recursive Radix-2 Decimation in Time FFT.
*/
{int x,Len2,Len4,Len8;
 Cmplx *Left,*Right;
 /* long */ bench_real Pow_r=1.0,Pow_i=0.0;
 /* long */ bench_real Nth_r,Nth_i;

Len2 = Len/2;
Len4 = Len/4;
Len8 = Len/8;

if (Len<=(L2_CACHE/sizeof(Cmplx))) {FwdFFT2_T(Data,Len);return;}
if (Len2 >= 2) {FwdRFFT2_T(Data,Len2);FwdRFFT2_T(Data+Len2,Len2);}

x=Log2(Len);
Nth_r=-2.0*SineTable[x+1]*SineTable[x+1];
Nth_i=SineTable[x];

Left=&Data[0];
Right=&Data[Len2];
RFFT2_TButterfly(0,POS,1.0,POS,0.0);
#ifdef RQUAD
for (x=1;x<Len8;x++)
#else
for (x=1;x<Len2;x++)
#endif
  {
    NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
    RFFT2_TButterfly(x,     POS,Pow_r,POS,Pow_i);
#ifdef RQUAD
    RFFT2_TButterfly(Len2-x,NEG,Pow_r,POS,Pow_i);
    RFFT2_TButterfly(Len4-x,POS,Pow_i,POS,Pow_r);
    RFFT2_TButterfly(Len4+x,NEG,Pow_i,POS,Pow_r);
#endif
  }
#ifdef RQUAD
if (Len8)
   { /* long */ bench_real sq=MY_SQRT2;
    RFFT2_TButterfly(Len8,     POS,sq,POS,sq);
    RFFT2_TButterfly(Len2-Len8,NEG,sq,POS,sq);
   }
if (Len4)
    RFFT2_TButterfly(Len4,POS,0,POS,1.0);
#endif
}

static void
FFTReOrder(Cmplx *Data, int Len)
{int Index,xednI,k;

xednI = 0;
for (Index = 0;Index < Len;Index++)
  {
   if (xednI > Index)
     {bench_real r,i;
      r=Data[xednI].r;
      Data[xednI].r = Data[Index].r;
      Data[Index].r = r;
      i=Data[xednI].i;
      Data[xednI].i = Data[Index].i;
      Data[Index].i = i;
     }
   k=Len/2;
   while ((k <= xednI) && (k >=1)) {xednI-=k;k/=2;}
   xednI+=k;
  }
}

static void
RealFFTWrapper(bench_real *RealData,int Len, int Dir)
{bench_real Pow_r, Pow_i, Nth_r, Nth_i, temp;
 int i, j;
 Cmplx *Data=(Cmplx*)RealData;
 int Half;
 bench_real NegHalf, PosHalf;

Len /= 2;       /* Len/2 Cmplx data */
Half = Len/2;
NegHalf = -0.5;
PosHalf = 0.5;
if (Dir > 0)
  {
    NegHalf = 0.5;
    PosHalf = -0.5;
    FFTReOrder(Data,Len);FwdRFFT2_T(Data,Len);
/*    FFTReOrder(Data,Len);RFFT2_T(Data,Len,Dir);
      RFFT2_F(Data,Len,Dir);FFTReOrder(Data,Len); */
  }

i=Log2(Len)+1;
Nth_r = -2.0*SineTable[i+1]*SineTable[i+1];
Nth_i = Dir*SineTable[i];
Pow_r = 1 + Nth_r;
Pow_i = Nth_i;
for (i = 1, j = Len - i; i < Half; i++, j--)
  {bench_real p1r,p1i,p2r,p2i;
   /* Seperate the two points from the jumbled points */
   p1r =     0.5 * (Data[i].r + Data[j].r);
   p2i = PosHalf * (Data[i].r - Data[j].r);
   p2r = NegHalf * (Data[i].i + Data[j].i);
   p1i =     0.5 * (Data[i].i - Data[j].i);
   /* Almost a standard Decimation in Time butterfly... */
   {bench_real tmp;
    tmp = p2r;
    p2r = p2r * Pow_r - p2i * Pow_i;
    p2i = p2i * Pow_r + tmp * Pow_i;
   }
   Data[i].r = p1r + p2r;
   Data[i].i = p1i + p2i;
   Data[j].r = p1r - p2r;
   Data[j].i = -(p1i - p2i); /* ... except this is negated */
   NextTrigPow(Pow_r,Pow_i,Nth_r,Nth_i);
  }

temp = Data[0].r;
Data[0].r = temp + Data[0].i;
Data[0].i = temp - Data[0].i;

if (Dir < 0)
  {
    Data[0].r /= 2.0;
    Data[0].i /= 2.0;
    FFTReOrder(Data,Len);RevRFFT2_T(Data,Len);
/*    FFTReOrder(Data,Len);RFFT2_T(Data,Len,Dir);
      RFFT2_F(Data,Len,Dir);FFTReOrder(Data,Len); */
  }
}


/*
** Perform a real value forward FFT
*/
void
Bloodworth_Q2_FwdRealFFT(bench_real *RealData,int Len)
{
RealFFTWrapper(RealData,Len,1);
}

/*
** Perform a real value reverse FFT
*/
void
Bloodworth_Q2_RevRealFFT(bench_real *RealData,int Len)
{int x;bench_real Inv=2.0/Len;
RealFFTWrapper(RealData,Len,-1);
for (x=0;x<Len;x++) {RealData[x] *= Inv;}
}

/*
** Perform a complex value forward FFT
*/
void
Bloodworth_Q2_FwdFFT(Cmplx *Data,int Len)
{
FFTReOrder(Data,Len);FwdRFFT2_T(Data,Len);
/* FFTReOrder(Data,Len);RFFT2_T(Data,Len,1);
   RFFT2_F(Data,Len,1);FFTReOrder(Data,Len); */
}

/*
** Perform a complex value reverse FFT
*/
void
Bloodworth_Q2_RevFFT(Cmplx *Data,int Len)
{int x;bench_real Inv=1.0/Len;
FFTReOrder(Data,Len);RevRFFT2_T(Data,Len);
/* FFTReOrder(Data,Len);RFFT2_T(Data,Len,-1);
   RFFT2_F(Data,Len,-1);FFTReOrder(Data,Len); */
for (x=0;x<Len;x++) {Data[x].r *= Inv;Data[x].i *= Inv;}
}


void
Bloodworth_Q2_InitFFT(unsigned int Len)
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

