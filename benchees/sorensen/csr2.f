CC=================================================================CC
CC                                                                 CC
CC  Subroutine CTFFTSR(X,Y,M,CT1,CT3,ST1,ST3,ITAB):                CC
CC      An in-place, split-radix complex FFT program               CC
CC      Decimation-in-frequency, cos/sin in third loop             CC
CC      and is looked-up in table. Tables CT1,CT3,ST1,ST3          CC
CC      have to have length>=N/8-1. The bit reverser uses partly   CC
CC      table lookup.                                              CC
CC                                                                 CC
CC  Input/output                                                   CC
CC      X    Array of real part of input/output (length >= N)      CC
CC      Y    Array of imaginary part of input/output (length >= N) CC
CC      M    Transform length is N=2**M                            CC
CC      CT1  Array of cos() table (length >= N/8-1)                CC
CC      CT3  Array of cos() table (length >= N/8-1)                CC
CC      ST1  Array of sin() table (length >= N/8-1)                CC
CC      ST3  Array of sin() table (length >= N/8-1)                CC
CC      ITAB Array of bitreversal indices (length >= sqrt(2*N)     CC
CC                                                                 CC
CC  Calls:                                                         CC
CC      CTSTAG                                                     CC
CC      and TINIT has to be called before this program!!           CC
CC                                                                 CC
CC  Author:                                                        CC
CC      H.V. Sorensen,   University of Pennsylvania,  Dec. 1984    CC
CC                       Arpa address: hvs@ee.upenn.edu            CC
CC  Modified:                                                      CC
CC      H.V. Sorensen,   University of Pennsylvania,  Jul. 1987    CC
CC                                                                 CC
CC  Reference:                                                     CC
CC      Sorensen, Heideman, Burrus :"On computing the split-radix  CC
CC      FFT", IEEE Tran. ASSP, Vol. ASSP-34, No. 1, pp. 152-156    CC
CC      Feb. 1986                                                  CC
CC      Mitra&Kaiser: "Digital Signal Processing Handbook, Chap.   CC
CC      8, page 491-610, John Wiley&Sons, 1993                     CC
CC                                                                 CC
CC      This program may be used and distributed freely as long    CC
CC      as this header is included                                 CC
CC                                                                 CC
CC=================================================================CC
	SUBROUTINE CTFFTSR(X,Y,M,CT1,CT3,ST1,ST3,ITAB)
	REAL X(1),Y(1),CT1(1),CT3(1),ST1(1),ST3(1)
        INTEGER ITAB(1)
	N = 2**M
C-------L shaped butterflies----------------------------------------C
        ITS = 1
	N2 = 2*N
	DO  10 K = 1, M-1
	    N2 = N2/2
	    N4 = N2/4
	    CALL CTSTAG(N,N2,N4,ITS,X(1),X(N4+1),X(2*N4+1),X(3*N4+1),
     $                              Y(1),Y(N4+1),Y(2*N4+1),Y(3*N4+1),
     $                              CT1,CT3,ST1,ST3)
            ITS = 2 * ITS
 10     CONTINUE
C-------Length two butterflies--------------------------------------C
	IS = 1
	ID = 4
 20         DO  30 I1 = IS,N,ID
		T1      = X(I1)
		X(I1)   = T1 + X(I1+1)
		X(I1+1) = T1 - X(I1+1)
		T1      = Y(I1)
		Y(I1)   = T1 + Y(I1+1)
		Y(I1+1) = T1 - Y(I1+1)
 30         CONTINUE
	    IS = 2*ID - 1
	    ID = 4*ID
	IF (IS.LT.N) GOTO 20
C-------Digit reverse counter---------------------------------------C
        M2 = M/2
        NBIT = 2**M2
	DO  50 K = 2, NBIT
            J0 = NBIT * ITAB(K) + 1
            I = K
            J = J0
            DO  40 L = 2, ITAB(K)+1
                T1   = X(I)
                X(I) = X(J)
                X(J) = T1
                T1   = Y(I)
                Y(I) = Y(J)
                Y(J) = T1
                I = I + NBIT
                J = J0 + ITAB(L)
 40         CONTINUE
 50     CONTINUE
	RETURN
	END
C===================================================================C
C   Subroutine CTSTAG - the work-horse of CTFFTSR                   C
C       Computes a stage of a length N split-radix transform        C
C   Author:                                                         C
C       H.V. Sorensen,   University of Pennsylvania,  Jul. 1987     C
C===================================================================C
        SUBROUTINE CTSTAG(N,N2,N4,ITS,X1,X2,X3,X4,Y1,Y2,Y3,Y4,
     $                    CT1,CT3,ST1,ST3)
	REAL X1(1),X2(1),X3(1),X4(1),Y1(1),Y2(1),Y3(1),Y4(1)
        REAL CT1(1),CT3(1),ST1(1),ST3(1)
	N8 = N4/2
C-------Zero butterfly----------------------------------------------C
	IS = 0
	ID = 2*N2
 10         DO  20  I1 = IS+1,N,ID
		T1     = X1(I1) - X3(I1)
		X1(I1) = X1(I1) + X3(I1)
                T2     = Y2(I1) - Y4(I1)
                Y2(I1) = Y2(I1) + Y4(I1)
                X3(I1) = T1     + T2
                T2     = T1     - T2
                T1     = X2(I1) - X4(I1)
                X2(I1) = X2(I1) + X4(I1)
                X4(I1) = T2
		T2     = Y1(I1) - Y3(I1)
                Y1(I1) = Y1(I1) + Y3(I1)
		Y3(I1) = T2     - T1
		Y4(I1) = T2     + T1
 20         CONTINUE
	    IS = 2*ID - N2
	    ID = 4*ID
	IF (IS .LT. N) GOTO 10
C
	IF  (N4-1) 100,100,30
C-------N/8 butterfly-----------------------------------------------C
 30     IS = 0
	ID = 2*N2
 40         DO  50 I1 = IS+1+N8,N,ID
                T1     = X1(I1) - X3(I1)
                X1(I1) = X1(I1) + X3(I1)
                T2     = X2(I1) - X4(I1)
                X2(I1) = X2(I1) + X4(I1)
                T3     = Y1(I1) - Y3(I1)
                Y1(I1) = Y1(I1) + Y3(I1)
                T4     = Y2(I1) - Y4(I1)
                Y2(I1) = Y2(I1) + Y4(I1)
                T5     = (T4 - T1)*0.707106778
                T1     = (T4 + T1)*0.707106778
                T4     = (T3 - T2)*0.707106778
                T2     = (T3 + T2)*0.707106778
		X3(I1) = T4 + T1
		Y3(I1) = T4 - T1
		X4(I1) = T5 + T2
		Y4(I1) = T5 - T2
 50         CONTINUE
	    IS = 2*ID - N2
	    ID = 4*ID
	IF (IS .LT. N-1) GOTO 40
C
	IF  (N8-1) 100,100,60
C-------General butterfly. Two at a time----------------------------C
 60     IS = 1
        ID = N2*2
 70         DO  90  I = IS, N, ID
                IT = 0
                JN = I + N4
		DO  80 J=1, N8-1
                    IT = IT+ITS
                    I1 = I+J
		    T1     = X1(I1) - X3(I1)
		    X1(I1) = X1(I1) + X3(I1)
		    T2     = X2(I1) - X4(I1)
		    X2(I1) = X2(I1) + X4(I1)
		    T3     = Y1(I1) - Y3(I1)
		    Y1(I1) = Y1(I1) + Y3(I1)
		    T4     = Y2(I1) - Y4(I1)
		    Y2(I1) = Y2(I1) + Y4(I1)
		    T5 = T1 - T4
		    T1 = T1 + T4
		    T4 = T2 - T3
		    T2 = T2 + T3
		    X3(I1) =  T1*CT1(IT) - T4*ST1(IT)
		    Y3(I1) = -T4*CT1(IT) - T1*ST1(IT)
		    X4(I1) =  T5*CT3(IT) + T2*ST3(IT)
		    Y4(I1) =  T2*CT3(IT) - T5*ST3(IT)
		    I2 = JN - J
		    T1     = X1(I2) - X3(I2)
		    X1(I2) = X1(I2) + X3(I2)
		    T2     = X2(I2) - X4(I2)
		    X2(I2) = X2(I2) + X4(I2)
		    T3     = Y1(I2) - Y3(I2)
		    Y1(I2) = Y1(I2) + Y3(I2)
		    T4     = Y2(I2) - Y4(I2)
		    Y2(I2) = Y2(I2) + Y4(I2)
		    T5 = T1 - T4
		    T1 = T1 + T4
		    T4 = T2 - T3
		    T2 = T2 + T3
		    X3(I2) =  T1*ST1(IT) - T4*CT1(IT)
		    Y3(I2) = -T4*ST1(IT) - T1*CT1(IT)
		    X4(I2) = -T5*ST3(IT) - T2*CT3(IT)
		    Y4(I2) = -T2*ST3(IT) + T5*CT3(IT)
 80             CONTINUE
 90         CONTINUE
            IS = 2*ID - N2 +1
            ID = 4*ID
	IF (IS .LT. N) GOTO 70
 100    RETURN
	END
CC=================================================================CC
CC  Subroutine TINIT:                                              CC
CC      Initialize SIN/COS and bit reversal tables                 CC
CC  Author:                                                        CC
CC      H.V. Sorensen,   University of Pennsylvania,  Jul. 1987    CC
CC=================================================================CC
	SUBROUTINE TINIT(M,CT1,CT3,ST1,ST3,ITAB)
	REAL CT1(1),CT3(1),ST1(1),ST3(1)
        INTEGER ITAB(1)
C-------Sin/Cos table-----------------------------------------------C
        N = 2**M
	ANG = 6.283185307179586/N
	DO  10 I=1,N/8-1
	    CT1(I) = COS(ANG*I)
	    CT3(I) = COS(ANG*I*3)
	    ST1(I) = SIN(ANG*I)
	    ST3(I) = SIN(ANG*I*3)
 10     CONTINUE
C-------Bit reversal table------------------------------------------C
        M2 = M/2
        NBIT = 2**M2
        IF (2*M2.NE.M) M2 = M2 + 1
	ITAB(1) = 0
	ITAB(2) = 1
        IMAX = 1
	DO  30 LBSS = 2, M2
            IMAX = 2 * IMAX
            DO  20 I = 1, IMAX
                ITAB(I)      = 2 * ITAB(I)
                ITAB(I+IMAX) = 1 + ITAB(I)
 20         CONTINUE
 30     CONTINUE
	RETURN
	END
