CC=================================================================CC
CC                                                                 CC
CC  Subroutine CFFTSR(X,Y,M):                                      CC
CC      An in-place, split-radix complex FFT program               CC
CC      Decimation-in-frequency, cos/sin in second loop            CC
CC      and is computed recursively                                CC
CC      The program is based on Tran ASSP Feb 1986 pp152-156       CC
CC                                                                 CC
CC  Input/output                                                   CC
CC      X    Array of real part of input/output (length >= N)      CC
CC      Y    Array of imaginary part of input/output (length >= N) CC
CC      M    Transform length is N=2**M                            CC
CC                                                                 CC
CC  Calls:                                                         CC
CC      CSTAGE,CBITREV                                             CC
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
CC      This program may be used and distributed freely as         CC
CC      as long as this header is included                         CC
CC                                                                 CC
CC=================================================================CC
      SUBROUTINE CFFTSR(X,Y,M)
      REAL X(1), Y(1)
      N = 2**M
C-----L shaped butterflies------------------------------------------C
      N2 = 2*N
      DO  10 K = 1, M-1
          N2 = N2/2
          N4 = N2/4
          CALL CSTAGE(N,N2,N4,X(1),X(N4+1),X(2*N4+1),X(3*N4+1),
     $         Y(1),Y(N4+1),Y(2*N4+1),Y(3*N4+1))
 10   CONTINUE
C-----Length two butterflies----------------------------------------C
      IS = 1
      ID = 4
 20       DO  30 I1 = IS,N,ID
              T1      = X(I1)
              X(I1)   = T1 + X(I1+1)
              X(I1+1) = T1 - X(I1+1)
              T1      = Y(I1)
              Y(I1)   = T1 + Y(I1+1)
              Y(I1+1) = T1 - Y(I1+1)
 30       CONTINUE
          IS = 2*ID - 1
          ID = 4*ID
      IF (IS.LT.N) GOTO 20
C-------Digit reverse counter---------------------------------------C
      CALL CBITREV(X,Y,M)
      RETURN
      END
C===================================================================C
C  Subroutine CSTAGE - the work-horse of the CFFTSR                 C
C       Computes a stage of a complex split-radix length N          C
C       transform.                                                  C
C  Author                                                           C
C       H.V. Sorensen,   University of Pennsylvania,  Jul. 1987     C
C===================================================================C
      SUBROUTINE CSTAGE(N,N2,N4,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
      REAL X1(1),X2(1),X3(1),X4(1),Y1(1),Y2(1),Y3(1),Y4(1)
      N8 = N4/2
C-------Zero butterfly----------------------------------------------C
      IS = 0
      ID = 2*N2
 10       DO  20  I1 = IS+1,N,ID
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
 20       CONTINUE
          IS = 2*ID - N2
          ID = 4*ID
      IF (IS .LT. N) GOTO 10
C
      IF  (N4-1) 100,100,30
C-------N/8 butterfly-----------------------------------------------C
 30   IS = 0
      ID = 2*N2
 40       DO  50 I1 = IS+1+N8,N,ID
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
 50       CONTINUE
          IS = 2*ID - N2
          ID = 4*ID
      IF (IS .LT. N-1) GOTO 40
C
      IF  (N8-1) 100,100,60
C-------General butterfly. Two at a time----------------------------C
 60   E  = 6.283185307179586/N2
      SS1 = SIN(E)
      SD1 = SS1
      SD3 = 3.*SD1-4.*SD1**3
      SS3 = SD3
      CC1 = COS(E)
      CD1 = CC1
      CD3 = 4.*CD1**3-3.*CD1
      CC3 = CD3
      DO  90  J = 2,N8
          IS = 0
          ID = 2*N2
          JN = N4 - 2*J + 2
 70           DO  80 I1=IS+J,N+J,ID
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
                  X3(I1) =  T1*CC1 - T4*SS1
                  Y3(I1) = -T4*CC1 - T1*SS1
                  X4(I1) =  T5*CC3 + T2*SS3
                  Y4(I1) =  T2*CC3 - T5*SS3
                  I2 = I1 + JN
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
                  X3(I2) =  T1*SS1 - T4*CC1
                  Y3(I2) = -T4*SS1 - T1*CC1
                  X4(I2) = -T5*SS3 - T2*CC3
                  Y4(I2) = -T2*SS3 + T5*CC3
 80           CONTINUE
              IS = 2*ID - N2
              ID = 4*ID
          IF (IS .LT. N) GOTO 70
C
          T1  = CC1*CD1 - SS1*SD1
          SS1 = CC1*SD1 + SS1*CD1
          CC1 = T1
          T3  = CC3*CD3 - SS3*SD3
          SS3 = CC3*SD3 + SS3*CD3
          CC3 = T3
 90   CONTINUE
 100  RETURN
      END
CC=================================================================CC
CC                                                                 CC
CC Subroutine CBITREV(X,Y,M):                                      CC
CC      Bitreverses the array X of length 2**M. It generates a     CC
CC      table ITAB (minimum length is SQRT(2**M) if M is even      CC
CC      or SQRT(2*2**M) if M is odd). ITAB need only be generated  CC
CC      once for a given transform length.                         CC
CC      The program uses the technique described by D. Evans       CC
CC      in Tran. ASSP Aug. 1987 pp1120-1125                        CC
CC                                                                 CC
CC Author:                                                         CC
CC      H.V. Sorensen,   University of Pennsylvania,  Aug. 1987    CC
CC                       Arpa address: hvs@ee.upenn.edu            CC
CC                                                                 CC
CC      This program may be used and distributed freely as long    CC
CC      as this header is included.                                CC
CC                                                                 CC
CC=================================================================CC
      SUBROUTINE CBITREV(X,Y,M)
      DIMENSION X(1),Y(1),ITAB(256)
C-------Initialization of ITAB array--------------------------------C
      M2 = M/2
      NBIT = 2**M2
      IF (2*M2.NE.M) M2 = M2 + 1
      ITAB(1) = 0
      ITAB(2) = 1
      IMAX = 1
      DO  10 LBSS = 2, M2
          IMAX = 2 * IMAX
          DO  10 I = 1, IMAX
              ITAB(I)      = 2 * ITAB(I)
              ITAB(I+IMAX) = 1 + ITAB(I)
 10   CONTINUE
C-----The actual bitreversal----------------------------------------C
      DO  20 K = 2, NBIT
          J0 = NBIT * ITAB(K) + 1
          I = K
          J = J0
          DO  20 L = 2, ITAB(K)+1
              T1   = X(I)
              X(I) = X(J)
              X(J) = T1
              T1   = Y(I)
              Y(I) = Y(J)
              Y(J) = T1
              I = I + NBIT
              J = J0 + ITAB(L)
 20   CONTINUE
      RETURN
      END
