C-----------------------------------------------------------------------
C This is not the original program!  
C
C This file was modified by Matteo Frigo to name the common block containing
C X and Y.  See morris.f.orig for the original program.
C
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C 
C MAIN PROGRAM: TIME-EFFICIENT RADIX-4 FAST FOURIER TRANSFORM
C AUTHOR:       L. ROBERT MORRIS
C               DEPARTMENT OF SYSTEMS AND COMPUTER ENGINEERING
C               CARLETON UNIVERSITY, OTTAWA, CANADA K1S 5B6
C
C INPUT:        THE ARRAYS "X,Y" CONTAINS THE Re, Im PARTS OF THE
C               DATA TO BE TRANSFORMED
C-----------------------------------------------------------------------
        SUBROUTINE FFT4(N,M)
        INTEGER KV1(5),KV2(5)
        COMMON /AA/X(1024),Y(1024),TR(2295)
        DATA KV1/1,4,16,64,256/
        DATA KV2/3,12,48,192,768/
        DATA C21/0.707106778/
        DATA CM/-1.414213563/
        N2=N
        KKINC = 9
C       ------------MAIN FFT LOOPS---------------------------
        DO 10 K=1, M
          N1 = N2
          N2 = N2/4
          N3 = N2 + N2
          N4 = N3 + N2
          JT = N2/2 + 1
C       ------------SPECIAL BUTTERFLY FOR W = 1 --------------
          DO 1 I=1, N, N1
             T1        = X(I)      + X(I + N3)
             R         = X(I)      - X(I + N3)
             T2        = X(I + N2) + X(I + N4)
             X(I)      = T1        + T2
             X(I + N3) = T1        - T2
C
             T1        = Y(I)      + Y(I + N3)
             S         = Y(I)      - Y(I + N3)
             T2        = Y(I + N2) + Y(I + N4)
             Y(I)      = T1        + T2
             Y(I + N3) = T1        - T2
C
             T1        = Y(I + N2) - Y(I + N4)
             T2        = X(I + N2) - X(I + N4)
C        
             X(I + N2) = R         + T1
             X(I + N4) = R         - T1
             Y(I + N2) = S         - T2
1            Y(I + N4) = S         + T2
C
          IF(K.EQ.M) GOTO 10
          KK = -8
C       ---------------MAIN BUTTERFLIES--------------------
          DO 20 J = 2, N2
             KK = KK + KKINC
             IF(J.EQ.JT) GOTO 50
C       ------------BUTTERFLIES WITH SAME W----------------
          DO 30 I = J, N, N1
             X1        = X(I)      + X(I + N3)
             X2        = X(I)      - X(I + N3)
             T         = X(I + N2) + X(I + N4)
             X(I)      = X1        + T
             X1        = X1        - T
C
             Y1        = Y(I)      + Y(I + N3)
             Y2        = Y(I)      - Y(I + N3)
             T         = Y(I + N2) + Y(I + N4)
             Y(I)      = Y1        + T
             Y1        = Y1        - T
C
             T         = (X1+Y1)*TR(KK + 3)
             X(I + N3) = TR(KK + 4)*Y1 + T
             Y(I + N3) = TR(KK + 5)*X1 + T
C
             T         = Y(I + N2) - Y(I + N4)
             X1        = X2        + T
             X2        = X2        - T
C
             T         = X(I + N2) - X(I + N4)
             Y1        = Y2        - T
             Y2        = Y2        + T
C
             T         = (X1+Y1)*TR(KK)
             X(I + N2) = TR(KK + 1)*Y1 + T
             Y(I + N2) = TR(KK + 2)*X1 + T
C
             T         = (X2+Y2)*TR(KK + 6)
             X(I + N4) = TR(KK + 7)*Y2 + T
30           Y(I + N4) = TR(KK + 8)*X2 + T
C
        GOTO 20
C       --------------SPECIAL BUTTERFLY FOR W = J ------------
50           DO 40 I = J, N, N1
             X1        = X(I)      + X(I + N3)
             X2        = X(I)      - X(I + N3)
             Y1        = Y(I)      + Y(I + N3)
             Y2        = Y(I)      - Y(I + N3)
C
             T         = X(I + N2) + X(I + N4)
             X(I)      = T         + X1
             Y(I+N3)   = T         - X1
C
             T         = Y(I + N2) + Y(I + N4)
             Y(I)      = Y1        + T
             X(I+N3)   = Y1        - T
C
             X1        = X(I + N2) - X(I + N4)
             Y1        = Y(I + N2) - Y(I + N4)
C
             T         = X2        + Y1
             X(I + N2) = (T        - X1  + Y2)*C21
             Y(I + N2) = T*CM      + X(I + N2)
C
             T         = X1        + Y2
             X(I + N4) =(T         - X2  + Y1)*C21
40           Y(I + N4) = T*CM      + X(I + N4) 
C
20      CONTINUE
        KKINC = KKINC + KKINC
10      KKINC = KKINC + KKINC
C       --------DIGIT REVERSE COUNTER--------------------
100     J = 1
        NM1 = N - 1
        DO 104 I = 1, NM1
           IF (I.GE.J) GOTO 101
           T    = X(J)
           X(J) = X(I)
           X(I) = T
           T    = Y(J)
           Y(J) = Y(I)
           Y(I) = T
101     INDEX = M
102     IF (KV2(INDEX).GE.J) GOTO 104
103        J = J - KV2(INDEX)
           INDEX = INDEX -1
        IF (KV2(INDEX).LT.J) GOTO 103
104        J = J + KV1(INDEX)
        RETURN
        END
C
        SUBROUTINE INIR4(M)
        COMMON /AA/X(1024),Y(1024),TR(2295)
C
        KK = 1
        N2 = 4**(M-1)
        P = 8.*ATAN(1.0)/FLOAT(4**M)
        DO 20 J = 2,N2
         TR(KK+3) =  COS((2*J-2)*P)
         TR(KK+4) =  SIN((2*J-2)*P)-TR(KK+3)
         TR(KK+5) = -SIN((2*J-2)*P)-TR(KK+3)          
         TR(KK)   =  COS((J-1)*P)
         TR(KK+1) =  SIN((J-1)*P)  -TR(KK)
         TR(KK+2) = -SIN((J-1)*P)  -TR(KK)          
         TR(KK+6) =  COS((3*J-3)*P)
         TR(KK+7) =  SIN((3*J-3)*P)-TR(KK+6)
         TR(KK+8) = -SIN((3*J-3)*P)-TR(KK+6)          
20       KK = KK + 9
        RETURN
        END
