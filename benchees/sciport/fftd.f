C-------------------------------------------------------------     ************
C                                                                     CRFFT2
C                                                                  ************
      SUBROUTINE CRFFT2(INIT,IX,N,CX,CWORK,CY)
C
C     THE STOCKHAM AUTO-SORT FFT
C
      INTEGER POPCNT
      DOUBLE COMPLEX CWORK(1),CX(1)
      DOUBLE PRECISION CY(1)
C
      L2 = POPCNT(N)
      IF (L2 .NE. 1) CALL ABORT('CRFFT2')
      N2 = N + 3
      NN = N/2
      IF (INIT .EQ. 0) GOTO 10
      IF (N .LT. 4) CALL ABORT('CRFFT2')
      CALL CRBLE1(NN,CWORK(N2))
      RETURN
  10  CONTINUE
      NS = N/4
      CALL CRFORM(IX,NS,NN,CX,CWORK(NN+1),CWORK(N2))
      CALL CROCK1(NS,CWORK,CWORK(NN+1))
      IF (IX .LE. 0) GOTO 50
      LS = 2
      NS = NS/2
  20  CONTINUE
      IF (NS .EQ. 1) GOTO 30
         CALL CROCK2(LS,NS,CWORK(NN+1),CWORK(1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      IF (NS .EQ. 1) GOTO 130
	 CALL CROCK2(LS,NS,CWORK(1),CWORK(NN+1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      GOTO 20
  30  CONTINUE
         CALL CROCK2(LS,NS,CY,CWORK(1),CWORK(N2))
      RETURN
 130  CONTINUE
         CALL CROCK2(LS,NS,CY,CWORK(NN+1),CWORK(N2))
      RETURN
  50  CONTINUE
         LS = 2
         NS = NS/2
  60  CONTINUE
      IF (NS .EQ. 1) GOTO 70
	 CALL CROCK3(LS,NS,CWORK(NN+1),CWORK(1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      IF (NS .EQ. 1) GOTO 120
         CALL CROCK3(LS,NS,CWORK(1),CWORK(NN+1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      GOTO 60
  70  CONTINUE
         CALL CROCK3(LS,NS,CY,CWORK(1),CWORK(N2))
      RETURN
 120  CONTINUE
         CALL CROCK3(LS,NS,CY,CWORK(NN+1),CWORK(N2))
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CRFORM
C                                                 ************
      SUBROUTINE CRFORM(IX,NS,NDIV2,CX,C,CH2)
C
      DOUBLE COMPLEX CX(1),WYK1,C(NS,2),WYK
      DOUBLE PRECISION CH2(NDIV2,2)
C
      IF (IX .GT. 0) GOTO 50
      K = NS + 1
      DO 10 I=1, NS
         WYK = CONJG(CX(NDIV2-I+2))
         C(I,1)= CX(I)+WYK + (CX(I) - WYK) * DCMPLX(CH2(I,2),CH2(I,1))
         WYK1 = CONJG(CX(NDIV2-K+2))
	 C(I,2)= CX(K)+WYK1+ (CX(K) -WYK1) * DCMPLX(CH2(K,2),CH2(K,1))
         K = K + 1
  10  CONTINUE
      RETURN
  50  CONTINUE
      K = NS + 1
      DO 20 I=1, NS
         WYK = CONJG(CX(NDIV2-I+2))
         C(I,1)= CX(I)+WYK + (CX(I) - WYK) * DCMPLX(-CH2(I,2),CH2(I,1))
         WYK1 = CONJG(CX(NDIV2-K+2))
         C(I,2)= CX(K)+WYK1 +(CX(K) -WYK1) * DCMPLX(-CH2(K,2),CH2(K,1))
         K = K + 1
  20  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CROCK1
C                                                 ************
      SUBROUTINE CROCK1(NS,C,CH)
C
      DOUBLE COMPLEX C(NS,2),CH(NS,2)
C
      DO 300 J=1, NS
         C(J,1) = CH(J,1) + CH(J,2)
         C(J,2) = CH(J,1) - CH(J,2)
 300  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CROCK2
C                                                 ************
      SUBROUTINE CROCK2(LS,NS,C,CH,CH2)
C
      DOUBLE COMPLEX WYK,C(NS,LS,2),CH(NS,2,LS)
      DOUBLE PRECISION CH2(2,NS,LS,2)
C
      IF (LS .GT. NS) GOTO 20
      DO 200 I=1, LS
	 DO 200 J=1, NS
            WYK = DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2)) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 200  CONTINUE
      RETURN
  20  CONTINUE
      DO 400 J=1, NS
         DO 400 I=1, LS
            WYK = DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2)) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 400  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CROCK3
C                                                 ************
      SUBROUTINE CROCK3(LS,NS,C,CH,CH2)
C
      DOUBLE COMPLEX WYK,C(NS,LS,2),CH(NS,2,LS)
      DOUBLE PRECISION CH2(2,NS,LS,2)
C
      IF (LS .GT. NS) GOTO 30
      DO 600 I=1, LS
         DO 600 J=1, NS
            WYK =CONJG(DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2))) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 600  CONTINUE
      RETURN
  30  CONTINUE
      DO 800 J=1, NS
         DO 800 I=1, LS
            WYK =CONJG(DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2))) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 800  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CRBLE1
C                                                 ************
      SUBROUTINE CRBLE1(NN,WORK)
C
      DOUBLE PRECISION WORK(NN,2),P2,TWOPI
      DATA TWOPI /6.28318530717958647692/
C
      N = 2 * NN
      P2 = TWOPI/N
      DO 10 I=1, NN
         WORK(I,1) = COS(P2 * (I-1))
         WORK(I,2) = SIN(P2 * (I-1))
  10  CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     RCFFT2
C                                                                  ************
      SUBROUTINE RCFFT2(INIT,IX,N,CX,CWORK,CY)
C
C     THE STOCKHAM AUTO-SORT FFT
C
      INTEGER POPCNT
      DOUBLE COMPLEX CWORK(1),CY(1),CX(1)
C
      L2 = POPCNT(N)
      IF (L2 .NE. 1) CALL ABORT('RCFFT2')
      N2 = N + 3
      NN = N/2
      IF (INIT .EQ. 0) GOTO 10
      IF (N .LT. 4) CALL ABORT('RCFFT2')
      CALL RABLE1(NN,CWORK(N2))
      RETURN
  10  CONTINUE
      NS = N/4
      CALL RTOCK1(NS,CWORK,CX)
      IF (IX .LE. 0) GOTO 50
      LS = 2
      NS = NS/2
  20  CONTINUE
         CALL RTOCK2(LS,NS,CWORK(NN+1),CWORK(1),CWORK(N2))
      IF (NS .EQ. 1) GOTO 30
         LS = LS + LS
         NS = NS/2
	 CALL RTOCK2(LS,NS,CWORK(1),CWORK(NN+1),CWORK(N2))
      IF (NS .EQ. 1) GOTO 130
         LS = LS + LS
         NS = NS/2
      GOTO 20
  30  CONTINUE
         CALL RCONV1(N,CY,CWORK(NN+1),CWORK(N2))
      RETURN
 130  CONTINUE
         CALL RCONV1(N,CY,CWORK,CWORK(N2))
      RETURN
  50  CONTINUE
         LS = 2
         NS = NS/2
  60  CONTINUE
	 CALL RTOCK3(LS,NS,CWORK(NN+1),CWORK(1),CWORK(N2))
      IF (NS .EQ. 1) GOTO 70
         LS = LS + LS
         NS = NS/2
         CALL RTOCK3(LS,NS,CWORK(1),CWORK(NN+1),CWORK(N2))
      IF (NS .EQ. 1) GOTO 120
         LS = LS + LS
         NS = NS/2
      GOTO 60
  70  CONTINUE
         CALL RCONV2(N,CY,CWORK(NN+1),CWORK(N2))
      RETURN
 120  CONTINUE
         CALL RCONV2(N,CY,CWORK,CWORK(N2))
      RETURN
      END

C-----------------------------------------------  ************
C                                                    RTOCK2
C                                                 ************
      SUBROUTINE RTOCK2(LS,NS,C,CH,CH2)
C
      DOUBLE COMPLEX WYK,C(NS,LS,2),CH(NS,2,LS)
      DOUBLE PRECISION CH2(2,NS,LS,2)
C
      IF (LS .GT. NS) GOTO 20
      DO 200 I=1, LS
         DO 200 J=1, NS
            WYK = DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2)) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 200  CONTINUE
      RETURN
  20  CONTINUE
      DO 400 J=1, NS
	 DO 400 I=1, LS
            WYK = DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2)) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 400  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    RTOCK1
C                                                 ************
      SUBROUTINE RTOCK1(NS,C,CH)
C
      DOUBLE COMPLEX C(NS,1,2),CH(NS,2,1)
C
      DO 300 J=1, NS
         C(J,1,1) = CH(J,1,1) + CH(J,2,1)
         C(J,1,2) = CH(J,1,1) - CH(J,2,1)
 300  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    RTOCK3
C                                                 ************
      SUBROUTINE RTOCK3(LS,NS,C,CH,CH2)
C
      DOUBLE COMPLEX WYK,C(NS,LS,2),CH(NS,2,LS)
      DOUBLE PRECISION CH2(2,NS,LS,2)
C
      IF (LS .GT. NS) GOTO 30
      DO 600 I=1, LS
         DO 600 J=1, NS
            WYK =CONJG(DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2))) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 600  CONTINUE
      RETURN
  30  CONTINUE
      DO 800 J=1, NS
         DO 800 I=1, LS
            WYK =CONJG(DCMPLX(CH2(1,1,I,1),CH2(1,1,I,2))) * CH(J,2,I)
            C(J,I,1) = CH(J,1,I) + WYK
            C(J,I,2) = CH(J,1,I) - WYK
 800  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    RABLE1
C                                                 ************
      SUBROUTINE RABLE1(NN,WORK)
C
      DOUBLE PRECISION WORK(NN,2),P2,TWOPI
      DATA TWOPI /6.28318530717958647692/
C
      N = 2 * NN
      P2 = TWOPI/N
      DO 10 I=1, NN
         WORK(I,1) = COS(P2 * (I-1))
         WORK(I,2) = SIN(P2 * (I-1))
  10  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    RCONV1
C                                                 ************
      SUBROUTINE RCONV1(N,CY,C,CH)
C
      DOUBLE COMPLEX CY(1)
      DOUBLE PRECISION CH(N/2,2),P(2,1),C(2,N/2)
      DOUBLE PRECISION X,Y,Z,Z1
C
      N2 = N/2
      P(1,1) = (C(1,1) + C(2,1)) * 2
      P(2,1) = (C(1,1) - C(2,1)) * 2
      CY(1) = DCMPLX(P(1,1),0.0D0)
      CY(N2+1) = DCMPLX(P(2,1),0.0D0)
      K = N2
      DO 10 I=2, N2
         X = C(1,I)+C(1,K)
         Y = C(2,I)+C(2,K)
         Z = C(1,I)-C(1,K)
         Z1= C(2,I)-C(2,K)
         P(1,1) = X + CH(I,1) * Y + CH(I,2) * Z
         P(2,1) = Z1 + CH(I,2) * Y - CH(I,1) * Z
         CY(I) = DCMPLX(P(1,1),P(2,1))
         K = K - 1
  10  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    RCONV2
C                                                 ************
      SUBROUTINE RCONV2(N,CY,C,CH)
C
      DOUBLE COMPLEX CY(1)
      DOUBLE PRECISION CH(N/2,2),P(2,1),C(2,N/2)
      DOUBLE PRECISION X,Y,Z,Z1
C
      N2 = N/2
      P(1,1) = (C(1,1) + C(2,1)) * 2
      P(2,1) = (C(1,1) - C(2,1)) * 2
      CY(1) = DCMPLX(P(1,1),0.0D0)
      CY(N2+1) = DCMPLX(P(2,1),0.0D0)
      K = N2
      DO 10 I=2, N2
         X = C(1,I)+C(1,K)
         Y = C(2,I)+C(2,K)
         Z = C(1,I)-C(1,K)
         Z1= C(2,I)-C(2,K)
	 P(1,1) = X + CH(I,1) * Y - CH(I,2) * Z
         P(2,1) = Z1 - CH(I,2) * Y - CH(I,1) * Z
         CY(I) = DCMPLX(P(1,1),P(2,1))
         K = K - 1
  10  CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CFFT2
C                                                                  ************
      SUBROUTINE CFFT2(INIT,IX,N,CX,CWORK,CY)
C
C     THE STOCKHAM AUTO-SORT FFT
C
      INTEGER POPCNT
      DOUBLE COMPLEX CX(1),CWORK(1),CY(1)
C
C     IS N IS THE POWER OF 2 ?
C
      L2 = POPCNT(N)
      IF (L2 .NE. 1) CALL ABORT('CFFT2 ')
      N2 = N + N + 1
C
C     WANT TO GENERATE THE TABLE FOR SINE AND COSINE, IF INIT .NE. 0.
C
      NS = N/2
      IF (INIT .EQ. 0) GOTO 10
C
C     IF N = 0 OR 2, RETURN.
C     IF N = 4,8 OR 16 THEN USING TABLE TO GENERATE.(CBALE1)
C     IF N > 16 THEN USING THE NEW METHOD TO GENERATE.(CBLE2)
C
      IF (N .LT. 4) CALL ABORT('CFFT2 ')
      CALL CABLE2(NS,CWORK(N2))
      RETURN
C
C     USING THE PREVIOUS TABLE TO PERFORM FOURIER TRANSFORMATION,
C     IF INIT = 0.
C
  10  CONTINUE
      CALL CTOCK1(NS,CWORK,CX)
      IF (IX .LE. 0) GOTO 50
C
C     THE FOLLOWING STATEMENTS FOR FOURIER ANALYSIS.
C            ^I.^E. IX > 0.
C
      LS = 2
      NS = NS/2
  20  CONTINUE
      IF (NS .EQ. 1) GOTO 120
         CALL CTOCK2(LS,NS,CWORK(N+1),CWORK(1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      IF (NS .EQ. 1) GOTO 30
         CALL CTOCK2(LS,NS,CWORK(1),CWORK(N+1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      GOTO 20
  30  CONTINUE
         CALL CTOCK2(LS,NS,CY,CWORK(N+1),CWORK(N2))
      RETURN
 120  CONTINUE
	 CALL CTOCK2(LS,NS,CY,CWORK(1),CWORK(N2))
      RETURN
C
C     THE FOLLOWING STATEMENTS FOR FOURIER SYNTHESIS.
C            ^I.^E. IX < 0.
C
  50  CONTINUE
         LS = 2
         NS = NS/2
  60  CONTINUE
      IF (NS .EQ. 1) GOTO 130
         CALL CTOCK3(LS,NS,CWORK(N+1),CWORK(1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      IF (NS .EQ. 1) GOTO 70
	 CALL CTOCK3(LS,NS,CWORK(1),CWORK(N+1),CWORK(N2))
         LS = LS + LS
         NS = NS/2
      GOTO 60
  70  CONTINUE
         CALL CTOCK3(LS,NS,CY,CWORK(N+1),CWORK(N2))
      RETURN
 130  CONTINUE
         CALL CTOCK3(LS,NS,CY,CWORK(1),CWORK(N2))
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CTOCK1
C                                                 ************
      SUBROUTINE CTOCK1(NS,C,CH)
C
      DOUBLE COMPLEX C(NS,1,2),CH(NS,2,1)
C
      DO 300 J=1, NS
	 C(J,1,1) = CH(J,1,1) + CH(J,2,1)
         C(J,1,2) = CH(J,1,1) - CH(J,2,1)
 300  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CTOCK2
C                                                 ************
      SUBROUTINE CTOCK2(LS,NS,C,CH,CH2)
C
      DOUBLE COMPLEX C(NS,LS,2),CH(NS,2,LS),CH2(NS,LS)
C
      IF (LS .GT. NS) GOTO 20
      DO 200 I=1, LS
         DO 200 J=1, NS
            C(J,I,1) = CH(J,1,I) + CH2(1,I) * CH(J,2,I)
	    C(J,I,2) = CH(J,1,I) - CH2(1,I) * CH(J,2,I)
 200  CONTINUE
      RETURN
  20  CONTINUE
      DO 400 J=1, NS
         DO 400 I=1, LS
            C(J,I,1) = CH(J,1,I) + CH2(1,I) * CH(J,2,I)
            C(J,I,2) = CH(J,1,I) - CH2(1,I) * CH(J,2,I)
 400  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CTOCK3
C                                                 ************
      SUBROUTINE CTOCK3(LS,NS,C,CH,CH2)
C
      DOUBLE COMPLEX C(NS,LS,2),CH(NS,2,LS),CH2(NS,LS)
C
      IF (LS .GT. NS) GOTO 30
      DO 600 I=1, LS
         DO 600 J=1, NS
            C(J,I,1) = CH(J,1,I) + CONJG(CH2(1,I)) * CH(J,2,I)
            C(J,I,2) = CH(J,1,I) - CONJG(CH2(1,I)) * CH(J,2,I)
 600  CONTINUE
      RETURN
  30  CONTINUE
      DO 800 J=1, NS
         DO 800 I=1, LS
	    C(J,I,1) = CH(J,1,I) + CONJG(CH2(1,I)) * CH(J,2,I)
            C(J,I,2) = CH(J,1,I) - CONJG(CH2(1,I)) * CH(J,2,I)
 800  CONTINUE
      RETURN
      END

C-----------------------------------------------  ************
C                                                    CABLE2
C                                                 ************
      SUBROUTINE CABLE2(NN,WORK)
C
      DOUBLE PRECISION WORK(2,NN),P2,TWOPI
      DATA TWOPI /6.28318530717958647692/
C
      N = 2 * NN
      P2 = TWOPI/N
      DO 10 I=1, NN
         WORK(1,I) = COS(P2 * (I-1))
         WORK(2,I) = SIN(P2 * (I-1))
  10  CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                      ABORT
C                                                                  ************
      SUBROUTINE ABORT(NME)
C
C     Routine to abort execution if N is not of the form 2**I.
C     Martin J. McBride.  2/27/86.
C     General Electric CRD, Information System Operation.
C
      CHARACTER*6 NME
      DOUBLE PRECISION DIV0,DEN,NUM

      PRINT*
      PRINT*,NME,' called with N not of the form N=2**I where I=>2.'
      PRINT*
      NUM = 1.0
      DEN = 0.0
      DIV0 = NUM/DEN

      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     POPCNT
C                                                                  ************
      INTEGER FUNCTION POPCNT(N)
C
C     Routine to determine how many bits of N are 1 (is N = 2**I).
C     Martin J. McBride.  3/4/86.
C     General Electric CRD, Information System Operation.
C
      INTEGER N,I,TMP,IB

      I = 0
      TMP = N
   10 IB = MOD(TMP,2)
      IF (IB .NE. 0) I = I + 1
      TMP = INT(TMP/2)
      IF (TMP .NE. 0) GO TO 10

      POPCNT = I

      RETURN
      END
