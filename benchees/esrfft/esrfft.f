C-------------------------------------------------------------C
C  An Extended Split-Radix DIF FFT                            C
C  Complex input and output in data arrays X and Y            C
C  Length is N = 2**M                                         C
C                                                             C
C  D. Takahashi, University of Tsukuba, Sep. 2002             C
C                                                             C
C  Reference: D. Takahashi, "An extended split-radix FFT      C
C             algorithm," IEEE Signal Processing Lett.,       C
C             vol. 8, pp. 145-147, May 2001.                  C
C-------------------------------------------------------------C
      SUBROUTINE ESRFFT(X,Y,IBETA,N,M,IV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(N),Y(N)
      INTEGER*4 IBETA(N/2)
      PARAMETER (C21=0.70710678118654752D0,
     +           TWOPI=6.28318530717958647D0)
C-------Beta vector table-------
      IBETA(1)=1
      IBETA(2)=0
      IBETA(3)=0
      IBETA(4)=0
      ID=1
      DO 30 K=3,M-1
        DO 20 J=1,4
          IS=(J+3)*ID+1
          DO 10 I=IS,IS+ID-1
            IBETA(I)=IBETA(I-IS+1)
   10     CONTINUE
   20   CONTINUE
        ID=2*ID
   30 CONTINUE
C
      IF (IV .EQ. -1) THEN
        DO 40 I=1,N
          Y(I)=-Y(I)
   40   CONTINUE
      END IF
C-------L shaped butterflies-------
      L=1
      N2=N
      DO 70 K=1,M-2
        N8=N2/8
        E=-TWOPI/DBLE(N2)
        A=0.0D0
        DO 60 J=1,N8
          A3=3.0D0*A
          A5=5.0D0*A
          A7=7.0D0*A
          CC1=DCOS(A)
          SS1=DSIN(A)
          CC3=DCOS(A3)
          SS3=DSIN(A3)
          CC5=DCOS(A5)
          SS5=DSIN(A5)
          CC7=DCOS(A7)
          SS7=DSIN(A7)
          A=DBLE(J)*E
          DO 50 I=1,L
            IF (IBETA(I) .EQ. 1) THEN
              I0=(I-1)*N2+J
              I1=I0+N8
              I2=I1+N8
              I3=I2+N8
              I4=I3+N8
              I5=I4+N8
              I6=I5+N8
              I7=I6+N8
              X0=X(I0)-X(I4)
              Y0=Y(I0)-Y(I4)
              X1=X(I1)-X(I5)
              Y1=Y(I1)-Y(I5)
              X2=Y(I2)-Y(I6)
              Y2=X(I6)-X(I2)
              X3=X(I3)-X(I7)
              Y3=Y(I3)-Y(I7)
              X(I0)=X(I0)+X(I4)
              Y(I0)=Y(I0)+Y(I4)
              X(I1)=X(I1)+X(I5)
              Y(I1)=Y(I1)+Y(I5)
              X(I2)=X(I2)+X(I6)
              Y(I2)=Y(I2)+Y(I6)
              X(I3)=X(I3)+X(I7)
              Y(I3)=Y(I3)+Y(I7)
              U0=X0+C21*(X1-X3)
              V0=Y0+C21*(Y1-Y3)
              U1=X0-C21*(X1-X3)
              V1=Y0-C21*(Y1-Y3)
              U2=X2+C21*(Y1+Y3)
              V2=Y2-C21*(X1+X3)
              U3=X2-C21*(Y1+Y3)
              V3=Y2+C21*(X1+X3)
              X(I4)=CC1*(U0+U2)-SS1*(V0+V2)
              Y(I4)=CC1*(V0+V2)+SS1*(U0+U2)
              X(I5)=CC5*(U1+U3)-SS5*(V1+V3)
              Y(I5)=CC5*(V1+V3)+SS5*(U1+U3)
              X(I6)=CC3*(U1-U3)-SS3*(V1-V3)
              Y(I6)=CC3*(V1-V3)+SS3*(U1-U3)
              X(I7)=CC7*(U0-U2)-SS7*(V0-V2)
              Y(I7)=CC7*(V0-V2)+SS7*(U0-U2)
            END IF
   50     CONTINUE
   60   CONTINUE
        L=2*L
        N2=N2/2
   70 CONTINUE
C-------Length four butterflies-------
      DO 80 I=1,N/4
        IF (IBETA(I) .EQ. 1) THEN
          I0=4*I-3
          I1=I0+1
          I2=I1+1
          I3=I2+1
          X0=X(I0)-X(I2)
          Y0=Y(I0)-Y(I2)
          X1=Y(I1)-Y(I3)
          Y1=X(I3)-X(I1)
          X(I0)=X(I0)+X(I2)
          Y(I0)=Y(I0)+Y(I2)
          X(I1)=X(I1)+X(I3)
          Y(I1)=Y(I1)+Y(I3)
          X(I2)=X0+X1
          Y(I2)=Y0+Y1
          X(I3)=X0-X1
          Y(I3)=Y0-Y1
        END IF
   80 CONTINUE
C-------Length two butterflies-------
      DO 90 I=1,N/2
        IF (IBETA(I) .EQ. 1) THEN
          I0=2*I-1
          I1=I0+1
          X0=X(I0)-X(I1)
          Y0=Y(I0)-Y(I1)
          X(I0)=X(I0)+X(I1)
          Y(I0)=Y(I0)+Y(I1)
          X(I1)=X0
          Y(I1)=Y0
        END IF
   90 CONTINUE
C-------Digit reverse counter-------
      J=1
      DO 110 I=1,N-1
        IF (I .LT. J) THEN
          X0=X(I)
          Y0=Y(I)
          X(I)=X(J)
          Y(I)=Y(J)
          X(J)=X0
          Y(J)=Y0
        END IF
        K=N/2
  100   IF (K .LT. J) THEN
          J=J-K
          K=K/2
          GO TO 100
        END IF
        J=J+K
  110 CONTINUE
C
      IF (IV .EQ. -1) THEN
        DO 120 I=1,N
          X(I)=X(I)/DBLE(N)
          Y(I)=-Y(I)/DBLE(N)
  120   CONTINUE
      END IF
      RETURN
      END
