      SUBROUTINE SFFTFU( X, Y, N, M, ITYPE )
c
c  Decimation-in-time radix-2 split-radix complex FFT
c
c  Arguments:
c     X - real part of data sequence (in/out)
c     Y - imag part of data sequence (in/out)
c     N,M - integers such that N = 2**M (in)
c     ITYPE - integer transform type (in)
c         ITYPE .ne. -1 --> forward transform
c         ITYPE .eq. -1 --> backward transform
c
c  The forward transform computes
c     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)
c
c  The backward transform computes
c     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)
c
c
c  Here is the original program header...
c
C-------------------------------------------------------------C
C  A Duhamel-Hollman Split-Radix DIT FFT                      C
C  Reference:  Electronics Letters, January 5, 1984           C
C  Complex input and output in data arrays X and Y            C
C  Length is N = 2**M                                         C
C                                                             C
C  H.V. Sorensen        Rice University         Dec 1984      C
C-------------------------------------------------------------C
c
c ... Scalar arguments ...
      INTEGER  N, M, ITYPE
c ... Array arguments ...
      REAL  X(*), Y(*)
c ... Local scalars ...
      INTEGER  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
      REAL  TWOPI, A, A3, E, XT, R1, R2, R3, S1, S2
      REAL  CC1, CC3, SS1, SS3
c ... Parameters ...
      PARAMETER  ( TWOPI = 6.283185307179586476925287 )
c ... Intrinsic functions ...
      INTRINSIC  SIN, COS
c
c ... Exe. statements ...
c
c ... Quick return ...
      IF ( N .EQ. 1 ) RETURN
c
c ... Conjugate if necessary ...
      IF ( ITYPE .EQ. -1 ) THEN
         DO 1, I = 1, N
            Y(I) = - Y(I)
 1       CONTINUE
      ENDIF
c
c ... Bit reversal permutation ...
 100  J = 1
      N1 = N - 1
      DO 104, I = 1, N1
	 IF ( I .GE. J ) GOTO 101
	 XT = X(J)
	 X(J) = X(I)
	 X(I) = XT
	 XT = Y(J)
	 Y(J) = Y(I)
	 Y(I) = XT
 101     K = N / 2
 102     IF ( K .GE. J ) GOTO 103
	 J = J - K
	 K = K / 2
	 GOTO 102
 103     J = J + K
 104  CONTINUE
c
c ... Length two transforms ...
      IS = 1
      ID = 4
 70   DO 60, I0 = IS, N, ID
	 I1 = I0 + 1
	 R1 = X(I0)
	 X(I0) = R1 + X(I1)
	 X(I1) = R1 - X(I1)
	 R1 = Y(I0)
	 Y(I0) = R1 + Y(I1)
	 Y(I1) = R1 - Y(I1)
 60   CONTINUE
      IS = 2 * ID - 1
      ID = 4 * ID
      IF ( IS .LT. N ) GOTO 70
c
c ... L shaped butterflies ...
      N2 = 2
      DO 10, K = 2, M
	 N2 = N2 * 2
	 N4 = N2 / 4
	 E = TWOPI / N2
	 A = 0.0
	 DO 20, J = 1, N4
	    A3 = 3 * A
	    CC1 = COS( A )
	    SS1 = SIN( A )
	    CC3 = COS( A3 )
	    SS3 = SIN( A3 )
	    A = J * E
	    IS = J
	    ID = 2 * N2
 40         DO 30, I0 = IS, N-1, ID
	       I1 = I0 + N4
	       I2 = I1 + N4
	       I3 = I2 + N4
	       R1 = X(I2) * CC1 + Y(I2) * SS1
	       S1 = Y(I2) * CC1 - X(I2) * SS1
	       R2 = X(I3) * CC3 + Y(I3) * SS3
	       S2 = Y(I3) * CC3 - X(I3) * SS3
	       R3 = R1 + R2
	       R2 = R1 - R2
	       R1 = S1 + S2
	       S2 = S1 - S2
	       X(I2) = X(I0) - R3
	       X(I0) = X(I0) + R3
	       X(I3) = X(I1) - S2
	       X(I1) = X(I1) + S2
	       Y(I2) = Y(I0) - R1
	       Y(I0) = Y(I0) + R1
	       Y(I3) = Y(I1) + R2
	       Y(I1) = Y(I1) - R2
 30         CONTINUE
	    IS = 2 * ID - N2 + J
	    ID = 4 * ID
	    IF ( IS .LT. N ) GOTO 40
 20      CONTINUE
 10   CONTINUE
c
c ... Conjugate and normalize if necessary ...
      IF ( ITYPE .EQ. -1 ) THEN
         DO 2, I = 1, N
	    Y(I) = - Y(I) / N
	    X(I) = X(I) / N
 2       CONTINUE
      ENDIF
c
      RETURN
c
c ... End of subroutine SFFTFU ...
c
      END
