c-------------------------------------------------------------c
c                                                             c
c  Subroutine sffteu( x, y, n, m, itype )                     c
c                                                             c
c  This routine is a slight modification of a complex split   c
c  radix FFT routine presented by C.S. Burrus.  The original  c
c  program header is shown below.                             c
c                                                             c
c  Arguments:                                                 c
c     x - real array containing real parts of transform       c
c              sequence (in/out)                              c
c     y - real array containing imag parts of transform       c
c              sequence (in/out)                              c
c     n - integer length of transform (in)                    c
c     m - integer such that n = 2**m  (in)                    c
c     itype - integer job specifier (in)                      c
c              itype .ne. -1 --> foward transform             c
c              itype .eq. -1 --> backward transform           c
c                                                             c
c  The forward transform computes                             c
c     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)             c
c                                                             c
c  The backward transform computes                            c
c     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)      c
c                                                             c
c                                                             c
c  Requires standard FORTRAN functions - sin, cos             c
c                                                             c
c  Steve Kifowit, 9 July 1997                                 c
c                                                             c
C-------------------------------------------------------------C
C  A Duhamel-Hollman Split-Radix DIF FFT                      C
C  Reference:  Electronics Letters, January 5, 1984           C
C  Complex input and output in data arrays X and Y            C
C  Length is N = 2**M                                         C
C                                                             C
C  C.S. Burrus          Rice University         Dec 1984      C
C-------------------------------------------------------------C
c
      SUBROUTINE SFFTEU( X, Y, N, M, ITYPE )
      INTEGER  N, M, ITYPE
      REAL  X(*), Y(*)
      INTEGER  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
      REAL  TWOPI, E, A, A3, CC1, SS1, CC3, SS3
      REAL  R1, R2, S1, S2, S3, XT
      INTRINSIC  SIN, COS
      PARAMETER  ( TWOPI = 6.2831853071795864769 )
c
      IF ( N .EQ. 1 ) RETURN
c
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 1, I = 1, N
	    Y(I) = - Y(I)
 1       CONTINUE
      ENDIF
c
      N2 = 2 * N
      DO 10, K = 1, M-1
	 N2 = N2 / 2
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
	       R1 = X(I0) - X(I2)
	       X(I0) = X(I0) + X(I2)
	       R2 = X(I1) - X(I3)
	       X(I1) = X(I1) + X(I3)
	       S1 = Y(I0) - Y(I2)
	       Y(I0) = Y(I0) + Y(I2)
	       S2 = Y(I1) - Y(I3)
	       Y(I1) = Y(I1) + Y(I3)
	       S3 = R1 - S2
	       R1 = R1 + S2
	       S2 = R2 - S1
	       R2 = R2 + S1
	       X(I2) = R1 * CC1 - S2 * SS1
	       Y(I2) = - S2 * CC1 - R1 * SS1
	       X(I3) = S3 * CC3 + R2 * SS3
	       Y(I3) = R2 * CC3 - S3 * SS3
 30         CONTINUE
	    IS = 2 * ID - N2 + J
	    ID = 4 * ID
	    IF ( IS .LT. N ) GOTO 40
 20      CONTINUE
 10   CONTINUE
c
C--------LAST STAGE, LENGTH-2 BUTTERFLY ----------------------C
c
      IS = 1
      ID = 4
 50   DO 60, I0 = IS, N, ID
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
      IF ( IS .LT. N ) GOTO 50
c
C-------BIT REVERSE COUNTER-----------------------------------C
c
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
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 2, I = 1, N
	    X(I) = X(I) / N
	    Y(I) = - Y(I) / N
 2       CONTINUE
      ENDIF
c
      RETURN
c
c ... End of subroutine SFFTEU ...
c
      END