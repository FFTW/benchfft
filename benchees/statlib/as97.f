      SUBROUTINE FORRT(X, M)
C
C     ALGORITHM AS 97  APPL. STATIST. (1976) VOL.25, NO. 2
C
C     Forward discrete Fourier transform in one dimension of real
C     data using complex transform subroutine FASTG.
C
C     X = array of real input data, type real, dimension M.
C     M = length of the transform, must be a power of 2.
C     The minimum length is 8, maximum 2**21.
C
C     The result is placed in X as described in the text of the paper.
C
C     Auxiliary routines required: SCRAG (or SCRAM) & FASTG from AS 83,
C     but with SCRAG modified as described on page 168 of the paper for
C     this algorithm.
C
      REAL X(M)
      DATA ZERO/0.0/, QUART/0.25/, HALF/0.5/, ONE/1.0/, ONE5/1.5/,
     *	TWO/2.0/, FOUR/4.0/
C
C     Check for valid transform size.
C
      II = 8
      DO 2 K = 3, 21
	IPOW = K
	IF (II .EQ. M) GO TO 3
	II = II * 2
    2 CONTINUE
C
C     If this point is reached, an illegal size was specified.
C
      RETURN
    3 PIE = FOUR * ATAN(ONE)
C
C     Separate odd and even parts into two halves.
C     First bit reverse the whole array of length M.
C
      CALL SCRAG(X, M, IPOW)
C
C     Next bit reverse the half arrays separately.
C
      N = M / 2
      JPOW = IPOW - 1
      CALL SCRAG(X, N, JPOW)
      CALL SCRAG(X(N+1), N, JPOW)
C
C     Faster alternative to the two lines above to SCRAM.
C     	CALL SCRAM(X, X(N+1), N, JPOW)
C
C     Now do the transform.
C
      CALL FASTG(X, X(N+1), N, 1)
C
C     Unscramble the transform results.
C
      CALL SCRAG(X, N, JPOW)
      CALL SCRAG(X(N+1), N, JPOW)
C
C     Faster alternative to the two lines above to SCRAM.
C     	CALL SCRAM(X, X(N+1), N, JPOW)
C
      NN = N / 2
C
C     Now unravel the result; first the special cases.
C
      Z = HALF * (X(1) + X(N+1))
      X(N+1) = HALF * (X(1) - X(N+1))
      X(1) = Z
      NN1 = NN + 1
      NN2 = NN1 + N
      X(NN1) = HALF * X(NN1)
      X(NN2) = -HALF * X(NN2)
      Z = PIE / N
      BCOS = -TWO * (SIN(Z / TWO) **2)
      BSIN = SIN(Z)
      UN = ONE
      VN = ZERO
      DO 4 K = 2, NN
	Z = UN * BCOS + VN * BSIN + UN
	VN = VN * BCOS - UN * BSIN + VN
	SAVE1 = ONE5 - HALF * (Z * Z + VN * VN)
	UN = Z * SAVE1
	VN = VN * SAVE1
	KI = N + K
	L = N + 2 - K
	LI = N + L
	AN = QUART * (X(K) + X(L))
	BN = QUART * (X(KI) - X(LI))
	CN = QUART * (X(KI) + X(LI))
	DN = QUART * (X(L) - X(K))
	XN = UN * CN - VN * DN
	YN = UN * DN + VN * CN
	X(K) = AN + XN
	X(KI) = BN + YN
	X(L) = AN - XN
	X(LI) = YN - BN
    4 CONTINUE
      RETURN
      END
C

      SUBROUTINE REVRT(X, M)
C
C     ALGORITHM AS 97.1  APPL. STATIST. (1976) VOL.25, NO. 2
C
C     Inverse discrete Fourier transform in one dimension of real
C     data using complex transform subroutine FASTG.
C
C     X = array of Fourier components as output from subroutine FORRT,
C         type real, dimension M.	
C     M = length of the inverse transform, must be a power of 2.
C     The minimum length is 8, maximum 2**21.
C
C     Auxiliary routines required: SCRAG & FASTG from AS 83, but
C     with SCRAG modified as described on page 168 of the paper for
C     this algorithm.
C
      REAL X(M)
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/, ONE5/1.5/,
     *	TWO/2.0/, FOUR/4.0/
C
C     Check for valid transform size.
C
      II = 8
      DO 2 K = 3, 21
	IPOW = K
	IF (II .EQ. M) GO TO 3
	II = II * 2
    2 CONTINUE
C
C     If this point is reached, an illegal size was specified.
C
      RETURN
    3 PIE = FOUR * ATAN(ONE)
      N = M / 2
      NN = N / 2
C
C     Undo the spectrum into that of two interleaved series.
C     First, the special cases.
C
      Z = X(1) + X(N+1)
      X(N+1) = X(1) - X(N+1)
      X(1) = Z
      NN1 = NN + 1
      NN2 = NN1 + N
      X(NN1) = TWO * X(NN1)
      X(NN2) = -TWO * X(NN2)
      Z = PIE / N
      BCOS = -TWO * (SIN(Z / TWO) **2)
      BSIN = SIN(Z)
      UN = ONE
      VN = ZERO
      DO 4 K = 2, NN
	Z = UN * BCOS + VN * BSIN + UN
	VN = VN * BCOS - UN * BSIN + VN
	SAVE1 = ONE5 - HALF * (Z * Z + VN * VN)
	UN = Z * SAVE1
	VN = VN * SAVE1
	KI = N + K
	L = N + 2 - K
	LI = N + L
	AN = X(K) + X(L)
	BN = X(KI) - X(LI)
	PN = X(K) - X(L)
	QN = X(KI) + X(LI)
	CN = UN * PN + VN * QN
	DN = UN * QN - VN * PN
	X(K) = AN - DN
	X(KI) = BN + CN
	X(L) = AN + DN
	X(LI) = CN - BN
    4 CONTINUE
C
C     Now do the inverse transform
C
      CALL FASTG(X, X(N+1), N, -1)
C
C     Now undo the order - the half arrays are already bit reversed;
C     bit reverse the whole array.
C
      CALL SCRAG(X, M, IPOW)
      RETURN
      END
C
	subroutine fastg(xreal, ximag, n, itype)
c
c       Algorithm AS 83.2 Appl. Statist. (1975) vol.24, no.1
c
c       Radix 4 complex discrete fast Fourier transform without
c       unscrambling, suitable for convolutions or other applications
c       which do not require unscrambling.   Called by subroutine
c       FASTF which also does the unscrambling.
c
	real    xreal(n), ximag(n)
	data    zero, half, one, one5, two, four
     +          /0.0,  0.5, 1.0,  1.5, 2.0,  4.0/
	pi = four * atan(one)
	ifaca = n / 4
	if (itype .eq. 0) return
	if (itype .gt. 0) go to 5
c
c       ITYPE < 0 indicates inverse transform required.
c       Calculate conjugate.
c
	do 4 k = 1, n
    4   ximag(k) = -ximag(k)
c
c       Following code is executed for IFACA = N/4, N/16, N/64, ...
c       until IFACA <= 1.
c
    5   ifcab = ifaca * 4
	z = pi / ifcab
	bcos = -two * sin(z)**2
	bsin = sin(two * z)
	cw1 = one
	sw1 = zero
	do 10 litla = 1, ifaca
	  do 8 i0 = litla, n, ifcab
	    i1 = i0 + ifaca
	    i2 = i1 + ifaca
	    i3 = i2 + ifaca
	    xs0 = xreal(i0) + xreal(i2)
	    xs1 = xreal(i0) - xreal(i2)
	    ys0 = ximag(i0) + ximag(i2)
	    ys1 = ximag(i0) - ximag(i2)
	    xs2 = xreal(i1) + xreal(i3)
	    xs3 = xreal(i1) - xreal(i3)
	    ys2 = ximag(i1) + ximag(i3)
	    ys3 = ximag(i1) - ximag(i3)
	    xreal(i0) = xs0 + xs2
	    ximag(i0) = ys0 + ys2
	    x1 = xs1 + ys3
	    y1 = ys1 - xs3
	    x2 = xs0 - xs2
	    y2 = ys0 - ys2
	    x3 = xs1 - ys3
	    y3 = ys1 + xs3
	    if (litla .eq. 1) then
	      xreal(i2) = x1
	      ximag(i2) = y1
	      xreal(i1) = x2
	      ximag(i1) = y2
	      xreal(i3) = x3
	      ximag(i3) = y3
	    else
	      xreal(i2) = x1 * cw1 + y1 * sw1
	      ximag(i2) = y1 * cw1 - x1 * sw1
	      xreal(i1) = x2 * cw2 + y2 * sw2
	      ximag(i1) = y2 * cw2 - x2 * sw2
	      xreal(i3) = x3 * cw3 + y3 * sw3
	      ximag(i3) = y3 * cw3 - x3 * sw3
	    end if
    8     continue
c
c       Calculate a new set of twiddle factors.
c
	  if (litla .lt. ifaca) then
	    z = cw1 * bcos - sw1 * bsin + cw1
	    sw1 = bcos * sw1 + bsin * cw1 + sw1
	    tempr = one5 - half * (z * z + sw1 * sw1)
	    cw1 = z * tempr
	    sw1 = sw1 * tempr
	    cw2 = cw1 * cw1 - sw1 * sw1
	    sw2 = two * cw1 * sw1
	    cw3 = cw1 * cw2 - sw1 * sw2
	    sw3 = cw1 * sw2 + cw2 * sw1
	  end if
   10   continue
	if (ifaca .le. 1) go to 14
c
c       Set up the transform split for the next stage.
c
	ifaca = ifaca / 4
	if (ifaca .gt. 0) go to 5
c
c       Radix 2 calculation, if needed.
c
	if (ifaca .lt. 0) return
	do 13 k = 1, n, 2
	  tempr = xreal(k) + xreal(k+1)
	  xreal(k+1) = xreal(k) - xreal(k+1)
	  xreal(k) = tempr
	  tempr = ximag(k) + ximag(k+1)
	  ximag(k+1) = ximag(k) - ximag(k+1)
	  ximag(k) = tempr
   13   continue
   14   if (itype .lt. 0) then
c
c       Inverse transform; conjugate the result.
c
	  do 16 k = 1, n
   16     ximag(k) = -ximag(k)
	  return
	end if
c
c       Forward transform
c
	z = one / n
	do 18 k = 1, n
	  xreal(k) = xreal(k) * z
	  ximag(k) = ximag(k) * z
   18   continue
c
	return
	end
c
c
c
	subroutine scrag(xreal, n, ipow)
c
c       Algorithm AS 83.3 Appl. Statist. (1975) vol.24, no.1
c ***   MODIFIED FOR USE WITH AS 97 ***
c
c       Subroutine for unscrambling FFT data.
c
	real    xreal(n)
	integer l(19)
	equivalence (l1,l(1)), (l2,l(2)), (l3,l(3)), (l4,l(4)),
     +          (l5,l(5)), (l6,l(6)), (l7,l(7)), (l8,l(8)), (l9,l(9)),
     +          (l10,l(10)), (l11,l(11)), (l12,l(12)), (l13,l(13)),
     +          (l14,l(14)), (l15,l(15)), (l16,l(16)), (l17,l(17)),
     +          (l18,l(18)), (l19,l(19))
c
	ii = 1
	itop = 2 ** (ipow - 1)
	i = 20 - ipow
	do 5 k = 1, i
    5   l(k) = ii
	l0 = ii
	i = i + 1
	do 6 k = i, 19
	  ii = ii * 2
	  l(k) = ii
    6   continue
c
	ii = 0
	do 9 j1 = 1, l1, l0
	  do 9 j2 = j1, l2, l1
	    do 9 j3 = j2, l3, l2
	      do 9 j4 = j3, l4, l3
		do 9 j5 = j4, l5, l4
		  do 9 j6 = j5, l6, l5
		    do 9 j7 = j6, l7, l6
		      do 9 j8 = j7, l8, l7
			do 9 j9 = j8, l9, l8
			  do 9 j10 = j9, l10, l9
			    do 9 j11 = j10, l11, l10
			      do 9 j12 = j11, l12, l11
				do 9 j13 = j12, l13, l12
				  do 9 j14 = j13, l14, l13
				    do 9 j15 = j14, l15, l14
				      do 9 j16 = j15, l16, l15
					do 9 j17 = j16, l17, l16
					  do 9 j18 = j17, l18, l17
					    do 9 j19 = j18, l19, l18
					      j20 = j19
					      do 9 i = 1, 2
						ii = ii + 1
						if (ii .lt. j20) then
c
c       J20 is the bit-reverse of II pairwise interchange.
c
						  tempr = xreal(ii)
						  xreal(ii) = xreal(j20)
						  xreal(j20) = tempr
						end if
						j20 = j20 + itop
    9   continue
c
	return
	end

