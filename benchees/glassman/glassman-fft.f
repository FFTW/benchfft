*=======================================================================
* Single Precision Complex Fast Fourier Transform
*
*  A subroutine to compute the discrete Fourier transform by the fastest
* available algorithm for any length input series.
*
* Reference:
*        Ferguson, W., (1979),   A Simple Derivation of Glassmans's
*          General N Fast Fourier Transform, MRC Tech. Summ. Rept. 2029,
*          Math. Res. Cent. U. of Wisconsin, Madison, Wis.
*
*  REFERENCES
*  ----------
*
* Routines called:
* SPCPFT
*
* Functions called:
* MOD
* FLOAT
*
* VAX extensions:
* DO WHILE
*=======================================================================

      SUBROUTINE SPCFFT(U,N,ISIGN,WORK,INTERP)

* VARIABLES
* ---------

      IMPLICIT NONE

      LOGICAL*1
     |  INU      ! Flag for calling SUBROUTINE SPCPFT( arguments ).
     | ,SCALE    ! .TRUE.=inverse transform -- .FALSE.=forward transform

      INTEGER*4
     |  A        ! After    \
     | ,B        ! Before    >  Factors of N.
     | ,C        ! Current  /
     | ,N        ! Length of the array to be transformed.
     | ,I        ! DO loop index.
     | ,ISIGN    ! sign of transform

      REAL*4
     |  INTERP   ! interpolation factor

      COMPLEX*8
     |  U(*)            !  Vector to be transformed
     | ,WORK(*)         !  Working storage.

*     Initialize parameters.

      A = 1
      B = N
      C = 1

      INU = .TRUE.

      IF (ISIGN.EQ.1) THEN

         SCALE = .TRUE.

      ELSE IF (ISIGN.EQ.-1) THEN

         SCALE = .FALSE.

      END IF

* Calculate Fourier transform by means of Glassman's algorithm

      DO WHILE ( B .GT. 1 )

         A = C * A

* Start of routine to get factors of N

         C = 2

* Increment C until it is an integer factor of B

         DO WHILE ( MOD(B,C) .NE. 0 )

                  C = C + 1

         END DO

* Calculate largest factor of B

         B = B / C


* Call Glassman's Fast Fourier Transform routine

         IF ( INU ) THEN

            CALL SPCPFT (A,B,C,U,WORK,ISIGN)

          ELSE

            CALL SPCPFT (A,B,C,WORK,U,ISIGN)

         END IF

* Set flag to swap input & output arrays to SPCPFT

         INU = ( .NOT. INU )

      END DO

* If odd number of calls to SPCPFT swap WORK into U

      IF ( .NOT. INU ) THEN

         DO I = 1, N
            U(I) = WORK(I)
         END DO

      END IF

* Scale the output for inverse Fourier transforms.

      IF ( SCALE ) THEN

         DO I = 1, N
            U(I) = U(I) / (N/INTERP)
         END DO

      END IF


* TERMINATION
* -----------

      RETURN
      END


*=======================================================================
* Single Precision Complex Prime Factor Transform
*
*  REFERENCES
*  ----------
*
* Calling routines:
* SPCFFT
*
* Subroutines called:
* -none-
*
* Functions called:
* CMLPX
* COS
* SIN
* FLOAT
*=======================================================================

      SUBROUTINE SPCPFT( A, B, C, UIN, UOUT, ISIGN )

* VARIABLES
* ---------

      IMPLICIT NONE

      INTEGER*4
     |  ISIGN           !  ISIGN of the Fourier transform.
     | ,A               !  After   \
     | ,B               !  Before   >  Factors of N.
     | ,C               !  Current /
     | ,IA              !  \
     | ,IB              !   \  DO loop indicies.
     | ,IC              !   /
     | ,JCR             !  /
     | ,JC              !  Dummy index.

      REAL*8
     |  ANGLE

      COMPLEX*8
     |  UIN(B,C,A)      !  Input vector.
     | ,UOUT(B,A,C)     !  Output vector.
     | ,DELTA           !  Fourier transform kernel.
     | ,OMEGA           !  Multiples of exp( i TWOPI ).
     | ,SUM             !  Dummy register for addition for UOUT(B,A,C)

* ALGORITHM
* ---------

* Initialize run time parameters.


      ANGLE =6.28318530717958 / FLOAT( A * C )
      OMEGA = CMPLX( 1.0, 0.0 )

* Check the ISIGN of the transform.

      IF( ISIGN .EQ. 1 ) THEN

         DELTA = CMPLX( DCOS(ANGLE), DSIN(ANGLE) )

      ELSE

         DELTA = CMPLX( DCOS(ANGLE), -DSIN(ANGLE) )

      END IF



* Do the computations.

      DO IC = 1, C

         DO IA = 1, A

            DO IB = 1, B

               SUM = UIN( IB, C, IA )

               DO JCR = 2, C

                  JC = C + 1 - JCR
                  SUM = UIN( IB, JC, IA ) + OMEGA * SUM

               END DO

               UOUT( IB, IA, IC ) = SUM

            END DO

            OMEGA = DELTA * OMEGA

         END DO

      END DO

* TERMINATION
* -----------

      RETURN
      END


#ifdef FFT_TEST

      program test

      implicit none

      integer*4 nmax
      parameter( nmax=10000 )

      complex*8 input(nmax)
      complex*8 scratch(nmax)

      real*4 a , b

      integer*4 i
      integer*4 n

      n = 0
 1000 format( 2e25.15 )
      do i=1, nmax

         read( 5 ,1000 ,end=10 ) input(i)
         n = n + 1

      end do

   10 continue

      call spcfft( input ,n ,0 ,scratch ,1.0 )

      do i=1 ,n

         write( 6 ,1000 ) real(input(i)) ,imag(input(i))

      end do

      call spcfft( input ,n ,1 ,scratch ,1.0 )

      do i=1 ,n

         write( 6 ,1000 ) real(input(i)) ,imag(input(i))

      end do

      end

#endif
#ifdef TEST

      program test

      implicit none

      real*4 data_in(10000)
      real*4 data_out(10000)

      integer*4 i
      integer*4 m
      integer*4 n


      read( 5 ,* ) m ,n

      do i=1,m
         read( 5 ,* ,end=10 ) data_in(i)
      end do

   10 continue

      write( 6 ,* ) '#input'
      do i=1 ,m
         write( 6 ,'(2f25.10)' ) real(i-1)/(m) ,data_in(i)
      end do

      call resample( data_in ,m ,data_out ,n )

      write( 6 ,* ) '#output'
      do i=1 ,n
         write( 6 ,'(2f25.10)' ) real(i-1)/(n) ,data_out(i)
      end do

      end

#endif
