         PROGRAM TEST
C
C***************** TESTING PROGRAM OF THE COMPLEX 3-DIM FFT ROUTINE
C                                       C3FFT
C
      INTEGER NPRNT,YPRNT
      PARAMETER(NPRNT=0,YPRNT=0)
      PARAMETER(MAXTEST=1)
      PARAMETER(NSIZE=151)
C
      INTEGER IDERR,FACERR,TBERR
      PARAMETER (IDERR=1,FACRR=2,TBERR=3)
      INTEGER SIZERR,DFFTERR,IFFTERR
      PARAMETER (SIZERR=10,DFFTERR=100,IFFTERR=200)
C
      PRINT '(/////)'
      PRINT '(6X,''TEST OF COMPLEX THREE DIMENSIONAL MIXED-RADIX FFT ROU
     $TINES'',//)'
      PRINT '(1X,''INPUT PARAMETERS ARE :'',/,
     $T6,''ID'',T14,''NL'',T22,''NM'',T30,''NN'',T38,''II''
     $,T46,''JJ'',T54,''KK'',T62,''IOPT'',T68,''IORD'')'
C
         DO 1 NTEST=1,MAXTEST
C
C***************** GENERATE MATRIX DIMENSIONS NL,NM,NN
C                  AND ALL THE OTHER INPUT PARAMETERS
C
           CALL GENERA(ID,NL,NM,NN,II,JJ,KK,IOPT,IORD,NSIZE)
           PRINT'(/9(I7,1X))',ID,NL,NM,NN,II,JJ,KK,IOPT,IORD
C
C
C***************** CALL TESTING ROUTINE
C
           CALL TESTFT(ID,NL,NM,NN,II,JJ,KK,IOPT,IORD,IEXIT,NPRNT)
C
           IF(IEXIT.NE.0)THEN
               PRINT '(//,''Fail- exit code '',I7)',IEXIT
             IF(IEXIT.EQ.IFFTERR.OR.IEXIT.EQ.DFFTERR.OR.
     $          IEXIT.EQ.TBERR                           )THEN
C
C***************** IN CASE OF ERROR, PRINT THE CONTENT OF DATA MATRIX
C                  ONLY IF SMALL ENOUGH (I.E. MAX(NL,NM) LESS THAN 32)
C
               IF(MAX(NL,NM).LE.32)THEN
                 CALL TESTFT(ID,NL,NM,NN,II,JJ,KK,IOPT,IORD,IEXIT,YPRNT)
               ENDIF
               STOP
             ENDIF
           ENDIF
C
1        CONTINUE
C
         END
C
      SUBROUTINE GENERA(ID,NL,NM,NN,II,JJ,KK,IOPT,IORD,NSIZE)
C
C***************** THIS ROUTINE PERFORMS A RANDOM GENERATION OF THE
C                  INPUT DATA DIMENSIONS NL,NM,NN AS WELL AS OF ALL
C                  THE CALLING PARAMETERS OF THE COMPLEX 3-DIM FFT
C                  ROUTINE C3FFT
C
      REAL MAXDIM
      PARAMETER (MAXDIM=512)
C
       L2MAX=(LOG(REAL(NSIZE*2)))/LOG(2.0)
       L3MAX=(LOG(REAL(NSIZE*3)))/LOG(3.0)
       L5MAX=(LOG(REAL(NSIZE*5)))/LOG(5.0)
C
10    CONTINUE
C
      L2=DRAND48()*L2MAX
      L3=DRAND48()*L3MAX
      L5=DRAND48()*L5MAX
C
      M2=DRAND48()*L2MAX
      M3=DRAND48()*L3MAX
      M5=DRAND48()*L5MAX
C
      N2=DRAND48()*L2MAX
      N3=DRAND48()*L3MAX
      N5=DRAND48()*L5MAX
C
      NL=2**L2*3**L3*5**L5
      IF(NL.EQ.1) NL=2
      NM=2**M2*3**M3*5**M5
      IF(NM.EQ.1) NM=2
      NN=2**N2*3**N3*5**N5
      IF(NN.EQ.1) NN=2
C
      NPL=DRAND48()*2
      ID=NL+NPL
C
C***************** GENERATED DIMENSIONS SHOULD NOT EXCEED
C                  MAXIMUM STORAGE CAPACITY
C
      IF(ID*NM*NN.GT.NSIZE*NSIZE*NSIZE) GO TO 10
      IF(NL.GT.MAXDIM.OR.NM.GT.MAXDIM.OR.NN.GT.MAXDIM) GO TO 10
C
      II=(2.0*DRAND48())*NL
      JJ=(2.0*DRAND48())*NM
      KK=(2.0*DRAND48())*NN
C
      IORD=DRAND48()*2
      IOPT=DRAND48()*2
C
C
C***************** THE PRESENT VERSION IS SUITED TO PUBLICATION
C                  THEREFORE THE INPUT PARAMETERS ARE FORCED TO
C                  ASSUME THE FOLLOWING VALUES..........
C
 1000 continue
        write(*,*) ' Please give matrix dimension'
        read(*,*) nl, nm, nn
        if( nl .le. 0 .or. nm .le. 0 .or. nn .le. 0 .or.
     $       nl*nm*nn .gt. (NSIZE-1)*NSIZE**2) THEN
           write(*,*)  NL,'is out of range -- try again.'
           goto 1000
        end if
           
      ID=Nl
      II=MAX(2,II)
      JJ=MAX(2,JJ)
      KK=MAX(2,KK)
      IOPT=1
      IORD=1
C
      RETURN
      END
C
      SUBROUTINE TESTFT(ID,NL,NM,NN,II,JJ,KK,IOPT,IORD,IEXIT,NPRINT)
C
C***************** THIS ROUTINE PERFORMS THE TESTS ON C3FFT
C
      PARAMETER (NSIZE=151)
      INTEGER SIZERR,DFFTERR,IFFTERR
      PARAMETER (SIZERR=10,DFFTERR=100,IFFTERR=200)
C
      INTEGER IDERR,FACERR,TBERR
      PARAMETER (IDERR=1,FACERR=2,TBERR=3)
C
      COMPLEX A(NSIZE*NSIZE*NSIZE),B((NSIZE+1)*NSIZE*NSIZE)
      INTEGER WL(4*NSIZE+14),WM(4*NSIZE*(NSIZE+1)+14),WN(4*NSIZE+14)
      INTEGER IWORK(NSIZE), cloops
      DIMENSION NDIMS(3)
C
      DATA EPS/1.E-7/
C
C***************** SET IEXIT TO 1 (MFFTP SHOULD RESET IT TO 0)
C
      IEXIT=1
C
C***************** TEST IF INPUT PARAMETERS EXCEED MAXIMUM
C                   WORK AREA SPACE; AN ERROR MESSAGE IS ISSUED
C                   AND EXECUTION IS STOPPED
C
      MXWORK=2*(NM+NM*ID*IOPT)
      IWSIZE=4*NSIZE*(NSIZE+1)+14
      IF(MXWORK.GT.IWSIZE) THEN
       IEXIT=SIZERR
       RETURN
      ENDIF
C
C***************** FILL INPUT MATRICES WITH
C                  MONOCHROMATIC SIGNALS
      DO 10 I=1,ID*NM*NN
         A(I)=0
10    CONTINUE
      NDIMS(1)=NL
      NDIMS(2)=NM
      NDIMS(3)=NN

      CALL FILL(A,ID,NL,NM,NN,II,JJ,KK)
C
C***************** THE PRESENT VERSION IS SUITED TO PUBLICATION
C                  THEREFORE THE INPUT MATRIX IS NOW PRINTED
C
      PRINT '(///,26X,''INPUT DATA MATRIX'')'
      IF(NPRINT.NE.0) CALL SHOWRES(A,ID,NL,NM,NN)
C
C***************** INITIALIZE FFT PACKAGE
C
C
C       CALL C3FFT(A,ID,NL,NM,NN,WL,WM,WN,IOPT,0,IORD,IWORK,IEXIT)
C       IF(IEXIT.NE.0)RETURN
C
C***************** CALCULATE DIRECT FOURIER TRANSFORM OF A
C     
       cloops=MAX(1,100000/(nl*nm*nn))
       IEXIT=0
       call stimer
       do itest=1, cloops
          ISIG=-1
C          CALL C3FFT(A,ID,NL,NM,NN,WL,WM,WN,IOPT,ISIG,IORD,IWORK,IEXIT)
          CALL FFT3D(A,B,NDIMS,ISIG)
          IF(IEXIT.NE.0)RETURN
C
C***************** INVERSE FFT OF A IS NOW COMPUTED
C
          ISIG=1
          CALL FFT3D(A,B,NDIMS,ISIG)
C          CALL C3FFT(A,ID,NL,NM,NN,WL,WM,WN,IOPT,ISIG,IORD,IWORK,IEXIT)
          IF(IEXIT.NE.0)RETURN
C     
C***************** NORMALIZE AFTER TWO TRANSFORMS
C
          CALL TMSCON(A,ID,NL,NM,NN,1.0/(REAL(NL*NM*NN)))
       end do
       call ptimer(dble(nl),dble(cloops) )

      CALL FTEST(A,ID,NL,NM,NN,II,JJ,KK)
      ERROR=XNORM(A,ID,NL,NM,NN)
C
       IF(ERROR.GT.EPS .OR. ERROR .NE. ERROR) THEN
         IEXIT=IFFTERR
         WRITE(*,*) 'Test failed, error =',error
         RETURN
       ENDIF
        IEXIT=0
C
C***************** THE TEST HAS BEEN NOW SATISFACTORILY
C                  COMPLETED. EVENTUALLY PRINT HERE
C                  HOW HAPPY YOU ARE.........
C
        END
C
      INTEGER FUNCTION IPERM(I,INDEX,NELEM)
C
C***************** THIS FUNCTION CALCULATES THE POSITION OF
C                  INDEXES WHEN ORDERING IS NOT REQUESTED
C
      INTEGER INDEX(-14:*)
      IPERM=INDEX(3*NELEM+I)
      END
C
       SUBROUTINE FILL (A,ID,NL,NM,NN,II,JJ,KK)
C
C***************** FILL INPUT MATRIX WITH MONOCHROMATIC
C                  SIGNAL OF FREQUENCY (II,JJ,KK)
C
       COMPLEX A(0:ID-1,0:NM-1,0:NN-1),CI,CJ,CK
C
       PI2=ATAN(1.)*8.
       PI2L=PI2/NL
       PI2M=PI2/NM
       PI2N=PI2/NN
C
          DO 10 I=0,NL-1
            CI=CMPLX(COS(PI2L*I*II),SIN(PI2L*I*II))
              DO 11 J=0,NM-1
                CJ=CMPLX(COS(PI2M*J*JJ),SIN(PI2M*J*JJ))
                  DO 12   K=0,NN-1
                   CK= CMPLX(COS(PI2N*K*KK),SIN(PI2N*K*KK))
C
                   A(I,J,K)=CI*CJ*CK
C
  12               CONTINUE
  11          CONTINUE
  10       CONTINUE
C
      END
       SUBROUTINE FTEST (A,ID,NL,NM,NN,II,JJ,KK)
C
C***************** TEST INPUT MATRIX AGAINST MONOCHROMATIC
C                  SIGNAL OF FREQUENCY (II,JJ,KK)
C
       COMPLEX A(0:ID-1,0:NM-1,0:NN-1),CI,CJ,CK
C
       PI2=ATAN(1.)*8.
       PI2L=PI2/NL
       PI2M=PI2/NM
       PI2N=PI2/NN
       EPS=1E-7
C
          DO 10 I=0,NL-1
            CI=CMPLX(COS(PI2L*I*II),SIN(PI2L*I*II))
              DO 11 J=0,NM-1
                CJ=CMPLX(COS(PI2M*J*JJ),SIN(PI2M*J*JJ))
                  DO 12   K=0,NN-1
                   CK= CI*CJ*CMPLX(COS(PI2N*K*KK),SIN(PI2N*K*KK))
C
D                   IF( ABS(A(I,J,K)-CK) .GT. EPS ) THEN
D                      WRITE(*,*) ' FAIL ',I,J,K,A(I,J,K),CK
D                   END IF
                   A(I,J,K)=A(I,J,K)-CK
C
  12               CONTINUE
  11          CONTINUE
  10       CONTINUE
C
      END
        SUBROUTINE SHOWRES(A,ID,NL,NM,NN)
C
C***************** PRINT THE MARIX ELEMENTS OF A
C
        COMPLEX A(ID,NM,NN)
C
        DO 100 II=1,NN
        WRITE(6,'(/)')
        DO 100 JJ=1,NL
 100    WRITE(6,'(16F7.3)') (A(JJ,K,II),K=1,NM)
        RETURN
        END
        SUBROUTINE TMSCON(A,ID,NL,NM,NN,X)
C
C***************** MULTIPLY MATRIX A BY THE REAL NUMBER X
C
        COMPLEX A(ID,NM,NN)
        REAL X

        DO 1 I=1,NN
        DO 1 J=1,NM
        DO 1 K=1,ID
          A(K,J,I)=A(K,J,I)*X
1       CONTINUE
        END
C
        FUNCTION XNORM(A,ID,NL,NM,NN)
C
C***************** CALCULATE MAXIMUM NORM OF MATRIX A
C
        COMPLEX A(ID,NM,NN)
        XNORM=0
        DO 1 N=1,NN
        DO 1 M=1,NM
        DO 1 L=1,NL
         XNORM=MAX(XNORM,ABS(A(L,M,N)))
1       CONTINUE
        END
C
        SUBROUTINE MINDEL(A,ID,NL,NM,NN,II,JJ,KK,X)
C
C***************** SUBTRACT THE REAL NUMBER X BY MATRIX A
C
        COMPLEX A(0:ID-1,0:NM-1,0:NN-1)
C
        A(II,JJ,KK)=A(II,JJ,KK)-X
C
        END
C
        SUBROUTINE PLSDEL(A,ID,NL,NM,NN,II,JJ,KK,X)
C
C***************** THE REAL NUMBER X IS ADDED TO MATRIX A
C
        COMPLEX A(0:ID-1,0:NM-1,0:NN-1)
C
        A(II,JJ,KK)=A(II,JJ,KK)+X
C
        END

