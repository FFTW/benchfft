C-----------------------------------------------------------------------
C This is not the original program!  
C
C This file was modified by Matteo Frigo to remove the assumption that
C local variables are SAVEd.  See dsp.tar.gz for the original program.
C
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
C MAIN PROGRAM: TIME-EFFICIENT RADIX-4 FAST FOURIER TRANSFORM
C AUTHOR:       L. ROBERT MORRIS
C               DEPARTMENT OF SYSTEMS ENGINEERING AND COMPUTING SCIENCE
C               CARLETON UNIVERSITY, OTTAWA, CANADA K1S 5B6
C
C INPUT:        THE ARRAY "A" CONTAINS THE DATA TO BE TRANSFORMED
C-----------------------------------------------------------------------
C
C     TEST PROGRAM FOR AUTOGEN RADIX-4 FFT
C
      DIMENSION A(2048),B(2048)
      COMMON /AA/A
C
      IOUTD=I1MACH(2)
C
C     COMPUTE DFT AND IDFT FOR N = 64, 256, AND 1024 COMPLEX POINTS
C
      DO 1 MM=3,5
      N=4**MM
      DO 2 J=1,N
      A(2*J-1)=UNI(0)
      A(2*J  )=UNI(0)
      B(2*J-1)=A(2*J-1)
 2    B(2*J  )=A(2*J)
C
C     FORWARD DFT
C
      CALL RADIX4(MM,1,-1)
C
      IF(MM.NE.3) GO TO 5
C
C     LIST DFT INPUT, OUTPUT FOR N = 64 ONLY
C
      WRITE(IOUTD,98)
      WRITE(IOUTD,100)
      DO 3 J=1,N
      WRITE(IOUTD,96) B(2*J-1),B(2*J),A(2*J-1),A(2*J)
 3    CONTINUE
C
C     INVERSE DFT
C
 5    CALL RADIX4(MM,0, 1)
C
C     LIST DFT INPUT, IDFT OUTPUT FOR N = 64 ONLY
C
      IF(MM.NE.3) GO TO 7
C
      WRITE(IOUTD,99)
      WRITE(IOUTD,100)
      DO 6 J=1,N
      WRITE(IOUTD,96) B(2*J-1),B(2*J),A(2*J-1),A(2*J)
 6    CONTINUE
C
C     CALCULATE RMS ERROR
C
 7    ERR=0.0
      DO 8 J=1,N
 8    ERR=ERR+(A(2*J-1)-B(2*J-1))**2+(A(2*J)-B(2*J))**2
      ERR=SQRT(ERR/FLOAT(N))
      WRITE(IOUTD,97) MM,ERR
 1    CONTINUE
C
 96   FORMAT(1X,4(F10.6,2X))
 97   FORMAT(1X,20H   RMS ERROR FOR M =,I2,4H IS ,E14.6/)
 98   FORMAT(1X,43H       DFT INPUT                DFT  OUTPUT/)
 99   FORMAT(1X,43H       DFT INPUT                IDFT OUTPUT/)
 100  FORMAT(1X,44H   REAL        IMAG       REAL          IMAG/)
      STOP
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE:  RADIX4
C COMPUTES FORWARD OR INVERSE COMPLEX DFT VIA RADIX-4 FFT.
C USES AUTOGEN TECHNIQUE TO YIELD TIME EFFICIENT PROGRAM.
C-----------------------------------------------------------------------
C
      SUBROUTINE RADIX4(MM,IFLAG,JFLAG)
C
C       INPUT:
C             MM = POWER OF 4 (I.E., N = 4**MM COMPLEX POINT TRANSFORM)
C             (MM.GE.2 AND MM.LE.5)
C
C          IFLAG = 1 ON FIRST PASS FOR GIVEN N
C                = 0 ON SUBSEQUENT PASSES FOR GIVEN N
C
C          JFLAG = -1 FOR FORWARD TRANSFORM
C                = +1 FOR INVERSE TRANSFORM
C
C       INPUT/OUTPUT:
C              A = ARRAY OF DIMENSIONS 2*N WITH REAL AND IMAGINARY PARTS
C                  OF DFT INPUT/OUTPUT IN ODD, EVEN ARRAY COMPONENTS.
C
C       FOR OPTIMAL TIME EFFICIENCY, COMMON IS USED TO PASS ARRAYS.
C       THIS MEANS THAT DIMENSIONS OF ARRAYS A, IX, AND T CAN BE
C       MODIFIED TO REFLECT MAXIMUM VALUE OF N = 4**MM TO BE USED. NOTE
C       THAT ARRAY "IX" IS ALSO DIMENSIONED IN SUBROUTINE "RAD4SB".
C
C    I.E.,    A(    ) IX(  ) T(   )
C
C       M =2      32     38     27
C       M<=3     128    144    135
C       M<=4     512    658    567
C       M<=5    2048   2996   2295
C
      DIMENSION A(2048),IX(2996),T(2295)
      DIMENSION NFAC(11),NP(209)
      COMMON NTYPL,KKP,INDEX,IXC
      COMMON /AA/A
      COMMON /XX/IX
      SAVE N, XP, NTOT, N2, N1TEST, N2TEST, N3TEST, NSPAN
      SAVE IBASE, ISN, INC, RAD, PI, C707, CM141, C383, C924, CM924
      SAVE C541, CM541, C131, CM131, NT, KS, KSPAN, JC, RADF, I, M, K
      SAVE NFAC, NP, S2, C2, S1, C1, S3, C3, T
C
C     CHECK FOR MM<2 OR MM>5
C
      IF(MM.LT.2.OR.MM.GT.5)STOP
C
C     INITIALIZE ON FIRST PASS """""""""""""""""""""""""""""""""""""""""
C
      IF(IFLAG.EQ.1) GO TO 9999
C
C     FAST FOURIER TRANSFORM START #####################################
C
 8885 KSPAN=2*4**MM
      IF(JFLAG.EQ.1) GO TO 8887
C
C     CONJUGATE DATA FOR FORWARD TRANSFORM
C
      DO 8886 J=2,N2,2
 8886 A(J)=-A(J)
      GO TO 8889
C
C     MULTIPLY DATA BY N**(-1) IF INVERSE TRANSFORM
C
 8887 DO 8888 J=1,N2,2
      A(J)=A(J)*XP
 8888 A(J+1)=A(J+1)*XP
 8889 I=3
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8),IT
C***********************************************************************
C
C                                   8 MULTIPLY BUTTERFLY
C
 1    KK=IX(I)
C
 11   K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
C
      AKP=A(KK)+A(K2)
      AKM=A(KK)-A(K2)
      AJP=A(K1)+A(K3)
      AJM=A(K1)-A(K3)
      A(KK)=AKP+AJP
C
      BKP=A(KK+1)+A(K2+1)
      BKM=A(KK+1)-A(K2+1)
      BJP=A(K1+1)+A(K3+1)
      BJM=A(K1+1)-A(K3+1)
      A(KK+1)=BKP+BJP
C
      BJP=BKP-BJP
C
      A(K2+1)=(AKP+BJP-AJP)*C707
      A(K2)=A(K2+1)+BJP*CM141
C
      BKP=BKM+AJM
      AKP=AKM-BJM
C
      AC0=(AKP+BKP)*C924
      A(K1+1)=AC0+AKP*CM541
      A(K1)  =AC0+BKP*CM131
C
      BKM=BKM-AJM
      AKM=AKM+BJM
C
      AC0=(AKM+BKM)*C383
      A(K3+1)=AC0+AKM*C541
      A(K3)  =AC0+BKM*CM131
C
      I=I+1
      KK=IX(I)
      IF (KK) 111,111,11
 111  I=I+2
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8), IT
C***********************************************************************
C
C                                   4 MULTIPLY BUTTERFLY
C
 2    KK=IX(I)
C
 22   K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
C
      AKP=A(KK)+A(K2)
      AKM=A(KK)-A(K2)
      AJP=A(K1)+A(K3)
      AJM=A(K1)-A(K3)
      A(KK)=AKP+AJP
C
      BKP=A(KK+1)+A(K2+1)
      BKM=A(KK+1)-A(K2+1)
      BJP=A(K1+1)+A(K3+1)
      BJM=A(K1+1)-A(K3+1)
      A(KK+1)=BKP+BJP
      A(K2)=-BKP+BJP
      A(K2+1)=AKP-AJP
C
      BKP=BKM+AJM
C
      A(K1+1)=(BKP+AKM-BJM)*C707
      A(K1)=A(K1+1)+BKP*CM141
C
      AKM=AKM+BJM
C
      A(K3+1)=(AKM+AJM-BKM)*C707
      A(K3)=A(K3+1)+AKM*CM141
C
      I=I+1
      KK=IX(I)
      IF (KK) 222,222,22
 222  I=I+2
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8), IT
C***********************************************************************
C
C                                   8 MULTIPLY BUTTERFLY
C
 3    KK=IX(I)
C
 33   K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
C
      AKP=A(KK)+A(K2)
      AKM=A(KK)-A(K2)
      AJP=A(K1)+A(K3)
      AJM=A(K1)-A(K3)
      A(KK)=AKP+AJP
C
      BKP=A(KK+1)+A(K2+1)
      BKM=A(KK+1)-A(K2+1)
      BJP=A(K1+1)+A(K3+1)
      BJM=A(K1+1)-A(K3+1)
      A(KK+1)=BKP+BJP
C
      AJP=AKP-AJP
C
      A(K2+1)=(AJP+BJP-BKP)*C707
      A(K2)=A(K2+1)+AJP*CM141
C
      BKP=BKM+AJM
      AKP=AKM-BJM
C
      AC0=(AKP+BKP)*C383
      A(K1+1)=AC0+AKP*C541
      A(K1)  =AC0+BKP*CM131
C
      BKM=BKM-AJM
      AKM=AKM+BJM
C
      AC0=(AKM+BKM)*CM924
      A(K3+1)=AC0+AKM*C541
      A(K3)  =AC0+BKM*C131
C
      I=I+1
      KK=IX(I)
      IF (KK) 333,333,33
 333  I=I+2
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8), IT
C***********************************************************************
C
C                                   GENERAL 9 MULTIPLY BUTTERFLY
C
 4    KK=IX(I)
C
 44   K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
C
      AKP=A(KK)+A(K2)
      AKM=A(KK)-A(K2)
      AJP=A(K1)+A(K3)
      AJM=A(K1)-A(K3)
      A(KK)=AKP+AJP
C
      BKP=A(KK+1)+A(K2+1)
      BKM=A(KK+1)-A(K2+1)
      BJP=A(K1+1)+A(K3+1)
      BJM=A(K1+1)-A(K3+1)
      A(KK+1)=BKP+BJP
C
      AJP=AKP-AJP
      BJP=BKP-BJP
C
      J=IX(I+1)
C
      AC0=(AJP+BJP)*T(J+8)
      A(K2+1)=AC0+AJP*T(J+6)
      A(K2)  =AC0+BJP*T(J+7)
C
      BKP=BKM+AJM
      AKP=AKM-BJM
C
      AC0=(AKP+BKP)*T(J+5)
      A(K1+1)=AC0+AKP*T(J+3)
      A(K1)  =AC0+BKP*T(J+4)
C
      BKM=BKM-AJM
      AKM=AKM+BJM
C
      AC0=(AKM+BKM)*T(J+2)
      A(K3+1)=AC0+AKM*T(J)
      A(K3)  =AC0+BKM*T(J+1)
C
      I=I+2
      KK=IX(I)
      IF (KK) 444,444,44
 444  I=I+2
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8), IT
C***********************************************************************
C
C                                   0 MULTIPLY BUTTERFLY
C
 5    KK=IX(I)
C
 55   K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
C
      AKP=A(KK)+A(K2)
      AKM=A(KK)-A(K2)
      AJP=A(K1)+A(K3)
      AJM=A(K1)-A(K3)
      A(KK)=AKP+AJP
      A(K2)=AKP-AJP
C
      BKP=A(KK+1)+A(K2+1)
      BKM=A(KK+1)-A(K2+1)
      BJP=A(K1+1)+A(K3+1)
      BJM=A(K1+1)-A(K3+1)
      A(KK+1)=BKP+BJP
      A(K2+1)=BKP-BJP
C
      A(K3+1)=BKM-AJM
      A(K1+1)=BKM+AJM
      A(K3)=AKM+BJM
      A(K1)=AKM-BJM
C
      I=I+1
      KK=IX(I)
      IF (KK) 555,555,55
 555  I=I+2
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8), IT
C***********************************************************************
C
C                                      OFFSET REDUCED
C
 6    KSPAN=KSPAN/4
      I=I+2
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8), IT
C***********************************************************************
C
C                                      BIT REVERSAL (SHUFFLING)
C
 7    IP1=IX(I)
 77   IP2=IX(I+1)
      T1=A(IP2)
      A(IP2)=A(IP1)
      A(IP1)=T1
      T1=A(IP2+1)
      A(IP2+1)=A(IP1+1)
      A(IP1+1)=T1
      I=I+2
      IP1=IX(I)
      IF (IP1) 777,777,77
 777  I=I+2
      IT=IX(I-1)
      GO TO (1,2,3,4,5,6,7,8), IT
C***********************************************************************
 8    IF(JFLAG.EQ.1) GO TO 888
C
C     CONJUGATE OUTPUT IF FORWARD TRANSFORM
C
      DO 88 J=2,N2,2
 88   A(J)=-A(J)
 888  RETURN
C
C     FAST FOURIER TRANSFORM ENDS ######################################
C
C     INITIALIZATION PHASE STARTS. DONE ONLY ONCE
C
 9999 IXC=1
      N=4**MM
      XP=N
      XP=1./XP
      NTOT=N
      N2=N*2
      NSPAN=N
      N1TEST=N/16
      N2TEST=N/8
      N3TEST=(3*N)/16
      NSPAN4=NSPAN/4
      IBASE=0
      ISN=1
      INC=ISN
      RAD=8.0*ATAN(1.0)
      PI=4.*ATAN(1.0)
      C707=SIN(PI/4.)
      CM141=-2.*C707
      C383=SIN(PI/8.)
      C924=COS(PI/8.)
      CM924=-C924
      C541=C924-C383
      CM541=-C541
      C131=C924+C383
      CM131=-C131
 10   NT=INC*NTOT
      KS=INC*NSPAN
      KSPAN=KS
      JC=KS/N
      RADF=RAD*FLOAT(JC)*.5
      I=0
C
C     DETERMINE THE FACTORS OF N
C     ALL FACTORS MUST BE 4 FOR THIS VERSION
C
      M=0
      K=N
 15   M=M+1
      NFAC(M)=4
      K=K/4
 20   IF(K-(K/4)*4.EQ.0) GO TO 15
      KT=1
      IF(N.GE.256) KT=2
      KSPAN0=KSPAN
      NTYPL=0
C
 100  NDELTA=KSPAN0/KSPAN
      INDEX=0
      SD=RADF/FLOAT(KSPAN)
      CD=2.0*SIN(SD)**2
      SD=SIN(SD+SD)
      KK=1
      I=I+1
C
C     TRANSFORM FOR A FACTOR OF 4
C
      KSPAN=KSPAN/4
      IX(IXC)=0
      IX(IXC+1)=6
      IXC=IXC+2
C
 410  C1=1.0
      S1=0.0
 420  K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
      IF(S1.EQ.0.0) GO TO 460
 430  IF(KSPAN.NE.NSPAN4) GO TO 431
      T(IBASE+5)=-(S1+C1)
      T(IBASE+6)=C1
      T(IBASE+4)=S1-C1
      T(IBASE+8)=-(S2+C2)
      T(IBASE+9)=C2
      T(IBASE+7)=S2-C2
      T(IBASE+2)=-(S3+C3)
      T(IBASE+3)=C3
      T(IBASE+1)=S3-C3
      IBASE=IBASE+9
C
 431  KKP=(KK-1)*2
      IF(INDEX.NE.N1TEST) GO TO 150
      CALL RAD4SB(1)
      GO TO 5035
 150  IF(INDEX.NE.N2TEST) GO TO 160
      CALL RAD4SB(2)
      GO TO 5035
 160  IF(INDEX.NE.N3TEST) GO TO 170
      CALL RAD4SB(3)
      GO TO 5035
 170  CALL RAD4SB(4)
 5035 KK=K3+KSPAN
      IF(KK.LE.NT) GO TO 420
 440  INDEX=INDEX+NDELTA
      C2=C1-(CD*C1+SD*S1)
      S1=(SD*C1-CD*S1)+S1
      C1=C2
      C2=C1*C1-S1*S1
      S2=C1*S1+C1*S1
      C3=C2*C1-S2*S1
      S3=C2*S1+S2*C1
      KK=KK-NT+JC
      IF(KK.LE.KSPAN) GO TO 420
      KK=KK-KSPAN+INC
      IF(KK.LE.JC) GO TO 410
      IF(KSPAN.EQ.JC) GO TO 800
      GO TO 100
 460  KKP=(KK-1)*2
      CALL RAD4SB(5)
 5050 KK=K3+KSPAN
      IF(KK.LE.NT) GO TO 420
      GO TO 440
C
 800  IX(IXC)=0
      IX(IXC+1)=7
      IXC=IXC+2
C
C     COMPUTE PARAMETERS TO PERMUTE THE RESULTS TO NORMAL ORDER
C     DONE IN TWO STEPS
C     PERMUTATION FOR SQUARE FACTORS OF N
C
      NP(1)=KS
      K=KT+KT+1
      IF(M.LT.K) K=K-1
      J=1
      NP(K+1)=JC
 810  NP(J+1)=NP(J)/NFAC(J)
      NP(K)=NP(K+1)*NFAC(J)
      J=J+1
      K=K-1
      IF(J.LT.K) GO TO 810
      K3=NP(K+1)
      KSPAN=NP(2)
      KK=JC+1
      K2=KSPAN+1
      J=1
C
C     PERMUTATION FOR SINGLE VARIATE TRANSFORM
C
 820  KKP=(KK-1)*2
      K2P=(K2-1)*2
      IX(IXC)=KKP+1
      IX(IXC+1)=K2P+1
      IXC=IXC+2
      KK=KK+INC
      K2=KSPAN+K2
      IF(K2.LT.KS) GO TO 820
 830  K2=K2-NP(J)
      J=J+1
      K2=NP(J+1)+K2
      IF(K2.GT.NP(J)) GO TO 830
      J=1
 840  IF(KK.LT.K2) GO TO 820
      KK=KK+INC
      K2=KSPAN+K2
      IF(K2.LT.KS) GO TO 840
      IF(KK.LT.KS) GO TO 830
      JC=K3
      IX(IXC)=0
      IX(IXC+1)=8
      GO TO 8885
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE:  RAD4SB
C USED BY SUBROUTINE RADIX4. NEVER DIRECTLY ACCESSED BY USER.
C-----------------------------------------------------------------------
C
      SUBROUTINE RAD4SB(NTYPE)
C
C     INPUT: NTYPE = TYPE OF BUTTERFLY INVOKED
C     OUTPUT: PARAMETERS USED BY SUBROUTINE RADIX4
C
      DIMENSION IX(2996)
      COMMON /XX/IX
      COMMON NTYPL,KKP,INDEX,IXC
      IF(NTYPE.EQ.NTYPL) GO TO 7
      IX(IXC)=0
      IX(IXC+1)=NTYPE
      IXC=IXC+2
      IF(NTYPE.NE.4) GO TO 4
      INDEXP=(INDEX-1)*9
      IX(IXC)=KKP+1
      IX(IXC+1)=INDEXP+1
      IXC=IXC+2
      GO TO 6
 4    IX(IXC)=KKP+1
      IXC=IXC+1
 6    NTYPL=NTYPE
      RETURN
 7    IF(NTYPE.NE.4) GO TO 8
      INDEXP=(INDEX-1)*9
      IX(IXC)=KKP+1
      IX(IXC+1)=INDEXP+1
      IXC=IXC+2
      RETURN
 8    IX(IXC)=KKP+1
      IXC=IXC+1
      RETURN
      END
