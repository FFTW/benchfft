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
C MAIN PROGRAM: TEST PROGRAM TO EXERCISE THE WFTA SUBROUTINE
C   THE TEST WAVEFORM IS A COMPLEX EXPONENTIAL A**I WHOSE
C   TRANSFORM IS KNOWN ANALYTICALLY TO BE (1 - A**N)/(1 - A*W**K).
C
C AUTHORS:
C   JAMES H. MCCLELLAN     AND     HAMID NAWAB
C   DEPARTMENT OF ELECTRICAL ENGINEERING AND COMPUTER SCIENCE
C   MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C   CAMBRIDGE, MASS.  02139
C
C INPUTS:
C   N-- TRANSFORM LENGTH. IT MUST BE FORMED AS THE PRODUCT OF
C       RELATIVELY PRIME INTEGERS FROM THE SET:
C           2,3,4,5,7,8,9,16
C   INVRS IS THE FLAG FOR FORWARD OR INVERSE TRANSFORM.
C           INVRS = 1 YIELDS INVERSE TRANSFORM
C           INVRS .NE. 1 GIVES FORWARD TRANSFORM
C   RAD AND PHI ARE THE MAGNITUDE AND ANGLE (AS A FRACTION OF
C       2*PI/N) OF THE COMPLEX EXPONENTIAL TEST SIGNAL.
C          SUGGESTION: RAD = 0.98, PHI = 0.5.
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION PI2,PIN,XN,XJ,XT
      DIMENSION XR(5040),XI(5040)
      COMPLEX CONE,CA,CAN,CNUM,CDEN
C
C   OUTPUT WILL BE PUNCHED
C
      IOUT=I1MACH(3)
      INPUT=I1MACH(1)
      CONE=CMPLX(1.0,0.0)
      PI2=8.0D0*DATAN(1.0D0)
 50   CONTINUE
      READ(INPUT,130)N
 130  FORMAT(I5)
      WRITE(IOUT,150) N
 150  FORMAT(10H LENGTH = ,I5)
      IF(N.LE.0 .OR. N.GT.5040) STOP
C
C   ENTER A 1 TO PERFORM THE INVERSE
C
      READ(INPUT,130) INVRS
C
C   ENTER MAGNITUDE AND ANGLE (IN FRACTION OF 2*PI/N)
C   AVOID MULTIPLES OF N FOR THE ANGLE IF THE RADIUS IS
C   CLOSE TO ONE.  SUGGESTION: RAD = 0.98, PHI = 0.5.
C
      READ(INPUT,160) RAD,PHI
 160  FORMAT(2F15.10)
      XN=FLOAT(N)
      PIN=PHI
      PIN=PIN*PI2/XN
C
C   GENERATE Z**J
C
      INIT=0
      DO 200 J=1,N
      AN=RAD**(J-1)
      XJ=J-1
      XJ=XJ*PIN
      XT=DCOS(XJ)
      XR(J)=XT
      XR(J)=XR(J)*AN
      XT=DSIN(XJ)
      XI(J)=XT
      XI(J)=XI(J)*AN
 200  CONTINUE
      CAN=CMPLX(XR(N),XI(N))
      CA=CMPLX(XR(2),XI(2))
      CAN=CAN*CA
C
C   PRINT FIRST 50 VALUES OF INPUT SEQUENCE
C
      MAX=50
      IF(N.LT.50)MAX=N
      WRITE(IOUT,300)(J,XR(J),XI(J),J=1,MAX)
C
C   CALL THE WINOGRAD FOURIER TRANSFORM ALGORITHM
C
      CALL WFTA(XR,XI,N,INVRS,INIT,IERR)
C
C   CHECK FOR ERROR RETURN
C
      IF(IERR.LT.0) WRITE(IOUT,250) IERR
 250  FORMAT(1X,5HERROR,I5)
      IF(IERR.LT.0) GO TO 50
C
C   PRINT FIRST 50 VALUES OF THE TRANSFORMED SEQUENCE
C
      WRITE(IOUT,300)(J,XR(J),XI(J),J=1,MAX)
 300  FORMAT(1X,3HJ =,I3,6HREAL =,E20.12,6HIMAG =,E20.12)
C
C   CALCULATE ABSOLUTE AND RELATIVE DEVIATIONS
C
      DEVABS=0.0
      DEVREL=0.0
      DEVMX1=0.0
      DEVMX2=0.0
      CNUM=CONE-CAN
      PIN=PI2/XN
      DO 350 J=1,N
      XJ=J-1
      XJ=-XJ*PIN
      IF(INVRS.EQ.1) XJ=-XJ
      TR=DCOS(XJ)
      TI=DSIN(XJ)
      CAN=CMPLX(TR,TI)
      CDEN=CONE-CA*CAN
      CDEN=CNUM/CDEN
C
C   TRUE VALUE OF THE TRANSFORM (1. - A**N)/(1. - A*W**K),
C   WHERE A = RAD*EXP(J*PHI*(2*PI/N)), W = EXP(-J*2*PI/N).
C   FOR THE INVERSE TRANSFORM THE COMPLEX EXPONENTIAL W
C   IS CONJUGATED.
C
      TR=REAL(CDEN)
      TI=AIMAG(CDEN)
      IF(INVRS.NE.1) GO TO 330
C
C   SCALE INVERSE TRANSFORM BY 1/N
C
      TR=TR/FLOAT(N)
      TI=TI/FLOAT(N)
 330  TR=XR(J)-TR
      TI=XI(J)-TI
      DEVABS=SQRT(TR*TR+TI*TI)
      XMAG=SQRT(XR(J)*XR(J)+XI(J)*XI(J))
      DEVREL=100.0*DEVABS/XMAG
      IF(DEVABS.LE.DEVMX1)GO TO 340
      DEVMX1=DEVABS
      LABS=J-1
 340  IF(DEVREL.LE.DEVMX2)GO TO 350
      DEVMX2=DEVREL
      LREL=J-1
 350  CONTINUE
C
C   PRINT THE ABSOLUTE AND RELATIVE DEVIATIONS TOGETHER
C   WITH THEIR LOCATIONS.
C
      WRITE(IOUT,380) DEVMX1,LABS,DEVMX2,LREL
 380  FORMAT(1X,21HABSOLUTE DEVIATION = ,E20.12,9H AT INDEX,I5/
     1 1X,21HRELATIVE DEVIATION = ,F11.7,8H PERCENT,1X,9H AT INDEX,I5)
      GO TO 50
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE: INISHL
C   THIS SUBROUTINE INITIALIZES THE WFTA ROUTINE FOR A GIVEN
C   VALUE OF THE TRANSFORM LENGTH N.  THE FACTORS OF N ARE
C   DETERMINED, THE MULTIPLICATION COEFFICIENTS ARE CALCULATED
C   AND STORED IN THE ARRAY COEF(.), THE INPUT AND OUTPUT
C   PERMUTATION VECTORS ARE COMPUTED AND STORED IN THE ARRAYS
C   INDX1(.) AND INDX2(.)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE INISHL(N,COEF,XR,XI,INDX1,INDX2,IERR)
      DIMENSION COEF(1),XR(1),XI(1)
      INTEGER S1,S2,S3,S4,INDX1(1),INDX2(1),P1
      DIMENSION CO3(3),CO4(4),CO8(8),CO9(11),CO16(18),CDA(18),CDB(11),
     1CDC(9),CDD(6)
      COMMON NA,NB,NC,ND,ND1,ND2,ND3,ND4
C
C   DATA STATEMENTS ASSIGN SHORT DFT COEFFICIENTS.
C
      DATA CO4(1),CO4(2),CO4(3),CO4(4)/4*1.0/
C
      DATA CDA(1),CDA(2),CDA(3),CDA(4),CDA(5),CDA(6),CDA(7),
     1 CDA(8),CDA(9),CDA(10),CDA(11),CDA(12),CDA(13),CDA(14),
     2 CDA(15),CDA(16),CDA(17),CDA(18)/18*1.0/
C
      DATA CDB(1),CDB(2),CDB(3),CDB(4),CDB(5),CDB(6),CDB(7),CDB(8),
     1 CDB(9),CDB(10),CDB(11)/11*1.0/
C
      DATA IONCE/1/
C
C   GET MULTIPLIER CONSTANTS
C
      IF(IONCE.NE.1) GO TO 20
      CALL CONST(CO3,CO8,CO16,CO9,CDC,CDD)
 20   IONCE=-1
C
C   FOLLOWING SEGMENT DETERMINES FACTORS OF N AND CHOOSES
C   THE APPROPRIATE SHORT DFT COEFFICIENTS.
C
      IOUT=I1MACH(2)
      IERR=0
      ND1=1
      NA=1
      NB=1
      ND2=1
      NC=1
      ND3=1
      ND=1
      ND4=1
      IF(N.LE.0 .OR. N.GT.5040) GO TO 190
      IF(16*(N/16).EQ.N) GO TO 30
      IF(8*(N/8).EQ.N) GO TO 40
      IF(4*(N/4).EQ.N) GO TO 50
      IF(2*(N/2).NE.N) GO TO 70
      ND1=2
      NA=2
      CDA(2)=1.0
      GO TO 70
 30   ND1=18
      NA=16
      DO 31 J=1,18
 31   CDA(J)=CO16(J)
      GO TO 70
 40   ND1=8
      NA=8
      DO 41 J=1,8
 41   CDA(J)=CO8(J)
      GO TO 70
 50   ND1=4
      NA=4
      DO 51 J=1,4
 51   CDA(J)=CO4(J)
 70   IF(3*(N/3).NE.N) GO TO 120
      IF(9*(N/9).EQ.N) GO TO 100
      ND2=3
      NB=3
      DO 71 J=1,3
 71   CDB(J)=CO3(J)
      GO TO 120
 100  ND2=11
      NB=9
      DO 110 J=1,11
 110  CDB(J)=CO9(J)
 120  IF(7*(N/7).NE.N) GO TO 160
      ND3=9
      NC=7
 160  IF(5*(N/5).NE.N) GO TO 190
      ND4=6
      ND=5
 190  M=NA*NB*NC*ND
      IF(M.EQ.N) GO TO 250
      WRITE(IOUT,210)
 210  FORMAT(21H THIS N DOES NOT WORK)
      IERR=-1
      RETURN
C
C   NEXT SEGMENT GENERATES THE DFT COEFFICIENTS BY
C   MULTIPLYING TOGETHER THE SHORT DFT COEFFICIENTS
C
 250  J=1
      DO 300 N4=1,ND4
      DO 300 N3=1,ND3
      DO 300 N2=1,ND2
      DO 300 N1=1,ND1
      COEF(J)=CDA(N1)*CDB(N2)*CDC(N3)*CDD(N4)
      J=J+1
 300  CONTINUE
C
C   FOLLOWING SEGMENT FORMS THE INPUT INDEXING VECTOR
C
      J=1
      NU=NB*NC*ND
      NV=NA*NC*ND
      NW=NA*NB*ND
      NY=NA*NB*NC
      K=1
      DO 440 N4=1,ND
      DO 430 N3=1,NC
      DO 420 N2=1,NB
      DO 410 N1=1,NA
 405  IF(K.LE.N) GO TO 408
      K=K-N
      GO TO 405
 408  INDX1(J)=K
      J=J+1
 410  K=K+NU
 420  K=K+NV
 430  K=K+NW
 440  K=K+NY
C
C   FOLLOWING SEGMENT FORMS THE OUTPUT INDEXING VECTOR
C
      M=1
      S1=0
      S2=0
      S3=0
      S4=0
      IF(NA.EQ.1) GO TO 530
 520  P1=M*NU-1
      IF((P1/NA)*NA.EQ.P1) GO TO 510
      M=M+1
      GO TO 520
 510  S1=P1+1
 530  IF(NB.EQ.1) GO TO 540
      M=1
 550  P1=M*NV-1
      IF((P1/NB)*NB.EQ.P1) GO TO 560
      M=M+1
      GO TO 550
 560  S2=P1+1
 540  IF(NC.EQ.1) GO TO 630
      M=1
 620  P1=M*NW-1
      IF((P1/NC)*NC.EQ.P1) GO TO 610
      M=M+1
      GO TO 620
 610  S3=P1+1
 630  IF(ND.EQ.1) GO TO 660
      M=1
 640  P1=M*NY-1
      IF((P1/ND)*ND.EQ.P1) GO TO 650
      M=M+1
      GO TO 640
 650  S4=P1+1
 660  J=1
      DO 810 N4=1,ND
      DO 810 N3=1,NC
      DO 810 N2=1,NB
      DO 810 N1=1,NA
      INDX2(J)=S1*(N1-1)+S2*(N2-1)+S3*(N3-1)+S4*(N4-1)+1
 900  IF(INDX2(J).LE.N) GO TO 910
      INDX2(J)=INDX2(J)-N
      GO TO 900
 910  J=J+1
 810  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE: CONST
C   COMPUTES THE MULTIPLIERS FOR THE VARIOUS MODULES
C-----------------------------------------------------------------------
C
      SUBROUTINE CONST(CO3,CO8,CO16,CO9,CDC,CDD)
      DOUBLE PRECISION DTHETA,DTWOPI,DSQ32,DSQ2
      DOUBLE PRECISION DCOS1,DCOS2,DCOS3,DCOS4
      DOUBLE PRECISION DSIN1,DSIN2,DSIN3,DSIN4
      DIMENSION CO3(3),CO8(8),CO16(18),CO9(11),CDC(9),CDD(6)
      DTWOPI=8.0D0*DATAN(1.0D0)
      DSQ32=DSQRT(0.75D0)
      DSQ2=DSQRT(0.5D0)
C
C   MULTIPLIERS FOR THE THREE POINT MODULE
C
      CO3(1)=1.0
      CO3(2)=-1.5
      CO3(3)=-DSQ32
C
C   MULTIPLIERS FOR THE FIVE POINT MODULE
C
      DTHETA=DTWOPI/5.0D0
      DCOS1=DCOS(DTHETA)
      DCOS2=DCOS(2.0D0*DTHETA)
      DSIN1=DSIN(DTHETA)
      DSIN2=DSIN(2.0D0*DTHETA)
      CDD(1)=1.0
      CDD(2)=-1.25
      CDD(3)=-DSIN1-DSIN2
      CDD(4)=0.5*(DCOS1-DCOS2)
      CDD(5)=DSIN1-DSIN2
      CDD(6)=DSIN2
C
C
C   MULTIPLIERS FOR THE SEVEN POINT MODULE
C
      DTHETA=DTWOPI/7.0D0
      DCOS1=DCOS(DTHETA)
      DCOS2=DCOS(2.0D0*DTHETA)
      DCOS3=DCOS(3.0D0*DTHETA)
      DSIN1=DSIN(DTHETA)
      DSIN2=DSIN(2.0D0*DTHETA)
      DSIN3=DSIN(3.0D0*DTHETA)
      CDC(1)=1.0
      CDC(2)=-7.0D0/6.0D0
      CDC(3)=-(DSIN1+DSIN2-DSIN3)/3.0D0
      CDC(4)=(DCOS1+DCOS2-2.0D0*DCOS3)/3.0D0
      CDC(5)=(2.0D0*DCOS1-DCOS2-DCOS3)/3.0D0
      CDC(6)=-(2.0D0*DSIN1-DSIN2+DSIN3)/3.0D0
      CDC(7)=-(DSIN1+DSIN2+2.0D0*DSIN3)/3.0D0
      CDC(8)=(DCOS1-2.0D0*DCOS2+DCOS3)/3.0D0
      CDC(9)=-(DSIN1-2.0D0*DSIN2-DSIN3)/3.0D0
C
C   MULTIPLIERS FOR THE EIGHT POINT MODULE
C
      CO8(1)=1.0
      CO8(2)=1.0
      CO8(3)=1.0
      CO8(4)=-1.0
      CO8(5)=1.0
      CO8(6)=-DSQ2
      CO8(7)=-1.0
      CO8(8)=DSQ2
C
C   MULTIPLIERS FOR THE NINE POINT MODULE
C
      DTHETA=DTWOPI/9.0D0
      DCOS1=DCOS(DTHETA)
      DCOS2=DCOS(2.0D0*DTHETA)
      DCOS4=DCOS(4.0D0*DTHETA)
      DSIN1=DSIN(DTHETA)
      DSIN2=DSIN(2.0D0*DTHETA)
      DSIN4=DSIN(4.0D0*DTHETA)
      CO9(1)=1.0
      CO9(2)=-1.5
      CO9(3)=-DSQ32
      CO9(4)=0.5
      CO9(5)=(2.0D0*DCOS1-DCOS2-DCOS4)/3.0D0
      CO9(6)=(DCOS1-2.0D0*DCOS2+DCOS4)/3.0D0
      CO9(7)=(DCOS1+DCOS2-2.0D0*DCOS4)/3.0D0
      CO9(8)=-(2.0D0*DSIN1+DSIN2-DSIN4)/3.0D0
      CO9(9)=-(DSIN1+2.0D0*DSIN2+DSIN4)/3.0D0
      CO9(10)=-(DSIN1-DSIN2-2.0D0*DSIN4)/3.0D0
      CO9(11)=-DSQ32
C
C   MULTIPLIERS FOR THE SIXTEEN POINT MODULE
C
      DTHETA=DTWOPI/16.0D0
      DCOS1=DCOS(DTHETA)
      DCOS3=DCOS(3.0D0*DTHETA)
      DSIN1=DSIN(DTHETA)
      DSIN3=DSIN(3.0D0*DTHETA)
      CO16(1)=1.0
      CO16(2)=1.0
      CO16(3)=1.0
      CO16(4)=-1.0
      CO16(5)=1.0
      CO16(6)=-DSQ2
      CO16(7)=-1.0
      CO16(8)=DSQ2
      CO16(9)=1.0
      CO16(10)=-(DSIN1-DSIN3)
      CO16(11)=-DSQ2
      CO16(12)=-CO16(10)
      CO16(13)=-1.0
      CO16(14)=-(DSIN1+DSIN3)
      CO16(15)=DSQ2
      CO16(16)=-CO16(14)
      CO16(17)=-DSIN3
      CO16(18)=DCOS3
      RETURN
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE: WFTA
C WINOGRAD FOURIER TRANSFORM ALGORITHM
C-----------------------------------------------------------------------
C
      SUBROUTINE WFTA(XR,XI,N,INVRS,INIT,IERR)
      DIMENSION XR(1),XI(1)
C
C INPUTS:
C   N-- TRANSFORM LENGTH.  MUST BE FORMED AS THE PRODUCT OF
C       RELATIVELY PRIME INTEGERS FROM THE SET:
C           2,3,4,5,7,8,9,16
C       THUS THE LARGEST POSSIBLE VALUE OF N IS 5040.
C   XR(.)-- ARRAY THAT HOLDS THE REAL PART OF THE DATA
C           TO BE TRANSFORMED.
C   XI(.)-- ARRAY THAT HOLDS THE IMAGINARY PART OF THE
C           DATA TO BE TRANSFORMED.
C   INVRS-- PARAMETER THAT FLAGS WHETHER OR NOT THE INVERSE
C           TRANSFORM IS TO BE CALCULATED.  A DIVISION BY N
C           IS INCLUDED IN THE INVERSE.
C           INVRS = 1 YIELDS INVERSE TRANSFORM
C           INVRS .NE. 1 GIVES FORWARD TRANSFORM
C   INIT-- PARAMETER THAT FLAGS WHETHER OR NOT THE PROGRAM
C          IS TO BE INITIALIZED FOR THIS VALUE OF N.  THE
C          INITIALIZATION IS PERFORMED ONLY ONCE IN ORDER TO
C          TO SPEED UP THE COMPUTATION ON SUCCEEDING CALLS
C          TO THE WFTA ROUTINE, WHEN N IS HELD FIXED.
C          INIT = 0 RESULTS IN INITIALIZATION.
C   IERR-- ERROR CODE THAT IS NEGATIVE WHEN THE WFTA
C          TERMINATES INCORRECTLY.
C           0 = SUCCESSFUL COMPLETION
C          -1 = THIS VALUE OF N DOES NOT FACTOR PROPERLY
C          -2 = AN INITIALIZATION HAS NOT BEEN DONE FOR
C               THIS VALUE OF N.
C
C
C   THE FOLLOWING TWO CARDS MAY BE CHANGED IF THE MAXIMUM
C   DESIRED TRANSFORM LENGTH IS LESS THAN 5040
C
C  *********************************************************************
      DIMENSION SR(10692),SI(10692),COEF(10692)
      INTEGER INDX1(5040),INDX2(5040)
      SAVE INDX1
      SAVE INDX2
C  *********************************************************************
C
      COMMON NA,NB,NC,ND,ND1,ND2,ND3,ND4
C
C   TEST FOR INITIAL RUN
C
      IF(INIT.EQ.0) CALL INISHL(N,COEF,XR,XI,INDX1,INDX2,IERR)
C
      IF(IERR.LT.0) RETURN
      M=NA*NB*NC*ND
      IF(M.EQ.N) GO TO 100
      IERR=-2
      RETURN
C
C  ERROR(-2)-- PROGRAM NOT INITIALIZED FOR THIS VALUE OF N
C
 100  NMULT=ND1*ND2*ND3*ND4
C
C   THE FOLLOWING CODE MAPS THE DATA ARRAYS XR AND XI TO
C   THE WORKING ARRAYS SR AND SI VIA THE MAPPING VECTOR
C   INDX1(.).  THE PERMUTATION OF THE DATA FOLLOWS THE
C   SINO CORRESPONDENCE OF THE CHINESE REMAINDER THEOREM.
C
      J=1
      K=1
      INC1=ND1-NA
      INC2=ND1*(ND2-NB)
      INC3=ND1*ND2*(ND3-NC)
      DO 140 N4=1,ND
      DO 130 N3=1,NC
      DO 120 N2=1,NB
      DO 110 N1=1,NA
      IND=INDX1(K)
      SR(J)=XR(IND)
      SI(J)=XI(IND)
      K=K+1
 110  J=J+1
 120  J=J+INC1
 130  J=J+INC2
 140  J=J+INC3
C
C   DO THE PRE-WEAVE MODULES
C
      CALL WEAVE1(SR,SI)
C
C   THE FOLLOWING LOOP PERFORMS ALL THE MULTIPLICATIONS OF THE
C   WINOGRAD FOURIER TRANSFORM ALGORITHM.  THE MULTIPLICATION
C   COEFFICIENTS ARE STORED ON THE INITIALIZATION PASS IN THE
C   ARRAY COEF(.).
C
      DO 200 J=1,NMULT
      SR(J)=SR(J)*COEF(J)
      SI(J)=SI(J)*COEF(J)
 200  CONTINUE
C
C   DO THE POST-WEAVE MODULES
C
      CALL WEAVE2(SR,SI)
C
C
C   THE FOLLOWING CODE MAPS THE WORKING ARRAYS SR AND SI
C   TO THE DATA ARRAYS XR AND XI VIA THE MAPPING VECTOR
C   INDX2(.).  THE PERMUTATION OF THE DATA FOLLOWS THE
C   CHINESE REMAINDER THEOREM.
C
      J=1
      K=1
      INC1=ND1-NA
      INC2=ND1*(ND2-NB)
      INC3=ND1*ND2*(ND3-NC)
C
C   CHECK FOR INVERSE
C
      IF(INVRS.EQ.1) GO TO 400
      DO 340 N4=1,ND
      DO 330 N3=1,NC
      DO 320 N2=1,NB
      DO 310 N1=1,NA
      KNDX=INDX2(K)
      XR(KNDX)=SR(J)
      XI(KNDX)=SI(J)
      K=K+1
 310  J=J+1
 320  J=J+INC1
 330  J=J+INC2
 340  J=J+INC3
      RETURN
C
C   DIFFERENT PERMUTATION FOR THE INVERSE
C
 400  FN=FLOAT(N)
      NP2=N+2
      INDX2(1)=N+1
      DO 440 N4=1,ND
      DO 430 N3=1,NC
      DO 420 N2=1,NB
      DO 410 N1=1,NA
      KNDX=NP2-INDX2(K)
      XR(KNDX)=SR(J)/FN
      XI(KNDX)=SI(J)/FN
      K=K+1
 410  J=J+1
 420  J=J+INC1
 430  J=J+INC2
 440  J=J+INC3
      RETURN
      END
C
C-----------------------------------------------------------------------
C SUBROUTINE: WEAVE1
C   THIS SUBROUTINE IMPLEMENTS THE DIFFERENT PRE-WEAVE
C   MODULES OF THE WFTA.  THE WORKING ARRAYS ARE SR AND SI.
C   THE ROUTINE CHECKS TO SEE WHICH FACTORS ARE PRESENT
C   IN THE TRANSFORM LENGTH N = NA*NB*NC*ND AND EXECUTES
C   THE PRE-WEAVE CODE FOR THESE FACTORS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE WEAVE1(SR,SI)
      COMMON NA,NB,NC,ND,ND1,ND2,ND3,ND4
      DIMENSION Q(8),T(16)
      DIMENSION SR(1),SI(1)
      IF(NA.EQ.1) GO TO 300
      IF(NA.NE.2) GO TO 800
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 2 POINT PRE-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=2*(ND2-NB)
      NLUP23=2*ND2*(ND3-NC)
      NBASE=1
      DO 240 N4=1,ND
      DO 230 N3=1,NC
      DO 220 N2=1,NB
      NR1=NBASE+1
      T0=SR(NBASE)+SR(NR1)
      SR(NR1)=SR(NBASE)-SR(NR1)
      SR(NBASE)=T0
      T0=SI(NBASE)+SI(NR1)
      SI(NR1)=SI(NBASE)-SI(NR1)
      SI(NBASE)=T0
 220  NBASE=NBASE+2
 230  NBASE=NBASE+NLUP2
 240  NBASE=NBASE+NLUP23
 800  IF(NA.NE.8) GO TO 1600
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 8 POINT PRE-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=8*(ND2-NB)
      NLUP23=8*ND2*(ND3-NC)
      NBASE=1
      DO 840 N4=1,ND
      DO 830 N3=1,NC
      DO 820 N2=1,NB
      NR1=NBASE+1
      NR2=NR1+1
      NR3=NR2+1
      NR4=NR3+1
      NR5=NR4+1
      NR6=NR5+1
      NR7=NR6+1
      T3=SR(NR3)+SR(NR7)
      T7=SR(NR3)-SR(NR7)
      T0=SR(NBASE)+SR(NR4)
      SR(NR4)=SR(NBASE)-SR(NR4)
      T1=SR(NR1)+SR(NR5)
      T5=SR(NR1)-SR(NR5)
      T2=SR(NR2)+SR(NR6)
      SR(NR6)=SR(NR2)-SR(NR6)
      SR(NBASE)=T0+T2
      SR(NR2)=T0-T2
      SR(NR1)=T1+T3
      SR(NR3)=T1-T3
      SR(NR5)=T5+T7
      SR(NR7)=T5-T7
      T3=SI(NR3)+SI(NR7)
      T7=SI(NR3)-SI(NR7)
      T0=SI(NBASE)+SI(NR4)
      SI(NR4)=SI(NBASE)-SI(NR4)
      T1=SI(NR1)+SI(NR5)
      T5=SI(NR1)-SI(NR5)
      T2=SI(NR2)+SI(NR6)
      SI(NR6)=SI(NR2)-SI(NR6)
      SI(NBASE)=T0+T2
      SI(NR2)=T0-T2
      SI(NR1)=T1+T3
      SI(NR3)=T1-T3
      SI(NR5)=T5+T7
      SI(NR7)=T5-T7
 820  NBASE=NBASE+8
 830  NBASE=NBASE+NLUP2
 840  NBASE=NBASE+NLUP23
 1600 IF(NA.NE.16) GO TO 300
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 16 POINT PRE-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=18*(ND2-NB)
      NLUP23=18*ND2*(ND3-NC)
      NBASE=1
      DO 1640 N4=1,ND
      DO 1630 N3=1,NC
      DO 1620 N2=1,NB
      NR1=NBASE+1
      NR2=NR1+1
      NR3=NR2+1
      NR4=NR3+1
      NR5=NR4+1
      NR6=NR5+1
      NR7=NR6+1
      NR8=NR7+1
      NR9=NR8+1
      NR10=NR9+1
      NR11=NR10+1
      NR12=NR11+1
      NR13=NR12+1
      NR14=NR13+1
      NR15=NR14+1
      NR16=NR15+1
      NR17=NR16+1
      JBASE=NBASE
      DO 1645 J=1,8
      T(J)=SR(JBASE)+SR(JBASE+8)
      T(J+8)=SR(JBASE)-SR(JBASE+8)
      JBASE=JBASE+1
 1645 CONTINUE
      DO 1650 J=1,4
      Q(J)=T(J)+T(J+4)
      Q(J+4)=T(J)-T(J+4)
 1650 CONTINUE
      SR(NBASE)=Q(1)+Q(3)
      SR(NR2)=Q(1)-Q(3)
      SR(NR1)=Q(2)+Q(4)
      SR(NR3)=Q(2)-Q(4)
      SR(NR5)=Q(6)+Q(8)
      SR(NR7)=Q(6)-Q(8)
      SR(NR4)=Q(5)
      SR(NR6)=Q(7)
      SR(NR8)=T(9)
      SR(NR9)=T(10)+T(16)
      SR(NR15)=T(10)-T(16)
      SR(NR13)=T(14)+T(12)
      SR(NR11)=T(14)-T(12)
      SR(NR17)=SR(NR11)+SR(NR15)
      SR(NR16)=SR(NR9)+SR(NR13)
      SR(NR10)=T(11)+T(15)
      SR(NR14)=T(11)-T(15)
      SR(NR12)=T(13)
      JBASE=NBASE
      DO 1745 J=1,8
      T(J)=SI(JBASE)+SI(JBASE+8)
      T(J+8)=SI(JBASE)-SI(JBASE+8)
      JBASE=JBASE+1
 1745 CONTINUE
      DO 1750 J=1,4
      Q(J)=T(J)+T(J+4)
      Q(J+4)=T(J)-T(J+4)
 1750 CONTINUE
      SI(NBASE)=Q(1)+Q(3)
      SI(NR2)=Q(1)-Q(3)
      SI(NR1)=Q(2)+Q(4)
      SI(NR3)=Q(2)-Q(4)
      SI(NR5)=Q(6)+Q(8)
      SI(NR7)=Q(6)-Q(8)
      SI(NR4)=Q(5)
      SI(NR6)=Q(7)
      SI(NR8)=T(9)
      SI(NR9)=T(10)+T(16)
      SI(NR15)=T(10)-T(16)
      SI(NR13)=T(14)+T(12)
      SI(NR11)=T(14)-T(12)
      SI(NR17)=SI(NR11)+SI(NR15)
      SI(NR16)=SI(NR9)+SI(NR13)
      SI(NR10)=T(11)+T(15)
      SI(NR14)=T(11)-T(15)
      SI(NR12)=T(13)
 1620 NBASE=NBASE+18
 1630 NBASE=NBASE+NLUP2
 1640 NBASE=NBASE+NLUP23
 300  IF(NB.EQ.1) GO TO 700
      IF(NB.NE.3) GO TO 900
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 3 POINT PRE-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=2*ND1
      NLUP23=3*ND1*(ND3-NC)
      NBASE=1
      NOFF=ND1
      DO 340 N4=1,ND
      DO 330 N3=1,NC
      DO 310 N2=1,ND1
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      T1=SR(NR1)+SR(NR2)
      SR(NBASE)=SR(NBASE)+T1
      SR(NR2)=SR(NR1)-SR(NR2)
      SR(NR1)=T1
      T1=SI(NR1)+SI(NR2)
      SI(NBASE)=SI(NBASE)+T1
      SI(NR2)=SI(NR1)-SI(NR2)
      SI(NR1)=T1
 310  NBASE=NBASE+1
 330  NBASE=NBASE+NLUP2
 340  NBASE=NBASE+NLUP23
 900  IF(NB.NE.9) GO TO 700
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 9 POINT PRE-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=10*ND1
      NLUP23=11*ND1*(ND3-NC)
      NBASE=1
      NOFF=ND1
      DO 940 N4=1,ND
      DO 930 N3=1,NC
      DO 910 N2=1,ND1
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      NR3=NR2+NOFF
      NR4=NR3+NOFF
      NR5=NR4+NOFF
      NR6=NR5+NOFF
      NR7=NR6+NOFF
      NR8=NR7+NOFF
      NR9=NR8+NOFF
      NR10=NR9+NOFF
      T3=SR(NR3)+SR(NR6)
      T6=SR(NR3)-SR(NR6)
      SR(NBASE)=SR(NBASE)+T3
      T7=SR(NR7)+SR(NR2)
      T2=SR(NR7)-SR(NR2)
      SR(NR2)=T6
      T1=SR(NR1)+SR(NR8)
      T8=SR(NR1)-SR(NR8)
      SR(NR1)=T3
      T4=SR(NR4)+SR(NR5)
      T5=SR(NR4)-SR(NR5)
      SR(NR3)=T1+T4+T7
      SR(NR4)=T1-T7
      SR(NR5)=T4-T1
      SR(NR6)=T7-T4
      SR(NR10)=T2+T5+T8
      SR(NR7)=T8-T2
      SR(NR8)=T5-T8
      SR(NR9)=T2-T5
      T3=SI(NR3)+SI(NR6)
      T6=SI(NR3)-SI(NR6)
      SI(NBASE)=SI(NBASE)+T3
      T7=SI(NR7)+SI(NR2)
      T2=SI(NR7)-SI(NR2)
      SI(NR2)=T6
      T1=SI(NR1)+SI(NR8)
      T8=SI(NR1)-SI(NR8)
      SI(NR1)=T3
      T4=SI(NR4)+SI(NR5)
      T5=SI(NR4)-SI(NR5)
      SI(NR3)=T1+T4+T7
      SI(NR4)=T1-T7
      SI(NR5)=T4-T1
      SI(NR6)=T7-T4
      SI(NR10)=T2+T5+T8
      SI(NR7)=T8-T2
      SI(NR8)=T5-T8
      SI(NR9)=T2-T5
 910  NBASE=NBASE+1
 930  NBASE=NBASE+NLUP2
 940  NBASE=NBASE+NLUP23
 700  IF(NC.NE.7) GO TO 500
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 7 POINT PRE-WEAVE MODULE
C
C **********************************************************************
C
      NOFF=ND1*ND2
      NBASE=1
      NLUP2=8*NOFF
      DO 740 N4=1,ND
      DO 710 N1=1,NOFF
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      NR3=NR2+NOFF
      NR4=NR3+NOFF
      NR5=NR4+NOFF
      NR6=NR5+NOFF
      NR7=NR6+NOFF
      NR8=NR7+NOFF
      T1=SR(NR1)+SR(NR6)
      T6=SR(NR1)-SR(NR6)
      T4=SR(NR4)+SR(NR3)
      T3=SR(NR4)-SR(NR3)
      T2=SR(NR2)+SR(NR5)
      T5=SR(NR2)-SR(NR5)
      SR(NR5)=T6-T3
      SR(NR2)=T5+T3+T6
      SR(NR6)=T5-T6
      SR(NR8)=T3-T5
      SR(NR3)=T2-T1
      SR(NR4)=T1-T4
      SR(NR7)=T4-T2
      T1=T1+T4+T2
      SR(NBASE)=SR(NBASE)+T1
      SR(NR1)=T1
      T1=SI(NR1)+SI(NR6)
      T6=SI(NR1)-SI(NR6)
      T4=SI(NR4)+SI(NR3)
      T3=SI(NR4)-SI(NR3)
      T2=SI(NR2)+SI(NR5)
      T5=SI(NR2)-SI(NR5)
      SI(NR5)=T6-T3
      SI(NR2)=T5+T3+T6
      SI(NR6)=T5-T6
      SI(NR8)=T3-T5
      SI(NR3)=T2-T1
      SI(NR4)=T1-T4
      SI(NR7)=T4-T2
      T1=T1+T4+T2
      SI(NBASE)=SI(NBASE)+T1
      SI(NR1)=T1
 710  NBASE=NBASE+1
 740  NBASE=NBASE+NLUP2
 500  IF(ND.NE.5) RETURN
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 5 POINT PRE-WEAVE MODULE
C
C **********************************************************************
C
      NOFF=ND1*ND2*ND3
      NBASE=1
      DO 510 N1=1,NOFF
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      NR3=NR2+NOFF
      NR4=NR3+NOFF
      NR5=NR4+NOFF
      T4=SR(NR1)-SR(NR4)
      T1=SR(NR1)+SR(NR4)
      T3=SR(NR3)+SR(NR2)
      T2=SR(NR3)-SR(NR2)
      SR(NR3)=T1-T3
      SR(NR1)=T1+T3
      SR(NBASE)=SR(NBASE)+SR(NR1)
      SR(NR5)=T2+T4
      SR(NR2)=T4
      SR(NR4)=T2
      T4=SI(NR1)-SI(NR4)
      T1=SI(NR1)+SI(NR4)
      T3=SI(NR3)+SI(NR2)
      T2=SI(NR3)-SI(NR2)
      SI(NR3)=T1-T3
      SI(NR1)=T1+T3
      SI(NBASE)=SI(NBASE)+SI(NR1)
      SI(NR5)=T2+T4
      SI(NR2)=T4
      SI(NR4)=T2
 510  NBASE=NBASE+1
      RETURN
      END
C
C-----------------------------------------------------------------------
C
C SUBROUTINE: WEAVE2
C   THIS SUBROUTINE IMPLEMENTS THE POST-WEAVE MODULES
C   OF THE WFTA.  THE WORKING ARRAYS ARE SR AND SI.
C   THE ROUTINE CHECKS TO SEE WHICH FACTORS ARE PRESENT
C   IN THE TRANSFORM LENGTH N = NA*NB*NC*ND AND EXECUTES
C   THE POST-WEAVE CODE FOR THESE FACTORS.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE WEAVE2(SR,SI)
      COMMON NA,NB,NC,ND,ND1,ND2,ND3,ND4
      DIMENSION SR(1),SI(1)
      DIMENSION Q(8),T(16)
      IF(ND.NE.5) GO TO 700
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 5 POINT POST-WEAVE MODULE
C
C **********************************************************************
C
      NOFF=ND1*ND2*ND3
      NBASE=1
      DO 510 N1=1,NOFF
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      NR3=NR2+NOFF
      NR4=NR3+NOFF
      NR5=NR4+NOFF
      T1=SR(NBASE)+SR(NR1)
      T3=T1-SR(NR3)
      T1=T1+SR(NR3)
      T4=SI(NR2)+SI(NR5)
      T2=SI(NR4)+SI(NR5)
      SR(NR1)=T1-T4
      SR4=T1+T4
      SR2=T3+T2
      SR(NR3)=T3-T2
      T1=SI(NBASE)+SI(NR1)
      T3=T1-SI(NR3)
      T1=T1+SI(NR3)
      T4=SR(NR2)+SR(NR5)
      T2=SR(NR4)+SR(NR5)
      SI(NR1)=T1+T4
      SI(NR4)=T1-T4
      SI(NR2)=T3-T2
      SI(NR3)=T3+T2
      SR(NR2)=SR2
      SR(NR4)=SR4
 510  NBASE=NBASE+1
 700  IF(NC.NE.7) GO TO 300
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 7 POINT POST-WEAVE MODULE
C
C **********************************************************************
C
      NOFF=ND1*ND2
      NBASE=1
      NLUP2=8*NOFF
      DO 740 N4=1,ND
      DO 710 N1=1,NOFF
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      NR3=NR2+NOFF
      NR4=NR3+NOFF
      NR5=NR4+NOFF
      NR6=NR5+NOFF
      NR7=NR6+NOFF
      NR8=NR7+NOFF
      T1=SR(NR1)+SR(NBASE)
      T2=T1-SR(NR3)-SR(NR4)
      T4=T1+SR(NR3)-SR(NR7)
      T1=T1+SR(NR4)+SR(NR7)
      T6=SI(NR2)+SI(NR5)+SI(NR8)
      T5=SI(NR2)-SI(NR5)-SI(NR6)
      T3=SI(NR2)+SI(NR6)-SI(NR8)
      SR(NR1)=T1-T6
      SR6=T1+T6
      SR2=T2-T5
      SR5=T2+T5
      SR(NR4)=T4-T3
      SR(NR3)=T4+T3
      T1=SI(NR1)+SI(NBASE)
      T2=T1-SI(NR3)-SI(NR4)
      T4=T1+SI(NR3)-SI(NR7)
      T1=T1+SI(NR4)+SI(NR7)
      T6=SR(NR2)+SR(NR5)+SR(NR8)
      T5=SR(NR2)-SR(NR5)-SR(NR6)
      T3=SR(NR2)+SR(NR6)-SR(NR8)
      SI(NR1)=T1+T6
      SI(NR6)=T1-T6
      SI(NR2)=T2+T5
      SI(NR5)=T2-T5
      SI(NR4)=T4+T3
      SI(NR3)=T4-T3
      SR(NR2)=SR2
      SR(NR5)=SR5
      SR(NR6)=SR6
 710  NBASE=NBASE+1
 740  NBASE=NBASE+NLUP2
 300  IF(NB.EQ.1) GO TO 400
      IF(NB.NE.3) GO TO 900
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 3 POINT POST-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=2*ND1
      NLUP23=3*ND1*(ND3-NC)
      NBASE=1
      NOFF=ND1
      DO 340 N5=1,ND
      DO 330 N4=1,NC
      DO 310 N2=1,ND1
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      T1=SR(NBASE)+SR(NR1)
      SR(NR1)=T1-SI(NR2)
      SR2=T1+SI(NR2)
      T1=SI(NBASE)+SI(NR1)
      SI(NR1)=T1+SR(NR2)
      SI(NR2)=T1-SR(NR2)
      SR(NR2)=SR2
 310  NBASE=NBASE+1
 330  NBASE=NBASE+NLUP2
 340  NBASE=NBASE+NLUP23
 900  IF(NB.NE.9) GO TO 400
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 9 POINT POST-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=10*ND1
      NLUP23=11*ND1*(ND3-NC)
      NBASE=1
      NOFF=ND1
      DO 940 N4=1,ND
      DO 930 N3=1,NC
      DO 910 N2=1,ND1
      NR1=NBASE+NOFF
      NR2=NR1+NOFF
      NR3=NR2+NOFF
      NR4=NR3+NOFF
      NR5=NR4+NOFF
      NR6=NR5+NOFF
      NR7=NR6+NOFF
      NR8=NR7+NOFF
      NR9=NR8+NOFF
      NR10=NR9+NOFF
      T3=SR(NBASE)-SR(NR3)
      T7=SR(NBASE)+SR(NR1)
      SR(NBASE)=SR(NBASE)+SR(NR3)+SR(NR3)
      T6=T3+SI(NR10)
      SR(NR3)=T3-SI(NR10)
      T4=T7+SR(NR5)-SR(NR6)
      T1=T7-SR(NR4)-SR(NR5)
      T7=T7+SR(NR4)+SR(NR6)
      SR(NR6)=T6
      T8=SI(NR2)-SI(NR7)-SI(NR8)
      T5=SI(NR2)+SI(NR8)-SI(NR9)
      T2=SI(NR2)+SI(NR7)+SI(NR9)
      SR(NR1)=T7-T2
      SR8=T7+T2
      SR(NR4)=T1-T8
      SR(NR5)=T1+T8
      SR7=T4-T5
      SR2=T4+T5
      T3=SI(NBASE)-SI(NR3)
      T7=SI(NBASE)+SI(NR1)
      SI(NBASE)=SI(NBASE)+SI(NR3)+SI(NR3)
      T6=T3-SR(NR10)
      SI(NR3)=T3+SR(NR10)
      T4=T7+SI(NR5)-SI(NR6)
      T1=T7-SI(NR4)-SI(NR5)
      T7=T7+SI(NR4)+SI(NR6)
      SI(NR6)=T6
      T8=SR(NR2)-SR(NR7)-SR(NR8)
      T5=SR(NR2)+SR(NR8)-SR(NR9)
      T2=SR(NR2)+SR(NR7)+SR(NR9)
      SI(NR1)=T7+T2
      SI(NR8)=T7-T2
      SI(NR4)=T1+T8
      SI(NR5)=T1-T8
      SI(NR7)=T4+T5
      SI(NR2)=T4-T5
      SR(NR2)=SR2
      SR(NR7)=SR7
      SR(NR8)=SR8
 910  NBASE=NBASE+1
 930  NBASE=NBASE+NLUP2
 940  NBASE=NBASE+NLUP23
 400  IF(NA.EQ.1) RETURN
      IF(NA.NE.4) GO TO 800
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 4 POINT POST-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=4*(ND2-NB)
      NLUP23=4*ND2*(ND3-NC)
      NBASE=1
      DO 440 N4=1,ND
      DO 430 N3=1,NC
      DO 420 N2=1,NB
      NR1=NBASE+1
      NR2=NR1+1
      NR3=NR2+1
      TR0=SR(NBASE)+SR(NR2)
      TR2=SR(NBASE)-SR(NR2)
      TR1=SR(NR1)+SR(NR3)
      TR3=SR(NR1)-SR(NR3)
      TI1=SI(NR1)+SI(NR3)
      TI3=SI(NR1)-SI(NR3)
      SR(NBASE)=TR0+TR1
      SR(NR2)=TR0-TR1
      SR(NR1)=TR2+TI3
      SR(NR3)=TR2-TI3
      TI0=SI(NBASE)+SI(NR2)
      TI2=SI(NBASE)-SI(NR2)
      SI(NBASE)=TI0+TI1
      SI(NR2)=TI0-TI1
      SI(NR1)=TI2-TR3
      SI(NR3)=TI2+TR3
 420  NBASE=NBASE+4
 430  NBASE=NBASE+NLUP2
 440  NBASE=NBASE+NLUP23
 800  IF(NA.NE.8) GO TO 1600
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 8 POINT POST-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=8*(ND2-NB)
      NLUP23=8*ND2*(ND3-NC)
      NBASE=1
      DO 840 N4=1,ND
      DO 830 N3=1,NC
      DO 820 N2=1,NB
      NR1=NBASE+1
      NR2=NR1+1
      NR3=NR2+1
      NR4=NR3+1
      NR5=NR4+1
      NR6=NR5+1
      NR7=NR6+1
      T1=SR(NBASE)-SR(NR1)
      SR(NBASE)=SR(NBASE)+SR(NR1)
      SR6=SR(NR2)+SI(NR3)
      SR(NR2)=SR(NR2)-SI(NR3)
      T4=SR(NR4)-SI(NR5)
      T5=SR(NR4)+SI(NR5)
      T6=SR(NR7)-SI(NR6)
      T7=SR(NR7)+SI(NR6)
      SR(NR4)=T1
      SR(NR1)=T4+T6
      SR3=T4-T6
      SR5=T5-T7
      SR(NR7)=T5+T7
      T1=SI(NBASE)-SI(NR1)
      SI(NBASE)=SI(NBASE)+SI(NR1)
      T3=SI(NR2)-SR(NR3)
      SI(NR2)=SI(NR2)+SR(NR3)
      T4=SI(NR4)+SR(NR5)
      T5=SI(NR4)-SR(NR5)
      SI(NR6)=T3
      T6=SR(NR6)+SI(NR7)
      T7=SR(NR6)-SI(NR7)
      SI(NR4)=T1
      SI(NR1)=T4+T6
      SI(NR3)=T4-T6
      SI(NR5)=T5+T7
      SI(NR7)=T5-T7
      SR(NR3)=SR3
      SR(NR5)=SR5
      SR(NR6)=SR6
 820  NBASE=NBASE+8
 830  NBASE=NBASE+NLUP2
 840  NBASE=NBASE+NLUP23
 1600 IF(NA.NE.16) RETURN
C
C **********************************************************************
C
C     THE FOLLOWING CODE IMPLEMENTS THE 16 POINT POST-WEAVE MODULE
C
C **********************************************************************
C
      NLUP2=18*(ND2-NB)
      NLUP23=18*ND2*(ND3-NC)
      NBASE=1
      DO 1640 N4=1,ND
      DO 1630 N3=1,NC
      DO 1620 N2=1,NB
      NR1=NBASE+1
      NR2=NR1+1
      NR3=NR2+1
      NR4=NR3+1
      NR5=NR4+1
      NR6=NR5+1
      NR7=NR6+1
      NR8=NR7+1
      NR9=NR8+1
      NR10=NR9+1
      NR11=NR10+1
      NR12=NR11+1
      NR13=NR12+1
      NR14=NR13+1
      NR15=NR14+1
      NR16=NR15+1
      NR17=NR16+1
      T(2)=SR(NBASE)-SR(NR1)
      SR(NBASE)=SR(NR1)+SR(NBASE)
      T(4)=SR(NR2)+SI(NR3)
      T(3)=SR(NR2)-SI(NR3)
      T(6)=SR(NR4)+SI(NR5)
      T(5)=SR(NR4)-SI(NR5)
      T(8)=-SI(NR6)-SR(NR7)
      T(7)=-SI(NR6)+SR(NR7)
      T(9)=SR(NR8)+SR(NR14)
      T(15)=SR(NR8)-SR(NR14)
      T(13)=-SI(NR10)-SI(NR12)
      T(11)=SI(NR10)-SI(NR12)
      T(16)=SR(NR15)-SR(NR17)
      T(12)=SR(NR11)-SR(NR17)
      T(10)=-SI(NR9)-SI(NR16)
      T(14)=-SI(NR16)+SI(NR13)
      SR(NR2)=T(5)+T(7)
      SR6=T(5)-T(7)
      SR10=T(6)+T(8)
      SR(NR14)=T(6)-T(8)
      Q(7)=T(9)+T(10)
      Q(8)=T(9)-T(10)
      Q(1)=T(11)+T(12)
      Q(2)=T(11)-T(12)
      Q(4)=T(14)+T(15)
      Q(5)=T(15)-T(14)
      Q(3)=T(13)+T(16)
      Q(6)=T(13)-T(16)
      SR(NR1)=Q(3)+Q(7)
      SR(NR7)=Q(7)-Q(3)
      SR9=Q(6)+Q(8)
      SR(NR15)=Q(8)-Q(6)
      SR5=Q(1)+Q(4)
      SR3=Q(4)-Q(1)
      SR13=Q(2)+Q(5)
      SR11=Q(5)-Q(2)
      SR(NR8)=T(2)
      SR(NR4)=T(3)
      SR12=T(4)
      T(2)=SI(NBASE)-SI(NR1)
      SI(NBASE)=SI(NR1)+SI(NBASE)
      T(4)=SI(NR2)-SR(NR3)
      T(3)=SI(NR2)+SR(NR3)
      T(6)=SI(NR4)-SR(NR5)
      T(5)=SI(NR4)+SR(NR5)
      T(8)=SR(NR6)-SI(NR7)
      T(7)=SR(NR6)+SI(NR7)
      T(9)=SI(NR8)+SI(NR14)
      T(15)=SI(NR8)-SI(NR14)
      T(13)=SR(NR10)+SR(NR12)
      T(11)=SR(NR12)-SR(NR10)
      T(16)=SI(NR15)-SI(NR17)
      T(12)=SI(NR11)-SI(NR17)
      T(10)=SR(NR9)+SR(NR16)
      T(14)=SR(NR16)-SR(NR13)
      SI(NR2)=T(5)+T(7)
      SI(NR6)=T(5)-T(7)
      SI(NR10)=T(6)+T(8)
      SI(NR14)=T(6)-T(8)
      Q(7)=T(9)+T(10)
      Q(8)=T(9)-T(10)
      Q(1)=T(11)+T(12)
      Q(2)=T(11)-T(12)
      Q(4)=T(14)+T(15)
      Q(5)=T(15)-T(14)
      Q(3)=T(13)+T(16)
      Q(6)=T(13)-T(16)
      SI(NR1)=Q(3)+Q(7)
      SI(NR7)=Q(7)-Q(3)
      SI(NR9)=Q(6)+Q(8)
      SI(NR15)=Q(8)-Q(6)
      SI(NR5)=Q(1)+Q(4)
      SI(NR3)=Q(4)-Q(1)
      SI(NR13)=Q(2)+Q(5)
      SI(NR11)=Q(5)-Q(2)
      SI(NR8)=T(2)
      SI(NR4)=T(3)
      SI(NR12)=T(4)
      SR(NR3)=SR3
      SR(NR5)=SR5
      SR(NR6)=SR6
      SR(NR9)=SR9
      SR(NR10)=SR10
      SR(NR11)=SR11
      SR(NR12)=SR12
      SR(NR13)=SR13
 1620 NBASE=NBASE+18
 1630 NBASE=NBASE+NLUP2
 1640 NBASE=NBASE+NLUP23
      RETURN
      END
