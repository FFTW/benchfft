      SUBROUTINE EZFFTB (N,R,AZERO,A,B,WSAVE)
      DIMENSION       R(1)       ,A(1)       ,B(1)       ,WSAVE(1)
      IF (N-2) 101,102,103
  101 R(1) = AZERO
      RETURN
  102 R(1) = AZERO+A(1)
      R(2) = AZERO-A(1)
      RETURN
  103 NS2 = (N-1)/2
      DO 104 I=1,NS2
         R(2*I) = .5*A(I)
         R(2*I+1) = -.5*B(I)
  104 CONTINUE
      R(1) = AZERO
      IF (MOD(N,2) .EQ. 0) R(N) = A(NS2+1)
      CALL RFFTB (N,R,WSAVE(N+1))
      RETURN
      END
