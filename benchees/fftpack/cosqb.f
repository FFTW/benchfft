      SUBROUTINE COSQB (N,X,WSAVE)
      DIMENSION       X(1)       ,WSAVE(1)
      DATA TSQRT2 /2.82842712474619/
      IF (N-2) 101,102,103
  101 X(1) = 4.*X(1)
      RETURN
  102 X1 = 4.*(X(1)+X(2))
      X(2) = TSQRT2*(X(1)-X(2))
      X(1) = X1
      RETURN
  103 CALL COSQB1 (N,X,WSAVE,WSAVE(N+1))
      RETURN
      END
