C
C      ________________________________________________________
C     |                                                        |
C     |          THE CONJUGATE FAST FOURIER TRANSFORM          |
C     |             NEW A SUB I = SUM FROM J=1 TO N            |
C     |    EXP(-2*PI*SQRT(-1)/N)**((I-1)*(J-1)) OLD A SUB J    |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --COEFFICIENTS                           |
C     |                                                        |
C     |         N     --NUMBER OF COEFFICIENTS                 |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         A     --TRANSFORMED COEFFICIENTS               |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ACOS,CEXP,CMPLX                  |
C     |________________________________________________________|
C
      SUBROUTINE FFC(A,N,W)
      COMPLEX A(1),W(1),S,T,U,V
      INTEGER D,E,F,G,H,I,J,K,L,M,N,NP,O,P(25)
      DATA P/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
     1 73,79,83,89,97/
      DATA NP/25/
      M = N
      F = 0
      U = CMPLX(0.,-2.*ACOS(-1.)/N)
C     --------------------------------------------------
C     |*** SOME COMPILERS USE ARCOS INSTEAD OF ACOS ***|
C     --------------------------------------------------
10    IF ( M .EQ. 1 ) GOTO 900
      DO 20 I = 1,NP
20         IF ( (M/P(I))*P(I) .EQ. M ) GOTO 30
      L = M
      GOTO 40
30    L = P(I)
40    O = M
      M = M/L
      V = CEXP(M*U)
      S = (1.,0.)
      H = 0
      IF ( F .EQ. 1 ) GOTO 470
      IF ( L .EQ. 2 ) GOTO 50
      IF ( L .EQ. 3 ) GOTO 170
      GOTO 290
50    GOTO (150,130,110,90),M
60    J = -H
70    I = H + 1
      H = H + M
      E = J + M
      DO 80 K = I,H
80         W(K) = A(J+K) + S*A(E+K)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 70
      IF ( H .LT. N ) GOTO 60
      F = 1
      GOTO 10
90    J = -H
100   H = H + 1
      E = J + M
      W(H) = A(J+H) + S*A(E+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 100
      IF ( H .LT. N ) GOTO 90
      F = 1
      GOTO 10
110   J = -H
120   H = H + 1
      E = J + M
      W(H) = A(J+H) + S*A(E+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 120
      IF ( H .LT. N ) GOTO 110
      F = 1
      GOTO 10
130   J = -H
140   H = H + 1
      E = J + M
      W(H) = A(J+H) + S*A(E+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 140
      IF ( H .LT. N ) GOTO 130
      F = 1
      GOTO 10
150   J = -H
160   H = H + 1
      E = J + M
      W(H) = A(J+H) + S*A(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 160
      IF ( H .LT. N ) GOTO 150
      F = 1
      GOTO 10
170   GOTO (270,250,230,210),M
180   J = -H
190   I = H + 1
      H = H + M
      E = J + M
      D = E + M
      T = S*S
      DO 200 K = I,H
200        W(K) = A(J+K) + S*A(E+K) + T*A(D+K)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 190
      IF ( H .LT. N ) GOTO 180
      F = 1
      GOTO 10
210   J = -H
220   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 220
      IF ( H .LT. N ) GOTO 210
      F = 1
      GOTO 10
230   J = -H
240   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 240
      IF ( H .LT. N ) GOTO 230
      F = 1
      GOTO 10
250   J = -H
260   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      H = H + 1
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 260
      IF ( H .LT. N ) GOTO 250
      F = 1
      GOTO 10
270   J = -H
280   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      W(H) = A(J+H) + S*A(E+H) + T*A(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 280
      IF ( H .LT. N ) GOTO 270
      F = 1
      GOTO 10
290   GOTO (440,410,380,350),M
300   J = -H
310   T = (1.,0.)
      I = H + 1
      H = H + M
      G = J + O
      DO 320 K = I,H
320        W(K) =  (0.,0.)
330   DO 340 K = I,H
340        W(K) = W(K) + T*A(J+K)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 330
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 310
      IF ( H .LT. N ) GOTO 300
      F = 1
      GOTO 10
350   J = -H
360   T = (1.,0.)
      I = H + 1
      E = I + 1
      D = E + 1
      H = H + M
      G = J + O
      W(I) = (0.,0.)
      W(E) = (0.,0.)
      W(D) = (0.,0.)
      W(H) = (0.,0.)
370   W(I) = W(I) + T*A(J+I)
      W(E) = W(E) + T*A(J+E)
      W(D) = W(D) + T*A(J+D)
      W(H) = W(H) + T*A(J+H)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 370
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 360
      IF ( H .LT. N ) GOTO 350
      F = 1
      GOTO 10
380   J = -H
390   T = (1.,0.)
      I = H + 1
      E = I + 1
      H = H + M
      G = J + O
      W(I) = (0.,0.)
      W(E) = (0.,0.)
      W(H) = (0.,0.)
400   W(I) = W(I) + T*A(J+I)
      W(E) = W(E) + T*A(J+E)
      W(H) = W(H) + T*A(J+H)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 400
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 390
      IF ( H .LT. N ) GOTO 380
      F = 1
      GOTO 10
410   J = -H
420   T = (1.,0.)
      I = H + 1
      H = H + M
      G = J + O
      W(I) = (0.,0.)
      W(H) = (0.,0.)
430   W(I) = W(I) + T*A(J+I)
      W(H) = W(H) + T*A(J+H)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 430
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 420
      IF ( H .LT. N ) GOTO 410
      F = 1
      GOTO 10
440   J = -H
450   T = (1.,0.)
      I = H + 1
      H = H + M
      G = J + O
      W(I) = (0.,0.)
460   W(I) = W(I) + T*A(J+I)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 460
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 450
      IF ( H .LT. N ) GOTO 440
      F = 1
      GOTO 10
470   IF ( L .EQ. 2 ) GOTO 480
      IF ( L .EQ. 3 ) GOTO 600
      GOTO 720
480   GOTO (580,560,540,520),M
490   J = -H
500   I = H + 1
      H = H + M
      E = J + M
      DO 510 K = I,H
510        A(K) = W(J+K) + S*W(E+K)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 500
      IF ( H .LT. N ) GOTO 490
      F = 0
      GOTO 10
520   J = -H
530   H = H + 1
      E = J + M
      A(H) = W(J+H) + S*W(E+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 530
      IF ( H .LT. N ) GOTO 520
      F = 0
      GOTO 10
540   J = -H
550   H = H + 1
      E = J + M
      A(H) = W(J+H) + S*W(E+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 550
      IF ( H .LT. N ) GOTO 540
      F = 0
      GOTO 10
560   J = -H
570   H = H + 1
      E = J + M
      A(H) = W(J+H) + S*W(E+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 570
      IF ( H .LT. N ) GOTO 560
      F = 0
      GOTO 10
580   J = -H
590   H = H + 1
      E = J + M
      A(H) = W(J+H) + S*W(E+H)
      J = E
      S = S*V
      IF ( J+H .LT. N ) GOTO 590
      IF ( H .LT. N ) GOTO 580
      F = 0
      GOTO 10
600   GOTO (700,680,660,640),M
610   J = -H
620   I = H + 1
      H = H + M
      E = J + M
      D = E + M
      T = S*S
      DO 630 K = I,H
630        A(K) = W(J+K) + S*W(E+K) + T*W(D+K)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 620
      IF ( H .LT. N ) GOTO 610
      F = 0
      GOTO 10
640   J = -H
650   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 650
      IF ( H .LT. N ) GOTO 640
      F = 0
      GOTO 10
660   J = -H
670   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 670
      IF ( H .LT. N ) GOTO 660
      F = 0
      GOTO 10
680   J = -H
690   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      H = H + 1
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 690
      IF ( H .LT. N ) GOTO 680
      F = 0
      GOTO 10
700   J = -H
710   H = H + 1
      E = J + M
      D = E + M
      T = S*S
      A(H) = W(J+H) + S*W(E+H) + T*W(D+H)
      J = D
      S = S*V
      IF ( J+H .LT. N ) GOTO 710
      IF ( H .LT. N ) GOTO 700
      F = 0
      GOTO 10
720   GOTO (870,840,810,780),M
730   J = -H
740   T = (1.,0.)
      I = H + 1
      H = H + M
      G = J + O
      DO 750 K = I,H
750        A(K) =  (0.,0.)
760   DO 770 K = I,H
770        A(K) = A(K) + T*W(J+K)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 760
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 740
      IF ( H .LT. N ) GOTO 730
      F = 0
      GOTO 10
780   J = -H
790   T = (1.,0.)
      I = H + 1
      E = I + 1
      D = E + 1
      H = H + M
      G = J + O
      A(I) = (0.,0.)
      A(E) = (0.,0.)
      A(D) = (0.,0.)
      A(H) = (0.,0.)
800   A(I) = A(I) + T*W(J+I)
      A(E) = A(E) + T*W(J+E)
      A(D) = A(D) + T*W(J+D)
      A(H) = A(H) + T*W(J+H)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 800
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 790
      IF ( H .LT. N ) GOTO 780
      F = 0
      GOTO 10
810   J = -H
820   T = (1.,0.)
      I = H + 1
      E = I + 1
      H = H + M
      G = J + O
      A(I) = (0.,0.)
      A(E) = (0.,0.)
      A(H) = (0.,0.)
830   A(I) = A(I) + T*W(J+I)
      A(E) = A(E) + T*W(J+E)
      A(H) = A(H) + T*W(J+H)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 830
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 820
      IF ( H .LT. N ) GOTO 810
      F = 0
      GOTO 10
840   J = -H
850   T = (1.,0.)
      I = H + 1
      H = H + M
      G = J + O
      A(I) = (0.,0.)
      A(H) = (0.,0.)
860   A(I) = A(I) + T*W(J+I)
      A(H) = A(H) + T*W(J+H)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 860
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 850
      IF ( H .LT. N ) GOTO 840
      F = 0
      GOTO 10
870   J = -H
880   T = (1.,0.)
      I = H + 1
      H = H + M
      G = J + O
      A(I) = (0.,0.)
890   A(I) = A(I) + T*W(J+I)
      T = T*S
      J = J + M
      IF ( J .LT. G ) GOTO 890
      J = J - M
      S = S*V
      IF ( J+H .LT. N ) GOTO 880
      IF ( H .LT. N ) GOTO 870
      F = 0
      GOTO 10
900   IF ( F .EQ. 0 ) RETURN
      DO 910 I = 1,N
910        A(I) = W(I)
      RETURN
      END
