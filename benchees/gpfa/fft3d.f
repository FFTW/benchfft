C
C   GPFFT3 Replacement for FFT3D for CASTEP
C
C  Calls Generalized Prime Factor Algorithm (Temperton, 1992).
C  Workspace array W is for compatibility and is not used.
C  
C  For Vector machines with interleaved memory, a stride of a
C  multiple of 2 can give memory-bank conflicts which degrage
C  performance by a factor of 2,3 or more.  Hence data of
C  dimension C(l,mln)  is "expanded" into an array C(l+1,m,n).
C  The expansion is done in-place, so C MUST BE DIMENSIONED
C  C(L+1,M,N) in calling prog.
C
      subroutine FFT3D(c,w,nn,is)
      dimension c(0:*),w(0:*)
      integer nn(3)
      logical lvecmc
      parameter (lvecmc = .FALSE. )

      if( lvecmc .AND. MOD(nn(1),4) .eq. 0 ) then
         call EXPAND(c,nn(1),nn(2)*nn(3),1)
         call GPF3D(c,nn(1)+1,nn,is)
         call EXPAND(c,nn(1),nn(2)*nn(3),-1)
      else
         call GPF3D(c,nn(1),nn,is)
      endif
      return
      end
