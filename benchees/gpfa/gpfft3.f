C
C GPF3D The basic 3D complex FFT.
C

C
C Arguments
C
C COMPLEX C(ID,NN(2),NN(3))
C ID first dimension of data array.  
C NN(3) Dimensions of 3D complex data.
C IS Forward (+1)/Reverse (-1)
C
      subroutine GPF3D(c,id,nn,is)
      parameter (NMAX=256)
      dimension c(*), nn(3), onn(3), trigs(NMAX,3)
      integer onn
      save onn,trigs
      data onn/0,0,0/,nthreads/1/


      do i=1,3
         if( onn(i) .ne. nn(i) ) then
            onn(i) = nn(i)
            call SETGPFA(trigs(1,i), nn(i))
         end if
      end do

C
C  The following three loops are prime candidates for parallelism
C  as they call independent multi-1D ffts.  Just set the variable
C  "nthreads" to the number of processors to use.
C
      lot=(nn(2)*nn(3)+nthreads-1)/nthreads
      do i=0,nthreads-1
         iofset=2*id*lot*i
         call GPFA(c(iofset+1),c(iofset+2),trigs(1,1),
     $        2,2*id,nn(1),MIN(lot,nn(2)*nn(3)-i*lot),-is)
      end do

      do i=0,nn(3)-1
         iofset=2*id*nn(2)*i
         call GPFA(c(iofset+1),c(iofset+2),trigs(1,2),
     $             2*id,2,nn(2),nn(1),-is)
      end do

      lot=(id*nn(2)+nthreads-1)/nthreads
      do i=0,nthreads-1
         iofset=2*lot*i
         call GPFA(c(iofset+1),c(iofset+2),trigs(1,3),
     $             2*id*nn(2),2,nn(3),MIN(lot,id*nn(2)-i*lot),-is)
      end do

      return
      end

      subroutine EXPAND(c,nl,nmn,is)
C
C  Stretch data stored in an array with a leading dimension of
C  nl=2**n into storage with a leading dimension of 2**n + 1.
C  This prevents memory-bank conflicts and therefore dramatically
C  improves performance on (vector) machines with interleaved memory.
C  Arrays are complex, so inital data of size 2*nl*nm*nn is increased
C  to 2*(nl+1)*nm*nn words.  Ensure array is dimensioned big enough
C  to take it!
C  In-place stretch can not be done in parallel as order of data move
C  is critical.  Vectorization does preserve order, if stride of -1
C  really does shift data top first.  May not be true for all vector
C  machines, so beware.
      dimension c(0:*)
      
      if( is .gt. 0 ) then
         do io=2*nmn-2,2,-2
            do ii=2*nl-1,0,-1
               c((ii+io*nl)+io) = c(ii+io*nl)
            end do
         end do
      else
         do io=2,2*nmn-2,2
            do ii=0,2*nl-1
               c(ii+io*nl) = c((ii+io*nl)+io)
            end do
         end do
      end if
      end

