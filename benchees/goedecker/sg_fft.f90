!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fft
!! NAME
!! sg_fft
!!
!! FUNCTION
!! Calculates the discrete fourier transform
!! ftarr(i1,i2,i3)=exp(ris*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) arr(j1,j2,j3)
!!
!! COPYRIGHT
!! Copyright by Stefan Goedecker, Ithaca, NY USA, July 14, 1993
!! Copyright (C) 1998-2001 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  arr(2,nd1,nd2,nd3)=input complex array with alternating real and imaginary
!!   elements; data resides in 2*n1*n2*n3 of this array, spread out.
!!  fftcache=size of the cache (kB)
!!  nd1,nd2,nd3=memory dimension of arr and ftarr
!!  n1,n2,n3=physical dimension of the transform
!!  ris=(real(dp)) sign of exponential in transform
!!
!! OUTPUT
!!  ftarr(2,nd1,nd2,nd3)=working space for transform and contains output
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  ndi must always be greater or equal to ni.  Recommended choice for nd1
!!  and nd2 is: ni for ni=odd or ni+1 for ni=even (hence 2*(ni/2)+1);
!!  nd3 should always be n3.  Note that choosing nd1 or nd2 larger than
!!  the recommended value can severely degrade efficiency of this routine.
!!  Avoiding even ndi for nd1 and nd2 avoids cache conflicts on cache machines.
!!  Each of n1,n2,n3 must be a
!!  product of the prime factors 2,3,5. If two ni s are equal
!!  it is recommended to place them behind each other.
!!  The largest any of these may be is set by parameter "mg" below.
!!  This fft is particularly efficient for cache architectures.
!!  Note that the meaning of fftcache has changed from the original
!!  ncache of SG (that was the maximum number of COMPLEX*16 in the cache)
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      ccfft
!!
!! SOURCE

 subroutine sg_fft(fftcache,nd1,nd2,nd3,n1,n2,n3,arr,ftarr,ris) 
 
 use basis_defs
 implicit none 

!Arguments ------------------------------------
 integer :: fftcache,nd1,nd2,nd3,n1,n2,n3
 real(dp) :: ris
 real(dp) :: arr(2,nd1,nd2,nd3),ftarr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!mfac sets maximum number of factors (5, 4, 3, or 2) which may be 
!contained within any n1, n2, or n3
!mg sets the maximum 1 dimensional fft length (any one of n1, n2, or n3)
 integer, parameter :: mfac=11,mg=2048
 integer :: ic,i2,nd13,nd23,n1i,n12,n2i,n23,n3i
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: trig(2,mg)
 character :: message*500

! *************************************************************************
 
!Check that dimension is not exceeded
 if (n1>mg.or.n2>mg.or.n3>mg) then
  write(message, '(a,a,a,a,3i10,a,i10,a)' ) ch10,&
&  ' sg_fft : BUG -',ch10,&
&  '  one of the dimensions n1,n2,n3=',n1,n2,n3,&
&  '  exceeds allowed dimension mg=',mg,ch10
!  call wrtout(06,message,'PERS')
!  call leave_new('PERS')
 endif
 
!transform along x direction
 call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 call sg_fftx(fftcache,mfac,mg,nd1,nd2,nd3,n2,n3,&
& arr,ftarr,trig,aft,now,bef,ris,ind,ic)
 
!transform along y direction
 if (n2/=n1)then
  call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 endif
 n1i=1 ; n3i=1
 call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,&
& ftarr,arr,trig,aft,now,bef,ris,ind,ic)
 
!transform along z direction
 if (n3/=n2)then
  call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 endif

!$OMP PARALLEL DO SHARED(aft,arr,bef,fftcache,ftarr,ind,ic)&
!$OMP&SHARED(nd1,nd2,nd3,now,n1,n2,ris,trig)&
!$OMP&PRIVATE(i2)
  do i2=1,n2
   call sg_fftz(fftcache,mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&   trig,aft,now,bef,ris,ind,ic)
  enddo
!$OMP END PARALLEL DO

 end subroutine
!!***
