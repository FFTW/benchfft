!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fftpad
!! NAME
!! sg_fftpad
!!
!! FUNCTION
!! Fast fourier transform. 
!! This is the zero-padding version of "fft".  See fft for comments.
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
!!  gbound(4*mgbound+4)=integer ranges of reciprocal space points on which data
!!   is nonzero, as needed to skip fft of arrays of 0 s
!!  mgbound=needed to dimension gbound
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
!!  See fft.f
!!
!! TODO
!!
!! PARENTS
!!      fourwf
!!
!! SOURCE

 subroutine sg_fftpad(fftcache,mgbound,nd1,nd2,nd3,n1,n2,n3,&
& arr,ftarr,ris,gbound)
 
 use basis_defs
 implicit none

!Arguments ------------------------------------
 integer :: fftcache,mgbound,nd1,nd2,nd3,n1,n2,n3
 integer :: gbound(2*mgbound+4)
 real(dp) :: ris
 real(dp) :: arr(2,nd1,nd2,nd3),ftarr(2,nd1,nd2,nd3)

!Local variables-------------------------------
!mfac sets maximum number of factors (5, 4, 3, or 2) which may be
!contained within any n1, n2, or n3
!mg sets the maximum 1 dimensional fft length (any one of n1, n2, or n3)
!xg : the signification of mg is changed with respect to fft3dp !!!
 integer, parameter :: mfac=11,mg=2048
 integer :: i1,i2,i3,g3min,g3max,ic,index,jx,jy,jyy,jzz,jz,nd13,&
& n1i,n2i,n23,n3i,n3p
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp), parameter :: tol=1.0d-12
 real(dp) :: trig(2,mg)
 character :: message*500

! *************************************************************************

! DEBUG
!write(6,*)' sg_fftpad : enter '
!if(.true.) stop
!return
! ENDDEBUG
 
!Check that dimension is not exceeded
 if (n1>mg.or.n2>mg.or.n3>mg) then
  write(message, '(a,a,a,a,3i10,a,i10,a)' ) ch10,&
&  ' fft: BUG -',ch10,&
&  '  one of the dimensions n1,n2,n3=',n1,n2,n3,&
&  '  exceeds allowed dimension mg=',mg,ch10
!  call wrtout(06,message,'PERS')
!  call leave_new('PERS')
 endif

 g3min=gbound(1)
 g3max=gbound(2)

! --------------------------------------------------------------------------

 if (abs(ris-1.d0)<tol) then
 
! Handle G -> r  transform (G sphere to fft box)
  
! Transform along x direction
  call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)

! Zero out the untransformed (0) data part of the work array
! -- at every (y,z) there are 0 s to be added to the ends of
! the x data so have to zero whole thing.
  ftarr(:,:,:,:)=0.0d0

  call sg_fftpx(fftcache,mfac,mg,mgbound,nd1,nd2,nd3,n2,n3,&
&  arr,ftarr,trig,aft,now,bef,ris,ind,ic,gbound)
 
! Transform along y direction in two regions of z
  if (n2/=n1)then
   call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
  endif
 
! First y transform: z=1..g3max+1
  n3p=g3max+1
  n1i=1 ; n3i=1
  call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3p,ftarr,arr,&
&  trig,aft,now,bef,ris,ind,ic)
 
! Zero out the untransformed (0) data part of the work array
! -- only need to zero specified ranges of z
  arr(:,:,:,n3p+1:g3min+n3)=0.0d0
 
! Second y transform: z=g3min+1..0 (wrapped around)
  n3p=-g3min
  if (n3p>0) then
   n3i=1+g3min+n3 ; n1i=1
   call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,ftarr,arr,&
&   trig,aft,now,bef,ris,ind,ic)
  endif
 
! Transform along z direction
  if (n3/=n2) then
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

 else
 
!*************************************************
 
! Handle r -> G transform (from fft box to G sphere)
 
! Transform along z direction
  call sg_ctrig(n3,trig,aft,bef,now,ris,ic,ind,mfac,mg)

!$OMP PARALLEL DO SHARED(aft,arr,bef,fftcache,ftarr,ind,ic)&
!$OMP&SHARED(nd1,nd2,nd3,now,n1,n2,ris,trig)&
!$OMP&PRIVATE(i2)
  do i2=1,n2
   call sg_fftz(fftcache,mfac,mg,nd1,nd2,nd3,n1,i2,i2,arr,ftarr,&
&   trig,aft,now,bef,ris,ind,ic)
  enddo
!$OMP END PARALLEL DO

! Transform along y direction in two regions of z
  if (n2/=n3) then
   call sg_ctrig(n2,trig,aft,bef,now,ris,ic,ind,mfac,mg)
  endif
 
! First y transform: z=1..g3max+1
  n3p=g3max+1
  n1i=1 ; n3i=1
  call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3p,ftarr,arr,&
&  trig,aft,now,bef,ris,ind,ic)
 
! Second y transform: z=g3min+1..0 (wrapped around)
  n3p=-g3min
  if (n3p>0) then
   n1i=1 ; n3i=1+g3min+n3
   call sg_ffty(fftcache,mfac,mg,nd1,nd2,nd3,n1i,n1,n3i,n3,ftarr,arr,&
&   trig,aft,now,bef,ris,ind,ic)
  endif

! Transform along x direction
  if (n1/=n2) then 
   call sg_ctrig(n1,trig,aft,bef,now,ris,ic,ind,mfac,mg)
  endif
 
! Zero out the untransformed (0) data part of the work array
! -- at every (y,z) there are 0 s to be added to the ends of
! the x data so have to zero whole thing.
 ftarr(:,:,:,:)=0.0d0

  call sg_fftpx(fftcache,mfac,mg,mgbound,nd1,nd2,nd3,n2,n3,&
&       arr,ftarr,trig,aft,now,bef,ris,ind,ic,gbound)
 
! Data is now ready to be extracted from fft box to sphere

 endif
 
! DEBUG
!write(6,*)' sg_fftpad : enter '
!write(6,*)allocated(arr)
!if(.true.) stop
! ENDDEBUG

 end subroutine
!!***
