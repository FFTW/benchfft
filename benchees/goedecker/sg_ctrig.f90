!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_ctrig
!! NAME
!! sg_ctrig
!!
!! FUNCTION
!! Precalculates trigonometric expressions and bitreversal key IND
!! (Stefan Goedecker lib)
!!
!! COPYRIGHT
!! Copyright (C) Stefan Goedecker, Ithaca, NY USA, July 14, 1993
!! Copyright (C) 1998-2001 ABINIT group (SG, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! ind(mg)
!! trig(2,mg) 
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This version of sg_ctrig produces cos and tan instead of sin and cos--
!! this allows for much greater efficiency on the superscalar architecture
!! of ibm rs6000 where floating point multiply and add (FMA) is used.
!!
!! TODO
!! Should describe arguments, and order them
!! Should suppress one-letter variables
!!
!! PARENTS
!!      fourdp,sg_fft,sg_fftpad,sg_fftrisc
!!
!! SOURCE

 subroutine sg_ctrig(n,trig,aft,bef,now,ris,ic,ind,mfac,mg)
 
 use basis_defs
 implicit none

!Arguments ------------------------------------
 integer :: ic,mfac,mg,n
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: ris
 real(dp) :: trig(2,mg)

!Local variables-------------------------------
 integer, save :: nextmx=4
 integer :: i,ii,inc,irep,l,j,k,next,nh
!"prime" is the set of radices coded elsewhere for fft
 integer,save :: prime(4)=(/5,4,3,2/)
 real(dp) :: angle,trigc,trigs,twopi
 character :: message*500

! *************************************************************************

!**Note**
!2*Pi must not be defined too accurately here or else
!cos(twopi/2) will be exactly 0 and sin/cos below will be
!infinite; if a small error is left in Pi, then sin/cos will
!be about 10**14 and later cos * (sin/cos) will be 1 to within
!about 10**(-14) and the fft routines will work
!The precision on sgi causes the algorithm to fail if 
!twopi is defined as 8.d0*atan(1.0d0).

 twopi=6.2831853071795867d0

 angle=ris*twopi/n 
!trig(1,0)=1.d0
!trig(2,0)=0.d0
 if (mod(n,2)==0) then
  nh=n/2
  trig(1,nh)=-1.d0
  trig(2,nh)=0.d0
  do i=1,nh-1
   trigc=cos(i*angle)
   trigs=sin(i*angle)
   trig(1,i)=trigc
   trig(2,i)=trigs/trigc
   trig(1,n-i)=trigc
   trig(2,n-i)=-trigs/trigc
  enddo
 else
  nh=(n-1)/2
  do i=1,nh
   trigc=cos(i*angle)
   trigs=sin(i*angle)
   trig(1,i)=trigc
   trig(2,i)=trigs/trigc
   trig(1,n-i)=trigc
   trig(2,n-i)=-trigs/trigc
  enddo
 endif

 ic=1
 aft(ic)=1
 bef(ic)=n
 next=1
 
!An infinite loop, with exit or cycle instructions
 do
  if( (bef(ic)/prime(next))*prime(next)<bef(ic) ) then
   next=next+1
   if (next<=nextmx) then
    cycle
   else
    now(ic)=bef(ic)
    bef(ic)=1
   endif
  else
   now(ic)=prime(next)
   bef(ic)=bef(ic)/prime(next)
  endif
  aft(ic+1)=aft(ic)
  now(ic+1)=now(ic)
  bef(ic+1)=bef(ic)
  ic=ic+1
  if (ic>mfac) then
   write(message, '(a,a,a,a,i10,a,i5,a)' ) ch10, &
&   ' sg_ctrig: BUG -',ch10,&
&   '  number of factors ic=',ic,ch10,&
&   '  exceeds dimensioned mfac=',mfac,ch10
!   call wrtout(06,message,'PERS')
!   call leave_new('PERS')
  endif
  if (bef(ic)/=1) then
   aft(ic)=aft(ic)*now(ic)
   cycle
  endif 
! If not cycled, exit
  exit
 enddo

 ic=ic-1
 
!DEBUG
!print*,'now',(now(i),i=1,ic)
!print*,'aft',(aft(i),i=1,ic)
!print*,'bef',(bef(i),i=1,ic)
!ENDDEBUG

 do i=1,n
  ind(i)=1
 enddo

 irep=1
 inc=n
 do l=ic,1,-1
  inc=inc/now(l)
  ii=0
  do k=1,1+(n-1)/(now(l)*irep)
   do j=0,now(l)-1
    do i=1,irep
     ii=ii+1
     ind(ii)=ind(ii)+j*inc
    enddo
   enddo
  enddo
  irep=irep*now(l)
 enddo

 if (irep/=n) then
  write(message, '(a,a,a,a,i10,a,i10)' ) ch10,&
&  ' sg_ctrig : BUG -',ch10,&
&  '  irep should equal n ; irep=',irep,' n=',n
!  call wrtout(06,message,'PERS')
!  call leave_new('PERS')
 endif

 if (inc/=1) then
  write(message, '(a,a,a,a,i10)' ) ch10,&
&  ' sg_ctrig : BUG -',ch10,&
&  '  inc should equal 1 in sg_ctrig; inc=',inc
!  call wrtout(06,message,'PERS')
!  call leave_new('PERS')
 endif

 end subroutine
!!***
