!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fftz
!! NAME
!! sg_fftz
!!
!! FUNCTION
!! This subroutine is called by the 3-dimensional fft to conduct the
!! "z" transforms for all x and y.
!!
!! COPYRIGHT
!! Copyright by Stefan Goedecker, Ithaca, NY USA, July 14, 1993
!! Copyright (C) 1998-2001 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fftcache=size of the cache (kB)
!!  mfac = maximum number of factors in 1D FFTs
!!  mg = maximum length of 1D FFTs
!!  nd1=first dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd2=second dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  nd3=third dimension of (complex) arrays z and zbr (treated as real within
!!   this subroutine)
!!  n1=actual length of x and y transforms
!!  n2=upper i2 index, used for blocking, put usual n2 for usual ffty
!!  n2i=lower i2 index, used for blocking : the do-loop will be i2=n2i,n2
!!   put to 1 for usual ffty
!!  z(2,nd1,nd2,nd3)=INPUT array; destroyed by transformation
!!  trig, aft, now, bef, ind=provided by previous call to ctrig
!!   Note that in this routine (and in ctrig) the values in array trig are
!!   actually cos and tan, not cos and sin.  Use of tan allows advantageous
!!   use of FMA on the ibm rs6000.
!!  ic=number of (radix) factors of x transform length (from ctrig)
!!  ris=sign of exponential in transform (should be 1 or -1; real)
!!
!! OUTPUT
!!  zbr(2,nd1,nd2,nd3)=OUTPUT transformed array; no scaling applied
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      fourdp,sg_fft,sg_fftpad
!!
!! SOURCE

 subroutine sg_fftz(fftcache,mfac,mg,nd1,nd2,nd3,n1,n2i,n2,z,zbr,&
&           trig,aft,now,bef,ris,ind,ic)

 use basis_defs
 implicit none

!Arguments ------------------------------------
 integer :: ic,fftcache,mfac,mg,nd1,nd2,nd3,n1,n2,n2i,n3
!Dimensions of aft, now, bef, ind, and trig should agree with
!those in subroutine ctrig.
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: ris
 real(dp) :: trig(2,mg),z(2,nd1,nd2,nd3),zbr(2,nd1,nd2,nd3)

!Local variables-------------------------------
 integer :: b_i,i,ia,ib,indx,i2,j,jj,lot,ma,mb,ntb
 real(dp),parameter :: &
& cos2=0.3090169943749474d0,&   !cos(2.d0*pi/5.d0)
& cos4=-0.8090169943749474d0,&  !cos(4.d0*pi/5.d0)
& sin42=0.6180339887498948d0    !sin(4.d0*pi/5.d0)/sin(2.d0*pi/5.d0)
 real(dp) :: bb,cr2,cr2s,cr3,cr3p,cr4,cr5,ct2,ct3,ct4,ct5,&
& r,r1,r2,r25,r3,r34,r4,r5,s,sin2,s1,s2,s25,s3,s34,s4,s5
 character message*500

! *************************************************************************

!n12 occurs as a loop index repeated below; do z transform while 
!looping over all n12 lines of data
 
!Direct transformation (to ic-1), bitreversal will be in second part
!of routine

 do i=1,ic-1
  ntb=now(i)*bef(i)
  b_i=bef(i)
 
! Treat radix 4
  if (now(i)==4) then
   ia=0
 
!  First step of radix 4
   do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,ia,ib,n1,n2i,n2,ntb,ris,z)
    do i2=n2i,n2
     do j=1,n1
      r4=z(1,j,i2,ia*ntb+3*b_i+ib)
      s4=z(2,j,i2,ia*ntb+3*b_i+ib)
      r3=z(1,j,i2,ia*ntb+2*b_i+ib)
      s3=z(2,j,i2,ia*ntb+2*b_i+ib)
      r2=z(1,j,i2,ia*ntb+b_i+ib)
      s2=z(2,j,i2,ia*ntb+b_i+ib)
      r1=z(1,j,i2,ia*ntb+ib)
      s1=z(2,j,i2,ia*ntb+ib)

      r=r1 + r3
      s=r2 + r4
      z(1,j,i2,ia*ntb+ib) = r + s
      z(1,j,i2,ia*ntb+2*b_i+ib) = r - s
      r=r1 - r3
      s=s2 - s4
      z(1,j,i2,ia*ntb+b_i+ib) = r - s*ris
      z(1,j,i2,ia*ntb+3*b_i+ib) = r + s*ris
      r=s1 + s3
      s=s2 + s4
      z(2,j,i2,ia*ntb+ib) = r + s
      z(2,j,i2,ia*ntb+2*b_i+ib) = r - s
      r=s1 - s3
      s=r2 - r4
      z(2,j,i2,ia*ntb+b_i+ib) = r + s*ris
      z(2,j,i2,ia*ntb+3*b_i+ib) = r - s*ris
     enddo ! j
    enddo ! i2 
!$OMP END PARALLEL DO
   enddo ! ib
 
!  Second step of radix 4
   do ia=1,aft(i)-1
    indx=ind(ia*4*b_i+1)-1
    indx=indx*b_i
    cr2=trig(1,indx)
    ct2=trig(2,indx)
    cr3=trig(1,2*indx)
    ct3=trig(2,2*indx)
    cr4=trig(1,3*indx)
    ct4=trig(2,3*indx)
    cr4=cr4/cr2
    cr2s=cr2*ris
    do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,cr2,cr3,cr4,ct2,cr2s,ct3,ct4,i,ia,ib,n1,n2i,n2,ntb,ris,z)
     do i2=n2i,n2
      do j=1,n1
       r4=z(1,j,i2,ia*ntb+3*b_i+ib) - &
&         z(2,j,i2,ia*ntb+3*b_i+ib)*ct4
       s4=z(1,j,i2,ia*ntb+3*b_i+ib)*ct4 + &
&         z(2,j,i2,ia*ntb+3*b_i+ib)
       r3=z(1,j,i2,ia*ntb+2*b_i+ib) - &
&         z(2,j,i2,ia*ntb+2*b_i+ib)*ct3
       s3=z(1,j,i2,ia*ntb+2*b_i+ib)*ct3 + &
&         z(2,j,i2,ia*ntb+2*b_i+ib)
       r2=z(1,j,i2,ia*ntb+b_i+ib) - &
&         z(2,j,i2,ia*ntb+b_i+ib)*ct2
       s2=z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&         z(2,j,i2,ia*ntb+b_i+ib)
       r1=z(1,j,i2,ia*ntb+ib)
       s1=z(2,j,i2,ia*ntb+ib)

       r=r1 + r3*cr3
       s=r2 + r4*cr4 
       z(1,j,i2,ia*ntb+ib) = r + s*cr2
       z(1,j,i2,ia*ntb+2*b_i+ib) = r - s*cr2
       r=r1 - r3*cr3
       s=s2 - s4*cr4 
       z(1,j,i2,ia*ntb+b_i+ib) = r - s*cr2s
       z(1,j,i2,ia*ntb+3*b_i+ib) = r + s*cr2s
       r=s1 + s3*cr3
       s=s2 + s4*cr4 
       z(2,j,i2,ia*ntb+ib) = r + s*cr2
       z(2,j,i2,ia*ntb+2*b_i+ib) = r - s*cr2
       r=s1 - s3*cr3
       s=r2 - r4*cr4 
       z(2,j,i2,ia*ntb+b_i+ib) = r + s*cr2s
       z(2,j,i2,ia*ntb+3*b_i+ib) = r - s*cr2s
       enddo ! j
      enddo ! i2 
!$OMP END PARALLEL DO
    enddo ! ib

   enddo ! ia
 
! Treat radix 2
  else if (now(i)==2) then
   ia=0
 
!  First step of radix 2
   do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,ia,ib,n1,n2,n2i,ntb,z)
    do i2=n2i,n2
     do j=1,n1
      r1=z(1,j,i2,ia*ntb+ib)
      s1=z(2,j,i2,ia*ntb+ib)
      r2=z(1,j,i2,ia*ntb+b_i+ib)
      s2=z(2,j,i2,ia*ntb+b_i+ib)
      z(1,j,i2,ia*ntb+ib) =  r2 + r1
      z(2,j,i2,ia*ntb+ib) =  s2 + s1
      z(1,j,i2,ia*ntb+b_i+ib) = -r2 + r1
      z(2,j,i2,ia*ntb+b_i+ib) = -s2 + s1
     enddo ! j
    enddo ! i2 
!$OMP END PARALLEL DO
   enddo ! ib
 
!  Second step of radix 2
   do ia=1,aft(i)-1
    indx=ind(ia*2*b_i+1)-1
    indx=indx*b_i
    cr2=trig(1,indx)
    ct2=trig(2,indx)
    do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,cr2,ct2,ia,ib,n1,n2,n2i,ntb,z)
     do i2=n2i,n2
      do j=1,n1
       r1=z(1,j,i2,ia*ntb+ib)
       s1=z(2,j,i2,ia*ntb+ib)
       r2=z(1,j,i2,ia*ntb+b_i+ib) - &
&         z(2,j,i2,ia*ntb+b_i+ib)*ct2
       s2=z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&         z(2,j,i2,ia*ntb+b_i+ib)
       z(1,j,i2,ia*ntb+ib) =  r2*cr2 + r1
       z(2,j,i2,ia*ntb+ib) =  s2*cr2 + s1
       z(1,j,i2,ia*ntb+b_i+ib) = -r2*cr2 + r1
       z(2,j,i2,ia*ntb+b_i+ib) = -s2*cr2 + s1
      enddo ! j
     enddo ! i2 
!$OMP END PARALLEL DO
    enddo ! ib

   enddo ! ia
 
! Treat radix 3
  else if (now(i)==3) then  
!  .5d0*sqrt(3.d0)=0.8660254037844387d0
   ia=0
   bb=ris*0.8660254037844387d0
 
!  First step of radix 3
   do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,b_i,ia,ib,n1,n2,n2i,ntb,z)
    do i2=n2i,n2
     do j=1,n1
      r1=z(1,j,i2,ia*ntb+ib)
      s1=z(2,j,i2,ia*ntb+ib)
      r2=z(1,j,i2,ia*ntb+b_i+ib)
      s2=z(2,j,i2,ia*ntb+b_i+ib)
      r3=z(1,j,i2,ia*ntb+2*b_i+ib)
      s3=z(2,j,i2,ia*ntb+2*b_i+ib)
      r=r2 + r3
      s=s2 + s3
      z(1,j,i2,ia*ntb+ib) = r + r1
      z(2,j,i2,ia*ntb+ib) = s + s1
      r1=r1 - r*.5d0
      s1=s1 - s*.5d0
      r2=r2-r3
      s2=s2-s3
      z(1,j,i2,ia*ntb+b_i+ib) = r1 - s2*bb
      z(2,j,i2,ia*ntb+b_i+ib) = s1 + r2*bb
      z(1,j,i2,ia*ntb+2*b_i+ib) = r1 + s2*bb
      z(2,j,i2,ia*ntb+2*b_i+ib) = s1 - r2*bb
     enddo ! j
    enddo ! i2 
!$OMP END PARALLEL DO
   enddo ! ib 
 
!  Second step of radix 3
   do ia=1,aft(i)-1
    indx=ind(ia*3*b_i+1)-1
    indx=indx*b_i
    cr2=trig(1,indx)
    ct2=trig(2,indx)
    cr3=trig(1,2*indx)
    ct3=trig(2,2*indx)
    cr2=cr2/cr3
    cr3p=.5d0*cr3
    bb=ris*cr3*0.8660254037844387d0
    do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,b_i,cr2,cr3,cr3p,ct2,ct3,ia,ib,n1,n2,n2i,ntb,z)
     do i2=n2i,n2
      do j=1,n1
       r1=z(1,j,i2,ia*ntb+ib)
       s1=z(2,j,i2,ia*ntb+ib)
       r2=z(1,j,i2,ia*ntb+b_i+ib) - &
&         z(2,j,i2,ia*ntb+b_i+ib)*ct2
       s2=z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&         z(2,j,i2,ia*ntb+b_i+ib)
       r3=z(1,j,i2,ia*ntb+2*b_i+ib) - &
&         z(2,j,i2,ia*ntb+2*b_i+ib)*ct3
       s3=z(1,j,i2,ia*ntb+2*b_i+ib)*ct3 + &
&         z(2,j,i2,ia*ntb+2*b_i+ib)
       r=cr2*r2 + r3
       s=cr2*s2 + s3
       z(1,j,i2,ia*ntb+ib) = r*cr3 + r1
       z(2,j,i2,ia*ntb+ib) = s*cr3 + s1
       r1=r1 - r*cr3p
       s1=s1 - s*cr3p
       r2=cr2*r2-r3
       s2=cr2*s2-s3
       z(1,j,i2,ia*ntb+b_i+ib) = r1 - s2*bb
       z(2,j,i2,ia*ntb+b_i+ib) = s1 + r2*bb
       z(1,j,i2,ia*ntb+2*b_i+ib) = r1 + s2*bb
       z(2,j,i2,ia*ntb+2*b_i+ib) = s1 - r2*bb
      enddo ! j
     enddo ! i2 
!$OMP END PARALLEL DO
    enddo ! ib

   enddo ! ia
 
! Treat radix 5
  else if (now(i)==5) then
   sin2=ris*0.9510565162951536d0
   ia=0
 
!  First step of radix 5
   do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,ia,ib,n1,n2,n2i,ntb,sin2,z)
    do i2=n2i,n2
     do j=1,n1
      r1=z(1,j,i2,ia*ntb+ib)
      s1=z(2,j,i2,ia*ntb+ib)
      r2=z(1,j,i2,ia*ntb+b_i+ib)
      s2=z(2,j,i2,ia*ntb+b_i+ib)
      r3=z(1,j,i2,ia*ntb+2*b_i+ib)
      s3=z(2,j,i2,ia*ntb+2*b_i+ib)
      r4=z(1,j,i2,ia*ntb+3*b_i+ib)
      s4=z(2,j,i2,ia*ntb+3*b_i+ib)
      r5=z(1,j,i2,ia*ntb+4*b_i+ib)
      s5=z(2,j,i2,ia*ntb+4*b_i+ib)
      r25 = r2 + r5
      r34 = r3 + r4
      s25 = s2 - s5
      s34 = s3 - s4
      z(1,j,i2,ia*ntb+ib) = r1 + r25 + r34
      r = r1 + cos2*r25 + cos4*r34
      s = s25 + sin42*s34
      z(1,j,i2,ia*ntb+b_i+ib) = r - sin2*s
      z(1,j,i2,ia*ntb+4*b_i+ib) = r + sin2*s
      r = r1 + cos4*r25 + cos2*r34
      s = sin42*s25 - s34
      z(1,j,i2,ia*ntb+2*b_i+ib) = r - sin2*s
      z(1,j,i2,ia*ntb+3*b_i+ib) = r + sin2*s
      r25 = r2 - r5
      r34 = r3 - r4
      s25 = s2 + s5
      s34 = s3 + s4
      z(2,j,i2,ia*ntb+ib) = s1 + s25 + s34
      r = s1 + cos2*s25 + cos4*s34
      s = r25 + sin42*r34
      z(2,j,i2,ia*ntb+b_i+ib) = r + sin2*s
      z(2,j,i2,ia*ntb+4*b_i+ib) = r - sin2*s
      r = s1 + cos4*s25 + cos2*s34
      s = sin42*r25 - r34
      z(2,j,i2,ia*ntb+2*b_i+ib) = r + sin2*s
      z(2,j,i2,ia*ntb+3*b_i+ib) = r - sin2*s
     enddo ! j
    enddo ! i2 
!$OMP END PARALLEL DO
   enddo ! ib
 
!  Second step of radix 5
   do ia=1,aft(i)-1
    indx=ind(ia*5*b_i+1)-1
    indx=indx*b_i
    cr2=trig(1,indx)
    ct2=trig(2,indx)
    cr3=trig(1,2*indx)
    ct3=trig(2,2*indx)
    cr4=trig(1,3*indx)
    ct4=trig(2,3*indx)
    cr5=trig(1,4*indx)
    ct5=trig(2,4*indx)
    do ib=1,b_i 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(b_i,cr2,cr3,cr4,cr5,ct2,ct3,ct4,ct5,ia,ib,n1,n2,n2i,ntb,sin2,z)
     do i2=n2i,n2
      do j=1,n1
       r1=z(1,j,i2,ia*ntb+ib)
       s1=z(2,j,i2,ia*ntb+ib)
       r2=cr2*(z(1,j,i2,ia*ntb+b_i+ib) - &
&              z(2,j,i2,ia*ntb+b_i+ib)*ct2)
       s2=cr2*(z(1,j,i2,ia*ntb+b_i+ib)*ct2 + &
&              z(2,j,i2,ia*ntb+b_i+ib))
       r3=cr3*(z(1,j,i2,ia*ntb+2*b_i+ib) - &
&              z(2,j,i2,ia*ntb+2*b_i+ib)*ct3)
       s3=cr3*(z(1,j,i2,ia*ntb+2*b_i+ib)*ct3 + &
&              z(2,j,i2,ia*ntb+2*b_i+ib))
       r4=z(1,j,i2,ia*ntb+3*b_i+ib) - &
&         z(2,j,i2,ia*ntb+3*b_i+ib)*ct4
       s4=z(1,j,i2,ia*ntb+3*b_i+ib)*ct4 + &
&         z(2,j,i2,ia*ntb+3*b_i+ib)
       r5=z(1,j,i2,ia*ntb+4*b_i+ib) - &
&         z(2,j,i2,ia*ntb+4*b_i+ib)*ct5
       s5=z(1,j,i2,ia*ntb+4*b_i+ib)*ct5 + &
&         z(2,j,i2,ia*ntb+4*b_i+ib)
       r25 = r2 + r5*cr5
       r34 = r3 + r4*cr4
       s25 = s2 - s5*cr5
       s34 = s3 - s4*cr4
       z(1,j,i2,ia*ntb+ib) = r1 + r25 + r34
       r = r1 + cos2*r25 + cos4*r34
       s = s25 + sin42*s34
       z(1,j,i2,ia*ntb+b_i+ib) = r - sin2*s
       z(1,j,i2,ia*ntb+4*b_i+ib) = r + sin2*s
       r = r1 + cos4*r25 + cos2*r34
       s = sin42*s25 - s34
       z(1,j,i2,ia*ntb+2*b_i+ib) = r - sin2*s
       z(1,j,i2,ia*ntb+3*b_i+ib) = r + sin2*s
       r25 = r2 - r5*cr5
       r34 = r3 - r4*cr4
       s25 = s2 + s5*cr5
       s34 = s3 + s4*cr4
       z(2,j,i2,ia*ntb+ib) = s1 + s25 + s34
       r = s1 + cos2*s25 + cos4*s34
       s = r25 + sin42*r34
       z(2,j,i2,ia*ntb+b_i+ib) = r + sin2*s
       z(2,j,i2,ia*ntb+4*b_i+ib) = r - sin2*s
       r = s1 + cos4*s25 + cos2*s34
       s = sin42*r25 - r34
       z(2,j,i2,ia*ntb+2*b_i+ib) = r + sin2*s
       z(2,j,i2,ia*ntb+3*b_i+ib) = r - sin2*s
      enddo ! j
     enddo ! i2 
!$OMP END PARALLEL DO
    enddo ! ib

   enddo ! ia
 
! All radices treated
  else
   write(message, '(a,a,a,a)' )ch10,&
&   ' sg_fftz : BUG -',ch10,&
&   '  called with factors other than 2, 3, and 5'
!   call wrtout(06,message,'PERS')
!   call leave_new('PERS')
  endif
 
!End of direct transformation
 enddo
 
! ------------------------------------------------------------
! bitreversal  (zbr is for z"bit-reversed")
 
!Treat radix 4
 if (now(ic)==4) then
  ia=0
 
! First step of radix 4
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,ind,n1,n2,n2i,ntb,ris,z,zbr)
  do i2=n2i,n2
   do j=1,n1
    r4=z(1,j,i2,ia*4+4)
    s4=z(2,j,i2,ia*4+4)
    r3=z(1,j,i2,ia*4+3)
    s3=z(2,j,i2,ia*4+3)
    r2=z(1,j,i2,ia*4+2)
    s2=z(2,j,i2,ia*4+2)
    r1=z(1,j,i2,ia*4+1)
    s1=z(2,j,i2,ia*4+1)

    r=r1 + r3
    s=r2 + r4
    zbr(1,j,i2,ind(ia*4+1)) = r + s
    zbr(1,j,i2,ind(ia*4+3)) = r - s
    r=r1 - r3
    s=s2 - s4
    zbr(1,j,i2,ind(ia*4+2)) = r - s*ris
    zbr(1,j,i2,ind(ia*4+4)) = r + s*ris
    r=s1 + s3
    s=s2 + s4
    zbr(2,j,i2,ind(ia*4+1)) = r + s
    zbr(2,j,i2,ind(ia*4+3)) = r - s
    r=s1 - s3
    s=r2 - r4
    zbr(2,j,i2,ind(ia*4+2)) = r + s*ris
    zbr(2,j,i2,ind(ia*4+4)) = r - s*ris
   enddo ! j
  enddo ! i2 
!$OMP END PARALLEL DO
 
! Second step of radix 4
  do ia=1,aft(ic)-1
   indx=ind(ia*4+1)-1
   cr2=trig(1,indx)
   ct2=trig(2,indx)
   cr3=trig(1,2*indx)
   ct3=trig(2,2*indx)
   cr4=trig(1,3*indx)
   ct4=trig(2,3*indx)
   cr4=cr4/cr2
   cr2s=cr2*ris 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,cr2,cr2s,cr3,cr4,ct2,ct3,ct4,ind,n1,n2,n2i,z,zbr)
   do i2=n2i,n2
    do j=1,n1
     r4=z(1,j,i2,ia*4+4) - z(2,j,i2,ia*4+4)*ct4
     s4=z(1,j,i2,ia*4+4)*ct4 + z(2,j,i2,ia*4+4)
     r3=z(1,j,i2,ia*4+3) - z(2,j,i2,ia*4+3)*ct3
     s3=z(1,j,i2,ia*4+3)*ct3 + z(2,j,i2,ia*4+3)
     r2=z(1,j,i2,ia*4+2) - z(2,j,i2,ia*4+2)*ct2
     s2=z(1,j,i2,ia*4+2)*ct2 + z(2,j,i2,ia*4+2)
     r1=z(1,j,i2,ia*4+1)
     s1=z(2,j,i2,ia*4+1)

     r=r1 + r3*cr3
     s=r2 + r4*cr4 
     zbr(1,j,i2,ind(ia*4+1)) = r + s*cr2
     zbr(1,j,i2,ind(ia*4+3)) = r - s*cr2
     r=r1 - r3*cr3
     s=s2 - s4*cr4 
     zbr(1,j,i2,ind(ia*4+2)) = r - s*cr2s
     zbr(1,j,i2,ind(ia*4+4)) = r + s*cr2s
     r=s1 + s3*cr3
     s=s2 + s4*cr4 
     zbr(2,j,i2,ind(ia*4+1)) = r + s*cr2
     zbr(2,j,i2,ind(ia*4+3)) = r - s*cr2
     r=s1 - s3*cr3
     s=r2 - r4*cr4 
     zbr(2,j,i2,ind(ia*4+2)) = r + s*cr2s
     zbr(2,j,i2,ind(ia*4+4)) = r - s*cr2s
    enddo ! j
   enddo ! i2 
!$OMP END PARALLEL DO

  enddo ! ia
 
!Treat radix 2
 else if (now(ic)==2) then
  ia=0
 
! First step of radix 2
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,ind,n1,n2,n2i,z,zbr)
  do i2=n2i,n2
   do j=1,n1
    r1=z(1,j,i2,ia*2+1)
    s1=z(2,j,i2,ia*2+1)
    r2=z(1,j,i2,ia*2+2)
    s2=z(2,j,i2,ia*2+2)
    zbr(1,j,i2,ind(ia*2+1)) =  r2 + r1
    zbr(2,j,i2,ind(ia*2+1)) =  s2 + s1
    zbr(1,j,i2,ind(ia*2+2)) = -r2 + r1
    zbr(2,j,i2,ind(ia*2+2)) = -s2 + s1
   enddo ! j
  enddo ! i2 
!$OMP END PARALLEL DO
 
! Second step of radix 2
  do ia=1,aft(ic)-1
   indx=ind(ia*2+1)-1
   cr2=trig(1,indx)
   ct2=trig(2,indx) 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(cr2,ct2,ia,ind,n1,n2,n2i,z,zbr)
   do i2=n2i,n2
    do j=1,n1
     r1=z(1,j,i2,ia*2+1)
     s1=z(2,j,i2,ia*2+1)
     r2=z(1,j,i2,ia*2+2) - z(2,j,i2,ia*2+2)*ct2
     s2=z(1,j,i2,ia*2+2)*ct2 + z(2,j,i2,ia*2+2)
     zbr(1,j,i2,ind(ia*2+1)) =  r2*cr2 + r1
     zbr(2,j,i2,ind(ia*2+1)) =  s2*cr2 + s1
     zbr(1,j,i2,ind(ia*2+2)) = -r2*cr2 + r1
     zbr(2,j,i2,ind(ia*2+2)) = -s2*cr2 + s1
    enddo ! j
   enddo ! i2 
!$OMP END PARALLEL DO
  enddo ! ia
 
!Treat radix 3
 else if (now(ic)==3) then  
! .5d0*sqrt(3.d0)=0.8660254037844387d0
  ia=0
  bb=ris*0.8660254037844387d0
 
! First step of radix 3
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,ia,ind,n1,n2,n2i,z,zbr)
  do i2=n2i,n2
   do j=1,n1
    r1=z(1,j,i2,ia*3+1)
    s1=z(2,j,i2,ia*3+1)
    r2=z(1,j,i2,ia*3+2)
    s2=z(2,j,i2,ia*3+2)
    r3=z(1,j,i2,ia*3+3)
    s3=z(2,j,i2,ia*3+3)
    r=r2 + r3
    s=s2 + s3
    zbr(1,j,i2,ind(ia*3+1)) = r + r1
    zbr(2,j,i2,ind(ia*3+1)) = s + s1
    r1=r1 - r*.5d0
    s1=s1 - s*.5d0
    r2=r2-r3
    s2=s2-s3
    zbr(1,j,i2,ind(ia*3+2)) = r1 - s2*bb
    zbr(2,j,i2,ind(ia*3+2)) = s1 + r2*bb
    zbr(1,j,i2,ind(ia*3+3)) = r1 + s2*bb
    zbr(2,j,i2,ind(ia*3+3)) = s1 - r2*bb
   enddo ! j
  enddo ! i2 
!$OMP END PARALLEL DO
 
! Second step of radix 3
  do ia=1,aft(ic)-1
   indx=ind(ia*3+1)-1
   cr2=trig(1,indx)
   ct2=trig(2,indx)
   cr3=trig(1,2*indx)
   ct3=trig(2,2*indx)
   cr2=cr2/cr3
   cr3p=.5d0*cr3
   bb=ris*cr3*0.8660254037844387d0 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(bb,cr2,cr3,cr3p,ct2,ct3,ia,ind,n1,n2,n2i,z,zbr)
   do i2=n2i,n2
    do j=1,n1
     r1=z(1,j,i2,ia*3+1)
     s1=z(2,j,i2,ia*3+1)
     r2=z(1,j,i2,ia*3+2) - z(2,j,i2,ia*3+2)*ct2
     s2=z(1,j,i2,ia*3+2)*ct2 + z(2,j,i2,ia*3+2)
     r3=z(1,j,i2,ia*3+3) - z(2,j,i2,ia*3+3)*ct3
     s3=z(1,j,i2,ia*3+3)*ct3 + z(2,j,i2,ia*3+3)
     r=cr2*r2 + r3
     s=cr2*s2 + s3
     zbr(1,j,i2,ind(ia*3+1)) = r*cr3 + r1
     zbr(2,j,i2,ind(ia*3+1)) = s*cr3 + s1
     r1=r1 - r*cr3p
     s1=s1 - s*cr3p
     r2=cr2*r2-r3
     s2=cr2*s2-s3
     zbr(1,j,i2,ind(ia*3+2)) = r1 - s2*bb
     zbr(2,j,i2,ind(ia*3+2)) = s1 + r2*bb
     zbr(1,j,i2,ind(ia*3+3)) = r1 + s2*bb
     zbr(2,j,i2,ind(ia*3+3)) = s1 - r2*bb
    enddo ! j
   enddo ! i2 
!$OMP END PARALLEL DO
  enddo ! ia
 
!Treat radix 5
 else if (now(ic)==5) then 
! sin(2.d0*pi/5.d0)
  sin2=ris*0.9510565162951536d0
  ia=0
 
! First step of radix 5
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(ia,ind,n1,n2,n2i,sin2,z,zbr)
  do i2=n2i,n2
   do j=1,n1
    r1=z(1,j,i2,ia*5+1)
    s1=z(2,j,i2,ia*5+1)
    r2=z(1,j,i2,ia*5+2)
    s2=z(2,j,i2,ia*5+2)
    r3=z(1,j,i2,ia*5+3)
    s3=z(2,j,i2,ia*5+3)
    r4=z(1,j,i2,ia*5+4)
    s4=z(2,j,i2,ia*5+4)
    r5=z(1,j,i2,ia*5+5)
    s5=z(2,j,i2,ia*5+5)
    r25 = r2 + r5
    r34 = r3 + r4
    s25 = s2 - s5
    s34 = s3 - s4
    zbr(1,j,i2,ind(ia*5+1)) = r1 + r25 + r34
    r = r1 + cos2*r25 + cos4*r34
    s = s25 + sin42*s34
    zbr(1,j,i2,ind(ia*5+2)) = r - sin2*s
    zbr(1,j,i2,ind(ia*5+5)) = r + sin2*s
    r = r1 + cos4*r25 + cos2*r34
    s = sin42*s25 - s34
    zbr(1,j,i2,ind(ia*5+3)) = r - sin2*s
    zbr(1,j,i2,ind(ia*5+4)) = r + sin2*s
    r25 = r2 - r5
    r34 = r3 - r4
    s25 = s2 + s5
    s34 = s3 + s4
    zbr(2,j,i2,ind(ia*5+1)) = s1 + s25 + s34
    r = s1 + cos2*s25 + cos4*s34
    s = r25 + sin42*r34
    zbr(2,j,i2,ind(ia*5+2)) = r + sin2*s
    zbr(2,j,i2,ind(ia*5+5)) = r - sin2*s
    r = s1 + cos4*s25 + cos2*s34
    s = sin42*r25 - r34
    zbr(2,j,i2,ind(ia*5+3)) = r + sin2*s
    zbr(2,j,i2,ind(ia*5+4)) = r - sin2*s
   enddo ! j
  enddo ! i2 
!$OMP END PARALLEL DO
 
! Second step of radix 5
  do ia=1,aft(ic)-1
   indx=ind(ia*5+1)-1
   cr2=trig(1,indx)
   ct2=trig(2,indx)
   cr3=trig(1,2*indx)
   ct3=trig(2,2*indx)
   cr4=trig(1,3*indx)
   ct4=trig(2,3*indx)
   cr5=trig(1,4*indx)
   ct5=trig(2,4*indx) 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(cr2,cr3,cr4,cr5,ct2,ct3,ct4,ct5,ia,ind,n1,n2,n2i,sin2,z,zbr) 
   do i2=n2i,n2
    do j=1,n1
     r1=z(1,j,i2,ia*5+1)
     s1=z(2,j,i2,ia*5+1)
     r2=cr2*(z(1,j,i2,ia*5+2) - z(2,j,i2,ia*5+2)*ct2)
     s2=cr2*(z(1,j,i2,ia*5+2)*ct2 + z(2,j,i2,ia*5+2))
     r3=cr3*(z(1,j,i2,ia*5+3) - z(2,j,i2,ia*5+3)*ct3)
     s3=cr3*(z(1,j,i2,ia*5+3)*ct3 + z(2,j,i2,ia*5+3))
     r4=z(1,j,i2,ia*5+4) - z(2,j,i2,ia*5+4)*ct4
     s4=z(1,j,i2,ia*5+4)*ct4 + z(2,j,i2,ia*5+4)
     r5=z(1,j,i2,ia*5+5) - z(2,j,i2,ia*5+5)*ct5
     s5=z(1,j,i2,ia*5+5)*ct5 + z(2,j,i2,ia*5+5)
     r25 = r2 + r5*cr5
     r34 = r3 + r4*cr4
     s25 = s2 - s5*cr5
     s34 = s3 - s4*cr4
     zbr(1,j,i2,ind(ia*5+1)) = r1 + r25 + r34
     r = r1 + cos2*r25 + cos4*r34
     s = s25 + sin42*s34
     zbr(1,j,i2,ind(ia*5+2)) = r - sin2*s
     zbr(1,j,i2,ind(ia*5+5)) = r + sin2*s
     r = r1 + cos4*r25 + cos2*r34
     s = sin42*s25 - s34
     zbr(1,j,i2,ind(ia*5+3)) = r - sin2*s
     zbr(1,j,i2,ind(ia*5+4)) = r + sin2*s
     r25 = r2 - r5*cr5
     r34 = r3 - r4*cr4
     s25 = s2 + s5*cr5
     s34 = s3 + s4*cr4
     zbr(2,j,i2,ind(ia*5+1)) = s1 + s25 + s34
     r = s1 + cos2*s25 + cos4*s34
     s = r25 + sin42*r34
     zbr(2,j,i2,ind(ia*5+2)) = r + sin2*s
     zbr(2,j,i2,ind(ia*5+5)) = r - sin2*s
     r = s1 + cos4*s25 + cos2*s34
     s = sin42*r25 - r34
     zbr(2,j,i2,ind(ia*5+3)) = r + sin2*s
     zbr(2,j,i2,ind(ia*5+4)) = r - sin2*s
    enddo ! j
   enddo ! i2 
!$OMP END PARALLEL DO
  enddo ! ia
 
!All radices treated
 else
  write(message, '(a,a,a,a)' )ch10,&
&  ' sg_fftz : BUG -',ch10,&
&  '  called with factors other than 2, 3, and 5'
!  call wrtout(06,message,'PERS')
!  call leave_new('PERS')
 endif

 end subroutine
!!***
