!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fftx
!! NAME
!! sg_fftx
!!
!! FUNCTION
!! This subroutine is called by the 3-dimensional fft to conduct the
!! "x" transforms for all y and z.
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
!!  n2,n3=actual length of y and z transforms
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
!! This routine blocks the x transforms
!! so that all transforms under consideration at one step fit within
!! the cache memory, which is crucial for optimal performance.
!! The blocking factor is set by parameter "fftcache" below, which should
!! be adjusted to be somewhat smaller (say 3/4) than the actual cache size
!! of the machine.
!!
!! TODO
!! Use latex for the equation above
!!
!! PARENTS
!!      fourdp,sg_fft,sg_fftrisc
!!
!! SOURCE

 subroutine sg_fftx(fftcache,mfac,mg,nd1,nd2,nd3,n2,n3,z,zbr,&
& trig,aft,now,bef,ris,ind,ic)
 
 use basis_defs
 implicit none

!Arguments ------------------------------------
 integer :: fftcache,ic,mfac,mg,nd1,nd2,nd3,n1,n2,n3
!Dimensions of aft, now, bef, ind, and trig should agree with
!those in subroutine ctrig.
 integer :: aft(mfac),bef(mfac),ind(mg),now(mfac)
 real(dp) :: ris
 real(dp) :: trig(2,mg),z(2,nd1,nd2,nd3),zbr(2,nd1,nd2,nd3)

!Local variables-------------------------------
 integer :: i,ia,ib,indx,i3,j,jj,lot,ma,mb,ntb
 real(dp),parameter :: &
& cos2=0.3090169943749474d0,&   !cos(2.d0*pi/5.d0)
& cos4=-0.8090169943749474d0,&  !cos(4.d0*pi/5.d0)
& sin42=0.6180339887498948d0    !sin(4.d0*pi/5.d0)/sin(2.d0*pi/5.d0)
 real(dp) :: bb,cr2,cr2s,cr3,cr3p,cr4,cr5,ct2,ct3,ct4,ct5,&
& factor,r,r1,r2,r25,r3,r34,r4,r5,s,sin2,s1,s2,s25,s3,s34,s4,s5
 character message*500 

! *************************************************************************

!Do x transforms in blocks of size "lot" which is set by how
!many x transform arrays (of size nd1 each) fit into the nominal
!cache size "fftcache".
 factor=0.75d0
 lot=(fftcache*factor*1000d0)/(nd1*8*2)

!XG : due to the dimension problems on the P6, I have slightly
!modified this part of the code, with an external loop
!on n3 ...
!Modifications are indicated explicitely, or
!are related to the increase of the number of dimensions of z and
!zbr ...

 factor=0.75d0
 lot=(fftcache*factor*1000d0)/(nd1*8*2)
 
!$OMP PARALLEL DO DEFAULT(PRIVATE)&
!$OMP&SHARED(aft,bef,ic,ind,lot,n2,n3,now,ris,trig,z,zbr)
 do i3=1,n3
  do jj=1,n2,lot 
!  end of modification
 
!  For each jj, ma and mb give starting and ending addresses for fft
!  ma starts where we left off after last block
   ma=jj 
!  mb runs to the end of the block or else to the end of the data
!  modified XG 980107
!  mb=min(jj+(lot-1),n23)
   mb=min(jj+(lot-1),n2) 
 
!  Run over all factors except the last (to ic-1), performing
!  x transform
 
!  Note: fortran should skip this loop if ic=1; beware "onetrip"
!  compiler option which forces each loop at least once
 
!  ------------------------------------------------------------------------
 
!  Direct transformation (to be followed by bit reversal)

   do i=1,ic-1
    ntb=now(i)*bef(i) 
!   radix 4
 
!   Treat radix 4
    if (now(i)==4) then
     ia=0
 
!    First step of factor 4
     do ib=1,bef(i)
      do j=ma,mb
       r4=z(1,ia*ntb+3*bef(i)+ib,j,i3)
       s4=z(2,ia*ntb+3*bef(i)+ib,j,i3)
       r3=z(1,ia*ntb+2*bef(i)+ib,j,i3)
       s3=z(2,ia*ntb+2*bef(i)+ib,j,i3)
       r2=z(1,ia*ntb+bef(i)+ib,j,i3)
       s2=z(2,ia*ntb+bef(i)+ib,j,i3)
       r1=z(1,ia*ntb+ib,j,i3)
       s1=z(2,ia*ntb+ib,j,i3)

       r=r1 + r3
       s=r2 + r4
       z(1,ia*ntb+ib,j,i3) = r + s
       z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - s
       r=r1 - r3
       s=s2 - s4
       z(1,ia*ntb+bef(i)+ib,j,i3) = r - s*ris
       z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + s*ris
       r=s1 + s3
       s=s2 + s4
       z(2,ia*ntb+ib,j,i3) = r + s
       z(2,ia*ntb+2*bef(i)+ib,j,i3) = r - s
       r=s1 - s3
       s=r2 - r4
       z(2,ia*ntb+bef(i)+ib,j,i3) = r + s*ris
       z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - s*ris
      enddo
     enddo
 
!    Second step of factor 4
     do ia=1,aft(i)-1
      indx=ind(ia*4*bef(i)+1)-1
      indx=indx*bef(i)
      cr2=trig(1,indx)
      ct2=trig(2,indx)
      cr3=trig(1,2*indx)
      ct3=trig(2,2*indx)
      cr4=trig(1,3*indx)
      ct4=trig(2,3*indx)
      cr4=cr4/cr2
      cr2s=cr2*ris
      do ib=1,bef(i)
       do j=ma,mb
        r4=z(1,ia*ntb+3*bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+3*bef(i)+ib,j,i3)*ct4
        s4=z(1,ia*ntb+3*bef(i)+ib,j,i3)*ct4 + &
&          z(2,ia*ntb+3*bef(i)+ib,j,i3)
        r3=z(1,ia*ntb+2*bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+2*bef(i)+ib,j,i3)*ct3
        s3=z(1,ia*ntb+2*bef(i)+ib,j,i3)*ct3 + &
&          z(2,ia*ntb+2*bef(i)+ib,j,i3)
        r2=z(1,ia*ntb+bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+bef(i)+ib,j,i3)*ct2
        s2=z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&          z(2,ia*ntb+bef(i)+ib,j,i3)
        r1=z(1,ia*ntb+ib,j,i3)
        s1=z(2,ia*ntb+ib,j,i3)

        r=r1 + r3*cr3
        s=r2 + r4*cr4 
        z(1,ia*ntb+ib,j,i3) = r + s*cr2
        z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - s*cr2
        r=r1 - r3*cr3
        s=s2 - s4*cr4 
        z(1,ia*ntb+bef(i)+ib,j,i3) = r - s*cr2s
        z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + s*cr2s
        r=s1 + s3*cr3
        s=s2 + s4*cr4 
        z(2,ia*ntb+ib,j,i3) = r + s*cr2
        z(2,ia*ntb+2*bef(i)+ib,j,i3) = r - s*cr2
        r=s1 - s3*cr3
        s=r2 - r4*cr4 
        z(2,ia*ntb+bef(i)+ib,j,i3) = r + s*cr2s
        z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - s*cr2s
       enddo
      enddo
     enddo
 
!   Treat radix 2
    else if (now(i)==2) then 
     ia=0
 
!    First step of factor 2
     do ib=1,bef(i)
      do j=ma,mb
       r1=z(1,ia*ntb+ib,j,i3)
       s1=z(2,ia*ntb+ib,j,i3)
       r2=z(1,ia*ntb+bef(i)+ib,j,i3)
       s2=z(2,ia*ntb+bef(i)+ib,j,i3)
       z(1,ia*ntb+ib,j,i3) =  r2 + r1
       z(2,ia*ntb+ib,j,i3) =  s2 + s1
       z(1,ia*ntb+bef(i)+ib,j,i3) = -r2 + r1
       z(2,ia*ntb+bef(i)+ib,j,i3) = -s2 + s1
      enddo
     enddo
 
!    Second step of factor 2
     do ia=1,aft(i)-1
      indx=ind(ia*2*bef(i)+1)-1
      indx=indx*bef(i)
      cr2=trig(1,indx)
      ct2=trig(2,indx)
      do ib=1,bef(i)
       do j=ma,mb
        r1=z(1,ia*ntb+ib,j,i3)
        s1=z(2,ia*ntb+ib,j,i3)
        r2=z(1,ia*ntb+bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+bef(i)+ib,j,i3)*ct2
        s2=z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&          z(2,ia*ntb+bef(i)+ib,j,i3)
        z(1,ia*ntb+ib,j,i3) =  r2*cr2 + r1
        z(2,ia*ntb+ib,j,i3) =  s2*cr2 + s1
        z(1,ia*ntb+bef(i)+ib,j,i3) = -r2*cr2 + r1
        z(2,ia*ntb+bef(i)+ib,j,i3) = -s2*cr2 + s1
       enddo
      enddo
     enddo
 
!   Treat radix 3
    else if (now(i)==3) then  
!    .5d0*sqrt(3.d0)=0.8660254037844387d0
     ia=0
     bb=ris*0.8660254037844387d0
 
!    First step of factor 3
     do ib=1,bef(i)
      do j=ma,mb
       r1=z(1,ia*ntb+ib,j,i3)
       s1=z(2,ia*ntb+ib,j,i3)
       r2=z(1,ia*ntb+bef(i)+ib,j,i3)
       s2=z(2,ia*ntb+bef(i)+ib,j,i3)
       r3=z(1,ia*ntb+2*bef(i)+ib,j,i3)
       s3=z(2,ia*ntb+2*bef(i)+ib,j,i3)
       r=r2 + r3
       s=s2 + s3
       z(1,ia*ntb+ib,j,i3) = r + r1
       z(2,ia*ntb+ib,j,i3) = s + s1
       r1=r1 - r*.5d0
       s1=s1 - s*.5d0
       r2=r2-r3
       s2=s2-s3
       z(1,ia*ntb+bef(i)+ib,j,i3) = r1 - s2*bb
       z(2,ia*ntb+bef(i)+ib,j,i3) = s1 + r2*bb
       z(1,ia*ntb+2*bef(i)+ib,j,i3) = r1 + s2*bb
       z(2,ia*ntb+2*bef(i)+ib,j,i3) = s1 - r2*bb
      enddo
     enddo
 
!    Second step of factor 3
     do ia=1,aft(i)-1
      indx=ind(ia*3*bef(i)+1)-1
      indx=indx*bef(i)
      cr2=trig(1,indx)
      ct2=trig(2,indx)
      cr3=trig(1,2*indx)
      ct3=trig(2,2*indx)
      cr2=cr2/cr3
      cr3p=.5d0*cr3
      bb=ris*cr3*0.8660254037844387d0
      do ib=1,bef(i)
       do j=ma,mb
        r1=z(1,ia*ntb+ib,j,i3)
        s1=z(2,ia*ntb+ib,j,i3)
        r2=z(1,ia*ntb+bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+bef(i)+ib,j,i3)*ct2
        s2=z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&          z(2,ia*ntb+bef(i)+ib,j,i3)
        r3=z(1,ia*ntb+2*bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+2*bef(i)+ib,j,i3)*ct3
        s3=z(1,ia*ntb+2*bef(i)+ib,j,i3)*ct3 + &
&          z(2,ia*ntb+2*bef(i)+ib,j,i3)
        r=cr2*r2 + r3
        s=cr2*s2 + s3
        z(1,ia*ntb+ib,j,i3) = r*cr3 + r1
        z(2,ia*ntb+ib,j,i3) = s*cr3 + s1
        r1=r1 - r*cr3p
        s1=s1 - s*cr3p
        r2=cr2*r2-r3
        s2=cr2*s2-s3
        z(1,ia*ntb+bef(i)+ib,j,i3) = r1 - s2*bb
        z(2,ia*ntb+bef(i)+ib,j,i3) = s1 + r2*bb
        z(1,ia*ntb+2*bef(i)+ib,j,i3) = r1 + s2*bb
        z(2,ia*ntb+2*bef(i)+ib,j,i3) = s1 - r2*bb
       enddo
      enddo
     enddo
 
!   Treat radix 5
    else if (now(i)==5) then 
!    sin(2.d0*pi/5.d0)
     sin2=ris*0.9510565162951536d0
     ia=0
 
!    First step of factor 5
     do ib=1,bef(i)
      do j=ma,mb
       r1=z(1,ia*ntb+ib,j,i3)
       s1=z(2,ia*ntb+ib,j,i3)
       r2=z(1,ia*ntb+bef(i)+ib,j,i3)
       s2=z(2,ia*ntb+bef(i)+ib,j,i3)
       r3=z(1,ia*ntb+2*bef(i)+ib,j,i3)
       s3=z(2,ia*ntb+2*bef(i)+ib,j,i3)
       r4=z(1,ia*ntb+3*bef(i)+ib,j,i3)
       s4=z(2,ia*ntb+3*bef(i)+ib,j,i3)
       r5=z(1,ia*ntb+4*bef(i)+ib,j,i3)
       s5=z(2,ia*ntb+4*bef(i)+ib,j,i3)
       r25 = r2 + r5
       r34 = r3 + r4
       s25 = s2 - s5
       s34 = s3 - s4
       z(1,ia*ntb+ib,j,i3) = r1 + r25 + r34
       r = r1 + cos2*r25 + cos4*r34
       s = s25 + sin42*s34
       z(1,ia*ntb+bef(i)+ib,j,i3) = r - sin2*s
       z(1,ia*ntb+4*bef(i)+ib,j,i3) = r + sin2*s
       r = r1 + cos4*r25 + cos2*r34
       s = sin42*s25 - s34
       z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - sin2*s
       z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + sin2*s
       r25 = r2 - r5
       r34 = r3 - r4
       s25 = s2 + s5
       s34 = s3 + s4
       z(2,ia*ntb+ib,j,i3) = s1 + s25 + s34
       r = s1 + cos2*s25 + cos4*s34
       s = r25 + sin42*r34
       z(2,ia*ntb+bef(i)+ib,j,i3) = r + sin2*s
       z(2,ia*ntb+4*bef(i)+ib,j,i3) = r - sin2*s
       r = s1 + cos4*s25 + cos2*s34
       s = sin42*r25 - r34
       z(2,ia*ntb+2*bef(i)+ib,j,i3) = r + sin2*s
       z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - sin2*s
      enddo
     enddo
 
!    Second step of factor 5
     do ia=1,aft(i)-1
      indx=ind(ia*5*bef(i)+1)-1
      indx=indx*bef(i)
      cr2=trig(1,indx)
      ct2=trig(2,indx)
      cr3=trig(1,2*indx)
      ct3=trig(2,2*indx)
      cr4=trig(1,3*indx)
      ct4=trig(2,3*indx)
      cr5=trig(1,4*indx)
      ct5=trig(2,4*indx)
      do ib=1,bef(i)
       do j=ma,mb
        r1=z(1,ia*ntb+ib,j,i3)
        s1=z(2,ia*ntb+ib,j,i3)
        r2=cr2*(z(1,ia*ntb+bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+bef(i)+ib,j,i3)*ct2)
        s2=cr2*(z(1,ia*ntb+bef(i)+ib,j,i3)*ct2 + &
&               z(2,ia*ntb+bef(i)+ib,j,i3))
        r3=cr3*(z(1,ia*ntb+2*bef(i)+ib,j,i3) - &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3)*ct3)
        s3=cr3*(z(1,ia*ntb+2*bef(i)+ib,j,i3)*ct3 + &
&               z(2,ia*ntb+2*bef(i)+ib,j,i3))
        r4=z(1,ia*ntb+3*bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+3*bef(i)+ib,j,i3)*ct4
        s4=z(1,ia*ntb+3*bef(i)+ib,j,i3)*ct4 + &
&          z(2,ia*ntb+3*bef(i)+ib,j,i3)
        r5=z(1,ia*ntb+4*bef(i)+ib,j,i3) - &
&          z(2,ia*ntb+4*bef(i)+ib,j,i3)*ct5
        s5=z(1,ia*ntb+4*bef(i)+ib,j,i3)*ct5 + &
&          z(2,ia*ntb+4*bef(i)+ib,j,i3)
        r25 = r2 + r5*cr5
        r34 = r3 + r4*cr4
        s25 = s2 - s5*cr5
        s34 = s3 - s4*cr4
        z(1,ia*ntb+ib,j,i3) = r1 + r25 + r34
        r = r1 + cos2*r25 + cos4*r34
        s = s25 + sin42*s34
        z(1,ia*ntb+bef(i)+ib,j,i3) = r - sin2*s
        z(1,ia*ntb+4*bef(i)+ib,j,i3) = r + sin2*s
        r = r1 + cos4*r25 + cos2*r34
        s = sin42*s25 - s34
        z(1,ia*ntb+2*bef(i)+ib,j,i3) = r - sin2*s
        z(1,ia*ntb+3*bef(i)+ib,j,i3) = r + sin2*s
        r25 = r2 - r5*cr5
        r34 = r3 - r4*cr4
        s25 = s2 + s5*cr5
        s34 = s3 + s4*cr4
        z(2,ia*ntb+ib,j,i3) = s1 + s25 + s34
        r = s1 + cos2*s25 + cos4*s34
        s = r25 + sin42*r34
        z(2,ia*ntb+bef(i)+ib,j,i3) = r + sin2*s
        z(2,ia*ntb+4*bef(i)+ib,j,i3) = r - sin2*s
        r = s1 + cos4*s25 + cos2*s34
        s = sin42*r25 - r34
        z(2,ia*ntb+2*bef(i)+ib,j,i3) = r + sin2*s
        z(2,ia*ntb+3*bef(i)+ib,j,i3) = r - sin2*s
       enddo
      enddo
     enddo

    else
 
!    All factors have been treated
     write(message, '(a,a,a,a)' ) ch10,&
&     ' sg_fftx : BUG -',ch10,&
&     '  called with factors other than 2, 3, and 5'
!     call wrtout(06,message,'PERS')
!     call leave_new('PERS')
    endif

   enddo
 
!  ---------------------------------------------------------------
 
!  bitreversal
 
!  Perform bit reversal on last factor of transformation
 
!  Treat factor 4
   if (now(ic)==4) then 
!   radix 4
    ia=0
  
!   First step of factor 4
    do j=ma,mb
     r4=z(1,ia*4+4,j,i3)
     s4=z(2,ia*4+4,j,i3)
     r3=z(1,ia*4+3,j,i3)
     s3=z(2,ia*4+3,j,i3)
     r2=z(1,ia*4+2,j,i3)
     s2=z(2,ia*4+2,j,i3)
     r1=z(1,ia*4+1,j,i3)
     s1=z(2,ia*4+1,j,i3)

     r=r1 + r3
     s=r2 + r4
     zbr(1,ind(ia*4+1),j,i3) = r + s
     zbr(1,ind(ia*4+3),j,i3) = r - s
     r=r1 - r3
     s=s2 - s4
     zbr(1,ind(ia*4+2),j,i3) = r - s*ris
     zbr(1,ind(ia*4+4),j,i3) = r + s*ris
     r=s1 + s3
     s=s2 + s4
     zbr(2,ind(ia*4+1),j,i3) = r + s
     zbr(2,ind(ia*4+3),j,i3) = r - s
     r=s1 - s3
     s=r2 - r4
     zbr(2,ind(ia*4+2),j,i3) = r + s*ris
     zbr(2,ind(ia*4+4),j,i3) = r - s*ris
    enddo
 
!   Second step of factor 4
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
     do j=ma,mb
      r4=z(1,ia*4+4,j,i3) - z(2,ia*4+4,j,i3)*ct4
      s4=z(1,ia*4+4,j,i3)*ct4 + z(2,ia*4+4,j,i3)
      r3=z(1,ia*4+3,j,i3) - z(2,ia*4+3,j,i3)*ct3
      s3=z(1,ia*4+3,j,i3)*ct3 + z(2,ia*4+3,j,i3)
      r2=z(1,ia*4+2,j,i3) - z(2,ia*4+2,j,i3)*ct2
      s2=z(1,ia*4+2,j,i3)*ct2 + z(2,ia*4+2,j,i3)
      r1=z(1,ia*4+1,j,i3)
      s1=z(2,ia*4+1,j,i3)

      r=r1 + r3*cr3
      s=r2 + r4*cr4 
      zbr(1,ind(ia*4+1),j,i3) = r + s*cr2
      zbr(1,ind(ia*4+3),j,i3) = r - s*cr2
      r=r1 - r3*cr3
      s=s2 - s4*cr4 
      zbr(1,ind(ia*4+2),j,i3) = r - s*cr2s
      zbr(1,ind(ia*4+4),j,i3) = r + s*cr2s
      r=s1 + s3*cr3
      s=s2 + s4*cr4 
      zbr(2,ind(ia*4+1),j,i3) = r + s*cr2
      zbr(2,ind(ia*4+3),j,i3) = r - s*cr2
      r=s1 - s3*cr3
      s=r2 - r4*cr4 
      zbr(2,ind(ia*4+2),j,i3) = r + s*cr2s
      zbr(2,ind(ia*4+4),j,i3) = r - s*cr2s
     enddo
    enddo
  
!  Treat factor 2
   else if (now(ic)==2) then 
!   radix 2
    ia=0
 
!   First step of factor 2
    do j=ma,mb
     r1=z(1,ia*2+1,j,i3)
     s1=z(2,ia*2+1,j,i3)
     r2=z(1,ia*2+2,j,i3)
     s2=z(2,ia*2+2,j,i3)
     zbr(1,ind(ia*2+1),j,i3) =  r2 + r1
     zbr(2,ind(ia*2+1),j,i3) =  s2 + s1
     zbr(1,ind(ia*2+2),j,i3) = -r2 + r1
     zbr(2,ind(ia*2+2),j,i3) = -s2 + s1
    enddo
 
!   Second step of factor 2
    do ia=1,aft(ic)-1
     indx=ind(ia*2+1)-1
     cr2=trig(1,indx)
     ct2=trig(2,indx)
     do j=ma,mb
      r1=z(1,ia*2+1,j,i3)
      s1=z(2,ia*2+1,j,i3)
      r2=z(1,ia*2+2,j,i3) - z(2,ia*2+2,j,i3)*ct2
      s2=z(1,ia*2+2,j,i3)*ct2 + z(2,ia*2+2,j,i3)
      zbr(1,ind(ia*2+1),j,i3) =  r2*cr2 + r1
      zbr(2,ind(ia*2+1),j,i3) =  s2*cr2 + s1
      zbr(1,ind(ia*2+2),j,i3) = -r2*cr2 + r1
      zbr(2,ind(ia*2+2),j,i3) = -s2*cr2 + s1
     enddo
    enddo
 
!  Treat factor 3
   else if (now(ic)==3) then  
!   radix 3
!   .5d0*sqrt(3.d0)=0.8660254037844387d0
    ia=0
    bb=ris*0.8660254037844387d0
 
!   First step of factor 3
    do j=ma,mb
     r1=z(1,ia*3+1,j,i3)
     s1=z(2,ia*3+1,j,i3)
     r2=z(1,ia*3+2,j,i3)
     s2=z(2,ia*3+2,j,i3)
     r3=z(1,ia*3+3,j,i3)
     s3=z(2,ia*3+3,j,i3)
     r=r2 + r3
     s=s2 + s3
     zbr(1,ind(ia*3+1),j,i3) = r + r1
     zbr(2,ind(ia*3+1),j,i3) = s + s1
     r1=r1 - r*.5d0
     s1=s1 - s*.5d0
     r2=r2-r3
     s2=s2-s3
     zbr(1,ind(ia*3+2),j,i3) = r1 - s2*bb
     zbr(2,ind(ia*3+2),j,i3) = s1 + r2*bb
     zbr(1,ind(ia*3+3),j,i3) = r1 + s2*bb
     zbr(2,ind(ia*3+3),j,i3) = s1 - r2*bb
    enddo
 
!   Second step of factor 3
    do ia=1,aft(ic)-1
     indx=ind(ia*3+1)-1
     cr2=trig(1,indx)
     ct2=trig(2,indx)
     cr3=trig(1,2*indx)
     ct3=trig(2,2*indx)
     cr2=cr2/cr3
     cr3p=.5d0*cr3
     bb=ris*cr3*0.8660254037844387d0
     do j=ma,mb
      r1=z(1,ia*3+1,j,i3)
      s1=z(2,ia*3+1,j,i3)
      r2=z(1,ia*3+2,j,i3) - z(2,ia*3+2,j,i3)*ct2
      s2=z(1,ia*3+2,j,i3)*ct2 + z(2,ia*3+2,j,i3)
      r3=z(1,ia*3+3,j,i3) - z(2,ia*3+3,j,i3)*ct3
      s3=z(1,ia*3+3,j,i3)*ct3 + z(2,ia*3+3,j,i3)
      r=cr2*r2 + r3
      s=cr2*s2 + s3
      zbr(1,ind(ia*3+1),j,i3) = r*cr3 + r1
      zbr(2,ind(ia*3+1),j,i3) = s*cr3 + s1
      r1=r1 - r*cr3p
      s1=s1 - s*cr3p
      r2=cr2*r2-r3
      s2=cr2*s2-s3
      zbr(1,ind(ia*3+2),j,i3) = r1 - s2*bb
      zbr(2,ind(ia*3+2),j,i3) = s1 + r2*bb
      zbr(1,ind(ia*3+3),j,i3) = r1 + s2*bb
      zbr(2,ind(ia*3+3),j,i3) = s1 - r2*bb
     enddo
    enddo
 
!  Treat factor 5
   else if (now(ic)==5) then 
!   radix 5
!   sin(2.d0*pi/5.d0)
    sin2=ris*0.9510565162951536d0
    ia=0
 
!   First step of factor 5
    do j=ma,mb
     r1=z(1,ia*5+1,j,i3)
     s1=z(2,ia*5+1,j,i3)
     r2=z(1,ia*5+2,j,i3)
     s2=z(2,ia*5+2,j,i3)
     r3=z(1,ia*5+3,j,i3)
     s3=z(2,ia*5+3,j,i3)
     r4=z(1,ia*5+4,j,i3)
     s4=z(2,ia*5+4,j,i3)
     r5=z(1,ia*5+5,j,i3)
     s5=z(2,ia*5+5,j,i3)
     r25 = r2 + r5
     r34 = r3 + r4
     s25 = s2 - s5
     s34 = s3 - s4
     zbr(1,ind(ia*5+1),j,i3) = r1 + r25 + r34
     r = r1 + cos2*r25 + cos4*r34
     s = s25 + sin42*s34
     zbr(1,ind(ia*5+2),j,i3) = r - sin2*s
     zbr(1,ind(ia*5+5),j,i3) = r + sin2*s
     r = r1 + cos4*r25 + cos2*r34
     s = sin42*s25 - s34
     zbr(1,ind(ia*5+3),j,i3) = r - sin2*s
     zbr(1,ind(ia*5+4),j,i3) = r + sin2*s
     r25 = r2 - r5
     r34 = r3 - r4
     s25 = s2 + s5
     s34 = s3 + s4
     zbr(2,ind(ia*5+1),j,i3) = s1 + s25 + s34
     r = s1 + cos2*s25 + cos4*s34
     s = r25 + sin42*r34
     zbr(2,ind(ia*5+2),j,i3) = r + sin2*s
     zbr(2,ind(ia*5+5),j,i3) = r - sin2*s
     r = s1 + cos4*s25 + cos2*s34
     s = sin42*r25 - r34
     zbr(2,ind(ia*5+3),j,i3) = r + sin2*s
     zbr(2,ind(ia*5+4),j,i3) = r - sin2*s
    enddo
 
!   Second step of factor 5
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
     do j=ma,mb
      r1=z(1,ia*5+1,j,i3)
      s1=z(2,ia*5+1,j,i3)
      r2=cr2*(z(1,ia*5+2,j,i3) - z(2,ia*5+2,j,i3)*ct2)
      s2=cr2*(z(1,ia*5+2,j,i3)*ct2 + z(2,ia*5+2,j,i3))
      r3=cr3*(z(1,ia*5+3,j,i3) - z(2,ia*5+3,j,i3)*ct3)
      s3=cr3*(z(1,ia*5+3,j,i3)*ct3 + z(2,ia*5+3,j,i3))
      r4=z(1,ia*5+4,j,i3) - z(2,ia*5+4,j,i3)*ct4
      s4=z(1,ia*5+4,j,i3)*ct4 + z(2,ia*5+4,j,i3)
      r5=z(1,ia*5+5,j,i3) - z(2,ia*5+5,j,i3)*ct5
      s5=z(1,ia*5+5,j,i3)*ct5 + z(2,ia*5+5,j,i3)
      r25 = r2 + r5*cr5
      r34 = r3 + r4*cr4
      s25 = s2 - s5*cr5
      s34 = s3 - s4*cr4
      zbr(1,ind(ia*5+1),j,i3) = r1 + r25 + r34
      r = r1 + cos2*r25 + cos4*r34
      s = s25 + sin42*s34
      zbr(1,ind(ia*5+2),j,i3) = r - sin2*s
      zbr(1,ind(ia*5+5),j,i3) = r + sin2*s
      r = r1 + cos4*r25 + cos2*r34
      s = sin42*s25 - s34
      zbr(1,ind(ia*5+3),j,i3) = r - sin2*s
      zbr(1,ind(ia*5+4),j,i3) = r + sin2*s
      r25 = r2 - r5*cr5
      r34 = r3 - r4*cr4
      s25 = s2 + s5*cr5
      s34 = s3 + s4*cr4
      zbr(2,ind(ia*5+1),j,i3) = s1 + s25 + s34
      r = s1 + cos2*s25 + cos4*s34
      s = r25 + sin42*r34
      zbr(2,ind(ia*5+2),j,i3) = r + sin2*s
      zbr(2,ind(ia*5+5),j,i3) = r - sin2*s
      r = s1 + cos4*s25 + cos2*s34
      s = sin42*r25 - r34
      zbr(2,ind(ia*5+3),j,i3) = r + sin2*s
      zbr(2,ind(ia*5+4),j,i3) = r - sin2*s
     enddo
    enddo

   else
 
!   All factors treated
    write(message, '(a,a,a,a)' )ch10,&
&    ' sg_fftx : BUG -',ch10,&
&    ' called with factors other than 2, 3, and 5'
!    call wrtout(06,message,'PERS')
!    call leave_new('PERS')
   endif
 
!  ---------------------------------------------------------------

  enddo ! do i3=1,n3
 enddo  ! do jj=1,n2,lot
!$OMP END PARALLEL DO

 end subroutine
!!***
