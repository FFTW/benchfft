!{\src2tex{textfont=tt}}
!!****f* ABINIT/sg_fftrisc
!! NAME
!! sg_fftrisc
!!
!! FUNCTION
!! Carry out fourier transforms between real and reciprocal (G) space,
!! for wavefunctions, contained in a sphere in reciprocal space,
!! in both directions. Also accomplish some post-processing.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2001 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .   
!!
!! NOTES
!! Specifically uses rather sophisticated algorithms, based on S Goedecker
!! routines, specialized for superscalar RISC architecture.
!! Zero padding : saves 7/12 execution time
!! Bi-dimensional data locality in most of the routine : cache reuse
!! For k-point (0 0 0) : takes advantage of symmetry of data.
!! Note however that no blocking is used, in both 1D z-transform
!! or subsequent 2D transform. This should be improved.
!!
!! INPUTS
!!  cplex= if 1 , denpot is real, if 2 , denpot is complex
!!     (cplex=2 only allowed for option=2 when istwf_k=1)
!!     one can also use cplex=0 if option=0 or option=3
!!  fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!  gboundin(2*mgfft+4)=sphere boundary info for reciprocal to real space
!!  gboundout(2*mgfft+4)=sphere boundary info for real to reciprocal space
!!  indpw_kin(4,npwin)=fft box index for reciprocal to real space
!!  indpw_kout(4,npwout)=fft box index for real to reciprocal space
!!  istwf_k=option parameter that describes the storage of wfs
!!  npwin=number of elements in fofgin array (for option 0, 1 and 2)
!!  npwout=number of elements in fofgout array (for option 2 and 3)
!!  mgfft=maximum size of 1D FFTs
!!  ngfft(1:3)=integer fft box dimensions
!!  n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!!  ngfft(8)=size of the cache (kB)
!!  option= if 0: do direct FFT
!!          if 1: do direct FFT, then sum the density
!!          if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!          if 3: do reverse FFT only  
!!  weight=weight to be used for the accumulation of the density in real space
!!          (needed only when option=1)    
!!
!! OPTIONS
!!  The different options are:
!!  - reciprocal to real space and output the result (when option=0),
!!  - reciprocal to real space and accumulate the density (when option=1) or
!!  - reciprocal to real space, apply the local potential to the wavefunction
!!    in real space and produce the result in reciprocal space (when option=2)
!!  - real space to reciprocal space (when option=3).
!!  option=0 IS NOT ALLOWED when istwf_k>2     
!!  option=3 IS NOT ALLOWED when istwf_k>=2     
!!
!! SIDE EFFECTS
!!  for option==0, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 fofr(2,n4,n5,n6) contains the Fourier Transform of fofgin;
!!                 no use of denpot, fofgout and npwout.
!!  for option==1, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*n4,n5,n6) contains the input density at input,
!!                 and the updated density at output;
!!                 fofr(2,n4,n5,n6) contains the Fourier transform of fofgin,
!!                 except in the case of the hp library subroutine;
!!                 no use of fofgout and npwout.
!!  for option==2, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                 denpot(cplex*n4,n5,n6) contains the input local potential;
!!                 fofgout(2,npwout) contains the output function;
!!                 fofr(2,n4,n5,n6) contains the Fourier transform of fofgin,
!!                 except in the case of the hp library subroutine.
!!  for option==3, fofr(2,n4,n5,n6) contains the real space wavefunction;
!!                 fofgout(2,npwout) contains its Fourier transform;
!!                 no use of fofgin and npwin.     
!!
!! TODO
!! Complete input and output list.
!!
!! PARENTS
!!      fourwf
!!
!! SOURCE
 
 subroutine sg_fftrisc(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
&  indpw_kin,indpw_kout,istwf_k,mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight)

 use basis_defs
 implicit none                              
 
!Arguments ------------------------------------
 integer :: cplex,istwf_k,mgfft,npwin,npwout,n4,n5,n6,option
 integer :: gboundin(2*mgfft+4),gboundout(2*mgfft+4),&
& indpw_kin(4,npwin),indpw_kout(4,npwout),ngfft(8)
 real(dp) :: weight
 real(dp) :: denpot(cplex*n4,n5,n6),fofgin(2,npwin),&
&  fofgout(2,npwout),fofr(2,n4,n5,n6)

!Local variables-------------------------------
 integer, parameter :: mfac=11,mg=2048
 integer, save :: ic1,ic2,ic3,ic4,ic5,ic6,&
& n1_save=0,n2_save=0,n3_save=0,n4_save=0,n5_save=0,n6_save=0
 integer :: fftcache,g2max,g2min,&
& ier,ifftig,ig,igb,igbmax,igb_inv,ii1,ii2,ii3,index,isign,&
& i1,i1inv,i1max,i2,i3,i3inv,lot,lotin,lotout,&
& mgb,nfft,ngbin,ngbout,nlot,nproc_omp,&
& n1,n1halfm,n1half1,n1i,n2,n2half1,n3,n4half1,n5half1
 integer,save :: &
& aft1(mfac),aft2(mfac),aft3(mfac),aft4(mfac),aft5(mfac),aft6(mfac),&
& bef1(mfac),bef2(mfac),bef3(mfac),bef4(mfac),bef5(mfac),bef6(mfac),&
& ind1(mg)  ,ind2(mg)  ,ind3(mg),  ind4(mg),  ind5(mg),  ind6(mg),&
& now1(mfac),now2(mfac),now3(mfac),now4(mfac),now5(mfac),now6(mfac)
 integer :: gbound_dum(4)
 real(dp) :: ai,ar,fraction,norm,phai,phar,ris,sqr,wkim,wkre
 real(dp), save :: trig1(2,mg),trig2(2,mg),trig3(2,mg),&
& trig4(2,mg),trig5(2,mg),trig6(2,mg),tsec(2)
 real(dp), allocatable :: pha1(:,:),pha2(:,:),pha3(:,:),&
& wk1d_a(:,:,:,:),wk1d_b(:,:,:,:),&
& wk2d_a(:,:,:,:),wk2d_b(:,:,:,:),wk2d_c(:,:,:,:),wk2d_d(:,:,:,:),&
& work1(:,:,:,:)
 character*500 :: message
 
! *************************************************************************
 
! DEBUG
!write(6,*)' sg_fftrisc : enter, istwf_k= ',istwf_k
! stop
! ENDDEBUG
!
 if(istwf_k>2 .and. option==0)then
  write(message, '(a,a,a,a,i4)' )ch10,&
&  ' sg_fftrisc : BUG -',ch10,&
&  '  It is not allowed to use option=0 with istwf_k=',istwf_k,'.'
!  call wrtout(6,message,'PERS')
!  call leave_new('PERS')
 endif

 if(istwf_k>=2 .and. option==3)then
  write(message, '(a,a,a,a,i4)' )ch10,&
&  ' sg_fftrisc : BUG -',ch10,&
&  '  It is not allowed to use option=3 with istwf_k=',istwf_k,'.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 endif

!For all other tests of validity of inputs, assume that they
!have been done in the calling routine

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3) ; nfft=n1*n2*n3
 fftcache=ngfft(8)

!Get the number of 1D transforms along z
 ngbin=maxval(indpw_kin(4,1:npwin))
 ngbout=maxval(indpw_kout(4,1:npwout))
!Define the dimension of the first work arrays, for 1D transforms along z ,
!taking into account the need to avoid the cache trashing
 mgb=max(ngbin,ngbout) 
 if(mod(mgb,2)/=1)mgb=mgb+1

!Initialise openmp, if needed
!$OMP PARALLEL
!$OMP SINGLE
 nproc_omp=1
!$OMP END SINGLE
!$OMP END PARALLEL

!For the treatment of the z transform, 
!one tries to use only a fraction of the cache, since the
!treatment of the array wk1d_a will not involve contiguous segments
 fraction=0.25
!First estimation of lot and nlot
 lot=(fftcache*fraction*1000)/(n3*8*2)+1
!Select the smallest integer multiple of nproc_omp, larger
!or equal to nlot. In this way, the cache size is not exhausted,
!and one takes care correctly of the number of processors.
!Treat separately the in and out cases
 nlot=(ngbin-1)/lot+1
 nlot=nproc_omp*((nlot-1)/nproc_omp+1)
 lotin=(ngbin-1)/nlot+1 
 nlot=(ngbout-1)/lot+1
 nlot=nproc_omp*((nlot-1)/nproc_omp+1)
 lotout=(ngbout-1)/nlot+1
!The next line impose only one lot. Usually, comment it.
!lotin=mgb ; lotout=mgb

!Compute auxiliary arrays needed for FFTs
 if(n1/=n1_save)then
  call sg_ctrig(n1,trig1,aft1,bef1,now1,1.0d0,ic1,ind1,mfac,mg)
  call sg_ctrig(n1,trig4,aft4,bef4,now4,-1.0d0,ic4,ind4,mfac,mg)
  n1_save=n1
 endif
 if(n2/=n2_save)then
  call sg_ctrig(n2,trig2,aft2,bef2,now2,1.0d0,ic2,ind2,mfac,mg)
  call sg_ctrig(n2,trig5,aft5,bef5,now5,-1.0d0,ic5,ind5,mfac,mg)
  n2_save=n2
 endif
 if(n3/=n3_save)then
  call sg_ctrig(n3,trig3,aft3,bef3,now3,1.0d0,ic3,ind3,mfac,mg)
  call sg_ctrig(n3,trig6,aft6,bef6,now6,-1.0d0,ic6,ind6,mfac,mg)
  n3_save=n3
 endif

!------------------------------------------------------------------
!Here, call general k-point code

 if(istwf_k==1)then

! Note that the z transform will appear as a y transform
  allocate (wk1d_a(2,mgb,n3,1),wk1d_b(2,mgb,n3,1))

  if(option/=3)then

!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(n3,ngbin,wk1d_a) 
   do i3=1,n3
    do igb=1,ngbin
     wk1d_a(1,igb,i3,1)=0.0d0
     wk1d_a(2,igb,i3,1)=0.0d0
    enddo
   enddo
!$OMP END PARALLEL DO

!  Insert fofgin into the work array
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(fofgin,indpw_kin,npwin,wk1d_a) 
   do ig=1,npwin
    igb=indpw_kin(4,ig) ; i3=indpw_kin(3,ig)
    wk1d_a(1,igb,i3,1)=fofgin(1,ig)
    wk1d_a(2,igb,i3,1)=fofgin(2,ig)
   enddo
!$OMP END PARALLEL DO

!  Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!  However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3,lotin,mgb)&
!$OMP&SHARED(ngbin,now3,n3,trig3,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
   do igb=1,ngbin,lotin
    igbmax=min(igb+lotin-1,ngbin)
!   Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!   However, due to special packing of data, use routine ffty
    call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_a,wk1d_b, &
&    trig3,aft3,now3,bef3,1.0d0,ind3,ic3)
   enddo
!$OMP END PARALLEL DO

  endif !  if(option/=3)

! Do-loop on the planes stacked in the z direction
!$OMP PARALLEL DEFAULT(PRIVATE) & 
!$OMP&SHARED(aft1,aft2,aft4,aft5,bef1,bef2,bef4,bef5,cplex,denpot) &
!$OMP&SHARED(fftcache,fofr,gboundin,gboundout)&
!$OMP&SHARED(ic1,ic2,ic4,ic5,ind1,ind2,ind4) &
!$OMP&SHARED(ind5,indpw_kin,indpw_kout,mgb,n1,n2,n3,n4,n5,ngbin) &
!$OMP&SHARED(ngbout,now1,now2,now4,now5,option,trig1,trig2,trig4,trig5) &
!$OMP&SHARED(weight,wk1d_a,wk1d_b) 

! Allocate two 2-dimensional work arrays
  allocate(wk2d_a(2,n4,n5,1),wk2d_b(2,n4,n5,1))
!$OMP DO
  do i3=1,n3

   if(option/=3)then
!   Zero the values on the current plane
!   wk2d_a(1:2,1:n1,1:n2,1)=0.0d0
    do i2=1,n2
     do i1=1,n1
      wk2d_a(1,i1,i2,1)=0.0d0
      wk2d_a(2,i1,i2,1)=0.0d0
     enddo
    enddo
!   Copy the data in the current plane
    do igb=1,ngbin
     i1=indpw_kin(1,igb) ; i2=indpw_kin(2,igb)
     wk2d_a(1,i1,i2,1)=wk1d_b(1,igb,i3,1)
     wk2d_a(2,i1,i2,1)=wk1d_b(2,igb,i3,1)
    enddo
!   Perform x transform, taking into account arrays of zeros
    g2min=gboundin(1) ; g2max=gboundin(2)
    if ( g2min+n2 >= g2max+2 ) then
     do i2=g2max+2,g2min+n2
      do i1=1,n1
       wk2d_b(1,i1,i2,1)=0.0d0
       wk2d_b(2,i1,i2,1)=0.0d0
      enddo
     enddo
    endif
    gbound_dum(1)=1 ; gbound_dum(2)=1 
    gbound_dum(3)=g2min ; gbound_dum(4)=g2max
    call sg_fftpx(fftcache,mfac,mg,0,n4,n5,1,n2,1,wk2d_a,wk2d_b,&
&    trig1,aft1,now1,bef1,1.0d0,ind1,ic1,gbound_dum)
!   Perform y transform
    n1i=1
    call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1,1,1,wk2d_b,wk2d_a, &
&    trig2,aft2,now2,bef2,1.0d0,ind2,ic2)
!   The wave function is now in real space, for the current plane
   endif

   if(option==0)then
!   Copy the transformed function at the right place
    do i2=1,n2
     do i1=1,n1
      fofr(1,i1,i2,i3)=wk2d_a(1,i1,i2,1)
      fofr(2,i1,i2,i3)=wk2d_a(2,i1,i2,1)
     enddo
    enddo
   endif

   if(option==1)then
!   Accumulate density
    do i2=1,n2
     do i1=1,n1
      denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&            weight*(wk2d_a(1,i1,i2,1)**2+wk2d_a(2,i1,i2,1)**2)
     enddo
    enddo
   endif

   if(option==2)then
!   Apply local potential
    if(cplex==1)then
     do i2=1,n2
      do i1=1,n1
       wk2d_a(1,i1,i2,1)=denpot(i1,i2,i3)*wk2d_a(1,i1,i2,1)
       wk2d_a(2,i1,i2,1)=denpot(i1,i2,i3)*wk2d_a(2,i1,i2,1)
      enddo
     enddo
    else
     do i2=1,n2
      do i1=1,n1
       wkre=wk2d_a(1,i1,i2,1)
       wkim=wk2d_a(2,i1,i2,1)
       wk2d_a(1,i1,i2,1)=denpot(2*i1-1,i2,i3)*wkre &
&                       -denpot(2*i1  ,i2,i3)*wkim
       wk2d_a(2,i1,i2,1)=denpot(2*i1-1,i2,i3)*wkim &
&                       +denpot(2*i1  ,i2,i3)*wkre
      enddo
     enddo
    endif
   endif

   if(option==3)then
!   Copy the function to be tranformed at the right place
    do i2=1,n2
     do i1=1,n1
      wk2d_a(1,i1,i2,1)=fofr(1,i1,i2,i3)
      wk2d_a(2,i1,i2,1)=fofr(2,i1,i2,i3)
     enddo
    enddo
   endif

   if(option==2 .or. option==3)then
!   Perform y transform
    n1i=1
    call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1,1,1,wk2d_a,wk2d_b, &
&    trig5,aft5,now5,bef5,-1.0d0,ind5,ic5)
!   Perform x transform, taking into account arrays of zeros
    gbound_dum(1)=1 ; gbound_dum(2)=1
    gbound_dum(3)=gboundout(1) ; gbound_dum(4)=gboundout(2)
    call sg_fftpx(fftcache,mfac,mg,0,n4,n5,1,n2,1,wk2d_b,wk2d_a,&
&   trig4,aft4,now4,bef4,-1.0d0,ind4,ic4,gbound_dum)
!   Copy the data from the current plane to wk1d_b
    do igb=1,ngbout
     i1=indpw_kout(1,igb) ; i2=indpw_kout(2,igb)
     wk1d_b(1,igb,i3,1)=wk2d_a(1,i1,i2,1)
     wk1d_b(2,igb,i3,1)=wk2d_a(2,i1,i2,1)
    enddo
   endif

! End loop on planes
  enddo
!$OMP END DO
  deallocate(wk2d_a,wk2d_b)
!$OMP END PARALLEL

  if(option==2 .or. option==3)then

!  Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!  However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft6,bef6,fftcache,ind6,ic6,lotout,mgb)&
!$OMP&SHARED(ngbout,now6,n3,trig6,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
   do igb=1,ngbout,lotout
    igbmax=min(igb+lotout-1,ngbout)
!   Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!   However, due to special packing of data, use routine ffty
    call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_b,wk1d_a, &
&    trig6,aft6,now6,bef6,-1.0d0,ind6,ic6)

   enddo
!$OMP END PARALLEL DO

!  Transfer the data in the output array, after normalization
   norm=1.d0/dble(nfft)
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(fofgout,indpw_kout,norm,npwout,wk1d_a) 
   do ig=1,npwout
    igb=indpw_kout(4,ig) ; i3=indpw_kout(3,ig)
    fofgout(1,ig)=wk1d_a(1,igb,i3,1)*norm
    fofgout(2,ig)=wk1d_a(2,igb,i3,1)*norm
   enddo
!$OMP END PARALLEL DO
  endif

  deallocate(wk1d_a,wk1d_b)

!End general k-point part
 endif

!------------------------------------------------------------------
!Here, use of time-reversal symmetry

 if(istwf_k>=2)then

  n1half1=n1/2+1 ; n1halfm=(n1+1)/2
  n2half1=n2/2+1
! n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
  n4half1=(n1half1/2)*2+1
  n5half1=(n2half1/2)*2+1
! Note that the z transform will appear as a y transform
  allocate (wk1d_a(2,mgb,n3,1),wk1d_b(2,mgb,n3,1))
  
  if(istwf_k/=2)then
   allocate(pha1(2,n1),pha2(2,n2),pha3(3,n3))
   do i1=1,n1
    pha1(1,i1)=cos(dble(i1-1)*pi/dble(n1))
    pha1(2,i1)=sin(dble(i1-1)*pi/dble(n1))
   enddo
   do i2=1,n2
    pha2(1,i2)=cos(dble(i2-1)*pi/dble(n2))
    pha2(2,i2)=sin(dble(i2-1)*pi/dble(n2))
   enddo
   do i3=1,n3
    pha3(1,i3)=cos(dble(i3-1)*pi/dble(n3))
    pha3(2,i3)=sin(dble(i3-1)*pi/dble(n3))
   enddo
  endif

  if(option/=3)then

!  Zero the components of wk1d_a
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(n3,ngbin,wk1d_a) 
   do i3=1,n3
    do igb=1,ngbin
     wk1d_a(1,igb,i3,1)=0.0d0
     wk1d_a(2,igb,i3,1)=0.0d0
    enddo
   enddo
!$OMP END PARALLEL DO

!  Insert fofgin into the work array
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(fofgin,indpw_kin,npwin,wk1d_a) 
   do ig=1,npwin
    igb=indpw_kin(4,ig) ; i3=indpw_kin(3,ig)
    wk1d_a(1,igb,i3,1)=fofgin(1,ig)
    wk1d_a(2,igb,i3,1)=fofgin(2,ig)
   enddo
!$OMP END PARALLEL DO

!  Must complete the i2=1 plane when $k_y \equiv 0$

!  Take care of i1=1 when $k_x \equiv 0$
   if(istwf_k==2)then
!   Take care of i1=1
    do i3=n3/2+1,n3
     i3inv=n3+2-i3
     wk1d_a(1,1,i3,1)= wk1d_a(1,1,i3inv,1)
     wk1d_a(2,1,i3,1)=-wk1d_a(2,1,i3inv,1)
    enddo
   else if(istwf_k==4)then
!   Take care of i1=1
    do i3=n3/2+1,n3
     i3inv=n3+1-i3
     wk1d_a(1,1,i3,1)= wk1d_a(1,1,i3inv,1)
     wk1d_a(2,1,i3,1)=-wk1d_a(2,1,i3inv,1)
    enddo
   endif

!  Now, take care of other i1 values, except i3==1 when $k_z \equiv 0$ 
   i1max=gboundin(4)+1
   if(istwf_k==2)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(i1max,n3,wk1d_a) 
    do igb=2,2*i1max-1
     igb_inv=2*i1max+1-igb
     do i3=n3/2+1,n3
      i3inv=n3+2-i3 
      wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
      wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
     enddo
    enddo
!$OMP END PARALLEL DO

   else if(istwf_k==3)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(i1max,n3,wk1d_a) 
    do igb=1,2*i1max
     igb_inv=2*i1max+1-igb
     do i3=n3/2+1,n3
      i3inv=n3+2-i3
      wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
      wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
     enddo
    enddo
!$OMP END PARALLEL DO

   else if(istwf_k==4)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(i1max,n3,wk1d_a) 
    do igb=2,2*i1max-1
     igb_inv=2*i1max+1-igb
     do i3=n3/2+1,n3
      i3inv=n3+1-i3
      wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
      wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
     enddo
    enddo
!$OMP END PARALLEL DO

   else if(istwf_k==5)then
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(i1max,n3,wk1d_a) 
    do igb=1,2*i1max
     igb_inv=2*i1max+1-igb
     do i3=n3/2+1,n3
      i3inv=n3+1-i3
      wk1d_a(1,igb,i3,1)= wk1d_a(1,igb_inv,i3inv,1)
      wk1d_a(2,igb,i3,1)=-wk1d_a(2,igb_inv,i3inv,1)
     enddo
    enddo
!$OMP END PARALLEL DO

   endif

!  Now, i3==1
   if(istwf_k==2)then
    do igb=2,i1max
     igb_inv=2*i1max+1-igb
     wk1d_a(1,igb_inv,1,1)= wk1d_a(1,igb,1,1)
     wk1d_a(2,igb_inv,1,1)=-wk1d_a(2,igb,1,1)
    enddo
   else if(istwf_k==3)then
    do igb=1,i1max
     igb_inv=2*i1max+1-igb
     wk1d_a(1,igb_inv,1,1)= wk1d_a(1,igb,1,1)
     wk1d_a(2,igb_inv,1,1)=-wk1d_a(2,igb,1,1)
    enddo
   endif

!  Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!  However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3,lotin,mgb)&
!$OMP&SHARED(ngbin,now3,n3,trig3,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
   do igb=1,ngbin,lotin
    igbmax=min(igb+lotin-1,ngbin)
!   Go from wk1d_a to wk1d_b, using 1D FFTs on the z direction
!   However, due to special packing of data, use routine ffty
    call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_a,wk1d_b, &
&    trig3,aft3,now3,bef3,1.0d0,ind3,ic3)
   enddo
!$OMP END PARALLEL DO


!  Change the phase if $k_z \neq 0$
   if(istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(ngbin,n3,pha3,wk1d_b) 
    do i3=1,n3
     phar=pha3(1,i3)
     phai=pha3(2,i3)
     do igb=1,ngbin
      ar=wk1d_b(1,igb,i3,1)
      ai=wk1d_b(2,igb,i3,1)
      wk1d_b(1,igb,i3,1)=phar*ar-phai*ai
      wk1d_b(2,igb,i3,1)=phai*ar+phar*ai
     enddo
    enddo
!$OMP END PARALLEL DO
   endif

  endif !  if(option/=3)

! Do-loop on the planes stacked in the z direction

!$OMP PARALLEL DEFAULT(PRIVATE) & 
!$OMP&SHARED(aft1,aft2,aft4,aft5,bef1,bef2,bef4,bef5,denpot) &
!$OMP&SHARED(fftcache,fofr,gboundin,ic1,ic2,ic4,ic5,ind1,ind2,ind4,ind5) &
!$OMP&SHARED(indpw_kin,indpw_kout,istwf_k,mgb,n1,n1half1) &
!$OMP&SHARED(n1halfm,n2,n2half1,n3,n4,n5,ngbin,ngbout) &
!$OMP&SHARED(now1,now2,now4,now5,option,pha1,pha2,trig1) &
!$OMP&SHARED(trig2,trig4,trig5,weight,wk1d_a,wk1d_b) 

! Allocate two 2-dimensional work arrays
  allocate(wk2d_a(2,n4,n5,1),wk2d_b(2,n4,n5,1))
  allocate(wk2d_c(2,2*n1halfm,n5,1),wk2d_d(2,2*n1halfm,n5,1))
!$OMP DO
  do i3=1,n3

!DEBUG
!  write(6,*)' sg_fftrisc : i3=',i3
!ENDDEBUG

   g2max=gboundin(2)

   if(option/=3)then
!   Zero the values on the current plane : need only from i2=1 to g2max+1
    do i2=1,g2max+1
     do i1=1,n1
      wk2d_a(1,i1,i2,1)=0.0d0
      wk2d_a(2,i1,i2,1)=0.0d0
     enddo
    enddo

!   Copy the data in the current plane
    do igb=1,ngbin
     i1=indpw_kin(1,igb) ; i2=indpw_kin(2,igb)
     wk2d_a(1,i1,i2,1)=wk1d_b(1,igb,i3,1)
     wk2d_a(2,i1,i2,1)=wk1d_b(2,igb,i3,1)
    enddo

!   Perform x transform, taking into account arrays of zeros
    call sg_fftx(fftcache,mfac,mg,n4,n5,1,g2max+1,1,wk2d_a,wk2d_b,&
&     trig1,aft1,now1,bef1,1.0d0,ind1,ic1)

!   Change the phase if $k_x \neq 0$
    if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9)then
     do i1=1,n1
      phar=pha1(1,i1)
      phai=pha1(2,i1)
      do i2=1,g2max+1
       ar=wk2d_b(1,i1,i2,1)
       ai=wk2d_b(2,i1,i2,1)
       wk2d_b(1,i1,i2,1)=phar*ar-phai*ai
       wk2d_b(2,i1,i2,1)=phai*ar+phar*ai
      enddo
     enddo
    endif

!   Compute symmetric and antisymmetric combinations
    if(istwf_k>=2 .and. istwf_k<=5)then
     do i1=1,n1half1-1
      wk2d_a(1,i1,1,1)=wk2d_b(1,2*i1-1,1,1)
      wk2d_a(2,i1,1,1)=wk2d_b(1,2*i1  ,1,1)
     enddo
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      wk2d_a(1,n1half1,1,1)=wk2d_b(1,n1,1,1)
      wk2d_a(2,n1half1,1,1)=0.0d0
     endif
     ii2=2
    else 
     ii2=1
    endif
    if( g2max+1 >= ii2)then
     do i2=ii2,g2max+1
      do i1=1,n1half1-1
       wk2d_a(1,i1,i2,1)=        wk2d_b(1,2*i1-1,i2,1)-wk2d_b(2,2*i1,i2,1)
       wk2d_a(2,i1,i2,1)=        wk2d_b(2,2*i1-1,i2,1)+wk2d_b(1,2*i1,i2,1)
       wk2d_a(1,i1,n2+ii2-i2,1)= wk2d_b(1,2*i1-1,i2,1)+wk2d_b(2,2*i1,i2,1)
       wk2d_a(2,i1,n2+ii2-i2,1)=-wk2d_b(2,2*i1-1,i2,1)+wk2d_b(1,2*i1,i2,1)
      enddo
      if((2*n1half1-2)/=n1)then
       wk2d_a(1,n1half1,i2,1)=        wk2d_b(1,n1,i2,1)
       wk2d_a(2,n1half1,i2,1)=        wk2d_b(2,n1,i2,1)
       wk2d_a(1,n1half1,n2+ii2-i2,1)= wk2d_b(1,n1,i2,1)
       wk2d_a(2,n1half1,n2+ii2-i2,1)=-wk2d_b(2,n1,i2,1)
      endif
     enddo
    endif
    if ( n2half1 >= g2max+2 ) then
     do i2=g2max+2,n2half1
      do i1=1,n1half1-1
       wk2d_a(1,i1,i2,1)=0.0d0
       wk2d_a(2,i1,i2,1)=0.0d0
       wk2d_a(1,i1,n2+ii2-i2,1)=0.0d0
       wk2d_a(2,i1,n2+ii2-i2,1)=0.0d0
      enddo
      if((2*n1half1-2)/=n1)then
       wk2d_a(1,n1half1,i2,1)=0.0d0
       wk2d_a(2,n1half1,i2,1)=0.0d0
       wk2d_a(1,n1half1,n2+ii2-i2,1)=0.0d0
       wk2d_a(2,n1half1,n2+ii2-i2,1)=0.0d0
      endif
     enddo
    endif

    n1i=1
    call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1halfm,1,1,wk2d_a,wk2d_b,&
&    trig2,aft2,now2,bef2,1.0d0,ind2,ic2)

!   Change the phase if $k_y \neq 0$
    if(istwf_k>=6 .and. istwf_k<=9)then
     do i2=1,n2
      phar=pha2(1,i2)
      phai=pha2(2,i2)
      do i1=1,n1halfm
       ar=wk2d_b(1,i1,i2,1)
       ai=wk2d_b(2,i1,i2,1)
       wk2d_b(1,i1,i2,1)= phar*ar-phai*ai
       wk2d_b(2,i1,i2,1)= phai*ar+phar*ai
      enddo
     enddo
    endif
 
   endif ! option/=3

!  The wave function is now in real space, for the current plane,
!  represented by REAL numbers, although packed in the complex array wk2d_b

   g2max=gboundin(2)

   if(option==0)then
!   This option is only permitted for istwf_k==2 (Gamma point)
!   Copy the transformed function at the right place
    do i2=1,n2
     do i1=1,n1half1-1
      fofr(1,2*i1-1,i2,i3)=wk2d_b(1,i1,i2,1)
      fofr(1,2*i1  ,i2,i3)=wk2d_b(2,i1,i2,1)
      fofr(2,2*i1-1,i2,i3)=0.0d0
      fofr(2,2*i1  ,i2,i3)=0.0d0
     enddo
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      fofr(1,n1,i2,i3)=wk2d_b(1,n1half1,i2,1)
      fofr(2,n1,i2,i3)=0.0d0
     endif
    enddo
   endif

   if(option==1)then
!   Accumulate density
    do i2=1,n2
     do i1=1,n1half1-1
      denpot(2*i1-1,i2,i3)=denpot(2*i1-1,i2,i3)+weight*wk2d_b(1,i1,i2,1)**2
      denpot(2*i1  ,i2,i3)=denpot(2*i1  ,i2,i3)+weight*wk2d_b(2,i1,i2,1)**2
     enddo
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      denpot(n1,i2,i3)=denpot(n1,i2,i3)+weight*wk2d_b(1,n1half1,i2,1)**2
     endif
    enddo
   endif

   if(option==2)then
!   Apply local potential
    do i2=1,n2
     do i1=1,n1half1-1
      wk2d_a(1,i1,i2,1)=denpot(2*i1-1,i2,i3)*wk2d_b(1,i1,i2,1)
      wk2d_a(2,i1,i2,1)=denpot(2*i1  ,i2,i3)*wk2d_b(2,i1,i2,1)
     enddo
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      wk2d_a(1,n1half1,i2,1)=denpot(n1,i2,i3)*wk2d_b(1,n1half1,i2,1)
      wk2d_a(2,n1half1,i2,1)=0.0d0
     endif
    enddo
   endif

   if(option==3)then
!   This option is only permitted for istwf_k==2 (Gamma point)
!   Copy the transformed function at the right place
    do i2=1,n2
     do i1=1,n1half1-1
      wk2d_b(1,i1,i2,1)=fofr(1,2*i1-1,i2,i3)
      wk2d_b(2,i1,i2,1)=fofr(1,2*i1  ,i2,i3)
     enddo
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      wk2d_b(1,n1half1,i2,1)=fofr(1,n1,i2,i3)
     endif
    enddo
   endif

   if(option==2 .or. option==3)then
!   Change the phase if $k_y \neq 0$
    if(istwf_k>=6 .and. istwf_k<=9)then
     do i2=1,n2
      phar=pha2(1,i2)
      phai=pha2(2,i2)
      do i1=1,n1halfm
       ar=wk2d_a(1,i1,i2,1)
       ai=wk2d_a(2,i1,i2,1)
       wk2d_a(1,i1,i2,1)= phar*ar+phai*ai
       wk2d_a(2,i1,i2,1)=-phai*ar+phar*ai
      enddo
     enddo
    endif

!   Perform y transform
    n1i=1
    call sg_ffty(fftcache,mfac,mg,n4,n5,1,n1i,n1halfm,1,1,wk2d_a,wk2d_b, &
&    trig5,aft5,now5,bef5,-1.0d0,ind5,ic5)

!   Decompose symmetric and antisymmetric parts
    if(istwf_k>=2 .and. istwf_k<=5)then
     do i1=1,n1halfm
      wk2d_c(1,2*i1-1,1,1)=wk2d_b(1,i1,1,1)
      wk2d_c(2,2*i1-1,1,1)=0.0d0
      wk2d_c(1,2*i1,1,1)=wk2d_b(2,i1,1,1)
      wk2d_c(2,2*i1,1,1)=0.0d0
     enddo
     ii2=2
    else 
     ii2=1
    endif
    do i2=ii2,g2max+1
     do i1=1,n1halfm
      wk2d_c(1,2*i1-1,i2,1)=(wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+ii2-i2,1))*0.5d0
      wk2d_c(2,2*i1-1,i2,1)=(wk2d_b(2,i1,i2,1)-wk2d_b(2,i1,n2+ii2-i2,1))*0.5d0
      wk2d_c(1,2*i1,i2,1)= ( wk2d_b(2,i1,i2,1)+wk2d_b(2,i1,n2+ii2-i2,1))*0.5d0
      wk2d_c(2,2*i1,i2,1)= (-wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+ii2-i2,1))*0.5d0
     enddo
    enddo

!   Change the phase if $k_x \neq 0$
    if(istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9 )then
     do i1=1,n1
      phar=pha1(1,i1)
      phai=pha1(2,i1)
      do i2=1,g2max+1
       ar=wk2d_c(1,i1,i2,1)
       ai=wk2d_c(2,i1,i2,1)
       wk2d_c(1,i1,i2,1)= phar*ar+phai*ai
       wk2d_c(2,i1,i2,1)=-phai*ar+phar*ai
      enddo
     enddo
    endif

!   Perform x transform : for y=1 to g2max+1, to benefit from zeros
    call sg_fftx(fftcache,mfac,mg,2*n1halfm,n5,1,g2max+1,1,wk2d_c,wk2d_d,&
&    trig4,aft4,now4,bef4,-1.0d0,ind4,ic4)

!   Copy the data from the current plane to wk1d_b
    do igb=1,ngbout
     i1=indpw_kout(1,igb) ; i2=indpw_kout(2,igb)
     wk1d_b(1,igb,i3,1)=wk2d_d(1,i1,i2,1)
     wk1d_b(2,igb,i3,1)=wk2d_d(2,i1,i2,1)
    enddo

   endif ! option==2 or 3

! End loop on planes
  enddo
!$OMP END DO
 deallocate(wk2d_a,wk2d_b,wk2d_c,wk2d_d)
!$OMP END PARALLEL

  if(option==2 .or. option==3)then

!  Change the phase if $k_z \neq 0$
   if(istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(ngbout,n3,pha3,wk1d_b) 
    do i3=1,n3
     phar=pha3(1,i3)
     phai=pha3(2,i3)
     do igb=1,ngbout
      ar=wk1d_b(1,igb,i3,1)
      ai=wk1d_b(2,igb,i3,1)
      wk1d_b(1,igb,i3,1)= phar*ar+phai*ai
      wk1d_b(2,igb,i3,1)=-phai*ar+phar*ai
     enddo
    enddo
!$OMP END PARALLEL DO
   endif

!  Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!  However, due to special packing of data, use routine ffty
!$OMP PARALLEL DO SHARED(aft6,bef6,fftcache,ind6,ic6,lotout,mgb)&
!$OMP&SHARED(ngbout,now6,n3,trig6,wk1d_a,wk1d_b)&
!$OMP&PRIVATE(igb,igbmax)
   do igb=1,ngbout,lotout
    igbmax=min(igb+lotout-1,ngbout)
!   Go from wk1d_b to wk1d_a, using 1D FFTs on the z direction
!   However, due to special packing of data, use routine ffty
    call sg_ffty(fftcache,mfac,mg,mgb,n3,1,igb,igbmax,1,1,wk1d_b,wk1d_a, &
&    trig6,aft6,now6,bef6,-1.0d0,ind6,ic6)

   enddo
!$OMP END PARALLEL DO

!  Transfer the data in the output array, after normalization
   norm=1.d0/dble(nfft)
!$OMP PARALLEL DO DEFAULT(PRIVATE) & 
!$OMP&SHARED(fofgout,indpw_kout,norm,npwout,wk1d_a) 
   do ig=1,npwout
    igb=indpw_kout(4,ig) ; i3=indpw_kout(3,ig)
    fofgout(1,ig)=wk1d_a(1,igb,i3,1)*norm
    fofgout(2,ig)=wk1d_a(2,igb,i3,1)*norm
   enddo
!$OMP END PARALLEL DO

  endif

  deallocate(wk1d_a,wk1d_b)

  if(istwf_k/=2)then
   deallocate(pha1,pha2,pha3)
  endif

!End time-reversal symmetry
 endif

!------------------------------------------------------------------

! DEBUG
! write(6,*)' sg_fftrisc : exit '
! stop
! ENDDEBUG

 end subroutine
!!***
