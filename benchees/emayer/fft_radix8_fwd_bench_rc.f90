!   Copyright 1997 by Ernst W. Mayer. This program may be used and
!   redistributed freely as long as this header is included. You may
!   modify this program freely, as long as any redistribution
!   contains the original header and a summary of any changes made.
!
!   This software is offered "as is," without warranty of any kind.
!   It is intended only for the private recreational use of our
!   audience. If it causes your CPU to go "poof," that's tough.
!
!   Any modifications of this program intended for distribution MUST
!   be put into the public domain, or provided freely to the original
!   author. Under no circumstances is this program, or any program
!   derived from it or making use of it in any way, to be used for
!   commercial activities of any kind without the express written consent
!   of the original author. (Read: you aren't allowed to make money off
!   it without my say-so.)
!
!   If you use this program or any program derived from it in work
!   leading to a publication, proper acknowledgement should be made in
!   said publication.
!
!   The author appreciates if users send bug reports, inquiries related
!   to the software, or reports of useful modifications or novel and
!   interesting applications to him at the following address:
!
!       Ernst Mayer
!       Dept. of Mechanical and Aerospace Engineering
!       Case Western Reserve University
!       10900 Euclid Avenue
!       Cleveland, OH 44106-7222 USA
!       E-mail: mayer@nigel.mae.cwru.edu
!           or  ewmayer@aol.com
!***************

!...
!   Subroutine which replaces a real vector A of length N by its forward DFT.
!   
!   Since A is real, only the positive-frequency half of the transform is
!   computed, owing to the conjugate-symmetry of the coefficients. The real
!   coefficients of the transform, F(1) and F(N/2+1), are stored as the real
!   and imaginary parts of the first complex element of the transformed array.
!   For details regarding this kind of storage scheme, cf. "Numerical Recipes
!   (in Fortran or C): the Art of Scientific Computation," 2nd Edition,
!   by Press, Flannery, Teukolsky and Vetterling (Cambridge Univ. Press).
!   
!   NOTE: N must be an integer power of 2 - the routine does not check this,
!         so the user must.
!   
!   Calls DLANCZOS_FWD.
!   
!	This version dated 27 February 1997.
!
!...
	subroutine dfft_fwd(a,n)
!...
	implicit none
	integer, intent(in) :: n
!* look for an extended precision real data type for sincos inits; if none available, use real*8
	integer, parameter :: r16=max(selected_real_kind(18,50),selected_real_kind(15,30))
	integer :: i,i1,i2,i3,i4,np3
	integer, save :: n2,n4
	real*8, intent(inout) :: a(n)
	real*8 :: h1r,h1i,h2r,h2i,t1,t2,t3,t4
	real*8, allocatable, save :: w1(:,:)
	real(kind=r16), save :: mt,one=1,pi,theta
	logical, save :: first_entry=.true.

!...initialize things upon first entry
	if(first_entry) then
	  first_entry=.false.
	  n2=n/2
	  n4=n/4
	  allocate(w1(2,n/4))
	  pi=acos(-one)
	  theta=pi/n2	!* = 2*pi/n
	  do i=2,n4	!* exp(2*pi*i*k/n), k=0, ... , n/4-1 are here
	    mt=(i-1)*theta
	    w1(1,i)=cos(mt)
	    w1(2,i)=sin(mt)
	  enddo
	endif

!...enter transform. Result returned in A.

	call dlanczos_fwd(a,n)

!...the rest is to unscramble the data...
	np3=n+3
	do i=2,n4
	  i1=2*i-1
	  i2=i1+1
	  i3=np3-i2
	  i4=i3+1
!...
	  t1=a(i1)
	  t2=a(i2)
	  t3=a(i3)
	  t4=a(i4)
!...Separate the two transforms
	  h1r= .5d0*(t1+t3)
	  h1i= .5d0*(t2-t4)
	  h2r= .5d0*(t2+t4)
	  h2i= .5d0*(t3-t1)
!...to obtain the true transform of the original real data...
	  t1=w1(1,i)
	  t2=w1(2,i)
	  t3=t1*h2r-t2*h2i
	  t4=t2*h2r+t1*h2i
	  a(i1)= h1r+t3
	  a(i2)= h1i+t4
	  a(i3)= h1r-t3
	  a(i4)= t4-h1i
	enddo
!...squeeze first and last datum together...
	h1r=a(1)
	a(1)=h1r+a(2)
	a(2)=h1r-a(2)
!...
	end subroutine dfft_fwd

!***************

	subroutine dlanczos_fwd(a,n)
!...
!   Subroutine which replaces a complex vector A of length N/2 by its forward DFT.
	implicit none
	integer, intent(in) :: n
!* if your compiler does not support real(kind=r16), use real*8 and keep an eye out for error warnings
	integer, parameter :: r16=max(selected_real_kind(18,50),selected_real_kind(15,30))
	integer :: i,ilo,ihi,incr,j,j1,j2,j3,j4,j5,j6,j7,j8,k,m,mm,nh,n2bit
	integer, save :: n2,n4,n8
	integer, allocatable, save :: index(:)
	real*8, intent(inout) :: a(n)
	real*8, allocatable, save :: b(:),w2(:,:),w3(:,:)
	real*8, save :: isqrt2
	real*8 :: tr1,tr2,ti1,ti2,tt,c1,s1,c2,s2,c3,s3
	real*8 :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
	real(kind=r16) :: mt,one=1,two=2,pi,theta,theta2,theta4,mt2,mt4
	logical, save :: first_entry=.true.

!...initialize things upon first entry
	if(first_entry) then
	  first_entry=.false.
	  allocate(b(n),index(n))
	  pi=acos(-one)
	  isqrt2=one/sqrt(two)
	  n2=n/2; n2bit=nint(log(1d0*n2)/log(2d0))
	  n4=n/4
	  n8=n/8
!...bit-reversal array is here...
	  do i=1,n
	    index(i)=i
	  enddo
	  j=1
	  do i=1,n,2
	    if(j > i) then
	      k=index(j)
	      index(j)=index(i)
	      index(i)=k
	      k=index(j+1)
	      index(j+1)=index(i+1)
	      index(i+1)=k
	    endif
	    nh=n/2
	    do
	    if( nh < 2 .or. j < nh ) EXIT
	      j=j-nh
	      nh=nh/2
	    enddo
	    j=j+nh
	  enddo
!...trig array is here: first, calculate the needed dimension...
	  j=1
	  mm=1
	  do
	  if( n2 <= mm ) EXIT
	    mm=8*mm
	    j=j+mm
	  enddo
	  j=3*j+1
	  allocate(w2(2,j))
!...then calculate the sincos data. First part is for the radix-8 loops.
	  j=1
	  mm=1
	  do
	  if( n2 <= mm ) EXIT
	    theta =pi/mm
	    theta2=theta /2
	    theta4=theta2/2
	    w2(1,j:j+2)=1d0
	    w2(2,j:j+2)=0d0
	    do m=1,mm
	      mt =m*theta
	      mt2=m*theta2
	      mt4=m*theta4
	      j=j+3
	      w2(1,j  )=cos(mt)
	      w2(2,j  )=sin(mt)
	      w2(1,j+1)=cos(mt2)
	      w2(2,j+1)=sin(mt2)
	      w2(1,j+2)=cos(mt4)
	      w2(2,j+2)=sin(mt4)
	    enddo
	    mm=8*mm
	  enddo
!...Now take care of the radix-2 and radix-4 wrappers.
!...n2 == 1 (mod 3) sincos data are here...
	  if(mod(n2bit,3)==1)then
	    allocate(w3(2,n4))
	    theta =pi/n4
	    w3(1,1)=1d0
	    w3(2,1)=0d0
	    do m=1,n4-1
	      mt =m*theta
	      w3(1,m+1)=cos(mt)
	      w3(2,m+1)=sin(mt)
	    enddo
	  elseif(mod(n2bit,3)==2)then
!...n2 == 2 (mod 3) sincos data are here.
	    allocate(w3(2,n2))
	    theta =pi/n8
	    theta2=theta/2
	    w3(1,1:2)=1d0
	    w3(2,1:2)=0d0
	    do m=1,n4-1
	      mt =m*theta
	      mt2=m*theta2
	      j=2*m+1
	      w3(1,j  )=cos(mt)
	      w3(2,j  )=sin(mt)
	      w3(1,j+1)=cos(mt2)
	      w3(2,j+1)=sin(mt2)
	    enddo
	  endif
	endif

!...Outer loop, which is executed log_2(n2)/3 times, is unrolled in this implementation.
!   We bit-reverse the vector on the first pass triplet. The pass numbering
!   reflects the passes which would be performed in a standard radix-2 algorithm.
!...Passes 1-3...
	  incr=n
	  ilo=1; ihi=incr
	    do i=ilo,ihi,16
!	get indices are here
	      j1=index(i  )
	      j2=index(i+2)
	      j3=index(i+4)
	      j4=index(i+6)
	      j5=index(i+8)
	      j6=index(i+10)
	      j7=index(i+12)
	      j8=index(i+14)
!	gather the needed data (8 64-byte complex, i.e. 16 64-byte reals)...
	      t1 =a(j1  )
	      t2 =a(j1+1)
	      t3 =a(j2)
	      t4 =a(j2+1)
	      t5 =a(j3)
	      t6 =a(j3+1)
	      t7 =a(j4)
	      t8 =a(j4+1)
	      t9 =a(j5)
	      t10=a(j5+1)
	      t11=a(j6)
	      t12=a(j6+1)
	      t13=a(j7)
	      t14=a(j7+1)
	      t15=a(j8)
	      t16=a(j8+1)
!	first get the 4 length-2 transforms...
	      tt=t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tt =t11
	      t11=t9 -tt
	      t9 =t9 +tt
	      tt =t12
	      t12=t10-tt
	      t10=t10+tt
	      tt =t15
	      t15=t13-tt
	      t13=t13+tt
	      tt =t16
	      t16=t14-tt
	      t14=t14+tt
!	combine to get the 2 length-4 transform...
	      tr1=t5
	      ti1=t6
	      tr2=t7
	      ti2=t8
	      t5=t1-tr1
	      t1=t1+tr1
	      t6=t2-ti1
	      t2=t2+ti1
	      t7=t3+ti2
	      t3=t3-ti2
	      t8=t4-tr2
	      t4=t4+tr2

	      tr1=t13
	      ti1=t14
	      tr2=t15
	      ti2=t16
	      t13=t9 -tr1
	      t9 =t9 +tr1
	      t14=t10-ti1
	      t10=t10+ti1
	      t15=t11+ti2
	      t11=t11-ti2
	      t16=t12-tr2
	      t12=t12+tr2
!	put indices are here
	      j1=ishft(i,-3)+1
	      j2=j1+n8
	      j3=j2+n8
	      j4=j3+n8
	      j5=j4+n8
	      j6=j5+n8
	      j7=j6+n8
	      j8=j7+n8
!	now combine the two half-transforms.
	      b(j1  )=t1+t9
	      b(j1+1)=t2+t10
	      b(j5  )=t1-t9
	      b(j5+1)=t2-t10

	      tr1=(t11-t12)*isqrt2	!* mpy by (1+i)/sqrt(2) is here...
	      ti1=(t11+t12)*isqrt2
	      b(j2  )=t3+tr1
	      b(j2+1)=t4+ti1
	      b(j6  )=t3-tr1
	      b(j6+1)=t4-ti1

	      b(j3  )=t5-t14		!* mpy by i is inlined here...
	      b(j3+1)=t6+t13
	      b(j7  )=t5+t14
	      b(j7+1)=t6-t13

	      tr1=(t15+t16)*isqrt2	!* mpy by (1-i)/sqrt(2) is here...
	      ti1=(t16-t15)*isqrt2
	      b(j4  )=t7-tr1		!* and get (i-1)/sqrt by flipping signs here.
	      b(j4+1)=t8-ti1
	      b(j8  )=t7+tr1
	      b(j8+1)=t8+ti1
	    enddo
	if(n2==8)then; a=b; RETURN; endif
!...Pass 4...
	incr=incr/8
	if(n2==16)then
	  ilo=1; ihi=incr
	  do m=1,8
	    c1=w3(1,m)
	    s1=w3(2,m)
	    do i=ilo,ihi,4
	      j1=ishft(i,-1)+1
	      j2=j1+n2
	      t1=b(i  )
	      t2=b(i+1)
	      t3=b(i+2)
	      t4=b(i+3)
	      tr1=c1*t3-s1*t4
	      ti1=c1*t4+s1*t3
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j2  )=t1-tr1
	      a(j2+1)=t2-ti1
	    enddo
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  RETURN
	endif
!...Passes 4 and 5...
	if(n2==32)then
	  ilo=1; ihi=incr
	  j=1
	  do m=1,8
	    c1=w3(1,j  )
	    s1=w3(2,j  )
	    c2=w3(1,j+1)
	    s2=w3(2,j+1)
	    do i=ilo,ihi,8
	      j1=ishft(i,-2)+1
	      j2=j1+n4
	      j3=j2+n4
	      j4=j3+n4
	      t1=b(i  )
	      t2=b(i+1)
	      t3=b(i+2)
	      t4=b(i+3)
	      t5=b(i+4)
	      t6=b(i+5)
	      t7=b(i+6)
	      t8=b(i+7)
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      tr2=c2*t7-s2*t8
	      ti2=c2*t8+s2*t7
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j2  )=t3-ti2
	      a(j2+1)=t4+tr2
	      a(j3  )=t1-tr1
	      a(j3+1)=t2-ti1
	      a(j4  )=t3+ti2
	      a(j4+1)=t4-tr2
	    enddo
	    j=j+2
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  RETURN
	endif
!...Passes 4-6...
	  j=4	!* 3*(1)+1
	  ilo=1; ihi=incr
	  do m=1,8
	    c1=w2(1,j  )
	    s1=w2(2,j  )
	    c2=w2(1,j+1)
	    s2=w2(2,j+1)
	    c3=w2(1,j+2)
	    s3=w2(2,j+2)
	    tr2=c3*isqrt2
	    ti2=s3*isqrt2
	    do i=ilo,ihi,16
!	put indices are here
	      j1=ishft(i,-3)+1
	      j2=j1+n8
	      j3=j2+n8
	      j4=j3+n8
	      j5=j4+n8
	      j6=j5+n8
	      j7=j6+n8
	      j8=j7+n8
!	gather the needed data (8 64-byte complex, i.e. 16 64-byte reals)...
	      t1 =b(i  )
	      t2 =b(i+1)
	      t3 =b(i+2)
	      t4 =b(i+3)
	      t5 =b(i+4)
	      t6 =b(i+5)
	      t7 =b(i+6)
	      t8 =b(i+7)
	      t9 =b(i+8)
	      t10=b(i+9)
	      t11=b(i+10)
	      t12=b(i+11)
	      t13=b(i+12)
	      t14=b(i+13)
	      t15=b(i+14)
	      t16=b(i+15)
!	first get the 4 length-2 transforms...
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tt =c1*t11-s1*t12
	      t12=c1*t12+s1*t11
	      t11=t9 -tt
	      t9 =t9 +tt
	      tt =t12
	      t12=t10-tt
	      t10=t10+tt
	      tt =c1*t15-s1*t16
	      t16=c1*t16+s1*t15
	      t15=t13-tt
	      t13=t13+tt
	      tt =t16
	      t16=t14-tt
	      t14=t14+tt
!	combine to get the 2 length-4 transform...
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      t5=t1-tr1
	      t1=t1+tr1
	      t6=t2-ti1
	      t2=t2+ti1
	      tr1=c2*t7-s2*t8
	      ti1=c2*t8+s2*t7
	      t7=t3+ti1
	      t3=t3-ti1
	      t8=t4-tr1
	      t4=t4+tr1
	      tr1=c2*t13-s2*t14
	      ti1=c2*t14+s2*t13
	      t13=t9 -tr1
	      t9 =t9 +tr1
	      t14=t10-ti1
	      t10=t10+ti1
	      tr1=c2*t15-s2*t16
	      ti1=c2*t16+s2*t15
	      t15=t11+ti1
	      t11=t11-ti1
	      t16=t12-tr1
	      t12=t12+tr1
!	now combine the two half-transforms
	      tr1=c3*t9 -s3*t10
	      ti1=c3*t10+s3*t9
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j5  )=t1-tr1
	      a(j5+1)=t2-ti1

	      tt =tr2*t11-ti2*t12
	      ti1=tr2*t12+ti2*t11
	      tr1=tt-ti1		!* mpy by (1+i)/sqrt(2) is here...
	      ti1=tt+ti1
	      a(j2  )=t3+tr1
	      a(j2+1)=t4+ti1
	      a(j6  )=t3-tr1
	      a(j6+1)=t4-ti1

	      tr1=c3*t13-s3*t14
	      ti1=c3*t14+s3*t13
	      a(j3  )=t5-ti1		!* mpy by i is inlined here...
	      a(j3+1)=t6+tr1
	      a(j7  )=t5+ti1
	      a(j7+1)=t6-tr1

	      tt =tr2*t15-ti2*t16
	      ti1=tr2*t16+ti2*t15
	      tr1=tt+ti1		!* mpy by (1-i)/sqrt(2) is here...
	      ti1=ti1-tt
	      a(j4  )=t7-tr1		!* and get (i-1)/sqrt by flipping signs here.
	      a(j4+1)=t8-ti1
	      a(j8  )=t7+tr1
	      a(j8+1)=t8+ti1
	    enddo
	    j=j+3
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	if(n2==64)RETURN
!...Pass 7...
	incr=incr/8
	if(n2==128)then
	  ilo=1; ihi=incr
	  do m=1,64
	    c1=w3(1,m)
	    s1=w3(2,m)
	    do i=ilo,ihi,4
	      j1=ishft(i,-1)+1
	      j2=j1+n2
	      t1=a(i  )
	      t2=a(i+1)
	      t3=a(i+2)
	      t4=a(i+3)
	      tr1=c1*t3-s1*t4
	      ti1=c1*t4+s1*t3
	      b(j1  )=t1+tr1
	      b(j1+1)=t2+ti1
	      b(j2  )=t1-tr1
	      b(j2+1)=t2-ti1
	    enddo
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  a=b; RETURN
	endif
!...Passes 7 and 8...
	if(n2==256)then
	  ilo=1; ihi=incr
	  j=1
	  do m=1,64
	    c1=w3(1,j  )
	    s1=w3(2,j  )
	    c2=w3(1,j+1)
	    s2=w3(2,j+1)
	    do i=ilo,ihi,8
	      j1=ishft(i,-2)+1
	      j2=j1+n4
	      j3=j2+n4
	      j4=j3+n4
	      t1=a(i  )
	      t2=a(i+1)
	      t3=a(i+2)
	      t4=a(i+3)
	      t5=a(i+4)
	      t6=a(i+5)
	      t7=a(i+6)
	      t8=a(i+7)
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      tr2=c2*t7-s2*t8
	      ti2=c2*t8+s2*t7
	      b(j1  )=t1+tr1
	      b(j1+1)=t2+ti1
	      b(j2  )=t3-ti2
	      b(j2+1)=t4+tr2
	      b(j3  )=t1-tr1
	      b(j3+1)=t2-ti1
	      b(j4  )=t3+ti2
	      b(j4+1)=t4-tr2
	    enddo
	    j=j+2
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  a=b; RETURN
	endif
!...Passes 7-9...
	  j=28	!* 3*(1+8)+1
	  ilo=1; ihi=incr
	  do m=1,64
	    c1=w2(1,j  )
	    s1=w2(2,j  )
	    c2=w2(1,j+1)
	    s2=w2(2,j+1)
	    c3=w2(1,j+2)
	    s3=w2(2,j+2)
	    tr2=c3*isqrt2
	    ti2=s3*isqrt2
	    do i=ilo,ihi,16
!	put indices are here
	      j1=ishft(i,-3)+1
	      j2=j1+n8
	      j3=j2+n8
	      j4=j3+n8
	      j5=j4+n8
	      j6=j5+n8
	      j7=j6+n8
	      j8=j7+n8
!	gather the needed data (8 64-byte complex, i.e. 16 64-byte reals)...
	      t1 =a(i  )
	      t2 =a(i+1)
	      t3 =a(i+2)
	      t4 =a(i+3)
	      t5 =a(i+4)
	      t6 =a(i+5)
	      t7 =a(i+6)
	      t8 =a(i+7)
	      t9 =a(i+8)
	      t10=a(i+9)
	      t11=a(i+10)
	      t12=a(i+11)
	      t13=a(i+12)
	      t14=a(i+13)
	      t15=a(i+14)
	      t16=a(i+15)
!	first get the 4 length-2 transforms...
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tt =c1*t11-s1*t12
	      t12=c1*t12+s1*t11
	      t11=t9 -tt
	      t9 =t9 +tt
	      tt =t12
	      t12=t10-tt
	      t10=t10+tt
	      tt =c1*t15-s1*t16
	      t16=c1*t16+s1*t15
	      t15=t13-tt
	      t13=t13+tt
	      tt =t16
	      t16=t14-tt
	      t14=t14+tt
!	combine to get the 2 length-4 transform...
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      t5=t1-tr1
	      t1=t1+tr1
	      t6=t2-ti1
	      t2=t2+ti1
	      tr1=c2*t7-s2*t8
	      ti1=c2*t8+s2*t7
	      t7=t3+ti1
	      t3=t3-ti1
	      t8=t4-tr1
	      t4=t4+tr1
	      tr1=c2*t13-s2*t14
	      ti1=c2*t14+s2*t13
	      t13=t9 -tr1
	      t9 =t9 +tr1
	      t14=t10-ti1
	      t10=t10+ti1
	      tr1=c2*t15-s2*t16
	      ti1=c2*t16+s2*t15
	      t15=t11+ti1
	      t11=t11-ti1
	      t16=t12-tr1
	      t12=t12+tr1
!	now combine the two half-transforms
	      tr1=c3*t9 -s3*t10
	      ti1=c3*t10+s3*t9
	      b(j1  )=t1+tr1
	      b(j1+1)=t2+ti1
	      b(j5  )=t1-tr1
	      b(j5+1)=t2-ti1

	      tt =tr2*t11-ti2*t12
	      ti1=tr2*t12+ti2*t11
	      tr1=tt-ti1		!* mpy by (1+i)/sqrt(2) is here...
	      ti1=tt+ti1
	      b(j2  )=t3+tr1
	      b(j2+1)=t4+ti1
	      b(j6  )=t3-tr1
	      b(j6+1)=t4-ti1

	      tr1=c3*t13-s3*t14
	      ti1=c3*t14+s3*t13
	      b(j3  )=t5-ti1		!* mpy by i is inlined here...
	      b(j3+1)=t6+tr1
	      b(j7  )=t5+ti1
	      b(j7+1)=t6-tr1

	      tt =tr2*t15-ti2*t16
	      ti1=tr2*t16+ti2*t15
	      tr1=tt+ti1		!* mpy by (1-i)/sqrt(2) is here...
	      ti1=ti1-tt
	      b(j4  )=t7-tr1		!* and get (i-1)/sqrt by flipping signs here.
	      b(j4+1)=t8-ti1
	      b(j8  )=t7+tr1
	      b(j8+1)=t8+ti1
	    enddo
	    j=j+3
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	if(n2==512)then; a=b; RETURN; endif
!...Pass 10...
	incr=incr/8
	if(n2==1024)then
	  ilo=1; ihi=incr
	  do m=1,512
	    c1=w3(1,m)
	    s1=w3(2,m)
	    do i=ilo,ihi,4
	      j1=ishft(i,-1)+1
	      j2=j1+n2
	      t1=b(i  )
	      t2=b(i+1)
	      t3=b(i+2)
	      t4=b(i+3)
	      tr1=c1*t3-s1*t4
	      ti1=c1*t4+s1*t3
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j2  )=t1-tr1
	      a(j2+1)=t2-ti1
	    enddo
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  RETURN
	endif
!...Passes 10 and 11...
	if(n2==2048)then
	  ilo=1; ihi=incr
	  j=1
	  do m=1,512
	    c1=w3(1,j  )
	    s1=w3(2,j  )
	    c2=w3(1,j+1)
	    s2=w3(2,j+1)
	    do i=ilo,ihi,8
	      j1=ishft(i,-2)+1
	      j2=j1+n4
	      j3=j2+n4
	      j4=j3+n4
	      t1=b(i  )
	      t2=b(i+1)
	      t3=b(i+2)
	      t4=b(i+3)
	      t5=b(i+4)
	      t6=b(i+5)
	      t7=b(i+6)
	      t8=b(i+7)
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      tr2=c2*t7-s2*t8
	      ti2=c2*t8+s2*t7
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j2  )=t3-ti2
	      a(j2+1)=t4+tr2
	      a(j3  )=t1-tr1
	      a(j3+1)=t2-ti1
	      a(j4  )=t3+ti2
	      a(j4+1)=t4-tr2
	    enddo
	    j=j+2
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  RETURN
	endif
!...Passes 10-12...
	  j=220	!* 3*(1+8+64)+1
	  ilo=1; ihi=incr
	  do m=1,512
	    c1=w2(1,j  )
	    s1=w2(2,j  )
	    c2=w2(1,j+1)
	    s2=w2(2,j+1)
	    c3=w2(1,j+2)
	    s3=w2(2,j+2)
	    tr2=c3*isqrt2
	    ti2=s3*isqrt2
	    do i=ilo,ihi,16
!	put indices are here
	      j1=ishft(i,-3)+1
	      j2=j1+n8
	      j3=j2+n8
	      j4=j3+n8
	      j5=j4+n8
	      j6=j5+n8
	      j7=j6+n8
	      j8=j7+n8
!	gather the needed data (8 64-byte complex, i.e. 16 64-byte reals)...
	      t1 =b(i  )
	      t2 =b(i+1)
	      t3 =b(i+2)
	      t4 =b(i+3)
	      t5 =b(i+4)
	      t6 =b(i+5)
	      t7 =b(i+6)
	      t8 =b(i+7)
	      t9 =b(i+8)
	      t10=b(i+9)
	      t11=b(i+10)
	      t12=b(i+11)
	      t13=b(i+12)
	      t14=b(i+13)
	      t15=b(i+14)
	      t16=b(i+15)
!	first get the 4 length-2 transforms...
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tt =c1*t11-s1*t12
	      t12=c1*t12+s1*t11
	      t11=t9 -tt
	      t9 =t9 +tt
	      tt =t12
	      t12=t10-tt
	      t10=t10+tt
	      tt =c1*t15-s1*t16
	      t16=c1*t16+s1*t15
	      t15=t13-tt
	      t13=t13+tt
	      tt =t16
	      t16=t14-tt
	      t14=t14+tt
!	combine to get the 2 length-4 transform...
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      t5=t1-tr1
	      t1=t1+tr1
	      t6=t2-ti1
	      t2=t2+ti1
	      tr1=c2*t7-s2*t8
	      ti1=c2*t8+s2*t7
	      t7=t3+ti1
	      t3=t3-ti1
	      t8=t4-tr1
	      t4=t4+tr1
	      tr1=c2*t13-s2*t14
	      ti1=c2*t14+s2*t13
	      t13=t9 -tr1
	      t9 =t9 +tr1
	      t14=t10-ti1
	      t10=t10+ti1
	      tr1=c2*t15-s2*t16
	      ti1=c2*t16+s2*t15
	      t15=t11+ti1
	      t11=t11-ti1
	      t16=t12-tr1
	      t12=t12+tr1
!	now combine the two half-transforms
	      tr1=c3*t9 -s3*t10
	      ti1=c3*t10+s3*t9
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j5  )=t1-tr1
	      a(j5+1)=t2-ti1

	      tt =tr2*t11-ti2*t12
	      ti1=tr2*t12+ti2*t11
	      tr1=tt-ti1		!* mpy by (1+i)/sqrt(2) is here...
	      ti1=tt+ti1
	      a(j2  )=t3+tr1
	      a(j2+1)=t4+ti1
	      a(j6  )=t3-tr1
	      a(j6+1)=t4-ti1

	      tr1=c3*t13-s3*t14
	      ti1=c3*t14+s3*t13
	      a(j3  )=t5-ti1		!* mpy by i is inlined here...
	      a(j3+1)=t6+tr1
	      a(j7  )=t5+ti1
	      a(j7+1)=t6-tr1

	      tt =tr2*t15-ti2*t16
	      ti1=tr2*t16+ti2*t15
	      tr1=tt+ti1		!* mpy by (1-i)/sqrt(2) is here...
	      ti1=ti1-tt
	      a(j4  )=t7-tr1		!* and get (i-1)/sqrt by flipping signs here.
	      a(j4+1)=t8-ti1
	      a(j8  )=t7+tr1
	      a(j8+1)=t8+ti1
	    enddo
	    j=j+3
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	if(n2==4096)RETURN
!...Pass 13...
	incr=incr/8
	if(n2==8192)then
	  ilo=1; ihi=incr
	  do m=1,4096
	    c1=w3(1,m)
	    s1=w3(2,m)
	    do i=ilo,ihi,4
	      j1=ishft(i,-1)+1
	      j2=j1+n2
	      t1=a(i  )
	      t2=a(i+1)
	      t3=a(i+2)
	      t4=a(i+3)
	      tr1=c1*t3-s1*t4
	      ti1=c1*t4+s1*t3
	      b(j1  )=t1+tr1
	      b(j1+1)=t2+ti1
	      b(j2  )=t1-tr1
	      b(j2+1)=t2-ti1
	    enddo
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  a=b; RETURN
	endif
!...Passes 13 and 14...
	if(n2==16384)then
	  ilo=1; ihi=incr
	  j=1
	  do m=1,4096
	    c1=w3(1,j  )
	    s1=w3(2,j  )
	    c2=w3(1,j+1)
	    s2=w3(2,j+1)
	    do i=ilo,ihi,8
	      j1=ishft(i,-2)+1
	      j2=j1+n4
	      j3=j2+n4
	      j4=j3+n4
	      t1=a(i  )
	      t2=a(i+1)
	      t3=a(i+2)
	      t4=a(i+3)
	      t5=a(i+4)
	      t6=a(i+5)
	      t7=a(i+6)
	      t8=a(i+7)
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      tr2=c2*t7-s2*t8
	      ti2=c2*t8+s2*t7
	      b(j1  )=t1+tr1
	      b(j1+1)=t2+ti1
	      b(j2  )=t3-ti2
	      b(j2+1)=t4+tr2
	      b(j3  )=t1-tr1
	      b(j3+1)=t2-ti1
	      b(j4  )=t3+ti2
	      b(j4+1)=t4-tr2
	    enddo
	    j=j+2
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  a=b; RETURN
	endif
!...Passes 13-15...
	  j=1756	!* 3*(1+8+64+512)+1
	  ilo=1; ihi=incr
	  do m=1,4096
	    c1=w2(1,j  )
	    s1=w2(2,j  )
	    c2=w2(1,j+1)
	    s2=w2(2,j+1)
	    c3=w2(1,j+2)
	    s3=w2(2,j+2)
	    tr2=c3*isqrt2
	    ti2=s3*isqrt2
	    do i=ilo,ihi,16
!	put indices are here
	      j1=ishft(i,-3)+1
	      j2=j1+n8
	      j3=j2+n8
	      j4=j3+n8
	      j5=j4+n8
	      j6=j5+n8
	      j7=j6+n8
	      j8=j7+n8
!	gather the needed data (8 64-byte complex, i.e. 16 64-byte reals)...
	      t1 =a(i  )
	      t2 =a(i+1)
	      t3 =a(i+2)
	      t4 =a(i+3)
	      t5 =a(i+4)
	      t6 =a(i+5)
	      t7 =a(i+6)
	      t8 =a(i+7)
	      t9 =a(i+8)
	      t10=a(i+9)
	      t11=a(i+10)
	      t12=a(i+11)
	      t13=a(i+12)
	      t14=a(i+13)
	      t15=a(i+14)
	      t16=a(i+15)
!	first get the 4 length-2 transforms...
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tt =c1*t11-s1*t12
	      t12=c1*t12+s1*t11
	      t11=t9 -tt
	      t9 =t9 +tt
	      tt =t12
	      t12=t10-tt
	      t10=t10+tt
	      tt =c1*t15-s1*t16
	      t16=c1*t16+s1*t15
	      t15=t13-tt
	      t13=t13+tt
	      tt =t16
	      t16=t14-tt
	      t14=t14+tt
!	combine to get the 2 length-4 transform...
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      t5=t1-tr1
	      t1=t1+tr1
	      t6=t2-ti1
	      t2=t2+ti1
	      tr1=c2*t7-s2*t8
	      ti1=c2*t8+s2*t7
	      t7=t3+ti1
	      t3=t3-ti1
	      t8=t4-tr1
	      t4=t4+tr1
	      tr1=c2*t13-s2*t14
	      ti1=c2*t14+s2*t13
	      t13=t9 -tr1
	      t9 =t9 +tr1
	      t14=t10-ti1
	      t10=t10+ti1
	      tr1=c2*t15-s2*t16
	      ti1=c2*t16+s2*t15
	      t15=t11+ti1
	      t11=t11-ti1
	      t16=t12-tr1
	      t12=t12+tr1
!	now combine the two half-transforms
	      tr1=c3*t9 -s3*t10
	      ti1=c3*t10+s3*t9
	      b(j1  )=t1+tr1
	      b(j1+1)=t2+ti1
	      b(j5  )=t1-tr1
	      b(j5+1)=t2-ti1

	      tt =tr2*t11-ti2*t12
	      ti1=tr2*t12+ti2*t11
	      tr1=tt-ti1		!* mpy by (1+i)/sqrt(2) is here...
	      ti1=tt+ti1
	      b(j2  )=t3+tr1
	      b(j2+1)=t4+ti1
	      b(j6  )=t3-tr1
	      b(j6+1)=t4-ti1

	      tr1=c3*t13-s3*t14
	      ti1=c3*t14+s3*t13
	      b(j3  )=t5-ti1		!* mpy by i is inlined here...
	      b(j3+1)=t6+tr1
	      b(j7  )=t5+ti1
	      b(j7+1)=t6-tr1

	      tt =tr2*t15-ti2*t16
	      ti1=tr2*t16+ti2*t15
	      tr1=tt+ti1		!* mpy by (1-i)/sqrt(2) is here...
	      ti1=ti1-tt
	      b(j4  )=t7-tr1		!* and get (i-1)/sqrt by flipping signs here.
	      b(j4+1)=t8-ti1
	      b(j8  )=t7+tr1
	      b(j8+1)=t8+ti1
	    enddo
	    j=j+3
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	if(n2==32768)then; a=b; RETURN; endif
!...Pass 16...
	incr=incr/8
	if(n2==65536)then
	  ilo=1; ihi=incr
	  do m=1,32768
	    c1=w3(1,m)
	    s1=w3(2,m)
	    do i=ilo,ihi,4
	      j1=ishft(i,-1)+1
	      j2=j1+n2
	      t1=b(i  )
	      t2=b(i+1)
	      t3=b(i+2)
	      t4=b(i+3)
	      tr1=c1*t3-s1*t4
	      ti1=c1*t4+s1*t3
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j2  )=t1-tr1
	      a(j2+1)=t2-ti1
	    enddo
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  RETURN
	endif
!...Passes 16 and 17...
	if(n2==131072)then
	  ilo=1; ihi=incr
	  j=1
	  do m=1,32768
	    c1=w3(1,j  )
	    s1=w3(2,j  )
	    c2=w3(1,j+1)
	    s2=w3(2,j+1)
	    do i=ilo,ihi,8
	      j1=ishft(i,-2)+1
	      j2=j1+n4
	      j3=j2+n4
	      j4=j3+n4
	      t1=b(i  )
	      t2=b(i+1)
	      t3=b(i+2)
	      t4=b(i+3)
	      t5=b(i+4)
	      t6=b(i+5)
	      t7=b(i+6)
	      t8=b(i+7)
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      tr2=c2*t7-s2*t8
	      ti2=c2*t8+s2*t7
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j2  )=t3-ti2
	      a(j2+1)=t4+tr2
	      a(j3  )=t1-tr1
	      a(j3+1)=t2-ti1
	      a(j4  )=t3+ti2
	      a(j4+1)=t4-tr2
	    enddo
	    j=j+2
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	  RETURN
	endif
!...Passes 16-18...
	  j=14044	!* 3*(1+8+64+512+4096)+1
	  ilo=1; ihi=incr
	  do m=1,32768
	    c1=w2(1,j  )
	    s1=w2(2,j  )
	    c2=w2(1,j+1)
	    s2=w2(2,j+1)
	    c3=w2(1,j+2)
	    s3=w2(2,j+2)
	    tr2=c3*isqrt2
	    ti2=s3*isqrt2
	    do i=ilo,ihi,16
!	put indices are here
	      j1=ishft(i,-3)+1
	      j2=j1+n8
	      j3=j2+n8
	      j4=j3+n8
	      j5=j4+n8
	      j6=j5+n8
	      j7=j6+n8
	      j8=j7+n8
!	gather the needed data (8 64-byte complex, i.e. 16 64-byte reals)...
	      t1 =b(i  )
	      t2 =b(i+1)
	      t3 =b(i+2)
	      t4 =b(i+3)
	      t5 =b(i+4)
	      t6 =b(i+5)
	      t7 =b(i+6)
	      t8 =b(i+7)
	      t9 =b(i+8)
	      t10=b(i+9)
	      t11=b(i+10)
	      t12=b(i+11)
	      t13=b(i+12)
	      t14=b(i+13)
	      t15=b(i+14)
	      t16=b(i+15)
!	first get the 4 length-2 transforms...
	      tt=c1*t3-s1*t4
	      t4=c1*t4+s1*t3
	      t3=t1-tt
	      t1=t1+tt
	      tt=t4
	      t4=t2-tt
	      t2=t2+tt
	      tt=c1*t7-s1*t8
	      t8=c1*t8+s1*t7
	      t7=t5-tt
	      t5=t5+tt
	      tt=t8
	      t8=t6-tt
	      t6=t6+tt
	      tt =c1*t11-s1*t12
	      t12=c1*t12+s1*t11
	      t11=t9 -tt
	      t9 =t9 +tt
	      tt =t12
	      t12=t10-tt
	      t10=t10+tt
	      tt =c1*t15-s1*t16
	      t16=c1*t16+s1*t15
	      t15=t13-tt
	      t13=t13+tt
	      tt =t16
	      t16=t14-tt
	      t14=t14+tt
!	combine to get the 2 length-4 transform...
	      tr1=c2*t5-s2*t6
	      ti1=c2*t6+s2*t5
	      t5=t1-tr1
	      t1=t1+tr1
	      t6=t2-ti1
	      t2=t2+ti1
	      tr1=c2*t7-s2*t8
	      ti1=c2*t8+s2*t7
	      t7=t3+ti1
	      t3=t3-ti1
	      t8=t4-tr1
	      t4=t4+tr1
	      tr1=c2*t13-s2*t14
	      ti1=c2*t14+s2*t13
	      t13=t9 -tr1
	      t9 =t9 +tr1
	      t14=t10-ti1
	      t10=t10+ti1
	      tr1=c2*t15-s2*t16
	      ti1=c2*t16+s2*t15
	      t15=t11+ti1
	      t11=t11-ti1
	      t16=t12-tr1
	      t12=t12+tr1
!	now combine the two half-transforms
	      tr1=c3*t9 -s3*t10
	      ti1=c3*t10+s3*t9
	      a(j1  )=t1+tr1
	      a(j1+1)=t2+ti1
	      a(j5  )=t1-tr1
	      a(j5+1)=t2-ti1

	      tt =tr2*t11-ti2*t12
	      ti1=tr2*t12+ti2*t11
	      tr1=tt-ti1		!* mpy by (1+i)/sqrt(2) is here...
	      ti1=tt+ti1
	      a(j2  )=t3+tr1
	      a(j2+1)=t4+ti1
	      a(j6  )=t3-tr1
	      a(j6+1)=t4-ti1

	      tr1=c3*t13-s3*t14
	      ti1=c3*t14+s3*t13
	      a(j3  )=t5-ti1		!* mpy by i is inlined here...
	      a(j3+1)=t6+tr1
	      a(j7  )=t5+ti1
	      a(j7+1)=t6-tr1

	      tt =tr2*t15-ti2*t16
	      ti1=tr2*t16+ti2*t15
	      tr1=tt+ti1		!* mpy by (1-i)/sqrt(2) is here...
	      ti1=ti1-tt
	      a(j4  )=t7-tr1		!* and get (i-1)/sqrt by flipping signs here.
	      a(j4+1)=t8-ti1
	      a(j8  )=t7+tr1
	      a(j8+1)=t8+ti1
	    enddo
	    j=j+3
	    ilo=ilo+incr; ihi=ihi+incr
	  enddo
	if(n2==262144)RETURN
!...
	print*,'FFT length not a power of 2 or out of range'
	STOP
!...
	end subroutine dlanczos_fwd

