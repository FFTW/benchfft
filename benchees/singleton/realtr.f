      subroutine realtr(a,b,n,isn)
c  if isn=1, this subroutine completes the fourier transform
c    of 2*n real data values, where the original data values are
c    stored alternately in arrays a and b, and are first
c    transformed by a complex fourier transform of dimension n.
c    the cosine coefficients are in a(1),a(2),...,a(n+1) and
c    the sine coefficients are in b(1),b(2),...,b(n+1).
c    a typical calling sequence is
c      call fft(a,b,n,n,n,1)
c      call realtr(a,b,n,1)
c    the results should be multiplied by 1.0/(2.0*n) to give the
c    usual scaling of coefficients.
c  if isn=-1, the inverse transformation is done, the first step
c    in evaluating a real fourier series.
c    a typical calling sequence is
c      call realtr(a,b,n,-1)
c      call fft(a,b,n,n,n,-1)
c    the results should be multiplied by 0.5 to give the usual
c    scaling, and the time domain results alternate in arrays a
c    and b, i.e. a(1),b(1),a(2),b(2),...,a(n),b(n).
c  with most fortran compilers the data can alternatively be
c    stored in a single complex array a, then the magnitude of isn
c    changed to two to give the correct indexing increment and a(2)
c    used to pass the initial address for the sequence of imaginary
c    values, e.g.
c      call fft(a,a(2),n,n,n,2)
c      call realtr(a,a(2),n,2)
c    in this case, the cosine and sine coefficients alternate in a.
c  by r. c. singleton, stanford research institute, sept. 1968
      dimension a(1),b(1)
      real im
      inc=iabs(isn)
      nk=n*inc+2
      nh=nk/2
      sd=2.0*atan(1.0)/float(n)
      cd=2.0*sin(sd)**2
      sd=sin(sd+sd)
      sn=0.0
      if(isn .lt. 0) go to 30
      cn=1.0
      a(nk-1)=a(1)
      b(nk-1)=b(1)
   10 do 20 j=1,nh,inc
      k=nk-j
      aa=a(j)+a(k)
      ab=a(j)-a(k)
      ba=b(j)+b(k)
      bb=b(j)-b(k)
      re=cn*ba+sn*ab
      im=sn*ba-cn*ab
      b(k)=im-bb
      b(j)=im+bb
      a(k)=aa-re
      a(j)=aa+re
      aa=cn-(cd*cn+sd*sn)
      sn=(sd*cn-cd*sn)+sn
      cn=2.0-(aa**2+sn**2)
      sn=cn*sn
   20 cn=cn*aa
      return
   30 cn=-1.0
      sd=-sd
      go to 10
      end
