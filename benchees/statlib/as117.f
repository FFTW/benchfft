        subroutine chrft(x, y, wr, wi, isize, lwork, itype, ifault)
c
c       Algorithm AS 117.1 Appl. Statist. (1977) vol.26, no.3
c
c       Code for complex discrete Fourier transform of a sequence of
c       general length by the Chirp-Z transform
c
        real    x(lwork), y(lwork), wr(lwork), wi(lwork)
        data    zero, one/0.0, 1.0/
c
c       Check that ISIZE and LWORK are valid
c
        ifault = 0
        if (isize .lt. 3) ifault = 1
        ii = 4
        do 1 k = 3, 20
          ii = ii * 2
          if (ii .eq. lwork) go to 3
          if (ii .gt. lwork) go to 2
    1   continue

    2   ifault = ifault + 2
    3   if (lwork .lt. 2 * isize) ifault = ifault + 4
        if (itype .eq. 0) ifault = ifault + 8
        if (ifault .ne. 0) return
c
c       Multiply by the Chirp function in the time domain
c
        call chrp(x, y, one, isize, lwork, itype)
        np1 = isize + 1
        do 4 nn = np1, lwork
          x(nn) = zero
          y(nn) = zero
    4   continue
c
c       Fourier transform the chirped series
c
        call fastf(x, y, lwork, 1)
c
c       Convolve by frequency domain multiplication with (wr, wi)
c
        do 5 nn = 1, lwork
          xwi = wi(nn)
          if (itype .lt. 0) xwi = -xwi
          z = x(nn) * wr(nn) - y(nn) * xwi
          y(nn) = x(nn) * xwi + y(nn) * wr(nn)
          x(nn) = z
    5   continue
c
c       Inverse Fourier transform
c
        call fastf(x, y, lwork, -1)
c
c       Multiply by Chirp function & gain correction
c
        gc = lwork
        if (itype .gt. 0) gc = gc/isize
        call chrp(x, y, gc, isize, lwork, itype)

        return
        end



        subroutine chfor(xray, wr, wi, jsize, lwork, mwork, ifault)
c
c       Algorithm AS 117.2 Appl. Statist. (1977) vol.26, no.3
c
c       Forward Chirp DFT for real even length input sequence
c
        real    xray(lwork), wr(mwork), wi(mwork)
        data    zero, quart, half, one, one5, two, four
     +          /0.0,  0.25,  0.5, 1.0,  1.5, 2.0,  4.0/
c
c       Check for valid JSIZE, LWORK & MWORK
c
        ifault = 0
        ii = 8
        do 1 k = 4, 21
          ii = ii * 2
          if (ii .eq. lwork) go to 3
          if (ii .gt. lwork) go to 2
    1   continue
    2   ifault = 2
    3   if (lwork .lt. 2 * jsize) ifault = ifault + 4
        if (lwork .ne. 2 * mwork) ifault = ifault + 16
        nn = jsize / 2
        if (2 * nn .ne. jsize) ifault = ifault + 32
        if (nn .lt. 3) ifault = ifault + 1
        if (ifault .ne. 0) return
        ll = mwork
c
c       Split XRAY into even and odd sequences of length n/2 placing
c       the even terms in the bottom half of XRAY as dummy real terms
c       and odd terms in the top half as dummy imaginary terms.
c
        do 4 i = 1, nn
          ii = 2 * i
          m = ll + i
          xray(m) = xray(ii)
          xray(i) = xray(ii-1)
    4   continue
c
c       Perform forward Chirp DFT on even and odd sequences.
c       IFA is not tested as possible faults have already been checked.
c
        call chrft(xray, xray(ll+1), wr, wi, nn, ll, 1, ifa)
c
c       Reorder the results according to output format.
c       First, the unique real terms.
c
        z = half * (xray(1) + xray(ll+1))
        xray(nn+1) = half * (xray(1) - xray(ll+1))
        xray(1) = z
c
c       If nn is even, terms at nn/2 + 1 and 3*nn/2 can be calculated
c
        nnn = nn / 2
        if (nn .ne. 2 * nnn) go to 5
        m = nnn + 1
        xray(m) = half * xray(m)
        mm = m + ll
        m = m + nn
        xray(m) = -half * xray(mm)
        go to 6
    5   nnn = nnn + 1
c
c       Set up trig functions for calculations of remaining
c       non-unique terms
c
    6   z = four * atan(one) / nn
        bcos = -two * (sin(z/two) ** 2)
        bsin = sin(z)
        un = one
        vn = zero
c
c       Calculate & place remaining terms in correct locations.
c
        do 7 i = 2, nnn
          z = un * bcos + un + vn * bsin
          vn = vn * bcos + vn - un * bsin
          save1 = one5 - half * (z * z + vn * vn)
          un = z * save1
          vn = vn * save1
          ii = nn + 2 - i
          m = ll + i
          mm = ll + ii
          an = quart * (xray(i) + xray(ii))
          bn = quart * (xray(m) - xray(mm))
          cn = quart * (xray(m) + xray(mm))
          dn = quart * (xray(i) - xray(ii))
          xray(i) = an + un * cn + vn * dn
          xray(ii) = two * an - xray(i)
          j = nn + i
          jj = nn + ii
          xray(j) = bn - un * dn + vn * cn
          xray(jj) = xray(j) - two * bn
    7   continue

        return
        end



        subroutine chrev(xray, wr, wi, jsize, lwork, mwork, ifault)
c
c       Algorithm AS 117.3 Appl. Statist. (1977) vol.26, no.3
c
c       Inverse Chirp DFT to give real even length sequence
c
        real    xray(lwork), wr(mwork), wi(mwork)
        data    zero, half, one, one5, two, four
     +          /0.0,  0.5, 1.0,  1.5, 2.0,  4.0/
c
c       Check for valid JSIZE, LWORK & MWORK
c
        ifault = 0
        ii = 8
        do 1 k = 4, 21
          ii = ii * 2
          if (ii .eq. lwork) go to 3
          if (ii .gt. lwork) go to 2
    1   continue
    2   ifault = 2
    3   if (lwork .lt. 2 * jsize) ifault = ifault + 4
        if (lwork .ne. 2 * mwork) ifault = ifault + 16
        nn = jsize / 2
        if (2 * nn .ne. jsize) ifault = ifault + 32
        if (nn .lt. 3) ifault = ifault + 1
        if (ifault .ne. 0) return
        ll = mwork
c
c       Reorder the spectrum; first the unique terms.
c
        z = xray(1) + xray(nn+1)
        xray(ll+1) = xray(1) - xray(nn+1)
        xray(1) = z
c
c       If nn is even then terms nn/2 + 1 and 3*nn/2 + 1 must be
c       reordered.
c
        nnn = nn / 2
        if (nn .ne. 2 * nnn) go to 4
        m = nnn + 1
        xray(m) = two * xray(m)
        mm = m + ll
        m = m + nn
        xray(mm) = -two * xray(m)
        go to 5
    4   nnn = nnn + 1
c
c       Set up trig functions for manipulation of remaining terms.
c
    5   z = four * atan(one) / nn
        bcos = -two * (sin(z/two) ** 2)
        bsin = sin(z)
        un = one
        vn = zero
c
c       Perform manipulation and reordering of remaining non-unique
c       terms.
c
        do 6 i = 2, nnn
          z = un * bcos + un + vn * bsin
          vn = vn * bcos + vn - un * bsin
          save1 = one5 - half * (z * z + vn * vn)
          un = z * save1
          vn = vn * save1
          ii = nn + 2 - i
          j = nn + i
          jj = nn + ii
          an = xray(i) + xray(ii)
          bn = xray(j) - xray(jj)
          pn = xray(i) - xray(ii)
          qn = xray(j) + xray(jj)
          cn = un * pn + vn * qn
          dn = vn * pn - un * qn
          xray(i) = an + dn
          xray(ii) = an - dn
          m = ll + i
          mm = ll + ii
          xray(m) = cn + bn
          xray(mm) = cn - bn
    6   continue
c
c       Do inverse Chirp DFT to give even and odd sequences.
c       IFA is not tested as possible faults have already been checked
c
        call chrft(xray, xray(ll+1), wr, wi, nn, ll, -1, ifa)
c
c       Interlace the results to produce the required output sequence.
c
        nnn = nn + 1
        ii = jsize + 2
        do 7 i = 1, nn
          j = nnn - i
          jj = ll + j
          m = ii - 2 * i
          mm = m - 1
          xray(m) = xray(jj)
          xray(mm) = xray(j)
    7   continue

        return
        end



      subroutine setwt(wr, wi, ksize, kwork, ifault)
c
c       Algorithm AS 117.4 Appl. Statist. (1977) vol.26, no.3
c
c       Subroutine to set up Fourier transformed Chirp function for
c       use by subroutine CHRFT.
c
      real    wr(kwork), wi(kwork)
      data    zero, one, four /0.0, 1.0, 4.0/
c
c       Check that KSIZE & KWORK are valid
c
      ifault = 0
      if (ksize .lt. 3) ifault = 1
      ii = 4
      do 1 k = 3, 20
        ii = ii * 2
        if (ii .eq. kwork) go to 3
        if (ii .gt. kwork) go to 2
    1   continue
    2   ifault = 2
    3   if (kwork .lt. 2 * ksize) ifault = ifault + 4
      if (ifault .ne. 0) return
      tc = four * atan(one) / ksize
c
c       Set up bottom segment of Chirp function
c
      do 4 nn = 1, ksize
        z = nn - 1
        z = z * z * tc
        wr(nn) = cos(z)
        wi(nn) = sin(z)
    4 continue
c
c       Clear the rest
c
      do 5 nn = ksize+1, kwork
        wr(nn) = zero
        wi(nn) = zero
    5 continue
c
c       Copy to the top segment
c
      do 6 nn = kwork-ksize+2, kwork
        ll = kwork - nn + 2
        wr(nn) = wr(ll)
        wi(nn) = wi(ll)
    6 continue
c
c       Fourier transform the Chirp function
c
      call fastf(wr, wi, kwork, 1)

      return
      end


      subroutine chrp(x, y, gc, isize, lwork, itype)
c
c       Algorithm AS 117.5 Appl. Statist. (1977) vol.26, no.3
c
c       Subroutine to multiply time series by Chirp function
c
      dimension x(lwork), y(lwork)
      data    one, four /1.0, 4.0/

      tc = four * atan(one) / isize
      do 3 nn = 1, isize
        z = nn - 1
        z = z * z * tc
        xwr = cos(z)
        xwi = -sin(z)
        if (itype .lt. 0) xwi = -xwi
        z = x(nn) * xwr - y(nn) * xwi
        y(nn) = (x(nn) * xwi + y(nn) * xwr) * gc
        x(nn) = z * gc
    3 continue

      return
      end
