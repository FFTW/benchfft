      subroutine stimer
      real tw, tc
      double precision tf,fputim
      real tarray(2)
      common /timer/ tw,tc,tf
      tw = secnds(0.0)
      tc = etime(tarray)
      return
      end

      subroutine ptimer(nl,cloops)
      real tw, tc
      double precision tf,nl,cloops,fputim
      real tarray(2)
      common /timer/ tw,tc,tf
      tw = secnds(tw)
      tc = etime(tarray) -tc
      tf = tf*1000./(2*cloops)
       tw = tw*1000./(2*cloops)
 220   tc = tc*1000./(2*cloops)
            print '(A13,'' n = '',F7.0,''    wall time = '',G12.2,
     *           '' cpu time = '', G12.2,
     *           '' fpu time = '', G12.2, '' ('',F7.0,'')'')',
     *           'GPFA' , nl, tw, tc, tf, cloops
      return
      end

      function drand48()
      drand48=drand(0)
      return
      end
