      FUNCTION dfridr(func,x,h,err)
      INTEGER NTAB
      REAL dfridr,err,h,x,func,CON,CON2,BIG,SAFE
      PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.E30,NTAB=10,SAFE=2.)
      EXTERNAL func
CU    USES func
      INTEGER i,j
      REAL errt,fac,hh,a(NTAB,NTAB)
      if(h.eq.0.) pause 'h must be nonzero in dfridr'
      hh=h
      a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)
      err=BIG
      do 12 i=2,NTAB
        hh=hh/CON
        a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
        fac=CON2
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dfridr=a(j,i)
          endif
11      continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.
