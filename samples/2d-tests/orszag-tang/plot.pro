; To be compared with
; http://www.astro.princeton.edu/~jstone/tests/orszag-tang/pagesource.html
rho = 25.0*reform(lnrho[l1:l2,m1:m2,3])/(36.*!pi)
tvlct,red,green,blue,/get
loadct,15
tv,rebin(bytscl(rho,min=0.,max=.75),2*[nx,ny])
tvlct,red,green,blue
end
