
function pc_meridional,field_in,rad,tht,$
                  plot=plot,amin=amin,amax=amax,ncolors=ncolors;,datadir=datadir
;
default,amin,0.
default,amax,1.
default,ncolors,256
;
pc_read_param, obj=param, dim=dim, datadir=datadir
;
field=field_in
;
nrad=n_elements(rad)
ntht=n_elements(tht)
;
rext=rad(nrad-1)
rint=rad(0)
;
if (keyword_set(plot)) then begin
    !p.multi=[0,2,1]
    contour,field,rad#sin(tht),rad#cos(tht),/fill,$
      lev=grange(amin,amax,ncolors),xs=3,ys=3,$
      title='Spherical plot',xtitle='!8r!x',ytitle='!8z!x'
endif
;
; limits of the cartesian grid
;
xxlim=minmax(rad#sin(tht))
x0=xxlim[0] & xn=xxlim[1]
nxc=nrad & xc=grange(x0,xn,Nxc)
;
zzlim=minmax(rad#cos(tht))
z0=zzlim[0] & zn=zzlim[1] & lpoint=nrad/2
zratio=rad(lpoint)/rad(0) & nzc=round(ntht*zratio)
zc=grange(z0,zn,Nzc)
;
; make the grid. use the resolution at lpoint in the z direction
;
fieldxz=fltarr(nxc,nzc)
fieldxz=fieldxz*0.
;
drad0=rad(1)-rad(0)
dtht=tht(1)-tht(0) & dtht1=1./dtht
;
rint=rad[0]
rext=rad[nrad-1]
;
for ix=0,nxc-1 do begin
  for iz=0,nzc-1 do begin
    radius = sqrt(xc[ix]^2+zc[iz]^2)
    theta=atan(xc[ix],zc[iz])
;
    if ((radius ge rint)  and (radius le rext)  $
                          and                   $
        (theta ge tht(0)) and (theta le tht(ntht-1))) then begin
;
       if (param.lequidist[0] eq -1) then begin
       ; usual indexing for equidistant grids
          ir1=floor((radius-rint)/drad0)
       endif else begin
       ; bisect to find the index in non-equidistant grids
       ; (imported from the fortran subroutine in the main code)
          jl = 0
          ju = n_elements(rad)-1
          while ((ju - jl) gt 1) do begin
             jm = (ju + jl) / 2
             if (radius gt rad[jm-1]) then begin
                jl = jm
             endif else begin
                ju = jm
             endelse
          endwhile
          if (radius - rad[jl-1] le rad[ju-1] - radius) then begin
             ir1 = jl-1
          endif else begin
             ir1 = ju-1
          endelse
       endelse

       ir2=ir1+1
       drad=rad(ir2)-rad(ir1) & drad1=1./drad

      if (ir2 eq Nrad) then ir2=Nrad-1
      if (ir1 lt    0) then begin
        print,'ir1 lt 0. radius,ir1=',radius,ir1
        stop
      endif
      if (ir2 gt Nrad) then begin
        print,'ir2 ge Nr. radius,ir2=',radius,ir2
        stop
      endif
      delr = radius - rad[ir1]
;
      it1=floor((theta-tht(0))*dtht1) & it2=it1+1
      if (it2 eq Ntht) then it2=Ntht-1
      if (it1 lt    0) then begin
        print,'it1 lt 0. theta,it1=',theta,it1
        stop
      endif
        
      if (it2 gt Ntht) then begin
        print,'it2 ge Nt. theta,it2=',theta,it2
        stop
      endif
      delt=theta-tht[it1]
;
            ;take closest r and tht
;
      p1=field(ir1,it1)&p2=field(ir2,it1)
      p3=field(ir1,it2)&p4=field(ir2,it2)
;
      fr=delr*drad1
      ft=delt*dtht1
;
      fieldxz[ix,iz] = fr*ft*(p1-p2-p3+p4) + fr*(p2-p1) + ft*(p3-p1) + p1
    endif else begin
      fieldxz[ix,iz]=0.
    endelse
  endfor
endfor
;
if (keyword_set(plot)) then begin
    contour,fieldxz,xc,zc,/fill,lev=grange(amin,amax,ncolors),/iso,xs=3,ys=3,$
      title='Cartesian plot',xtitle='X',ytitle='Z'
endif
;
fc={field:fieldxz,$
   xc:    xc    ,$
   zc:    zc     }
;
return,fc
;
end
