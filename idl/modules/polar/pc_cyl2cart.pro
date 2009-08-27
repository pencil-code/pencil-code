
function pc_cyl2cart,field_in,rad,phi,perc_expand=perc_expand,time=time,$
                  plot=plot,amin=amin,amax=amax,ncolors=ncolors,$
                  keep_resolution=keep_resolution   

;;;
;;; pc_cyl2cart.pro
;;;
;;; Author: Wladimir Lyra
;;; Date: 2008/03/15
;;;
;;; Description:
;;;   Transform a 2d variable from Cylindrical to Cartesian
;;;   coordinates, for visualization purposes. It returns a
;;;   structure with the variable in Cartesian coordinates 
;;;   (field), and the cartesian axes xc and yc.
;;;   
;;;   
;;; Usage:
;;; 
;;;     IDL> rad=ff.x
;;; 
;;; In Cylindrical coordinates
;;;
;;;     IDL> phi=ff.y
;;;     IDL> midplane=reform(ff.lnrho(*,*,npoint))	
;;;
;;; In Spherical coordinates
;;;
;;;     IDL> phi=ff.z
;;;     IDL> midplane=reform(ff.lnrho(*,mpoint,*))
;;; 
;;; Transform to Cartesian
;;; 
;;;     IDL> fc=pc_cyl2cart(midplane,rad,phi) 
;;;   or
;;;     IDL> fc=pc_cyl2cart(midplane,rad,phi,perc_expand=10)
;;;   then
;;;     IDL> tmp=minmax(fc.field)
;;;     IDL> contour, fc.field,fc.xc,fc.yc,/fill,lev=grange(tmp(0),tmp(1),256),/iso,xs=3,ys=3
;;;   use
;;;     IDL> help, fc, /STRUCT
;;;   for a full list of slots
;;; 
;;; 
;;;   or for a direct plot,
;;;
;;;     IDL> fc=pc_cyl2cart(midplane,rad,phi,amin=0.0,amax=2.0,ncolors=256,/plot) 
;;;    
;
default,perc_expand,5
default,amin,0.
default,amax,1.
default,ncolors,256
;
field=field_in
;
nr=n_elements(rad)
nphi=n_elements(phi)
;
rext=rad(nr-1)
rint=rad(0)
Lr=rext-rint
;
xn=(1+perc_expand/100.)*rext & x0=-xn
yn=xn                        & y0=-yn
Lx=xn-x0
;
fac=Lx/Lr
;
nxc=round(fac*nr) & nyc=nxc
;
xc=grange(x0,xn,Nxc)
yc=grange(y0,yn,Nyc)
;
fieldxy=fltarr(Nxc,Nyc)
fieldxy=fieldxy*0.
;
drad=rad(1)-rad(0) & drad1=1./drad
dphi=phi(1)-phi(0) & dphi1=1./dphi
;
if (keyword_set(plot)) then begin
    !p.multi=[0,2,1]
    contour,field,rad,phi,/fill,lev=grange(amin,amax,ncolors),xs=3,ys=3,$
      title='Cylindrical plot',xtitle='!8r!x',ytitle=textoidl("\phi")
endif
;
;  Transform to Cartesian with Bilinear interpolation
;
for ix=0,nxc-1 do begin
    for iy=0,nyc-1 do begin
        radius = sqrt(xc[ix]^2+yc[iy]^2)
        if ((radius ge rint) and (radius le rext)) then begin
            ir1=floor((radius-rint)*drad1) & ir2=ir1+1
            delr = radius - rad[ir1]
;
            if (ir1 lt  0) then begin
                print,'ir1 lt 0. ir1=',ir1
                stop
            endif
            if (ir2 ge Nr) then begin
                print,'ir2 ge Nr. ir2=',ir2
                stop
            endif
;
            phis=atan(yc[iy],xc[ix])
            ;take closest r and phi
            ip1=floor((phis-phi(0))*dphi1) & ip2=ip1+1
            if (ip1 eq -1) then begin
                ip1=nphi-1
                delp=phis - phi(ip1) + 2*!pi
            endif else begin
                delp=phis - phi(ip1)
            endelse
            if (ip2 eq nphi) then ip2=1
;
            p1=field(ir1,ip1)&p2=field(ir2,ip1)
            p3=field(ir1,ip2)&p4=field(ir2,ip2)
;
            fr=delr*drad1
            fp=delp*dphi1
;
            fieldxy[ix,iy] = fr*fp*(p1-p2-p3+p4) + fr*(p2-p1) + fp*(p3-p1) + p1
        endif else begin
            fieldxy[ix,iy]=0.
        endelse
    endfor
endfor
;
if (keyword_set(plot)) then begin
    contour,fieldxy,xc,yc,/fill,lev=grange(amin,amax,ncolors),/iso,xs=3,ys=3,$
      title='Cartesian plot',xtitle='X',ytitle='Y'
endif
;
; Rebin xy to the original resolution rphi
;
if (keyword_set(keep_resolution)) then begin
    field=rebin(fieldxy,nr,nphi)
    x=rebin(xc,nr) & y=rebin(yc,nphi)
endif else begin
    field=fieldxy
    x=xc & y=yc
endelse
;
fc={field:field,$
   xc:    x    ,$
   yc:    y      }
;
return,fc
;
end
