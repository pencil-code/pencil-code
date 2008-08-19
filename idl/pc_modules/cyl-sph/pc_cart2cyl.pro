
function pc_cart2cyl,field_in,x,y,rint,perc_expand=perc_expand,time=time,planet=planet,plot=plot,$
                  amin=amin,amax=amax,ncol=ncol,shift=shift,minimum_floor=minimum_floor

;;;
;;; pc_cart2cyl.pro
;;;
;;; Author: Wladimir Lyra
;;; Date: 2008/03/15
;;;
;;; Description:
;;;   Transform a 2d variable from Cartesian to cylindrical
;;;   coordinates, for visualization purposes. In addition to
;;;   the coordinates, one has to provide rint, the start of the
;;;   constructed cylindrical grid, since the Cartesian goes 
;;;   all the way to the origin.
;;;   
;;;   
;;; Usage:
;;;     IDL> fc=pc_cart2cyl(ff.lnrho(*,*,npoint),ff.x,ff.y,rint) 
;;;   or
;;;     IDL> fc=pc_cart2cyl(ff.lnrho(*,*,npoint),ff.x,ff.y,rint,perc_expand=10)
;;;   or
;;;     IDL> fc=pc_cart2cyl(ff.lnrho(*,*,npoint),ff.x,ff.y,rint,perc_expand=10,/planet) 
;;;   then
;;;     IDL> tmp=minmax(fc.field)
;;;     IDL> contour, fc.field,fc.xc,fc.yc,/fill,lev=grange(tmp(0),tmp(1),256),/iso,xs=3,ys=3
;;;   use
;;;     IDL> help, fc, /STRUCT
;;;   for a full list of slots
;;;
;;;
;;;   or for a direct plot,
;;;     IDL> fc=pc_cart2cyl(ff.lnrho,ff.x,ff.y,rint,time=ff.t,/planet,amin=0.0,amax=2.0,ncol=256,/plot,shift=0) 
;;;    


default,perc_expand,5
default,time,0.
default,amin,0.
default,amax,1.
default,ncol,256
default,shift,-!pi
default,minimum_floor,0.

field=field_in

nx=n_elements(x)
ny=n_elements(y)

rext=x(nx-1)

nnc= nx > ny

nr=nnc       & nphi=nr

rad=grange(rint,rext,Nr)
phi=grange(-!pi,!pi,Nphi)

fieldrp=fltarr(Nr,Nphi)
fieldrp=fieldrp*0.

dx=x(1)-x(0)
dy=y(1)-y(0)
dx1=1./dx
dy1=1./dy

for ir=0,Nr-1 do begin
    for iphi=0,Nphi-1 do begin

        distx = rad[ir]*cos(phi[iphi]) - x[0]
        disty = rad[ir]*sin(phi[iphi]) - y[0]

        ix1 = floor(distx*dx1) & ix2 = ceil(distx*dx1)

        iy1 = floor(disty*dy1) & iy2 = ceil(disty*dy1)

        if (ix1 lt 0 ) then ix1=0
        if (iy1 lt 0 ) then iy1=0
        if (ix2 ge Nx) then ix2=Nx-1
        if (iy2 ge Ny) then iy2=Ny-1
        
        dudxdo   = (field[ix2,iy1] - field[ix1,iy1])*dx1
        dudxup   = (field[ix2,iy2] - field[ix1,iy2])*dx1

        delx = rad[ir]*cos(phi[iphi]) - x[ix1]

        fieldx_do = field[ix1,iy1] + dudxdo*delx
        fieldx_up = field[ix1,iy2] + dudxup*delx
        fielddudy = (fieldx_up - fieldx_do)*dy1

        dely = rad[ir]*sin(phi[iphi]) - y[iy1]

        fieldrp[ir,iphi] = fieldx_do + fielddudy*dely

    endfor
endfor

field=rebin(fieldrp,nx,ny)
rad=rebin(rad,nx)
phi=rebin(phi,ny)


fc={field:field,$
   rad:    rad    ,$
   phi:    phi      }

return,fc


end
