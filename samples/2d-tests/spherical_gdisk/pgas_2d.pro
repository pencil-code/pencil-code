;
; $Id: pgas_2d.pro,v 1.1 2008-03-21 23:13:06 wlyra Exp $
;

!p.charsize=3.0
!p.multi=[0,2,4]

pc_read_var,obj=ff,dim=dim
pc_read_var,obj=ff0,ivar=0

pc_read_param,obj=par
rad=ff.x & tht=ff.y & phi=ff.z & uu=ff.uu
if (par.ldensity_nolog eq 0) then begin
    rho=exp(ff.lnrho)
    rho0=exp(ff0.lnrho)
endif else begin
    rho=ff.lnrho
    rho0=ff0.lnrho
endelse


rint = par.r_int
rext = par.r_ext
lr = par.Lxyz(0)/2.

Mx=dim.mxgrid & My=dim.mygrid & Mz=dim.mzgrid
Nx=dim.nxgrid & Ny=dim.nygrid & Nz=dim.nzgrid

print,'Mx,My,Mz=',Mx,My,Mz

loadct,3

lpoint = fix(Mx*0.5)
mpoint = fix(My*0.5)
npoint = fix(Mz*0.5)

print,'xp,yp,np'
print,rad(lpoint),tht(mpoint),phi(npoint)

l1=dim.l1 & l2=dim.l2
m1=dim.m1 & m2=dim.m2
n1=dim.n1 & n2=dim.n2

n=npoint

urad   = reform(uu(*,*,*,0)) 
utht   = reform(uu(*,*,*,1)) 
uphi   = reform(uu(*,*,*,2)) 

lev=grange(0.,2.,100)

contour,reform(rho(*,*,npoint)),rad,tht/!pi*180-90,/fill,lev=lev, $
  xtitle='r',ytitle='tht (degrees)',title='Density - Polar',xs=3,ys=3

;aa=minmax(ff.shock(*,*,npoint))
;levv=grange(aa(0),aa(1),20)
;loadct,5
;contour,reform(ff.shock(*,*,npoint)),rad,tht/!pi*180-90,/fill,lev=levv, $
;  xtitle='r',ytitle='tht (degrees)',title='shock',xs=3,ys=3
;loadct,3
ratio=rho/rho0

contour,reform(ratio(*,*,npoint)),rad,tht/!pi*180-90,/fill,lev=lev, $
  xtitle='r',ytitle='tht (degrees)',title='Density - Polar',xs=3,ys=3

plot,rad,reform(rho(*,mpoint,npoint))$
  ,xstyle=3,ystyle=3,title="Rad Density",yr=[0,2]
plot,tht,reform(rho(lpoint,*,npoint))$
  ,/ylog,title="Polar Density",ys=3,xs=3

plot,rad,reform(urad(*,mpoint,npoint))$
  ,xstyle=3,ystyle=3,title="Rad Vel",yr=[0,2]
plot,rad,reform(utht(lpoint,*,npoint))$
  ,title="Polar Vel",ys=3,xs=3
plot,rad,reform(uphi(*,mpoint,npoint))$
  ,xstyle=3,ystyle=3,title="Azimuthal Vel",yr=[0,2]


!p.multi=0


END
