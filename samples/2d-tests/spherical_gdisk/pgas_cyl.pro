;
; $Id: pgas_cyl.pro,v 1.1 2008-03-21 23:13:06 wlyra Exp $
;

!p.charsize=2.0
!p.multi=[0,2,2]

pc_read_var,obj=ff,dim=dim
;pc_read_var,obj=df,varfile='dvar.dat'

pc_read_param,obj=par
rad=ff.x & tht=ff.y & phi=ff.z 
rho=ff.lnrho & uu=ff.uu
;drho=df.lnrho & duu=df.uu



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

urad   = reform(uu(*,mpoint,*,0)) 
utht   = reform(uu(*,mpoint,*,1)) 
uphi   = reform(uu(*,mpoint,*,2)) 

;durad   = reform(duu(*,mpoint,*,0)) 
;dutht   = reform(duu(*,mpoint,*,1)) 
;duphi   = reform(duu(*,mpoint,*,2)) 

rho  = reform(rho(*,mpoint,*))

lev=grange(-1,1,100)

contour,alog10(rho),rad,phi,/fill,lev=lev, $
  xtitle='r',ytitle='phi',title='Density Map',xs=3,ys=3

plot,rad,urad(*,npoint)
plot,rad,utht(*,npoint)
plot,rad,uphi(*,npoint)




!p.multi=0


END
