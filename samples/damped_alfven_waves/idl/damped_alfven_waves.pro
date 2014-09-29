; $Id$
;+
; NAME:
;	DAMPED_ALFVEN_WAVES
;
; PURPOSE:
;	This procedure compares the analytical solution of a damped
;	Alfven wave with the numerical solution calculated by the Pencil
;	Code.
;
; CATEGORY:
;	The Pencil Code - Code Test
;
; CALLING SEQUENCE:
;	DAMPED_ALFVEN_WAVES [, /DOUBLE] [, VARFILE=string]
;
; KEYWORDS:
;	DOUBLE:	Set this keyword if the numerical data is in double
;		precision.
;	VARFILE:	Set this keyword to a string containing the name
;		of the snapshot.  If this keyword is omitted, the
;		snapshot 'var.dat' will be read.
;
; EXAMPLE:
;	DAMPED_ALFVEN_WAVES, /DOUBLE, VARFILE='VAR10'
;
; MODIFICATION HISTORY:
; 	Written by:	Chao-Chin Yang, August 22, 2013.
;-

pro damped_alfven_waves, double=double, varfile=varfile
  compile_opt idl2

; Read the data.
  pc_read_dim, obj=dim, /quiet
  pc_read_param, obj=par1, /quiet
  pc_read_param, obj=par2, /param2, /quiet
  pc_read_var, obj=f, varfile=varfile, /bb, /trimall, /quiet

; Calculate relevant parameters.
  ki = par1.init_k0
  nu = par2.nu
  eta = par2.eta
  omega = transpose(ki) # par2.b_ext / sqrt(par1.mu0 * par1.rho0)
  ot = omega * f.t
  decay = exp(-0.5 * total(ki^2, double=double) * (nu + eta) * f.t)
  du = [par1.init_amp0[0], par1.init_amp0[1], par1.init_amp0[2]] * decay
  db = du * sqrt(par1.mu0 * par1.rho0)
  du2 = total(du^2, double=double)
  db2 = total(db^2, double=double)

; Allocate memory.
  if keyword_set(double) then begin
    uex = dblarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    uey = dblarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    uez = dblarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    bex = dblarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    bey = dblarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    bez = dblarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
  endif else begin
    uex = fltarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    uey = fltarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    uez = fltarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    bex = fltarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    bey = fltarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
    bez = fltarr(dim.nxgrid,dim.nygrid,dim.nzgrid)
  endelse

; Find the exact solution for the velocity and magnetic field.
  for k = 0, dim.nzgrid - 1 do begin
    for j = 0, dim.nygrid - 1 do begin
      phase = ki[0] * f.x + (ki[1] * (f.y)[j] + ki[2] * (f.z)[k] + ot)[0]
      a = sin(phase)
      uex[*,j,k] = du[0] * a
      uey[*,j,k] = du[1] * a
      uez[*,j,k] = du[2] * a
      bex[*,j,k] = db[0] * a
      bey[*,j,k] = db[1] * a
      bez[*,j,k] = db[2] * a
    endfor
  endfor

; Compare the solutions.
  ix = fix(dim.nxgrid * randomu(seed))
  iy = fix(dim.nygrid * randomu(seed))
  iz = fix(dim.nzgrid * randomu(seed))
  bmax = max([max(bex), max(bey), max(bez)])
  plot, f.x, bex[*,iy,iz], yrange=[-bmax,bmax]
  oplot, f.x, bey[*,iy,iz], color=100
  oplot, f.x, bez[*,iy,iz], color=200
  oplot, f.x, f.bb[*,iy,iz,0], psym=1
  oplot, f.x, f.bb[*,iy,iz,1], psym=1, color=100
  oplot, f.x, f.bb[*,iy,iz,2], psym=1, color=200

; Evaluate the errors.
  vol = (f.dx / par1.lxyz[0]) * (f.dy / par1.lxyz[1]) * (f.dz / par1.lxyz[2])
  err2 = ((f.uu[*,*,*,0] - uex)^2 + (f.uu[*,*,*,1] - uey)^2 + (f.uu[*,*,*,2] - uez)^2) / du2 + $
         ((f.bb[*,*,*,0] - bex)^2 + (f.bb[*,*,*,1] - bey)^2 + (f.bb[*,*,*,2] - bez)^2) / db2
  two_norm = sqrt(vol * total(err2, double=double))
  max_norm = sqrt(max(err2))
  print, 'dx = ', f.dx / par1.lxyz[0]
  print, 'Two-norm error = ', two_norm
  print, 'Max-norm error = ', max_norm

end
