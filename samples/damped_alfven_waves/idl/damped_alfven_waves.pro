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
; 	Written by:	Chao-Chin Yang, October 17, 2012.
;-

pro damped_alfven_waves, double=double, varfile=varfile
  compile_opt idl2

; Read the data.
  pc_read_var, obj=f, varfile=varfile, /bb, /trimall, /quiet
  pc_read_dim, obj=dim, /quiet
  pc_read_param, obj=par1, /quiet
  pc_read_param, obj=par2, /param2, /quiet

; Calculate relevant parameters.
  ki = [par1.kx_ux[0], par1.ky_ux[0], par1.kz_ux[0]]
  omega = transpose(ki) # par2.b_ext / sqrt(par1.mu0 * par1.rho0)
  ot = omega * f.t
  decay = exp(-0.5 * total(ki^2, double=double) * (par2.nu + par2.eta) * f.t)
  du = [par1.ampl_ux[0], par1.ampl_uy[0], par1.ampl_uz[0]] * decay
  db = du * sqrt(par1.mu0 * par1.rho0)

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
      phase = ki[0] * f.x + (ki[1] * f.y[j] + ki[2] * (f.z)[k] + ot)[0]
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
  plot, f.x, bex[*,iy,iz]
  oplot, f.x, bey[*,iy,iz], color=100
  oplot, f.x, bez[*,iy,iz], color=200
  oplot, f.x, f.bb[*,iy,iz,0], psym=1
  oplot, f.x, f.bb[*,iy,iz,1], psym=1, color=100
  oplot, f.x, f.bb[*,iy,iz,2], psym=1, color=200

; Find the 2-norm error.
  vol = f.dx * f.dy * f.dz
  error_uu = [sqrt(vol * total((f.uu[*,*,*,0] - uex)^2, double=double)), $
              sqrt(vol * total((f.uu[*,*,*,1] - uey)^2, double=double)), $
              sqrt(vol * total((f.uu[*,*,*,2] - uez)^2, double=double))]
  error_bb = [sqrt(vol * total((f.bb[*,*,*,0] - bex)^2, double=double)), $
              sqrt(vol * total((f.bb[*,*,*,1] - bey)^2, double=double)), $
              sqrt(vol * total((f.bb[*,*,*,2] - bez)^2, double=double))]
  print, 'dx = ', f.dx
  print, 'uu error = ', error_uu
  print, 'bb error = ', error_bb

end
