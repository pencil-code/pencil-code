;;;;;;;;;;;;;;;;;;;
;;;  start.pro  ;;;
;;;;;;;;;;;;;;;;;;;
;;;
;;; Initialise coordinate arrays, detect precision and dimensions.
;;; Typically run only once before running `r.pro' and other
;;; plotting/analysing scripts.
;;; $Id$

common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
common cdat_coords, coord_system
;forward_function safe_get_tag
;
;  Compile the derivative routines for data that have ghost zones
;  For analysis purposes, one may want to use other routines (called
;  xder_6th, yder_6th, ..., zder2_6th in Axel's idl/pro repo).
;
@xder_6th_ghost
@yder_6th_ghost
@zder_6th_ghost
@xder2_6th_ghost
@yder2_6th_ghost
@zder2_6th_ghost
;
;  The following avoids a mysterious bug when using esrg_legend later
;  (box size was wrong, because lenstr(['a','b']) would be wrong,
;  because xyout would write all letters onto one posision) ..
;
@lenstr
;
;; Flag for controlling inverse verbosity
;; quiet=0: print everything (default)
;;       1: print mildly interesting messages
;;       2: print interesting, but irrelevant info (status of reading, etc.)
;;       3: print important info
;;       4: print warnings about inconsistencies and unexpected behavior
;;       5: don't print anything (to the extent this is possible)
default, quiet, 0
default, proc, 0
default, varfile, 'var.dat'

;
;  Read the grid, dimensions and startup parameters
;
pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=par

l12=dim.l1+indgen(dim.nx)
m12=dim.m1+indgen(dim.ny)
n12=dim.n1+indgen(dim.nz)

;  Read the positions of variables in f from index.pro
;  Can't just use `@data/index', as the data directory may have a different name
;  [We used to concatenate all lines in a Perl line and execute this, but
;   by now on some systems this blows the limit on the number of commands
;   that can be concatenated with '&' (445 on x86-Linux, IDL6.0; see also
;   the comment to `maxtags' above)]
openr, 1, datadir+'/index.pro', ERROR=err
if (err ne 0) then begin
  free_lun, in_file
  message, !err_string
endif
;
line = ''
lineno = 0
while (not eof(1)) do begin
  readf, 1, line
  lineno = lineno + 1
  if (execute(line) ne 1) then $
      message, /INFO, $
               'There was a problem with index.pro, line ' $
               + strtrim(lineno,2) + ': ' + line
endwhile
close, 1
;
if (quiet le 2) then print,'nname=',nname

  ;; Abbreviate some frequently used parameters
  x0=par.xyz0[0] & y0=par.xyz0[1] & z0=par.xyz0[2]
  Lx=par.Lxyz[0] & Ly=par.Lxyz[1] & Lz=par.Lxyz[2]
  unit_system=par.unit_system
  unit_length=par.unit_length
  unit_velocity=par.unit_velocity
  unit_density=par.unit_density
  unit_temperature=par.unit_temperature
  ;
  default, STRUCT=par, ['leos_ionization','leos_fixed_ionization'],  0L
  default, STRUCT=par, ['leos_temperature_ionization'],  0L
  default, STRUCT=par, 'lequidist', [-1L, -1L, -1L]
  lequidist = par.lequidist
  lhydro    = par.lhydro
  ldensity  = par.ldensity
  lentropy  = par.lentropy
  ltemperature = par.ltemperature
  lmagnetic = par.lmagnetic
  lradiation= par.lradiation
  leos_ionization=par.leos_ionization
  leos_temperature_ionization=par.leos_temperature_ionization
  leos_fixed_ionization=par.leos_fixed_ionization
  ;lvisc_shock=par.lvisc_shock
  ;lvisc_hyper3=par.lvisc_hyper3
  lpscalar  = par.lpscalar
  ldustvelocity = par.ldustvelocity
  ldustdensity = par.ldustdensity
  lforcing  = par.lforcing
  lshear    = par.lshear
  coord_system  = par.coord_system

  if (ldensity and (lentropy or ltemperature)) then begin
     if (not (leos_ionization or leos_temperature_ionization)) then begin
      cs0=par.cs0 & rho0=par.rho0
      gamma=par.gamma & gamma_m1=gamma-1.
      cs20 = cs0^2 & lnrho0 = alog(rho0)
    endif
  endif

  if (lentropy) then begin
    mpoly0=par.mpoly0 & mpoly1=par.mpoly1
    mpoly2=par.mpoly2 & isothtop=par.isothtop
  endif

if (quiet le 2) then begin
  print, 'To get index arrays for accessing'
  print, '  scalar[lmn12] \equiv scalar[l1:l2,m1:m2,n1:n2] and'
  print, '  vector[lmn123] \equiv vector[l1:l2,m1:m2,n1:n2,*] type'
  print, 'lmn12 = l1+spread(indgen(nx),[1,2],[ny,nz]) + mx*(m1+spread(indgen(ny),[0,2],[nx,nz])) + mx*my*(n1+spread(indgen(nz),[0,1],[nx,ny]))'
  print, 'lmn123 = spread(lmn12,3,3) + (mx*my*mz)*spread(indgen(3),[0,1,2],[nx,ny,nz])'
endif
;
if (quiet le 2) then print, '..done'
;
started=1
END
