;;;;;;;;;;;;;;;;;;;
;;;  start.pro  ;;;
;;;;;;;;;;;;;;;;;;;
;;;
;;; Initialise coordinate arrays, detect precision and dimensions.
;;; Typically run only once before running `r.pro' and other
;;; plotting/analysing scripts.
;;; $Id$

function param
  COMPILE_OPT IDL2,HIDDEN 
; Dummy to keep IDL from complaining. The real param() routine will be
; compiled below
end

common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
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
if (n_elements(datatopdir) eq 0) then datatopdir=pc_get_datadir()
default, varfile, 'var.dat'
default, dimfile, 'dim.dat'
datadir = datatopdir+'/proc'+str(proc)
; Directory for temporary output by the IDL scripts. Defaults to
; data/, but can be overwritten in case you don't have write access to
; data/ (e.g. when working in another user's directory or on CD data):
default, tmpdir, datatopdir
;; Max # of structure tags in one execute statement; if we use more,
;; we get `% Program code area full' error.
;; Typically 398 or 569, but stupid Solaris doesn't tolerate more than 127.
default, maxtags, 120
;
;  Read the dimensions and precision (single or double) from dim.dat
;
mx=0L & my=0L & mz=0L & nvar=0L & naux=0L & nglobal=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
close,1
openr,1,datadir+'/'+dimfile
readf,1,mx,my,mz,nvar,naux,nglobal
readf,1,prec
readf,1,nghostx,nghosty,nghostz
close,1
;
mw=mx*my*mz  ;(this must be calculated; its not in dim.dat)
prec = (strtrim(prec,2))        ; drop leading zeros
prec = strmid(prec,0,1)
if ((prec eq 'S') or (prec eq 's')) then begin ; single precision
  one = 1.e0
  nl2idl_d_opt = ''
endif else if ((prec eq 'D') or (prec eq 'd')) then begin ; double precision
  one = 1.D0
  nl2idl_d_opt = '-d'
endif else begin
  if (quiet le 4) then print, "prec = `", prec, "' makes no sense to me"
  STOP
endelse
zero = 0*one

;  Read the positions of variables in f from index.pro
;  Can't just use `@data/index', as the data directory may have a different name
;  [We used to concatenate all lines in a Perl line and execute this, but
;   by now on some systems this blows the limit on the number of commands
;   that can be concatenated with '&' (445 on x86-Linux, IDL6.0; see also
;   the comment to `maxtags' above)]
openr, 1, datatopdir+'/index.pro', ERROR=err
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

;
;  Read grid
;
t=zero
x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
dx=zero &  dy=zero &  dz=zero & dxyz=zero
gridfile=datadir+'/'+'grid.dat'
dummy=file_search(gridfile, COUNT=cgrid)
dummy2=zero
if (cgrid gt 0) then begin
  if (quiet le 2) then print, 'Reading grid.dat..'
  openr,1, gridfile, /F77
  readu,1, t,x,y,z
  readu,1, dx,dy,dz
  readu,1, dummy2,dummy2,dummy2 ; instead of Lx,Ly,Lz
  point_lun,-1,pos              ; to read nonequidistant stuff below 
  close,1
endif else begin
  if (quiet le 4) then print, 'Warning: cannot find file ', gridfile
endelse
;
;print,'calculating xx,yy,zz (comment this out if there is not enough memory)'
;xx = spread(x, [1,2], [my,mz])
;yy = spread(y, [0,2], [mx,mz])
;zz = spread(z, [0,1], [mx,my])
;
;  set boundary values for physical (sub)domain
;  l12 etc are useful as array indices, e.g.
;    plot_3d_vect, bb[l12,m12,nz/2,*],x[l12],y[m12]
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
l1=nghostx & l2=mx-nghostx-1 & l12=l1+indgen(nx)
m1=nghosty & m2=my-nghosty-1 & m12=m1+indgen(ny)
n1=nghostz & n2=mz-nghostz-1 & n12=n1+indgen(nz)
;
;  Read startup parameters
;
pfile = datatopdir+'/param.nml'
dummy = file_search(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  if (quiet le 2) then print, 'Reading param.nml..'
  spawn, 'for d in . $TMPDIR $TMP /tmp /var/tmp; do if [ -d $d -a -w $d ]; then echo $d; fi; done', /SH, result
  if (strlen(result[0])) le 0 then begin
    message, "Can't find writeable directory for temporary files"
  endif else begin
    tmpdir = result[0]
  endelse
  tmpfile = tmpdir+'/param.pro'
  ;; Write content of param.nml to temporary file:
  spawn, '$PENCIL_HOME/bin/nl2idl '+nl2idl_d_opt+' -m ' + $
      datatopdir+'/param.nml > ' + tmpfile , result
  ;; Compile that file. Should be easy, but is incredibly awkward, as
  ;; there is no way in IDL to compile a given file at run-time
  ;; outside the command line:
  ;; Save old path and pwd
  _path = !path
  cd, tmpdir, CURRENT=_pwd
  !path = '.:'
  resolve_routine, 'param', /IS_FUNCTION
  ;; Restore old path and pwd
  !path = _path & cd, _pwd
  ;; Delete temporary file
  ; file_delete, tmpfile      ; not in IDL <= 5.3
  spawn, 'rm -f '+tmpfile, /SH
  par = param()

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
  lgravz    = par.lgravz
  lgravr    = par.lgravr
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
  lradiation_fld = par.lradiation_fld
  coord_system  = par.coord_system

  ;
  ;  Read coefficients for nonequidistant grid
  ;
  dx_1=fltarr(mx)*zero & dy_1=fltarr(my)*zero & dz_1=fltarr(mz)*zero
  dx_tilde=fltarr(mx)*zero& dy_tilde=fltarr(my)*zero& dz_tilde=fltarr(mz)*zero
  if (any(lequidist eq 0)) then begin
    openr,1,gridfile,/F77
    point_lun,1,pos
    readu,1, dx_1,     dy_1,     dz_1
    readu,1, dx_tilde, dy_tilde, dz_tilde
    close,1
  endif else begin
    ;
    ;  Ensure we don't use these values
    ;
    dx_1 = dx_1*!values.f_nan
    dy_1 = dy_1*!values.f_nan
    dz_1 = dz_1*!values.f_nan
    dx_tilde = dx_tilde*!values.f_nan
    dy_tilde = dy_tilde*!values.f_nan
    dz_tilde = dz_tilde*!values.f_nan
  endelse
  ;
  if (ldensity and (lentropy or ltemperature)) then begin
     if (not (leos_ionization or leos_temperature_ionization)) then begin
      cs0=par.cs0 & rho0=par.rho0
      gamma=par.gamma & gamma_m1=gamma-1.
      cs20 = cs0^2 & lnrho0 = alog(rho0)
    endif
  endif
  ;
  if (lgravz) then begin
    z1=par.z1 & z2=par.z2
    zref=par.zref
    gravz=par.gravz
    ztop=z[n2] & z3=ztop
  endif
  ;
  if (lentropy) then begin
    mpoly0=par.mpoly0 & mpoly1=par.mpoly1
    mpoly2=par.mpoly2 & isothtop=par.isothtop
  endif
endif else begin
  if (quiet le 4) then print, 'Warning: cannot find file ', pfile
endelse

if (quiet le 2) then begin
  print, 'To get index arrays for accessing'
  print, '  scalar[lmn12] \equiv scalar[l1:l2,m1:m2,n1:n2] and'
  print, '  vector[lmn123] \equiv vector[l1:l2,m1:m2,n1:n2,*] type'
  print,'lmn12 = l1+spread(indgen(nx),[1,2],[ny,nz]) + mx*(m1+spread(indgen(ny),[0,2],[nx,nz])) + mx*my*(n1+spread(indgen(nz),[0,1],[nx,ny]))'
  print, 'lmn123 = spread(lmn12,3,3) + (mx*my*mz)*spread(indgen(3),[0,1,2],[nx,ny,nz])'
endif
;
if (quiet le 2) then print, '..done'
;
started=1
END
