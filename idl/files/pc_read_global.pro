;
; $Id: pc_read_global.pro,v 1.5 2006-06-02 19:10:17 joishi Exp $
;
;   Read global variable from file.
;  
pro pc_read_global, gvar, varfile=varfile, datadir=datadir, $
    dim=dim, param=param, TRIMALL=TRIMALL, QUIET=QUIET,     $
    SWAP_ENDIAN=SWAP_ENDIAN
COMPILE_OPT IDL2,HIDDEN
  common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
  common pc_precision, zero, one
; Default data directory

default, datadir, 'data'
default, varfile, 'np.dat'

; Get necessary dimensions, inheriting QUIET
if (n_elements(dim) eq 0) then pc_read_dim,object=dim,datadir=datadir,proc=proc,quiet=quiet
if (n_elements(param) eq 0) then pc_read_param,object=param,dim=dim,datadir=datadir,QUIET=QUIET 

if (n_elements(proc) eq 1) then begin
  procdim=dim
endif else begin
  pc_read_dim,object=procdim,datadir=datadir,proc=0,QUIET=QUIET
endelse

; and check precision is set
pc_set_precision,dim=dim,quiet=quiet

if keyword_set(TRIMALL) then TRIMXYZ=1L

nx=dim.nx
ny=dim.ny
nz=dim.nz
nw=dim.nx*dim.ny*dim.nz
mx=dim.mx
my=dim.my
mz=dim.mz
mvar=dim.mvar
precision=dim.precision
mxloc=procdim.mx
myloc=procdim.my
mzloc=procdim.mz

if (n_elements(proc) eq 1) then ncpus=1 else ncpus = dim.nprocx*dim.nprocy*dim.nprocz

;
; Initialize / set default returns for ALL variables
;

t=zero
x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
dx=zero &  dy=zero &  dz=zero & deltay=zero

if (n_elements(proc) ne 1) then begin
  xloc=fltarr(procdim.mx)*one & yloc=fltarr(procdim.my)*one & zloc=fltarr(procdim.mz)*one
endif
;
;  Read data
;
if (varfile eq 'np.dat') then begin
  varname='np'
  vartype='s'
endif else if (varfile eq 'uup.dat') then begin
  varname='uup'
  vartype='v'
endif

if (n_elements(proc) eq 1) then begin
endif else begin
  if (vartype eq 's') then begin
    res=varname+'=fltarr(mx,my,mz)*one'
  endif else if (vartype eq 'v') then begin
    res=varname+'=fltarr(mx,my,mz,3)*one'
  endif else begin
    print, 'pc_read_global: no such variable type vartype=', vartype
    stop
  endelse
  if (execute(res) ne 1) then message, 'ERROR: could not define variable array'
endelse



; Get a unit number
GET_LUN, file

for i=0,ncpus-1 do begin 
  if (n_elements(proc) eq 1) then begin
    ; Build the full path and filename
    filename=datadir+'/proc'+str(proc)+'/'+varfile 
  endif else begin
    filename=datadir+'/proc'+str(i)+'/'+varfile 
      pc_read_dim, object=procdim, datadir=datadir, proc=i, QUIET=QUIET 
  endelse
; Check for existence
  dummy=findfile(filename, COUNT=countfile)
  if (not countfile gt 0) then begin
    message, 'ERROR: cannot find file '+ filename
  endif
;
  close, file
  openr, file, filename, /F77,SWAP_ENDIAN=SWAP_ENDIAN

  if (execute('readu, file, '+varname) ne 1) then $
      message, 'Error reading: ' + 'readu, file, '+varname
;
  if (n_elements(proc) eq 1) then begin
    if (param.lshear) then begin
      readu, file, t, x, y, z, dx, dy, dz, deltay
    endif else begin
      readu, file, t, x, y, z, dx, dy, dz
    endelse
  endif else begin
    if (param.lshear) then begin
      readu, file, t, xloc, yloc, zloc, dx, dy, dz, deltay
    endif else begin
      readu, file, t, xloc, yloc, zloc, dx, dy, dz
    endelse
;
;  Don't overwrite ghost zones of processor to the left (and
;  accordingly in y and z direction makes a difference on the
;  diagonals)
;
;    if (procdim.ipx eq 0L) then begin
;      i0x=0L
;      i1x=i0x+procdim.mx-1L
;      i0xloc=0L 
;      i1xloc=procdim.mx-1L
;    endif else begin
;      i0x=procdim.ipx*procdim.nx+procdim.nghostx 
;      i1x=i0x+procdim.mx-1L-procdim.nghostx
;      i0xloc=procdim.nghostx & i1xloc=procdim.mx-1L
;    endelse
;    
;    if (procdim.ipy eq 0L) then begin
;      i0y=0L
;      i1y=i0y+procdim.my-1L
;      i0yloc=0L 
;      i1yloc=procdim.my-1L
;    endif else begin
;      i0y=procdim.ipy*procdim.ny+procdim.nghosty 
;      i1y=i0y+procdim.my-1L-procdim.nghosty
;      i0yloc=procdim.nghosty 
;      i1yloc=procdim.my-1L
;    endelse
    
;    if (procdim.ipz eq 0L) then begin
;      i0z=0L
;      i1z=i0z+procdim.mz-1L
;      i0zloc=0L 
;      i1zloc=procdim.mz-1L
;    endif else begin
;      i0z=procdim.ipz*procdim.nz+procdim.nghostz 
;      i1z=i0z+procdim.mz-1L-procdim.nghostz
;      i0zloc=procdim.nghostz 
;      i1zloc=procdim.mz-1L
;    endelse
;
;    x[i0x:i1x] = xloc[i0xloc:i1xloc]
;    y[i0y:i1y] = yloc[i0yloc:i1yloc]
;    z[i0z:i1z] = zloc[i0zloc:i1zloc]
;
;    
;    cmd =   varname $
;          + "[i0x:i1x,i0y:i1y,i0z:i1z,*,*]=" $
;          + varname+'_loc' $
;          +"[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*,*]"
;    if (execute(cmd) ne 1) then $
;        message, 'Error combining data for ' + varname
;
    close,file
    FREE_LUN,file
  endelse
endfor

; Tidy memory a little
if (n_elements(proc) ne 1) then begin
  undefine,xloc
  undefine,yloc
  undefine,zloc
endif

if (keyword_set(TRIMXYZ)) then begin
  x=x[dim.l1:dim.l2]
  y=y[dim.m1:dim.m2]
  z=z[dim.n1:dim.n2]
endif

if keyword_set(TRIMALL) then begin
  res = varname+'=pc_noghost('+varname+')'
  if (execute(res) ne 1) then print, 'pc_read_global: problems with trimming'
endif

; If requested print a summary
;if keyword_set(STATS) or (not (keyword_set(NOSTATS) or keyword_set(QUIET))) then begin
;  pc_object_stats,object,dim=dim,QUIET=QUIET
  print,' t = ', t
;endif

res='gvar='+varname
if (execute(res) ne 1) then print, 'pc_read_global: problem'


end
