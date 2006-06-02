; $Id: pc_read_grid.pro,v 1.14 2006-06-02 19:10:17 joishi Exp $
;
;   Read grid.dat
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2006-06-02 19:10:17 $
;  $Revision: 1.14 $
;
;  27-nov-02/tony: coded 
;
;  
pro pc_read_grid,object=object, dim=dim, param=param, $
                 TRIMXYZ=TRIMXYZ, $
                 datadir=datadir,proc=proc,PRINT=PRINT,QUIET=QUIET,HELP=HELP, $
                 SWAP_ENDIAN=SWAP_ENDIAN
COMPILE_OPT IDL2,HIDDEN
  common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
  common pc_precision, zero, one
; If no meaningful parameters are given show some help!
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_read_grid, t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, $                      "
    print, "              datadir=datadir, proc=proc, $                                   "
    print, "              /PRINT, /QUIET, /HELP                                           "
    print, "                                                                              "
    print, "Returns the grid position arrays and grid deltas of a Pencil-Code run. For a  "
    print, "specific processor. Returns zeros and empty in all variables on failure.      "
    print, "                                                                              "
    print, "  datadir: specify root data directory. Default: './data'          [string]   "
    print, "     proc: specify a processor to get the data from. Default is 0  [integer]  "
    print, "                                                                              "
    print, "   object: optional structure in which to return all the above     [structure]"
    print, "                                                                              "
    print, "   /PRINT: instruction to print all variables to standard output              "
    print, "   /QUIET: instruction not to print any 'helpful' information                 "
    print, "    /HELP: display this usage information, and exit                           "
    return
  ENDIF

; Default data directory

default, datadir, 'data'

; Get necessary dimensions, inheriting QUIET
if n_elements(dim) eq 0 then  $
     pc_read_dim,object=dim,datadir=datadir,proc=proc,QUIET=QUIET 
if n_elements(param) eq 0 then  $
     pc_read_param,object=param,datadir=datadir,QUIET=QUIET 

ncpus=dim.nprocx*dim.nprocy*dim.nprocz

; Set mx,my,mz in common block for derivative routines
mx=dim.mx
my=dim.my
mz=dim.mz
lequidist=safe_get_tag(param,'lequidist',default=[1,1,1])

; and check pc_precision is set!
pc_set_precision,dim=dim,QUIET=QUIET

;
; Initialize / set default returns for ALL variables
;
t=zero
x=fltarr(dim.mx)*one & y=fltarr(dim.my)*one & z=fltarr(dim.mz)*one
dx=zero &  dy=zero &  dz=zero 
Lx=zero &  Ly=zero &  Lz=zero 
dx_1=fltarr(dim.mx)*one & dy_1=fltarr(dim.my)*one & dz_1=fltarr(dim.mz)*one
dx_tilde=fltarr(dim.mx)*one & dy_tilde=fltarr(dim.my)*one & dz_tilde=fltarr(dim.mz)*one

; Get a unit number
GET_LUN, file

for i=0,ncpus-1 do begin
  ; Build the full path and filename
  if n_elements(proc) ne 0 then begin
    filename=datadir+'/proc'+str(proc)+'/grid.dat'   ; Read processor box dimensions
  endif else begin
    filename=datadir+'/proc'+str(i)+'/grid.dat'   ; Read processor box dimensions
    pc_read_dim,object=procdim,datadir=datadir,proc=i,QUIET=QUIET 
    xloc=fltarr(procdim.mx)*one  
    yloc=fltarr(procdim.my)*one 
    zloc=fltarr(procdim.mz)*one
    dx_1loc=fltarr(procdim.mx)*one  
    dy_1loc=fltarr(procdim.my)*one 
    dz_1loc=fltarr(procdim.mz)*one
    dx_tildeloc=fltarr(procdim.mx)*one  
    dy_tildeloc=fltarr(procdim.my)*one 
    dz_tildeloc=fltarr(procdim.mz)*one
    ;
    ;  Don't overwrite ghost zones of processor to the left (and
    ;  accordingly in y and z direction makes a difference on the
    ;  diagonals)
    ;
    if (procdim.ipx eq 0L) then begin
      i0x=0L
      i1x=i0x+procdim.mx-1L
      i0xloc=0L 
      i1xloc=procdim.mx-1L
    endif else begin
      i0x=procdim.ipx*procdim.nx+procdim.nghostx 
      i1x=i0x+procdim.mx-1L-procdim.nghostx
      i0xloc=procdim.nghostx & i1xloc=procdim.mx-1L
    endelse
    ;
    if (procdim.ipy eq 0L) then begin
      i0y=0L
      i1y=i0y+procdim.my-1L
      i0yloc=0L 
      i1yloc=procdim.my-1L
    endif else begin
      i0y=procdim.ipy*procdim.ny+procdim.nghosty 
      i1y=i0y+procdim.my-1L-procdim.nghosty
      i0yloc=procdim.nghosty 
      i1yloc=procdim.my-1L
    endelse
    ;
    if (procdim.ipz eq 0L) then begin
      i0z=0L
      i1z=i0z+procdim.mz-1L
      i0zloc=0L 
      i1zloc=procdim.mz-1L
    endif else begin
      i0z=procdim.ipz*procdim.nz+procdim.nghostz 
      i1z=i0z+procdim.mz-1L-procdim.nghostz
      i0zloc=procdim.nghostz 
      i1zloc=procdim.mz-1L
    endelse
  endelse
  ; Check for existance and read the data
  dummy=findfile(filename, COUNT=countfile)
  if (not countfile gt 0) then begin
    FREE_LUN,file
    message, 'ERROR: cannot find file '+ filename
  endif

  IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' , filename , '...'

  openr,file,filename,/F77,SWAP_ENDIAN=SWAP_ENDIAN
    
  if n_elements(proc) ne 0 then begin
      readu,file, t,x,y,z
      readu,file, dx,dy,dz
      found_Lxyz=0
      found_grid_der=0
      ON_IOERROR, missing
      found_Lxyz=1
      readu,file, Lx,Ly,Lz
      readu,file, dx_1,dy_1,dz_1
      readu,file, dx_tilde,dy_tilde,dz_tilde
  endif else begin
      readu,file, t,xloc,yloc,zloc
      x[i0x:i1x] = xloc[i0xloc:i1xloc]
      y[i0y:i1y] = yloc[i0yloc:i1yloc]
      z[i0z:i1z] = zloc[i0zloc:i1zloc]
      
      readu,file, dx,dy,dz
      found_Lxyz=0
      found_grid_der=0
      ON_IOERROR, missing
      found_Lxyz=1
      readu,file, Lx,Ly,Lz

      readu,file,xloc,yloc,zloc
      dx_1[i0x:i1x] = xloc[i0xloc:i1xloc]
      dy_1[i0y:i1y] = yloc[i0yloc:i1yloc]
      dz_1[i0z:i1z] = zloc[i0zloc:i1zloc]

      readu,file,xloc,yloc,zloc
      dx_tilde[i0x:i1x] = xloc[i0xloc:i1xloc]
      dy_tilde[i0y:i1y] = yloc[i0yloc:i1yloc]
      dz_tilde[i0z:i1z] = zloc[i0zloc:i1zloc]
  endelse
  found_grid_der=1    

missing:
  ON_IOERROR, Null

  close,file 

endfor
FREE_LUN,file

if (keyword_set(TRIMXYZ)) then begin
  x=x[dim.l1:dim.l2]
  y=y[dim.m1:dim.m2]
  z=z[dim.n1:dim.n2]
endif
; Build structure of all the variables
if (found_Lxyz and found_grid_der) then begin
  object = CREATE_STRUCT(name="pc_read_grid_" + $
                          str((size(x))[1]) + '_' + $
                          str((size(y))[1]) + '_' + $
                          str((size(z))[1]), $
                     ['t','x','y','z','dx','dy','dz', $
                      'Lx','Ly','Lz', $
                      'dx_1','dy_1','dz_1', $
                      'dx_tilde','dy_tilde','dz_tilde'], $
                       t,x,y,z,dx,dy,dz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
endif else if (found_Lxyz) then begin
  object = CREATE_STRUCT(name="pc_read_grid_" + $
                          str((size(x))[1]) + '_' + $
                          str((size(y))[1]) + '_' + $
                          str((size(z))[1]), $
                     ['t','x','y','z','dx','dy','dz', $
                      'Lx','Ly','Lz'], $
                       t,x,y,z,dx,dy,dz,Lx,Ly,Lz)

  dx_1=zero  & dy_1=zero  & dz_1=zero
  dx_tilde=zero & dy_tilde=zero & dz_tilde=zero
endif else if (found_grid_der) then begin
  object = CREATE_STRUCT(name="pc_read_grid_" + $
                          str((size(x))[1]) + '_' + $
                          str((size(y))[1]) + '_' + $
                          str((size(z))[1]), $
                     ['t','x','y','z','dx','dy','dz', $
                      'dx_1','dy_1','dz_1', $
                      'dx_tilde','dy_tilde','dz_tilde'], $
                       t,x,y,z,dx,dy,dz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
endif


; If requested print a summary
fmt = '(A,4G15.6)'
if keyword_set(PRINT) then begin
  if n_elements(proc) eq 0 then begin
    print, FORMAT='(A,I2,A)', 'For all processors calculation domain:'
  endif else begin
    print, FORMAT='(A,I2,A)', 'For processor ',proc,' calculation domain:'
  endelse
  print, '             t = ', t
  print, 'min(x), max(x) = ',min(x),', ',max(x)
  print, 'min(y), max(y) = ',min(y),', ',max(y)
  print, 'min(z), max(z) = ',min(z),', ',max(z)
  print, '    dx, dy, dz = ' , dx , ', ' , dy , ', ' , dz
endif

end


