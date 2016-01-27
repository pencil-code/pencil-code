;;
;; $Id$
;;
;;   Read grid.dat
;;
;;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;;  $Date: 2008-06-12 13:26:15 $
;;  $Revision: 1.17 $
;;
;;  27-nov-02/tony: coded
;;   2-oct-14/MR: keyword parameter down added for use with downsampled data
;;  27-jan-16/MR: added check for FORTRAN consistency og grid data + automatic endian swapping if necessary
;;
pro pc_read_grid, object=object, dim=dim, param=param, trimxyz=trimxyz, $
    datadir=datadir, proc=proc, print=print, quiet=quiet, help=help, $
    swap_endian=swap_endian, allprocs=allprocs, reduced=reduced, down=down
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
  common cdat_coords,coord_system
;
; Show some help!
;
  if (keyword_set(help)) then begin
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
    print, " allprocs: specify a collectively written grid file. Default is -1 [integer]  "
    print, "           if allprocs is equal to -1, allprocs will be chosen automatically. "
    print, " /reduced: use a reduced dataset grid file.                                   "
    print, "                                                                              "
    print, "   object: optional structure in which to return all the above     [structure]"
    print, "                                                                              "
    print, "   /PRINT: instruction to print all variables to standard output              "
    print, "   /QUIET: instruction not to print any 'helpful' information                 "
    print, "    /HELP: display this usage information, and exit                           "
    return
  ENDIF
;
; Default data directory
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
;
; Default filename
;
if (not keyword_set(down)) then $
  gridfile='grid.dat' $
else $
  gridfile='grid_down.dat'
;
; Default allprocs.
;
default, allprocs, -1
if (allprocs eq -1) then begin
  allprocs=0
  if (file_test(datadir+'/allprocs/'+gridfile) and (n_elements(proc) eq 0)) then allprocs=1
end
;
; Check if allprocs is consistent with proc.
;
if ((allprocs gt 0) and (n_elements(proc) ne 0)) then $
    message, "pc_read_grid: 'allprocs' and 'proc' cannot be set both."
;
; Check if reduced keyword is set.
;
if (keyword_set(reduced) and (n_elements(proc) ne 0)) then $
    message, "pc_read_grid: /reduced and 'proc' cannot be set both."
;
; Get necessary dimensions.
;
if (n_elements(dim) eq 0) then $
    pc_read_dim,object=dim,datadir=datadir,proc=proc,reduced=reduced,QUIET=QUIET, down=down
if (n_elements(param) eq 0) then $
    pc_read_param,object=param,datadir=datadir,QUIET=QUIET

if ((allprocs gt 0) or keyword_set (reduced)) then begin
  ncpus=1
endif else begin
  ncpus=dim.nprocx*dim.nprocy*dim.nprocz
endelse
;
; Set mx,my,mz in common block for derivative routines
;
nx=dim.nx
ny=dim.ny
nz=dim.nz
nw=nx*ny*nz
mx=dim.mx
my=dim.my
mz=dim.mz
mw=mx*my*mz
l1=dim.l1
l2=dim.l2
m1=dim.m1
m2=dim.m2
n1=dim.n1
n2=dim.n2
nghostx=dim.nghostx
nghosty=dim.nghosty
nghostz=dim.nghostz

lequidist=safe_get_tag(param,'lequidist',default=[1,1,1]) ne 0
lperi=param.lperi ne 0
ldegenerated=[nx,ny,nz] eq 1
;
;  Set coordinate system.
;
coord_system=param.coord_system
;
;  Check pc_precision is set!
;
pc_set_precision, dim=dim, quiet=quiet
;
; Initialize / set default returns for ALL variables
;
t=zero
x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
dx=zero & dy=zero & dz=zero
Lx=zero & Ly=zero & Lz=zero
dx_1=fltarr(mx)*one & dy_1=fltarr(my)*one & dz_1=fltarr(mz)*one
dx_tilde=fltarr(mx)*one & dy_tilde=fltarr(my)*one & dz_tilde=fltarr(mz)*one
;
; Get a unit number
;
get_lun, file

for i=0,ncpus-1 do begin
  ; Build the full path and filename
  if (keyword_set (reduced)) then begin
    filename=datadir+'/reduced/'+gridfile
  end else if (allprocs gt 0) then begin
    filename=datadir+'/allprocs/'+gridfile
  endif else if (n_elements(proc) ne 0) then begin
    filename=datadir+'/proc'+str(proc)+'/'+gridfile
  endif else begin
    filename=datadir+'/proc'+str(i)+'/'+gridfile
    ; Read processor box dimensions
    pc_read_dim,object=procdim,datadir=datadir,proc=i,QUIET=QUIET, down=down
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
  if (not file_test(filename)) then begin
    FREE_LUN,file
    message, 'ERROR: cannot find file '+ filename
  endif

  check=check_ftn_consistency(filename)
  if check eq -1 then begin
    print, 'File "'+trim(filename)+'" corrupted!'
    return
  endif else $
   if check eq 1 then begin
      print, 'Try to read file "'+trim(filename)+'" with swapped endian!'
      swap_endian=1-swap_endian
    endif
 
  if (not keyword_set(quiet)) then print, 'Reading ' , filename , '...'

  openr,file,filename,/F77,SWAP_ENDIAN=swap_endian

  if ((allprocs gt 0) or (n_elements(proc) ne 0) or keyword_set(reduced)) then begin
    readu,file, t,x,y,z
    readu,file, dx,dy,dz
    found_Lxyz=0
    found_grid_der=0
    on_ioerror, missing
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
    on_ioerror, missing
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
  on_ioerror, Null

  close,file

endfor
;
;
;
free_lun,file
;
;  Trim ghost zones of coordinate arrays.
;
if (keyword_set(trimxyz)) then begin
  x=x[l1:l2]
  y=y[m1:m2]
  z=z[n1:n2]
  dx_1=dx_1[l1:l2]
  dy_1=dy_1[m1:m2]
  dz_1=dz_1[n1:n2]
  dx_tilde=dx_tilde[l1:l2]
  dy_tilde=dy_tilde[m1:m2]
  dz_tilde=dz_tilde[n1:n2]
endif
;
;  Build structure of all the variables
;
if (found_Lxyz and found_grid_der) then begin
  if (keyword_set(trimxyz)) then begin
    Ox = x[0] - lperi[0] * 0.5 / dx_1[0]
    Oy = y[0] - lperi[1] * 0.5 / dy_1[0]
    Oz = z[0] - lperi[2] * 0.5 / dz_1[0]
  endif else begin
    Ox = x[nghostx] - lperi[0] * 0.5 / dx_1[nghostx]
    Oy = y[nghosty] - lperi[1] * 0.5 / dy_1[nghosty]
    Oz = z[nghostz] - lperi[2] * 0.5 / dz_1[nghostz]
  endelse
  object = create_struct(name="pc_read_grid_" + $
      str((size(x))[1]) + '_' + $
      str((size(y))[1]) + '_' + $
      str((size(z))[1]), $
      ['t','x','y','z','dx','dy','dz','Ox','Oy','Oz','Lx','Ly','Lz', $
       'dx_1','dy_1','dz_1','dx_tilde','dy_tilde','dz_tilde', $
       'lequidist','lperi','ldegenerated'], $
      t,x,y,z,dx,dy,dz,Ox,Oy,Oz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde, $
        lequidist,lperi,ldegenerated)
endif else if (found_Lxyz) then begin
  object = create_struct(name="pc_read_grid_" + $
      str((size(x))[1]) + '_' + $
      str((size(y))[1]) + '_' + $
      str((size(z))[1]), $
      ['t','x','y','z','dx','dy','dz','Lx','Ly','Lz'], $
      t,x,y,z,dx,dy,dz,Lx,Ly,Lz)
  dx_1=zero  & dy_1=zero  & dz_1=zero
  dx_tilde=zero & dy_tilde=zero & dz_tilde=zero
endif else if (found_grid_der) then begin
  object = create_struct(name="pc_read_grid_" + $
      str((size(x))[1]) + '_' + $
      str((size(y))[1]) + '_' + $
      str((size(z))[1]), $
      ['t','x','y','z','dx','dy','dz','dx_1','dy_1','dz_1', $
       'dx_tilde','dy_tilde','dz_tilde'], $
      t,x,y,z,dx,dy,dz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
endif
;
; If requested print a summary
;
fmt = '(A,4G15.6)'
if (keyword_set(print)) then begin
  if (n_elements(proc) eq 0) then begin
    print, FORMAT='(A,I2,A)', 'For all processors calculation domain:'
  endif else if (keyword_set(reduced)) then begin
    print, FORMAT='(A,I2,A)', 'For reduced calculation domain:'
  endif else begin
    print, FORMAT='(A,I2,A)', 'For processor ',proc,' calculation domain:'
  endelse
  print, '             t = ', t
  print, 'min(x), max(x) = ',min(x),', ',max(x)
  print, 'min(y), max(y) = ',min(y),', ',max(y)
  print, 'min(z), max(z) = ',min(z),', ',max(z)
  print, '    dx, dy, dz = ' , dx , ', ' , dy , ', ' , dz
  print, '   periodicity = ' , lperi
  print, '  degeneration = ' , ldegenerated
endif
;
end
