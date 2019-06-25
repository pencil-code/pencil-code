;+
;
; 14-Aug-2006/dintrans: coded
;
; NAME:
;       PC_READ_PENCILED
;
; PURPOSE:
;       Read a data/proc*/*.dat file which has been written using 
;       the 'ouput_pencil' subroutine (e.g. hcond.dat, glhc.dat, gg.dat,...)
;       Read both scalar and vector fields.
;
;       Works for one or all processors.
;
; CATEGORY:
;       Pencil Code, File I/O
;
; CALLING SEQUENCE:
;       PC_READ_PENCILED, datafile=datafile, field=field, vec=vec, $
;          TRIMALL=TRIMALL, dim=dim, datadir=datadir, HELP=HELP
;
; KEYWORD PARAMETERS:    
;  datafile: name for the data file [string]
;    datadir: specify the root data directory. Default is './data' [string]
;     field: the returned array; it can be a scalar (the default) 
;            or a vector [array]
;       vec: put vec=1 if you want to read a vector field (vec=0 by
;            default) [integer]
;   TRIMALL: remove ghost points from the field that is returned
;       dim: return the object 'dim' if not set up before [structure]
;      /HELP: display this usage information, and exit
;
;-
pro pc_read_penciled,datafile=datafile,field=field, vec=vec, $
TRIMALL=TRIMALL,dim=dim,HELP=HELP,datadir=datadir
;
IF (keyword_set(HELP)) THEN BEGIN
  doc_library,'pc_read_penciled'
  return
ENDIF
datadir = pc_get_datadir(datadir)
;
; Load HDF5 varfile if requested or available.
;
  if (strmid (datafile, strlen(datafile)-3) eq '.h5') then begin
    message, "pc_read_penciled: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    t = pc_read ('time', file=datafile, datadir=datadir+'/allprocs')
    quantities = h5_content ('/', number=num_quantities)
    found = where (strmid (quantities, 0, 4) eq 'data', num_quantities)
    if (num_quantities eq 0) then message, "pc_read_penciled: ERROR: no valid quantities found!"
    quantities = quantities[found]
help, quantities
print, quantities
    field = pc_read (quantities)
    h5_close_file
    return
  end
;
; - assume by default that we want to read a scalar
if (n_elements(vec) eq 0)   then vec=0
if (datafile eq 'glhc.dat') then vec=1
if (datafile eq 'gg.dat')   then vec=1
;
if (n_elements(dim) eq 0) then pc_read_dim,obj=dim,datadir=datadir
;
;  read local sizes
;
mxloc=0L & myloc=0L & mzloc=0L
close,1
openr,1,datadir+'/proc0/dim.dat'
readf,1,mxloc,myloc,mzloc
close,1
if (vec eq 0) then begin
  floc=fltarr(mxloc,myloc,mzloc)
  field=fltarr(dim.mx,dim.my,dim.mz)
endif else begin
  floc=fltarr(mxloc,myloc,mzloc,3)
  field=fltarr(dim.mx,dim.my,dim.mz,3)
endelse
;
nghostx=dim.nghostx & nghosty=dim.nghosty & nghostz=dim.nghostz
nxloc=mxloc-2*nghostx
nyloc=myloc-2*nghosty
nzloc=mzloc-2*nghostz
;
ncpus=dim.nprocx*dim.nprocy*dim.nprocz
dimfile='dim.dat'
; -- loop on all processors
for i=0,ncpus-1 do begin        ; read data from individual files
  tag='proc'+str(i)
  localdir=datadir+'/proc'+strtrim(i,2)
  ; read processor position
  dummy=''
  ipx=0L &ipy=0L &ipz=0L
  close,1
  openr,1,localdir+'/'+dimfile
  readf,1, dummy
  readf,1, dummy
  readf,1, dummy
  readf,1, ipx,ipy,ipz
  close,1
;
  if (ipx eq 0) then begin
    i0x=ipx*nxloc & i1x=i0x+mxloc-1
    i0xloc=0 & i1xloc=mxloc-1
  endif else begin
    i0x=ipx*nxloc+nghostx & i1x=i0x+mxloc-1-nghostx
    i0xloc=nghostx & i1xloc=mxloc-1
  endelse
;
  if (ipy eq 0) then begin
    i0y=ipy*nyloc & i1y=i0y+myloc-1
    i0yloc=0 & i1yloc=myloc-1
  endif else begin
    i0y=ipy*nyloc+nghosty & i1y=i0y+myloc-1-nghosty
    i0yloc=nghosty & i1yloc=myloc-1
  endelse
;
  if (ipz eq 0) then begin
    i0z=ipz*nzloc & i1z=i0z+mzloc-1
    i0zloc=0 & i1zloc=mzloc-1
  endif else begin
    i0z=ipz*nzloc+nghostz & i1z=i0z+mzloc-1-nghostz
    i0zloc=nghostz & i1zloc=mzloc-1
  endelse
; -- put the local data in the global array
  close,1
  print,'Reading '+localdir+'/'+datafile
  openr,1,localdir+'/'+datafile,/f77
  readu,1,floc
  close,1
  if (vec eq 0) then begin
    field[i0x:i1x,i0y:i1y,i0z:i1z]= $
      floc[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc]
  endif else begin
    field[i0x:i1x,i0y:i1y,i0z:i1z,*]= $
      floc[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]
  endelse
endfor
; -- trimall
if keyword_set(TRIMALL) then begin
  if (vec eq 0) then $
    field=reform(field[dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2]) else $
    field=reform(field[dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2,*])
endif
;
end
