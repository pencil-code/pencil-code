;
; $Id: pc_read_pdim.pro 9839 2008-09-05 07:24:02Z ajohan $
;+
;  Read block domain decomposition from file.
;-
pro pc_read_block, object=object, datadir=datadir, varfile=varfile, $
    proc=proc, swap=swap, quiet=quiet, single=single, help=help
compile_opt IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
  if (keyword_set(help)) then begin
    doc_library, 'pc_read_block'
    return
  endif
;
; Default settings.
;
default, varfile, 'blocks.dat'
default, proc, 0
default, swap, 0
default, quiet, 0
default, single, 0
;
; Default data directory.
;
datadir = pc_get_datadir(datadir)
;
; Read dimension.
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
;
; Read block dimensions from file.
;
pc_read_bdim, obj=bdim, /quiet
nxb=bdim.nxb
nyb=bdim.nyb
nzb=bdim.nzb
mxb=bdim.mxb
myb=bdim.myb
mzb=bdim.mzb
nbx=bdim.nbrickx/dim.nprocx
nby=bdim.nbricky/dim.nprocy
nbz=bdim.nbrickz/dim.nprocz
nbricks=nbx*nby*nbz
nblockmax=bdim.nblockmax
;
; Initialize all variables.
;
t=zero
nblock_loc=0L
nproc_parent=0L
nproc_foster=0L
iproc_foster_brick=lonarr(nbricks)
;
; Check for existence and read the data.
;
filename=datadir+'/proc'+strtrim(proc,2)+'/'+varfile
if (file_test(filename)) then begin
  if ( not keyword_set(quiet) ) then print, 'Reading ' + filename + '...'
  get_lun, file
  openr, file, filename, /f77, swap_endian=swap
  readu, file, t
  readu, file, nblock_loc, nproc_parent, nproc_foster
  readu, file, iproc_foster_brick
  if (nblock_loc gt 0) then begin
    iproc_parent_block=lonarr(nblock_loc)
    ibrick_parent_block=lonarr(nblock_loc)
    xb=make_array(mxb,nblock_loc, type=type_idl)
    yb=make_array(myb,nblock_loc, type=type_idl)
    zb=make_array(mzb,nblock_loc, type=type_idl)
    dx1b=make_array(mxb,nblock_loc, type=type_idl)
    dy1b=make_array(myb,nblock_loc, type=type_idl)
    dz1b=make_array(mzb,nblock_loc, type=type_idl)
    dVol1xb=make_array(mxb,nblock_loc, type=type_idl)
    dVol1yb=make_array(myb,nblock_loc, type=type_idl)
    dVol1zb=make_array(mzb,nblock_loc, type=type_idl)
    readu, file, iproc_parent_block
    readu, file, ibrick_parent_block
    readu, file, xb
    readu, file, yb
    readu, file, zb
    readu, file, dx1b
    readu, file, dy1b
    readu, file, dz1b
    readu, file, dVol1xb
    readu, file, dVol1yb
    readu, file, dVol1zb
    if (single and precision eq 'D') then begin
      xb=float(xb)
      yb=float(yb)
      zb=float(zb)
      dx1b=float(dx1b)
      dy1b=float(dy1b)
      dz1b=float(dz1b)
      dVol1xb=float(dVol1xb)
      dVol1yb=float(dVol1yb)
      dVol1zb=float(dVol1zb)
    endif
  endif else begin
    dummy=0L
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    iproc_parent_block=-1
    ibrick_parent_block=-1
    xb=single ? -1. : -one
    yb=xb & zb=xb & dx1b=xb & dy1b=xb & dz1b=xb
    dVol1xb=xb & dVol1yb=xb & dVol1zb=xb
  endelse
  if (nproc_parent gt 0) then begin
    iproc_parent_list=lonarr(nproc_parent)
    readu, file, iproc_parent_list
  endif else begin
    iproc_parent_list=-1
    dummy=0L
    readu, file, dummy
  endelse
  if (nproc_foster gt 0) then begin
    iproc_foster_list=lonarr(nproc_foster)
    readu, file, iproc_foster_list
  endif else begin
    iproc_foster_list=-1
    dummy=0L
    readu, file, dummy
  endelse
  free_lun,file
endif else begin
  message, 'ERROR: cannot find file ' + filename
endelse
if single then t=float(t)
;
; Build structure of all the variables
;
if nblock_loc gt 0 then $
  object = create_struct(name=objectname, $
      ['t','nblock_loc','nproc_parent','nproc_foster','iproc_foster_brick', $
       'iproc_parent_block','ibrick_parent_block','xb','yb','zb', $
       'dx1b','dy1b','dz1b','dVol1xb','dVol1yb','dVol1zb', $
       'iproc_parent_list','iproc_foster_list'], $
       t,nblock_loc,nproc_parent,nproc_foster,iproc_foster_brick, $
       iproc_parent_block,ibrick_parent_block,xb,yb,zb, $
       dx1b,dy1b,dz1b,dVol1xb,dVol1yb,dVol1zb, $
       iproc_parent_list,iproc_foster_list)
;
end
