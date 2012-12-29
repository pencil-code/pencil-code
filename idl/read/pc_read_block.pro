;
; $Id: pc_read_pdim.pro 9839 2008-09-05 07:24:02Z ajohan $
;
;  Read block domain decomposition from file.
;
pro pc_read_block, object=object, datadir=datadir, proc=proc, $
    swap=swap, quiet=quiet
compile_opt IDL2,HIDDEN
COMMON pc_precision, zero, one
;
; Default settings.
;
default, proc, 0
default, swap, 0
default, quiet, 0
;
; Default data directory.
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
;
; Set precision.
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
pc_set_precision, dim=dim, /quiet
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
t=0.0d
nblock_loc=0L
nproc_parent=0L
nproc_foster=0L
iproc_foster_brick=lonarr(nbricks)
;
; Check for existence and read the data.
;
filename=datadir+'/proc'+strtrim(proc,2)+'/blocks.dat'
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
    xb=fltarr(mxb,nblock_loc)*one
    yb=fltarr(myb,nblock_loc)*one
    zb=fltarr(mzb,nblock_loc)*one
    dx1b=fltarr(mxb,nblock_loc)*one
    dy1b=fltarr(myb,nblock_loc)*one
    dz1b=fltarr(mzb,nblock_loc)*one
    dVol1xb=fltarr(mxb,nblock_loc)*one
    dVol1yb=fltarr(myb,nblock_loc)*one
    dVol1zb=fltarr(mzb,nblock_loc)*one
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
  endif else begin
    dummy=0L
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    iproc_parent_block=-1
    ibrick_parent_block=-1
    xb=-1.0*one
    yb=-1.0*one
    zb=-1.0*one
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
  free_lun,file
  message, 'ERROR: cannot find file ' + filename
endelse
;
; Build structure of all the variables
;
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
