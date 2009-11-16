;
; $Id: pc_read_pdim.pro 9839 2008-09-05 07:24:02Z ajohan $
;
;  Read block domain decomposition from file.
;
pro pc_read_block, object=object, datadir=datadir, proc=proc, $
    print=print, quiet=quiet
compile_opt IDL2,HIDDEN
COMMON pc_precision, zero, one
;
; Default processor.
;
default, proc, 0
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
pc_read_bdim, obj=bdim
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
filename=datadir+'/proc'+strtrim(proc,2)+'/block.dat'
dummy=findfile(filename, count=found)
if (found gt 0) then begin
  if ( not keyword_set(quiet) ) then print, 'Reading ' + filename + '...'
  get_lun, file
  openr, file, filename, /f77
  readu, file, t
  readu, file, nblock_loc, nproc_parent, nproc_foster
  readu, file, iproc_foster_brick
  if (nblock_loc gt 0) then begin
    iproc_parent_block=lonarr(nblock_loc)
    ibrick_parent_block=lonarr(nblock_loc)
    xb=fltarr(mxb,nblock_loc)*one
    yb=fltarr(myb,nblock_loc)*one
    zb=fltarr(mzb,nblock_loc)*one
    readu, file, iproc_parent_block
    readu, file, ibrick_parent_block
    readu, file, xb
    readu, file, yb
    readu, file, zb
  endif else begin
    dummy=0L
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
    readu, file, dummy
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
     'iproc_parent_list','iproc_foster_list'], $
     t,nblock_loc,nproc_parent,nproc_foster,iproc_foster_brick, $
     iproc_parent_block,ibrick_parent_block,xb,yb,zb, $
     iproc_parent_list,iproc_foster_list)
;
end
