;
;  $Id$
;
;  Same as pc_varcontent, but for global variables.
;
FUNCTION pc_varcontent_global, datadir=datadir, dim=dim, $
    run2D=run2D, param=param, quiet=quiet, scalar=scalar
COMPILE_OPT IDL2,HIDDEN
;
;  Read grid dimensions, input parameters and location of datadir.
;
if (n_elements(dim) eq 0) then pc_read_dim,obj=dim,datadir=datadir,quiet=quiet
if (n_elements(param) eq 0) then pc_read_param,obj=param,datadir=datadir, $
    dim=dim,quiet=quiet
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
; 
;  Read the positions of variables in f.
;
cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),     "\n"'+"' "+datadir+'/index.pro'
spawn, cmd, result
res = flatten_strings(result)
index=strsplit(res,'&',/extract)
nindex=n_elements(index)
for i=0,nindex-1 do begin
  if (execute(index[i]) ne 1) then $
  message, 'pc_varcontent_global: there was a problem with index.pro', /info
endfor
;
;  Make an array of structures in which to store their descriptions.
;  Index zero is kept as a dummy entry.
;
varcontent=replicate({varcontent_all, $
    variable:'UNKNOWN', $
    idlvar:'dummy', $
    idlinit:'fltarr(mx,my,mz)*one', $
    idlvarloc:'dummy_loc', $
    idlinitloc:'fltarr(mxloc,myloc,mzloc)*one', $
    skip:0}, dim.mglobal+1)
;
;  Predefine some variable types used regularly
;
if (keyword_set(run2D)) then begin
;
;  For 2-D runs with lwrite_2d=T. Data has been written by the code without
;  ghost zones in the missing direction. We add ghost zones here anyway so
;  that the array can be treated exactly like 3-D data.
;  
  if (dim.nx eq 1) then begin
;  2-D run in (y,z) plane.
    INIT_3VECTOR     = 'fltarr(mx,my,mz,3)*one'
    INIT_3VECTOR_LOC = 'fltarr(myloc,mzloc,3)*one'
    INIT_SCALAR      = 'fltarr(mx,my,mz)*one'
    INIT_SCALAR_LOC  = 'fltarr(myloc,mzloc)*one'
  endif else if (dim.ny eq 1) then begin
;  2-D run in (x,z) plane.
    INIT_3VECTOR     = 'fltarr(mx,my,mz,3)*one'
    INIT_3VECTOR_LOC = 'fltarr(mxloc,mzloc,3)*one'
    INIT_SCALAR      = 'fltarr(mx,my,mz)*one'
    INIT_SCALAR_LOC  = 'fltarr(mxloc,mzloc)*one'
  endif else begin
;  2-D run in (x,y) plane.
    INIT_3VECTOR     = 'fltarr(mx,my,mz,3)*one'
    INIT_3VECTOR_LOC = 'fltarr(mxloc,myloc,3)*one'
    INIT_SCALAR      = 'fltarr(mx,my,mz)*one'
    INIT_SCALAR_LOC  = 'fltarr(mxloc,myloc)*one'
  endelse
endif else begin
;
;  Regular 3-D run.
;
  INIT_3VECTOR     = 'fltarr(mx,my,mz,3)*one'
  INIT_3VECTOR_LOC = 'fltarr(mxloc,myloc,mzloc,3)*one'
  INIT_SCALAR      = 'fltarr(mx,my,mz)*one'
  INIT_SCALAR_LOC  = 'fltarr(mxloc,myloc,mzloc)*one'
endelse
;
; For EVERY POSSIBLE variable in a var file, store a
; description of the variable in an indexed array of structures
; where the indexes line up with those in the saved f array
;
; Any variable not stored should have iXXXXXX set to zero
; and will only update the dummy index zero entry
;
default, iglobal_bx_ext, 0
if (iglobal_bx_ext ne 0) then iglobal_bx_ext=iglobal_bx_ext-dim.mvar-dim.maux
varcontent[iglobal_bx_ext].variable   = 'Magnetic field (bx_ext)'
varcontent[iglobal_bx_ext].idlvar     = 'bx_ext'
varcontent[iglobal_bx_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_bx_ext].idlvarloc  = 'bx_ext_loc'
varcontent[iglobal_bx_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_by_ext, 0
if (iglobal_by_ext ne 0) then iglobal_by_ext=iglobal_by_ext-dim.mvar-dim.maux
varcontent[iglobal_by_ext].variable   = 'Magnetic field (by_ext)'
varcontent[iglobal_by_ext].idlvar     = 'by_ext'
varcontent[iglobal_by_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_by_ext].idlvarloc  = 'by_ext_loc'
varcontent[iglobal_by_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_bz_ext, 0
if (iglobal_bz_ext ne 0) then iglobal_bz_ext=iglobal_bz_ext-dim.mvar-dim.maux
varcontent[iglobal_bz_ext].variable   = 'Magnetic field (bz_ext)'
varcontent[iglobal_bz_ext].idlvar     = 'bz_ext'
varcontent[iglobal_bz_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_bz_ext].idlvarloc  = 'bz_ext_loc'
varcontent[iglobal_bz_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_jx_ext, 0
if (iglobal_jx_ext ne 0) then iglobal_jx_ext=iglobal_jx_ext-dim.mvar-dim.maux
varcontent[iglobal_jx_ext].variable   = 'Current density (jx_ext)'
varcontent[iglobal_jx_ext].idlvar     = 'jx_ext'
varcontent[iglobal_jx_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_jx_ext].idlvarloc  = 'jx_ext_loc'
varcontent[iglobal_jx_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_jy_ext, 0
if (iglobal_jy_ext ne 0) then iglobal_jy_ext=iglobal_jy_ext-dim.mvar-dim.maux
varcontent[iglobal_jy_ext].variable   = 'Current density (jy_ext)'
varcontent[iglobal_jy_ext].idlvar     = 'jy_ext'
varcontent[iglobal_jy_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_jy_ext].idlvarloc  = 'jy_ext_loc'
varcontent[iglobal_jy_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_jz_ext, 0
if (iglobal_jz_ext ne 0) then iglobal_jz_ext=iglobal_jz_ext-dim.mvar-dim.maux
varcontent[iglobal_jz_ext].variable   = 'Current density (jz_ext)'
varcontent[iglobal_jz_ext].idlvar     = 'jz_ext'
varcontent[iglobal_jz_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_jz_ext].idlvarloc  = 'jz_ext_loc'
varcontent[iglobal_jz_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_ex_ext, 0
if (iglobal_ex_ext ne 0) then iglobal_ex_ext=iglobal_ex_ext-dim.mvar-dim.maux
varcontent[iglobal_ex_ext].variable   = 'Electromotive force (ex_ext)'
varcontent[iglobal_ex_ext].idlvar     = 'ex_ext'
varcontent[iglobal_ex_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_ex_ext].idlvarloc  = 'ex_ext_loc'
varcontent[iglobal_ex_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_ey_ext, 0
if (iglobal_ey_ext ne 0) then iglobal_ey_ext=iglobal_ey_ext-dim.mvar-dim.maux
varcontent[iglobal_ey_ext].variable   = 'Electromotive force (ey_ext)'
varcontent[iglobal_ey_ext].idlvar     = 'ey_ext'
varcontent[iglobal_ey_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_ey_ext].idlvarloc  = 'ey_ext_loc'
varcontent[iglobal_ey_ext].idlinitloc = INIT_SCALAR_LOC
;
default, iglobal_ez_ext, 0
if (iglobal_ez_ext ne 0) then iglobal_ez_ext=iglobal_ez_ext-dim.mvar-dim.maux
varcontent[iglobal_ez_ext].variable   = 'Electromotive force (ez_ext)'
varcontent[iglobal_ez_ext].idlvar     = 'ez_ext'
varcontent[iglobal_ez_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_ez_ext].idlvarloc  = 'ez_ext_loc'
varcontent[iglobal_ez_ext].idlinitloc = INIT_SCALAR_LOC
;
default, ics2, 0
if (ics2 ne 0) then ics2=ics2-dim.mvar-dim.maux
varcontent[ics2].variable   = 'Sound speed'
varcontent[ics2].idlvar     = 'cs2'
varcontent[ics2].idlinit    = INIT_SCALAR
varcontent[ics2].idlvarloc  = 'cs2_loc'
varcontent[ics2].idlinitloc = INIT_SCALAR_LOC
;
default, iglnTT, 0
if (iglnTT ne 0) then iglnTT=iglnTT-dim.mvar-dim.maux
varcontent[iglnTT].variable   = 'Gradient of logarithmic temperature'
varcontent[iglnTT].idlvar     = 'glnTT'
varcontent[iglnTT].idlinit    = INIT_3VECTOR
varcontent[iglnTT].idlvarloc  = 'lnTT_loc'
varcontent[iglnTT].idlinitloc = INIT_3VECTOR_LOC
varcontent[iglnTT].skip       = 2
;
default, igg, 0
if (igg ne 0) then igg=igg-dim.mvar-dim.maux
varcontent[igg].variable   = 'Gravitational acceleration'
varcontent[igg].idlvar     = 'gg'
varcontent[igg].idlinit    = INIT_3VECTOR
varcontent[igg].idlvarloc  = 'gg_loc'
varcontent[igg].idlinitloc = INIT_3VECTOR_LOC
varcontent[igg].skip       = 2
;
default, iglobal_gg, 0
if (iglobal_gg ne 0) then iglobal_gg=iglobal_gg-dim.mvar-dim.maux
varcontent[iglobal_gg].variable   = 'Gravitational acceleration'
varcontent[iglobal_gg].idlvar     = 'gg'
varcontent[iglobal_gg].idlinit    = INIT_3VECTOR
varcontent[iglobal_gg].idlvarloc  = 'gg_loc'
varcontent[iglobal_gg].idlinitloc = INIT_3VECTOR_LOC
varcontent[iglobal_gg].skip       = 2
;
;  Zero out default definition in case it has been set by mistake.
;
varcontent[0].variable = 'UNKNOWN'
varcontent[0].idlvar   = 'UNKNOWN'
varcontent[0].idlinit  = '0.0'
varcontent[0].skip     = 0
;
return, varcontent
;
END
