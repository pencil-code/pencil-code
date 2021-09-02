;
;  $Id$
;+
;  Same as pc_varcontent, but for global variables.
;-
FUNCTION pc_varcontent_global, datadir=datadir, dim=dim, help=help, $
    run2D=run2D, param=param, quiet=quiet, scalar=scalar, single=single
COMPILE_OPT IDL2,HIDDEN
;
  if (keyword_set(help)) then begin
    doc_library, 'pc_varcontent_global'
    return
  endif

default, single, 0
;
;  Read grid dimensions, input parameters and location of datadir.
;
if (n_elements(dim) eq 0) then pc_read_dim,obj=dim,datadir=datadir,quiet=quiet
if (n_elements(param) eq 0) then pc_read_param,obj=param,datadir=datadir, $
    dim=dim,quiet=quiet, single=single
datadir = pc_get_datadir(datadir)
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
;  Predefine some variable types used regularly
;

INIT_3VECTOR = single ? 'fltarr(mx,my,mz,3)' : 'make_array(mx,my,mz,3, type=type_idl)'
INIT_SCALAR  = single ? 'fltarr(mx,my,mz)' : 'make_array(mx,my,mz, type=type_idl)'
if (keyword_set(run2D)) then begin
;
;  For 2-D runs with lwrite_2d=T. Data has been written by the code without
;  ghost zones in the missing direction. We add ghost zones here anyway so
;  that the array can be treated exactly like 3-D data.
;  
  if (dim.nx eq 1) then begin
;  2-D run in (y,z) plane.
    INIT_3VECTOR_LOC = 'make_array(myloc,mzloc,3, type=type_idl)'
    INIT_SCALAR_LOC  = 'make_array(myloc,mzloc, type=type_idl)'
  endif else if (dim.ny eq 1) then begin
;  2-D run in (x,z) plane.
    INIT_3VECTOR_LOC = 'make_array(mxloc,mzloc,3, type=type_idl)'
    INIT_SCALAR_LOC  = 'make_array(mxloc,mzloc, type=type_idl)'
  endif else begin
;  2-D run in (x,y) plane.
    INIT_3VECTOR_LOC = 'make_array(mxloc,myloc,3, type=type_idl)'
    INIT_SCALAR_LOC  = 'make_array(mxloc,myloc, type=type_idl)'
  endelse
endif else begin
;
;  Regular 3-D run.
;
  INIT_3VECTOR_LOC = 'make_array(mxloc,myloc,mzloc,3, type=type_idl)'
  INIT_SCALAR_LOC  = 'make_array(mxloc,myloc,mzloc, type=type_idl)'
endelse
;
;  Make an array of structures in which to store their descriptions.
;  Index zero is kept as a dummy entry.
;
varcontent=replicate({varcontent_all, $
    variable:'UNKNOWN', $
    idlvar:'dummy', $
    idlinit:INIT_SCALAR, $
    idlvarloc:'dummy_loc', $
    idlinitloc:'make_array(mxloc,myloc,mzloc, type=type_idl)', $
    skip:0}, dim.mglobal+1)
;
; For EVERY POSSIBLE variable in a var file, store a
; description of the variable in an indexed array of structures
; where the indexes line up with those in the saved f array
;
; Any variable not stored should have iXXXXXX set to zero
; and will only update the dummy index zero entry
;
default, iglobal_bx_ext, 0
if (iglobal_bx_ext gt 0) then begin
  iglobal_bx_ext=iglobal_bx_ext-dim.mvar-dim.maux
  varcontent[iglobal_bx_ext].variable   = 'Magnetic field (bx_ext)'
  varcontent[iglobal_bx_ext].idlvar     = 'bx_ext'
  varcontent[iglobal_bx_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_bx_ext].idlvarloc  = 'bx_ext_loc'
  varcontent[iglobal_bx_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_by_ext, 0
if (iglobal_by_ext gt 0) then begin
  iglobal_by_ext=iglobal_by_ext-dim.mvar-dim.maux
  varcontent[iglobal_by_ext].variable   = 'Magnetic field (by_ext)'
  varcontent[iglobal_by_ext].idlvar     = 'by_ext'
  varcontent[iglobal_by_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_by_ext].idlvarloc  = 'by_ext_loc'
  varcontent[iglobal_by_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_bz_ext, 0
if (iglobal_bz_ext gt 0) then begin
  iglobal_bz_ext=iglobal_bz_ext-dim.mvar-dim.maux
  varcontent[iglobal_bz_ext].variable   = 'Magnetic field (bz_ext)'
  varcontent[iglobal_bz_ext].idlvar     = 'bz_ext'
  varcontent[iglobal_bz_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_bz_ext].idlvarloc  = 'bz_ext_loc'
  varcontent[iglobal_bz_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_ax_ext, 0
if (iglobal_ax_ext gt 0) then begin
  iglobal_ax_ext=iglobal_ax_ext-dim.mvar-dim.maux
  varcontent[iglobal_ax_ext].variable   = 'Magnetic vector potential (ax_ext)'
  varcontent[iglobal_ax_ext].idlvar     = 'ax_ext'
  varcontent[iglobal_ax_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_ax_ext].idlvarloc  = 'ax_ext_loc'
  varcontent[iglobal_ax_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_ay_ext, 0
if (iglobal_ay_ext gt 0) then begin
  iglobal_ay_ext=iglobal_ay_ext-dim.mvar-dim.maux
  varcontent[iglobal_ay_ext].variable   = 'Magnetic vector potential (ay_ext)'
  varcontent[iglobal_ay_ext].idlvar     = 'ay_ext'
  varcontent[iglobal_ay_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_ay_ext].idlvarloc  = 'ay_ext_loc'
  varcontent[iglobal_ay_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_az_ext, 0
if (iglobal_az_ext gt 0) then begin
  iglobal_az_ext=iglobal_az_ext-dim.mvar-dim.maux
  varcontent[iglobal_az_ext].variable   = 'Magnetic vector potential (az_ext)'
  varcontent[iglobal_az_ext].idlvar     = 'az_ext'
  varcontent[iglobal_az_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_az_ext].idlvarloc  = 'az_ext_loc'
  varcontent[iglobal_az_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_jx_ext, 0
if (iglobal_jx_ext gt 0) then begin
  iglobal_jx_ext=iglobal_jx_ext-dim.mvar-dim.maux
  varcontent[iglobal_jx_ext].variable   = 'Current density (jx_ext)'
  varcontent[iglobal_jx_ext].idlvar     = 'jx_ext'
  varcontent[iglobal_jx_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_jx_ext].idlvarloc  = 'jx_ext_loc'
  varcontent[iglobal_jx_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_jy_ext, 0
if (iglobal_jy_ext gt 0) then begin
  iglobal_jy_ext=iglobal_jy_ext-dim.mvar-dim.maux
  varcontent[iglobal_jy_ext].variable   = 'Current density (jy_ext)'
  varcontent[iglobal_jy_ext].idlvar     = 'jy_ext'
  varcontent[iglobal_jy_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_jy_ext].idlvarloc  = 'jy_ext_loc'
  varcontent[iglobal_jy_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_jz_ext, 0
if (iglobal_jz_ext gt 0) then begin
  iglobal_jz_ext=iglobal_jz_ext-dim.mvar-dim.maux
  varcontent[iglobal_jz_ext].variable   = 'Current density (jz_ext)'
  varcontent[iglobal_jz_ext].idlvar     = 'jz_ext'
  varcontent[iglobal_jz_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_jz_ext].idlvarloc  = 'jz_ext_loc'
  varcontent[iglobal_jz_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_ex_ext, 0
if (iglobal_ex_ext gt 0) then begin
  iglobal_ex_ext=iglobal_ex_ext-dim.mvar-dim.maux
  varcontent[iglobal_ex_ext].variable   = 'Electromotive force (ex_ext)'
  varcontent[iglobal_ex_ext].idlvar     = 'ex_ext'
  varcontent[iglobal_ex_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_ex_ext].idlvarloc  = 'ex_ext_loc'
  varcontent[iglobal_ex_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_ey_ext, 0
if (iglobal_ey_ext gt 0) then begin
  iglobal_ey_ext=iglobal_ey_ext-dim.mvar-dim.maux
  varcontent[iglobal_ey_ext].variable   = 'Electromotive force (ey_ext)'
  varcontent[iglobal_ey_ext].idlvar     = 'ey_ext'
  varcontent[iglobal_ey_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_ey_ext].idlvarloc  = 'ey_ext_loc'
  varcontent[iglobal_ey_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglobal_ez_ext, 0
if (iglobal_ez_ext gt 0) then begin
  iglobal_ez_ext=iglobal_ez_ext-dim.mvar-dim.maux
  varcontent[iglobal_ez_ext].variable   = 'Electromotive force (ez_ext)'
  varcontent[iglobal_ez_ext].idlvar     = 'ez_ext'
  varcontent[iglobal_ez_ext].idlinit    = INIT_SCALAR
  varcontent[iglobal_ez_ext].idlvarloc  = 'ez_ext_loc'
  varcontent[iglobal_ez_ext].idlinitloc = INIT_SCALAR_LOC
endif
;
default, ics2, 0
if (ics2 gt 0) then begin
  ics2=ics2-dim.mvar-dim.maux
  varcontent[ics2].variable   = 'Sound speed'
  varcontent[ics2].idlvar     = 'cs2'
  varcontent[ics2].idlinit    = INIT_SCALAR
  varcontent[ics2].idlvarloc  = 'cs2_loc'
  varcontent[ics2].idlinitloc = INIT_SCALAR_LOC
endif
;
default, iglnTT, 0
default, iglnTx, 0
default, iglnTy, 0
default, iglnTz, 0
if ((iglnTT gt 0) or ((iglnTx gt 0) and (iglnTy gt 0) and (iglnTz gt 0))) then begin
  iglnTT=max ([iglnTT,iglnTx])-dim.mvar-dim.maux
  varcontent[iglnTT].variable   = 'Gradient of logarithmic temperature'
  varcontent[iglnTT].idlvar     = 'glnTT'
  varcontent[iglnTT].idlinit    = INIT_3VECTOR
  varcontent[iglnTT].idlvarloc  = 'lnTT_loc'
  varcontent[iglnTT].idlinitloc = INIT_3VECTOR_LOC
  varcontent[iglnTT].skip       = 2
endif
;
default, igg, 0
default, iglobal_ggx, 0
default, iglobal_ggy, 0
default, iglobal_ggz, 0
if ((igg gt 0) or ((iglobal_ggx gt 0) and (iglobal_ggy gt 0) and (iglobal_ggz gt 0))) then begin
  igg=max ([igg,iglobal_ggx])-dim.mvar-dim.maux
  varcontent[igg].variable   = 'Gravitational acceleration'
  varcontent[igg].idlvar     = 'gg'
  varcontent[igg].idlinit    = INIT_3VECTOR
  varcontent[igg].idlvarloc  = 'gg_loc'
  varcontent[igg].idlinitloc = INIT_3VECTOR_LOC
  varcontent[igg].skip       = 2
endif
;
default, iglobal_gg, 0
if (iglobal_gg gt 0) then begin
  iglobal_gg=iglobal_gg-dim.mvar-dim.maux
  varcontent[iglobal_gg].variable   = 'Gravitational acceleration'
  varcontent[iglobal_gg].idlvar     = 'gg'
  varcontent[iglobal_gg].idlinit    = INIT_3VECTOR
  varcontent[iglobal_gg].idlvarloc  = 'gg_loc'
  varcontent[iglobal_gg].idlinitloc = INIT_3VECTOR_LOC
  varcontent[iglobal_gg].skip       = 2
endif
;
;  Remove empty entry at position 0:
;
if (n_elements(varcontent) gt 1) then varcontent=varcontent[1:*] else varcontent=0 
;
return, varcontent
;
END
