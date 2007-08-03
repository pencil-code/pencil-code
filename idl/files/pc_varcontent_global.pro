;  $Id: pc_varcontent_global.pro,v 1.2 2007-08-03 09:53:27 ajohan Exp $
FUNCTION pc_varcontent_global,datadir=datadir,dim=dim, $
                       param=param,quiet=quiet,scalar=scalar,run2D=run2D
COMPILE_OPT IDL2,HIDDEN
; 
;  Read the positions of global variables in f
;  Can't just use `@data/index', as the data directory may have a different name
;
if (n_elements(dim) eq 0) then pc_read_dim,obj=dim,datadir=datadir,quiet=quiet
if (n_elements(param) eq 0) then pc_read_param,obj=param,datadir=datadir, $
    dim=dim,quiet=quiet

if (not keyword_set(datadir)) then datadir=pc_get_datadir()
cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),"\n"'+"' "+datadir+'/index.pro'
spawn, cmd, result
res = flatten_strings(result) 

if (execute(res) ne 1) then $
    message, 'There was a problem with index.pro', /INFO

; Make an array of structures in which to store their descriptions
; index zero is kept as a dummy entry.
varcontent=REPLICATE({varcontent_all, variable:'UNKNOWN', $ 
                                      idlvar:'dummy', $
                                      idlinit:'fltarr(mx,my,mz)*one', $
                                      idlvarloc:'dummy_loc', $
                                      idlinitloc:'fltarr(mxloc,myloc,mzloc)*one', $
                                      skip:0},dim.mglobal+1)
;Predefine some variable types used regularly
if (not keyword_set(run2D)) then begin
  ; classical 3D-run (x,y,z)
  INIT_3VECTOR     = 'fltarr(mx,my,mz,3)*one'
  INIT_3VECTOR_LOC = 'fltarr(mxloc,myloc,mzloc,3)*one'
  INIT_SCALAR      = 'fltarr(mx,my,mz)*one'
  INIT_SCALAR_LOC  = 'fltarr(mxloc,myloc,mzloc)*one'
endif else begin
  if (dim.ny eq 1) then begin
    ; 2D_run in plane (x,z)
    INIT_3VECTOR     = 'fltarr(mx,mz,3)*one'
    INIT_3VECTOR_LOC = 'fltarr(mxloc,mzloc,3)*one'
    INIT_SCALAR      = 'fltarr(mx,mz)*one'
    INIT_SCALAR_LOC  = 'fltarr(mxloc,mzloc)*one'
  endif else begin
    ; 2D_run in plane (x,y)
    INIT_3VECTOR     = 'fltarr(mx,my,3)*one'
    INIT_3VECTOR_LOC = 'fltarr(mxloc,myloc,3)*one'
    INIT_SCALAR      = 'fltarr(mx,my)*one'
    INIT_SCALAR_LOC  = 'fltarr(mxloc,myloc)*one'
  endelse
endelse

default, iglobal_bx_ext, 0
default, iglobal_by_ext, 0
default, iglobal_bz_ext, 0
default, iglobal_jx_ext, 0
default, iglobal_jy_ext, 0
default, iglobal_jz_ext, 0
default, iglobal_ex_ext, 0
default, iglobal_ey_ext, 0
default, iglobal_ez_ext, 0

; For EVERY POSSIBLE variable in a var file, store a
; description of the variable in an indexed array of structures
; where the indexes line up with those in the saved f array

; Any variable not stored should have iXXXXXX set to zero
; and will only update the dummy index zero entry

if (iglobal_bx_ext ne 0) then iglobal_bx_ext=iglobal_bx_ext-dim.mvar-dim.maux
if (iglobal_by_ext ne 0) then iglobal_by_ext=iglobal_by_ext-dim.mvar-dim.maux
if (iglobal_bz_ext ne 0) then iglobal_bz_ext=iglobal_bz_ext-dim.mvar-dim.maux
if (iglobal_jx_ext ne 0) then iglobal_jx_ext=iglobal_jx_ext-dim.mvar-dim.maux
if (iglobal_jy_ext ne 0) then iglobal_jy_ext=iglobal_jy_ext-dim.mvar-dim.maux
if (iglobal_jz_ext ne 0) then iglobal_jz_ext=iglobal_jz_ext-dim.mvar-dim.maux
if (iglobal_ex_ext ne 0) then iglobal_ex_ext=iglobal_ex_ext-dim.mvar-dim.maux
if (iglobal_ey_ext ne 0) then iglobal_ey_ext=iglobal_ey_ext-dim.mvar-dim.maux
if (iglobal_ez_ext ne 0) then iglobal_ez_ext=iglobal_ez_ext-dim.mvar-dim.maux

varcontent[iglobal_bx_ext].variable   = 'Magnetic field (bx_ext)'
varcontent[iglobal_bx_ext].idlvar     = 'bx_ext'
varcontent[iglobal_bx_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_bx_ext].idlvarloc  = 'bx_ext_loc'
varcontent[iglobal_bx_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_by_ext].variable   = 'Magnetic field (by_ext)'
varcontent[iglobal_by_ext].idlvar     = 'by_ext'
varcontent[iglobal_by_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_by_ext].idlvarloc  = 'by_ext_loc'
varcontent[iglobal_by_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_bz_ext].variable   = 'Magnetic field (bz_ext)'
varcontent[iglobal_bz_ext].idlvar     = 'bz_ext'
varcontent[iglobal_bz_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_bz_ext].idlvarloc  = 'bz_ext_loc'
varcontent[iglobal_bz_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_jx_ext].variable   = 'Current density (jx_ext)'
varcontent[iglobal_jx_ext].idlvar     = 'jx_ext'
varcontent[iglobal_jx_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_jx_ext].idlvarloc  = 'jx_ext_loc'
varcontent[iglobal_jx_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_jy_ext].variable   = 'Current density (jy_ext)'
varcontent[iglobal_jy_ext].idlvar     = 'jy_ext'
varcontent[iglobal_jy_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_jy_ext].idlvarloc  = 'jy_ext_loc'
varcontent[iglobal_jy_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_jz_ext].variable   = 'Current density (jz_ext)'
varcontent[iglobal_jz_ext].idlvar     = 'jz_ext'
varcontent[iglobal_jz_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_jz_ext].idlvarloc  = 'jz_ext_loc'
varcontent[iglobal_jz_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_ex_ext].variable   = 'Electromotive force (ex_ext)'
varcontent[iglobal_ex_ext].idlvar     = 'ex_ext'
varcontent[iglobal_ex_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_ex_ext].idlvarloc  = 'ex_ext_loc'
varcontent[iglobal_ex_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_ey_ext].variable   = 'Electromotive force (ey_ext)'
varcontent[iglobal_ey_ext].idlvar     = 'ey_ext'
varcontent[iglobal_ey_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_ey_ext].idlvarloc  = 'ey_ext_loc'
varcontent[iglobal_ey_ext].idlinitloc = INIT_SCALAR_LOC

varcontent[iglobal_ez_ext].variable   = 'Electromotive force (ez_ext)'
varcontent[iglobal_ez_ext].idlvar     = 'ez_ext'
varcontent[iglobal_ez_ext].idlinit    = INIT_SCALAR
varcontent[iglobal_ez_ext].idlvarloc  = 'ez_ext_loc'
varcontent[iglobal_ez_ext].idlinitloc = INIT_SCALAR_LOC

; ZERO out default 'should never be used' definition
; will have been filled in where i?????? has not been
; set above.

varcontent[0].variable = 'UNKNOWN'
varcontent[0].idlvar   = 'UNKNOWN'
varcontent[0].idlinit  = '0.'
varcontent[0].skip  = 0

return, varcontent

END
