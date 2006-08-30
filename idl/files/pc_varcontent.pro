;  $Id: pc_varcontent.pro,v 1.31 2006-08-30 13:28:37 dintrans Exp $
FUNCTION pc_varcontent,datadir=datadir,dim=dim, $
                       param=param,quiet=quiet,scalar=scalar,run2D=run2D
COMPILE_OPT IDL2,HIDDEN
;
;
;
default,ifcr,0

ishock=0
ishock_perp=0
inp=0
ivpxsum=0
irhop=0
igpotselfx=0 & igpotselfy=0 & igpotselfz=0
ipsi_real=0
ipsi_imag=0
ikapparho=0
; 
;  Read the positions of variables in f
;  Can't just use `@data/index', as the data directory may have a different name
;
if (n_elements(dim) eq 0) then pc_read_dim,obj=dim,datadir=datadir,quiet=quiet
if (n_elements(param) eq 0) then pc_read_param,obj=param,datadir=datadir, $
    dim=dim,quiet=quiet

default,datadir,'data'
cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),"\n"'+"' "+datadir+'/index.pro'
spawn, cmd, result
res = flatten_strings(result) 

;res=''
;get_lun,indexfile
;openr,indexfile,datadir+'/index.pro'
;repeat begin
;readf,indexfile,res
if (execute(res) ne 1) then $
    message, 'There was a problem with index.pro', /INFO
;endrep until eof(indexfile)

;close,indexfile
;free_lun,indexfile

;
; VARCONTENT STRUCTURE DESCRIPTION
;
; variable (string)
;   Human readable name for the variable
;
; idlvar (string)
;   Name of the variable (usually in the IDL global namespace)
;   in which the variables data will be stored
;
; idlinit (string)
;   IDL command to initialise the storage variable ready for
;   reading in from a file
;
; idlvarloc (string)
;   As idlvar but used when two grid sizes are used eg. global mesh
;   and processor mesh (local -> loc). And used in processes such
;   as in rall.pro.  Eg. uses mesh sizes of (mxloc,myloc,mzloc)
;
; idlinitloc (string)
;   Again as idlinit but used when two mesh sizes are required at once.
;   see idlvarloc

; How many variables are expected to be stored in the var file?


if (param.lwrite_aux ne 0) then totalvars=dim.mvar+dim.maux else totalvars=dim.mvar

; Make an array of structures in which to store their descriptions
; index zero is kept as a dummy entry.
varcontent=REPLICATE({varcontent_all, variable:'UNKNOWN', $ 
                                      idlvar:'dummy', $
                                      idlinit:'fltarr(mx,my,mz)*one', $
                                      idlvarloc:'dummy_loc', $
                                      idlinitloc:'fltarr(mxloc,myloc,mzloc)*one', $
                                      skip:0},totalvars+1)
;for i=1L,totalvars do begin
;  varcontent[i].idlvar='dummy'+str(i)
;endfor
;
; Declare ALL variables that MAY OCCUR
;

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

; For EVERY POSSIBLE variable in a var file, store a
; description of the variable in an indexed array of structures
; where the indexes line up with those in the saved f array

; Any variable not stored should have iXXXXXX set to zero
; and will only update the dummy index zero entry

; DO mvar VARIABLES FIRST

varcontent[iuu].variable   = 'Velocity (uu)'
varcontent[iuu].idlvar     = 'uu'
varcontent[iuu].idlinit    = INIT_3VECTOR
varcontent[iuu].idlvarloc  = 'uu_loc'
varcontent[iuu].idlinitloc = INIT_3VECTOR_LOC
varcontent[iuu].skip       = 2

varcontent[ilnrho].variable   = 'Log density (lnrho)'
varcontent[ilnrho].idlvar     = 'lnrho'
varcontent[ilnrho].idlinit    = INIT_SCALAR
varcontent[ilnrho].idlvarloc  = 'lnrho_loc'
varcontent[ilnrho].idlinitloc = INIT_SCALAR_LOC

varcontent[iss].variable = 'Entropy (ss)'
varcontent[iss].idlvar   = 'ss'
varcontent[iss].idlinit    = INIT_SCALAR
varcontent[iss].idlvarloc= 'ss_loc'
varcontent[iss].idlinitloc = INIT_SCALAR_LOC

varcontent[iaa].variable = 'Magnetic vector potential (aa)'
varcontent[iaa].idlvar   = 'aa'
varcontent[iaa].idlinit    = INIT_3VECTOR
varcontent[iaa].idlvarloc= 'aa_loc'
varcontent[iaa].idlinitloc = INIT_3VECTOR_LOC
varcontent[iaa].skip  = 2

default,iaatest,0
varcontent[iaatest].variable = 'Testfield vector potential (aatest)'
varcontent[iaatest].idlvar   = 'aatest'
varcontent[iaatest].idlinit    = 'fltarr(mx,my,mz,3*9)*one'
varcontent[iaatest].idlvarloc= 'aatest_loc'
varcontent[iaatest].idlinitloc = 'fltarr(mxloc,myloc,mzloc,3*9)*one'
varcontent[iaatest].skip  = 26+9

varcontent[ifx].variable = 'Radiation vector ?something? (ff)'
varcontent[ifx].idlvar   = 'ff'
varcontent[ifx].idlinit    = INIT_3VECTOR
varcontent[ifx].idlvarloc= 'ff_loc'
varcontent[ifx].idlinitloc = INIT_3VECTOR_LOC
varcontent[ifx].skip  = 2

varcontent[ie].variable = 'Radiation scalar ?something? (ee)'
varcontent[ie].idlvar   = 'ee'
varcontent[ie].idlinit    = INIT_SCALAR
varcontent[ie].idlvarloc= 'ee_loc'
varcontent[ie].idlinitloc = INIT_SCALAR_LOC

varcontent[icc].variable = 'Passive scalar (cc)'
varcontent[icc].idlvar   = 'cc'
varcontent[icc].idlinit    = INIT_SCALAR
varcontent[icc].idlvarloc= 'cc_loc'
varcontent[icc].idlinitloc = INIT_SCALAR_LOC

varcontent[ilncc].variable = 'Log passive scalar (lncc)'
varcontent[ilncc].idlvar   = 'lncc'
varcontent[ilncc].idlinit    = INIT_SCALAR
varcontent[ilncc].idlvarloc= 'lncc_loc'
varcontent[ilncc].idlinitloc = INIT_SCALAR_LOC

default,iXX_chiral,0
varcontent[iXX_chiral].variable = 'XX_chiral'
varcontent[iXX_chiral].idlvar   = 'XX_chiral'
varcontent[iXX_chiral].idlinit    = INIT_SCALAR
varcontent[iXX_chiral].idlvarloc= 'XX_chiral_loc'
varcontent[iXX_chiral].idlinitloc = INIT_SCALAR_LOC

default,iYY_chiral,0
varcontent[iYY_chiral].variable = 'YY_chiral'
varcontent[iYY_chiral].idlvar   = 'YY_chiral'
varcontent[iYY_chiral].idlinit    = INIT_SCALAR
varcontent[iYY_chiral].idlvarloc= 'YY_chiral_loc'
varcontent[iYY_chiral].idlinitloc = INIT_SCALAR_LOC

default,ispecial,0
varcontent[ispecial].variable = 'special'
varcontent[ispecial].idlvar   = 'special'
varcontent[ispecial].idlinit    = INIT_SCALAR
varcontent[ispecial].idlvarloc= 'special_loc'
varcontent[ispecial].idlinitloc = INIT_SCALAR_LOC

varcontent[iecr].variable = 'Cosmic ray energy density (ecr)'
varcontent[iecr].idlvar   = 'ecr'
varcontent[iecr].idlinit    = INIT_SCALAR
varcontent[iecr].idlvarloc= 'ecr_loc'
varcontent[iecr].idlinitloc = INIT_SCALAR_LOC

varcontent[ifcr].variable = 'Cosmic ray energy flux (fcr)'
varcontent[ifcr].idlvar   = 'fcr'
varcontent[ifcr].idlinit    = INIT_3VECTOR
varcontent[ifcr].idlvarloc= 'fcr_loc'
varcontent[ifcr].idlinitloc = INIT_3VECTOR_LOC
varcontent[ifcr].skip  = 2

varcontent[ipsi_real].variable = 'Wave function (real part)'
varcontent[ipsi_real].idlvar   = 'psi_real'
varcontent[ipsi_real].idlinit    = INIT_SCALAR
varcontent[ipsi_real].idlvarloc= 'psi_real_loc'
varcontent[ipsi_real].idlinitloc = INIT_SCALAR_LOC

varcontent[ipsi_imag].variable = 'Wave function (imaginary part)'
varcontent[ipsi_imag].idlvar   = 'psi_imag'
varcontent[ipsi_imag].idlinit    = INIT_SCALAR
varcontent[ipsi_imag].idlvarloc= 'psi_imag_loc'
varcontent[ipsi_imag].idlinitloc = INIT_SCALAR_LOC


dustcount=n_elements(iuud) 
if (dustcount gt 0L) then begin
  if (keyword_set(scalar)) then begin
    for i=0,dustcount-1 do begin
     istr=strcompress(string(i),/remove_all)
     varcontent[iuud[i]].variable = 'Dust velocity  (uud['+istr+'])'
     varcontent[iuud[i]].idlvar   = 'uud'+istr
     varcontent[iuud[i]].idlinit  = 'fltarr(mx,my,mz,3)*one' 
     varcontent[iuud[i]].idlvarloc= 'uud'+istr+'_loc'
     varcontent[iuud[i]].idlinitloc = 'fltarr(mxloc,myloc,mzloc,3)*one'
    endfor
  endif else begin
    varcontent[iuud[0]].variable = 'Dust velocity  (uud)'
    varcontent[iuud[0]].idlvar   = 'uud'
    varcontent[iuud[0]].idlinit  = 'fltarr(mx,my,mz,3,'+str(dustcount)+')*one' 
    varcontent[iuud[0]].idlvarloc= 'uud_loc'
    varcontent[iuud[0]].idlinitloc = 'fltarr(mxloc,myloc,mzloc,3,'+str(dustcount)+')*one'
    varcontent[iuud[0]].skip     = (dustcount * 3) - 1
  endelse
endif

dustcount=n_elements(ind)
if (dustcount gt 0L) then begin
  if (keyword_set(scalar)) then begin
    for i=0,dustcount-1 do begin
     istr=strcompress(string(i),/remove_all)
     varcontent[ind[0]].variable = 'Dust number density (nd'+istr+')'
     varcontent[ind[0]].idlvar   = 'nd'+istr
     varcontent[ind[0]].idlinit  = 'fltarr(mx,my,mz)*one' 
     varcontent[ind[0]].idlvarloc= 'nd'+istr+'_loc'
     varcontent[ind[0]].idlinitloc = 'fltarr(mxloc,myloc,mzloc)*one'
    endfor
  endif else begin
    varcontent[ind[0]].variable = 'Dust number density (nd)'
    varcontent[ind[0]].idlvar   = 'nd'
    varcontent[ind[0]].idlinit  = 'fltarr(mx,my,mz,'+str(dustcount)+')*one' 
    varcontent[ind[0]].idlvarloc= 'nd_loc'
    varcontent[ind[0]].idlinitloc = 'fltarr(mxloc,myloc,mzloc,'+str(dustcount)+')*one'
    varcontent[ind[0]].skip     = dustcount - 1
  endelse
endif

dustcount=n_elements(imd)
if (dustcount gt 0L) then begin
  if (keyword_set(scalar)) then begin
    for i=0,dustcount-1 do begin
     istr=strcompress(string(i),/remove_all)
      varcontent[imd[i]].variable = 'Dust density (md'+istr+')'
      varcontent[imd[i]].idlvar   = 'md'+istr
      varcontent[imd[i]].idlinit  = 'fltarr(mx,my,mz)*one' 
      varcontent[imd[i]].idlvarloc= 'md'+istr+'_loc'
      varcontent[imd[i]].idlinitloc = 'fltarr(mxloc,myloc,mzloc)*one'
    endfor
  endif else begin
    varcontent[imd[0]].variable = 'Dust density (md)'
    varcontent[imd[0]].idlvar   = 'md'
    varcontent[imd[0]].idlinit  = 'fltarr(mx,my,mz,'+str(dustcount)+')*one' 
    varcontent[imd[0]].idlvarloc= 'md_loc'
    varcontent[imd[0]].idlinitloc = 'fltarr(mxloc,myloc,mzloc,'+str(dustcount)+')*one'
    varcontent[imd[0]].skip     = dustcount - 1
  endelse
endif

dustcount=n_elements(imi)
if (dustcount gt 0L) then begin
  if (keyword_set(scalar)) then begin
    for i=0,dustcount-1 do begin
      istr=strcompress(string(i),/remove_all)
      varcontent[imi[i]].variable = 'Ice density (mi'+istr+')'
      varcontent[imi[i]].idlvar   = 'mi'+istr
      varcontent[imi[i]].idlinit  = 'fltarr(mx,my,mz)*one' 
      varcontent[imi[i]].idlvarloc= 'mi'+istr+'_loc'
      varcontent[imi[i]].idlinitloc = 'fltarr(mxloc,myloc,mzloc)*one'
    endfor
  endif else begin
    varcontent[imi[0]].variable = 'Ice density (mi)'
    varcontent[imi[0]].idlvar   = 'mi'
    varcontent[imi[0]].idlinit  = 'fltarr(mx,my,mz,'+str(dustcount)+')*one' 
    varcontent[imi[0]].idlvarloc= 'mi_loc'
    varcontent[imi[0]].idlinitloc = 'fltarr(mxloc,myloc,mzloc,'+str(dustcount)+')*one'
    varcontent[imi[0]].skip     = dustcount - 1
  endelse
endif

varcontent[igg].variable = 'Gravitational acceleration (gg)'
varcontent[igg].idlvar   = 'gg'
varcontent[igg].idlinit    = INIT_3VECTOR
varcontent[igg].idlvarloc= 'gg_loc'
varcontent[igg].idlinitloc = INIT_3VECTOR_LOC
varcontent[igg].skip     = 2


; Special condition as can be maux or mvar variable
if ((ilnTT le dim.mvar) or (param.lwrite_aux ne 0)) then begin
    varcontent[ilnTT].variable   = 'Log temperature (lnTT)'
    varcontent[ilnTT].idlvar     = 'lnTT'
    varcontent[ilnTT].idlinit    = INIT_SCALAR
    varcontent[ilnTT].idlvarloc  = 'lnTT_loc'
    varcontent[ilnTT].idlinitloc = INIT_SCALAR_LOC
end

; THEN DO maux VARIABLES 
; ** ONLY IF THEY HAVE BEEN SAVED **
if (param.lwrite_aux ne 0) then begin
  varcontent[iQrad].variable = 'Radiative heating rate (Qrad)'
  varcontent[iQrad].idlvar   = 'Qrad'
  varcontent[iQrad].idlinit    = INIT_SCALAR
  varcontent[iQrad].idlvarloc= 'Qrad_loc'
  varcontent[iQrad].idlinitloc = INIT_SCALAR_LOC

  varcontent[ikapparho].variable = 'Opacity (kapparho)'
  varcontent[ikapparho].idlvar   = 'kapparho'
  varcontent[ikapparho].idlinit    = INIT_SCALAR
  varcontent[ikapparho].idlvarloc= 'kapparho_loc'
  varcontent[ikapparho].idlinitloc = INIT_SCALAR_LOC

  varcontent[iFrad].variable = 'Radiative flux (Frad)'
  varcontent[iFrad].idlvar   = 'Frad'
  varcontent[iFrad].idlinit    = INIT_3VECTOR
  varcontent[iFrad].idlvarloc= 'Frad_loc'
  varcontent[iFrad].idlinitloc = INIT_3VECTOR_LOC
  varcontent[iFrad].skip  = 2

  varcontent[iyH].variable   = 'Hydrogen ionization fraction (yH)'
  varcontent[iyH].idlvar     = 'yH'
  varcontent[iyH].idlinit    = INIT_SCALAR
  varcontent[iyH].idlvarloc  = 'yH_loc'
  varcontent[iyH].idlinitloc = INIT_SCALAR_LOC

  varcontent[ishock].variable = 'Shock Profile (shock)'
  varcontent[ishock].idlvar   = 'shock'
  varcontent[ishock].idlinit    = INIT_SCALAR
  varcontent[ishock].idlvarloc= 'shock_loc'
  varcontent[ishock].idlinitloc = INIT_SCALAR_LOC

  varcontent[ishock_perp].variable = 'B-Perp Shock Profile (shock_perp)'
  varcontent[ishock_perp].idlvar   = 'shock_perp'
  varcontent[ishock_perp].idlinit    = INIT_SCALAR
  varcontent[ishock_perp].idlvarloc= 'shock_perp_loc'
  varcontent[ishock_perp].idlinitloc = INIT_SCALAR_LOC

  default,icooling,0
  varcontent[icooling].variable = 'Cooling Term (cooling)'
  varcontent[icooling].idlvar   = 'cooling'
  varcontent[icooling].idlinit    = INIT_SCALAR
  varcontent[icooling].idlvarloc= 'cooling_loc'
  varcontent[icooling].idlinitloc = INIT_SCALAR_LOC

  default,icooling2,0
  varcontent[icooling2].variable = 'Applied Cooling Term (cooling)'
  varcontent[icooling2].idlvar   = 'cooling2'
  varcontent[icooling2].idlinit    = INIT_SCALAR
  varcontent[icooling2].idlvarloc= 'cooling2_loc'
  varcontent[icooling2].idlinitloc = INIT_SCALAR_LOC

  varcontent[inp].variable   = 'Particle number (np)'
  varcontent[inp].idlvar     = 'np'
  varcontent[inp].idlinit    = INIT_SCALAR
  varcontent[inp].idlvarloc  = 'np_loc'
  varcontent[inp].idlinitloc = INIT_SCALAR_LOC

  varcontent[irhop].variable   = 'Particle mass density (rhop)'
  varcontent[irhop].idlvar     = 'rhop'
  varcontent[irhop].idlinit    = INIT_SCALAR
  varcontent[irhop].idlvarloc  = 'rhop_loc'
  varcontent[irhop].idlinitloc = INIT_SCALAR_LOC

  varcontent[ivpxsum].variable   = 'Sum of particle velocities (vvpsum)'
  varcontent[ivpxsum].idlvar     = 'vvpsum'
  varcontent[ivpxsum].idlinit    = INIT_3VECTOR
  varcontent[ivpxsum].idlvarloc  = 'vvpsum_loc'
  varcontent[ivpxsum].idlinitloc = INIT_3VECTOR_LOC
  varcontent[ivpxsum].skip       = 2

  varcontent[ipotself].variable   = 'Self gravity potential'
  varcontent[ipotself].idlvar     = 'potself'
  varcontent[ipotself].idlinit    = INIT_SCALAR
  varcontent[ipotself].idlvarloc  = 'potself_loc'
  varcontent[ipotself].idlinitloc = INIT_SCALAR_LOC

  varcontent[igpotselfx].variable   = 'Gradient of self gravity potential'
  varcontent[igpotselfx].idlvar     = 'gpotself'
  varcontent[igpotselfx].idlinit    = INIT_3VECTOR
  varcontent[igpotselfx].idlvarloc  = 'gpotself_loc'
  varcontent[igpotselfx].idlinitloc = INIT_3VECTOR_LOC
  varcontent[igpotselfx].skip       = 2
endif

if keyword_set(scalar) then begin
  for i = 1, totalvars do begin
    if (varcontent[i].skip eq 2) then begin
      varcontent[i+2].variable  = varcontent[i].variable + ' 3rd component' 
      varcontent[i+1].variable  = varcontent[i].variable + ' 2nd component' 
      varcontent[i  ].variable  = varcontent[i].variable + ' 1st component' 
      varcontent[i+2].idlvar    = varcontent[i].idlvar + '3' 
      varcontent[i+1].idlvar    = varcontent[i].idlvar + '2' 
      varcontent[i  ].idlvar    = varcontent[i].idlvar + '1' 
      varcontent[i+2].idlvarloc = varcontent[i].idlvarloc + '3' 
      varcontent[i+1].idlvarloc = varcontent[i].idlvarloc + '2' 
      varcontent[i  ].idlvarloc = varcontent[i].idlvarloc + '1' 
      varcontent[i:i+2].idlinit = INIT_SCALAR
      varcontent[i:i+2].idlinitloc = INIT_SCALAR_LOC
      varcontent[i:i+2].skip    = 0
      i=i+2
    endif   
  endfor
endif

; ZERO out default 'should never be used' definition
; will have been filled in where i?????? has not been
; set above.
varcontent[0].variable = 'UNKNOWN'
varcontent[0].idlvar   = 'UNKNOWN'
varcontent[0].idlinit  = '0.'
varcontent[0].skip  = 0

return,varcontent

END
