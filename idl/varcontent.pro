if (par.lwrite_aux ne 0) then totalvars=nvar+naux else totalvars=nvar
varcontent=REPLICATE({varcontent_all, variable:'UNKNOWN', $ 
                                      idlvar:'dummy', $
                                      idlinit:'fltarr(mx,my,mz)*one', $
                                      idlvarloc:'dummy_loc', $
                                      idlinitloc:'fltarr(mxloc,myloc,mzloc)*one', $
                                      skip:0},totalvars+1)

;
; Declare ALL variables that may occur
;
INIT_3VECTOR     = 'fltarr(mx,my,mz,3)*one'
INIT_3VECTOR_LOC = 'fltarr(mxloc,myloc,mzloc,3)*one'
INIT_SCALAR     = 'fltarr(mx,my,mz)*one'
INIT_SCALAR_LOC = 'fltarr(mxloc,myloc,mzloc)*one'

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

varcontent[ient].variable = 'Entropy (ss)'
varcontent[ient].idlvar   = 'ss'
varcontent[ient].idlinit    = INIT_SCALAR
varcontent[ient].idlvarloc= 'ss_loc'
varcontent[ient].idlinitloc = INIT_SCALAR_LOC

varcontent[iaa].variable = 'Magnetic vector potential (aa)'
varcontent[iaa].idlvar   = 'aa'
varcontent[iaa].idlinit    = INIT_3VECTOR
varcontent[iaa].idlvarloc= 'aa_loc'
varcontent[iaa].idlinitloc = INIT_3VECTOR_LOC
varcontent[iaa].skip  = 2

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

varcontent[ilncc].variable = 'Log passive scalar (lncc)'
varcontent[ilncc].idlvar   = 'lncc'
varcontent[ilncc].idlinit    = INIT_SCALAR
varcontent[ilncc].idlvarloc= 'lncc_loc'
varcontent[ilncc].idlinitloc = INIT_SCALAR_LOC

varcontent[iuud].variable = 'Dust velocity (uud)'
varcontent[iuud].idlvar   = 'uud'
varcontent[iuud].idlinit    = INIT_3VECTOR
varcontent[iuud].idlvarloc= 'uud_loc'
varcontent[iuud].idlinitloc = INIT_3VECTOR_LOC
varcontent[iuud].skip     = 2

varcontent[ilnrhod].variable = 'Dust log density (lnrhod)'
varcontent[ilnrhod].idlvar   = 'lnrhod'
varcontent[ilnrhod].idlinit    = INIT_SCALAR
varcontent[ilnrhod].idlvarloc= 'lnrhod_loc'
varcontent[ilnrhod].idlinitloc = INIT_SCALAR_LOC

if (par.lwrite_aux ne 0) then begin
    varcontent[iQrad].variable = 'Radiation (Qrad)'
    varcontent[iQrad].idlvar   = 'Qrad'
    varcontent[iQrad].idlinit    = INIT_SCALAR
    varcontent[iQrad].idlvarloc= 'Qrad_loc'
    varcontent[iQrad].idlinitloc = INIT_SCALAR_LOC
    
    ;varcontent[iSrad].variable = 'Radiation (Srad)'
    ;varcontent[iSrad].idlvar   = 'Srad'
    ;varcontent[iSrad].idlinit    = INIT_SCALAR
    ;varcontent[iSrad].idlvarloc= 'Srad_loc'
    ;varcontent[iSrad].idlinitloc = INIT_SCALAR_LOC
    
    ;varcontent[ikappa].variable = 'Radiation (kappa)'
    ;varcontent[ikappa].idlvar   = 'kappa'
    ;varcontent[ikappa].idlinit    = INIT_SCALAR
    ;varcontent[ikappa].idlvarloc= 'kappa_loc'
    ;varcontent[ikappa].idlinitloc = INIT_SCALAR_LOC
    
; May need special condition as can be maux or mvar variable?
    varcontent[iTT].variable   = 'Temperature (TT)'
    varcontent[iTT].idlvar     = 'TT'
    varcontent[iTT].idlinit    = INIT_SCALAR
    varcontent[iTT].idlvarloc  = 'TT_loc'
    varcontent[iTT].idlinitloc = INIT_SCALAR_LOC

    varcontent[iyH].variable   = 'Hydrogen ionization fraction (yH)'
    varcontent[iyH].idlvar     = 'yH'
    varcontent[iyH].idlinit    = INIT_SCALAR
    varcontent[iyH].idlvarloc  = 'yH_loc'
    varcontent[iyH].idlinitloc = INIT_SCALAR_LOC

    varcontent[ishock].variable = 'Shock characteristic (nu_shock)'
    varcontent[ishock].idlvar   = 'nu_shock'
    varcontent[ishock].idlinit    = INIT_SCALAR
    varcontent[ishock].idlvarloc= 'nu_shock_loc'
    varcontent[ishock].idlinitloc = INIT_SCALAR_LOC
end

varcontent[0].variable = 'UNKNOWN'
varcontent[0].idlvar   = 'dummy'
varcontent[0].idlinit  = '0.'
varcontent[0].skip  = 0
