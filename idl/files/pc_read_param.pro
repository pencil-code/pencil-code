; $Id: pc_read_param.pro,v 1.4 2004-05-05 17:10:31 mee Exp $
;
;   Read param.nml
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2004-05-05 17:10:31 $
;  $Revision: 1.4 $
;
;  27-nov-02/tony: coded mostly from Wolgang's start.pro
;
;  REQUIRES: external 'nl2idl' perl script (WD)
;  
;  
function param
COMPILE_OPT IDL2,HIDDEN
; Dummy to keep IDL from complaining. The real param() routine will be
; compiled below
end
;
pro pc_read_param,lhydro=lhydro,ldensity=ldensity,lgravz=lgravz,lgravr=lgravr,lentropy=lentropy, $
                  lmagnetic=lmagnetic, lradiation=lradiation,lpscalar=lpscalar, lforcing=lforcing, $
                  lshear=lshear, $
                  cs0=cs0, cs20=cs20, rho0=rho0, lnrho0=lnrho0, gamma=gamma, gamma1=gamma1, $
                  z1=z1, z2=z2, z3=z3, zref=zref, ztop=ztop, gravz=gravz, $
                  mpoly0=mpoly0, mpoly1=mpoly1, mpoly2=mpoly2, isothtop=isothtop, $
                  object=object, $
                  datadir=datadir,PRINT=PRINT,QUIET=QUIET,HELP=HELP
COMPILE_OPT IDL2,HIDDEN
  COMMON pc_precision, zero, one
; If no meaningful parameters are given show some help!
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_read_param, lhydro=lhydro, ldensity=ldensity, lgravz=lgravz, lgravr=lgravr,            " 
    print, "               lentropy=lentropy, lmagnetic=lmagnetic, lradiation=lradiation,             "
    print, "               lpscalar=lpscalar, lforcing=lforcing, lshear=lshear,                       "
    print, "               cs0=cs0, cs20=cs20, rho0=rho0, lnrho0=lnrho0, gamma=gamma, gamma1=gamma1,  "
    print, "               z1=z1, z2=z2, z3=z3, zref=zref, ztop=ztop, gravz=gravz,                    "
    print, "               mpoly0=mpoly0, mpoly1=mpoly1, mpoly2=mpoly2, isothtop=isothtop,            "
    print, "               object=object,                                                             "
    print, "               datadir=datadir, proc=proc,                                                "
    print, "               /PRINT, /QUIET, /HELP                                                       "
    print, "                                                                                           "
    print, "Returns the parameters of a Pencil-Code run. Returns zeros and empty in all variables on   "
    print, "fail ure.                                                                                  "
    print, "                                                                                           "
    print, "    datadir: specify the root data directory. Default is './data'                  [string]"
    print, ""
    print, "     lhydro: boolean indicating teh presence of Hydro module                         [long]"
    print, "   ldensity: boolean indicating the presence of Density module                       [long]"
    print, "     lgravz: boolean indicating the presence fo Gravity (gravz) module               [long]"
    print, "     lgravr: boolean indicating the presence fo Gravity (gravr) module               [long]"
    print, "   lentropy: boolean indicating the presence fo Entropy module                       [long]"
    print, "  lmagnetic: boolean indicating the presence fo Magnetic module                      [long]"
    print, " lradiation: boolean indicating the presence fo Radiation module                     [long]"
    print, "   lpscalar: boolean indicating the presence fo Pscalar module                       [long]"
    print, "   lforcing: boolean indicating the presence fo Forcing module                       [long]"
    print, "     lshear: boolean indicating the presence of the Shear module                     [long]"
    print, "        cs0: unit sound speed (code units)                                         [single]"
    print, "       cs20: cs0 squared                                                           [single]"
    print, "       rho0: unit density (code units)                                             [single]"
    print, "     lnrho0: natural log of rho0                                                   [single]"
    print, "      gamma: polytropic index for ideal gas = c_p / c_v                            [single]"
    print, "     gamma1: gamma - 1                                                             [single]"
    print, "         z1:                                                                       [single]"
    print, "         z2:                                                                       [single]"
    print, "         z3:                                                                       [single]"
    print, "       zref:                                                                       [single]"
    print, "       ztop:                                                                       [single]"
    print, "      gravz:                                                                       [single]"
    print, "     mpoly0:                                                                       [single]"
    print, "     mpoly1:                                                                       [single]"
    print, "     mpoly2:                                                                       [single]"
    print, "   isothtop:                                                                       [single]"
    print, ""
    print, "   object: optional structure in which to return all the above as tags          [structure] "
    print, ""
    print, "   /PRINT: instruction to print all variables to standard output                            "
    print, "   /QUIET: instruction not to print any 'helpful' information                               "
    print, "    /HELP: display this usage information, and exit                                         "
    return
  ENDIF

; Default data directory

default, datadir, 'data'

pc_set_precision   ; check precision is set
;
; Initialize / set default returns for ALL variables
;
lhydro    = 0L
ldensity  = 0L
lgravz    = 0L
lgravr    = 0L
lentropy  = 0L
lmagnetic = 0L
lradiation= 0L
lpscalar  = 0L
lforcing  = 0L
lshear    = 0L
cs0       = zero
rho0      = zero
gamma     = zero
gamma1    = zero
cs20      = zero
lnrho     = zero
z1        = zero
z2        = zero
zref      = zero
gravz     = zero
ztop      = zero
z3        = zero
mpoly0    = zero
mpoly1    = zero
mpoly2    = zero
isothtop  = zero
; Get a unit number
; GET_LUN, file   DON'T NEED TO OPEN A FILE!!

; Build the full path and filename
filename=datadir+'/param.nml'   

; Check for existance and read the data
dummy = findfile(filename, COUNT=found)
if (found gt 0) then begin
    IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' + filename + '...'
    spawn, 'for d in . $TMPDIR $TMP /tmp /var/tmp; do if [ -d $d -a -w $d ]; then echo $d; fi; done', /SH, result
    if (strlen(result[0])) le 0 then begin
      message, "Can't find writeable directory for temporary files"
    endif else begin
      tmpdir = result[0]
    endelse
    tmpfile = tmpdir+'/param.pro'
    ;; Write content of param.nml to temporary file:
    spawn, '$PENCIL_HOME/bin/nl2idl -m '+filename+'> ' $
         + tmpfile , result
    ;; Compile that file. Should be easy, but is incredibly awkward, as
    ;; there is no way in IDL to compile a given file at run-time
    ;; outside the command line:
    ;; Save old path and pwd
    _path = !path
    cd, tmpdir, CURRENT=_pwd
    !path = '.:'
    resolve_routine, 'param', /IS_FUNCTION
    ;; Restore old path and pwd
    !path = _path & cd, _pwd
    ;; Delete temporary file
    ; file_delete, tmpfile      ; not in IDL <= 5.3
    spawn, 'rm -f '+tmpfile, /SH
    object = param()

    x0=object.xyz0[0] & y0=object.xyz0[1] & z0=object.xyz0[2]
    Lx=object.Lxyz[0] & Ly=object.Lxyz[1] & Lz=object.Lxyz[2]
                                ;
    lhydro    = object.lhydro
    ldensity  = object.ldensity
    lgravz    = object.lgravz
    lgravr    = object.lgravr
    lentropy  = object.lentropy
    lmagnetic = object.lmagnetic
    lradiation= object.lradiation
    lpscalar  = object.lpscalar
    lforcing  = object.lforcing
    lshear    = object.lshear
    
                                ;
    if (ldensity) then begin
        cs0=object.cs0 & rho0=object.rho0
        gamma=object.gamma & gamma1=gamma-1.
        cs20 = cs0^2 & lnrho0 = alog(rho0)
    endif
                                ;
    if (lgravz) then begin
        z1=object.z1 & z2=object.z2
        zref=object.zref
        gravz=object.gravz
        ;ztop=object.ztop & z3=object.z3
    endif
                                ;
    if (lentropy) then begin
        mpoly0=object.mpoly0 & mpoly1=object.mpoly1
        mpoly2=object.mpoly2 & isothtop=object.isothtop
    endif
endif else begin
    message, 'Warning: cannot find file '+ filename
endelse

; If requested print a summary
if keyword_set(PRINT) then begin
    print, 'For GLOBAL calculation domain:'
    print, '    NO SUMMARY INFORMATION CONFIGURED - edit pc_read_param.pro'
endif

end


