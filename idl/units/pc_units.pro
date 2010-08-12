; $Id$
;
;  Read defined units from params and return them along with a bunch
;  of calculated units in a structure.
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;
;  26-jul-04/tony: coded 
;  
pro pc_units, object=object, symbols=symbols, param=param, dim=dim, $
                   datadir=datadir,QUIET=QUIET,HELP=HELP
;
COMPILE_OPT IDL2,HIDDEN
; 
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_units, object=object,                                                                   "
    print, "               datadir=datadir, proc=proc,                                                 "
    print, "               /PRINT, /QUIET, /HELP                                                       "
    print, "                                                                                           "
    print, "Returns a structure containing the effective units for various quantities in a neat        "
    print, "structure.                                                                                 "
    print, "                                                                                           "
    print, "    datadir: specify the root data directory. Default is './data'                  [string]"
    print, ""
    print, "   object: structure in which to return all the units                          [structure] "
    print, ""
    print, "   /QUIET: instruction not to print any 'helpful' information                              "
    print, "    /HELP: display this usage information, and exit                                        "
    return
  ENDIF
;
; Default data directory
;
pc_check_math,location='before entry to pc_units'
default, datadir, 'data'
if (n_elements(dim) eq 0) then pc_read_dim, datadir=datadir, object=dim, $
    quiet=quiet
if (n_elements(param) eq 0) then pc_read_param, datadir=datadir, object=param, $
    dim=dim,quiet=quiet
;
length=param.unit_length*1D0
temperature=param.unit_temperature*1D0
density=param.unit_density*1D0
velocity=param.unit_velocity*1D0
magnetic=param.unit_magnetic*1D0
;
if param.unit_system eq "cgs" then begin
;
  object=create_struct(['temperature',         $
                        'density',             $
                        'length',              $
                        'velocity',            $
                        'time',                $
                        'energy',              $
                        'specific_energy',     $
                        'magnetic_field'],     $
                        temperature,           $
                        density,               $
                        length,                $
                        velocity,              $
                        length/velocity,       $
                        1D0*density*velocity^2*length^3, $
                        velocity^2,            $
                        magnetic               $
                      )
  pc_check_math,location='pc_units - cgs unit calculation'
  tex=texsyms()
  symbols=create_struct(['temperature',        $
                         'density',            $
                         'length',             $
                         'velocity',           $
                         'time',               $
                         'energy',             $
                         'specific_energy',    $
                         'magnetic_field'],    $
                         'T',                  $
                         tex.varrho,           $
                         'cm',                 $
                         'cm/s',               $
                         's',                  $
                         'ergs',               $
                         'ergs/g',             $
                         'G'                   $
                      )
  
end else if param.unit_system eq "SI" then begin
;
  object=create_struct(['temperature',         $
                        'density',             $
                        'length',              $
                        'velocity',            $
                        'time',                $
                        'energy',              $
                        'specific_energy',     $
                        'magnetic_field'],     $
                        temperature,           $
                        density,               $
                        length,                $
                        velocity,              $
                        length/velocity,       $
                        1D0*density*velocity^2*length^3, $
                        velocity^2,            $
                        magnetic               $
                      )
  pc_check_math,location='pc_units - cgs unit calculation'
  tex=texsyms()
  symbols=create_struct(['temperature',        $
                         'density',            $
                         'length',             $
                         'velocity',           $
                         'time',               $
                         'energy',             $
                         'specific_energy',    $
                         'magnetic_field'],    $
                         'K',                  $
                         tex.varrho,           $
                         'm',                  $
                         's',                  $
                         'm/s',                $
                         'J',                  $
                         'J/kg',               $
                         'T'                   $
                      )
;
end else begin
;
  print,"pc_units: Unit system tranformations for unit_system=",param.unit_system," are not implemented."
;
end
;
end


