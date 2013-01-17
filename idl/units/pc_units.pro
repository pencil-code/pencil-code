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
              datadir=datadir, QUIET=QUIET, HELP=HELP
;
COMPILE_OPT IDL2,HIDDEN
; 
  if (keyword_set(HELP)) then begin
    print, "Usage:"
    print, ""
    print, "pc_units, object=object, datadir=datadir, proc=proc, /PRINT, /QUIET, /HELP"
    print, ""
    print, "Returns a structure containing the effective units for various quantities."
    print, ""
    print, "   object: structure that contains some commonly used units."
    print, "  symbols: structure that contains the unit symbols as strings."
    print, "  datadir: string that specifies the data directory."
    print, ""
    print, "   /QUIET: instruction not to print any 'helpful' information."
    print, "    /HELP: display this usage information, and exit."
    return
  end
;
; Default data directory
;
  pc_check_math,location='before entry to pc_units'
  if (not keyword_set(datadir)) then datadir = pc_get_datadir()
  if (n_elements(param) eq 0) then $
      pc_read_param, object=param, datadir=datadir, dim=dim, quiet=quiet
;
  temperature = double(param.unit_temperature)
  density = double(param.unit_density)
  length = double(param.unit_length)
  velocity = double(param.unit_velocity)
  if (any (strmatch (tag_names (param), "magnetic", /fold_case))) then magnetic = double(param.unit_magnetic)
;
  if (param.unit_system eq "cgs") then begin
;
    default, magnetic, sqrt(4*double (!pi)/param.mu0 * density) * velocity
    object = { temperature:temperature, $
               density:density, $
               mass:param.unit_density*param.unit_length^3, $
               length:length, $
               velocity:velocity, $
               time:length/velocity, $
               energy:density*velocity^2*length^3, $
               specific_energy:velocity^2, $
               magnetic_field:magnetic, $
               current_density:!Values.D_NaN } ; needs definition
    pc_check_math,location='pc_units - cgs unit calculation'
    tex=texsyms()
    symbols = { temperature:'T', $
                density:tex.varrho, $
                mass:'g', $
                length:'cm', $
                velocity:'cm/s', $
                time:'s', $
                energy:'ergs', $
                specific_energy:'ergs/g', $
                magnetic_field:'G', $
                current_density:'-undefined-' } ; needs definition
;
  end else if (param.unit_system eq "SI") then begin
;
    default, magnetic, sqrt (4*double (!pi)*1e-7/param.mu0 * density) * velocity
    object = { temperature:temperature, $
               density:density, $
               mass:param.unit_density*param.unit_length^3, $
               length:length, $
               velocity:velocity, $
               time:length/velocity, $
               energy:density*velocity^2*length^3, $
               specific_energy:velocity^2, $
               magnetic_field:magnetic, $
               current_density:velocity*sqrt(param.mu0*0.25/double(!Pi)*1.e7*density)/length }
    pc_check_math,location='pc_units - SI unit calculation'
    tex=texsyms()
    symbols = { temperature:'K', $
                density:tex.varrho, $
                mass:'kg', $
                length:'m', $
                velocity:'m/s', $
                time:'s', $
                energy:'J', $
                specific_energy:'J/kg', $
                magnetic_field:'T', $
                current_density:'A/m^2' }
;
  end else begin
;
    print,"pc_units: Unit system tranformations for unit_system='",param.unit_system,"' are not implemented."
;
  end
;
end


