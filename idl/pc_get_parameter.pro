;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_get_parameter.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Returns parameters in Pencil-units, after checking for their existence.
;   Checks first in run.in and then in start.in parameter lists.
;
;  Parameters:
;   * param       Name of the parameter to be returned or "" for setting up.
;   * label       Optional error label to be printed, if parameter is unfound.
;   * missing     Optional label of missing parameter to be printed, if unfound.
;   * datadir     Given data directory for loading the parameter structures.
;   * run_param   'run_pars' will be cached, and loaded if necessary.
;   * start_param 'start_pars' will be cached, and loaded if necessary.
;
;  Examples: (in ascending order of efficiency)
;  ============================================
;
;   Use cached parameter structures that are initially loaded once:
;   IDL> print, pc_get_parameter ('nu'[, datadir=datadir])
;   IDL> print, pc_get_parameter ('eta')
;   IDL> print, pc_get_parameter ('mu0')
;
;   Use given parameter structures and cache them:
;   IDL> print, pc_get_parameter ('nu', start_param=start_param, run_param=run_param)
;   IDL> print, pc_get_parameter ('eta')
;   IDL> print, pc_get_parameter ('mu0')
;


; Cleanup parameter cache, if requested, and return selected parameter.
function pc_get_parameter_cleanup, param, cleanup=cleanup

	common pc_get_parameter_common, start_par, run_par

	if (keyword_set (cleanup)) then begin
		undefine, start_par
		undefine, run_par
	endif

	return, param
end


; Get a parameter and look for alternatives.
function pc_get_parameter, param, label=label, missing=missing, dim=dim, datadir=datadir, start_param=start_param, run_param=run_param, cleanup=cleanup

	common pc_get_parameter_common, start_par, run_par

	if (keyword_set (start_param)) then start_par = start_param
	if (keyword_set (run_param)) then run_par = run_param
	if (not keyword_set (start_par)) then pc_read_param, obj=start_par, dim=dim, datadir=datadir, /quiet
	if (not keyword_set (run_par)) then pc_read_param, obj=run_par, dim=dim, datadir=datadir, /param2, /quiet
	start_param = start_par
	run_param = run_par

	if (param eq '') then return, pc_get_parameter_cleanup (!Values.D_NaN, cleanup=cleanup)

	; run.in parameters
	run_names = tag_names (run_par)
	if (strcmp (param, 'K_Spitzer', /fold_case)) then begin
		if (any (run_names eq "K_SPITZER")) then return, run_par.K_spitzer
		if (any (run_names eq "KGPARA")) then return, run_par.Kgpara
		if (any (run_names eq "KPARA")) then return, run_par.Kpara
		if (not keyword_set (missing)) then missing = "'K_Spitzer' or 'KPARA'"
	end
	if (strcmp (param, 'K_sat', /fold_case)) then begin
		if (any (run_names eq "KSAT")) then return, run_par.Ksat
		if (not keyword_set (missing)) then missing = "'Ksat'"
	end
	index = where (run_names eq strupcase (param))
	if (index[0] ge 0) then return, pc_get_parameter_cleanup (run_par.(index[0]), cleanup=cleanup)

	; start.in parameters
	start_names = tag_names (start_par)
	index = where (start_names eq strupcase (param))
	if (index[0] ge 0) then return, pc_get_parameter_cleanup (start_par.(index[0]), cleanup=cleanup)

	; Cleanup parameter cache, if requested
	dummy = pc_get_parameter_cleanup ('', cleanup=cleanup)

	; Some additional useful constants
	if (strcmp (param, 'q_electron', /fold_case)) then return, 1.6021766d-19 ; Electron charge [A*s]
	if (strcmp (param, 'm_electron', /fold_case)) then return, 9.109383d-31 ; Electron mass [kg]
	if (strcmp (param, 'm_proton', /fold_case)) then return, 1.6726218d-27 ; Proton mass [kg]
	if (strcmp (param, 'c', /fold_case)) then return, 299792458.d0 ; Speed of light [m/s]
	if (strcmp (param, 'mu0_SI', /fold_case)) then return, 4.0 * double (!Pi) * 1.d-7 ; Magnetic vacuum permeability in SI units [V*s/(A*m)]
	if (strcmp (param, 'pi', /fold_case)) then return, double (!Pi) ; Precise value of pi

	message = "find"
	if (keyword_set (label)) then message = "compute '"+label+"' without"
	if (not keyword_set (missing)) then missing = param
	print, "ERROR: Can't "+message+" parameter "+missing

	return, !Values.D_NaN
end

