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


; Generate parameter abbreviations.
function pc_generate_parameter_abbreviation, param, label=label

	common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
	common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
	common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated

	if (any (strcmp (param, ['mu0_SI', 'mu0_4_pi'], /fold_case))) then begin
		mu0 = pc_get_parameter ('mu0', label=label)
		unit_magnetic = pc_get_parameter ('unit_magnetic', label=label)
		unit_density = pc_get_parameter ('unit_density', label=label)
		unit_velocity = pc_get_parameter ('unit_velocity', label=label)
		return, mu0 * unit_magnetic^2/(unit_density*unit_velocity^2) ; Magnetic vacuum permeability [SI: 4*pi*10^-7 V*s/(A*m)]
	end
	if (any (strcmp (param, ['epsilon_0', 'epsilon0', 'epsilon0_SI'], /fold_case))) then begin
		mu0_SI = pc_get_parameter ('mu0_SI', label=label)
		c = pc_get_parameter ('c', label=label)
		return, 1 / (mu0_SI * c^2) ; Dielectrical vacuum permittivity [SI: A*s/(V*m)]
	end
	if (strcmp (param, 'sigma_total', /fold_case)) then begin
		mu0 = pc_get_parameter ('mu0', label=label)
		eta_total = pc_get_parameter ('eta_total', label=label)
		return, 1 / (mu0 * eta_total) ; Electric conductivity
	end
	if (strcmp (param, 'sigma_SI', /fold_case)) then begin
		mu0_SI = pc_get_parameter ('mu0_SI', label=label)
		eta_SI = pc_get_parameter ('eta_SI', label=label)
		return, 1 / (mu0_SI * eta_SI) ; Electric conductivity [SI: 1/(Ohm*m)]
	end
	if (strcmp (param, 'g_SI', /fold_case)) then begin
		unit_length = pc_get_parameter ('unit_length', label=label)
		unit_time = pc_get_parameter ('unit_time', label=label)
		g = pc_get_parameter ('g_ref', label=label)
		return, g * unit_length/(unit_time^2) ; Gravitational acceleration [SI: m/(s^2)]
	end
	if (strcmp (param, 'g_total', /fold_case)) then begin
		unit_length = pc_get_parameter ('unit_length', label=label)
		unit_time = pc_get_parameter ('unit_time', label=label)
		g_SI = pc_get_parameter ('g_SI', label=label)
		gravx_profile = pc_get_parameter ('gravx_profile', label=label)
		gravy_profile = pc_get_parameter ('gravy_profile', label=label)
		gravz_profile = pc_get_parameter ('gravz_profile', label=label)
		if (size (gravx_profile, /type) ne 7) then return, !Values.D_NaN
		g_total = dblarr (nx,ny,nz,3)
		if (strcmp (gravx_profile, 'const', /fold_case)) then begin
			g_ref = pc_get_parameter ('gravx', label=label)
			g_total[*,*,*,0] += g_ref * unit_length/(unit_time^2)
		end else if ((gravx_profile ne '') and not strcmp (gravx_profile, 'zero', /fold_case)) then begin
			print, "WARNING: The gravity x-profile '"+gravx_profile+"' is not yet implemented."
			stop
		end
		if (strcmp (gravy_profile, 'const', /fold_case)) then begin
			g_ref = pc_get_parameter ('gravy', label=label)
			g_total[*,*,*,1] += g_ref * unit_length/(unit_time^2)
		end else if ((gravy_profile ne '') and not strcmp (gravy_profile, 'zero', /fold_case)) then begin
			print, "WARNING: The gravity y-profile '"+gravy_profile+"' is not yet implemented."
			stop
		end
		if (strcmp (gravz_profile, 'solid_sphere', /fold_case)) then begin
			zref = pc_get_parameter ('zref', label=label) * unit_length
			sphere_rad = pc_get_parameter ('sphere_rad', label=label) * unit_length
			z_prof = z[n1:n2] * unit_length
			g_prof = sign (z_prof - zref + 2*sphere_rad) * sphere_rad^2 / (z_prof - zref + sphere_rad)^2
			inside = where (z_prof - zref + sphere_rad lt sphere_rad, num)
			if (num ge 1) then g_prof[inside] = (z_prof[inside] - zref + sphere_rad) / sphere_rad
			g_total[*,*,*,2] += g_SI * spread (g_prof, [0, 1], [nx, ny])
		end else if (strcmp (gravz_profile, 'const', /fold_case)) then begin
			g_ref = pc_get_parameter ('gravz', label=label)
			g_total[*,*,*,2] += g_ref * unit_length/(unit_time^2)
		end else if ((gravz_profile ne '') and not strcmp (gravz_profile, 'zero', /fold_case)) then begin
			print, "WARNING: The gravity z-profile '"+gravz_profile+"' is not yet implemented."
			stop
		end
		return, g_total
	end
	if (strcmp (param, 'eta_total', /fold_case)) then begin
		resistivities = pc_get_parameter ('iresistivity', label=label)
		eta = pc_get_parameter ('eta', label=label)
		eta_total = 0.0
		eta_found = 0
		for pos = 0, n_elements (resistivities)-1 do begin
			resistivity = resistivities[pos]
			if (any (strcmp (resistivity, ['shock',''], /fold_case))) then continue
			eta_found++
			if (strcmp (resistivity, 'eta-const', /fold_case)) then begin
				eta_total += eta
			end else if (any (strcmp (resistivity, 'eta-zdep', /fold_case))) then begin
				zdep_profile = pc_get_parameter ('zdep_profile', label=label)
				eta_z0 = pc_get_parameter ('eta_z0', label=label)
				eta_jump = pc_get_parameter ('eta_jump', label=label)
				eta_zwidth = pc_get_parameter ('eta_zwidth', label=label)
				if (strcmp (zdep_profile, 'cubic_step', /fold_case)) then begin
					if (eta_zwidth eq 0.0) then eta_zwidth = 5 * mean (1.0 / dz_1)
					eta_total += spread (eta + eta * (eta_jump - 1.) * cubic_step (z[n1:n2], eta_z0, -eta_zwidth), [0, 1], [nx, ny])
				end else begin
					print, "WARNING: The zdep_profile '"+zdep_profile+"' is not yet implemented."
					eta_found--
				end
			end else begin
				print, "WARNING: The resistivity '"+resistivity+"' is not yet implemented."
				eta_found--
			end
		end
		if (eta_found le 0) then eta_total = !Values.D_NaN
		return, eta_total
	end
	if (strcmp (param, 'eta_SI', /fold_case)) then begin
		unit_length = pc_get_parameter ('unit_length', label=label)
		unit_velocity = pc_get_parameter ('unit_velocity', label=label)
		eta_total = pc_get_parameter ('eta_total', label=label)
		return, eta_total * (unit_length*unit_velocity) ; Magnetic resistivity [SI: m^2/s]
	end
	if (strcmp (param, 'cp_SI', /fold_case)) then begin
		unit_velocity = pc_get_parameter ('unit_velocity', label=label)
		unit_temperature = pc_get_parameter ('unit_temperature', label=label)
		cp = pc_get_parameter ('cp', label=label)
		return, cp * unit_velocity^2/unit_temperature ; Specific heat capacity [SI: m^2/(s^2*K)]
	end
	if (strcmp (param, 'degrees_of_freedom', /fold_case)) then begin
		return, 3 ; default: f=3
	end
	if (strcmp (param, 'isentropic_exponent', /fold_case)) then begin
		DOF = pc_get_parameter ('degrees_of_freedom', label=label)
		return, (DOF + 2.0d0) / DOF ; Isentropic exponent (default: 5/3, f=3) [-]
	end
	if (strcmp (param, 'kappa_ideal_3', /fold_case)) then begin
		return, 5/3.0 ; Isentropic exponent for an ideal atomic gas (f=3) [-]
	end
	if (strcmp (param, 'kappa_ideal_5', /fold_case)) then begin
		return, 7/5.0 ; Isentropic exponent for an ideal bi-atomic fixed-molecular gas (f=5) [-]
	end
	if (strcmp (param, 'kappa_ideal_6', /fold_case)) then begin
		return, 8/6.0 ; Isentropic exponent for an ideal tri-atomic fixed-molecular gas (f=6) [-]
	end
	if (strcmp (param, 'kappa_ideal_7', /fold_case)) then begin
		return, 9/7.0 ; Isentropic exponent for an ideal tri-atomic flexible-molecular gas (f=7) [-]
	end

	return, !Values.D_NaN
end


; Get a parameter and look for alternatives.
function pc_get_parameter, param, label=label, missing=missing, dim=dim, datadir=datadir, start_param=start_param, run_param=run_param, cleanup=cleanup

	common pc_get_parameter_common, start_par, run_par

	if (keyword_set (start_param)) then start_par = start_param
	if (keyword_set (run_param)) then run_par = run_param
	if (n_elements (start_par) eq 0) then pc_read_param, obj=start_par, dim=dim, datadir=datadir, /quiet
	if (n_elements (run_par) eq 0) then pc_read_param, obj=run_par, dim=dim, datadir=datadir, /param2, /quiet
	start_param = start_par
	run_param = run_par

	if (param eq '') then return, pc_get_parameter_cleanup (!Values.D_NaN, cleanup=cleanup)

	; run.in parameters
	run_names = tag_names (run_par)
	if (strcmp (param, 'K_Spitzer', /fold_case)) then begin
		if (any (run_names eq "K_SPITZER")) then return, run_par.K_spitzer
		if (any (run_names eq "KPARA")) then return, run_par.Kpara
		if (any (run_names eq "KGPARA")) then return, run_par.Kgpara ; only for backwards compatibility, should be removed some day, because this belongs to the entropy module.
		if (not keyword_set (missing)) then missing = "'K_Spitzer' or 'KPARA'"
	end
	if (strcmp (param, 'K_sat', /fold_case)) then begin
		if (any (run_names eq "KSAT")) then return, run_par.Ksat
		if (not keyword_set (missing)) then missing = "'Ksat'"
	end
	index = where (run_names eq strupcase (param), found)
	if (found ge 1) then return, pc_get_parameter_cleanup (run_par.(index[0]), cleanup=cleanup)

	; start.in parameters
	start_names = tag_names (start_par)
	index = where (start_names eq strupcase (param), found)
	if (found ge 1) then return, pc_get_parameter_cleanup (start_par.(index[0]), cleanup=cleanup)

	; Some additional useful parameter abbreviations
	abbreviation = pc_generate_parameter_abbreviation (param, label=label)

	; Cleanup parameter cache, if requested
	dummy = pc_get_parameter_cleanup ('', cleanup=cleanup)

	; Return parameter abbreviation, if existent
	if (not finite (abbreviation[0], /NaN)) then return, abbreviation

	; Some additional physical constants
	if (strcmp (param, 'AU', /fold_case)) then return, 149597870700.d0 ; Astronomical unit [m]
	if (strcmp (param, 'c', /fold_case)) then return, 299792458.d0 ; Speed of light [m/s]
	if (strcmp (param, 'G', /fold_case)) then return, 6.674d-11 ; Gravitational constant [N*m^2/kg^2]
	if (strcmp (param, 'h_Planck', /fold_case)) then return, 6.62606957d-34 ; Planck constant [J/s]
	if (strcmp (param, 'k_Boltzmann', /fold_case)) then return, 1.3806488d-23 ; Boltzmann constant [J/K]
	if (strcmp (param, 'm_Earth', /fold_case)) then return, 5.974d24 ; Earth mass [kg]
	if (strcmp (param, 'm_electron', /fold_case)) then return, 9.109383d-31 ; Electron mass [kg]
	if (strcmp (param, 'm_neutron', /fold_case)) then return, 1.6749274d-27 ; Neutron mass [kg]
	if (strcmp (param, 'm_proton', /fold_case)) then return, 1.6726218d-27 ; Proton mass [kg]
	if (strcmp (param, 'm_Sun', /fold_case)) then return, 1.989d30 ; Sun mass [kg]
	if (strcmp (param, 'N_Avogadro', /fold_case)) then return, 6.02214129d23 ; Avogadro number [1/mol]
	if (strcmp (param, 'pi', /fold_case)) then return, !DPi ; Precise value of pi
	if (strcmp (param, 'q_electron', /fold_case)) then return, 1.6021766d-19 ; Electron charge [A*s]
	if (strcmp (param, 'R', /fold_case)) then return, 8.314462145d0 ; universal gas constant [J/K/mol] = k_Boltzman * N_Avogadro
	if (strcmp (param, 'R_Earth', /fold_case)) then return, 6367456d0 ; Earth average radius [m]
	if (strcmp (param, 'R_Sun', /fold_case)) then return, 696342.d3 ; Sun radius [m]
	if (strcmp (param, 'u', /fold_case)) then return, 1.660538921d-27 ; Atomic mass unit [kg]

	; Some additional units
	if (strcmp (param, 'unit_time', /fold_case)) then return, pc_get_parameter ('unit_length', label=label) / pc_get_parameter ('unit_velocity', label=label)
	if (strcmp (param, 'unit_energy', /fold_case)) then return, pc_get_parameter ('unit_density', label=label) * pc_get_parameter ('unit_velocity', label=label)^2 / pc_get_parameter ('unit_length', label=label)^3

	; Some additional mathematical constants
	if (strcmp (param, 'e', /fold_case)) then return, 2.718281828459045235d0 ; Euler constnat

	; Non-existent parameter
	message = "find"
	if (keyword_set (label)) then message = "compute '"+label+"' without"
	if (not keyword_set (missing)) then missing = param
	print, "ERROR: Can't "+message+" parameter "+missing

	return, !Values.D_NaN
end

