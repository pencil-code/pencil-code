;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_get_quantity.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Calculation of physical quantities. If the unit-structure is given,
;   this unit system is chosen, otherwise the result is in Pencil-units.
;
;  Parameters:
;   * quantity       Name of physical quantity to be computed.
;                    (If an array is given, a structure is returned.)
;   * vars           Data array or data structure as load by pc_read_*.
;   * index          Indices or tags of the given variables inside data.
;   * /ghost         If set, ghost cells are included in the returned data.
;   * /cache         If set, a chache is used to optimize computation.
;   * /cleanup       If set, the chache is being freed after computation.
;
;   Label            Description
;  ===============================================================
;   Temp             temperature
;   u_abs            absolute velocity
;   u_z              velocity z-component
;   rho              density
;   log_rho          decatic logarithm of density
;   ln_rho           natural logarithm of density
;   n_rho            particle density
;   P_therm          thermal pressure
;   HR_ohm           volumetric Ohmic heating rate
;   HR_viscous       volumetric viscous heating rate
;   B_z              magnetic field z-component
;   Spitzer_q        absolute value of Spitzer heat flux vector
;   HR_ohm           volumetric Ohmic heating rate
;   j_abs            current density
;   species_5        number density of 5th chemical species
;   [...]            more are listed in "pc_check_quantities.pro":
;                    IDL> help, pc_check_quantities (/all), /str
;
;  Examples: (in ascending order of efficiency)
;  ============================================
;
;  * Using 'pc_read_var': (NOT RECOMMENDED)
;
;   Load varfile and calculate separate quantities, simplest version:
;   IDL> pc_read_var, obj=vars
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', vars)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', vars)
;   IDL> B_z = pc_get_quantity ('B_z', vars)
;   IDL> tvscl, HR_viscous[*,*,20]
;
;   Load varfile and calculate an array of quantities:
;   IDL> pc_read_var, obj=vars
;   IDL> result = pc_get_quantity (['HR_viscous', 'HR_ohm', 'B_z'], vars)
;   IDL> tvscl, result.HR_viscous[*,*,20]
;
;  * Using 'pc_read_var_raw': (HIGHLY RECOMMENDED)
;
;   Load varfile and calculate one quantity only:
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', 'var.dat')
;
;   Load varfile and calculate separate quantities, using a data array:
;   IDL> pc_read_var_raw, obj=var, tags=tags
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', var, tags)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', var, tags)
;   IDL> B_z = pc_get_quantity ('B_z', var, tags)
;   IDL> tvscl, HR_viscous[*,*,20]
;
;   Load varfile and calculate an array of quantities: (RECOMMENDED FOR CLI)
;   IDL> pc_read_var_raw, obj=var, tags=tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param
;   IDL> result = pc_get_quantity (['HR_viscous', 'HR_ohm', 'B_z'], var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param)
;   IDL> tvscl, result.HR_viscous[*,*,20]
;
;   Load varfile and separately calculate quantities, using the cache manually: (RECOMMENDED FOR SCRIPTS)
;   IDL> pc_read_var_raw, obj=var, tags=tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
;   IDL> B_z = pc_get_quantity ('B_z', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache, /cleanup)
;   IDL> tvscl, HR_viscous[*,*,20]
;
;  * Using 'pc_read_slice_raw': (RECOMMENDED FOR 2D CUTS)
;
;   Load 2D-slice and calculate separate quantities, not using cache:
;   IDL> pc_read_slice_raw, obj=slice, tags=tags, cut_z=20, slice_dim=dim
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', slice, tags, dim=dim)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', slice, tags, dim=dim)
;   IDL> B_z = pc_get_quantity ('B_z', slice, tags, dim=dim)
;   IDL> tvscl, HR_viscous
;
;  * Using 'pc_read_subvol_raw': (RECOMMENDED FOR 3D SUBVOLUMES)
;
;   Load 3D-subvolume from 'VAR123' and calculate separate quantities, not using the cache:
;   IDL> pc_read_subvol_raw, obj=var, varfile='VAR123', tags=tags, xs=16, xe=47, ys=16, ye=47, zs=0, ze=31, /addghosts, sub_dim=dim, sub_grid=grid
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', var, tags, dim=dim, grid=grid)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', var, tags, dim=dim, grid=grid)
;   IDL> B_z = pc_get_quantity ('B_z', var, tags, dim=dim, grid=grid)
;   IDL> tvscl, HR_viscous[*,*,20]
;
;   Load 3D-subvolume from 'var.dat' and calculate separate quantities, using the cache manually:
;   IDL> pc_read_subvol_raw, obj=var, tags=tags, xs=16, xe=47, ys=16, ye=47, zs=0, ze=31, /addghosts, sub_dim=dim, sub_grid=grid, start_param=start_param, run_param=run_param
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
;   IDL> B_z = pc_get_quantity ('B_z', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache, /cleanup)
;   IDL> tvscl, HR_viscous[*,*,20]
;


;  Known issues:
;  =============
;  none


; Computation of physical quantities.
; =============================================================================
; PLEASE ADD MORE PHYSICAL QUANTITIES IN THIS FUNCTION.
; And update the availability and dependency list in "pc_check_quantities.pro".
; =============================================================================
function pc_compute_quantity, vars, index, quantity, ghost=ghost

	common quantity_cache, uu, rho, grad_rho, n_rho, Temp, grad_Temp, P_therm, grad_P_therm, bb, B_2, jj, EMF, Poynting, Poynting_j, Poynting_u, F_Lorentz
	common quantity_params, sources, l1, l2, m1, m2, n1, n2, nx, ny, nz, unit, start_par, run_par, alias
	common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
	common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated

	default, ghost, 0
	gl1 = l1
	gl2 = l2
	gm1 = m1
	gm2 = m2
	gn1 = n1
	gn2 = n2
	gnx = nx
	gny = ny
	gnz = nz
	if (keyword_set (ghost)) then begin
		gl1 = 0
		gl2 += nghostx
		gm1 = 0
		gm2 += nghosty
		gn1 = 0
		gn2 += nghostz
		gnx += 2 * nghostx
		gny += 2 * nghosty
		gnz += 2 * nghostz
	end

	if (strcmp (quantity, 'u', /fold_case)) then begin
		; Velocity [m / s]
		if (n_elements (uu) eq 0) then uu = vars[gl1:gl2,gm1:gm2,gn1:gn2,index.uu] * unit.velocity
		return, uu
	end
	if (strcmp (quantity, 'u_x', /fold_case)) then begin
		; Velocity x-component
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u', ghost=ghost)
		return, uu[*,*,*,0]
	end
	if (strcmp (quantity, 'u_y', /fold_case)) then begin
		; Velocity y-component
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u', ghost=ghost)
		return, uu[*,*,*,1]
	end
	if (strcmp (quantity, 'u_z', /fold_case)) then begin
		; Velocity z-component
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u', ghost=ghost)
		return, uu[*,*,*,2]
	end
	if (strcmp (quantity, 'u_abs', /fold_case)) then begin
		; Absolute value of the velocity
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u', ghost=ghost)
		if (n_elements (u_abs) eq 0) then u_abs = sqrt (dot2 (uu))
		return, u_abs
	end
	if (strcmp (quantity, 'grad_u', /fold_case)) then begin
		; Velocity gradient
		if (n_elements (grad_u) eq 0) then grad_u = (grad (sqrt (dot2 (vars[*,*,*,index.uu]))))[l1:l2,m1:m2,n1:n2,*] * unit.velocity
		return, grad_u
	end
	if (strcmp (quantity, 'grad_u_abs', /fold_case)) then begin
		; Velocity gradient absolute value
		if (n_elements (grad_u) eq 0) then grad_u = pc_compute_quantity (vars, index, 'grad_u', ghost=ghost)
		return, sqrt (dot2 (grad_u))
	end

	if (strcmp (quantity, 'E_therm', /fold_case)) then begin
		; Thermal energy in one unit volume [J] = specific thermal energy [J/m^3]
		gamma = pc_get_parameter ('isentropic_exponent', label=quantity)
		cp_SI = pc_get_parameter ('cp_SI', label=quantity)
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp', ghost=ghost)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		return, rho * Temp * cp_SI / gamma
	end
	if (strcmp (quantity, 'E_kin', /fold_case)) then begin
		; Kinetic energy in one unit volume [J] = kinetic energy density [J/m^3]
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u', ghost=ghost)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		return, 0.5 * dot2 (uu)^2 * rho
	end

	if (strcmp (quantity, 'Temp', /fold_case)) then begin
		; Temperature [K]
		if (n_elements (Temp) eq 0) then begin
			if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
				Temp = exp (vars[gl1:gl2,gm1:gm2,gn1:gn2,index.lnTT]) * unit.temperature
			end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
				Temp = vars[gl1:gl2,gm1:gm2,gn1:gn2,index.TT] * unit.temperature
			end else if (any (strcmp (sources, 'ss', /fold_case))) then begin
				cp = pc_get_parameter ('cp', label=quantity)
				cp_SI = pc_get_parameter ('cp_SI', label=quantity)
				cs0 = pc_get_parameter ('cs0', label=quantity)
				gamma = pc_get_parameter ('gamma', label=quantity)
				if (gamma eq 1.0) then tmp = 1.0 else tmp = gamma - 1.0
				ln_Temp_0 = alog (cs0^2 / cp / tmp) + alog (unit.temperature)
				ln_rho_0 = alog (pc_get_parameter ('rho0', label=quantity)) + alog (unit.density)
				S = pc_compute_quantity (vars, index, 'S', ghost=ghost)
				ln_rho = pc_compute_quantity (vars, index, 'ln_rho', ghost=ghost) - ln_rho_0
				Temp = exp (ln_Temp_0 + gamma/cp_SI * S + (gamma-1) * ln_rho)
			end
		end
		return, Temp
	end
	if (strcmp (quantity, 'grad_Temp', /fold_case)) then begin
		; Gradient of temperature
		if (n_elements (grad_Temp) eq 0) then begin
			if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
				grad_Temp = (grad (exp (vars[*,*,*,index.lnTT])))[l1:l2,m1:m2,n1:n2,*] * unit.temperature / unit.length
			end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
				grad_Temp = (grad (vars[*,*,*,index.TT]))[l1:l2,m1:m2,n1:n2,*] * unit.temperature / unit.length
			end else if (any (strcmp (sources, 'ss', /fold_case))) then begin
				cp = pc_get_parameter ('cp', label=quantity)
				cp_SI = pc_get_parameter ('cp_SI', label=quantity)
				cs0 = pc_get_parameter ('cs0', label=quantity)
				gamma = pc_get_parameter ('gamma', label=quantity)
				if (gamma eq 1.0) then tmp = 1.0 else tmp = gamma - 1.0
				ln_Temp_0 = alog (cs0^2 / cp / tmp) + alog (unit.temperature)
				ln_rho_0 = alog (pc_get_parameter ('rho0', label=quantity)) + alog (unit.density)
				S = pc_compute_quantity (vars, index, 'S', /ghost)
				ln_rho = pc_compute_quantity (vars, index, 'ln_rho', /ghost) - ln_rho_0
				Temp_with_ghosts = exp (ln_Temp_0 + gamma/cp_SI * S + (gamma-1) * ln_rho)
				; *** WORK HERE: this is maybe better computed with a direct derivation from the entropy S
				grad_Temp = (grad (Temp_with_ghosts))[l1:l2,m1:m2,n1:n2,*] / unit.length
			end
		end
		return, grad_Temp
	end
	if (strcmp (quantity, 'grad_Temp_abs', /fold_case)) then begin
		; Absolute value of temperature gradient
		if (n_elements (grad_Temp) eq 0) then grad_Temp = pc_compute_quantity (vars, index, 'grad_Temp')
		return, sqrt (dot2 (grad_Temp))
	end
	if (strcmp (quantity, 'log_Temp', /fold_case)) then begin
		; Logarithmic temperature
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.lnTT] / alog (10.0) + alog10 (unit.temperature)
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			return, alog10 (vars[gl1:gl2,gm1:gm2,gn1:gn2,index.TT]) + alog10 (unit.temperature)
		end else if (any (strcmp (sources, 'ss', /fold_case))) then begin
			cp = pc_get_parameter ('cp', label=quantity)
			cp_SI = pc_get_parameter ('cp_SI', label=quantity)
			cs0 = pc_get_parameter ('cs0', label=quantity)
			gamma = pc_get_parameter ('gamma', label=quantity)
			if (gamma eq 1.0) then tmp = 1.0 else tmp = gamma - 1.0
			ln_Temp_0 = alog (cs0^2 / cp / tmp) + alog (unit.temperature)
			ln_rho_0 = alog (pc_get_parameter ('rho0', label=quantity)) + alog (unit.density)
			S = pc_compute_quantity (vars, index, 'S', ghost=ghost)
			ln_rho = pc_compute_quantity (vars, index, 'ln_rho', ghost=ghost) - ln_rho_0
			return, ln_Temp_0 + gamma/cp_SI * S + (gamma-1) * ln_rho + alog10 (unit.temperature)
		end
	end
	if (strcmp (quantity, 'ln_Temp', /fold_case)) then begin
		; Natural logarithmic temperature
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.lnTT] + alog (unit.temperature)
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			return, alog (vars[gl1:gl2,gm1:gm2,gn1:gn2,index.TT]) + alog (unit.temperature)
		end else if (any (strcmp (sources, 'ss', /fold_case))) then begin
			cp = pc_get_parameter ('cp', label=quantity)
			cp_SI = pc_get_parameter ('cp_SI', label=quantity)
			cs0 = pc_get_parameter ('cs0', label=quantity)
			gamma = pc_get_parameter ('gamma', label=quantity)
			if (gamma eq 1.0) then tmp = 1.0 else tmp = gamma - 1.0
			ln_Temp_0 = alog (cs0^2 / cp / tmp) + alog (unit.temperature)
			ln_rho_0 = alog (pc_get_parameter ('rho0', label=quantity)) + alog (unit.density)
			S = pc_compute_quantity (vars, index, 'S', ghost=ghost)
			ln_rho = pc_compute_quantity (vars, index, 'ln_rho', ghost=ghost) - ln_rho_0
			return, ln_Temp_0 + gamma/cp_SI * S + (gamma-1) * ln_rho + alog (unit.temperature)
		end
	end

	if (strcmp (quantity, 'S', /fold_case)) then begin
		; Entropy
		if (any (strcmp (sources, 'ss', /fold_case))) then begin
			return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.ss] * (unit.velocity^2 * unit.mass / unit.temperature)
		end else begin
			cp = pc_get_parameter ('cp', label=quantity)
			cs0 = pc_get_parameter ('cs0', label=quantity)
			gamma = pc_get_parameter ('gamma', label=quantity)
			if (gamma eq 1.0) then tmp = 1.0 else tmp = gamma - 1.0
			ln_Temp_0 = alog (cs0^2 / cp / tmp) + alog (unit.temperature)
			ln_rho_0 = alog (pc_get_parameter ('rho0', label=quantity)) + alog (unit.density)
			ln_Temp = pc_compute_quantity (vars, index, 'ln_Temp', ghost=ghost) - ln_Temp_0
			ln_rho = pc_compute_quantity (vars, index, 'ln_rho', ghost=ghost) - ln_rho_0
			return, cp/gamma * (ln_Temp - (gamma-1) * ln_rho)
		end
	end

	if (strcmp (quantity, 'q_abs', /fold_case)) then begin
		; Absolute value of the heat flux density vector q [W / m^2] = [kg / s^3]
		chi = pc_get_parameter ('chi', label=quantity)
		return, chi * sqrt (dot2 (pc_compute_quantity (vars, index, 'grad_Temp'))) * unit.density / (unit.length^2 * unit.time^3 * unit.temperature)
	end
	if (strcmp (quantity, 'q_sat', /fold_case)) then begin
		; Absolute value of the saturation heat flux density vector q [W / m^2] = [kg / s^3]
		K_sat = pc_get_parameter ('K_sat', label=quantity)
		if (K_sat le 0.0) then K_sat = 1.0
		mu = pc_get_parameter ('mu', label=quantity)
		m_e = pc_get_parameter ('m_electron', label=quantity)
		m_p = pc_get_parameter ('m_proton', label=quantity)
		k_B = pc_get_parameter ('k_Boltzmann', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp', ghost=ghost)
		return, K_sat * 1.5 * sqrt (3 / m_e) / (m_e + m_p) * k_B^1.5 * mu * rho * Temp^1.5
		; This should correspond to the calculation in the solar_corona module, but is untested:
		; return, K_sat * sqrt (Temp / dot2 (pc_compute_quantity (vars, index, 'grad_Temp'))) * (7.28e7 * unit.density * unit.velocity^3 / unit.length * sqrt (unit.temperature))
	end
	if (strcmp (quantity, 'Spitzer_K_parallel', /fold_case)) then begin
		; Field-aligned Spitzer heat flux coefficient, not including T^2.5 [W / (m * K^3.5)] = [kg * m / (s^3 * K^3.5)]
		K_Spitzer = pc_get_parameter ('K_Spitzer', label=quantity)
		; TODO: check units consistency
		return, K_Spitzer * (unit.density * unit.velocity^3 * unit.length / unit.temperature)
	end
	if (strcmp (quantity, 'Spitzer_q', /fold_case)) then begin
		; Absolute value of the Spitzer heat flux density vector q [W / m^2] = [kg / s^3]
		Spitzer_K_parallel = pc_compute_quantity (vars, index, 'Spitzer_K_parallel')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		return, Spitzer_K_parallel * Temp^2.5 * sqrt (dot2 (pc_compute_quantity (vars, index, 'grad_Temp')))
	end
	if (strcmp (quantity, 'Spitzer_q_parallel', /fold_case)) then begin
		; Absolute value of the field-aligned Spitzer heat flux density vector q_parallel [W / m^2] = [kg / s^3]
		Spitzer_K_parallel = pc_compute_quantity (vars, index, 'Spitzer_K_parallel')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		return, Spitzer_K_parallel * Temp^2.5 * dot (pc_compute_quantity (vars, index, 'grad_Temp'), bb) / sqrt (B_2)
	end
	if (strcmp (quantity, 'Spitzer_q_perpendicular', /fold_case)) then begin
		; Absolute value of the field-perpendicular Spitzer heat flux density vector q_perpendicular [W / m^2] = [kg / s^3]
		return, pc_compute_quantity (vars, index, 'Spitzer_q_parallel') * pc_compute_quantity (vars, index, 'Spitzer_ratio')
	end
	if (strcmp (quantity, 'Spitzer_dt', /fold_case)) then begin
		; Spitzer heat flux timestep [s]
		cp_SI = pc_get_parameter ('cp_SI', label=quantity)
		gamma = pc_get_parameter ('gamma', label=quantity)
		cdtv = pc_get_parameter ('cdtv', label=quantity)
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		Spitzer_K_parallel = pc_compute_quantity (vars, index, 'Spitzer_K_parallel')
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			grad_ln_Temp = (grad (vars[*,*,*,index.lnTT]))[l1:l2,m1:m2,n1:n2,*]
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			grad_ln_Temp = (grad (alog (vars[*,*,*,index.TT])))[l1:l2,m1:m2,n1:n2,*]
		end
		dt = cdtv / (gamma * cp_SI * Spitzer_K_parallel) * rho * B_2 * sqrt (dot2 (grad_ln_Temp)) / (Temp^2.5 * abs (dot (bb, grad_ln_Temp)))
		; Iterate through the z-direction
		dx_inv = pc_compute_quantity (vars, index, 'inv_dx')
		dy_inv = pc_compute_quantity (vars, index, 'inv_dy')
		dz_inv = pc_compute_quantity (vars, index, 'inv_dz')
		if (all (lequidist[0:1])) then begin
			dxy_inv = dx_inv[0]^2 + dy_inv[0]^2
		end else begin
			dxy_inv = spread (dx_inv^2, 1, ny) + spread (dy_inv^2, 0, nx)
		end
		for pz = 0, nz - 1 do dt[*,*,pz] /= dxy_inv + dz_inv[pz]^2
		return, dt
	end
	if (strcmp (quantity, 'Spitzer_ratio', /fold_case)) then begin
		; Ratio of perpendicular to parallel Spitzer heat conduction coefficients
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		if (n_elements (n_rho) eq 0) then n_rho = pc_compute_quantity (vars, index, 'n_rho')
		return, 2.e-31 * n_rho^2 / (Temp^3 * B_2) ; [Solar MHD, E. Priest (1982/1984), p. 86]
	end
	if (strcmp (quantity, 'Spitzer_q_ratio', /fold_case)) then begin
		; Ratio of saturation heat flux to Spitzer heat flux
		return, pc_compute_quantity (vars, index, 'q_sat') / pc_compute_quantity (vars, index, 'Spitzer_q')
	end
	if (strcmp (quantity, 'Spitzer_collision_frequency_e', /fold_case)) then begin
		; Spitzer electron collision frequency [1 / s]
		mu = pc_get_parameter ('mu', label=quantity)
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		if (n_elements (n_rho) eq 0) then n_rho = pc_compute_quantity (vars, index, 'n_rho')
		return, (mu / 266d3) * n_rho * Temp^(-1.5) * pc_compute_quantity (vars, index, 'Coulomb_logarithm')
	end
	if (strcmp (quantity, 'Spitzer_conductivity', /fold_case)) then begin
		; Spitzer conductivity in a two-component WKB plasma [A / (V*m)]
		m_electron = pc_get_parameter ('m_electron', label=quantity)
		q_electron = pc_get_parameter ('q_electron', label=quantity)
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp', ghost=ghost)
		return, (266d3 * q_electron^2 / m_electron) * Temp^1.5 / pc_compute_quantity (vars, index, 'Coulomb_logarithm', ghost=ghost)
	end
	if (strcmp (quantity, 'Spitzer_mag_diffusivity', /fold_case)) then begin
		; Spitzer magnetic diffusivity [m^2 / s]
		mu0_SI_inv = 1 / pc_get_parameter ('mu0_SI', label=quantity)
		return, mu0_SI_inv / pc_compute_quantity (vars, index, 'Spitzer_conductivity', ghost=ghost)
	end
	if (strcmp (quantity, 'Coulomb_logarithm', /fold_case)) then begin
		Spitzer_K_parallel = pc_compute_quantity (vars, index, 'Spitzer_K_parallel', ghost=ghost)
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp', ghost=ghost)
		return, (1.8d-10 / Spitzer_K_parallel) * Temp^2.5 ; Coulomb logarithm [Solar MHD, E. Priest (1982/1984), p. 79, 86, eq. 2.34]
	end

	if (strcmp (quantity, 'collision_frequency_e', /fold_case)) then begin
		; Electron collision frequency [1 / s]
		c = pc_get_parameter ('c', label=quantity)
		m_electron = pc_get_parameter ('m_electron', label=quantity)
		q_electron = pc_get_parameter ('q_electron', label=quantity)
		return, (q_electron / (m_electron * c)) * pc_compute_quantity (vars, index, 'B_abs')
	end
	if (strcmp (quantity, 'WKB_conductivity', /fold_case)) then begin
		; Electrical conductivity in a two-component WKB plasma [A / (V*m)]
		c = pc_get_parameter ('c', label=quantity)
		q_electron = pc_get_parameter ('q_electron', label=quantity)
		if (n_elements (n_rho) eq 0) then n_rho = pc_compute_quantity (vars, index, 'n_rho')
		return, (q_electron * c) * n_rho / pc_compute_quantity (vars, index, 'B_abs')
	end
	if (strcmp (quantity, 'WKB_mag_diffusivity', /fold_case)) then begin
		; Magnetic diffusivity [m^2 / s]
		mu0_SI_inv = 1 / pc_get_parameter ('mu0_SI', label=quantity)
		return, mu0_SI_inv / pc_compute_quantity (vars, index, 'WKB_conductivity')
	end

	if (strcmp (quantity, 'rho', /fold_case)) then begin
		; Density [kg / m^3]
		if (n_elements (rho) eq 0) then begin
			if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
				rho = exp (vars[gl1:gl2,gm1:gm2,gn1:gn2,index.lnrho]) * unit.density
			end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
				rho = vars[gl1:gl2,gm1:gm2,gn1:gn2,index.rho] * unit.density
			end
		end
		return, rho
	end
	if (strcmp (quantity, 'grad_rho', /fold_case)) then begin
		; Gradient of density [kg / m^4]
		if (n_elements (grad_rho) eq 0) then begin
			if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
				grad_rho = (grad (exp (vars[*,*,*,index.lnrho])))[l1:l2,m1:m2,n1:n2,*] * unit.density / unit.length
			end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
				grad_rho = (grad (vars[*,*,*,index.rho]))[l1:l2,m1:m2,n1:n2,*] * unit.density / unit.length
			end
		end
		return, grad_rho
	end
	if (strcmp (quantity, 'log_rho', /fold_case)) then begin
		; Logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.lnrho] / alog (10.0) + alog10 (unit.density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			return, alog10 (vars[gl1:gl2,gm1:gm2,gn1:gn2,index.rho]) + alog10 (unit.density)
		end
	end
	if (strcmp (quantity, 'ln_rho', /fold_case)) then begin
		; Natural logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.lnrho] + alog (unit.density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			return, alog (vars[gl1:gl2,gm1:gm2,gn1:gn2,index.rho]) + alog (unit.density)
		end
	end
	if (strcmp (quantity, 'n_rho', /fold_case)) then begin
		; Particle density [1 / m^3]
		mu = pc_get_parameter ('mu', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		if (n_elements (n_rho) eq 0) then n_rho = rho / (pc_get_parameter ('m_proton', label=quantity) * mu)
		return, n_rho
	end

	if (strcmp (quantity, 'P_therm', /fold_case)) then begin
		; Thermal pressure [N / m^2]
		cp_SI = pc_get_parameter ('cp_SI', label=quantity)
		gamma = pc_get_parameter ('gamma', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp', ghost=ghost)
		if (n_elements (P_therm) eq 0) then P_therm = cp_SI * (gamma - 1.0) / gamma * rho * Temp
		return, P_therm
	end
	if (strcmp (quantity, 'grad_P_therm', /fold_case)) then begin
		; Gradient of thermal pressure [N / m^3]
		cp_SI = pc_get_parameter ('cp_SI', label=quantity)
		gamma = pc_get_parameter ('gamma', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		if (n_elements (grad_rho) eq 0) then grad_rho = pc_compute_quantity (vars, index, 'grad_rho')
		if (n_elements (grad_Temp) eq 0) then grad_Temp = pc_compute_quantity (vars, index, 'grad_Temp')
		if (n_elements (grad_P_therm) eq 0) then begin
			fact = cp_SI * (gamma - 1.0) / gamma
			grad_P_therm = grad_rho
			for pa = 0, 2 do grad_P_therm[*,*,*,pa] = fact * (grad_rho[*,*,*,pa] * Temp + rho * grad_Temp[*,*,*,pa])
		end
		return, grad_P_therm
	end
	if (strcmp (quantity, 'grad_P_therm_abs', /fold_case)) then begin
		; Absolute value of thermal pressure gradient [N / m^3]
		if (n_elements (grad_P_therm) eq 0) then grad_P_therm = pc_compute_quantity (vars, index, 'grad_P_therm')
		return, sqrt (dot2 (grad_P_therm))
	end
	if (strcmp (quantity, 'c_sound', /fold_case)) then begin
		; Sound speed [m / s]
		kappa = pc_get_parameter ('isentropic_exponent', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		if (n_elements (P_therm) eq 0) then P_therm = pc_compute_quantity (vars, index, 'P_therm', ghost=ghost)
		return, sqrt (kappa * P_therm / rho)
	end
	if (strcmp (quantity, 'H_P_therm_x', /fold_case)) then begin
		; Scaling height of thermal pressure x-component [m]
		if (n_elements (P_therm) eq 0) then P_therm = pc_compute_quantity (vars, index, 'P_therm')
		if (n_elements (grad_P_therm) eq 0) then grad_P_therm = pc_compute_quantity (vars, index, 'grad_P_therm')
		dP_therm_dx = grad_P_therm[*,*,*,0]
		return, -(P_therm / dP_therm_dx)
	end
	if (strcmp (quantity, 'H_P_therm_y', /fold_case)) then begin
		; Scaling height of thermal pressure y-component [m]
		if (n_elements (P_therm) eq 0) then P_therm = pc_compute_quantity (vars, index, 'P_therm')
		if (n_elements (grad_P_therm) eq 0) then grad_P_therm = pc_compute_quantity (vars, index, 'grad_P_therm')
		dP_therm_dy = grad_P_therm[*,*,*,1]
		return, -(P_therm / dP_therm_dy)
	end
	if (strcmp (quantity, 'H_P_therm_z', /fold_case)) then begin
		; Scaling height of thermal pressure z-component [m]
		if (n_elements (P_therm) eq 0) then P_therm = pc_compute_quantity (vars, index, 'P_therm')
		if (n_elements (grad_P_therm) eq 0) then grad_P_therm = pc_compute_quantity (vars, index, 'grad_P_therm')
		dP_therm_dz = grad_P_therm[*,*,*,2]
		return, -(P_therm / dP_therm_dz)
	end

	if (strcmp (quantity, 'a_grav', /fold_case)) then begin
		; Gravity acceleration vector
		return, pc_get_parameter ('g_total', label=quantity)
	end
	if (strcmp (quantity, 'a_grav_abs', /fold_case)) then begin
		; Gravity acceleration absolute value
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'a_grav')))
	end
	if (strcmp (quantity, 'a_grav_x', /fold_case)) then begin
		; Gravity acceleration x-component
		return, (pc_compute_quantity (vars, index, 'a_grav'))[*,*,*,0]
	end
	if (strcmp (quantity, 'a_grav_y', /fold_case)) then begin
		; Gravity acceleration y-component
		return, (pc_compute_quantity (vars, index, 'a_grav'))[*,*,*,1]
	end
	if (strcmp (quantity, 'a_grav_z', /fold_case)) then begin
		; Gravity acceleration z-component
		return, (pc_compute_quantity (vars, index, 'a_grav'))[*,*,*,2]
	end
	if (strcmp (quantity, 'F_grav', /fold_case)) then begin
		; Gravity force vector
		a_grav = pc_compute_quantity (vars, index, 'a_grav')
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		return, a_grav * spread_scalar_to_components (rho)
	end
	if (strcmp (quantity, 'F_grav_abs', /fold_case)) then begin
		; Gravity force absolute value
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'F_grav')))
	end
	if (strcmp (quantity, 'F_grav_x', /fold_case)) then begin
		; Gravity force x-component
		return, (pc_compute_quantity (vars, index, 'F_grav'))[*,*,*,0]
	end
	if (strcmp (quantity, 'F_grav_y', /fold_case)) then begin
		; Gravity force y-component
		return, (pc_compute_quantity (vars, index, 'F_grav'))[*,*,*,1]
	end
	if (strcmp (quantity, 'F_grav_z', /fold_case)) then begin
		; Gravity force z-component
		return, (pc_compute_quantity (vars, index, 'F_grav'))[*,*,*,2]
	end
	if (strcmp (quantity, 'rho_u', /fold_case)) then begin
		; Impulse density
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		return, spread_scalar_to_components (rho) * pc_compute_quantity (vars, index, 'u', ghost=ghost)
	end
	if (strcmp (quantity, 'rho_u_abs', /fold_case)) then begin
		; Impulse density absolute value
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'rho_u', ghost=ghost)))
	end
	if (strcmp (quantity, 'rho_u_x', /fold_case)) then begin
		; Impulse density x-component
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		return, rho * pc_compute_quantity (vars, index, 'u_x', ghost=ghost)
	end
	if (strcmp (quantity, 'rho_u_y', /fold_case)) then begin
		; Impulse density y-component
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		return, rho * pc_compute_quantity (vars, index, 'u_y', ghost=ghost)
	end
	if (strcmp (quantity, 'rho_u_z', /fold_case)) then begin
		; Impulse density z-component
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho', ghost=ghost)
		return, rho * pc_compute_quantity (vars, index, 'u_z', ghost=ghost)
	end

	if (strcmp (quantity, 'c_Alfven', /fold_case)) then begin
		; Alfvén velocity
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		B_abs = pc_compute_quantity (vars, index, 'B_abs')
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		return, B_abs / sqrt (mu0_SI * rho)
	end
	if (strcmp (quantity, 'c_Alfven_inv', /fold_case)) then begin
		; Inverse of the Alfvén velocity
		return, 1.0 / pc_compute_quantity (vars, index, 'c_Alfven')
	end
	if (strcmp (quantity, 'rho_c', /fold_case)) then begin
		; Minimum density for an Alfvén speed below the speed of light
		cdtv = pc_get_parameter ('cdtv', label=quantity)
		c = pc_get_parameter ('c', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		return, B_2 / (2 * mu0_SI * (c * cdtv)^2)
	end
	if (strcmp (quantity, 'rho_c_ratio', /fold_case)) then begin
		; Ratio of density to minimum density for an Alfvén speed below the speed of light
		return, pc_compute_quantity (vars, index, 'rho') / pc_compute_quantity (vars, index, 'rho_c')
	end

	if (strcmp (quantity, 'HR_viscous', /fold_case)) then begin
		; Viscous heating rate [W / m^3] = [kg / m^3] * [m / s]^3 / [m]
		nu = pc_get_parameter ('nu', label=quantity)
		u_xx = (xder (vars[*,*,*,index.ux]))[l1:l2,m1:m2,n1:n2]
		u_xy = (yder (vars[*,*,*,index.ux]))[l1:l2,m1:m2,n1:n2]
		u_xz = (zder (vars[*,*,*,index.ux]))[l1:l2,m1:m2,n1:n2]
		u_yx = (xder (vars[*,*,*,index.uy]))[l1:l2,m1:m2,n1:n2]
		u_yy = (yder (vars[*,*,*,index.uy]))[l1:l2,m1:m2,n1:n2]
		u_yz = (zder (vars[*,*,*,index.uy]))[l1:l2,m1:m2,n1:n2]
		u_zx = (xder (vars[*,*,*,index.uz]))[l1:l2,m1:m2,n1:n2]
		u_zy = (yder (vars[*,*,*,index.uz]))[l1:l2,m1:m2,n1:n2]
		u_zz = (zder (vars[*,*,*,index.uz]))[l1:l2,m1:m2,n1:n2]
		div_u3 = (u_xx + u_yy + u_zz) / 3.0
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		return, nu * rho * ( 2*((u_xx - div_u3)^2 + (u_yy - div_u3)^2 + (u_zz - div_u3)^2) + (u_xy + u_yx)^2 + (u_xz + u_zx)^2 + (u_yz + u_zy)^2 ) * unit.velocity^3 / unit.length
	end
	if (strcmp (quantity, 'HR_viscous_particle', /fold_case)) then begin
		; Viscous heating rate per particle [W]
		if (n_elements (n_rho) eq 0) then n_rho = pc_compute_quantity (vars, index, 'n_rho')
		return, pc_compute_quantity (vars, index, 'HR_viscous') / n_rho
	end
	if (strcmp (quantity, 'Rn_viscous', /fold_case)) then begin
		; Viscous mesh Reynolds number
		nu = pc_get_parameter ('nu', label=quantity)
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		Rx = spread (pc_compute_quantity (vars, index, 'dx'), [1,2], [ny,nz]) * abs (uu[*,*,*,0])
		Ry = spread (pc_compute_quantity (vars, index, 'dy'), [0,2], [nx,nz]) * abs (uu[*,*,*,1])
		Rz = spread (pc_compute_quantity (vars, index, 'dz'), [0,1], [nx,ny]) * abs (uu[*,*,*,2])
		return, ((Rx > Ry) > Rz) / (nu * unit.length^2/unit.time)
	end
	if (strcmp (quantity, 'Heat_cool_compression', /fold_case)) then begin
		; Heating/cooling due to compression  ;[J/m^3]
		gamma = pc_get_parameter ('isentropic_exponent', label=quantity)
		cp_SI = pc_get_parameter ('cp_SI', label=quantity)
		lnTT  = pc_compute_quantity (vars, index, 'ln_Temp')
		gamma_m1 = gamma-1
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (any (strcmp (sources, 'uu', /fold_case))) then begin
			divu = (div(vars[*,*,*,index.uu]))
			divu = divu  * unit.velocity  / unit.length
		end
		tmp = -(gamma_m1 * divu[l1:l2,m1:m2,n1:n2])
		return,(rho * (exp(lnTT + tmp) - exp(lnTT)) * cp_SI)/gamma
	end
	if (strcmp (quantity, 'Heat_cool_visc', /fold_case)) then begin
		; Heating/cooling due to viscous heating  ;[J/m^3]
		nu    = (pc_get_parameter ('nu', label=quantity))*(unit.length^2/unit.time)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		cp_SI = pc_get_parameter ('cp_SI', label=quantity)
		gamma = pc_get_parameter ('isentropic_exponent', label=quantity)
		lnTT  = pc_compute_quantity (vars, index, 'ln_Temp')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		if (any (strcmp (sources, 'uu', /fold_case))) then begin
			uij_PC = (gij(vars[*,*,*,index.uu]))
			uij = (uij_PC) * unit.velocity/unit.length
			divu_PC = (div(vars[*,*,*,index.uu]))
			divu = (divu_PC) * unit.velocity/unit.length
		end
		oo = make_array ([ (size(uij, /dim))[0:2], 3], type=size (uij, /type))
		oo[*,*,*,0]=uij[*,*,*,2,1]-uij[*,*,*,1,2]
		oo[*,*,*,1]=uij[*,*,*,0,2]-uij[*,*,*,2,0]
		oo[*,*,*,2]=uij[*,*,*,1,0]-uij[*,*,*,0,1]
		o2   = oo[*,*,*,0]^2+oo[*,*,*,1]^2+oo[*,*,*,2]^2
		visc_heat = nu * o2[l1:l2,m1:m2,n1:n2]

		sij = make_array (size(uij, /dim), type=size (uij, /type))
		for j=0,2 do begin
			sij[*,*,*,j,j]=uij[*,*,*,j,j]-(1./3.)*divu
			for i=j+1,2 do begin
				sij[*,*,*,i,j]=.5*(uij[*,*,*,i,j]+uij[*,*,*,j,i])
				sij[*,*,*,j,i]=sij[*,*,*,i,j]
			endfor
		endfor
		sij2 = sij[*,*,*,1,1]^2
		for i = 1, 2 do begin
			sij2 = sij2 + sij[*,*,*,i,i]^2
			for j = 0, i-1 do begin
				sij2 = sij2 + 2 * sij[*,*,*,i,j]^2
			endfor
		endfor

		tmp = (gamma/cp_SI) * (1/Temp) * (visc_heat + 2 * nu * sij2[l1:l2,m1:m2,n1:n2])
		return, (rho * (exp(lnTT + tmp) - exp(lnTT)) * cp_SI)/gamma
	end
	if (any (strcmp (quantity, ['A', 'A_contour'], /fold_case))) then begin
		; Magnetic vector potential [T * m]
		return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.aa] * (unit.magnetic_field*unit.length)
	end
	if (strcmp (quantity, 'A_abs', /fold_case)) then begin
		; Magnetic vector potential [T * m]
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'A')))
	end
	if (strcmp (quantity, 'A_x', /fold_case)) then begin
		; Magnetic vector potential x-component
		return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.ax] * (unit.magnetic_field*unit.length)
	end
	if (strcmp (quantity, 'A_y', /fold_case)) then begin
		; Magnetic vector potential y-component
		return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.ay] * (unit.magnetic_field*unit.length)
	end
	if (strcmp (quantity, 'A_z', /fold_case)) then begin
		; Magnetic vector potential z-component
		return, vars[gl1:gl2,gm1:gm2,gn1:gn2,index.az] * (unit.magnetic_field*unit.length)
	end
	if (strcmp (quantity, 'div_A', /fold_case)) then begin
		; Divergence of the magnetic vector potential [T]
		return, (div (vars[*,*,*,index.aa]))[l1:l2,m1:m2,n1:n2] * unit.magnetic_field
	end

	if (strcmp (quantity, 'B', /fold_case)) then begin
		; Magnetic field vector [Tesla]
		if (any (strcmp (sources, 'bb', /fold_case))) then begin
			if (n_elements (bb) eq 0) then bb = vars[gl1:gl2,gm1:gm2,gn1:gn2,index.bb] * unit.magnetic_field
		end else begin
			if (n_elements (bb) eq 0) then bb = (curl (vars[*,*,*,index.aa]))[l1:l2,m1:m2,n1:n2,*] * unit.magnetic_field
		end
		return, bb
	end
	if (strcmp (quantity, 'B_abs', /fold_case)) then begin
		; Magnetic field strengh [T]
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		return, sqrt (B_2)
	end
	if (strcmp (quantity, 'B_2', /fold_case)) then begin
		; Magnetic field strengh squared [T^2]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (B_2) eq 0) then B_2 = dot2 (bb)
		return, B_2
	end
	if (strcmp (quantity, 'B_x', /fold_case)) then begin
		; Magnetic field x-component [T]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, bb[*,*,*,0]
	end
	if (strcmp (quantity, 'B_y', /fold_case)) then begin
		; Magnetic field y-component [T]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, bb[*,*,*,1]
	end
	if (strcmp (quantity, 'B_z', /fold_case)) then begin
		; Magnetic field z-component [T]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, bb[*,*,*,2]
	end
	if (strcmp (quantity, 'dB_dx', /fold_case)) then begin
		; Magnetic field x-derivative [T / m]
		dB_dx = dblarr (nx, ny, nz, 3, /nozero)
		dB_dx[*,*,*,*] = unit.magnetic_field / unit.length
		dB_dx[*,*,*,0] *= (xderyder (vars[*,*,*,index.az]) - xderzder (vars[*,*,*,index.ay]))[l1:l2,m1:m2,n1:n2]
		dB_dx[*,*,*,1] *= (xder2    (vars[*,*,*,index.az]) - xderzder (vars[*,*,*,index.ax]))[l1:l2,m1:m2,n1:n2]
		dB_dx[*,*,*,2] *= (xder2    (vars[*,*,*,index.ay]) - xderyder (vars[*,*,*,index.ax]))[l1:l2,m1:m2,n1:n2]
		return, dB_dx
	end
	if (strcmp (quantity, 'dB_dy', /fold_case)) then begin
		; Magnetic field y-derivative [T / m]
		dB_dy = dblarr (nx, ny, nz, 3, /nozero)
		dB_dy[*,*,*,*] = unit.magnetic_field / unit.length
		dB_dy[*,*,*,0] *= (yder2    (vars[*,*,*,index.az]) - yderzder (vars[*,*,*,index.ay]))[l1:l2,m1:m2,n1:n2]
		dB_dy[*,*,*,1] *= (yderxder (vars[*,*,*,index.az]) - yderzder (vars[*,*,*,index.ax]))[l1:l2,m1:m2,n1:n2]
		dB_dy[*,*,*,2] *= (yderxder (vars[*,*,*,index.ay]) - yder2    (vars[*,*,*,index.ax]))[l1:l2,m1:m2,n1:n2]
		return, dB_dy
	end
	if (strcmp (quantity, 'dB_dz', /fold_case)) then begin
		; Magnetic field z-derivative [T / m]
		dB_dz = dblarr (nx, ny, nz, 3, /nozero)
		dB_dz[*,*,*,*] = unit.magnetic_field / unit.length
		dB_dz[*,*,*,0] *= (zderyder (vars[*,*,*,index.az]) - zder2    (vars[*,*,*,index.ay]))[l1:l2,m1:m2,n1:n2]
		dB_dz[*,*,*,1] *= (zderxder (vars[*,*,*,index.az]) - zder2    (vars[*,*,*,index.ax]))[l1:l2,m1:m2,n1:n2]
		dB_dz[*,*,*,2] *= (zderxder (vars[*,*,*,index.ay]) - zderyder (vars[*,*,*,index.ax]))[l1:l2,m1:m2,n1:n2]
		return, dB_dz
	end
	if (strcmp (quantity, 'grad_B', /fold_case)) then begin
		; Magnetic field gradient [T / m]
		if (n_elements (grad_B) eq 0) then grad_B = (gradcurl (pc_compute_quantity (vars, index, 'A', /ghost)))[l1:l2,m1:m2,n1:n2,*] / unit.length^2
		return, grad_B
	end
	if (strcmp (quantity, 'grad_B_abs', /fold_case)) then begin
		; Magnetic field gradient absolute value
		if (n_elements (grad_B) eq 0) then grad_B = pc_compute_quantity (vars, index, 'grad_B')
		return, sqrt (dot2 (grad_B))
	end

	if (strcmp (quantity, 'E', /fold_case)) then begin
		; Electric field [V / m]
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		if (n_elements (EMF) eq 0) then EMF = -cross (uu, bb)
		sigma_SI_inv = 1.0 / pc_get_parameter ('sigma_SI', label=quantity)
		E = EMF
		for pa = 0, 2 do E[*,*,*,pa] += sigma_SI_inv * jj[*,*,*,pa]
		return, E
	end
	if (strcmp (quantity, 'E_abs', /fold_case)) then begin
		; Electric field strengh [V / m]
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'E')))
	end
	if (strcmp (quantity, 'EMF', /fold_case)) then begin
		; Electro motive force [V / m]
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (EMF) eq 0) then EMF = -cross (uu, bb)
		return, EMF
	end
	if (strcmp (quantity, 'EMF_abs', /fold_case)) then begin
		; Electro motive force strength [V / m]
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'EMF')))
	end
	if (strcmp (quantity, 'EMF_x', /fold_case)) then begin
		; Electro motive force x [V / m]
		return, (pc_compute_quantity (vars, index, 'EMF'))[*,*,*,0]
	end
	if (strcmp (quantity, 'EMF_y', /fold_case)) then begin
		; Electro motive force y [V / m]
		return, (pc_compute_quantity (vars, index, 'EMF'))[*,*,*,1]
	end
	if (strcmp (quantity, 'EMF_z', /fold_case)) then begin
		; Electro motive force z [V / m]
		return, (pc_compute_quantity (vars, index, 'EMF'))[*,*,*,2]
	end
	if (strcmp (quantity, 'E_j', /fold_case)) then begin
		; Current electric field [V / m]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		sigma_SI_inv = 1.0 / pc_get_parameter ('sigma_SI', label=quantity)
		E = jj
		for pa = 0, 2 do E[*,*,*,pa] *= sigma_SI_inv
		return, E
	end
	if (strcmp (quantity, 'E_j_abs', /fold_case)) then begin
		; Current electric field strength [V / m]
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'E_j')))
	end
	if (strcmp (quantity, 'E_j_x', /fold_case)) then begin
		; Current electric field x [V / m]
		return, (pc_compute_quantity (vars, index, 'E_j'))[*,*,*,0]
	end
	if (strcmp (quantity, 'E_j_y', /fold_case)) then begin
		; Current electric field y [V / m]
		return, (pc_compute_quantity (vars, index, 'E_j'))[*,*,*,1]
	end
	if (strcmp (quantity, 'E_j_z', /fold_case)) then begin
		; Current electric field z [V / m]
		return, (pc_compute_quantity (vars, index, 'E_j'))[*,*,*,2]
	end
	if (strcmp (quantity, 'E_x', /fold_case)) then begin
		; Electric field x-component [V / m]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		sigma_SI_inv = 1.0 / pc_get_parameter ('sigma_SI', label=quantity)
		return, sigma_SI_inv * jj[*,*,*,0]
	end
	if (strcmp (quantity, 'E_y', /fold_case)) then begin
		; Electric field y-component [V / m]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		sigma_SI_inv = 1.0 / pc_get_parameter ('sigma_SI', label=quantity)
		return, sigma_SI_inv * jj[*,*,*,1]
	end
	if (strcmp (quantity, 'E_z', /fold_case)) then begin
		; Electric field z-component [V / m]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		sigma_SI_inv = 1.0 / pc_get_parameter ('sigma_SI', label=quantity)
		return, sigma_SI_inv * jj[*,*,*,2]
	end
	if (strcmp (quantity, 'E_parallel', /fold_case)) then begin
		; Electric field strengh parallel to the magnetic field [V / m]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		return, dot (pc_compute_quantity (vars, index, 'E'), bb) / sqrt (B_2)
	end
	if (strcmp (quantity, 'E_perpendicular', /fold_case)) then begin
		; Electric field strengh perpendicular to the magnetic field [V / m]
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'E')) - pc_compute_quantity (vars, index, 'E_parallel')^2)
	end
	if (strcmp (quantity, 'grad_E_abs', /fold_case)) then begin
		; Gradient of electric field strenght [V / m^2]
		return, !Values.D_NaN ; not yet implemented...

		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		c = 1.0 / pc_get_parameter ('c', label=quantity)
		sigma_SI_inv = 1.0 / pc_get_parameter ('sigma_SI', label=quantity)
		E = -(1.0/c) * cross (uu, bb)
		for pa = 0, 2 do E[*,*,*,pa] += sigma_SI_inv * jj[*,*,*,pa]
		return, E

		return, grad (pc_compute_quantity (vars, index, 'E_abs'))
	end
	if (strcmp (quantity, 'grad_E_abs_abs', /fold_case)) then begin
		; Absolute value of electric field strength gradient [V / m^2]
		return, sqrt (dot2 (pc_compute_quantity (vars, index, 'grad_E_abs')))
	end
	if (strcmp (quantity, 'beta', /fold_case)) then begin
		; Plasma beta
		if (n_elements (P_therm) eq 0) then P_therm = pc_compute_quantity (vars, index, 'P_therm')
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		return, 2 * mu0_SI * P_therm / B_2
	end
	if (strcmp (quantity, 'rho_mag', /fold_case)) then begin
		; Magnetic energy density [J / m^3]
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		return, B_2 / (2.0 * mu0_SI)
	end
	if (strcmp (quantity, 'Rn_mag', /fold_case)) then begin
		; Magnetic mesh Reynolds number of velocities perpendicular to the magnetic field
		eta = pc_get_parameter ('eta_total', label=quantity)
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (B_2) eq 0) then B_2 = pc_compute_quantity (vars, index, 'B_2')
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		if (n_elements (u_abs) eq 0) then u_abs = pc_compute_quantity (vars, index, 'u_abs')
		fact = 1.0 / (u_abs * B_2)
		BxUxB = cross (cross (bb, uu), bb)
		Rx = spread (pc_compute_quantity (vars, index, 'dx'), [1,2], [ny,nz]) * fact * uu[*,*,*,0] * BxUxB[*,*,*,0]
		Ry = spread (pc_compute_quantity (vars, index, 'dy'), [0,2], [nx,nz]) * fact * uu[*,*,*,1] * BxUxB[*,*,*,1]
		Rz = spread (pc_compute_quantity (vars, index, 'dz'), [0,1], [nx,ny]) * fact * uu[*,*,*,2] * BxUxB[*,*,*,2]
		return, ((Rx > Ry) > Rz) / (eta * unit.length^2/unit.time)
	end

	if (strcmp (quantity, 'forcing', /fold_case)) then begin
		; Forcing function [kg * m / s^2]
		if (n_elements (ff) eq 0) then ff = (vars[*,*,*,index.fx:index.fz])[gl1:gl2,gm1:gm2,gn1:gn2,*] * unit.mass * unit.length / unit.time^2
		return, ff
	end
	if (strcmp (quantity, 'forcing_abs', /fold_case)) then begin
		; Absolute value of the forcing function [kg * m / s^2]
		if (n_elements (ff) eq 0) then ff = pc_compute_quantity (vars, index, 'forcing', ghost=ghost)
		return, sqrt (dot2 (ff))
	end
	if (strcmp (quantity, 'forcing_x', /fold_case)) then begin
		; Forcing function x-component [kg * m / s^2]
		if (n_elements (ff) eq 0) then ff = pc_compute_quantity (vars, index, 'forcing', ghost=ghost)
		return, ff[*,*,*,0]
	end
	if (strcmp (quantity, 'forcing_y', /fold_case)) then begin
		; Forcing function y-component [kg * m / s^2]
		if (n_elements (ff) eq 0) then ff = pc_compute_quantity (vars, index, 'forcing', ghost=ghost)
		return, ff[*,*,*,1]
	end
	if (strcmp (quantity, 'forcing_z', /fold_case)) then begin
		; Forcing function z-component [kg * m / s^2]
		if (n_elements (ff) eq 0) then ff = pc_compute_quantity (vars, index, 'forcing', ghost=ghost)
		return, ff[*,*,*,2]
	end

	if (strcmp (quantity, 'j', /fold_case)) then begin
		; Current density [A / m^2]
		if (n_elements (jj) eq 0) then jj = (curlcurl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] * unit.current_density
		return, jj
	end
	if (strcmp (quantity, 'j_abs', /fold_case)) then begin
		; Absolute value of the current density [A / m^2]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, sqrt (dot2 (jj))
	end
	if (strcmp (quantity, 'j_x', /fold_case)) then begin
		; Current density x-component [A / m^2]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, jj[*,*,*,0]
	end
	if (strcmp (quantity, 'j_y', /fold_case)) then begin
		; Current density y-component [A / m^2]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, jj[*,*,*,1]
	end
	if (strcmp (quantity, 'j_z', /fold_case)) then begin
		; Current density z-component [A / m^2]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, jj[*,*,*,2]
	end
	if (strcmp (quantity, 'eta_j', /fold_case)) then begin
		; Current density times eta [A / s]
		eta = pc_get_parameter ('eta_total', label=quantity) * unit.length*unit.velocity
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		if (size (eta, /n_dimensions) eq 0) then return, eta * jj
		eta_j = jj
		eta_j[*,*,*,0] *= eta
		eta_j[*,*,*,1] *= eta
		eta_j[*,*,*,2] *= eta
		return, eta_j
	end

	if (strcmp (quantity, 'H_mag', /fold_case)) then begin
		; Magnetic field helicity density [T^2 * m]
		aa = pc_compute_quantity (vars, index, 'A')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, dot (aa, bb)
	end
	if (strcmp (quantity, 'H_mag_pos', /fold_case)) then begin
		; Magnetic field helicity density (positive part) [T^2 * m]
		H_mag_pos = pc_compute_quantity (vars, index, 'H_mag') > 0.0
		return, H_mag_pos
	end
	if (strcmp (quantity, 'H_mag_neg', /fold_case)) then begin
		; Magnetic field helicity density (negative part) [T^2 * m]
		H_mag_neg = (-pc_compute_quantity (vars, index, 'H_mag')) > 0.0
		return, H_mag_neg
	end
	if (strcmp (quantity, 'H_j', /fold_case)) then begin
		; Electric current helicity [A * T * m]
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, mu0_SI * dot (jj, bb)
	end
	if (strcmp (quantity, 'dH_mag_dt', /fold_case)) then begin
		; Change rate of magnetic field helicity density [T^2 * m / s =?= A * T / s]
		eta = pc_get_parameter ('eta_total', label=quantity) * unit.length*unit.velocity
		H_j = pc_compute_quantity (vars, index, 'H_j')
		return, -2 * eta * H_j
	end

	if (strcmp (quantity, 'F_Lorentz', /fold_case)) then begin
		; Lorentz force [N]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		if (n_elements (F_Lorentz) eq 0) then F_Lorentz = cross (jj, bb)
		return, F_Lorentz
	end
	if (strcmp (quantity, 'F_Lorentz_x', /fold_case)) then begin
		; X-component of the Lorentz force [N]
		if (n_elements (F_Lorentz) eq 0) then F_Lorentz = pc_compute_quantity (vars, index, 'F_Lorentz')
		return, F_Lorentz[*,*,*,0]
	end
	if (strcmp (quantity, 'F_Lorentz_y', /fold_case)) then begin
		; Y-component of the Lorentz force [N]
		if (n_elements (F_Lorentz) eq 0) then F_Lorentz = pc_compute_quantity (vars, index, 'F_Lorentz')
		return, F_Lorentz[*,*,*,1]
	end
	if (strcmp (quantity, 'F_Lorentz_z', /fold_case)) then begin
		; Z-component of the Lorentz force [N]
		if (n_elements (F_Lorentz) eq 0) then F_Lorentz = pc_compute_quantity (vars, index, 'F_Lorentz')
		return, F_Lorentz[*,*,*,2]
	end
	if (strcmp (quantity, 'F_Lorentz_abs', /fold_case)) then begin
		; Absolute value of the Lorentz force [N]
		if (n_elements (F_Lorentz) eq 0) then F_Lorentz = pc_compute_quantity (vars, index, 'F_Lorentz')
		return, sqrt (dot2 (F_Lorentz))
	end
	if (strcmp (quantity, 'W_Lorentz', /fold_case)) then begin
		; Work done by the Lorentz force [J]
		if (n_elements (F_Lorentz) eq 0) then F_Lorentz = pc_compute_quantity (vars, index, 'F_Lorentz')
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		return, dot (uu, F_Lorentz)
	end
	if (strcmp (quantity, 'Lorentz_angle', /fold_case)) then begin
		; Angle between the current density and the magnetic field vectors [°]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		B_abs = pc_compute_quantity (vars, index, 'B_abs')
		j_abs = pc_compute_quantity (vars, index, 'j_abs')
		return, acos (dot (jj, bb) / (j_abs * B_abs)) * (180 / !DPi)
	end
	if (strcmp (quantity, 'Lorentz_angle_deviation', /fold_case)) then begin
		; Deviation of the angle (j,B) from 0° or 180° with values in [-90°,90°] where negative is anti-parallel
		Lorentz_angle = pc_compute_quantity (vars, index, 'Lorentz_angle')
		return, ((Lorentz_angle + 90) mod 180) - 90
	end

	if (strcmp (quantity, 'HR_ohm', /fold_case)) then begin
		; Ohming heating rate [W / m^3] = [kg/m^3 * (m/s)^3 / m]
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		eta_SI = pc_get_parameter ('eta_SI', label=quantity)
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, eta_SI * mu0_SI * dot2 (jj)
	end
	if (strcmp (quantity, 'HR_ohm_particle', /fold_case)) then begin
		; Ohming heating rate per particle [W]
		if (n_elements (n_rho) eq 0) then n_rho = pc_compute_quantity (vars, index, 'n_rho')
		return, pc_compute_quantity (vars, index, 'HR_ohm') / n_rho
	end

	if (strcmp (quantity, 'Poynting_j', /fold_case)) then begin
		; current Poynting flux vector [W / m^2]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (Poynting_j) eq 0) then Poynting_j = cross (pc_compute_quantity (vars, index, 'eta_j'), bb)
		return, Poynting_j
	end
	if (strcmp (quantity, 'Poynting_u', /fold_case)) then begin
		; velocity Poynting flux vector [W / m^2]
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (Poynting_u) eq 0) then Poynting_u = cross (cross (uu, bb), bb) / (-mu0_SI)
		return, Poynting_u
	end
	if (strcmp (quantity, 'Poynting', /fold_case)) then begin
		; Poynting flux vector [W / m^2]
		if (n_elements (Poynting) eq 0) then begin
			if ((n_elements (Poynting_j) gt 0) and (n_elements (Poynting_u) gt 0)) then begin
				Poynting = Poynting_j + Poynting_u
			end else begin
				mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
				if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
				if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
				eta_j = pc_compute_quantity (vars, index, 'eta_j')
				Poynting = cross ((eta_j - cross (uu, bb)/mu0_SI), bb)
			end
		end
		return, Poynting
	end
	if (strcmp (quantity, 'Poynting_j_x', /fold_case)) then begin
		; X-component of the current Poynting flux vector [W / m^2]
		if (n_elements (Poynting_j) eq 0) then Poynting_j = pc_compute_quantity (vars, index, 'Poynting_j')
		return, Poynting_j[*,*,*,0]
	end
	if (strcmp (quantity, 'Poynting_j_y', /fold_case)) then begin
		; Y-component of the current Poynting flux vector [W / m^2]
		if (n_elements (Poynting_j) eq 0) then Poynting_j = pc_compute_quantity (vars, index, 'Poynting_j')
		return, Poynting_j[*,*,*,1]
	end
	if (strcmp (quantity, 'Poynting_j_z', /fold_case)) then begin
		; Z-component of the current Poynting flux vector [W / m^2]
		if (n_elements (Poynting_j) eq 0) then Poynting_j = pc_compute_quantity (vars, index, 'Poynting_j')
		return, Poynting_j[*,*,*,2]
	end
	if (strcmp (quantity, 'Poynting_j_abs', /fold_case)) then begin
		; Absolute value of the current Poynting flux [W / m^2]
		if (n_elements (Poynting_j) eq 0) then Poynting_j = pc_compute_quantity (vars, index, 'Poynting_j')
		return, sqrt (dot2 (Poynting_j))
	end
	if (strcmp (quantity, 'Poynting_u_x', /fold_case)) then begin
		; X-component of the velocity Poynting flux vector [W / m^2]
		if (n_elements (Poynting_u) eq 0) then Poynting_u = pc_compute_quantity (vars, index, 'Poynting_u')
		return, Poynting_u[*,*,*,0]
	end
	if (strcmp (quantity, 'Poynting_u_y', /fold_case)) then begin
		; Y-component of the velocity Poynting flux vector [W / m^2]
		if (n_elements (Poynting_u) eq 0) then Poynting_u = pc_compute_quantity (vars, index, 'Poynting_u')
		return, Poynting_u[*,*,*,1]
	end
	if (strcmp (quantity, 'Poynting_u_z', /fold_case)) then begin
		; Z-component of the velocity Poynting flux vector [W / m^2]
		if (n_elements (Poynting_u) eq 0) then Poynting_u = pc_compute_quantity (vars, index, 'Poynting_u')
		return, Poynting_u[*,*,*,2]
	end
	if (strcmp (quantity, 'Poynting_u_abs', /fold_case)) then begin
		; Absolute value of the velocity Poynting flux [W / m^2]
		if (n_elements (Poynting_u) eq 0) then Poynting_u = pc_compute_quantity (vars, index, 'Poynting_u')
		return, sqrt (dot2 (Poynting_u))
	end
	if (strcmp (quantity, 'Poynting_x', /fold_case)) then begin
		; X-component of the Poynting flux vector [W / m^2]
		if (n_elements (Poynting) eq 0) then Poynting = pc_compute_quantity (vars, index, 'Poynting')
		return, Poynting[*,*,*,0]
	end
	if (strcmp (quantity, 'Poynting_y', /fold_case)) then begin
		; Y-component of the Poynting flux vector [W / m^2]
		if (n_elements (Poynting) eq 0) then Poynting = pc_compute_quantity (vars, index, 'Poynting')
		return, Poynting[*,*,*,1]
	end
	if (strcmp (quantity, 'Poynting_z', /fold_case)) then begin
		; Z-component of the Poynting flux vector [W / m^2]
		if (n_elements (Poynting) eq 0) then Poynting = pc_compute_quantity (vars, index, 'Poynting')
		return, Poynting[*,*,*,2]
	end
	if (strcmp (quantity, 'Poynting_abs', /fold_case)) then begin
		; Absolute value of the Poynting flux [W / m^2]
		if (n_elements (Poynting) eq 0) then Poynting = pc_compute_quantity (vars, index, 'Poynting')
		return, sqrt (dot2 (Poynting))
	end

	species = (stregex (quantity, '^species_([0-9]+)$', /subexpr, /extract, /fold_case))[1]
	if (species ne '') then begin
		result = execute ('species_index = index.yy'+species)
		if (not result) then begin
			print, "ERROR: Unknown species 'yy"+species+"'"
			return, !Values.D_NaN
		end
		return, vars[gl1:gl2,gm1:gm2,gn1:gn2,species_index] * 100
	end

	; Check for Pencil Code alias names
	if (n_elements (alias) eq 0) then alias = pc_check_quantities (sources=sources, /aliases)
	pos = find_tag (alias, quantity)
	if (pos ge 0) then return, pc_compute_quantity (vars, index, alias.(pos))

	; Timestamp
	if (strcmp (quantity, 'time', /fold_case)) then return, index.time * unit.time

	; Coordinates
	if (strcmp (quantity, 'x', /fold_case)) then return, x[gl1:gl2] * unit.length
	if (strcmp (quantity, 'y', /fold_case)) then return, y[gm1:gm2] * unit.length
	if (strcmp (quantity, 'z', /fold_case)) then return, z[gn1:gn2] * unit.length

	; Grid distances
	if (strcmp (quantity, 'dx', /fold_case)) then return, 1.0 / dx_1[gl1:gl2] * unit.length
	if (strcmp (quantity, 'dy', /fold_case)) then return, 1.0 / dy_1[gm1:gm2] * unit.length
	if (strcmp (quantity, 'dz', /fold_case)) then return, 1.0 / dz_1[gn1:gn2] * unit.length

	; Inverse grid distances
	if (strcmp (quantity, 'inv_dx', /fold_case)) then return, dx_1[gl1:gl2] / unit.length
	if (strcmp (quantity, 'inv_dy', /fold_case)) then return, dy_1[gm1:gm2] / unit.length
	if (strcmp (quantity, 'inv_dz', /fold_case)) then return, dz_1[gn1:gn2] / unit.length

	; Grid volume
	if (strcmp (quantity, 'dV', /fold_case)) then begin
		dx = pc_compute_quantity (vars, index, 'dx')
		dy = pc_compute_quantity (vars, index, 'dy')
		dz = pc_compute_quantity (vars, index, 'dz')
		if (all (lequidist[0:1])) then begin
			dV = dx[0] * dy[0] * dz[0]
		end else begin
			dV = spread (dx, [1,2], [gny,gnz]) * spread (dy, [0,2], [gnx,gnz]) * spread (dz, [0,1], [gnx,gny])
		end
		return, dV
	end

	; Box size
	if (strcmp (quantity, 'size_x', /fold_case)) then return, (x[l2]-x[l1] + lperi[0] * mean (1.0 / dx_1[[l1,l2]])) * unit.length
	if (strcmp (quantity, 'size_y', /fold_case)) then return, (y[m2]-y[m1] + lperi[1] * mean (1.0 / dy_1[[m1,m2]])) * unit.length
	if (strcmp (quantity, 'size_z', /fold_case)) then return, (z[n2]-z[n1] + lperi[2] * mean (1.0 / dz_1[[n1,n2]])) * unit.length

	; Box origin
	if (strcmp (quantity, 'origin_x', /fold_case)) then return, (x[l1] - lperi[0] * 0.5 / dx_1[l1]) * unit.length
	if (strcmp (quantity, 'origin_y', /fold_case)) then return, (y[m1] - lperi[1] * 0.5 / dy_1[m1]) * unit.length
	if (strcmp (quantity, 'origin_z', /fold_case)) then return, (z[n1] - lperi[2] * 0.5 / dz_1[n1]) * unit.length

	print, "ERROR: Unknown quantity '"+quantity+"'"
	return, !Values.D_NaN
end


; Clean up cache for computation of physical quantities.
pro pc_quantity_cache_cleanup

	common quantity_cache, uu, rho, grad_rho, n_rho, Temp, grad_Temp, P_therm, grad_P_therm, bb, B_2, jj, EMF, Poynting, Poynting_j, Poynting_u, F_Lorentz
	common quantity_params, sources, l1, l2, m1, m2, n1, n2, nx, ny, nz, unit, start_par, run_par, alias

	undefine, uu
	undefine, rho
	undefine, grad_rho
	undefine, n_rho
	undefine, Temp
	undefine, grad_Temp
	undefine, P_therm
	undefine, grad_P_therm
	undefine, bb
	undefine, B_2
	undefine, jj
	undefine, EMF
	undefine, Poynting
	undefine, Poynting_j
	undefine, Poynting_u
	undefine, F_Lorentz

	undefine, sources

	undefine, l1
	undefine, l2
	undefine, m1
	undefine, m2
	undefine, n1
	undefine, n2
	undefine, nx
	undefine, ny
	undefine, nz

	undefine, unit
	undefine, start_par
	undefine, run_par
end


; Calculation of physical quantities.
function pc_get_quantity, quantity, vars, index, varfile=varfile, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, datadir=datadir, ghost=ghost, cache=cache, cleanup=cleanup, verbose=verbose

	common quantity_params, sources, l1, l2, m1, m2, n1, n2, nx, ny, nz, unit, start_par, run_par, alias
	common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
	common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated

	if (keyword_set (ghost) and keyword_set (cache)) then message, "pc_get_quantity: keywords 'ghost' and 'cache' can not be used together."
	if (keyword_set (verbose)) then quiet = 0 else quiet = 1
	if (keyword_set (cleanup) and not keyword_set (cache)) then begin
		pc_quantity_cache_cleanup
		if (size (quantity, /type) eq 0) then return, !Values.D_NaN
	end

	if (n_elements (quantity) eq 0) then quantity = ""
	if (not any (quantity ne "") or ((n_elements (vars) eq 0) and (n_elements (varfile) eq 0)) or ((n_elements (index) eq 0) and (size (vars, /type) ne 8) and (size (vars, /type) ne 7) and (size (varfile, /type) ne 7))) then begin
		; Print usage
		print, "USAGE:"
		print, "======"
		print, "* using var-files:"
		print, "-----------------------"
		print, "HR = pc_get_quantity ('HR_viscous', 'var.dat')"
		print, ""
		print, "* using var-structures:"
		print, "-----------------------"
		print, "pc_read_var, obj=vars"
		print, "HR = pc_get_quantity ('HR_viscous', vars)"
		print, ""
		print, "* using var-arrays:"
		print, "-------------------"
		print, "pc_read_var_raw, obj=var, tags=tags"
		print, "HR = pc_get_quantity ('HR_viscous', var, tags)"
		print, ""
		print, "* using 2D-slices:"
		print, "------------------"
		print, "pc_read_slice_raw, obj=var, tags=tags, cut_x=20"
		print, "HR = pc_get_quantity ('HR_viscous', var, tags)"
		print, ""
		print, "* to get a list of available quantities:"
		print, "----------------------------------------"
		print, "help, pc_check_quantities (/all), /str"
		print, ""
		if (not any (quantity ne "")) then print, "ERROR: no quantity selected"
		if (n_elements (vars) eq 0) then print, "ERROR: no data source given"
		if (n_elements (index) eq 0) then print, "ERROR: data source has no associated index structure"
		return, -1
	end

	if (size (vars, /type) eq 7) then varfile = vars
	if (size (varfile, /type) eq 7) then begin
		pc_read_var_raw, obj=vars, varfile=varfile, tags=index, datadir=datadir, dim=dim, grid=grid, start_param=start_param, run_param=run_param, quiet=quiet
	end

	; Setup 'start.in' and 'run.in' parameters
	dummy = pc_get_parameter ('', start_param=start_param, run_param=run_param, dim=dim, datadir=datadir)
	start_par = start_param
	run_par = run_param

	; Set default units
	if (n_elements (units) eq 0) then pc_units, obj=units, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	unit = units

	if (size (vars, /type) eq 8) then begin
		; Need to have a valid varcontent
		if (size (index, /type) eq 8) then varcontent = index
		if (n_elements (varcontent) eq 0) then varcontent = pc_varcontent (datadir=datadir, dim=dim, param=start_param)
		; Create array out of given structure and pass recursively computed results
		array = pc_convert_vars_struct (vars, varcontent, index)
		return, pc_get_quantity (quantity, array, index, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, datadir=datadir, cache=cache, cleanup=cleanup)
	end

	sources = tag_names (index)

	if (n_elements (dim) eq 0) then begin
		; Check consistency of dimensions
		if (((size (vars))[1] ne mx) or ((size (vars))[2] ne my) or ((size (vars))[3] ne mz)) then begin
			print, "ERROR: Data doesn't fit to the loaded dim structure, please pass the corresponding dim structure as parameter."
			return, -1
		end
		pc_read_dim, obj=glob_dim, datadir=datadir, quiet=quiet
		l1 = glob_dim.nghostx
		l2 = mx - 1 - glob_dim.nghostx
		m1 = glob_dim.nghosty
		m2 = my - 1 - glob_dim.nghosty
		n1 = glob_dim.nghostz
		n2 = mz - 1 - glob_dim.nghostz
		nx = mx - 2*glob_dim.nghostx
		ny = my - 2*glob_dim.nghosty
		nz = mz - 2*glob_dim.nghostz
	end else begin
		; Set dimensions in common block for derivative routines
		mx = dim.mx
		my = dim.my
		mz = dim.mz
		nw = dim.nx * dim.ny * dim.nz
		l1 = dim.l1
		l2 = dim.l2
		m1 = dim.m1
		m2 = dim.m2
		n1 = dim.n1
		n2 = dim.n2
		nx = mx - 2*dim.nghostx
		ny = my - 2*dim.nghosty
		nz = mz - 2*dim.nghostz
		if (((size (vars))[1] ne mx) or ((size (vars))[2] ne my) or ((size (vars))[3] ne mz)) then begin
			print, "ERROR: Data doesn't fit to the given dim structure."
			return, -1
		end
	end

	if (n_elements (grid) eq 0) then begin
		; Check consistency of grid
		if (((size (x))[1] ne (size (vars))[1]) or ((size (y))[1] ne (size (vars))[2]) or ((size (z))[1] ne (size (vars))[3])) then begin
			print, "ERROR: Data doesn't fit to the loaded grid structure, please pass the corresponding grid structure as parameter."
			return, -1
		end
	end else begin
		; Set grid in common block for derivative routines
		x = grid.x
		y = grid.y
		z = grid.z
		dx = grid.dx
		dy = grid.dy
		dz = grid.dz
		dx_1 = grid.dx_1
		dy_1 = grid.dy_1
		dz_1 = grid.dz_1
		dx_tilde = grid.dx_tilde
		dy_tilde = grid.dy_tilde
		dz_tilde = grid.dz_tilde
		lequidist = grid.lequidist
		lperi = grid.lperi
		ldegenerated = grid.ldegenerated
		if (((size (x))[1] ne (size (vars))[1]) or ((size (y))[1] ne (size (vars))[2]) or ((size (z))[1] ne (size (vars))[3])) then begin
			print, "Data doesn't fit to the given grid structure."
			return, -1
		end
	end

	; Check availability of requested quantities
	check = quantity
	avail = pc_check_quantities (check=check, sources=sources, /indices, /warn)
	if (not any (avail ge 0)) then return, -1

	if (n_elements (quantity) gt 1) then begin
		; Iterate through availabe quantities
		num = n_elements (avail)
		result = create_struct (quantity[avail[0]], pc_compute_quantity (vars, index, quantity[avail[0]], ghost=ghost))
		if (num gt 1) then begin
			for pos = 1, num-1 do begin
				result = create_struct (result, quantity[avail[pos]], pc_compute_quantity (vars, index, quantity[avail[pos]], ghost=ghost))
			end
		end
	end else if (n_elements (quantity) eq 1) then begin
		; Compute requested quantity:
		result = pc_compute_quantity (vars, index, quantity, ghost=ghost)
	end else begin
		result = !Values.D_NaN
	end

	if (not keyword_set (cache) or keyword_set (cleanup)) then pc_quantity_cache_cleanup

	return, result

end

