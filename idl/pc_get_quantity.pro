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
;   * /cache         If activated, a chache is used to optimize computation.
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
;   P                thermal pressure
;   HR_ohm           volumetric Ohmic heating rate
;   HR_viscous       volumetric viscous heating rate
;   B_z              magnetic field z-component
;   Spitzer_q        absolute value of Spitzer heat flux vector
;   HR_ohm           volumetric Ohmic heating rate
;   j_abs            current density
;   [...]            more are listed in "pc_check_quantities.pro":
;                    IDL> help, pc_check_quantities (/all), /str
;
;  Examples: (in ascending order of efficiency)
;  ============================================
;
;  * Using 'pc_read_var':
;
;   Load varfile and calculate separate quantities, simplest version:
;   IDL> pc_read_var, obj=vars
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', vars)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', vars)
;   IDL> B_z = pc_get_quantity ('B_z', vars)
;   IDL> tvscl, HR_viscous[*,*,20]
;
;   Load varfile and calculate an array of quantities: (RECOMMENDED)
;   IDL> pc_read_var, obj=vars
;   IDL> result = pc_get_quantity (['HR_viscous', 'HR_ohm', 'B_z'], vars)
;   IDL> tvscl, result.HR_viscous[*,*,20]
;
;  * Using 'pc_read_var_raw':
;
;   Load varfile and calculate separate quantities, using a data array:
;   IDL> pc_read_var_raw, obj=var, tags=tags
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', var, tags)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', var, tags)
;   IDL> B_z = pc_get_quantity ('B_z', var, tags)
;   IDL> tvscl, HR_viscous[*,*,20]
;
;   Load varfile and calculate an array of quantities: (RECOMMENDED)
;   IDL> pc_read_var_raw, obj=var, tags=tags, dim=dim, grid=grid, param=param, par2=run_param
;   IDL> result = pc_get_quantity (['HR_viscous', 'HR_ohm', 'B_z'], var, tags, dim=dim, grid=grid, param=param, run_param=run_param)
;   IDL> tvscl, result.HR_viscous[*,*,20]
;
;   Load varfile and separately calculate quantities, using the cache manually:
;   IDL> pc_read_var_raw, obj=var, tags=tags, dim=dim, grid=grid, param=param, par2=run_param
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', var, tags, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', var, tags, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
;   IDL> B_z = pc_get_quantity ('B_z', var, tags, dim=dim, grid=grid, param=param, run_param=run_param, /cache, /cleanup)
;   IDL> tvscl, HR_viscous[*,*,20]
;
;  * Using 'pc_read_slice_raw':
;
;   Load 2D-slice and calculate separate quantities, using the cache manually:
;   IDL> pc_read_slice_raw, obj=slice, tags=tags, cut_z=20, slice_dim=dim, grid=grid, param=param, par2=run_param
;   IDL> HR_viscous = pc_get_quantity ('HR_viscous', slice, tags, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
;   IDL> HR_ohm = pc_get_quantity ('HR_ohm', slice, tags, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
;   IDL> B_z = pc_get_quantity ('B_z', slice, tags, dim=dim, grid=grid, param=param, run_param=run_param, /cache, /cleanup)
;   IDL> tvscl, HR_viscous
;


;  Known issues:
;  =============
;  'q_sat' is untested
;  'rho_c' is untested, only available in SI


; Computation of physical quantities.
; PLEASE ADD MORE PHYSICAL QUANTITIES IN THIS FUNCTION.
; And update the availability and dependency list in "pc_check_quantities.pro".
function pc_compute_quantity, vars, index, quantity

	common quantitiy_cache, uu, rho, grad_rho, n_rho, Temp, grad_Temp, grad_P_therm, bb, jj
	common quantitiy_params, sources, l1, l2, m1, m2, n1, n2, nx, ny, nz, unit, start_par, run_par, alias
	common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0
	common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated

	if (strcmp (quantity, 'u', /fold_case)) then begin
		; Velocity
		if (n_elements (uu) eq 0) then uu = vars[l1:l2,m1:m2,n1:n2,index.uu] * unit.velocity
		return, uu
	end
	if (strcmp (quantity, 'u_x', /fold_case)) then begin
		; Velocity x-component
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		return, uu[*,*,*,0]
	end
	if (strcmp (quantity, 'u_y', /fold_case)) then begin
		; Velocity y-component
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		return, uu[*,*,*,1]
	end
	if (strcmp (quantity, 'u_z', /fold_case)) then begin
		; Velocity z-component
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		return, uu[*,*,*,2]
	end
	if (strcmp (quantity, 'u_abs', /fold_case)) then begin
		; Absolute value of the velocity
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		return, sqrt (dot2 (uu))
	end

	if (strcmp (quantity, 'Temp', /fold_case)) then begin
		; Temperature
		if (n_elements (Temp) eq 0) then begin
			if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
				Temp = exp (vars[l1:l2,m1:m2,n1:n2,index.lnTT]) * unit.temperature
			end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
				Temp = vars[l1:l2,m1:m2,n1:n2,index.TT] * unit.temperature
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
			return, vars[l1:l2,m1:m2,n1:n2,index.lnTT] / alog (10.0) + alog10 (unit.temperature)
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			return, alog10 (vars[l1:l2,m1:m2,n1:n2,index.TT]) + alog10 (unit.temperature)
		end
	end
	if (strcmp (quantity, 'ln_Temp', /fold_case)) then begin
		; Natural logarithmic temperature
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			return, vars[l1:l2,m1:m2,n1:n2,index.lnTT] + alog (unit.temperature)
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			return, alog (vars[l1:l2,m1:m2,n1:n2,index.TT]) + alog (unit.temperature)
		end
	end

	if (strcmp (quantity, 'q_sat', /fold_case)) then begin
		; Absolute value of the saturation heat flux density vector q [W/m^2] = [kg/s^3]
		K_sat = pc_get_parameter ('K_sat', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		return, K_sat * sqrt (Temp / dot2 (pc_compute_quantity (vars, index, 'grad_Temp'))) * (7.28e7 * unit.density * unit.velocity^3 / unit.length * sqrt (unit.temperature))
	end
	if (strcmp (quantity, 'Spitzer_q', /fold_case)) then begin
		; Absolute value of the Spitzer heat flux density vector q [W/m^2] = [kg/s^3]
		K_spitzer = pc_get_parameter ('K_spitzer', label=quantity)
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		return, K_spitzer * Temp^2.5 * sqrt (dot2 (pc_compute_quantity (vars, index, 'grad_Temp'))) * (unit.density * unit.velocity^3 / unit.temperature^3.5 * unit.length)
	end
	if (strcmp (quantity, 'Spitzer_dt', /fold_case)) then begin
		; Spitzer heat flux timestep [s]
		K_spitzer = pc_get_parameter ('K_spitzer', label=quantity)
		cp = pc_get_parameter ('cp', label=quantity)
		gamma = pc_get_parameter ('gamma', label=quantity)
		cdtv = pc_get_parameter ('cdtv', label=quantity)
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			grad_ln_Temp = (grad (vars[*,*,*,index.lnTT]))[l1:l2,m1:m2,n1:n2,*]
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			grad_ln_Temp = (grad (alog (vars[*,*,*,index.TT])))[l1:l2,m1:m2,n1:n2,*]
		end
		dt = cdtv / (gamma * cp * K_spitzer) * rho * sqrt (dot2 (bb) * dot2 (grad_ln_Temp)) / (Temp^2.5 * abs (dot (bb, grad_ln_Temp))) * (unit.time * unit.temperature^2.5 / unit.density)
		; The z-direction may have a non-uniform gird, but not the x- and y-direction
		dxy_1 = dx_1[0]^2 + dy_1[0]^2
		for pz = 0, nz - 1 do dt[*,*,pz] /= dxy_1 + dz_1[pz]^2
		return, dt
	end
	if (strcmp (quantity, 'Spitzer_ratio', /fold_case)) then begin
		; Ratio of perpendicular to parallel Spitzer heat conduction coefficients
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (n_rho) eq 0) then n_rho = pc_compute_quantity (vars, index, 'n_rho')
		return, 2.e-31 * n_rho^2 / (Temp^3 * dot2 (bb)) ; [Solar MHD, E. Priest (1982/1984), p. 86]
	end

	if (strcmp (quantity, 'rho', /fold_case)) then begin
		; Density
		if (n_elements (rho) eq 0) then begin
			if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
				rho = exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * unit.density
			end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
				rho = vars[l1:l2,m1:m2,n1:n2,index.rho] * unit.density
			end
		end
		return, rho
	end
	if (strcmp (quantity, 'grad_rho', /fold_case)) then begin
		; Gradient of density
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
			return, vars[l1:l2,m1:m2,n1:n2,index.lnrho] / alog (10.0) + alog10 (unit.density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			return, alog10 (vars[l1:l2,m1:m2,n1:n2,index.rho]) + alog10 (unit.density)
		end
	end
	if (strcmp (quantity, 'ln_rho', /fold_case)) then begin
		; Natural logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			return, vars[l1:l2,m1:m2,n1:n2,index.lnrho] + alog (unit.density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			return, alog (vars[l1:l2,m1:m2,n1:n2,index.rho]) + alog (unit.density)
		end
	end
	if (strcmp (quantity, 'n_rho', /fold_case)) then begin
		; Particle density
		mu = pc_get_parameter ('mu', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (n_rho) eq 0) then n_rho = rho / (pc_get_parameter ('m_proton', label=quantity) * mu)
		return, n_rho
	end

	if (strcmp (quantity, 'P_therm', /fold_case)) then begin
		; Thermal pressure
		cp = pc_get_parameter ('cp', label=quantity)
		gamma = pc_get_parameter ('gamma', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		return, cp * (gamma - 1.0) / gamma * rho * Temp * unit.density * unit.velocity^2
	end
	if (strcmp (quantity, 'grad_P_therm', /fold_case)) then begin
		; Gradient of thermal pressure
		cp = pc_get_parameter ('cp', label=quantity)
		gamma = pc_get_parameter ('gamma', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (Temp) eq 0) then Temp = pc_compute_quantity (vars, index, 'Temp')
		if (n_elements (grad_rho) eq 0) then grad_rho = pc_compute_quantity (vars, index, 'grad_rho')
		if (n_elements (grad_Temp) eq 0) then grad_Temp = pc_compute_quantity (vars, index, 'grad_Temp')
		if (n_elements (grad_P_therm) eq 0) then begin
			fact = cp * (gamma - 1.0) / gamma * unit.density * unit.temperature / unit.length
			grad_P_therm = grad_rho
			for pa = 0, 2 do grad_P_therm[*,*,*,pa] = fact * (grad_rho[*,*,*,pa] * Temp + rho * grad_Temp[*,*,*,pa])
		end
		return, grad_P_therm
	end
	if (strcmp (quantity, 'grad_P_therm_abs', /fold_case)) then begin
		; Absolute value of thermal pressure gradient
		if (n_elements (grad_P_therm) eq 0) then grad_P_therm = pc_compute_quantity (vars, index, 'grad_P_therm')
		return, sqrt (dot2 (grad_P_therm))
	end

	if (strcmp (quantity, 'rho_u_z', /fold_case)) then begin
		; Impulse density z-component
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		return, rho * pc_compute_quantity (vars, index, 'u_z')
	end

	if (strcmp (quantity, 'rho_c', /fold_case)) then begin
		; Minimum density for an Alfvén speed below the speed of light
		cdtv = pc_get_parameter ('cdtv', label=quantity)
		if (n_elements (rho) eq 0) then rho = pc_compute_quantity (vars, index, 'rho')
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'bb')
		mu0_SI = pc_get_parameter ('mu0_SI', label=quantity)
		return, rho - dot2 (bb) / (2 * mu0_SI * (pc_get_parameter ('c', label=quantity) * cdtv)^2)
	end

	if (strcmp (quantity, 'HR_viscous', /fold_case)) then begin
		; Viscous heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
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

	if (any (strcmp (quantity, ['A', 'A_contour'], /fold_case))) then begin
		; Magnetic vector potential
		return, vars[l1:l2,m1:m2,n1:n2,index.aa] * unit.magnetic_field
	end
	if (strcmp (quantity, 'A_x', /fold_case)) then begin
		; Magnetic vector potential x-component
		return, vars[l1:l2,m1:m2,n1:n2,index.ax] * unit.magnetic_field
	end
	if (strcmp (quantity, 'A_y', /fold_case)) then begin
		; Magnetic vector potential y-component
		return, vars[l1:l2,m1:m2,n1:n2,index.ay] * unit.magnetic_field
	end
	if (strcmp (quantity, 'A_z', /fold_case)) then begin
		; Magnetic vector potential z-component
		return, vars[l1:l2,m1:m2,n1:n2,index.az] * unit.magnetic_field
	end

	if (strcmp (quantity, 'B', /fold_case)) then begin
		; Magnetic field vector
		if (n_elements (bb) eq 0) then bb = (curl (vars[*,*,*,index.aa]))[l1:l2,m1:m2,n1:n2,*] * unit.magnetic_field
		return, bb
	end
	if (strcmp (quantity, 'B_x', /fold_case)) then begin
		; Magnetic field x-component
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, bb[*,*,*,0]
	end
	if (strcmp (quantity, 'B_y', /fold_case)) then begin
		; Magnetic field y-component
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, bb[*,*,*,1]
	end
	if (strcmp (quantity, 'B_z', /fold_case)) then begin
		; Magnetic field z-component
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, bb[*,*,*,2]
	end
	if (strcmp (quantity, 'rho_mag', /fold_case)) then begin
		; Magnetic energy density [WORK HERE: unfinished, currently only computes B^2]
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		return, dot2 (bb)
	end
	if (strcmp (quantity, 'Rn_mag', /fold_case)) then begin
		; Magnetic mesh Reynolds number of velocities perpendicular to the magnetic field
		eta = pc_get_parameter ('eta', label=quantity)
		if (n_elements (bb) eq 0) then bb = pc_compute_quantity (vars, index, 'B')
		if (n_elements (uu) eq 0) then uu = pc_compute_quantity (vars, index, 'u')
		fact = 1.0 / (sqrt (dot2 (uu)) * dot2 (bb))
		BxUxB = cross (cross (bb, uu), bb)
		Rx = spread (pc_compute_quantity (vars, index, 'dx'), [1,2], [ny,nz]) * fact * uu[*,*,*,0] * BxUxB[*,*,*,0]
		Ry = spread (pc_compute_quantity (vars, index, 'dy'), [0,2], [nx,nz]) * fact * uu[*,*,*,1] * BxUxB[*,*,*,1]
		Rz = spread (pc_compute_quantity (vars, index, 'dz'), [0,1], [nx,ny]) * fact * uu[*,*,*,2] * BxUxB[*,*,*,2]
		return, ((Rx > Ry) > Rz) / (eta * unit.length^2/unit.time)
	end

	if (strcmp (quantity, 'j', /fold_case)) then begin
		; Current density
		mu0 = pc_get_parameter ('mu0', label=quantity)
		if (n_elements (jj) eq 0) then jj = (curlcurl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] / mu0 * unit.current_density
		return, jj
	end

	if (strcmp (quantity, 'HR_ohm', /fold_case)) then begin
		; Ohming heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
		mu0 = pc_get_parameter ('mu0', label=quantity)
		eta = pc_get_parameter ('eta', label=quantity)
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, eta * mu0 * dot2 (jj / unit.current_density) * unit.density * unit.velocity^3 / unit.length
	end
	if (strcmp (quantity, 'HR_ohm_particle', /fold_case)) then begin
		; Ohming heating rate per particle [W]
		if (n_elements (n_rho) eq 0) then n_rho = pc_compute_quantity (vars, index, 'n_rho')
		return, pc_compute_quantity (vars, index, 'HR_ohm') / n_rho
	end

	if (strcmp (quantity, 'j_abs', /fold_case)) then begin
		; Current density [A / m^2]
		if (n_elements (jj) eq 0) then jj = pc_compute_quantity (vars, index, 'j')
		return, sqrt (dot2 (jj))
	end

	; Check for Pencil Code alias names
	if (n_elements (alias) eq 0) then alias = pc_check_quantities (/alias)
	tags = strlowcase (tag_names (alias))
	pos = where (tags eq strlowcase (quantity))
	if (any (pos ge 0)) then return, pc_compute_quantity (vars, index, alias.(pos))

	; Timestamp
	if (strcmp (quantity, 'time', /fold_case)) then return, index.time * unit.time

	; Coordinates
	if (strcmp (quantity, 'x', /fold_case)) then return, x[l1:l2] * unit.length
	if (strcmp (quantity, 'y', /fold_case)) then return, y[m1:m2] * unit.length
	if (strcmp (quantity, 'z', /fold_case)) then return, z[n1:n2] * unit.length

	; Grid distances
	if (strcmp (quantity, 'dx', /fold_case)) then return, 1.0 / dx_1[l1:l2] * unit.length
	if (strcmp (quantity, 'dy', /fold_case)) then return, 1.0 / dy_1[m1:m2] * unit.length
	if (strcmp (quantity, 'dz', /fold_case)) then return, 1.0 / dz_1[n1:n2] * unit.length

	; Inverse grid distances
	if (strcmp (quantity, 'inv_dx', /fold_case)) then return, dx_1[l1:l2] / unit.length
	if (strcmp (quantity, 'inv_dy', /fold_case)) then return, dy_1[m1:m2] / unit.length
	if (strcmp (quantity, 'inv_dz', /fold_case)) then return, dz_1[n1:n2] / unit.length

	print, "ERROR: Unknown quantity '"+quantity+"'"
	return, !Values.D_NaN
end


; Clean up cache for computation of physical quantities.
pro pc_quantity_cache_cleanup

	common quantitiy_cache, uu, rho, grad_rho, n_rho, Temp, grad_Temp, grad_P_therm, bb, jj
	common quantitiy_params, sources, l1, l2, m1, m2, n1, n2, nx, ny, nz, unit, start_par, run_par, alias

	undefine, uu
	undefine, rho
	undefine, grad_rho
	undefine, n_rho
	undefine, Temp
	undefine, grad_Temp
	undefine, grad_P_therm
	undefine, bb
	undefine, jj

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
function pc_get_quantity, quantity, vars, index, units=units, dim=dim, grid=grid, param=param, run_param=run_param, datadir=datadir, cache=cache, cleanup=cleanup

	common quantitiy_cache, uu, rho, grad_rho, n_rho, Temp, grad_Temp, grad_P_therm, bb, jj
	common quantitiy_params, sources, l1, l2, m1, m2, n1, n2, nx, ny, nz, unit, start_par, run_par, alias
	common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0
	common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated

	if (keyword_set (cleanup) and not keyword_set (cache)) then pc_quantity_cache_cleanup

	if (n_elements (quantity) eq 0) then quantity = ""
	if (not any (quantity ne "") or (n_elements (vars) eq 0) or ((n_elements (index) eq 0) and (size (vars, /type) ne 8))) then begin
		; Print usage
		print, "USAGE:"
		print, "======"
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

	; Setup 'start.in' and 'run.in' parameters
	dummy = pc_get_parameter ('', start_param=param, run_param=run_param, dim=dim, datadir=datadir)
	start_par = param
	lequidist = safe_get_tag (param, 'lequidist', default=[1,1,1])
	run_par = run_param

	; Set default units
	if (n_elements (units) eq 0) then pc_units, obj=units, datadir=datadir, dim=dim, param=param, /quiet
	unit = units

	if (size (vars, /type) eq 8) then begin
		; Need to have a valid varcontent
		if (size (index, /type) eq 8) then varcontent = index
		if (n_elements (varcontent) eq 0) then varcontent = pc_varcontent (datadir=datadir, dim=dim, param=param)
		; Create array out of given structure and pass recursively computed results
		array = pc_convert_vars_struct (vars, varcontent, index)
		return, pc_get_quantity (quantity, array, index, units=units, dim=dim, grid=grid, param=param, run_param=run_param, datadir=datadir, cache=cache, cleanup=cleanup)
	end

	sources = tag_names (index)

	if (n_elements (dim) eq 0) then begin
		; Check consistency of dimensions
		if (((size (vars))[1] ne mx) or ((size (vars))[2] ne my) or ((size (vars))[3] ne mz)) then begin
			print, "ERROR: Data doesn't fit to the loaded dim structure, please pass the corresponding dim structure as parameter."
			return, -1
		end
		pc_read_dim, obj=glob_dim, datadir=datadir, /quiet
		l1 = glob_dim.nprocx
		l2 = mx - 1 - glob_dim.nprocx
		m1 = glob_dim.nprocy
		m2 = my - 1 - glob_dim.nprocy
		n1 = glob_dim.nprocz
		n2 = mz - 1 - glob_dim.nprocz
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
		result = create_struct (quantity[avail[0]], pc_compute_quantity (vars, index, quantity[avail[0]]))
		if (num gt 1) then begin
			for pos = 1, num-1 do begin
				result = create_struct (result, quantity[avail[pos]], pc_compute_quantity (vars, index, quantity[avail[pos]]))
			end
		end
	end else if (n_elements (quantity) eq 1) then begin
		; Compute requested quantity:
		result = pc_compute_quantity (vars, index, quantity)
	end else begin
		result = !Values.D_NaN
	end

	if (not keyword_set (cache) or keyword_set (cleanup)) then pc_quantity_cache_cleanup

	return, result

end

