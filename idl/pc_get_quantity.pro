;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_get_quantity.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Calculation of physical quantities. If the unit-structure is given,
;;;   this unit system is chosen, otherwise the result is in Pencil-units.
;;;
;;;  Available physical quantities are:
;;;
;;;   Label            Description
;;;  =============================================================
;;;   u_abs            absolute velocity
;;;   uu               velocity
;;;   Temp             temperature
;;;   Spitzer_q        absolute value of Spitzer heat flux vector
;;;   ln_rho           natural logarithm of density
;;;   log_rho          decatic logarithm of density
;;;   rho              density
;;;   P                pressure
;;;   HR_viscous       volumetric viscous heating rate
;;;   B                magnetic field
;;;   HR_ohm           volumetric Ohmic heating rate
;;;   j                current density
;;;
;;;  Examples: (in ascending order of efficiency)
;;;  ============================================
;;;
;;;   Load the most recent varfile and calculate viscous heating rate using a data structure:
;;;   IDL> pc_read_var, obj=vars, dim=dim, grid=grid, param=param, par2=run_param
;;;   IDL> HR_viscous = pc_get_quantity (vars, tags, 'HR_viscous', dim=dim, grid=grid, param=param, run_param=run_param)
;;;   IDL> tvscl, HR_viscous[*,*,20]
;;;
;;;   Load the most recent varfile and calculate viscous heating rate using a data array:
;;;   IDL> pc_read_var_raw, obj=var, tags=tags, dim=dim, grid=grid, param=param, par2=run_param
;;;   IDL> HR_viscous = pc_get_quantity (var, tags, 'HR_viscous', dim=dim, grid=grid, param=param, run_param=run_param)
;;;   IDL> tvscl, HR_viscous[*,*,20]
;;;
;;;   Load only a slice from the most recent varfile and calculate Ohmic heating rate:
;;;   IDL> pc_read_slice_raw, obj=slice, tags=tags, cut_z=20, dim=dim, grid=grid, param=param, par2=run_param
;;;   IDL> HR_ohm = pc_get_quantity (slice, tags, 'HR_ohm', dim=dim, grid=grid, param=param, run_param=run_param)
;;;   IDL> tvscl, HR_ohm
;;;

; Calculation of physical quantities.
function pc_get_quantity, vars, index, quantity, unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, datadir=datadir, cache=cache

	common quantitiy_cache, uu, rho, TT, bb, jj
	common cdat, x, y, z, mx, my, mz, nw
	common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated

	if (not keyword_set (cache)) then begin
		undefine, uu
		undefine, rho
		undefine, TT
		undefine, bb
		undefine, jj
	end

	if (size (vars, /type) eq 8) then begin
		; Create array out of given structure and pass recursively computed results
		array = pc_convert_vars_struct (vars, index, tags)
		return, pc_get_quantity (array, tags, quantity, unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, cache=cache)
	end

	sources = tag_names (index)

	if (n_elements (unit) eq 0) then begin
		; Set default units
		print, "WARNING: using unity as unit."
		unit = { length:1, velocity:1, time:1, temperature:1, density:1, mass:1, magnetic_field:1, current_density:1 }
	end

	; Default data directory
	if (not keyword_set (datadir)) then datadir = pc_get_datadir()

	if (n_elements (dim) eq 0) then begin
		; Check consistency of dimensions
		if (((size (vars))[1] ne mx) or ((size (vars))[2] ne my) or ((size (vars))[3] ne mz)) then $
			message, "Data doesn't fit to the loaded dim structure, please pass the corresponding dim structure as parameter."
		pc_read_dim, obj=glob_dim, datadir=datadir, /quiet
		l1 = glob_dim.nprocx
		l2 = mx - 1 - glob_dim.nprocx
		m1 = glob_dim.nprocy
		m2 = my - 1 - glob_dim.nprocy
		n1 = glob_dim.nprocz
		n2 = mz - 1 - glob_dim.nprocz
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
		if (((size (vars))[1] ne mx) or ((size (vars))[2] ne my) or ((size (vars))[3] ne mz)) then $
			message, "Data doesn't fit to the given dim structure."
	end

	if (n_elements (param) eq 0) then begin
		; Load 'start.in' parameters
		pc_read_param, obj=param, dim=dim, datadir=datadir, /quiet
	end

	if (n_elements (run_param) eq 0) then begin
		; Load 'run.in' parameters
		pc_read_param, obj=run_param, dim=dim, datadir=datadir, /param2, /quiet
	end

	if (n_elements (grid) eq 0) then begin
		; Check consistency of grid
		if (((size (x))[1] ne (size (vars))[1]) or ((size (y))[1] ne (size (vars))[2]) or ((size (z))[1] ne (size (vars))[3])) then $
			message, "Data doesn't fit to the loaded grid structure, please pass the corresponding grid structure as parameter."
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
		if (((size (x))[1] ne (size (vars))[1]) or ((size (y))[1] ne (size (vars))[2]) or ((size (z))[1] ne (size (vars))[3])) then $
			message, "Data doesn't fit to the given grid structure."
	end

	if (n_elements (param) gt 0) then begin
		lequidist = safe_get_tag (param, 'lequidist', default=[1,1,1])
	end

	; Check availability of requested quantities
	avail = pc_check_quantities (check=quantity, sources=sources, /indices, /warn)
	if (not any (avail ge 0)) then return, -1

	if (n_elements (quantity) gt 1) then begin
		; Iterate through availabe quantities
		num = n_elements (avail)
		result = create_struct (quantity[avail[0]], pc_get_quantity (vars, index, quantity[avail[0]], unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache))
		if (num gt 1) then begin
			for pos = 1, num-1 do begin
				result = create_struct (result, quantity[pos], pc_get_quantity (vars, index, quantity[pos], unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache))
			end
		end
		return, result
	end

	; Compute requested quantity:

	if (strcmp (quantity, 'u', /fold_case)) then begin
		; Velocity
		if (n_elements (uu) eq 0) then begin
			uu = vars[l1:l2,m1:m2,n1:n2,index.ux:index.uz] * unit.velocity
		end
		return, uu
	end
	if (strcmp (quantity, 'u_x', /fold_case)) then begin
		; Velocity x-component
		if (n_elements (uu) eq 0) then begin
			return, vars[l1:l2,m1:m2,n1:n2,index.ux] * unit.velocity
		end else begin
			return, uu[*,*,*,0]
		end
	end
	if (strcmp (quantity, 'u_y', /fold_case)) then begin
		; Velocity y-component
		if (n_elements (uu) eq 0) then begin
			return, vars[l1:l2,m1:m2,n1:n2,index.uy] * unit.velocity
		end else begin
			return, uu[*,*,*,1]
		end
	end
	if (strcmp (quantity, 'u_z', /fold_case)) then begin
		; Velocity z-component
		if (n_elements (uu) eq 0) then begin
			return, vars[l1:l2,m1:m2,n1:n2,index.uz] * unit.velocity
		end else begin
			return, uu[*,*,*,2]
		end
	end
	if (strcmp (quantity, 'u_abs', /fold_case)) then begin
		; Absolute velocity
		if (n_elements (uu) eq 0) then uu = pc_get_quantity (vars, index, 'uu', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
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

	if (strcmp (quantity, 'Spitzer_q', /fold_case)) then begin
		; Absolute value of the Spitzer heat flux density vector q [W/m^2] = [kg/s^3]
		if (not any (tag_names (run_param) eq "K_SPITZER")) then message, "Can't compute '"+quantity+"' without parameter 'K_SPITZER'"
		if (n_elements (Temp) eq 0) then Temp = pc_get_quantity (vars, index, 'Temp', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			grad_Temp = (grad (exp (vars[*,*,*,index.lnTT])))[l1:l2,m1:m2,n1:n2,*]
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			grad_Temp = (grad (vars[*,*,*,index.TT]))[l1:l2,m1:m2,n1:n2,*]
		end
		return, run_param.K_spitzer * Temp^2.5 * sqrt (dot2 (grad_temp)) * (unit.density * unit.velocity^3 / unit.temperature^2.5)
	end
	if (strcmp (quantity, 'Spitzer_dt', /fold_case)) then begin
		; Spitzer heat flux timestep [s]
		if (not any (tag_names (run_param) eq "K_SPITZER")) then message, "Can't compute '"+quantity+"' without parameter 'K_SPITZER'"
		if (n_elements (bb) eq 0) then bb = pc_get_quantity (vars, index, 'B', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			grad_ln_Temp = (grad (vars[*,*,*,index.lnTT]))[l1:l2,m1:m2,n1:n2,*]
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			grad_ln_Temp = (grad (alog (vars[*,*,*,index.TT])))[l1:l2,m1:m2,n1:n2,*]
		end
		dt = (run_param.K_spitzer / run_param.cdtv) * abs (dot (bb, grad_ln_Temp)) / sqrt (dot2 (bb * grad_ln_Temp)) * unit.time
		; The z-direction may have a non-uniform gird, but not the x- and y-direction
		dxy_1 = dx_1[0]^2 + dy_[0]^2
		for pz = 0, coord.nz - 1 do dt[*,*,pz] *= dxy_1 + dz_1[pz]^2
		return, dt
	end
	if (strcmp (quantity, 'Spitzer_ratio', /fold_case)) then begin
		; Ratio of perpendicular and parallel Spitzer heat conduction coefficients
		if (n_elements (bb) eq 0) then bb = pc_get_quantity (vars, index, 'B', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		if (n_elements (n_rho) eq 0) then n_rho = pc_get_quantity (vars, index, 'n_rho', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, 2.e-31 * n_rho^2 / (varsets[i].Temp^3 * dot2 (bb))
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
		if (n_elements (rho) eq 0) then rho = pc_get_quantity (vars, index, 'rho', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		if (n_elements (n_rho) eq 0) then begin
			m_p = 1.6726218e-27 ; Mass of a proton [kg]
			n_rho = rho / (m_p * param.mu)
		end
		return, n_rho
	end

	if (strcmp (quantity, 'P', /fold_case)) then begin
		; Thermal pressure
		if (not any (tag_names (run_param) eq "CP") or not any (tag_names (run_param) eq "GAMMA")) then message, "Can't compute '"+quantity+"' without parameter 'CP' or 'GAMMA'"
		if (n_elements (rho) eq 0) then rho = pc_get_quantity (vars, index, 'rho', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		if (n_elements (Temp) eq 0) then Temp = pc_get_quantity (vars, index, 'Temp', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, param.cp * (param.gamma - 1.0) / param.gamma * rho * Temp * unit.density * unit.velocity^2
	end

	if (strcmp (quantity, 'rho_u_z', /fold_case)) then begin
		; Impulse density z-component
		if (n_elements (rho) eq 0) then rho = pc_get_quantity (vars, index, 'rho', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, rho * pc_get_quantity (vars, index, 'u_z', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
	end

	if (strcmp (quantity, 'HR_viscous', /fold_case)) then begin
		; Viscous heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
		if (not any (tag_names (run_param) eq "NU")) then message, "Can't compute '"+quantity+"' without parameter 'NU'"
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
		if (n_elements (rho) eq 0) then rho = pc_get_quantity (vars, index, 'rho', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, run_param.nu * rho * ( 2*((u_xx - div_u3)^2 + (u_yy - div_u3)^2 + (u_zz - div_u3)^2) + (u_xy + u_yx)^2 + (u_xz + u_zx)^2 + (u_yz + u_zy)^2 ) * unit.density * unit.velocity^3 / unit.length
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

	; Magnetic field
	if (strcmp (quantity, 'B', /fold_case)) then begin
		if (n_elements (bb) eq 0) then begin
			bb = (curl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] * unit.magnetic_field
		end
		return, bb
	end
	if (strcmp (quantity, 'B_x', /fold_case)) then begin
		; Magnetic field x-component
		if (n_elements (bb) eq 0) then bb = pc_get_quantity (vars, index, 'B', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, bb[*,*,*,0] * unit.magnetic_field
	end
	if (strcmp (quantity, 'B_y', /fold_case)) then begin
		; Magnetic field y-component
		if (n_elements (bb) eq 0) then bb = pc_get_quantity (vars, index, 'B', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, bb[*,*,*,1] * unit.magnetic_field
	end
	if (strcmp (quantity, 'B_z', /fold_case)) then begin
		; Magnetic field z-component
		if (n_elements (bb) eq 0) then bb = pc_get_quantity (vars, index, 'B', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, bb[*,*,*,2] * unit.magnetic_field
	end
	if (strcmp (quantity, 'rho_mag', /fold_case)) then begin
		; Magnetic energy density [WORK HERE: unfinished, currently only computes B^2]
		if (n_elements (bb) eq 0) then bb = pc_get_quantity (vars, index, 'B', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, dot2 (bb)
	end

	; Current density
	if (strcmp (quantity, 'jj', /fold_case)) then begin
		if (n_elements (jj) eq 0) then begin
			jj = (curlcurl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] / param.mu0 * unit.current_density
		end
		return, jj
	end

	; Ohming heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
	if (strcmp (quantity, 'HR_ohm', /fold_case)) then begin
		if (not any (tag_names (run_param) eq "ETA")) then message, "Can't compute '"+quantity+"' without parameter 'ETA'"
		if (n_elements (jj) eq 0) then jj = pc_get_quantity (vars, index, 'jj', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, run_param.eta * param.mu0 * dot2 (jj / unit.current_density) * unit.density * unit.velocity^3 / unit.length
	end

	; Current density [A / m^2]
	if (strcmp (quantity, 'j_abs', /fold_case)) then begin
		if (n_elements (jj) eq 0) then jj = pc_get_quantity (vars, index, 'jj', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param, /cache)
		return, sqrt (dot2 (jj))
	end

	message, "Unknown quantity '"+quantity+"'"
	return, !Values.D_NaN

end

