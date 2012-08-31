;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_destretch.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Destretches a given 2D-data array to an equidistant grid by interpolation.
;
;  Parameters:
;   * data           2D-data array.
;   * grid           Original grid spacing.
;   * target         Target grid spacing.
;   * dimension      Dimension to be destretched.
;
;  Example:
;  ========
;
;   IDL> pc_read_slice_raw, obj=var, tags=tags, grid=grid, cut_y=500
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> dz = 1.0 / grid.dz_1
;   IDL> target_dz = mean (dz)
;   IDL> Temp = pc_destretch (Temp, dz, target_dz)
;   IDL> tvscl, alog10 (Temp)
;


; Destretch data from non-uniform to uniform grid.
function pc_destretch, data, grid, target

	if ((n_elements (data) eq 0) or (n_elements (grid) eq 0) or (n_elements (target) eq 0)) then begin
		; Print usage
		print, "USAGE:"
		print, "======"
		print, "pc_destretch, data, grid, target"
		print, ""
		return, -1
	end

	if (not keyword_set (datadir)) then datadir = pc_get_datadir()
	if (n_elements (param) eq 0) then pc_read_param, obj=param, dim=dim, datadir=datadir, /quiet
	if (n_elements (run_param) eq 0) then pc_read_param, obj=run_param, dim=dim, datadir=datadir, /param2, /quiet
	if (n_elements (units) eq 0) then begin
		pc_units, obj=unit, datadir=datadir, dim=dim, param=param, /quiet
		mu0_SI = 4.0 * !Pi * 1.e-7
		unit_current_density = unit.velocity * sqrt (param.mu0 / mu0_SI * unit.density) / unit.length
		units = { length:unit.length, default_length:1, default_length_str:'m', velocity:unit.velocity, default_velocity:1, default_velocity_str:'m/s', time:unit.time, default_time:1, default_time_str:'s', temperature:unit.temperature, default_temperature:1, default_temperature_str:'K', density:unit.density, default_density:1, default_density_str:'kg/m^3', mass:unit.density*unit.length^3, default_mass:1, default_mass_str:'kg', magnetic_field:unit.magnetic_field, default_magnetic_field:1, default_magnetic_field_str:'Tesla', current_density:unit_current_density, default_current_density:1, default_current_density_str:'A/m^2' }
	end

	lequidist = safe_get_tag (param, 'lequidist', default=[1,1,1])

	if (n_elements (target) eq 1) then begin
		tmp = target
		target = grid
		target[*] = tmp
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

