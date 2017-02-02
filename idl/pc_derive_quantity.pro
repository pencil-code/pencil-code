;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_derive_quantity.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Derive physical quantities. If the unit-structure is given,
;   this unit system is chosen, otherwise the result is in Pencil-units.
;
;  Parameters:
;   * derivative     Derivative direction: 'x', 'y', 'z', 'xz', 'yx', etc.
;   * quantity       Quantity passed on to 'pc_get_quantity', see there.
;   * vars           Data array or data structure as load by pc_read_*.
;   * index          Indices or tags of the given variables inside data.
;
;   
;  Examples: (in ascending order of efficiency)
;  ============================================
;
;  * Using 'pc_read_var_raw': (HIGHLY RECOMMENDED)
;
;   Load varfile and calculate separate quantities, using a data array:
;   IDL> pc_read_var_raw, obj=var, tags=tags
;   IDL> T = pc_get_quantity ('Temp', var, tags)
;   IDL> dT_dz = pc_derive_quantity ('z', 'Temp', var, tags)
;   IDL> d2T_dx_dy = pc_derive_quantity ('xy', 'Temp', var, tags)
;
;   Load varfile and separately calculate quantities, using the cache: (RECOMMENDED FOR SCRIPTS)
;   IDL> pc_read_var_raw, obj=var, tags=tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param
;   IDL> T = pc_get_quantity ('Temp', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
;   IDL> rho = pc_get_quantity ('rho', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
;   IDL> P_therm = pc_get_quantity ('P_therm', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache, /cleanup)
;   IDL> dT_dz = pc_derive_quantity ('z', 'Temp', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param)
;   IDL> d2rho_dx_dy = pc_derive_quantity ('xy', 'rho', var, tags, dim=dim, grid=grid, start_param=start_param, run_param=run_param)
;


; Derivation of physical quantities.
function pc_derive_quantity, derivative, quantity, vars, index, varfile=varfile, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, datadir=datadir, cleanup=cleanup, verbose=verbose

	common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz

	if (size (quantity, /type) eq 7) then begin
		label = "'"+quantity+"'"
		default, vars, "var.dat"
		quantity = pc_get_quantity (quantity, vars, index, varfile=varfile, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, datadir=datadir, /ghost, cleanup=cleanup, verbose=verbose)
	end else begin
		label = 'the given data'
	end

	if (size (quantity, /n_dimensions) lt 3) then message, "pc_derive_quantity: "+label+" is not a derivable 3D quantity."
	if (not all (size (quantity, /dimensions) eq [mx,my,mz])) then message, "pc_derive_quantity: "+label+" can not be derived because it lacks the ghost cells or it is already a derivatvie."

	derivative = strtrim (derivative, 2)
	num_derivatives = strlen (derivative)
	high_level_1st = [ 'div', 'grad', 'curl' ]
	high_level_2nd = [ 'curl_curl', 'grad_curl', 'delta', 'laplace' ]
	if (any (strcmp (derivative, high_level_1st, /fold_case))) then num_derivatives = 1
	if (any (strcmp (derivative, high_level_2nd, /fold_case))) then num_derivatives = 2
	if ((num_derivatives lt 1) or (num_derivatives gt 2)) then message, "pc_derive_quantity: can only compute first and second derivatives."

	unit_length = pc_get_parameter ('unit_length', label=label, dim=dim, datadir=datadir, start_param=start_param, run_param=run_param)

	l1 = nghostx
	l2 = mx - nghostx - 1
	m1 = nghosty
	m2 = my - nghosty - 1
	n1 = nghostz
	n2 = mz - nghostz - 1

	nx = mx - 2 * nghostx
	ny = my - 2 * nghosty
	nz = mz - 2 * nghostz
	na = 1

	vector_quantity = [ 'grad', 'curl', 'curl_curl', 'grad_curl' ]
	if (any (strcmp (derivative, vector_quantity, /fold_case))) then na = 3

	result = dblarr (nx, ny, nz, na, /nozero)
	if (num_derivatives eq 2) then begin
		result[*,*,*,*] = 1.0 / (unit_length^2)
	end else begin
		result[*,*,*,*] = 1.0 / unit_length
	end

	case (derivative) of
		'x': result *= (xder (quantity))[l1:l2,m1:m2,n1:n2]
		'y': result *= (yder (quantity))[l1:l2,m1:m2,n1:n2]
		'z': result *= (zder (quantity))[l1:l2,m1:m2,n1:n2]
		'xx': result *= (xder2    (quantity))[l1:l2,m1:m2,n1:n2]
		'xy': result *= (xderyder (quantity))[l1:l2,m1:m2,n1:n2]
		'xz': result *= (xderzder (quantity))[l1:l2,m1:m2,n1:n2]
		'yx': result *= (yderxder (quantity))[l1:l2,m1:m2,n1:n2]
		'yy': result *= (yder2    (quantity))[l1:l2,m1:m2,n1:n2]
		'yz': result *= (yderzder (quantity))[l1:l2,m1:m2,n1:n2]
		'zx': result *= (zderxder (quantity))[l1:l2,m1:m2,n1:n2]
		'zy': result *= (zderyder (quantity))[l1:l2,m1:m2,n1:n2]
		'zz': result *= (zder2    (quantity))[l1:l2,m1:m2,n1:n2]
		'div': result *= (div (quantity))[l1:l2,m1:m2,n1:n2]
		'grad': result *= (grad (quantity))[l1:l2,m1:m2,n1:n2,*]
		'curl': result *= (curl (quantity))[l1:l2,m1:m2,n1:n2,*]
		'curl_curl': result *= (curlcurl (quantity))[l1:l2,m1:m2,n1:n2,*]
		'grad_curl': result *= (gradcurl (quantity))[l1:l2,m1:m2,n1:n2,*]
		'delta': result *= (del2 (quantity))[l1:l2,m1:m2,n1:n2]
		'laplace': result *= (del2 (quantity))[l1:l2,m1:m2,n1:n2]
		else: message, "pc_derive_quantity: derivative '"+derivative+"' is unknown."
	endcase

	return, result

end

