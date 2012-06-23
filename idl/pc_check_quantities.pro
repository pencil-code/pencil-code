;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_check_quantities.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Returns a list of available quantities accepted by 'pc_get_quantity'
;;;   and depending on the optionally given varcontent.
;;;   A list can also be checked and only the valid quantities are returned.
;;;
;;;  Parameters:
;;;   * check		list of quantities to be checked for availability
;;;   * sources		array of varcontent IDL names or a varcontent structure
;;;   * /all		return all available quantities, without checking
;;;   * /warn		warn about missing dependencies in the varcontent
;;;   * /indices	list of indices, instead of the quantities itself
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

; Check if a dependency is fulfilled.
function dependency_ok, tag, depend, sources

	if (size (depend, /type) ne 8) then begin
		; Iterate through array of alternative sources without dependencies
		num = n_elements (tag)
		or_flags = bytarr (num)
		for pos = 0, num-1 do or_flags[pos] = any (strcmp (tag[pos], sources, /fold_case))
		return, any (or_flags)
	end

	index = where (strcmp (tag, tag_names (depend), /fold_case))
	if (any (index eq -1)) then return, dependency_ok (tag, -1, sources)
	dependency = depend.(index)

	if (size (dependency, /type) eq 8) then begin
		; Iterate through structure of alternative sources
		num = n_elements (dependency)
		or_flags = bytarr (num)
		for pos = 0, num-1 do begin
			if (strcmp (tag, (tag_names (dependency))[pos], /fold_case)) then begin
				or_flags[pos] = dependency_ok (dependency.(pos), -1, sources)
			end else begin
				or_flags[pos] = dependency_ok (dependency.(pos), depend, sources)
			end
		end
		return, any (or_flags)
	end

	num = n_elements (dependency)
	if (num gt 0) then begin
		; Iterate through array of mandatory sources
		and_flags = bytarr (num)
		for pos = 0, num-1 do and_flags[pos] = dependency_ok (dependency[pos], depend, sources)
		return, all (and_flags)
	end

	; Check dependency agains sources
	return, any (strcmp (dependency, sources, /fold_case))
end

; Return available quantities.
function pc_check_quantities, check=check, sources=sources, datadir=datadir, dim=dim, param=param, all=all, available=avail, aliases=aliases, vectorfields=vectorfields, warn=warn, indices=indices

	; List of available quantities.
	available = { $
		Temp:'temperature', $
		ln_Temp:'ln temperature', $
		log_Temp:'log temperature', $
		grad_Temp_abs:'grad temperature', $
		j_abs:'current density', $
		HR_ohm:'Ohmic heating rate', $
		HR_viscous:'viscous heating rate', $
		rho_mag:'magnetic energy', $
		A_x:'magnetic vector potential x', $
		A_y:'magnetic vector potential y', $
		A_z:'magnetic vector potential z', $
		B_x:'magnetic field x', $
		B_y:'magnetic field y', $
		B_z:'magnetic field z', $
		u_x:'velocity x', $
		u_y:'velocity y', $
		u_z:'velocity z', $
		u_abs:'velocity', $
		P_therm:'thermal pressure', $
		grad_P_therm_abs:'grad thermal pressure', $
		rho_u_z:'impulse density z', $
		Rn_visc:'viscous Rn', $
		Rn_mag:'magnetic Rn', $
		q_sat:'saturation heatflux', $
		Spitzer_q:'Spitzer heatflux', $
		Spitzer_dt:'Spitzer timestep', $
		Spitzer_ratio:'Spitzer perp./par. heatflux', $
		rho_c:'minimum density (Alfven < c)', $
		rho:'density', $
		ln_rho:'ln density', $
		log_rho:'log density', $
		n_rho:'particle density' $
	}

	; List of code variable aliases.
	alias = { $
		TT:'Temp', $
		uu:'u', $
		AA:'A', $
		lnrho:'ln_rho', $
		lnTT:'ln_Temp' $
	}

	; List of available overplot quantities.
	available_vectorfields = { $
		u:'velocities', $
		j:'current density', $
		A:'magnetic vector potential', $
		A_contour:'fieldlines', $
		B:'magnetic field', $
		grad_Temp:'grad temperature', $
		grad_P_therm:'grad thermal pressure' $
	}

	; List of dependencies.
	depend = { $
		Temp:{ Temp:['lnTT', 'TT'] }, $
		grad_Temp:'Temp', $
		grad_Temp_abs:'grad_Temp', $
		ln_Temp:'Temp', $
		log_Temp:'Temp', $
		A:'aa', $
		A_contour:'A', $
		B:'A', $
		j:'A', $
		j_abs:'j', $
		HR_ohm:'j', $
		HR_viscous:['u', 'rho'], $
		rho_mag:'B', $
		A_x:'A', $
		A_y:'A', $
		A_z:'A', $
		B_x:'B', $
		B_y:'B', $
		B_z:'B', $
		u:'uu', $
		u_x:'u', $
		u_y:'u', $
		u_z:'u', $
		u_abs:'u', $
		P_therm:['Temp', 'rho'], $
		grad_P_therm:'P_therm', $
		grad_P_therm_abs:'grad_P_therm', $
		rho_u_z:['u', 'rho'], $
		Rn_visc:'u', $
		Rn_mag:['u','B'], $
		q_sat:['Temp', 'rho'], $
		Spitzer_q:['Temp'], $
		Spitzer_dt:['Temp', 'rho', 'B'], $
		Spitzer_ratio:['Temp', 'B', 'n_rho'], $
		rho_c:['rho', 'B'], $
		rho:{ rho:['lnrho', 'rho'] }, $
		ln_rho:'rho', $
		log_rho:'rho', $
		n_rho:'rho' $
	}

	; Fill default values
	if (keyword_set (all)) then return, create_struct (available, alias)
	if (keyword_set (avail)) then return, available
	if (keyword_set (aliases)) then return, alias
	if (keyword_set (vectorfields)) then return, available_vectorfields
	if (not keyword_set (check)) then check = create_struct (available, alias)
	if (not keyword_set (sources)) then sources = pc_varcontent (datadir=datadir, dim=dim, param=param, /quiet)

	if (size (sources, /type) eq 8) then begin
		; Create array of IDL names out of given varcontent structure
		sources = sources.idlvar
		sources = sources[where (sources ne "dummy")]
	end

	if (size (check, /type) eq 7) then begin
		; Create structure out of given array
		names = check
		num = n_elements (names)
		if (num ge 1) then check = create_struct (names[0], names[0])
		if (num ge 2) then for pos = 1, num-1 do check = create_struct (check, names[pos], names[pos])
	end

	; Perform check and build list of available quantities, depending on the actual varcontent
	num = n_tags (check)
	list = ""
	pos_list = -1
	num_list = 0
	avail = create_struct (available, available_vectorfields)
	avail_list = tag_names (avail)
	alias_list = tag_names (alias)
	tags = tag_names (check)
	for pos = 0, num-1 do begin
		tag = tags[pos]
		index = where (avail_list eq tag)
		if (index lt 0) then begin
			index = where (alias_list eq tag)
			if (index ge 0) then begin
				tag = alias.(index)
				index = where (strcmp (avail_list, tag, /fold_case))
			end
		end
		if (index ge 0) then begin
			if (dependency_ok (tag, depend, sources)) then begin
				label = avail.(index)
				if (num_list eq 0) then begin
					list = create_struct (tag, label)
					pos_list = [ pos ]
				end else begin
					list = create_struct (list, tag, label)
					pos_list = [ pos_list, pos ]
				end
				num_list++
			end else if (keyword_set (warn)) then begin
				print, "WARNING: dependency '"+tag+"' not available."
			end
		end else if (keyword_set (warn)) then begin
			print, "WARNING: '"+tag+"' is not in availablility list."
		end
	end

	if (keyword_set (indices)) then return, pos_list

	return, list

end

