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
function pc_check_quantities, check=check, sources=sources, datadir=datadir, dim=dim, param=param, all=all, warn=warn, indices=indices

	; List of available quantities.
	available = { $
		Temp:'temperature', $
		ln_Temp:'ln temperature', $
		log_Temp:'log temperature', $
		j:'currentdensity', $
		HR_ohm:'Ohmic heating rate', $
		HR_viscous:'viscous heating rate', $
		rho_mag:'magnetic energy', $
		bx:'magnetic field x', $
		by:'magnetic field y', $
		bz:'magnetic field z', $
		u_abs:'velocity', $
		u_x:'velocity x', $
		u_y:'velocity y', $
		u_z:'velocity z', $
		P_therm:'thermal pressure', $
		rho_u_z:'impulse density z', $
		Rn_visc:'viscous Rn', $
		Rn_mag:'magnetic Rn', $
		rho_c:'minimum density', $
		Spitzer_q:'Spitzer heatflux', $
		Spitzer_dt:'Spitzer timestep', $
		rho:'density', $
		ln_rho:'ln density', $
		log_rho:'log density' $
	}

	; List of dependencies.
	depend = { $
		Temp:{ Temp:['lnTT', 'TT'] }, $
		ln_Temp:'Temp', $
		log_Temp:'Temp', $
		B:'aa', $
		j:'B', $
		HR_ohm:'B', $
		HR_viscous:'rho', $
		rho_mag:'B', $
		bx:'B', $
		by:'B', $
		bz:'B', $
		u:'uu', $
		u_abs:'u', $
		u_x:'u', $
		u_y:'u', $
		u_z:'u', $
		P_therm:['Temp', 'rho'], $
		rho_u_z:['u', 'rho'], $
		Rn_visc:'u', $
		Rn_mag:'u', $
		rho_c:'rho', $
		Spitzer_q:'Temp', $
		Spitzer_dt:'Temp', $
		rho:{ rho:['lnrho', 'rho'] }, $
		ln_rho:'rho', $
		log_rho:'rho' $
	}

	; Fill default values
	if (keyword_set (all)) then return, available
	if (not keyword_set (check)) then check = available
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
	for pos = 0, num-1 do begin
		tag = (tag_names (check))[pos]
		index = where (tag_names (available) eq tag)
		if (index ge 0) then begin
			if (dependency_ok (tag, depend, sources)) then begin
				label = available.(index)
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

