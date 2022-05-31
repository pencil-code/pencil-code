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
function dependency_ok, tag, depend, sources, ALT=ALT

	; Check for dependencies
	if (all (tag eq "")) then return, 1

	if (size (depend, /type) ne 8) then begin
		; Iterate through array of alternative sources
		num = n_elements (tag)
		or_flags = bytarr (num)
		for pos = 0, num-1 do or_flags[pos] = any (strcmp (tag[pos], sources, /fold_case))
		return, any (or_flags)
	end

	if (keyword_set (ALT) and (n_elements (tag) gt 1)) then begin
		; Iterate through array of alternative dependencies
		num = n_elements (tag)
		or_flags = bytarr (num)
		for pos = 0, num-1 do or_flags[pos] = dependency_ok (tag[pos], depend, sources)
		return, any (or_flags)
	end

	; If no dependency is found, check against sources
	index = find_tag (depend, tag)
	if (index eq -1) then return, dependency_ok (tag, -1, sources)

	dependency = depend.(index)

	if (size (dependency, /type) eq 8) then begin
		; Iterate through structure of alternative sources
		num = n_elements (dependency)
		or_flags = bytarr (num)
		for pos = 0, num-1 do begin
			alternative = dependency.(pos)
			if (find_tag (dependency, tag) eq pos) then begin
				or_flags[pos] = dependency_ok (alternative, -1, sources)
			end else begin
				or_flags[pos] = dependency_ok (alternative, depend, sources, /ALT)
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

	; Check dependency against sources
	return, any (strcmp (dependency, sources, /fold_case))
end

; Return available quantities.
function pc_check_quantities, check=check, sources=sources, datadir=datadir, dim=dim, param=param, all=all, available=avail, aliases=aliases, additionals=add_quant, vectorfields=vectorfields, warn=warn, indices=indices

	; Get available sources
	if (not keyword_set (sources)) then sources = pc_varcontent (datadir=datadir, dim=dim, param=param, /quiet)
	if (size (sources, /type) eq 8) then begin
		; Create array of IDL names out of given varcontent structure
		sources = sources.idlvar
		sources = sources[where (sources ne "dummy")]
	end
	num_sources = n_elements (sources)
	num_species = 0
	for pos = 0, num_sources - 1 do begin
		species = (stregex (sources[pos], '^YY([0-9]+)$', /subexpr, /extract, /fold_case))[1]
		if (species ne '') then begin
			num_species += 1
			if (num_species eq 1) then begin
				species_available = create_struct ('species_'+species, 'chemical species '+species+' mass-%')
				species_depend = create_struct ('species_'+species, 'YY'+species)
			end else begin
				species_available = create_struct (species_available, 'species_'+species, 'chemical species '+species+' mass-%')
				species_depend = create_struct (species_depend, 'species_'+species, 'YY'+species)
			end
		end
	end

	; List of available quantities
	available = { $
		Temp:'temperature', $
		ln_Temp:'ln temperature', $
		log_Temp:'log temperature', $
		grad_Temp_abs:'grad temperature', $
		S:'entropy', $
		j_abs:'current density', $
		j_x:'current density x', $
		j_y:'current density y', $
		j_z:'current density z', $
		F_Lorentz_x:'Lorentz force x', $
		F_Lorentz_y:'Lorentz force y', $
		F_Lorentz_z:'Lorentz force z', $
		F_Lorentz_abs:'Lorentz force', $
		W_Lorentz:'Lorentz work', $
		Lorentz_angle:'Lorentz angle', $
		Lorentz_angle_deviation:'Lorentz angle deviation', $
		HR_ohm:'Ohmic heating rate', $
		HR_ohm_particle:'Ohmic heating rate / particle', $
		HR_viscous:'viscous heating rate', $
		HR_viscous_particle:'viscous heating rate / particle', $
		A_abs:'magnetic vector potential', $
		A_x:'magnetic vector potential x', $
		A_y:'magnetic vector potential y', $
		A_z:'magnetic vector potential z', $
		div_A:'divergence of magnetic vector potential', $
		B_abs:'magnetic field strength', $
		B_2:'magnetic field squared', $
		B_x:'magnetic field x', $
		B_y:'magnetic field y', $
		B_z:'magnetic field z', $
		grad_B_abs:'magnetic field gradient', $
		EMF_abs:'electro motive force strength', $
		EMF_x:'electro motive force x', $
		EMF_y:'electro motive force y', $
		EMF_z:'electro motive force z', $
		E_j_abs:'current electric field strength', $
		E_j_x:'current electric field x', $
		E_j_y:'current electric field y', $
		E_j_z:'current electric field z', $
		E_x:'electric field x', $
		E_y:'electric field y', $
		E_z:'electric field z', $
		E_parallel:'field-aligned electric field', $
		E_perpendicular:'field-perpendicular electric field', $
		H_mag:'magnetic field helicity', $
		H_mag_pos:'positive magnetic field helicity', $
		H_mag_neg:'negative magnetic field helicity', $
		H_j:'electric current helicity', $
		dH_mag_dt:'change rate of magnetic field helicity', $
;		grad_E_abs_abs:'grad electric field strength', $
		beta:'plasma beta', $
		rho_mag:'magnetic energy density', $
		Poynting_j_x:'current Poynting flux x', $
		Poynting_j_y:'current Poynting flux y', $
		Poynting_j_z:'current Poynting flux z', $
		Poynting_j_abs:'current Poynting flux', $
		Poynting_u_x:'velocity Poynting flux x', $
		Poynting_u_y:'velocity Poynting flux y', $
		Poynting_u_z:'velocity Poynting flux z', $
		Poynting_u_abs:'velocity Poynting flux', $
		Poynting_x:'Poynting flux x', $
		Poynting_y:'Poynting flux y', $
		Poynting_z:'Poynting flux z', $
		Poynting_abs:'Poynting flux', $
		u_x:'velocity x', $
		u_y:'velocity y', $
		u_z:'velocity z', $
		u_abs:'velocity', $
		grad_u_abs:'velocity gradient', $
		E_therm:'thermal energy', $
		E_kin:'kinetic energy', $
		P_therm:'thermal pressure', $
		grad_P_therm_abs:'grad thermal pressure', $
		c_sound:'sound speed', $
		H_P_therm_x:'thermal pressure scaling height x', $
		H_P_therm_y:'thermal pressure scaling height y', $
		H_P_therm_z:'thermal pressure scaling height z', $
		a_grav_abs:'gravity acceleration', $
		a_grav_x:'gravity acceleration x', $
		a_grav_y:'gravity acceleration y', $
		a_grav_z:'gravity acceleration z', $
		F_grav_abs:'gravity force', $
		F_grav_x:'gravity force x', $
		F_grav_y:'gravity force y', $
		F_grav_z:'gravity force z', $
		rho_u_abs:'impulse density', $
		rho_u_x:'impulse density x', $
		rho_u_y:'impulse density y', $
		rho_u_z:'impulse density z', $
		Rn_viscous:'viscous Reynolds number', $
		Rn_mag:'magnetic Reynolds number', $
		forcing_x:'forcing function x', $
		forcing_y:'forcing function y', $
		forcing_z:'forcing function z', $
		forcing_abs:'forcing function', $
		q_abs:'isotropic heatflux', $
		q_sat:'saturation heatflux', $
		Spitzer_q:'Spitzer heatflux', $
		Spitzer_q_parallel:'field-aligned Spitzer heatflux', $
		Spitzer_q_perpendicular:'field-perpendicular Spitzer heatflux', $
		Spitzer_dt:'Spitzer timestep', $
		Spitzer_ratio:'Spitzer perp./par. heatflux', $
		Spitzer_q_ratio:'saturation/Spitzer heatflux', $
		Spitzer_collision_frequency_e:'Spitzer electron collision frequency', $
		Spitzer_conductivity:'Spitzer conductivity', $
		Spitzer_mag_diffusivity:'Spitzer magnetic diffusivity', $
		Coulomb_logarithm:'Coulomb logarithm', $
		collision_frequency_e:'electron collision frequency', $
		WKB_conductivity:'WKB electric conductivity', $
		WKB_mag_diffusivity:'WKB magnetic diffusivity', $
		c_Alfven:'Alfven velocity', $
		c_Alfven_inv:'inverse Alfven velocity', $
		rho_c:'minimum density (Alfven speed < c)', $
		rho_c_ratio:'density/min. Alfven density', $
		rho:'density', $
		ln_rho:'ln density', $
		log_rho:'log density', $
		n_rho:'particle density',$
		Heat_cool_compression: 'Heating / cooling due to compression',$
		Heat_cool_visc: 'Heating / cooling due to viscous heating',$
		Heat_cool_conduction: 'Heating / cooling due to heat conduction',$
		Heat_cool_rad_loss: 'Heating / cooling due to radiative losses' $
	}

	; List of code variable aliases.
	alias = { $
		t:'time', $
		TT:'Temp', $
		uu:'u', $
		AA:'A', $
		lnrho:'ln_rho', $
		lnTT:'ln_Temp' $
	}

	; List of available vector field quantities.
	available_vectorfields = { $
		u:'velocities', $
		grad_u:'velocity gradient', $
		j:'current density', $
		F_Lorentz:'Lorentz force', $
		Poynting_j:'current Poynting flux', $
		Poynting_u:'velocity Poynting flux', $
		Poynting:'Poynting flux', $
		A:'magnetic vector potential', $
		A_contour:'fieldlines', $
		B:'magnetic field', $
		grad_B:'magnetic field gradient', $
		dB_dx:'magnetic field x-derivative', $
		dB_dy:'magnetic field y-derivative', $
		dB_dz:'magnetic field z-derivative', $
		E:'electric field', $
		EMF:'electro motive force', $
		E_j:'current electric field', $
;		grad_E_abs:'grad electric field strength', $
		grad_Temp:'grad temperature', $
		grad_rho:'grad density', $
		grad_P_therm:'grad thermal pressure', $
		a_grav:'gravity acceleration', $
		F_grav:'gravity force', $
		forcing:'forcing function', $
		rho_u:'impulse density' $
	}

	; Additional quantities without dependencies.
	additional = { $
		time:'timestamp', $
		x:'x coordinates', $
		y:'y coordinates', $
		z:'z coordinates', $
		dx:'grid distance x', $
		dy:'grid distance y', $
		dz:'grid distance z', $
		dV:'grid cell volume', $
		inv_dx:'inverse grid distance x', $
		inv_dy:'inverse grid distance y', $
		inv_dz:'inverse grid distance z', $
		size_x:'box size x', $
		size_y:'box size y', $
		size_z:'box size z', $
		origin_x:'origin x', $
		origin_y:'origin y', $
		origin_z:'origin z', $
		Spitzer_K_parallel:'field-aligned Spitzer coefficient' $
	}

	; List of dependencies.
	; Arrays list a set of mandatory dependencies (see 'HR_viscous').
	; The elements of structures are all mandatory dependencies,
	; while contained arrays list alternative quantities (see 'Temp').
	; If a structure element has the same name as the quantity itself,
	; the contained elements are required to be data sources (see 'TT').
	depend = { $
		TT:{ TT:['lnTT', 'TT'] }, $
		Temp:{ Temp_alternatives:['TT', 'S_rho'] }, $
		S:{ S_alternatives:['ss', 'TT_rho'] }, $
		grad_Temp:'Temp', $
		grad_Temp_abs:'grad_Temp', $
		ln_Temp:'Temp', $
		log_Temp:'Temp', $
		A:'aa', $
		A_contour:'A', $
		B:{ B_alternatives:['A', 'bb'] }, $
		j:'A', $
		j_abs:'j', $
		j_x:'j', $
		j_y:'j', $
		j_z:'j', $
		F_Lorentz:['j', 'B'], $
		F_Lorentz_x:'F_Lorentz', $
		F_Lorentz_y:'F_Lorentz', $
		F_Lorentz_z:'F_Lorentz', $
		F_Lorentz_abs:'F_Lorentz', $
		W_Lorentz:['u', 'F_Lorentz'], $
		Lorentz_angle:['j', 'B'], $
		Lorentz_angle_deviation:'Lorentz_angle', $
		HR_ohm:'j', $
		HR_ohm_particle:['HR_ohm', 'n_rho'], $
		HR_viscous:['u', 'rho'], $
		HR_viscous_particle:['HR_viscous', 'n_rho'], $
		A_abs:'A', $
		A_x:'A', $
		A_y:'A', $
		A_z:'A', $
		div_A:'A', $
		B_abs:'B', $
		B_2:'B', $
		B_x:'B', $
		B_y:'B', $
		B_z:'B', $
		dB_dx:'A', $
		dB_dy:'A', $
		dB_dz:'A', $
		grad_B:'A', $
		grad_B_abs:'grad_B', $
		E:['u','A'], $
		E_abs:'E', $
		EMF:['u','A'], $
		EMF_abs:'EMF', $
		EMF_x:'EMF', $
		EMF_y:'EMF', $
		EMF_z:'EMF', $
		E_j:'j', $
		E_j_abs:'E_j', $
		E_j_x:'E_j', $
		E_j_y:'E_j', $
		E_j_z:'E_j', $
		E_x:'E', $
		E_y:'E', $
		E_z:'E', $
		E_parallel:'E', $
		E_perpendicular:'E', $
		grad_E_abs:'E', $
		grad_E_abs_abs:'E', $
		H_mag:['A','B'], $
		H_mag_pos:'H_mag', $
		H_mag_neg:'H_mag', $
		H_j:['B','j'], $
		dH_mag_dt:['H_j'], $
		beta:['P_therm', 'B_2'], $
		rho_mag:'B_2', $
		Poynting:['u', 'B', 'j'], $
		Poynting_j:['j', 'B'], $
		Poynting_j_x:'Poynting_j', $
		Poynting_j_y:'Poynting_j', $
		Poynting_j_z:'Poynting_j', $
		Poynting_j_abs:'Poynting_j', $
		Poynting_u:['u', 'B'], $
		Poynting_u_x:'Poynting_u', $
		Poynting_u_y:'Poynting_u', $
		Poynting_u_z:'Poynting_u', $
		Poynting_u_abs:'Poynting_u', $
		Poynting_x:'Poynting', $
		Poynting_y:'Poynting', $
		Poynting_z:'Poynting', $
		Poynting_abs:'Poynting', $
		u:'uu', $
		u_x:'u', $
		u_y:'u', $
		u_z:'u', $
		u_abs:'u', $
		grad_u:'u_abs', $
		grad_u_abs:'u_abs', $
		E_therm:['Temp','n_rho'], $
		E_kin:['u','rho'], $
		P_therm:['Temp', 'rho'], $
		grad_P_therm:['P_therm','grad_Temp'], $
		grad_P_therm_abs:'grad_P_therm', $
		c_sound:['P_therm', 'rho'], $
		H_P_therm_x:['P_therm','grad_P_therm'], $
		H_P_therm_y:['P_therm','grad_P_therm'], $
		H_P_therm_z:['P_therm','grad_P_therm'], $
		a_grav:'', $
		a_grav_abs:'a_grav', $
		a_grav_x:'a_grav', $
		a_grav_y:'a_grav', $
		a_grav_z:'a_grav', $
		F_grav:['rho', 'a_grav'], $
		F_grav_abs:'F_grav', $
		F_grav_x:'F_grav', $
		F_grav_y:'F_grav', $
		F_grav_z:'F_grav', $
		rho_u:['u', 'rho'], $
		rho_u_abs:'rho_u', $
		rho_u_x:['u_x', 'rho'], $
		rho_u_y:['u_y', 'rho'], $
		rho_u_z:['u_z', 'rho'], $
		Rn_viscous:'u', $
		Rn_mag:['u','B'], $
		forcing:'ff', $
		forcing_x:'forcing', $
		forcing_y:'forcing', $
		forcing_z:'forcing', $
		forcing_abs:'forcing', $
		q_abs:['Temp'], $
		q_sat:['Temp', 'rho'], $
		Spitzer_K_parallel:'Temp', $
		Spitzer_q:['Temp'], $
		Spitzer_q_parallel:['Temp', 'B'], $
		Spitzer_q_perpendicular:['Spitzer_q_parallel', 'Spitzer_ratio'], $
		Spitzer_dt:['Temp', 'rho', 'B'], $
		Spitzer_ratio:['Temp', 'B', 'n_rho'], $
		Spitzer_q_ratio:['q_sat', 'Spitzer_q'], $
		Spitzer_collision_frequency_e:['Temp', 'n_rho', 'Coulomb_logarithm'], $
		Spitzer_conductivity:['Temp', 'Coulomb_logarithm'], $
		Spitzer_mag_diffusivity:'Spitzer_conductivity', $
		Coulomb_logarithm:'Temp', $
		collision_frequency_e:['B'], $
		WKB_conductivity:['n_rho', 'collision_frequency_e'], $
		WKB_mag_diffusivity:'WKB_conductivity', $
		c_Alfven:['rho', 'B'], $
		c_Alfven_inv:['c_Alfven'], $
		rho_c:['rho', 'B'], $
		rho_c_ratio:['rho', 'rho_c'], $
		rho:{ rho:['lnrho', 'rho'] }, $
		grad_rho:'rho', $
		ln_rho:'rho', $
		log_rho:'rho', $
		n_rho:'rho', $
		; Virtual combined dependencies:
		TT_rho:['TT', 'rho'], $
		S_rho:['S', 'rho'], $
		; Additional quantities with dependencies:
		eta_j:'j', $
		; Additional quantities without dependencies:
		time:'', $
		x:'', $
		y:'', $
		z:'', $
		dx:'', $
		dy:'', $
		dz:'', $
		dV:['dx', 'dy', 'dz'], $
		inv_dx:'dx', $
		inv_dy:'dy', $
		inv_dz:'dz', $
		size_x:'', $
		size_y:'', $
		size_z:'', $
		origin_x:'', $
		origin_y:'', $
		origin_z:'' ,$
		Heat_cool_compression:'uu',$
		Heat_cool_visc:['Temp', 'uu'],$
		Heat_cool_conduction:['Temp', 'rho'],$
		Heat_cool_rad_loss:['Temp', 'rho'] $
	}

	if (num_species ge 1) then begin
		available = create_struct (available, species_available)
		depend = create_struct (depend, species_depend)
	end

	; Return requested listings
	if (keyword_set (all)) then return, create_struct (available, available_vectorfields, alias, additional)
	if (keyword_set (avail)) then return, available
	if (keyword_set (aliases)) then return, alias
	if (keyword_set (add_quant)) then return, additional
	if (keyword_set (vectorfields)) then return, available_vectorfields

	; Fill default values
	if (not keyword_set (check)) then check = create_struct (available)

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
	tags = tag_names (check)
	for pos = 0, num-1 do begin
		tag = tags[pos]
		index = find_tag (avail, tag)
		if (index lt 0) then begin
			index = find_tag (alias, tag)
			if (index ge 0) then begin
				tag = alias.(index)
				index = find_tag (avail, tag)
			end
		end
		if (index lt 0) then begin
			index = find_tag (additional, tag)
			if (index ge 0) then label = additional.(index)
		end else begin
			label = avail.(index)
		end
		if (index ge 0) then begin
			if (dependency_ok (tag, depend, sources)) then begin
				if (num_list eq 0) then begin
					list = create_struct (tag, label)
					pos_list = [ pos ]
				end else begin
					list = create_struct (list, tag, label)
					pos_list = [ pos_list, pos ]
				end
				num_list++
			end else if (keyword_set (warn)) then begin
				print, "WARNING: dependency '"+tag+"' is not available."
			end
		end else if (keyword_set (warn)) then begin
			print, "WARNING: '"+tag+"' is not in the availability list."
		end
	end

	if (keyword_set (indices)) then return, pos_list

	return, list

end
