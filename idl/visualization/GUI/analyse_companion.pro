;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   analyse_companion.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Framework for precalculation and comparision of output in pencil units.
;;;   Companion procedures needed by 'analyse.pro'.
;;;
;;;  To do:
;;;   Add more comments


; Prepares the varset
pro prepare_varset, num, units, coords, varset, overset

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources

	unit = units
	coord = { x:coords.x*unit.length, y:coords.y*unit.length, z:coords.z*unit.length }

	varfiles = { title:"-", loaded:0, number:0, precalc_done:0 }
	varfiles = replicate (varfiles, num)

	varsets = replicate (varset, num)
	oversets = replicate (overset, num)
end


; Precalculates a data set and loads data, if necessary
pro precalc, i, number=number, varfile=varfile, show_aver=show_aver, vars=vars

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources

	; Default settings
	default, show_aver, 0
	default, number, i

	if (varfiles[i].number le 0) then varfiles[i].number = number

	if (varfiles[i].loaded eq 0) then begin
		default, varfile, "var.dat"
		if (n_elements (vars) eq 0) then begin
			print, 'Reading: ', varfile, ' ... please wait!'
			pc_read_var, varfile=varfile, object=vars, /quiet
			sources = tag_names (vars)
		end
		varfiles[i].title = varfile
		varfiles[i].loaded = 1
		varfiles[i].precalc_done = 0
	end

	if (varfiles[i].precalc_done eq 0) then begin
		print, 'Calculation started...'
		precalc_data, number, vars
		varfiles[i].precalc_done = 1
		print, 'Ready.'
	end

	if (show_aver) then draw_averages, number
end


; Precalculates a data set
pro precalc_data, i, vars

	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources

	tags = tag_names (varsets[i])

	; Compute all desired quantities from available source data
	if (any (strcmp (sources, 'uu', /fold_case))) then begin
		if (any (strcmp (tags, 'u_abs', /fold_case))) then begin
			; Absolute velocity
			varsets[i].u_abs = sqrt (dot2 (vars.uu)) * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_x', /fold_case))) then begin
			; Velocity x-component
			varsets[i].u_x = vars.uu[*,*,*,0] * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_y', /fold_case))) then begin
			; Velocity y-component
			varsets[i].u_y = vars.uu[*,*,*,1] * unit.velocity / unit.default_velocity
		end
		if (any (strcmp (tags, 'u_z', /fold_case))) then begin
			; Velocity z-component
			varsets[i].u_z = vars.uu[*,*,*,2] * unit.velocity / unit.default_velocity
		end
	end
	if (any (strcmp (tags, 'Temp', /fold_case))) then begin
		; Temperature
		if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
			varsets[i].Temp = exp (vars.lnTT) * unit.temperature
		end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
			varsets[i].Temp = vars.TT * unit.temperature
		end
	end
	if (any (strcmp (sources, 'aa', /fold_case))) then begin
		if (any (strcmp (tags, 'Ax', /fold_case))) then begin
			; Magnetic vector potential x-component
			varsets[i].Ax = vars.aa[*,*,*,0]
		end
		if (any (strcmp (tags, 'Ay', /fold_case))) then begin
			; Magnetic vector potential y-component
			varsets[i].Ay = vars.aa[*,*,*,1]
		end
		if (any (strcmp (tags, 'Az', /fold_case))) then begin
			; Magnetic vector potential z-component
			varsets[i].Az = vars.aa[*,*,*,2]
		end
		; Magnetic field
		bb = curl (vars.aa) / unit.length
		if (any (strcmp (tags, 'bx', /fold_case))) then begin
			; Magnetic field x-component
			varsets[i].bx = bb[*,*,*,0]
		end
		if (any (strcmp (tags, 'by', /fold_case))) then begin
			; Magnetic field y-component
			varsets[i].by = bb[*,*,*,1]
		end
		if (any (strcmp (tags, 'bz', /fold_case))) then begin
			; Magnetic field z-component
			varsets[i].bz = bb[*,*,*,2]
		end
		if (any (strcmp (tags, 'rho_mag', /fold_case))) then begin
			; Magnetic energy density
			varsets[i].rho_mag = dot2 (bb)
		end
		if (any (strcmp (tags, 'j', /fold_case))) then begin
			; Current density
			varsets[i].j = sqrt (sqrt (dot2 (curlcurl (vars.aa))) / unit.length^2)
		end
	end
	if (any (strcmp (tags, 'ln_rho', /fold_case))) then begin
		; Natural logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].ln_rho = alog (exp (vars.lnrho) * unit.density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].ln_rho = alog (vars.rho * unit.density)
		end
	end else if (any (strcmp (tags, 'log_rho', /fold_case))) then begin
		; Logarithmic density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].log_rho = alog10 (exp (vars.lnrho) * unit.density)
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].log_rho = alog10 (vars.rho * unit.density)
		end
	end else if (any (strcmp (tags, 'rho', /fold_case))) then begin
		; Density
		if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
			varsets[i].rho = exp (vars.lnrho) * unit.density
		end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
			varsets[i].rho = vars.rho * unit.density
		end
	end

	over_tags = tag_names (oversets[i])
	if (any (strcmp (sources, 'uu', /fold_case))) then begin
		if (any (strcmp (over_tags, 'u', /fold_case))) then begin
			; Velocity overplot
			oversets[i].u = float (vars.uu * unit.velocity / unit.default_velocity)
		end
	end
	if (any (strcmp (sources, 'aa', /fold_case))) then begin
		if (any (strcmp (over_tags, 'b', /fold_case))) then begin
			; Magnetic field overplot
			oversets[i].b = float (bb)
		end
		if (any (strcmp (over_tags, 'a_contour', /fold_case))) then begin
			; Magnetic field lines overplot
			oversets[i].a_contour = float (vars.aa)
		end
	end
end


; Dummy routine
pro analyse_companion

	analyse_companion_loaded = 1
end

