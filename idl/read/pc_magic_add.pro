; Add magic variables to an existing variable data structure.

pro pc_magic_add, vars, varcontent, bb=bb, jj=jj, oo=oo, TT=TT, pp=pp, global=global, processor=proc, dim=dim, datadir=datadir, start_param=start_param, run_param=run_param

	if (keyword_set (global)) then begin
		; add global parameters to requested quantities (like external magnetic field to 'bb')
		pc_read_global, obj=gg, proc=proc, param=start_param, proc=proc, dim=dim, datadir=datadir, /quiet
		global_names = strlowcase (tag_names (gg))
	endif

	if (keyword_set (bb)) then begin
		; magnetic-field vector
		vars = create_struct (vars, 'bb', curl (vars.aa))
		if (keyword_set (global)) then begin
			if (max (where (global_names eq 'bx_ext')) ge 0) then vars.bb[*,*,*,0] += gg.bx_ext
			if (max (where (global_names eq 'by_ext')) ge 0) then vars.bb[*,*,*,1] += gg.by_ext
			if (max (where (global_names eq 'bz_ext')) ge 0) then vars.bb[*,*,*,2] += gg.bz_ext
		endif
	endif

	if (keyword_set (jj)) then begin
		; current density
		vars = create_struct (vars, 'jj', curlcurl (vars.aa))
		if (keyword_set (global)) then begin
			if (max (where (global_names eq 'jx_ext')) ge 0) then vars.jj[*,*,*,0] += gg.jx_ext
			if (max (where (global_names eq 'jy_ext')) ge 0) then vars.jj[*,*,*,1] += gg.jy_ext
			if (max (where (global_names eq 'jz_ext')) ge 0) then vars.bb[*,*,*,2] += gg.jz_ext
		endif
	endif

	if (keyword_set (oo)) then begin
		; vorticity
		vars = create_struct (vars, 'oo', curl (vars.uu))
	endif

	if (keyword_set (TT)) then begin
		; temperature
		if (size (lionization, /type) eq 0) then begin
			lionization = safe_get_tag (param, 'lionization', default=safe_get_tag (param, 'leos_ionization', default=0))
			lionization_fixed = safe_get_tag (param, 'lionization_fixed', default=safe_get_tag (param, 'leos_ionizationi_fixed', default=0))
		endif
		if (lionization and not lionization_fixed) then begin
			vars = create_struct (vars, 'TT', exp (vars.lnTT))
		endif else begin
			if (param.ldensity_nolog) then begin
				vars = create_struct (vars, 'TT', pc_eoscalc (vars.rho, vars.ss, /tt, /rho_ss, dim=dim, param=start_param, datadir=datadir))
			endif else begin
				vars = create_struct (vars, 'TT', pc_eoscalc (vars.lnrho, vars.ss, /tt, /lnrho_ss, dim=dim, param=start_param, datadir=datadir))
			endelse
		endelse
	endif

	if (keyword_set (pp)) then begin
		; pressure
		if (size (lionization, /type) eq 0) then begin
			lionization = safe_get_tag (param, 'lionization', default=safe_get_tag (param, 'leos_ionization', default=0))
			lionization_fixed = safe_get_tag (param, 'lionization_fixed', default=safe_get_tag (param, 'leos_ionizationi_fixed', default=0))
		endif
		if (lionization and not lionization_fixed) then begin
			if (param.ldensity_nolog) then begin
				vars = create_struct (vars, 'pp', pc_eoscalc (vars.rho, vars.lnTT, /pp, /rho_lnTT, dim=dim, param=param, datadir=datadir))
			endif else begin
				vars = create_struct (vars, 'pp', pc_eoscalc(vars.lnrho, vars.lnTT, /pp, /lnrho_lnTT, dim=dim, param=param, datadir=datadir))
			endelse
		endif else begin
			if (param.ldensity_nolog) then begin
				vars = create_struct (vars, 'pp', pc_eoscalc (vars.rho, vars.ss, /pp, /rho_ss, dim=dim, param=start_param, datadir=datadir))
			endif else begin
				vars = create_struct (vars, 'pp', pc_eoscalc (vars.lnrho, vars.ss, /pp, /lnrho_ss, dim=dim, param=start_param, datadir=datadir))
			endelse
		endelse
	endif

end

