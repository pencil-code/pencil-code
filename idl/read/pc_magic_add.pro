; Add magic variables to an existing variable data structure.

pro pc_magic_add, vars, bb=bb, jj=jj, oo=oo, TT=TT, pp=pp, global=global, processor=proc, dim=dim, datadir=datadir, start_param=param

	if (keyword_set (global)) then begin
		; add global parameters to requested quantities (like external magnetic field to 'bb')
		pc_read_global, obj=gg, proc=proc, param=param, proc=proc, dim=dim, datadir=datadir, /quiet
	endif

	if (keyword_set (bb) and not tag_exists (vars, 'bb')) then begin
		; magnetic-field vector
		if (tag_exists (vars, 'aa')) then begin
			vars = create_struct (vars, 'bb', curl (vars.aa))
		endif else begin
			message, "'pc_magic_add': no known source to compute 'bb'!", /info
		endelse
		if (keyword_set (global)) then begin
			if (tag_exists (global, 'bx_ext')) then vars.bb[*,*,*,0] += gg.bx_ext
			if (tag_exists (global, 'by_ext')) then vars.bb[*,*,*,1] += gg.by_ext
			if (tag_exists (global, 'bz_ext')) then vars.bb[*,*,*,2] += gg.bz_ext
		endif
	endif

	if (keyword_set (jj) and not tag_exists (vars, 'jj')) then begin
		; current density
		if (tag_exists (vars, 'aa')) then begin
			vars = create_struct (vars, 'jj', curlcurl (vars.aa))
		endif else if (tag_exists (vars, 'bb')) then begin
			vars = create_struct (vars, 'jj', curl (vars.bb))
		endif else begin
			message, "'pc_magic_add': no known source to compute 'jj'!", /info
		endelse
		if (keyword_set (global)) then begin
			if (tag_exists (global, 'jx_ext')) then vars.jj[*,*,*,0] += gg.jx_ext
			if (tag_exists (global, 'jy_ext')) then vars.jj[*,*,*,1] += gg.jy_ext
			if (tag_exists (global, 'jz_ext')) then vars.bb[*,*,*,2] += gg.jz_ext
		endif
	endif

	if (keyword_set (oo) and not tag_exists (vars, 'oo')) then begin
		; vorticity
		vars = create_struct (vars, 'oo', curl (vars.uu))
	endif

	if (keyword_set (TT) and not tag_exists (vars, 'TT')) then begin
		; temperature
		if (tag_exists (vars, 'lnTT')) then begin
			vars = create_struct (vars, 'TT', exp (vars.lnTT))
		endif else if (tag_exists (vars, 'rho') and tag_exists (vars, 'ss')) then begin
				vars = create_struct (vars, 'TT', pc_eoscalc (vars.rho, vars.ss, /tt, /rho_ss, dim=dim, param=param, datadir=datadir))
		endif else if (tag_exists (vars, 'lnrho') and tag_exists (vars, 'ss')) then begin
				vars = create_struct (vars, 'TT', pc_eoscalc (vars.lnrho, vars.ss, /tt, /lnrho_ss, dim=dim, param=param, datadir=datadir))
		endif else begin
			message, "'pc_magic_add': no known source to compute 'TT'!", /info
		endelse
	endif

	if (keyword_set (pp) and not tag_exists (vars, 'pp')) then begin
		; pressure
		if (tag_exists (vars, 'lnTT')) then begin
			if (tag_exists (vars, 'rho')) then begin
				vars = create_struct (vars, 'pp', pc_eoscalc (vars.rho, vars.lnTT, /pp, /rho_lnTT, dim=dim, param=param, datadir=datadir))
			endif else if (tag_exists (vars, 'lnrho')) then begin
				vars = create_struct (vars, 'pp', pc_eoscalc(vars.lnrho, vars.lnTT, /pp, /lnrho_lnTT, dim=dim, param=param, datadir=datadir))
			endif else begin
				message, "'pc_magic_add': no known source to compute 'pp'! from 'lnTT'", /info
			endelse
		endif else if (tag_exists (vars, 'ss')) then begin
			if (tag_exists (vars, 'rho')) then begin
				vars = create_struct (vars, 'pp', pc_eoscalc (vars.rho, vars.ss, /pp, /rho_ss, dim=dim, param=param, datadir=datadir))
			endif else if (tag_exists (vars, 'lnrho')) then begin
				vars = create_struct (vars, 'pp', pc_eoscalc (vars.lnrho, vars.ss, /pp, /lnrho_ss, dim=dim, param=param, datadir=datadir))
			endif else begin
				message, "'pc_magic_add': no known source to compute 'pp' from 'ss'!", /info
			endelse
		endif else begin
			message, "'pc_magic_add': no known source to compute 'pp'!", /info
		endelse
	endif

end

