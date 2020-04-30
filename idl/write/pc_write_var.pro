; Writes a data snapshot in the HDF5 format

pro pc_write_var, varfile, obj, tags=tags, group=group, time=time, append_list=append_list, truncate_list=truncate_list, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, run_param=run_param, quiet=quiet

	datadir = pc_get_datadir (datadir)
	default, group, 'data'
	if (not keyword_set (unit)) then pc_units, obj=unit, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	filename = varfile

	; case-insensitve replacements for dataset names
	replace = { lntt:'lnTT', tt:'TT' }
	search = strlowcase (tag_names (replace))
	num_replace = n_elements (search)

	is_structure = (size (obj, /type) eq 8)

	if (strmid (group, strlen (group)-1) eq '/') then group = strmid (group, 0, strlen (group)-1)
	if (strmid (filename, strlen (filename)-4) eq '.dat') then filename = strmid (filename, 0, strlen (filename)-4)
	if (strmid (filename, strlen (filename)-3) ne '.h5') then filename += '.h5'
	h5_open_file, datadir+'/allprocs/'+filename, /write, /truncate
	if (keyword_set (group)) then h5_create_group, group

	if (is_structure) then begin
		; write from var structure (pc_read_var)
		if (size (varcontent, /type) eq 0) then begin
			varcontent = pc_varcontent (datadir=datadir, dim=dim, param=start_param, par2=run_param, quiet=quiet)
		end
		labels = strlowcase (varcontent[*].idlvar)
		labels = labels[where (labels ne 'dummy')]
	end else begin
		; write from var array (pc_read_var_raw)
		labels = strlowcase (tag_names (tags))
	end

	num_content = n_elements (labels)
	for i = 0, num_content-1 do begin
		label = labels[i]
		if (is_structure) then begin
			if (label eq 'dummy') then continue
		end else begin
			if (size (tags.(i), /n_dimensions) ne 0) then continue
			if (label eq 'time') then continue
		end
		for j = 0, num_replace-1 do begin
			if (label eq search[j]) then label = replace.(j)
		end

		if (is_structure) then begin
			pos = (where (label eq strlowcase (tag_names (obj))))[0]
			if (size (obj.(pos), /n_dimensions) eq 4) then begin
				num_dims = (size (obj.(pos), /dimensions))[3]
				components = [ 'x', 'y', 'z' ]
				for comp = 0, num_dims-1 do begin
					if (num_dims eq 3) then begin
						comp_label = strmid (label, 0, strlen (label)-1) + components[comp]
					end else begin
						comp_label = label + str (comp+1)
					end
					h5_write, group+'/'+comp_label, reform ((obj.(pos))[*,*,*,comp])
				end
			end else begin
				h5_write, group+'/'+label, obj.(pos)
			end
		end else begin
			h5_write, group+'/'+label, obj[*,*,*,tags.(i)]
		end
	end

	if (is_structure and (size (time, /type) eq 0)) then begin
		t_pos = (where ('t' eq strlowcase (tag_names (obj))))[0]
		if (t_pos ge 0) then time = double (obj.(t_pos))
	end

	if (size (time, /type) ne 0) then h5_write, 'time', time

	pc_write_grid, grid=grid, filename='allprocs/'+filename, /append, datadir=datadir, dim=dim, unit=unit, start_param=start_param, quiet=quiet
	h5_close_file

	; update list file
	list_file = datadir+'/allprocs/varN.list'

	if (keyword_set (truncate_list)) then begin
		file_delete, list_file, /allow_nonexistent
		lists = file_search (datadir+'/*proc*/'+'varN.list')
		if (keyword_set (lists)) then file_delete, lists
	end

	if (keyword_set (append_list) and (filename ne 'var.h5')) then begin
		openw, lun, list_file, /get_lun, /append
		printf, lun, filename
		close, lun
		free_lun, lun
	end

end

