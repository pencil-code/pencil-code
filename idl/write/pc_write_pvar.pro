; Writes a data snapshot in the HDF5 format

pro pc_write_pvar, varfile, obj, tags=tags, group=group, time=time, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, quiet=quiet

	datadir = pc_get_datadir (datadir)
	default, group, 'part'
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	if (not keyword_set (unit)) then pc_units, obj=unit, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	filename = varfile

	; case-insensitve replacements for dataset names
	replace = { t:'', x:'', y:'', z:'', dx:'', dy:'', dz:'', distribution:'', npar_found:'', ipar:'ID', xx:'_p', vv:'vp_', vv_cart:'vp__cart' }
	search = strlowcase (tag_names (replace))
	num_replace = n_elements (search)

	if (strmid (group, strlen (group)-1) eq '/') then group = strmid (group, 0, strlen (group)-1)
	if (strmid (filename, strlen (filename)-4) eq '.dat') then filename = strmid (filename, 0, strlen (filename)-4)
	if (strmid (filename, strlen (filename)-3) ne '.h5') then filename += '.h5'
	h5_open_file, datadir+'/allprocs/'+filename, /write, /truncate
	if (keyword_set (group)) then h5_create_group, group

	; write from pvar structure (pc_read_pvar)
	labels = strlowcase (tag_names (obj))
	tags = labels

	num_content = n_elements (labels)
	for i = 0, num_content-1 do begin
		label = labels[i]
		tag_name = label
		for j = 0, num_replace-1 do begin
			if (label eq search[j]) then label = replace.(j)
		end
		if (label eq '') then continue

		pos = (where (tag_name eq strlowcase (tag_names (obj))))[0]
		if (size (obj.(pos), /n_dimensions) eq 2) then begin
			num_dims = (size (obj.(pos), /dimensions))[1]
			components = [ 'x', 'y', 'z' ]
			for comp = 0, num_dims-1 do begin
				if (num_dims eq 3) then begin
					space = strpos (label, '_')
					if (space ge 0) then begin
						comp_label = strmid (label, 0, space) + components[comp] + strmid (label, space+1)
					end else begin
						comp_label = strmid (label, 0, strlen (label)-1) + components[comp]
					end
				end else begin
					comp_label = label + str (comp+1)
				end
				h5_write, group+'/'+comp_label, reform ((obj.(pos))[*,comp])
			end
		end else begin
			h5_write, group+'/'+label, obj.(pos)
		end
	end

	if (size (time, /type) eq 0) then begin
		t_pos = (where ('t' eq strlowcase (tag_names (obj))))[0]
		if (t_pos ge 0) then time = double (obj.(t_pos))
	end

	if (size (time, /type) ne 0) then h5_write, 'time', time

	bounds_x = grid.ox + dindgen (dim.nprocx + 1) / dim.nprocx * grid.Lx
	bounds_y = grid.oy + dindgen (dim.nprocy + 1) / dim.nprocy * grid.Ly
	bounds_z = grid.oz + dindgen (dim.nprocz + 1) / dim.nprocz * grid.Lz
	if (tag_exists (obj, 'distribution')) then begin
		; generate a distribution of particles per processor
		ncpus = dim.nprocx * dim.nprocy * dim.nprocz
		distribution = replicate (obj.npar_found / ncpus, ncpus)
		; correct unequal number of particles per processor
		distribution[0] -= total (distribution) - obj.npar_found
	end else begin
		distribution = obj.npar_found
	end
	h5_create_group, 'proc'
	h5_write, 'proc/bounds_x', bounds_x
	h5_write, 'proc/bounds_y', bounds_y
	h5_write, 'proc/bounds_z', bounds_z
	h5_write, 'proc/distribution', distribution

	pc_write_grid, grid=grid, filename='allprocs/'+filename, /append, datadir=datadir, dim=dim, unit=unit, start_param=start_param, quiet=quiet
	h5_close_file

end

