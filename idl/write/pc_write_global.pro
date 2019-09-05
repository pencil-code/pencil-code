; Writes a global data snapshot in the HDF5 format

pro pc_write_global, varfile, obj, time=time, group=group, datadir=datadir, dim=dim, grid=grid, start_param=start_param, quiet=quiet

	datadir = pc_get_datadir (datadir)
	default, group, 'data'
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet

	if ((size (time, /type) eq 0) and tag_exists (obj, 't')) then time = obj.t

	labels = strlowcase (tag_names (obj))
	num_labels = n_elements (labels)
	h5_open_file, datadir+'/'+varfile, /write, /truncate
	h5_create_group, group
	for pos = 0, num_labels-1 do begin
		label = labels[pos]
		dims = size (obj.(pos), /dimensions)
		num_dims = n_elements (dims)
		if (num_dims eq 4) then begin
			if (dims[3] eq 3) then begin
				components = [ 'x', 'y', 'z' ]
			end else begin
				components = str (lindgen (dims[3]+1))
			end
			for comp = 0, dims[3]-1 do begin
				h5_write, group+'/global_'+label+components[comp], reform (obj.(pos)[*,*,*,comp], dims[0:2])
			end
		end else begin
			h5_write, group+'/global_'+label, obj.(pos)
		end
	end
	if (size (time, /type) eq 0) then h5_write, 'time', time
	h5_close_file

end

