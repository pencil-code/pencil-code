; Writes 1D averages in the HDF5 format

pro pc_write_1d_aver, varfile, obj, datadir=datadir, dim=dim, grid=grid, start_param=start_param, quiet=quiet

	datadir = pc_get_datadir (datadir)
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet

	if (not file_test (datadir+'/averages', /directory)) then file_mkdir, datadir+'/averages'

	num_files = n_elements (obj.t)
	h5_open_file, datadir+'/averages/'+averages[aver]+'.h5', /write, /truncate
	labels = strlowcase (tag_names (obj))
	num_labels = n_elements (labels)
	for file = 0, num_files-1 do begin
		group = str (file)
		h5_create_group, group
		h5_write, group+'/time', reform (obj.t[file])
		for pos = 0, num_labels-1 do begin
			label = labels[pos]
			if (any (label eq [ 't', 'last', 'pos', 'num_quantities', 'labels', 'x', 'y', 'z' ])) then continue
			h5_write, group+'/'+label, reform (obj.(pos)[*,file])
		end
	end
	h5_write, 'last', num_files-1
	h5_close_file

end

