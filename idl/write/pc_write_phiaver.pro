; Writes phi averages in the HDF5 format

pro pc_write_phiaver, obj, iteration, truncate=truncate, last=last, datadir=datadir, dim=dim, run_param=run_param, quiet=quiet

	datadir = pc_get_datadir (datadir)
	default, truncate, 0
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (run_param)) then pc_read_param, obj=run_param, /run_param, quiet=quiet

	varfile = 'phi.h5'
	h5_open_file, datadir+'/averages/'+varfile, /write, truncate=truncate

	labels = strlowcase (tag_names (obj))
	num_labels = n_elements (labels)
	group = str (iteration)
	h5_create_group, group
	h5_write, group+'/time', obj.t
	for pos = 0, num_labels-1 do begin
		label = labels[pos]
		if (any (label eq [ 't', 'nvars', 'z', 'labels' ])) then continue
		h5_write, group+'/'+label, obj.(pos)
	end

	if (keyword_set (truncate) or not h5_contains ('r')) then begin
		h5_write, 'r', obj.rcyl
		h5_write, 'dr', 2 * run_param.xyz1[0] / dim.nxgrid
	end

	if (size (last, /type) eq 0) then begin
		if (h5_contains ('last')) then begin
			last = h5_read ('last')
			if (last lt iteration) then begin
				last = iteration
				h5_write, 'last', last
			end
		end else begin
			h5_write, 'last', iteration
		end
	end else begin
		if (last lt iteration) then last = iteration
		h5_write, 'last', last
	end

	h5_close_file

end

