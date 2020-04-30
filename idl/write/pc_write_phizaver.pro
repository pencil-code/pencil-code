; Writes phi-z averages in the HDF5 format

pro pc_write_phizaver, obj, datadir=datadir, quiet=quiet

	datadir = pc_get_datadir (datadir)

	varfile = 'phi_z.h5'
	h5_open_file, datadir+'/averages/'+varfile, /write, /truncate
	num_files = n_elements (obj.t)
	labels = strlowcase (tag_names (obj))
	num_labels = n_elements (labels)
	for file = 0, num_files-1 do begin
		group = str (file)
		h5_create_group, group
		h5_write, group+'/time', reform (obj.t[file])
		for pos = 0, num_labels-1 do begin
			label = labels[pos]
			if (any (label eq [ 't', 'rcyl', 'last', 'pos', 'nvars', 'labels' ])) then continue
			h5_write, group+'/'+label, reform (obj.(pos)[*,file])
		end
	end
	h5_write, 'r', obj.rcyl
	h5_write, 'last', num_files-1
	h5_close_file

end

