; Writes a time-averages snapshot in the HDF5 format

pro pc_write_tavg, varfile, obj, group=group, append_list=append_list, truncate_list=truncate_list, datadir=datadir, dim=dim, grid=grid, start_param=start_param, quiet=quiet

	datadir = pc_get_datadir (datadir)
	default, group, 'data'
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet

	if (not file_test (datadir+'/averages', /directory)) then file_mkdir, datadir+'/averages'

	h5_open_file, datadir+'/averages/'+varfile, /write, /truncate
	labels = strlowcase (tag_names (obj))
	num_labels = n_elements (labels)
	h5_create_group, group
	for pos = 0, num_labels-1 do begin
		label = labels[pos]
		dims = size (obj.(pos), /dimensions)
		num_dims = n_elements (dims)
		if (num_dims eq 4) then begin
			if (dims[3] eq 3) then begin
				label = strmid (label, 0, strlen (label)-1)
				components = [ 'x', 'y', 'z' ]
			end else begin
				components = str (lindgen (dims[3]+1))
			end
			for comp = 0, dims[3]-1 do begin
				h5_write, group+'/'+label+components[comp], reform (obj.(pos)[*,*,*,comp], dims[0:2])
			end
		end else begin
			h5_write, group+'/'+label, obj.(pos)
		end
	end
	h5_write, 'time', time
	h5_close_file

	if (keyword_set (truncate_list)) then begin
		file_delete, list_file, /allow_nonexistent
		lists = file_search (datadir+'/*proc*/'+'tavgN.list')
		if (keyword_set (lists)) then file_delete, lists
	end

	if (keyword_set (append_list) and (varfile ne 'timeavg.dat')) then begin
		list_file = datadir+'/averages/tavgN.list'
		openw, lun, list_file, /get_lun, /append
		printf, lun, varfile
		close, lun
		free_lun, lun
	end

end

