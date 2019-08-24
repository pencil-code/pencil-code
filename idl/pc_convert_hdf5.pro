; Converts all snapshots of a run to the HDF5 format

pro pc_convert_hdf5, all=all, old=old, datadir=datadir, delete=delete

	datadir = pc_get_datadir (datadir)
	if (file_test (datadir+'/allprocs/var.dat')) then begin
		procdir = datadir+'/allprocs/'
	end else begin
		procdir = datadir+'/proc0/'
	end

	replace = { lntt:'lnTT', tt:'TT' }
	search = strlowcase (tag_names (replace))
	num_replace = n_elements (search)

	varfiles = 'var.dat'
	if (keyword_set (old)) then varfiles = 'VAR[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'VAR[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	varfiles = strmid (varfiles, strlen (procdir))
	print, varfiles

	num_files = n_elements (varfiles)
	for pos = 0, num_files-1 do begin
		varfile = varfiles[pos]
		if ((varfile eq '') or (strmid (varfile, strlen(varfile)-3) eq '.h5')) then continue
		pc_read_var_raw, obj=var, tags=tags, varfile=varfile, datadir=datadir
		labels = strlowcase (tag_names (tags))
		num_content = n_elements (labels)
		if (varfile eq 'var.dat') then varfile = 'var'
		h5_open_file, datadir+'/allprocs/'+varfile+'.h5', /write, /truncate
		h5_create_group, 'data'
		for i = 0, num_content-1 do begin
			if (size (tags.(i), /n_dimensions) ne 0) then continue
			label = labels[i]
			for j = 0, num_replace-1 do begin
				if (label eq search[j]) then label = replace.(j)
			end
			if (label eq 'time') then begin
				group = '/'
			end else begin
				group = 'data/'
			end
			h5_write, group+label, var[*,*,*,tags.(i)]
		end
		h5_close_file
	end

END

