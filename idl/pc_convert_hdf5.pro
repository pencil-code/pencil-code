; Converts all snapshots of a run to the HDF5 format

pro pc_convert_hdf5, all=all, old=old, delete=delete, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, run_param=run_param

	datadir = pc_get_datadir (datadir)

	if (not keyword_set (unit)) then pc_units, obj=unit, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet

	if (file_test (datadir+'/allprocs/var.dat')) then begin
		procdir = datadir+'/allprocs/'
	end else begin
		procdir = datadir+'/proc0/'
	end

	; MHD snapshots
	varfiles = 'var.dat'
	if (keyword_set (old) and not keyword_set (all)) then varfiles = 'VAR[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'VAR[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	varfiles = strmid (varfiles, strlen (procdir))

	num_files = n_elements (varfiles)
	for pos = 0, num_files-1 do begin
		varfile = varfiles[pos]
		if ((varfile eq '') or (strmid (varfile, strlen (varfile)-3) eq '.h5')) then continue
		pc_read_var_raw, obj=data, tags=tags, varfile=varfile, time=time, datadir=datadir, dim=dim, grid=grid, start_param=start_param, run_param=run_param
		pc_write_var, varfile, data, tags=tags, time=time, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, run_param=run_param
	end

	; particle snapshots
	varfiles = 'pvar.dat'
	if (keyword_set (old) and not keyword_set (all)) then varfiles = 'PVAR[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'PVAR[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	varfiles = strmid (varfiles, strlen (procdir))

	num_files = n_elements (varfiles)
	for pos = 0, num_files-1 do begin
		varfile = varfiles[pos]
		if ((varfile eq '') or (strmid (varfile, strlen (varfile)-3) eq '.h5')) then continue
		pc_read_pvar, obj=data, varfile=varfile, datadir=datadir
		pc_write_pvar, varfile, data, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, run_param=run_param
	end
END

