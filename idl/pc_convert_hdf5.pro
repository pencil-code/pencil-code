; Converts all snapshots of a run to the HDF5 format

pro pc_convert_hdf5, all=all, old=old, delete=delete, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, run_param=run_param, quiet=quiet

	datadir = pc_get_datadir (datadir)

	if (not keyword_set (unit)) then pc_units, obj=unit, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	if (not keyword_set (run_param)) then pc_read_param, obj=run_param, /run_param, quiet=quiet

	if (keyword_set (all)) then old = all

	if (file_test (datadir+'/allprocs/var.dat')) then begin
		procdir = datadir+'/allprocs/'
	end else begin
		procdir = datadir+'/proc0/'
	end

	; create necessary directories
	if (not file_test (datadir+'/allprocs', /directory)) then file_mkdir, datadir+'/allprocs'
	if (not file_test (datadir+'/averages', /directory)) then file_mkdir, datadir+'/averages'
	if (not file_test (datadir+'/slices', /directory)) then file_mkdir, datadir+'/slices'

	; MHD snapshots
	varfiles = 'var.dat'
	if (keyword_set (old) and not keyword_set (all)) then varfiles = 'VAR[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'VAR[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	varfiles = strmid (varfiles, strlen (procdir))
	varfiles = varfiles[sort (varfiles)]
	num_files = n_elements (varfiles)
	for pos = 0, num_files-1 do begin
		varfile = varfiles[pos]
		if ((varfile eq '') or (strmid (varfile, strlen (varfile)-3) eq '.h5')) then continue
		pc_read_var_raw, obj=data, tags=tags, varfile=varfile, time=time, datadir=datadir, dim=dim, grid=grid, start_param=start_param, run_param=run_param, quiet=quiet
		truncate = keyword_set (delete) and (varfile eq 'VAR0')
		pc_write_var, varfile, data, tags=tags, time=time, append=delete, truncate=truncate, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, run_param=run_param, quiet=quiet
		varfile = varfiles[pos]
		if (keyword_set (delete)) then file_delete, file_search (datadir+'/*proc*/'+varfile)
	end

	; global variables
	varfile = procdir+'global.dat'
	if (file_test (varfile)) then begin
		pc_read_var_time, time=time, datadir=datadir, param=start_param, quiet=quiet
		pc_read_global, obj=data, datadir=datadir, dim=dim, grid=grid, param=start_param, quiet=quiet
		pc_write_global, 'global.h5', data, time=time, datadir=datadir, dim=dim, grid=grid, start_param=start_param, quiet=quiet
		if (keyword_set (delete)) then file_delete, varfile
	end

	; particle snapshots
	varfiles = 'pvar.dat'
	if (keyword_set (old) and not keyword_set (all)) then varfiles = 'PVAR[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'PVAR[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	if (keyword_set (varfiles)) then begin
		varfiles = strmid (varfiles, strlen (procdir))
		num_files = n_elements (varfiles)
		for pos = 0, num_files-1 do begin
			varfile = varfiles[pos]
			if ((varfile eq '') or (strmid (varfile, strlen (varfile)-3) eq '.h5')) then continue
			pc_read_pvar, obj=data, varfile=varfile, datadir=datadir, quiet=quiet
			pc_write_pvar, varfile, data, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param
		end
		if (keyword_set (delete)) then file_delete, file_search (datadir+'/*proc*/'+varfiles)
	end

	; qvar snapshots
	varfiles = 'qvar.dat'
	if (keyword_set (old) and not keyword_set (all)) then varfiles = 'QVAR[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'QVAR[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	if (keyword_set (varfiles)) then begin
		varfiles = strmid (varfiles, strlen (procdir))
		num_files = n_elements (varfiles)
		for pos = 0, num_files-1 do begin
			varfile = varfiles[pos]
			if ((varfile eq '') or (strmid (varfile, strlen (varfile)-3) eq '.h5')) then continue
			pc_read_qvar, obj=data, varfile=varfile, datadir=datadir, quiet=quiet
			pc_write_qvar, varfile, data, datadir=datadir
		end
		if (keyword_set (delete)) then file_delete, file_search (datadir+'/*proc*/'+varfile)
	end

	; stalker particle snapshots
	varfile = 'particles_stalker.dat'
	if ((keyword_set (old) or keyword_set (all)) and file_test (procdir+varfile)) then begin
		pc_read_pstalk, obj=data, datadir=datadir, quiet=quiet
		pc_write_pstalk, data, datadir=datadir, dim=dim, quiet=quiet
		if (keyword_set (delete)) then file_delete, file_search (datadir+'/*proc*/'+varfile)
	end

	; grid file
	pc_write_grid, grid=grid, dim=dim, unit=unit, start_param=start_param

	; time series
	varfile = datadir+'/time_series.dat'
	if (file_test (varfile)) then begin
		pc_read_ts, obj=obj, datadir=datadir, it=iterations, time=time, dt=dt, num=num_labels, labels=labels, quiet=quiet
		h5_open_file, datadir+'/time_series.h5', /write, /truncate
		num_iterations = n_elements (iterations)
		for pos = 0, num_labels-1 do begin
			label = labels[pos]
			if (label eq 'it') then continue
			h5_create_group, label
			for pos_it = 0, num_iterations-1 do begin
				it = str (iterations[pos_it])
				h5_write, label+'/'+it, (obj.(pos))[pos_it]
			end
		end
		h5_write, 'last', iterations[num_iterations-1]
		h5_write, 'step', run_param.it1
		h5_close_file
	end

	; video files
	varfiles = file_search (datadir+'/*proc*/slice_*.*')
	if (file_test ('video.in') and keyword_set (varfiles) and (keyword_set (old) or keyword_set (all))) then begin
		if (not file_test ('src/read_all_videofiles.x')) then spawn, 'make "read_all_videofiles"'
		if (not file_test ('src/read_all_videofiles.x')) then spawn, 'pc_build "read_all_videofiles"'
		if (not file_test ('src/read_all_videofiles.x')) then message, 'Can not build "read_all_videofiles"!'
		spawn, 'src/read_all_videofiles.x'
		varfiles = strmid (varfiles, transpose (strpos (varfiles, '/slice_', /reverse_search) + 7))
		varfiles = varfiles[sort (varfiles)]
		num_files = n_elements (varfiles)
		last = ''
		last_field = ''
		indices = round (reform ((read_ascii (datadir+'/slice_position.dat')).field1[1,*]))
		for pos = 0, num_files-1 do begin
			varfile = varfiles[pos]
			if ((varfile eq last) or (varfile eq 'position.dat')) then continue
			last = varfile
			field = strmid (varfile, 0, strpos (varfile, '.'))
			plane = strmid (varfile, strpos (varfile, '.')+1)
			if (field ne last_field) then begin
				pc_read_video, obj=obj, field=field, /old_format
				num_steps = n_elements (obj.t)
			end
			case (plane) of
				'xy' : begin & index = indices[0] & coord = grid.z[index-1] & end
				'xy2': begin & index = indices[1] & coord = grid.z[index-1] & end
				'xy3': begin & index = indices[2] & coord = grid.z[index-1] & end
				'xy4': begin & index = indices[3] & coord = grid.z[index-1] & end
				'xz' : begin & index = indices[4] & coord = grid.y[index-1] & end
				'xz2': begin & index = indices[5] & coord = grid.y[index-1] & end
				'yz' : begin & index = indices[6] & coord = grid.x[index-1] & end
				else : begin & message, 'Unknown plane "'+plane+'"!' & end
			end
			last_field = field
			h5_open_file, datadir+'/slices/'+field+'_'+plane+'.h5', /write, /truncate
			for step = 1, num_steps do begin
				plane_pos = where (plane eq strlowcase (tag_names (obj)))
				if (plane_pos[0] eq -1) then begin & help, obj, /str & print, (strlowcase (tag_names (obj))+'|') & message, 'Plane "'+plane+'" not found!' & end
				group = str (step)
				h5_create_group, group
				h5_write, group+'/time', obj.t[step-1]
				h5_write, group+'/data', reform ((obj.(plane_pos))[*,*,step-1])
				h5_write, group+'/coordinate', index
				h5_write, group+'/position', coord
			end
			h5_write, 'last', num_steps
			h5_close_file
		end
		if (keyword_set (delete)) then begin
			file_delete, file_search (datadir+[ '/*proc*', '' ]+'/slice_*.*')
		end
	end

	; phi averages
	varfiles = file_search (datadir+'/averages/PHIAVG[0-9]*')
	if (keyword_set (varfiles)) then begin
		num_files = n_elements (varfiles)
		numbers = long (strmid (varfiles, strlen (datadir+'/averages/PHIAVG')))
		varfiles = varfiles[sort (numbers)]
		numbers = numbers[sort (numbers)]
		for file = 0, num_files-1 do begin
			varfile = varfiles[file]
			data = pc_read_phiavg (varfile, datadir=datadir, /old_format)
			pc_write_phiaver, data, numbers[file]-1, truncate=(file eq 0), last=(num_files-1), datadir=datadir, dim=dim, run_param=run_param, quiet=quiet
		end
		if (keyword_set (delete)) then begin
			file_delete, file_search (datadir+'/averages/'+[ 'PHIAVG[0-9]*', 'phiavg.files', 'phiavg.list' ])
		end
	end

	; phi-z averages
	varfile = datadir+'/phizaverages.dat'
	if (file_test (varfile)) then begin
		pc_read_phizaver, obj=data, datadir=datadir, quiet=quiet
		pc_write_phizaver, data, datadir=datadir, quiet=quiet
		if (keyword_set (delete)) then file_delete, varfile
	end

	; 1D averages
	averages = [ 'xy', 'xz', 'yz' ]
	directions = [ 'z', 'y', 'x' ]
	num_aver = n_elements (averages)
	for aver = 0, num_aver-1 do begin
		varfile = datadir+'/'+averages[aver]+'averages.dat'
		if (not file_test (varfile)) then continue
		pc_read_1d_aver, directions[aver], obj=data, datadir=datadir, dim=dim, grid=grid, quiet=quiet
		pc_write_1d_aver, averages[aver]+'.h5', data, datadir=datadir, dim=dim, grid=grid, start_param=start_param, quiet=quiet
		if (keyword_set (delete)) then file_delete, varfile
	end

	; 2D averages
	averages = [ 'y', 'z' ]
	num_aver = n_elements (averages)
	for aver = 0, num_aver-1 do begin
		varfile = procdir+averages[aver]+'averages.dat'
		if (not file_test (varfile)) then continue
		varfile = strmid (varfile, strlen (procdir))
		pc_read_2d_aver, averages[aver], obj=data, datadir=datadir, dim=dim, grid=grid, quiet=quiet
		pc_write_2d_aver, averages[aver]+'.h5', data, datadir=datadir, dim=dim, grid=grid, start_param=start_param, quiet=quiet
		if (keyword_set (delete)) then file_delete, file_search (datadir+'/*proc*/'+varfile)
	end

	; time averages
	varfiles = 'timeavg.dat'
	if (keyword_set (old) and not keyword_set (all)) then varfiles = 'TAVG[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'TAVG[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	if (keyword_set (varfiles)) then begin
		varfiles = strmid (varfiles, strlen (procdir))
		varfiles = varfiles[sort (varfiles)]
		num_files = n_elements (varfiles)
		tavg_vc = pc_varcontent (datadir=datadir, dim=dim, param=start_param, par2=run_param)
		fields = read_ascii (procdir+'tavgN.list')
		times = reform (fields.field1[1,*])
		num_times = n_elements (times)
		for file = 0, num_files-1 do begin
			varfile = varfiles[file]
			if ((varfile eq '') or (strmid (varfile, strlen (varfile)-3) eq '.h5')) then continue
			pc_read_global, obj=data, varfile=varfile, datadir=datadir, dim=dim, grid=grid, param=start_param, varcontent=tavg_vc, quiet=quiet
			if (varfile eq 'timeavg.dat') then begin
				h5_file = 'timeavg.h5'
				pc_read_var_time, time=time, datadir=datadir, param=start_param, quiet=quiet
			end else begin
				h5_file = varfile+'.h5'
				time = times[long (strmid (varfile, 4)) - 1]
			end
			truncate = keyword_set (delete) and (varfile eq 'TAVG1')
			pc_write_tavg, h5_file, data, append=delete, truncate=truncate, datadir=datadir, dim=dim, grid=grid, start_param=start_param, quiet=quiet
			if (keyword_set (delete)) then file_delete, file_search (datadir+'/*proc*/'+varfile)
		end
	end

	; cleanup datadir
	if (keyword_set (delete)) then begin
		if (file_test (datadir+'/proc0')) then file_delete, file_search (datadir+'/proc[0-9]*'), /recursive
		file_delete, datadir+'/allprocs/grid.dat', /allow_nonexistent
		file_delete, datadir+'/allprocs/dim.dat', /allow_nonexistent
		file_delete, datadir+'/dim.dat', /allow_nonexistent
		;file_delete, datadir+'/reduced', /allow_nonexistent, /recursive
	end

END

