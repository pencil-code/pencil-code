; Converts all snapshots of a run to the HDF5 format

pro pc_convert_hdf5, all=all, old=old, delete=delete, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param, run_param=run_param

	datadir = pc_get_datadir (datadir)

	if (not keyword_set (unit)) then pc_units, obj=unit, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, datadir=datadir, dim=dim, param=start_param, quiet=quiet
	if (not keyword_set (run_param)) then pc_read_param, obj=run_param, /run_param, quiet=quiet

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

	; global variables
	varfile = procdir+'global.dat'
	if (file_test (varfile)) then begin
		pc_read_var_time, time=time, datadir=datadir, param=start_param
		pc_read_global, obj=data, datadir=datadir, dim=dim, grid=grid, param=start_param
		labels = strlowcase (tag_names (data))
		num_labels = n_elements (labels)
		h5_open_file, datadir+'/global.h5', /write, /truncate
		h5_create_group, 'data'
		for pos = 0, num_labels-1 do begin
			label = labels[pos]
			dims = size (data.(pos), /dimensions)
			num_dims = n_elements (dims)
			if (num_dims eq 4) then begin
				if (dims[3] eq 3) then begin
					components = [ 'x', 'y', 'z' ]
				end else begin
					components = str (lindgen (dims[3]+1))
				end
				for comp = 0, dims[3]-1 do begin
					h5_write, 'data/global_'+label+components[comp], reform (data.(pos)[*,*,*,comp], dims[0:2])
				end
			end else begin
				h5_write, 'data/global_'+label, data.(pos)
			end
		end
		h5_write, 'time', time
		h5_close_file
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
			pc_read_pvar, obj=data, varfile=varfile, datadir=datadir
			pc_write_pvar, varfile, data, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param
		end
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
			pc_read_qvar, obj=data, varfile=varfile, datadir=datadir
			pc_write_qvar, varfile, data, datadir=datadir
		end
	end

	; stalker particle snapshots
	varfile = 'particles_stalker.dat'
	if ((keyword_set (old) or keyword_set (all)) and file_test (procdir+varfile)) then begin
		pc_read_pstalk, obj=data, datadir=datadir
		num_files = n_elements (data.t)
		num_particles = n_elements (data.ipar)
		num_procs = dim.nprocx * dim.nprocy * dim.nprocz
		distribution = replicate (num_particles / num_procs, num_procs)
		if (num_procs ge 2) then distribution[num_procs-1] += num_particles - total (distribution)
		for pos = 0, num_files-1 do begin
			h5_open_file, datadir+'/allprocs/PSTALK'+str(pos)+'.h5', /write, /truncate
			h5_write, 'time', data.t[pos]
			h5_create_group, 'proc'
			h5_write, 'proc/distribution', distribution
			h5_create_group, 'stalker'
			h5_write, 'stalker/ID', data.ipar
			h5_write, 'stalker/ap', reform (data.ap[*,pos])
			h5_write, 'stalker/npswarm', reform (data.npswarm[*,pos])
			h5_write, 'stalker/rho', reform (data.rho[*,pos])
			h5_write, 'stalker/ux', reform (data.ux[*,pos])
			h5_write, 'stalker/uy', reform (data.uy[*,pos])
			h5_write, 'stalker/uz', reform (data.uz[*,pos])
			h5_write, 'stalker/vpx', reform (data.vpx[*,pos])
			h5_write, 'stalker/vpy', reform (data.vpy[*,pos])
			h5_write, 'stalker/vpz', reform (data.vpz[*,pos])
			h5_write, 'stalker/xp', reform (data.xp[*,pos])
			h5_write, 'stalker/yp', reform (data.yp[*,pos])
			h5_write, 'stalker/zp', reform (data.zp[*,pos])
			h5_close_file
		end
	end

	; grid file
	h5_open_file, datadir+'/grid.h5', /write, /truncate
	h5_create_group, 'grid'
	h5_write, 'grid/Lx', grid.Lx
	h5_write, 'grid/Ly', grid.Ly
	h5_write, 'grid/Lz', grid.Lz
	h5_write, 'grid/x', grid.x
	h5_write, 'grid/y', grid.y
	h5_write, 'grid/z', grid.z
	h5_write, 'grid/dx', grid.dx
	h5_write, 'grid/dy', grid.dy
	h5_write, 'grid/dz', grid.dz
	h5_write, 'grid/dx_1', grid.dx_1
	h5_write, 'grid/dy_1', grid.dy_1
	h5_write, 'grid/dz_1', grid.dz_1
	h5_write, 'grid/dx_tilde', grid.dx_tilde
	h5_write, 'grid/dy_tilde', grid.dy_tilde
	h5_write, 'grid/dz_tilde', grid.dz_tilde
	h5_create_group, 'settings'
	h5_write, 'settings/l1', dim.l1
	h5_write, 'settings/l2', dim.l2
	h5_write, 'settings/m1', dim.m1
	h5_write, 'settings/m2', dim.m2
	h5_write, 'settings/n1', dim.n1
	h5_write, 'settings/n2', dim.n2
	h5_write, 'settings/nx', dim.nxgrid
	h5_write, 'settings/ny', dim.nygrid
	h5_write, 'settings/nz', dim.nzgrid
	h5_write, 'settings/mx', dim.mxgrid
	h5_write, 'settings/my', dim.mygrid
	h5_write, 'settings/mz', dim.mzgrid
	h5_write, 'settings/nghost', dim.nghostx
	h5_write, 'settings/nprocx', dim.nprocx
	h5_write, 'settings/nprocy', dim.nprocy
	h5_write, 'settings/nprocz', dim.nprocz
	h5_write, 'settings/mvar', dim.mvar
	h5_write, 'settings/maux', dim.maux
	h5_write, 'settings/mglobal', dim.mglobal
	h5_write, 'settings/precision', dim.precision
	h5_write, 'settings/version', 0
	h5_create_group, 'unit'
	h5_write, 'unit/density', unit.density
	h5_write, 'unit/energy', unit.energy
	h5_write, 'unit/flux', unit.energy / (unit.length^2 * unit.time)
	h5_write, 'unit/length', unit.length
	h5_write, 'unit/magnetic', unit.magnetic_field
	h5_write, 'unit/mass', unit.mass
	h5_write, 'unit/system', unit.system
	h5_write, 'unit/temperature', unit.temperature
	h5_write, 'unit/time', unit.time
	h5_write, 'unit/velocity', unit.velocity
	h5_close_file

	; time series
	varfile = datadir+'/time_series.dat'
	if (file_test (varfile)) then begin
		pc_read_ts, obj=obj, datadir=datadir, it=iterations, time=time, dt=dt, num=num_labels, labels=labels
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
	varfiles = file_search (datadir+'/proc*/slice_*.*')
	if (keyword_set (varfiles) and (keyword_set (old) or keyword_set (all))) then begin
		if (not file_test ('src/read_all_videofiles.x')) then spawn, 'pc_build "read_all_videofiles"'
		if (not file_test ('src/read_all_videofiles.x')) then message, 'Can not build "read_all_videofiles"!'
		spawn, 'src/read_all_videofiles.x'
		varfiles = strmid (varfiles, transpose (strpos (varfiles, '/slice_', /reverse_search) + 7))
		varfiles = varfiles[sort (varfiles)]
		num_files = n_elements (varfiles)
		last = ''
		last_field = ''
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
				'xy': begin & index = run_param.iz & coord = grid.z[index-1] & end
				'xy2': begin & index = run_param.iz2 & coord = grid.z[index-1] & end
				'xy3': begin & index = run_param.iz3 & coord = grid.z[index-1] & end
				'xy4': begin & index = run_param.iz4 & coord = grid.z[index-1] & end
				'xz': begin & index = run_param.iy & coord = grid.y[index-1] & end
				'yz': begin & index = run_param.ix & coord = grid.x[index-1] & end
				else: begin & message, 'Unknown plane "'+plane+'"!' & end
			end
			last_field = field
			if (not file_test (datadir+'/slices', /directory)) then file_mkdir, datadir+'/slices'
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
	end

	; phi averages
	varfiles = file_search (datadir+'/averages/PHIAVG[0-9]*')
	if (keyword_set (varfiles)) then begin
		num_files = n_elements (varfiles)
		numbers = long (strmid (varfiles, strlen (datadir+'/averages/PHIAVG')))
		varfiles = varfiles[sort (numbers)]
		numbers = numbers[sort (numbers)]
		h5_open_file, datadir+'/averages/phi.h5', /write, /truncate
		for file = 0, num_files-1 do begin
			varfile = varfiles[file]
			data = pc_read_phiavg (varfile, datadir=datadir, /old_format)
			labels = strlowcase (tag_names (data))
			num_labels = n_elements (labels)
			group = str (numbers[file]-1)
			h5_create_group, group
			h5_write, group+'/time', data.t
			for pos = 0, num_labels-1 do begin
				label = labels[pos]
				if (any (label eq [ 't', 'nvars', 'z', 'labels' ])) then continue
				h5_write, group+'/'+label, data.(pos)
			end
		end
		h5_write, 'r', data.rcyl
		h5_write, 'dr', 2 * run_param.xyz1[0] / dim.nxgrid
		h5_write, 'last', num_files-1
		h5_close_file
	end

	; phi-z averages
	varfile = datadir+'/phizaverages.dat'
	if (file_test (varfile)) then begin
		pc_read_phizaver, obj=data, datadir=datadir
		num_files = n_elements (data.t)
		h5_open_file, datadir+'/averages/phi_z.h5', /write, /truncate
		labels = strlowcase (tag_names (data))
		num_labels = n_elements (labels)
		for file = 0, num_files-1 do begin
			group = str (file)
			h5_create_group, group
			h5_write, group+'/time', reform (data.t[file])
			for pos = 0, num_labels-1 do begin
				label = labels[pos]
				if (any (label eq [ 't', 'rcyl', 'last', 'pos', 'nvars', 'labels' ])) then continue
				h5_write, group+'/'+label, reform (data.(pos)[*,file])
			end
		end
		h5_write, 'r', data.rcyl
		h5_write, 'last', num_files-1
		h5_close_file
		if (keyword_set (delete)) then file_delete, datadir+varfile
	end

END

