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
		pc_write_pvar, varfile, data, datadir=datadir, dim=dim, grid=grid, unit=unit, start_param=start_param
	end

	; qvar snapshots
	varfiles = 'qvar.dat'
	if (keyword_set (old) and not keyword_set (all)) then varfiles = 'QVAR[0-9]*'
	if (keyword_set (all)) then varfiles = [ varfiles, 'QVAR[0-9]*' ]
	varfiles = file_search (procdir+varfiles)
	varfiles = strmid (varfiles, strlen (procdir))

	num_files = n_elements (varfiles)
	for pos = 0, num_files-1 do begin
		varfile = varfiles[pos]
		if ((varfile eq '') or (strmid (varfile, strlen (varfile)-3) eq '.h5')) then continue
		pc_read_qvar, obj=data, varfile=varfile, datadir=datadir
		pc_write_qvar, varfile, data, datadir=datadir
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

END

