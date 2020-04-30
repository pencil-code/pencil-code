; Writes a stalker particles file in the HDF5 format

pro pc_write_pstalk, obj, datadir=datadir, dim=dim, quiet=quiet

	datadir = pc_get_datadir (datadir)
	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet

	num_files = n_elements (obj.t)
	num_particles = n_elements (obj.ipar)
	num_procs = dim.nprocx * dim.nprocy * dim.nprocz
	distribution = replicate (num_particles / num_procs, num_procs)
	if (num_procs ge 2) then distribution[num_procs-1] += num_particles - total (distribution)

	for pos = 0, num_files-1 do begin
		varfile = 'PSTALK'+str (pos)+'.h5'
		h5_open_file, datadir+'/allprocs/'+varfile, /write, /truncate
		h5_write, 'time', obj.t[pos]

		h5_create_group, 'proc'
		h5_write, 'proc/distribution', distribution

		h5_create_group, 'stalker'
		h5_write, 'stalker/ID', obj.ipar
		h5_write, 'stalker/ap', reform (obj.ap[*,pos])
		h5_write, 'stalker/npswarm', reform (obj.npswarm[*,pos])
		h5_write, 'stalker/rho', reform (obj.rho[*,pos])
		h5_write, 'stalker/ux', reform (obj.ux[*,pos])
		h5_write, 'stalker/uy', reform (obj.uy[*,pos])
		h5_write, 'stalker/uz', reform (obj.uz[*,pos])
		h5_write, 'stalker/vpx', reform (obj.vpx[*,pos])
		h5_write, 'stalker/vpy', reform (obj.vpy[*,pos])
		h5_write, 'stalker/vpz', reform (obj.vpz[*,pos])
		h5_write, 'stalker/xp', reform (obj.xp[*,pos])
		h5_write, 'stalker/yp', reform (obj.yp[*,pos])
		h5_write, 'stalker/zp', reform (obj.zp[*,pos])

		h5_close_file
	end

end

