; Default settings
default, varfile, 'var.dat'
default, datadir, pc_get_datadir()


pc_read_dim, obj=dim, /quiet
default, nghost_x, dim.nghostx
default, nghost_y, dim.nghosty
default, nghost_z, dim.nghostz
nx = dim.mx - 2*nghost_x
ny = dim.my - 2*nghost_y
nz = dim.mz - 2*nghost_z
lmn12 = nghost_x+spread(indgen(nx),[1,2],[ny,nz]) + dim.mx*(nghost_y+spread(indgen(ny),[0,2],[nx,nz])) + dim.mx*dim.my*(nghost_z+spread(indgen(nz),[0,1],[nx,ny]))

procdir = datadir+"/proc0/"
file_struct = file_info (procdir+varfile)
if (file_struct.exists eq 0) then begin
	procdir = datadir+"/allprocs/"
	file_struct = file_info (procdir+varfile)
	if (file_struct.exists eq 0) then begin
		print, "No '"+varfile+"' file found."
		stop
	endif
endif

pc_units, obj=unit
pc_read_grid, obj=grid, /quiet
pc_read_var, varfile=varfile, object=vars, /quiet

varfile_loaded = varfile


END

