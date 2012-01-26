;
;  Loads a varfile dataset, including dim, grid and units.
;  The nx, ny, nz, and lmn12 variables are set.
;
;  This file is intended to be included from scripts.
;  Please refer to 'pc_load_dataset.pro' to see an example.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Default settings:
default, varfile, 'var.dat'
default, datadir, pc_get_datadir()
;
pc_units, obj=unit, datadir=datadir
pc_read_grid, obj=grid, dim=dim, param=param, datadir=datadir, allprocs=allprocs
pc_read_var, varfile=varfile, datadir=datadir, object=vars, param=param, par2=run_param, dim=dim, allprocs=allprocs, /quiet
;
default, nghost_x, dim.nghostx
default, nghost_y, dim.nghosty
default, nghost_z, dim.nghostz
nx = dim.mx - 2*nghost_x
ny = dim.my - 2*nghost_y
nz = dim.mz - 2*nghost_z
lmn12 = nghost_x+spread(indgen(nx),[1,2],[ny,nz]) + dim.mx*(nghost_y+spread(indgen(ny),[0,2],[nx,nz])) + dim.mx*dim.my*(nghost_z+spread(indgen(nz),[0,1],[nx,ny]))
;
procdir = datadir+"/proc0/"
if (not file_test (procdir+varfile)) then begin
  procdir = datadir+"/allprocs/"
  if (not file_test (procdir+varfile)) then begin
    print, "No '"+varfile+"' file found."
    stop
  endif
endif
;
varfile_loaded = varfile
;
end
