;
;  Loads a varfile dataset, including dim, grid and units.
;  The nx, ny, nz, and lmn12 variables are set, optionally.
;
;  To use this code in interactive mode or as include file in scripts,
;  please try '.r pc_load_dataset' or use '@pc_load_dataset.inc'.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; This line is executed on a '.r pc_load_dataset':
@pc_load_dataset.inc
;
;
; Same code as procedure, using the include file:
;
pro pc_load_dataset, varfile, datadir=datadir, procdir=procdir, vars=vars, dim=dim, grid=grid, unit=unit, param=param, run_param=run_param, nx=nx, ny=ny, nz=nz, nghost_x=nghost_x, nghost_y=nghost_y, nghost_z=nghost_z, lmn12=lmn12
;
  default, varfile, 'var.dat'
  default, varfile_loaded, ''
  if (varfile ne varfile_loaded) then begin
    @pc_load_dataset.inc
;
end
