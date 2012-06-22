;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_gui.pro      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Framework for precalculation and comparision of Pencil VAR* files.
;;;   Calls 'cmp_cslice_cache' for visualisation of a full 3D dataset.
;;;   These routines are intended for usage with Euclidian coordinates.
;;;   Non-equidistant grid coordinates are in principle supported,
;;;   but will be displayed as if they were stretched to an equidistand grid.
;;;
;;;   To run the Graphical User Interface (GUI), please go to a simulation
;;;   directory, open IDL there, and type ".r pc_gui".
;;;
;;;   Optional settings that can be done before starting the GUI:
;;;   IDL> scaling = (0,+oo]            ; magnification factor
;;;   IDL> datadir = "my_data_dir"      ; alternative data directory
;;;   IDL> varfile = "VAR123"           ; default is "var.dat"
;;;   IDL> cut_y = 511                  ; only read an xz-slice at y-pos. 511
;;;   IDL> default_length = 1           ; default length display unit
;;;   IDL> default_length_str = '...'   ; default length string
;;;   IDL> default_velocity = 1         ; default velocity display unit
;;;   IDL> default_velocity_str = '...' ; default velocity string
;;;   IDL> default_density = 1          ; default density display unit
;;;   IDL> default_density_str = '...'  ; default density string
;;;   IDL> default_mass = 1             ; default mass display unit
;;;   IDL> default_mass_str = '...'     ; default mass string
;;;
;;;   The GUI can be closed, but the data stays in memory. Then, the scaling
;;;   parameter can be changed and the GUI can be started again, without the
;;;   need to reload all the data.
;;;
;;;   See the settings section below to select physical quantities for display.
;;;
;;;   At first startup, time series analysis windows are displayed. There,
;;;   only those quantities can be analysed that are listed in 'print.in'.
;;;

; Compile accompanying functions and routines:
@pc_gui_companion
resolve_routine, "cmp_cslice_cache", /COMPILE_FULL_FILE, /NO_RECOMPILE


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Settings that can be changed by the user:

;;;
;;; Physical quantities to be visualized
;;; Available quantities can be found and defined in 'pc_get_quantity'.
;;;
default, quantities, { $
	Temp:'temperature', $
	j_abs:'current density', $
;	B_x:'magnetic field x', $
;	B_y:'magnetic field y', $
	B_z:'magnetic field z', $
	u_abs:'velocity', $
;	u_x:'velocity x', $
;	u_y:'velocity y', $
	u_z:'velocity z', $
	P_therm:'thermal pressure', $
	grad_P_therm_abs:'thermal pressure gradient', $
;	Rn_visc:'viscous Rn', $
;	Rn_mag:'magnetic Rn', $
	log_rho:'logarithmic density' $
}


;;;
;;; Quantities to be overplotted
;;; Available quantities can be found and defined in 'pc_get_quantity'.
;;;
default, overplot_quantities, { $
	u:'velocities', $
;	grad_Temp:'temperature gradient', $
;	grad_P_therm:'thermal pressure gradient', $
;	B:'magnetic_field', $
	A_contour:'fieldlines' $
}


;;;
;;; Preferred units for display
;;;
default, default_length             , 1.e6
default, default_length_str         , 'Mm'
default, default_time               , 1.0
default, default_time_str           , 's'
default, default_velocity           , 1.e3
default, default_velocity_str       , 'km/s'
default, default_density            , 1.0
default, default_density_str        , 'kg/m^3'
default, default_mass               , 1.0
default, default_mass_str           , 'kg'
default, default_magnetic_field     , 1e-4
default, default_magnetic_field_str , 'Gauß'


;;;
;;; Initial varfile
;;;
default, varfile, 'var.dat'
default, crashfile, 'crash.dat'
default, pattern, 'VAR[0-9]*'


;;;
;;; Default data directory
;;;
default, datadir, pc_get_datadir()


;;;
;;; Default technical parameters
;;;
default, cut_x, -1
default, cut_y, -1
default, cut_z, -1
default, data_reduction, 1.0
if (n_elements (data_reduction) eq 1) then data_reduction = replicate (data_reduction, 3)
if (any (data_reduction lt 1.0)) then begin
	print, "Reset data reduction factor to 1.0, must not be smaller."
	print, "(Type .c to continue.)"
	stop
	data_reduction = data_reduction > 1.0
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MAIN PROGRAM:

default, pc_gui_loaded, 0

if (not pc_gui_loaded) then BEGIN

	default, addfile, crashfile

	if (n_elements (nghost) gt 0) then begin
		nghost_x = nghost
		nghost_y = nghost
		nghost_z = nghost
	end

	pc_read_dim, obj=orig_dim, datadir=datadir, /quiet
	pc_read_param, obj=param, dim=orig_dim, datadir=datadir, /quiet
	pc_read_param, obj=run_param, /param2, dim=orig_dim, datadir=datadir, /quiet
	pc_units, obj=unit, datadir=datadir
	mu0_SI = 4.0 * !Pi * 1.e-7
	unit_current_density = unit.velocity * sqrt (param.mu0 / mu0_SI * unit.density) / unit.length
	units = { length:unit.length, default_length:1e6, default_length_str:'Mm', velocity:unit.velocity, default_velocity:1e3, default_velocity_str:'km/s', time:unit.time, default_time:1, default_time_str:'s', temperature:unit.temperature, default_temperature:1, default_temperature_str:'K', density:unit.density, default_density:1, default_density_str:'kg/m^3', mass:unit.density*unit.length^3, default_mass:1, default_mass_str:'kg', magnetic_field:unit.magnetic_field, default_magnetic_field:1e-4, default_magnetic_field_str:'Gauß', current_density:unit_current_density, default_current_density:1, default_current_density_str:'A/m^2' }

	pc_select_files, files=files, num_selected=num_files, pattern=pattern, varfile=varfile, addfile=addfile, datadir=datadir, allprocs=allprocs, procdir=procdir, units=units, param=start_param, run_param=run_param, varcontent=varcontent, var_list=var_list, quantities=quantities, overplots=overplot_quantities, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z, dim=orig_dim
	if ((num_files le 0) or (n_elements (quantities) le 0)) then stop

	if (total([cut_x, cut_y, cut_z] < 0) ge -2) then begin
		pc_read_slice_raw, varfile=varfile, var_list=var_list, dim=dim, grid=grid, datadir=datadir, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z, /trim, /quiet
	end else begin
		dim = orig_dim
		pc_read_grid, obj=grid, dim=dim, datadir=datadir, allprocs=allprocs, /trim, /quiet
	end
	default, nghost_x, dim.nghostx
	default, nghost_y, dim.nghosty
	default, nghost_z, dim.nghostz
	nx = dim.mx - 2*nghost_x
	ny = dim.my - 2*nghost_y
	nz = dim.mz - 2*nghost_z
	disp_size_x = round ((dim.mx - 2*dim.nghostx) / data_reduction[0]) > 1
	disp_size_y = round ((dim.my - 2*dim.nghosty) / data_reduction[1]) > 1
	disp_size_z = round ((dim.mz - 2*dim.nghostz) / data_reduction[2]) > 1

	coords = { x:congrid (grid.x, disp_size_x, 1, 1, /center, /interp)*unit.length/default_length, y:congrid (grid.y, disp_size_y, 1, 1, /center, /interp)*unit.length/default_length, z:congrid (grid.z, disp_size_z, 1, 1, /center, /interp)*unit.length/default_length, dx:congrid (1.0/grid.dx_1, disp_size_x, 1, 1, /center, /interp)*unit.length, dy:congrid (1.0/grid.dy_1, disp_size_y, 1, 1, /center, /interp)*unit.length, dz:congrid (1.0/grid.dz_1, disp_size_z, 1, 1, /center, /interp)*unit.length, nx:disp_size_x, ny:disp_size_y, nz:disp_size_z, l1:dim.nghostx, l2:dim.mx-dim.nghostx-1, m1:dim.nghosty, m2:dim.my-dim.nghosty-1, n1:dim.nghostz, n2:dim.mz-dim.nghostz-1 }


	print, "Allocating memory..."
	dummy = dindgen (coords.nx, coords.ny, coords.nz)
	dummy_3D = findgen (coords.nx, coords.ny, coords.nz, 3)

	; Create varset dummy
	exec_str = "varset = { "
	for i = 0, n_tags (quantities) - 1 do begin
		if (i gt 0) then exec_str += ", "
		exec_str += (tag_names (quantities))[i]+":dummy"
	end
	exec_str += " }"
	res = execute (exec_str)
	if (not res) then begin
		print, "Could not create varset dummy!"
		stop
	end

	; Create overplot dummy
	exec_str = "overplot = { "
	for i = 0, n_tags (overplot_quantities) - 1 do begin
		if (i gt 0) then exec_str += ", "
		exec_str += (tag_names (overplot_quantities))[i]+":dummy_3D"
	end
	exec_str += " }"
	res = execute (exec_str)
	if (not res) then begin
		print, "Could not create overplot dummy!"
		stop
	end

	dummy = 0
	dummy_3D = 0
	print, "...finished."


	pc_gui_prepare_varset, num_files, units, coords, varset, overplot, datadir, param, run_param, var_list

	; Precalculate selected timesteps
	for i = 1, num_files do begin
		pc_gui_precalc, i-1, varfile=files[num_files-i], datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z
	end

	; Mark completition of preparational work
	pc_gui_loaded = 1

END

; scaling factor for visualisation
default, scaling, fix (256.0 / max ([coords.nx, coords.ny, coords.nz]))
if (n_elements (scaling) eq 1) then if (scaling le 0) then scaling = 1


cmp_cslice_cache, quantities, limits=limits, scaling=scaling, coords=coords, overplots=overplot_quantities

window, 0, xsize=8, ysize=8, retain=2
!P.MULTI = [0, 1, 1]
wdelete

end

