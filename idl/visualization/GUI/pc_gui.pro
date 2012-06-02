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
;;;
; Available quantities for visualization are:
; 'Temp', 'rho'             ; temperature and density
; 'ln_rho', 'log_rho'       ; logarithmic densities
; 'ux', 'uy', 'uz', 'u_abs' ; velocity components and the absolute value
; 'Ax', 'Ay', 'Az'          ; vector potential components
; 'Bx', 'By', 'Bz'          ; magnetic field components
; 'rho_mag'                 ; magnetic energy density
; 'j'                       ; absolute value of the current density
; 'HR_ohm'                  ; ohmic heating rate per volume (eta*mu0*j^2)
; 'rho_u_z'                 ; vertical component of the impulse density
; 'P_therm'                 ; thermal pressure
; 'Rn_visc'                 ; viscous mesh Reynolds number
; 'Rn_mag'                  ; magnetic mesh Reynolds number
; (more quantities can be defined in 'precalc_data', see pc_gui_companion.pro)
default, quantities, { $
	temperature:'Temp', $
	currentdensity:'j', $
	ohmic_heating_rate:'HR_ohm', $
	viscous_heating_rate:'HR_viscous', $
;	magnetic_energy:'rho_mag', $
;	magnetic_field_x:'bx', $
;	magnetic_field_y:'by', $
	magnetic_field_z:'bz', $
	velocity:'u_abs', $
;	velocity_x:'u_x', $
;	velocity_y:'u_y', $
	velocity_z:'u_z', $
;	thermal_pressure:'P_therm', $
	impulse_density_z:'rho_u_z', $
;	viscous_Rn:'Rn_visc', $
;	magnetic_Rn:'Rn_mag', $
;	minimum_density:'rho_c', $
	Spitzer_heatflux:'Spitzer_q', $
	Spitzer_timestep:'Spitzer_dt', $
	logarithmic_density:'log_rho' $
	}


;;;
;;; Quantities to be overplotted (calculated in 'precalc_data')
;;;
; Available quantities for overplotting are:
; 'b', 'a_contour', 'grad_P_therm', and 'u'
default, overplot_quantities, { $
;	magnetic_field:'b', $
	fieldlines:'a_contour', $
;	thermal_pressure_gradient:'grad_P_therm', $
	velocities:'u' $
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

	pc_units, obj=unit, datadir=datadir
	units = { velocity:unit.velocity, time:unit.time, temperature:unit.temperature, length:unit.length, density:unit.density, mass:unit.density*unit.length^3, magnetic_field:unit.magnetic_field, default_length:default_length, default_time:default_time, default_velocity:default_velocity, default_density:default_density, default_mass:default_mass, default_magnetic_field:default_magnetic_field, default_length_str:default_length_str, default_time_str:default_time_str, default_velocity_str:default_velocity_str, default_density_str:default_density_str, default_mass_str:default_mass_str, default_magnetic_field_str:default_magnetic_field_str }
	pc_read_dim, obj=orig_dim, datadir=datadir, /quiet
	pc_read_param, obj=param, dim=orig_dim, datadir=datadir, /quiet
	pc_read_param, obj=run_param, /param2, dim=orig_dim, datadir=datadir, /quiet

	pc_select_files, files=files, num_selected=num_files, pattern=pattern, varfile=varfile, addfile=addfile, datadir=datadir, allprocs=allprocs, procdir=procdir, units=units, param=start_param, run_param=run_param
	if (num_files eq 0) then stop

	if (total([cut_x, cut_y, cut_z] < 0) ge -2) then begin
		pc_read_slice_raw, varfile=varfile, dim=dim, grid=grid, datadir=datadir, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z, /trim, /quiet
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
		exec_str += quantities.(i)+":dummy"
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
		exec_str += overplot_quantities.(i)+":dummy_3D"
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


	prepare_varset, num_files, units, coords, varset, overplot, datadir, param, run_param

	; Precalculate selected timesteps
	for i = 1, num_files do begin
		precalc, i-1, varfile=files[num_files-i], datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z
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

