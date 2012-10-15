;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_gui.pro      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Framework for precalculation and comparision of Pencil VAR* files.
;;;   Calls 'cmp_cslice_cache' for visualisation of a full 3D dataset.
;;;   These routines are intended for usage with Euclidian coordinates.
;;;   Non-equidistant grids can be displayed as-is or de-stretched.
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
;;; Load file with user-defined default settings.
@pc_gui_settings

;;; Default data directory
default, datadir, pc_get_datadir()

;;; Default minimum size of the data display
default, min_display_size, 256

;;; Default technical parameters
default, cut_x, -1
default, cut_y, -1
default, cut_z, -1
default, reduced, 0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MAIN PROGRAM:

default, pc_gui_loaded, 0

if (not pc_gui_loaded) then BEGIN

	default, addfile, crashfile

	pc_read_dim, obj=orig_dim, datadir=datadir, reduced=reduced, /quiet
	pc_read_param, obj=param, dim=orig_dim, datadir=datadir, /quiet
	pc_read_param, obj=run_param, /param2, dim=orig_dim, datadir=datadir, /quiet
	pc_units, obj=unit, datadir=datadir, param=param, dim=orig_dim, /quiet
	unit = create_struct (unit, display_units)

	; Scaling factor for visualisation
	default, scaling, fix (min_display_size / max ([orig_dim.nx, orig_dim.ny, orig_dim.nz]))
	if (n_elements (scaling) eq 1) then if (scaling le 0) then scaling = 1

	pc_select_files, files=files, num_selected=num_files, pattern=pattern, varfile=varfile, addfile=addfile, datadir=datadir, allprocs=allprocs, reduced=reduced, procdir=procdir, unit=unit, param=start_param, run_param=run_param, varcontent=varcontent, var_list=var_list, quantities=quantities, overplots=overplot_quantities, cut_x=cut_x, cut_y=cut_y, cut_z=cut_z, xs=xs, xe=xe, ys=ys, ye=ye, zs=zs, ze=ze, dim=orig_dim, scaling=scaling
	if ((num_files le 0) or (n_elements (quantities) le 0)) then stop

	if ((xe-xs lt orig_dim.nx-1) or (ye-ys lt orig_dim.ny-1) or (ze-zs lt orig_dim.nz-1)) then begin
		pc_read_subvol_raw, varfile=files[0], var_list=['none'], dim=orig_dim, sub_dim=dim, sub_grid=grid, datadir=datadir, xs=xs, xe=xe, ys=ys, ye=ye, zs=zs, ze=ze, allprocs=allprocs, reduced=reduced, /addghosts, /trim, /quiet
	end else begin
		dim = orig_dim
		pc_read_grid, obj=grid, dim=dim, datadir=datadir, allprocs=allprocs, reduced=reduced, /trim, /quiet
	end

	coords = { $
		x:grid.x * unit.length/unit.default_length, $
		y:grid.y * unit.length/unit.default_length, $
		z:grid.z * unit.length/unit.default_length, $
		dx:1.0/grid.dx_1 * unit.length, $
		dy:1.0/grid.dy_1 * unit.length, $
		dz:1.0/grid.dz_1 * unit.length, $
		nx:dim.nx, ny:dim.ny, nz:dim.nz, $
		orig_nx:orig_dim.nx, orig_ny:orig_dim.ny, orig_nz:orig_dim.nz, $
		x_off:xs, y_off:ys, z_off:zs, $
		l1:dim.nghostx, l2:dim.mx-dim.nghostx-1, $
		m1:dim.nghosty, m2:dim.my-dim.nghosty-1, $
		n1:dim.nghostz, n2:dim.mz-dim.nghostz-1, $
		lequidist:grid.lequidist, lperi:grid.lperi, ldegenerated:grid.ldegenerated }


	print, "Allocating memory..."
	dummy = dindgen (dim.nx, dim.ny, dim.nz)
	dummy_3D = findgen (dim.nx, dim.ny, dim.nz, 3)

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


	pc_gui_prepare_varset, num_files, unit, coords, varset, overplot, datadir, param, run_param, var_list

	; Precalculate selected timesteps
	for i = 1, num_files do begin
		pc_gui_precalc, i-1, varfile=files[num_files-i], datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs, reduced=reduced, xs=xs, xe=xe, ys=ys, ye=ye, zs=zs, ze=ze
	end

	; Mark completition of preparational work
	pc_gui_loaded = 1

END


cmp_cslice_cache, quantities, limits=limits, scaling=scaling, coords=coords, overplots=overplot_quantities

window, 0, xsize=8, ysize=8, retain=2
!P.MULTI = [0, 1, 1]
wdelete

end

