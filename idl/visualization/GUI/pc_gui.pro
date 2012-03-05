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
;	ohmic_heating_rate:'HR_ohm', $
;	magnetic_energy:'rho_mag', $
;	magnetic_field_x:'bx', $
;	magnetic_field_y:'by', $
	magnetic_field_z:'bz', $
;	spitzer_ratio:'spitzer_ratio', $
	velocity:'u_abs', $
;	velocity_x:'u_x', $
;	velocity_y:'u_y', $
	velocity_z:'u_z', $
;	thermal_pressure:'P_therm', $
	impulse_density_z:'rho_u_z', $
;	viscous_Rn:'Rn_visc', $
;	magnetic_Rn:'Rn_mag', $
;	minimum_density:'rho_c', $
;	Spitzer_heatflux:'Spitzer_q', $
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


;;;
;;; Default data directory
;;;
default, datadir, pc_get_datadir()


;;;
;;; Default technical parameters
;;;
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

	allprocs = 1
	procdir = datadir+"/allprocs/"
	default, load_varfile, 1
	if (load_varfile and not file_test (procdir+varfile)) then begin
		allprocs = 0
		procdir = datadir+"/proc0/"
		if (not file_test (procdir+varfile)) then begin
			print, "No '"+varfile+"' file found."
			load_varfile = 0
			stop
		endif
	endif

	default, addfile, crashfile
	if ((addfile ne varfile) and file_test (procdir+addfile)) then begin
		print, "A '"+addfile+"' exists, do you want to load it instead of '"+varfile+"'? / as additional file?"
		repeat begin
			answer = "y"
			read, answer, format="(A)", prompt="(Y)es / (N)o / (A)dditional: "
		end until (any (strcmp (answer, ['n', 'y', 'a'], /fold_case)))
		if (strcmp (answer, 'y', /fold_case)) then varfile = addfile
		if (not strcmp (answer, 'a', /fold_case)) then addfile = ""
	end else begin
		addfile = ""
	end
	if (addfile) then num_additional = 1 else num_additional = 0

	if (n_elements (nghost) gt 0) then begin
		nghost_x = nghost
		nghost_y = nghost
		nghost_z = nghost
	end
	pc_read_dim, obj=dim, datadir=datadir, /quiet
	default, nghost_x, dim.nghostx
	default, nghost_y, dim.nghosty
	default, nghost_z, dim.nghostz
	nx = dim.mx - 2*nghost_x
	ny = dim.my - 2*nghost_y
	nz = dim.mz - 2*nghost_z
	disp_size_x = round ((dim.mx - 2*dim.nghostx) / data_reduction[0]) > 1
	disp_size_y = round ((dim.my - 2*dim.nghosty) / data_reduction[1]) > 1
	disp_size_z = round ((dim.mz - 2*dim.nghostz) / data_reduction[2]) > 1

	subdomains = dim.nprocx * dim.nprocy * dim.nprocz
	ghosts = 2*nghost_x*(dim.nprocx-1)*dim.mygrid*dim.mzgrid + $
                 2*nghost_y*(dim.nprocy-1)*(dim.mxgrid-2*nghost_y*(dim.nprocy-1))*dim.mzgrid + $
                 2*nghost_z*(dim.nprocz-1)*(dim.mxgrid-2*nghost_x*(dim.nprocx-1))*(dim.mygrid-2*nghost_y*(dim.nprocy-1))
	correction = 1.0 - ghosts / double (dim.mxgrid*dim.mygrid*dim.mzgrid)
	file_struct = file_info (procdir+varfile)
	gb_per_file = (file_struct.size * subdomains * correction) / 1024. / 1024. / 1024.

	snapfiles = file_search (procdir, "VAR[0-9]*")
	num_snapshots = n_elements (snapfiles)
	for pos = 0, num_snapshots - 1 do begin
		snapfiles[pos] = strmid (snapfiles[pos], strpos (snapfiles[pos], "VAR"))
	end
	snapshots = strarr (1)
	for len = min (strlen (snapfiles)), max (strlen (snapfiles)) do begin
		indices = where (strlen (snapfiles) eq len)
		if (n_elements (indices) eq 1) then if (indices eq -1) then continue
		sub = snapfiles[indices]
		reorder = sort (sub)
		snapshots = [ snapshots, sub[reorder] ]
	end
	snapshots = snapshots[1:*]

	files_total = num_snapshots
	skipping = 0
	stepping = 1
	if (num_snapshots gt 0) then begin
		print, ""
		print, "============================================================================="
		print, "There are > ", strtrim (num_snapshots, 2), " < snapshot files available."
		print, "(This corresponds to ", strtrim (round (num_snapshots * gb_per_file * 10) / 10., 2), " GB.)"
		if ((stepping eq 1) and (skipping eq 0)) then begin
			print, "'"+procdir+varfile+"' will be read anyways."
			print, "Do you want to load additional files into the cache?"
			repeat begin
				answer = "n"
				read, answer, format="(A)", prompt="Read (A)ll / (S)elected / (N)o additional files: "
			end until (any (strcmp (answer, ['n', 'a', 's'], /fold_case)))
			if (strcmp (answer, 'n', /fold_case)) then begin
				stepping = 0
				skipping = num_snapshots
				files_total = 0
			end
			print, "Available snapshots: ", snapshots
			if (strcmp (answer, 's', /fold_case)) then begin
				if (num_snapshots gt 1) then begin
					print, "How many files do you want to skip at start?"
					repeat begin
						read, skipping, format="(I)", prompt="(0..."+strtrim (num_snapshots-1, 2)+"): "
					end until ((skipping ge 0) and (skipping le num_snapshots-1))
				end
				if ((num_snapshots-skipping) gt 1) then begin
					print, "Please enter a stepping for reading files:"
					print, "(0=do not read any more files, 1=each file, 2=every 2nd, ...)"
					repeat begin
						read, stepping, format="(I)", prompt="(0..."+strtrim (num_snapshots-skipping, 2)+"): "
					end until ((stepping ge 0) and (stepping le num_snapshots-skipping))
				end
				max_files = floor ((num_snapshots-1-skipping)/stepping) + 1
				if ((max_files gt 1) and (stepping ge 1)) then begin
					print, "How many files do you want to read in total?"
					repeat begin
						read, files_total, format="(I)", prompt="(0=all, 1..."+strtrim (max_files, 2)+"): "
					end until ((files_total ge 0) and (files_total le max_files))
					if (files_total eq 0) then files_total = max_files
					if (files_total lt max_files) then begin
						print, "You selected to load less files than availabe."
						addquestion = ""
						if (num_additional gt 0) then addquestion = " and '"+addfile+"'"
						print, "Do you want to skip reading of '"+varfile+"'"+addquestion+"?"
						repeat begin
							answer = "n"
							read, answer, format="(A)", prompt="(Y)es / (N)o: "
						end until (any (strcmp (answer, ['y', 'n'], /fold_case)))
						if (strcmp (answer, 'y', /fold_case)) then begin
							load_varfile = 0
							addfile = ""
							num_additional = 0
						end
					end
				end
				if (stepping eq 0) then begin
					files_total = 1
					stepping = 1
				end
			end
		end
	end

	num_selected = num_snapshots
	if (skipping ge 1) then num_selected -= skipping
	if (stepping ge 1) then num_selected = (num_selected - 1) / stepping + 1 else num_selected = 0
	if (num_selected lt 0) then num_selected = 0
	if (files_total gt num_selected) then files_total = num_selected
	if (files_total lt num_selected) then num_selected = files_total
	ignore_end = num_snapshots - skipping - (files_total * stepping)

	print, ""
	print, "Selected snapshots: skipping the first ", strtrim (skipping, 2), " with stepping=", strtrim (stepping, 2)
	print, "(This corresponds to ", strtrim (num_selected * gb_per_file, 2), " GB = ", strtrim (num_selected, 2), " files.)"
	print, ""


	pc_units, obj=unit, datadir=datadir
	units = { velocity:unit.velocity, time:unit.time, temperature:unit.temperature, length:unit.length, density:unit.density, mass:unit.density*unit.length^3, magnetic_field:unit.magnetic_field, default_length:default_length, default_time:default_time, default_velocity:default_velocity, default_density:default_density, default_mass:default_mass, default_magnetic_field:default_magnetic_field, default_length_str:default_length_str, default_time_str:default_time_str, default_velocity_str:default_velocity_str, default_density_str:default_density_str, default_mass_str:default_mass_str, default_magnetic_field_str:default_magnetic_field_str }
	pc_read_grid, obj=grid, dim=dim, datadir=datadir, allprocs=allprocs, /trim, /quiet
	pc_read_param, obj=param, dim=dim, datadir=datadir, /quiet
	pc_read_param, obj=run_param, /param2, dim=dim, datadir=datadir, /quiet

	coords = { x:congrid (grid.x, disp_size_x, 1, 1, /center, /interp)*unit.length/default_length, y:congrid (grid.y, disp_size_y, 1, 1, /center, /interp)*unit.length/default_length, z:congrid (grid.z, disp_size_z, 1, 1, /center, /interp)*unit.length/default_length, dx:congrid (1.0/grid.dx_1, disp_size_x, 1, 1, /center, /interp)*unit.length, dy:congrid (1.0/grid.dy_1, disp_size_y, 1, 1, /center, /interp)*unit.length, dz:congrid (1.0/grid.dz_1, disp_size_z, 1, 1, /center, /interp)*unit.length, nx:disp_size_x, ny:disp_size_y, nz:disp_size_z, l1:dim.nghostx, l2:dim.mx-dim.nghostx-1, m1:dim.nghosty, m2:dim.my-dim.nghosty-1, n1:dim.nghostz, n2:dim.mz-dim.nghostz-1 }

	if ((n_elements (dt) le 0) and file_test (datadir+"/time_series.dat")) then pc_read_ts, obj=ts, datadir=datadir, /quiet
	pc_show_ts, ts, units=units, param=param, run_param=run_param


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


	prepare_varset, load_varfile+num_additional+num_selected, units, coords, varset, overplot, datadir, param, run_param

	if (addfile) then begin
		; Precalculate additional timestep
		precalc, 0, varfile=addfile, datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs
	end

	if (load_varfile) then begin
		; Precalculate initial timestep
		precalc, num_additional, varfile=varfile, datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs
	end

	if (num_selected gt 0) then begin
		; Precalculate first selected timestep
		precalc, load_varfile+num_additional+num_selected-1, varfile=snapshots[skipping], datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs
		if (num_selected gt 1) then begin
			; Precalculate last selected timestep
			pos_last = skipping + (num_selected-1)*stepping
			precalc, load_varfile+num_additional, varfile=snapshots[pos_last], datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs
			if (num_selected gt 2) then begin
				for i = 2, num_selected-1 do begin
					; Precalculate selected timesteps
					pos = skipping + (i-1)*stepping
					precalc, load_varfile+num_additional+num_selected-i, varfile=snapshots[pos], datadir=datadir, dim=dim, param=param, run_param=run_param, varcontent=varcontent, allprocs=allprocs
				end
			end
		end
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

