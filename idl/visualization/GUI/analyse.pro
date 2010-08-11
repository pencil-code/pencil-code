;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   analyse.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Framework for precalculation and comparision of output in pencil units.
;;;   Calls cmp_cslice_cache for visualisation of full 3D data.
;;;
;;;  To do:
;;;   Add more comments

@analyse_companion


;;; Settings:

; Quantities to be visualized:
quantities = { temperature:'Temp', currentdensity:'j',            $
               magnetic_energy:'rho_mag', magnetic_field_z:'bz',  $
               velocity:'u_abs', velocity_z:'u_z',                $
               logarithmic_density:'ln_rho' }

; Available quantities are:
; 'Temp', 'rho', 'ln_rho'   ; temperature, density and logarithmic density
; 'ux', 'uy', 'uz', 'u_abs' ; velocity components and the absolute value
; 'Ax', 'Ay', 'Az'          ; vector potential components
; 'Bx', 'By', 'Bz'          ; magnetic field components
; 'rho_mag'                 ; magnetic energy density
; 'j'                       ; absolute value of the current density
; (more quantities can be defined in 'precalc_data', see analyse_companion.pro)

; Quantities to be overplotted (calculated in 'precalc_data'):
overplot_quantities = { magnetic_field:'b', velocities:'u' }

; Available quantities for overplotting are:
; 'b' and 'u'

; Preferred units for display
default_length        = 1.e6
default_length_str    = 'Mm'
default_velocity      = 1.e3
default_velocity_str  = 'km/s'

; initial varfile
default, varfile, 'var.dat'

; default data directory
default, datadir, './data/'

; stepping for varfiles
default, stepping, 1

; skipping of first n varfiles
default, skipping, 0


default, analyse_loaded, 0

if (not analyse_loaded) then BEGIN

	if (n_elements (nghost) gt 0) then begin
		nghost_x = nghost
		nghost_y = nghost
		nghost_z = nghost
	end
	pc_read_dim, obj=dim, /quiet
	default, nghost_x, dim.nghostx
	default, nghost_y, dim.nghosty
	default, nghost_z, dim.nghostz
	nx = dim.mx - 2*nghost_x
	ny = dim.my - 2*nghost_y
	nz = dim.mz - 2*nghost_z
	lmn12 = nghost_x+spread(indgen(nx),[1,2],[ny,nz]) + dim.mx*(nghost_y+spread(indgen(ny),[0,2],[nx,nz])) + dim.mx*dim.my*(nghost_z+spread(indgen(nz),[0,1],[nx,ny]))

	pc_units, obj=unit

	file_struct = file_info (datadir+"/proc0/var.dat")
	subdomains = n_elements (file_search (datadir, "../proc*"))
	gb_per_file = (file_struct.size * subdomains) / 1024. / 1024. / 1024.

	snapshots = file_search (datadir+"/proc0/", "VAR?")
	add = file_search (datadir+"/proc0/", "VAR??")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	add = file_search (datadir+"/proc0/", "VAR???")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	add = file_search (datadir+"/proc0/", "VAR????")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	add = file_search (datadir+"/proc0/", "VAR?????")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	num_snapshots = n_elements (snapshots)
	files_total = num_snapshots
	if (num_snapshots gt 0) then begin
		print, ""
		print, "There are ", num_snapshots, " snapshot files available."
		print, "(This corresponds to ", (round (num_snapshots * gb_per_file * 10) / 10.), " GB.)"
		if ((stepping eq 1) and (skipping eq 0)) then begin
			repeat begin
				print, "Which files do you want to load into the cache?"
				answer = "n"
				read, answer, format="(A)", prompt="Read (A)ll / (S)elected / (N)o files : "
			end until (any (strcmp (answer, ['n', 'a', 's'], /fold_case)))
			if (strcmp (answer, 'n', /fold_case)) then begin
				stepping = 0
				skipping = num_snapshots
				files_total = 0
			end
			if (strcmp (answer, 's', /fold_case)) then begin
				print, "Please enter a stepping for reading files:"
				print, "(0=do not read any files, 1=each file, 2=every 2nd, ...)"
				read, stepping, format="(I)", prompt="(0 - number of snapshots) : "
				if (stepping le 0) then begin
					stepping = 0
					skipping = num_snapshots
					files_total = 0
				end else begin
					print, "How many files do you want to skip at start ?"
					read, skipping, format="(I)", prompt="(0 - number of snapshots) : "
					print, "How many files do you want to read in total ?"
					read, files_total, format="(I)", prompt="(0 - number of remaining files) : "
				end
			end
		end
	end

	if (skipping lt 0) then skipping = 0
	if (files_total le 0) then begin
		files_total = 0
		stepping = 0
	end
	if (stepping lt 0) then stepping = 1

	num_selected = num_snapshots
	if (skipping ge 1) then num_selected -= skipping
	if (stepping ge 1) then num_selected = (num_selected - 1) / stepping + 1 else num_snapshots = 0
	if (num_selected lt 0) then num_selected = 0
	if (files_total gt num_selected) then files_total = num_selected
	if (files_total lt num_selected) then num_selected = files_total
	ignore_end = num_snapshots - skipping - (files_total * stepping)

	print, ""
	print, "Selected snapshots: skipping the first ", skipping, " with stepping=", stepping
	print, "(This corresponds to ", (num_selected * gb_per_file), " GB = ", num_selected, " files.)"
	print, ""


	time_series = file_search (datadir, "time_series.dat")
	if ((n_elements (dt) le 0) and (strlen (time_series[0]) gt 0)) then pc_read_ts, obj=ts, /quiet
	if (n_elements (ts) gt 0) then begin
		window, 1, xsize=1000, ysize=800, title = 'time series analysis', retain=2
		!P.MULTI = [0, 2, 2]

		tags = tag_names (ts)
		y_minmax = minmax (ts.dt)
		if (any (strcmp (tags, 'dtu', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtu])
		if (any (strcmp (tags, 'dtv', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtv])
		if (any (strcmp (tags, 'dtnu', /fold_case)))   then y_minmax = minmax ([y_minmax, ts.dtnu])
		if (any (strcmp (tags, 'dtb', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtb])
		if (any (strcmp (tags, 'dteta', /fold_case)))  then y_minmax = minmax ([y_minmax, ts.dteta])
		if (any (strcmp (tags, 'dtc', /fold_case)))    then y_minmax = minmax ([y_minmax, ts.dtc])
		if (any (strcmp (tags, 'dtchi', /fold_case)))  then y_minmax = minmax ([y_minmax, ts.dtchi])
		if (any (strcmp (tags, 'dtchi2', /fold_case))) then y_minmax = minmax ([y_minmax, ts.dtchi2])

		print, "starting values:"
		print, "dt    :", ts.dt[0]
		plot, ts.dt, title = 'dt', yrange=y_minmax, /yl
		plot, ts.t, ts.dt, title = 'dt(tt) u{-t} v{-p} nu{.v} b{.r} eta{-g} c{.y} chi{-.b} chi2{-.o} [s]', yrange=y_minmax, /yl
		if (any (strcmp (tags, 'dtu', /fold_case))) then begin
			oplot, ts.t, ts.dtu, linestyle=2, color=11061000
			print, "dtu   :", ts.dtu[0]
		end
		if (any (strcmp (tags, 'dtv', /fold_case))) then begin
			oplot, ts.t, ts.dtv, linestyle=2, color=128255200
			print, "dtv   :", ts.dtv[0]
		end
		if (any (strcmp (tags, 'dtnu', /fold_case))) then begin
			oplot, ts.t, ts.dtnu, linestyle=1, color=128000128
			print, "dtnu  :", ts.dtnu[0]
		end
		if (any (strcmp (tags, 'dtb', /fold_case))) then begin
			oplot, ts.t, ts.dtb, linestyle=1, color=200
			print, "dtb   :", ts.dtb[0]
		end
		if (any (strcmp (tags, 'dteta', /fold_case))) then begin
			oplot, ts.t, ts.dteta, linestyle=2, color=220200200
			print, "dteta :", ts.dteta[0]
		end
		if (any (strcmp (tags, 'dtc', /fold_case))) then begin
			oplot, ts.t, ts.dtc, linestyle=1, color=61695
			print, "dtc   :", ts.dtc[0]
		end
		if (any (strcmp (tags, 'dtchi', /fold_case))) then begin
			oplot, ts.t, ts.dtchi, linestyle=3, color=115100200
			print, "dtchi :", ts.dtchi[0]
		end
		if (any (strcmp (tags, 'dtchi2', /fold_case))) then begin
			oplot, ts.t, ts.dtchi2, linestyle=3, color=41215
			print, "dtchi2:", ts.dtchi2[0]
		end
		max_subplots = 2
		num_subplots = 0
		if (any (strcmp (tags, 'TTmax', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			Temp_max = ts.TTmax * unit.temperature
			plot, ts.t, Temp_max, title = 'Temp_max(tt) [K]', /yl
		end
		if (any (strcmp (tags, 'umax', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			u_max = ts.umax * unit.velocity / default_velocity
			plot, ts.t, u_max, title = 'u_max(tt) ['+default_velocity_str+']'
		end
		if (any (strcmp (tags, 'rhomin', /fold_case)) and (num_subplots lt max_subplots)) then begin
			num_subplots += 1
			rho_min = ts.rhomin * unit.density
			plot, ts.t, rho_min, title = 'rho_min(tt)', /yl
		end
	end

	resolve_routine, "cmp_cslice_cache", /COMPILE_FULL_FILE, /NO_RECOMPILE

	units = { velocity:unit.velocity, temperature:unit.temperature, length:unit.length, density:unit.density, default_length:default_length, default_velocity:default_velocity, default_length_str:default_length_str, default_velocity_str:default_velocity_str }
	pc_read_grid, obj=grid, /trim, /quiet
	coords = { x:grid.x/default_length, y:grid.y/default_length, z:grid.z/default_length }
	dummy = dindgen (dim.mx, dim.my, dim.mz)
	dummy_3D = findgen (dim.mx, dim.my, dim.mz, 3)

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

	prepare_varset, num_selected+1, units, coords, varset, overplot

	; Precalculate initial timestep
	precalc, 0, varfile=varfile

	if (num_selected gt 0) then begin
		for i = 1, num_selected do begin
			; Precalculate selected timesteps
			pos = (num_snapshots - ignore_end) - (i - 1) * stepping - 1
			precalc, i, varfile=strmid (snapshots[pos], strpos (snapshots[pos], "VAR"))
		end
	end

	; Mark completition of preparational work
	analyse_loaded = 1

END

; scaling factor for visualisation
default, scaling, fix (256 / max ([dim.nx, dim.ny, dim.nz]))
if (n_elements (scaling) eq 1) then if (scaling lt 1) then scaling = 1


cmp_cslice_cache, quantities, lmn12, scaling=scaling, overplots=overplot_quantities

window, 0, xsize=8, ysize=8, retain=2
!P.MULTI = [0, 1, 1]
wdelete

end

