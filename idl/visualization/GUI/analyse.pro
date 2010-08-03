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

var_source = ['uu', 'lnrho', 'lnTT', 'aa']

; Quantities to be visualized (calculated in 'precalc_data'):
quantities = { temperature:'Temp', currentdensity:'j',            $
               magnetic_energy:'rho_mag', magnetic_field_z:'bz',  $
               velocity:'u_abs', velocity_z:'u_z',                $
               logarithmic_density:'ln_rho' }

; Quantities to be overplotted (calculated in 'precalc_data'):
overplot_quantities = { magnetic_field:'b', velocities:'u' }

; stepping for varfiles
default, stepping, 1

; skipping of first n varfiles
default, skipping, 0

; scaling factor for visualisation
default, scaling, fix (256 / max ([nx, ny, nz]))
if (n_elements (scaling) eq 1) then if (scaling lt 1) then scaling = 1


default, analyse_loaded, 0

if (not analyse_loaded) then BEGIN

	lmn12 = l1+spread(indgen(nx),[1,2],[ny,nz]) + mx*(m1+spread(indgen(ny),[0,2],[nx,nz])) + mx*my*(n1+spread(indgen(nz),[0,1],[nx,ny]))

	time_series = file_search (datadir+"/../", "time_series.dat")
	if ((n_elements (dt) le 0) and (strlen (time_series[0]) gt 0)) then begin
		@time_series

	default, tt, [0, 0]
	default, dt, [0, 0]
	default, umax, [0, 0]
	default, TTmax, [0, 0]
	default, rhomin, [0, 0]

	u_max = umax * unit_velocity * 1e-3
	Temp_max = TTmax * unit_temperature
	rho_min = rhomin * unit_density


	file_struct = file_info (datadir+"/var.dat")
	subdomains = n_elements (file_search (datadir, "../proc*"))
	gb_per_file = (file_struct.size * subdomains) / 1024. / 1024. / 1024.

	snapshots = file_search (datadir, "VAR?")
	add = file_search (datadir, "VAR??")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	add = file_search (datadir, "VAR???")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	add = file_search (datadir, "VAR????")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	add = file_search (datadir, "VAR?????")
	if (strlen (add[0]) gt 0) then snapshots = [ snapshots, add ]
	num_snapshots = n_elements (snapshots)
	files_total = num_snapshots
	if (num_snapshots gt 0) then begin
		print, ""
		print, "There are ", num_snapshots, " snapshot files available."
		print, "(This corresponds to ", (round (num_snapshots * gb_per_file * 10) / 10.), " GB.)"
		if ((stepping eq 1) and (skipping eq 0)) then begin
			print, "Do you want to read all those files into the cache ?"
			answer = "-"
			read, answer, format="(A)", prompt="(y/n) : "
			if (answer ne "y") then begin
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


	window, 1, xsize=1000, ysize=800, title = 'time series analysis', retain=2
	!P.MULTI = [0, 2, 2]

	default, dtu, dt[0:1]
	default, dtnu, dt[0:1]
	default, dtb, dt[0:1]
	default, dteta, dt[0:1]
	default, dtc, dt[0:1]
	default, dtchi, dt[0:1]
	default, dtchi2, dt[0:1]

	print, "starting values:"
	print, "dt    :", dt[0]
	print, "dtu   :", dtu[0]
	print, "dtnu  :", dtnu[0]
	print, "dtb   :", dtb[0]
	print, "dteta :", dteta[0]
	print, "dtc   :", dtc[0]
	print, "dtchi :", dtchi[0]
	print, "dtchi2:", dtchi2[0]

	plot, dt, title = 'dt', yrange=[min([dt,dtb,dteta,dtchi]),max([dt,dtb,dteta,dtchi])], /yl
	plot, tt, dt, title = 'dt(tt) u{-t} nu{.v} b{.r} eta{-g} c{.y} chi{-.b} chi2{-.o} [s]', yrange=[min([dt,dtu,dtnu,dtb,dteta,dtc,dtchi,dtchi2]),max([dt,dtu,dtnu,dtb,dteta,dtc,dtchi,dtchi2])], /yl
	oplot, tt, dtu, linestyle=2, color=11061000
	oplot, tt, dtnu, linestyle=1, color=128000128
	oplot, tt, dtb, linestyle=1, color=200
	oplot, tt, dteta, linestyle=2, color=220200200
	oplot, tt, dtc, linestyle=1, color=61695
	oplot, tt, dtchi, linestyle=3, color=115100200
	oplot, tt, dtchi2, linestyle=3, color=41215
	plot, tt, Temp_max, title = 'Temp_max(tt) [K]', /yl
;	plot, tt, u_max, title = 'u_max(tt)'
	plot, tt, rho_min, title = 'rho_min(tt)', /yl


	resolve_routine, "cmp_cslice_cache", /COMPILE_FULL_FILE, /NO_RECOMPILE

	units = { velocity:unit_velocity, temperature:unit_temperature, length:unit_length, density:unit_density }
	Mm_SI = 1.e6
	coords = { x:reform (xx[l1:l2,m1,n1])/Mm_SI, y:reform (yy[l1,m1:m2,n1])/Mm_SI, z:reform (zz[l1,m1,n1:n2])/Mm_SI }
	dummy = dindgen (mx, my, mz)
	dummy_3D = findgen (mx, my, mz, 3)

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

	prepare_varset, num_selected+1, units, coords, varset, overplot, var_source

	; Create variable source structure
	exec_str = "vars = { "
	for i = 0, n_elements (var_source) - 1 do begin
		if (i gt 0) then exec_str += ", "
		exec_str += var_source[i]+":"+var_source[i]
	end
	exec_str += " }"
	res = execute (exec_str)
	if (not res) then begin
		print, "Could not create variable source structure!"
		stop
	end

	; Precalculate initial timestep
	precalc, 0, varfile='var.dat', vars=vars

	if (num_selected gt 0) then begin
		for i = 1, num_selected do begin
			; Precalculate selected timesteps
			pos = (num_snapshots - ignore_end) - (i - 1) * stepping - 1
			precalc, i, show_aver=1, varfile=strmid (snapshots[pos], strpos (snapshots[pos], "VAR"))
		end
	end

	; Mark completition of preparational work
	analyse_loaded = 1

END

cmp_cslice_cache, quantities, lmn12, scaling=scaling, overplots=overplot_quantities

window, 0, xsize=8, ysize=8, retain=2
!P.MULTI = [0, 1, 1]
wdelete

end

