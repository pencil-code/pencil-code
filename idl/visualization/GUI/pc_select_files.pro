;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_select_files.pro      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   GUI for selecting a subset of available snapshot files.
;;;   The returned array contains a list of selected filenames.
;;;   Optional parameters are:
;;;   * pattern (contains the file-search pattern, default: "VAR[0-9]*")
;;;   * varfile (contains the varfile loaded by default, default: "var.dat")
;;;   * addfile (contains an additional filename, default: "crash.dat")
;;;   * datadir (contains the datadir, default: pc_get_datadir)
;;;   * allprocs (contains the IO-strategy parameter, default: automatic)
;;;   * procdir (contains procdir based on the chosen IO-strategy)
;;;   If an optional parameter is given as undefined, the default is returned.
;;;
;;;   Examples:
;;;   IDL> pc_select_files, obj=files
;;;   IDL> pc_select_files, obj=files, pattern="PVAR1[5-9]*"
;;;   IDL> pc_select_files, obj=files, varfile="myvar.dat", datadir=datadir
;;;   IDL> pc_select_files, obj=files, datadir="mydata", procdir=procdir
;;;   IDL> pc_select_files, obj=files, addfile="crash.dat", /allprocs
;;;   IDL> pc_select_files, obj=files, pattern="VAR[1-9]*", allprocs=allprocs
;;;


; Event handling of file dialog window
pro select_files_event, event

	common select_files_gui_common, b_var, b_add, c_list, i_skip, i_step, f_gb
	common select_files_common, num_files, selected, num_selected, var_selected, add_selected, skipping, stepping, data_dir, units, run_par, start_par, gb_per_file

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1

	SWITCH eventval of
	'ALL': begin
		WIDGET_CONTROL, c_list, SET_LIST_SELECT = indgen (num_files)
		num_selected = num_files
		WIDGET_CONTROL, f_gb, SET_VALUE = gb_per_file*(num_selected+var_selected+add_selected)
		break
	end
	'NONE': begin
		WIDGET_CONTROL, c_list, SET_LIST_SELECT = -1
		num_selected = 0
		WIDGET_CONTROL, f_gb, SET_VALUE = gb_per_file*(var_selected+add_selected)
		break
	end
	'SKIP': begin
		skipping = event.value
		if (skipping ge num_files) then begin
			skipping = num_files - 1
			WIDGET_CONTROL, i_skip, SET_VALUE = skipping
		end
		if (stepping gt num_files - skipping) then begin
			stepping = num_files - skipping
			WIDGET_CONTROL, i_step, SET_VALUE = stepping
		end
		break
	end
	'STEP': begin
		stepping = event.value
		if (stepping gt num_files - skipping) then begin
			stepping = num_files - skipping
			WIDGET_CONTROL, i_step, SET_VALUE = stepping
		end
		break
	end
	'APPLY': begin
		selected = indgen (num_files)
		selected = selected[skipping:*:stepping]
		WIDGET_CONTROL, c_list, SET_LIST_SELECT = selected
	end
	'LIST': begin
		selected = WIDGET_INFO (c_list, /LIST_SELECT)
		if (any (selected ne -1)) then num_selected = n_elements (selected) else num_selected = 0
		WIDGET_CONTROL, f_gb, SET_VALUE = gb_per_file*(num_selected+var_selected+add_selected)
		break
	end
	'ADD': begin
		add_selected = event.select
		WIDGET_CONTROL, f_gb, SET_VALUE = gb_per_file*(num_selected+var_selected+add_selected)
		break
	end
	'VAR': begin
		var_selected = event.select
		WIDGET_CONTROL, f_gb, SET_VALUE = gb_per_file*(num_selected+var_selected+add_selected)
		break
	end
	'SHOW_TIME': begin
		pc_show_ts, obj=ts, units=units, param=start_par, run_param=run_par, datadir=data_dir
		break
	end
	'OK': begin
		selected = WIDGET_INFO (c_list, /LIST_SELECT)
		quit = event.top
		break
	end
	'CANCEL': begin
		selected = -1
		add_selected = 0
		var_selected = 0
		quit = event.top
		break
	end
	endswitch

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


pro pc_select_files, files=files, pattern=pattern, varfile=varfile, addfile=addfile, datadir=datadir, allprocs=allprocs, procdir=procdir, units=units_struct, param=param, run_param=run_param

	common select_files_gui_common, b_var, b_add, c_list, i_skip, i_step, f_gb
	common select_files_common, num_files, selected, num_selected, var_selected, add_selected, skipping, stepping, data_dir, units, run_par, start_par, gb_per_file

	; Default settings
	if (not keyword_set (pattern)) then pattern = "VAR[0-9]*"
	if (not keyword_set (varfile)) then varfile = "var.dat"
	if (not keyword_set (addfile)) then addfile = "crash.dat"
	if (not keyword_set (datadir)) then datadir = pc_get_datadir ()
	default, skipping, 1
	default, stepping, 2

	if (datadir eq "") then datadir = "."
	data_dir = datadir
	if (keyword_set (units_struct)) then units = units_struct
	if (n_elements (units) le 0) then pc_units, obj=units, datadir=datadir
	if (keyword_set (param)) then start_par = param
	if (keyword_set (run_param)) then run_par = run_param

	if (n_elements (allprocs) eq 0) then begin
		allprocs = 1
		procdir = datadir+"/allprocs/"
		if (not file_test (procdir+varfile)) then begin
			allprocs = 0
			procdir = datadir+"/proc0/"
			if (file_test (datadir+"/proc1/")) then begin
				if (not file_test (datadir+"/proc1/"+varfile)) then allprocs = 2
			endif
		endif
	end else begin
		procdir = datadir+"/proc0/"
		if (allprocs eq 1) then procdir = datadir+"/allprocs/"
	end

	file_struct = file_info (procdir+varfile)
	gb_per_file = file_struct.size / 1073741824.

	files = file_search (procdir, pattern)
	num_files = n_elements (files)
	for pos = 0, num_files - 1 do begin
		files[pos] = strmid (files[pos], strpos (procdir, "/", /REVERSE_SEARCH) - 1)
	end
	sorted = [ "" ]
	for len = min (strlen (files)), max (strlen (files)) do begin
		indices = where (strlen (files) eq len)
		if (n_elements (indices) eq 1) then if (indices eq -1) then continue
		sub = files[indices]
		reorder = sort (sub)
		sorted = [ sorted, sub[reorder] ]
	end
	files = sorted[1:*]
	sorted = 1

	selected = -1
	num_selected = 0
	if ((addfile ne "") and file_test (procdir+addfile)) then add_selected = 1 else add_selected = 0
	if ((varfile ne "") and file_test (procdir+varfile)) then var_selected = 1 else var_selected = 0


	MOTHER	= WIDGET_BASE (title='PC file-dialog', EVENT_PRO=select_files_event)
	BASE	= WIDGET_BASE (MOTHER, /row)

	CTRL	= WIDGET_BASE (BASE, /col)

	tmp	= WIDGET_LABEL (CTRL, value='File Selection:', frame=0)
	EDIT	= WIDGET_BASE (CTRL, /row, /align_center)
	tmp	= WIDGET_BUTTON (EDIT, xsize=60, value='ALL', uvalue='ALL')
	tmp	= WIDGET_BUTTON (EDIT, xsize=60, value='NONE', uvalue='NONE')

	EDIT	= WIDGET_BASE (CTRL, /col, xsize=100, /align_center, frame=1)
	i_skip	= CW_FIELD (EDIT, title='skip first', /column, uvalue='SKIP', value=skipping, /integer, /return_events, xsize=8)
	i_step	= CW_FIELD (EDIT, title='stepping', /column, uvalue='STEP', value=stepping, /integer, /return_events, xsize=8)
	tmp	= WIDGET_BUTTON (EDIT, xsize=80, value='APPLY', uvalue='APPLY')

	XTRA	= WIDGET_BASE (CTRL, /col, /align_center)
	tmp	= WIDGET_LABEL (XTRA, value='Analysis:', frame=0)
	BUT	= WIDGET_BASE (XTRA, /row, /align_center)
	timeser	= WIDGET_BUTTON (BUT, value='show timeseries', uvalue='SHOW_TIME')

	tmp	= CW_FIELD (XTRA, title='GB per file', /column, value=gb_per_file, /float)
	f_gb	= CW_FIELD (XTRA, title='Total GB selected', /column, value=gb_per_file*(add_selected+var_selected), /float)

	BUT	= WIDGET_BASE (XTRA, /row, /align_center, frame=1)
	tmp	= WIDGET_BUTTON (BUT, xsize=60, value='CANCEL', uvalue='CANCEL')
	tmp	= WIDGET_BUTTON (BUT, xsize=60, value='OK', uvalue='OK')

	SEL	= WIDGET_BASE (BASE, /col)

	tmp	= WIDGET_LABEL (SEL, value='Available snapshots:', frame=0)
	c_list	= WIDGET_LIST (SEL, value=files, uvalue='LIST', YSIZE=num_files<24, /multiple)

	if (var_selected eq 1) then begin
		pc_read_var_time, t=var_time, varfile=varfile, datadir=datadir, allprocs=allprocs
		b_var	= CW_BGROUP (SEL, varfile+' ('+strtrim (var_time*units.time, 2)+' s)', /nonexcl, uvalue='VAR', set_value=1)
	end else if (varfile ne "") then begin
		b_var	= WIDGET_LABEL (SEL, value='"'+varfile+'" not found', frame=0)
	end else begin
		b_var	= WIDGET_LABEL (SEL, value='No varfile selected', frame=0)
	end

	if (addfile eq varfile) then begin
		b_add	= WIDGET_LABEL (SEL, value='Additional file is identical to "'+varfile+'"', frame=0)
	end else if (add_selected eq 1) then begin
		pc_read_var_time, t=add_time, varfile=addfile, datadir=datadir, allprocs=allprocs
		b_add	= CW_BGROUP (SEL, addfile+' ('+strtrim (add_time*units.time, 2)+' s)', /nonexcl, uvalue='ADD', set_value=1)
	end else begin
		b_add	= WIDGET_LABEL (SEL, value='No "'+addfile+'" found', frame=0)
	end

	WIDGET_CONTROL, MOTHER, /REALIZE
	wimg = !d.window

	WIDGET_CONTROL, BASE

	XMANAGER, "select_files", MOTHER


	if (any (selected ne -1)) then files = files[selected] else files = [ "" ]
	if (var_selected) then files = [ files, varfile ]
	if (add_selected) then files = [ files, addfile ]
	indices = where (files ne "")
	if (any (indices ne -1)) then begin
		files = files[indices]
		num_selected = n_elements (files)
	end else begin
		files = -1
		num_selected = 0
	end

	print, ""
	print, "Selected snapshots:"
	if (num_selected le 0) then begin
		print, "none"
	end else begin
		print, files
		print, "This corresponds to ", strtrim (num_selected * gb_per_file, 2), " GB = ", strtrim (num_selected, 2), " files"
	end
	print, ""

end

