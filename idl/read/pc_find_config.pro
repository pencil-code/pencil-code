; +
; NAME:
;       PC_FIND_CONFIG
;
; PURPOSE:
;       Automatically find VAR files and set configuration parameters like allprocs.
;
; CATEGORY:
;       Pencil Code, File I/O
;
; CALLING SEQUENCE:
;       pc_find_config, varfile, datadir=datadir, procdir=procdir, dim=dim, allprocs=allprocs, reduced=reduced, swap_endian=swap_endian, f77=f77, additional=additional, marker=marker, mvar_io=mvar_io, start_param=param
;
; KEYWORD PARAMETERS:
;     datadir: Specifies the root data directory. Default: './data'.   [string]
;     procdir: Name of the var file directory. Default: '/proc0/'.     [string]
;         dim: Dimensions structure.                                   [structure]
;    allprocs: Load data from the allprocs directory.                  [integer]
;  additional: Additional data block position (like time, grid, etc.). [long]
;      marker: Record marker for F2003 stream access files.            [string]
;     mvar_io: Number of IO variables.                                 [integer]
; start_param: Start parameters from "start.in".                       [structure]
;         f77: F77 record marker activation parameter inside VAR files.
;
;    /reduced: Load reduced collective varfiles.
; /check_file: Check VAR file for consistency.
;/swap_endian: Swap endianess while reading.
;
; ALLPROCS VALUES:
;       0 : "data/proc[0:N-1]/"        distributed VAR files from "io_dist"
;       1 : "data/allprocs/"           collective monolithic VAR file from "io_collect" or "pc_collect"
;       2 : "data/proc[0:NZ-1]*NX*NY/" xy-collective VAR files from "io_collect_xy"
;       3 : "data/proc[0:NZ-1]*NX*NY/" xy-collective VAR files from "io_collect_xy_f2003"
;
; EXAMPLE:
;       Automatically find VAR files and set all unset parameters:
;       pc_find_config, varfile, datadir=datadir, procdir=procdir, dim=dim, allprocs=allprocs, f77=f77, additional=additional, marker=marker, mvar_io=mvar_io, start_param=param
;
; MODIFICATION HISTORY:
;       $Id$
;       Written by: PABourdin, 06-Sep-2015
;
;-
pro pc_find_config, varfile, datadir=datadir, procdir=procdir, dim=dim, allprocs=allprocs, reduced=reduced, swap_endian=swap_endian, f77=f77, additional=additional, marker=marker, mvar_io=mvar_io, start_param=start_param, check_file=check_file

COMPILE_OPT IDL2,HIDDEN

	common pc_precision, zero, one

	; defaults
	if (not keyword_set (datadir)) then datadir = pc_get_datadir()
	if (keyword_set (reduced)) then allprocs = 1
	default, check_file, 0
	dead_beef = 'DEADBEEF'

	; determine file type (set allprocs and procdir)
	if (keyword_set (reduced)) then begin
		if (n_elements (allprocs) gt 0) then begin
			if (allprocs ne 1) then begin
				print, "ERROR: if 'reduced' is set, 'allprocs=1' must be chosen."
				print, "Type '.c' to continue..."
				stop
			end
		end
		allprocs = 1
	end
	if (size (allprocs, /type) eq 0) then begin
		procdir = datadir+'/allprocs/'
		if (file_test (procdir+varfile) or file_test (procdir+'VAR*')) then begin
			allprocs = 1
		end else begin
			procdir = datadir+'/proc0/'
			allprocs = 0
			pc_read_dim, object=pd, datadir=datadir, proc=0, /quiet
			if (size (start_param, /type) ne 8) then pc_read_param, object=start_param, dim=pd, datadir=datadir, /quiet
			nprocxy = pd.nprocx * pd.nprocy
			num_var = pd.mvar
			if (start_param.lwrite_aux) then num_var += pd.maux
			if ((nprocxy le 1) or not file_test (datadir+'/proc1/') or (not file_test (datadir+'/proc1/'+varfile) and not file_test (datadir+'/proc1/VAR*'))) then begin
				fstat = file_info (procdir+varfile)
				if (pd.precision eq 'D') then real_bytes=8 else real_bytes=4
				data_end = long64(pd.mxgrid)*long64(pd.mygrid)*long64(pd.mz)*long64(num_var*real_bytes)
				file_size = data_end + real_bytes * long64((1+1)+(pd.mxgrid+pd.mygrid+pd.mzgrid+3+1))
				f2003_size = data_end + 2*strlen (dead_beef) + real_bytes * long64(1+3*(pd.mxgrid+pd.mygrid+pd.mzgrid)+3)
				if (fstat.size lt data_end) then begin
					allprocs = 0
				end else if (fstat.size eq f2003_size) then begin
					allprocs = 3
					check_file = 1
				end else if ((fstat.size eq file_size) or (fstat.size eq file_size+8)) then begin
					; direct access or data record with F77 markers
					allprocs = 2
				end 
			end
		end
	end else begin
		procdir = datadir+'/proc0/'
		if (allprocs eq 1) then procdir = datadir+'/allprocs/'
	end
	if (keyword_set (reduced) and (allprocs eq 1)) then procdir = datadir+'/reduced/'

	; get dimensions
	if (size (dim, /type) ne 8) then begin
		if (allprocs eq 1) then begin
			pc_read_dim, object=dim, datadir=datadir, reduced=reduced, /quiet
		end else if (allprocs ge 2) then begin
			pc_read_dim, object=dim, datadir=datadir, proc=0, /quiet
			dim.nx = dim.nxgrid
			dim.ny = dim.nygrid
			dim.mx = dim.mxgrid
			dim.my = dim.mygrid
			dim.mw = dim.mx * dim.my * dim.mz
		end else begin
			pc_read_dim, object=dim, datadir=datadir, proc=0, /quiet
		end
	end

	; get "start.in" parameters
	if (size (start_param, /type) ne 8) then pc_read_param, object=start_param, dim=dim, datadir=datadir, /quiet

	; get number of variables in VAR file
	mvar_io = dim.mvar
	if (start_param.lwrite_aux) then mvar_io += dim.maux

	; set precision
	pc_set_precision, dim=dim, /quiet
	if (dim.precision eq 'D') then real_bytes=8 else real_bytes=4

	; set other parameters according to allprocs
	marker = ''
	marker_bytes = 4
	if (allprocs eq 0) then begin
		f77 = 1
		additional = long64(dim.mx)*long64(dim.my)*long64(dim.mz)*long64(mvar_io*real_bytes)+long64(2*marker_bytes)
	end else if (allprocs eq 1) then begin
		f77 = 1
		additional = long64(dim.mxgrid)*long64(dim.mygrid)*long64(dim.mzgrid)*long64(mvar_io*real_bytes)
	end else if (allprocs eq 2) then begin
		f77 = 1
		additional = long64(dim.mxgrid)*long64(dim.mygrid)*long64(dim.mz)*long64(mvar_io*real_bytes)+long64(2*marker_bytes)
	end else if (allprocs eq 3) then begin
		marker = dead_beef
		f77 = 0
		additional = long64(dim.mxgrid)*long64(dim.mygrid)*long64(dim.mz)*long64(mvar_io*real_bytes)
		; check file integrity
		record_marker = marker
		if (check_file) then begin
			openr, file, procdir+varfile, swap_endian=swap_endian, /get_lun
			point_lun, file, additional
			readu, file, record_marker
			if (marker ne record_marker) then message, "ERROR: '"+filename+"' has invalid record marker ('"+record_marker+"')."
			close, file
			free_lun, file
		end
		additional += strlen (marker)
	end
end
