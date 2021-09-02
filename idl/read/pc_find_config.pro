;+
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
;       pc_find_config, varfile, datadir=datadir, procdir=procdir, dim=dim, allprocs=allprocs, f77=f77, additional=additional, marker=marker, start_param=start_param
;
; MODIFICATION HISTORY:
;       $Id$
;       Written by: PABourdin, 06-Sep-2015
;
;-
pro pc_find_config, varfile, datadir=datadir, procdir=procdir, dim=dim, procdim=procdim, allprocs=allprocs, reduced=reduced, lcollect_xy=lcollect_xy, swap_endian=swap_endian, f77=f77, additional=additional, marker=marker, mvar_io=mvar_io, start_param=start_param, check_file=check_file, help=help

COMPILE_OPT IDL2,HIDDEN

	common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
	common cdat_coords, coord_system

        if (keyword_set(help)) then begin
          doc_library, 'pc_find_config'
          return
        endif

	; defaults
	default, reduced, 0
	datadir = pc_get_datadir(datadir)
	if (keyword_set (reduced)) then allprocs = 1
	default, check_file, 0
	dead_beef = 'DEADBEEF'

        pc_set_precision, dim=dim, datadir=datadir, /quiet
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
			pc_read_dim, object=procdim, datadir=datadir, proc=0, /quiet
			if (not is_defined(start_param)) then pc_read_param, object=start_param, dim=procdim, datadir=datadir, /quiet
			nprocxy = procdim.nprocx * procdim.nprocy
			num_var = procdim.mvar
			if (start_param.lwrite_aux) then num_var += procdim.maux
			if ((nprocxy le 1) or not file_test (datadir+'/proc1/') or (not file_test (datadir+'/proc1/'+varfile) and not file_test (datadir+'/proc1/VAR*'))) then begin
				fstat = file_info (procdir+varfile)
				data_end = long64(procdim.mxgrid)*long64(procdim.mygrid)*long64(procdim.mz)*long64(num_var*data_bytes)
				file_size = data_end + data_bytes * long64((1+1)+(procdim.mxgrid+procdim.mygrid+procdim.mzgrid+3+1))
				f2003_size = data_end + 2*strlen (dead_beef) + data_bytes * long64(1+3*(procdim.mxgrid+procdim.mygrid+procdim.mzgrid)+3)
				if (fstat.size lt data_end) then begin
					allprocs = 0
				end else if (fstat.size eq f2003_size) then begin
					allprocs = 3
					check_file = 1
				end else if ((fstat.size eq file_size) or (fstat.size eq file_size+8)) then begin
					; direct access or data record with F77 markers
					if file_test (datadir+'/param2.nml') then begin 
                                          	pc_read_param, obj=run_param, /run_param, /quiet
 					  	if (has_tag (run_param, 'lcollective_IO')) then begin
 							if run_param.lcollective_IO then allprocs = 2 else allprocs = 0
                                                endif
					endif else if (has_tag (start_param, 'lcollective_IO')) then begin
						if (start_param.lcollective_IO) then allprocs = 2 else allprocs = 0
					end else if (file_test (datadir+'/allprocs/grid.dat')) then begin
						allprocs = 2
					end else begin
						allprocs = 0
					end
				end 
			end
		end
	end else begin
		procdir = datadir+'/proc0/'
		if (allprocs eq 1) then procdir = datadir+'/allprocs/'
	end
	if (keyword_set (reduced) and (allprocs eq 1)) then procdir = datadir+'/reduced/'

	; get global dimensions
	if (size(dim, /type) ne 8) then $
		pc_read_dim, object=dim, datadir=datadir, reduced=reduced, /quiet

	; get "start.in" parameters
	if (not is_defined(start_param)) then pc_read_param, object=start_param, dim=dim, datadir=datadir, /quiet

	; set coordinate system
	coord_system = start_param.coord_system

	; get number of variables in VAR file
	mvar_io = dim.mvar
	if (start_param.lwrite_aux) then mvar_io += dim.maux

	; set also procdim
	if (not is_defined(procdim)) then begin
		if (allprocs eq 1) then procdim = dim else pc_read_dim, object=procdim, datadir=datadir, proc=0, /quiet
        endif
	; set other parameters according to allprocs
	marker = ''
	marker_bytes = 4
	if (allprocs eq 0) then begin
		f77 = 1
		additional = long64(procdim.mx)*long64(procdim.my)*long64(procdim.mz)*long64(mvar_io*data_bytes)+long64(2*marker_bytes)
	end else if (allprocs eq 1) then begin
		f77 = 1
		additional = long64(dim.mxgrid)*long64(dim.mygrid)*long64(dim.mzgrid)*long64(mvar_io*data_bytes)
	end else if (allprocs eq 2) then begin
		f77 = 1
		additional = long64(procdim.mxgrid)*long64(procdim.mygrid)*long64(procdim.mz)*long64(mvar_io*data_bytes)+long64(2*marker_bytes)
	end else if (allprocs eq 3) then begin
		marker = dead_beef
		f77 = 0
		additional = long64(procdim.mxgrid)*long64(procdim.mygrid)*long64(procdim.mz)*long64(mvar_io*data_bytes)+strlen(marker)
	end

	; set collective IO flags
	lcollect_xy = (allprocs ge 2) and (allprocs le 3)

	; check file integrity
	if (check_file) then begin
		t = zero
		x = make_array (dim.mx, type=type_idl)
		y = make_array (dim.my, type=type_idl)
		z = make_array (dim.mz, type=type_idl)
		dx = zero
		dy = zero
		dz = zero
		deltay = zero
		openr, file, procdir+varfile, f77=f77, swap_endian=swap_endian, /get_lun
		if (allprocs eq 0) then begin
			point_lun, file, additional
			if (start_param.lshear) then begin
				readu, file, t, x, y, z, dx, dy, dz, deltay
			end else begin
				readu, file, t, x, y, z, dx, dy, dz
			end
		end else if (allprocs eq 1) then begin
			point_lun, file, additional
			readu, file, t, x, y, z, dx, dy, dz
		end else if (allprocs eq 2) then begin
			point_lun, file, additional
			readu, file, t
		end else if (allprocs eq 3) then begin
			record_marker = marker
			point_lun, file, additional - strlen (marker)
			readu, file, record_marker
			if (marker ne record_marker) then message, "ERROR: '"+filename+"' has invalid record marker ('"+record_marker+"')."
		end
		close, file
		free_lun, file
	end
end
