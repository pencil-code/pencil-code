; $Id$
;
; NAME:
;       PC_WRITE_VAR_RAW
;
; PURPOSE:
;       Write a 'var.dat', or other VAR files in an efficient way!
;
;       Needs one array to write into a snapshot (var) file usable by
;       a Pencil Code run, and a structure containing the variable labels.
;       Works for distributed and collective IO modules.
;       Persistent variables need an additional treatment and should anyways be
;       stored separately setting the flag "lseparate_persist=T" in 'start.in'.
;
; CATEGORY:
;       Pencil Code, File I/O
;
; CALLING SEQUENCE:
;       pc_write_var_raw, object, varfile=varfile, datadir=datadir, dim=dim, /allprocs, /quiet
; KEYWORD PARAMETERS:
;     object: Object in which the data is given.                     [4D-array]
;       time: Time of the snapshot. Default: 0.0d0.                  [float or double]
;    datadir: Specifies the root data directory. Default: './data'.  [string]
;        dim: Dimensions structure.                                  [structure]
;      param: Parameters structure.                                  [structure]
;    varfile: Name of the var file. Default: 'var.dat'.              [string]
;   allprocs: Write distributed (0) or collective (1 or 2) varfiles. [integer]
;     /quiet: Suppress any information messages and summary statistics.
;
; EXAMPLES:
;       pc_read_var_raw, obj=vars, tags=tags, varfile='VAR123'     ;; read from data/proc*/VAR123
;       vars[*,*,*,tags.ux:tags:uz] = 0.0                          ;; alter some data
;       pc_write_var_raw, vars, varfile='VAR0'                     ;; write with t=0 to data/proc*/VAR0
;
;       pc_read_var_raw, obj=vars, tags=tags, time=time, /allprocs ;; read from data/allprocs/var.dat
;       vars[*,*,*,tags.ux:tags:uz] = 0.0                          ;; alter some data
;       pc_write_var_raw, vars, time=time, /allprocs               ;; write to data/allprocs/var.dat
;
; MODIFICATION HISTORY:
;       Adapted from: pc_read_var_raw.pro, 6th July 2013 (Bourdin.KIS)

pro pc_write_var_raw, vars, time=time, dim=dim, grid=grid, param=param, datadir=datadir, varfile=varfile, allprocs=allprocs, f77=f77, swap_endian=swap_endian, quiet=quiet

	common pc_precision, zero, one

	default, datadir, './data'
	default, varfile, 'var.dat'

	if (not keyword_set (dim)) then pc_read_dim, obj=dim, datadir=datadir
	if (not keyword_set (param)) then pc_read_param, obj=param, dim=dim, datadir=datadir, /quiet
	if (not keyword_set (grid)) then pc_read_grid, obj=grid, dim=dim, param=param, datadir=datadir, /quiet
	if (size (vars, /type) eq 8) then begin
		; Need to have a valid varcontent
		varcontent = pc_varcontent (datadir=datadir, dim=dim, param=param)
		; Create array out of given structure and pass recursively computed results
		vars = pc_convert_vars_struct (vars, varcontent, index)
	end

	default, time, zero
	default, allprocs, 0
	if (keyword_set (allprocs)) then default, f77, 0
	default, f77, 1

	nx = dim.nx
	ny = dim.ny
	nz = dim.nz
	nprocx = dim.nprocx
	nprocy = dim.nprocy
	nprocz = dim.nprocz
	nghostx = dim.nghostx
	nghosty = dim.nghosty
	nghostz = dim.nghostz

	if (allprocs ge 1) then begin
		last_proc_x = 0
		last_proc_y = 0
	end else begin
		last_proc_x = nprocx - 1
		last_proc_y = nprocy - 1
	end
	if (allprocs eq 1) then begin
		last_proc_z = 0
	end else begin
		last_proc_z = nprocz - 1
	end

	for ipx = 0, last_proc_x do begin
		for ipy = 0, last_proc_y do begin
			for ipz = 0, last_proc_z do begin

				iproc = ipx + ipy*nprocx + ipz*nprocx*nprocy
				if (allprocs eq 1) then begin
					procdir = 'allprocs'
				end else begin
					procdir = 'proc' + strtrim (string (iproc), 2)
				end
				filename = datadir + '/' + procdir + '/' + varfile

				if (allprocs ge 1) then begin
					xs = 0
					xe = dim.mx - 1
					ys = 0
					ye = dim.my - 1
				end else begin
					pnx = nx / nprocx
					pny = ny / nprocy
					xs = ipx * pnx
					xe = (ipx+1) * pnx + 2*nghostx - 1
					ys = ipy * pny
					ye = (ipy+1) * pny + 2*nghosty - 1
				end
				if (allprocs eq 1) then begin
					zs = 0
					ze = dim.mz - 1
				end else begin
					pnz = nz / nprocz
					zs = ipz * pnz
					ze = (ipz+1) * pnz + 2*nghostz - 1
				end

				if (not keyword_set (quiet)) then print, 'Writing: '+filename

				; write the actual data
				openw, lun, filename, f77=f77, swap_endian=swap_endian, /get_lun
				writeu, lun, vars[xs:xe,ys:ye,zs:ze,*]
				close, lun

				openw, lun, filename, /f77, swap_endian=swap_endian, /append
				if (allprocs eq 1) then begin
					; collectively written files
					writeu, lun, time, grid.x, grid.y, grid.z, grid.dx, grid.dy, grid.dz
				end else if (allprocs eq 2) then begin
					; xy-collectively written files for each ipz-layer
					writeu, lun, time
					if (iproc eq 0) then writeu, lun, grid.x, grid.y, grid.z, grid.dx, grid.dy, grid.dz
				end else begin
					; distributed files
					writeu, lun, time, grid.x[xs:xe], grid.y[ys:ye], grid.z[zs:ze], grid.dx, grid.dy, grid.dz
				end
				close, lun
				free_lun, lun

			endfor
		endfor
	endfor

end

