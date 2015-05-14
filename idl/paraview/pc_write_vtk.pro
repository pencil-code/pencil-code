;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_write_vtk.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Write VTK file from a data directory or data array.
;
;  Parameters:
;   * data           data array or data structure as load by pc_read_*.
;   * index          indices or tags of the given variables inside data (needed only if 'data' is an array).
;   * filename       filename of the VTK file (without '.vtk' suffix, optional).
;   * datadir        data directory (optional).
;   * grid           grid structure (optional).
;   * dim            dimensions structure (optional).
;   * selected       selected tag names of quantities for export (string array, optional).
;   * quiet          be quiet (optional).
;
;  History:
;   PABourdin, 13.05.2015: coded.
;
;  Examples: (in ascending order of efficiency)
;  ============================================
;
;  * Using 'pc_read_var': (NOT RECOMMENDED)
;
;   Read data from a VAR file and export all content to VTK:
;   IDL> pc_read_var, obj=vars, varfile='var.dat', grid=grid, dim=dim
;   IDL> pc_write_vtk, vars, grid=grid, dim=dim
;
;   Read data from a VAR file and export only selected quantities to VTK:
;   IDL> pc_read_var, obj=vars, varfile='var.dat', grid=grid, dim=dim
;   IDL> pc_write_vtk, vars, selected=['uu','ax','lnrho'], grid=grid, dim=dim
;
;  * Using 'pc_read_var_raw': (RECOMMENDED)
;
;   Read data from a VAR file efficiently and export selected quantities to VTK:
;   IDL> pc_read_var_raw, obj=var, tags=tags, varfile='var.dat', grid=grid, dim=dim
;   IDL> pc_write_vtk, var, tags, selected=['uu','ax','lnrho'], grid=grid, dim=dim
;
;  * Using 'pc_read_subvol_raw': (SUBVOLUME)
;
;   Read a subvolume from a VAR file efficiently and export selected quantities to VTK:
;   IDL> pc_read_subvol_raw, obj=var, tags=tags, varfile='var.dat', xs=3, xe=18, ys=3, ye=34, zs=11, ze=18, sub_grid=grid, sub_dim=dim
;   IDL> pc_write_vtk, var, tags, selected=['uu','ax','lnrho'], grid=grid, dim=dim

pro pc_write_vtk, data, index, filename=filename, datadir=datadir, grid=grid, dim=dim, selected=selected

	common pc_precision, zero, one

	; Default values.
	if (not keyword_set (filename)) then filename = 'export'
	filename += '.vtk'

	; Consistency check.
	if (size (data, /type) eq 0) then message, "ERROR: some 'data' must be given."
	if ((size (data, /type) ne 8) and (size (index, /type) ne 8)) then message, "ERROR: if 'data' is a data array, 'index' must be a structure."

	; Default data directory.
	if (not keyword_set (datadir)) then datadir = pc_get_datadir()

	; Get necessary dimensions quietly.
	if (size (dim, /type) ne 8) then pc_read_dim, object=dim, datadir=datadir, /quiet
	if (size (grid, /type) ne 8) then pc_read_grid, object=grid, dim=dim, datadir=datadir, /quiet
	if (dim.precision eq 'D') then data_type = 'double' else data_type = 'float'

	; Local shorthand for some parameters.
	precision = dim.precision
	nxgrid = dim.nxgrid
	nygrid = dim.nygrid
	nzgrid = dim.nzgrid
	nwgrid = nxgrid * nygrid * nzgrid
	data_type = strlowcase (size (zero, /tname))
	dimensions = [ nxgrid, nygrid, nzgrid ]
	num_dim = total (dimensions gt 1, /int)

	; Open the VTK file for write.
	if (not keyword_set (quiet)) then print, 'Writing "'+strtrim (filename)+'"...'
	openw, lun, file, /get_lun, /swap_if_little_endian

	; Write the header information.
	printf, lun, '# vtk DataFile Version 2.0'
	printf, lun, 'Pencil Code snapshot'
	printf, lun, 'BINARY'
	if (all (lequidist)) then begin
		printf, lun, 'DATASET STRUCTURED_POINTS'
		printf, lun, 'DIMENSIONS ', dimensions
		printf, lun, 'ORIGIN ', [ grid.Ox, grid.Oy, grid.Oz ]
		printf, lun, 'SPACING ', [ grid.dx, grid.dy, grid.dz ]
	endif else begin
		printf, lun, 'DATASET RECTILINEAR_GRID'
		printf, lun, 'DIMENSIONS ', dimensions
		printf, lun, 'X_COORDINATES ', nxgrid, ' ', data_type
		writeu, lun, x
		printf, lun, 'Y_COORDINATES ', nygrid, ' ', data_type
		writeu, lun, y
		printf, lun, 'Z_COORDINATES ', nzgrid, ' ', data_type
		writeu, lun, z
	endelse
	printf, lun, 'POINT_DATA ', nwgrid

	; Find quantities to be written.
	if (size (data, /type) eq 8) then tags = tag_names (data) else tags = tag_names (index)
	num_tags = n_elements (tags)

	; Write selected data quantities.
	skip = [ 'X', 'Y', 'Z', 'DX', 'DY', 'DZ', 'T', 'TIME', 'DELTAY' ]
	for i = 0, num_tags - 1 do begin
		if (any (tags[i] eq skip)) then continue
		if (keyword_set (selected)) then if (not any (tags[i] eq selected)) then continue
		if (size (data, /type) eq 8) then begin
			num_dims = size (data.(i), /n_dimensions)
		endif else begin
			num_dims = size (index.(i), /n_dimensions)
		endelse
		if (num_dims eq 1) then begin
			; Write scalar field.
			print, 'SCALARS ', strlowcase (tags[i]), ' ', data_type, '...'
			printf, lun, 'SCALARS ', strlowcase (tags[i]), ' ', data_type
			printf, lun, 'LOOKUP_TABLE default'
			if (size (data, /type) eq 8) then begin
				writeu, lun, swap_endian (data.(i), /swap_if_big_endian)
			endif else begin
				writeu, lun, swap_endian (data[*,*,*,i], /swap_if_big_endian)
			endelse
		endif else if (num_dims eq 3) then begin
			; Write vector field.
			print, 'VECTORS ', strlowcase (tags[i]), ' ', data_type
			printf, lun, 'VECTORS ', strlowcase (tags[i]), ' ', data_type
			if (size (data, /type) eq 8) then begin
				writeu, lun, swap_endian (transpose (reform (data.(i), [nwgrid, num_dims])), /swap_if_big_endian)
			endif else begin
				writeu, lun, swap_endian (transpose (reform (data[*,*,*,i:i+num_dims-1], [nwgrid, num_dims])), /swap_if_big_endian)
			endelse
			if (not keyword_set (selected)) then i += 3
		endif else begin
			message, "ERROR: data with unrecognized dimension."
		endelse
	endfor

	; Clean up.
	if (not keyword_set (quiet)) then print, 'Ready.'
	close, lun
	free_lun, lun

end

