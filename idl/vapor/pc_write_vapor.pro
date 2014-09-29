; Program to convert a directory of pencil files into a vapor dataset.
; Arguments are:
;	datadir  = absolute path to directory where pencil data is stored.
;	nprocs   = number of procN subdirectories (N goes from 0 to nprocs -1).
;	vdf_file = absolute path, including filename, of VDF-file for the result.
;	           The VDF-file should specify names for the variables, in the
;	           same order as they appear in the Pencil "var.dat" files.
;	varfile  = a string (or array of strings) containing the varfiles to be read.
;	           (Default: "var.dat", if nothing is set.)
;       ivar     = varfile numbers, or a min-max range of varfile numbers to be read.
;	 
; 	The chunks associated with the different processors are not assumed to be
;	the same size, however, all chunks at a given x-coordinate have the same x-thickness,
;	and likewize for the chunks at a given y or z coordinate
;
; $Id$

pro pc_write_vapor, vdf_file=vdf_file, varfile=varfile, datadir=datadir, $
	proc=proc, ivar=ivar, variables=variables, proc=proc, scalar=scalar, $
	varcontent=varcontent, additional=additional, _extra=_extra

	pc_read_dim, obj=dim, proc=proc
	if (size (proc, /type) ne 0) then nprocs = 1 else nprocs = dim.nprocx*dim.nprocy*dim.nprocz

	if (n_elements (ivar) eq 2) then begin
		range = minmax (ivar)
		varfile = "VAR"+strtrim (indgen (range[1]-range[0]+1)+range[0], 2)
	end else if (n_elements (ivar) ge 1) then begin
		varfile = "VAR"+strtrim (ivar, 2)
	end
	num_files = n_elements (varfile)

	default, datadir, 'data'
	default, vdf_file, datadir+'/var.vdf'

	default, varcontent, pc_varcontent (datadir=datadir, dim=dim, param=param, quiet=quiet, scalar=scalar)
	filevars = varcontent[where ((varcontent[*].idlvar ne 'dummy'))].idlvar
	if (size (variables, /type) eq 0) then variables = filevars
	if (keyword_set (additional)) then variables = [ variables, additional ]

	; Usage: vdfcreate [options] filename
	;    -dimensio arg0 : Volume dimensions (NXxNYxNZ)
	;    -numts    arg0 : Number of timesteps
	;    -bs       arg0 : Internal storage blocking factor (NXxNYxNZ)
	;    -level    arg0 : Maximum refinement level. 0 => no refinement
	;    -nfilter  arg0 : Number of wavelet filter coefficients
	;    -nlifting arg0 : Number of wavelet lifting coefficients
	;    -comment  arg0 : Top-level comment
	;    -gridtype arg0 : Data grid type (regular|streched|block_amr)
	;    -coordsys arg0 : Top-level comment (cartesian|spherical)
	;    -extents  arg0 : Domain extents in user coordinates
	;    -varnames arg0 : Colon delimited list of variable names
	;    -mtkcompa      : Force compatibility with older mtk files
	;    -help          : Print this message and exit
	dim_str = strtrim (dim.nx, 2) + 'x' + strtrim (dim.ny, 2) + 'x' + strtrim (dim.nz, 2)
	vdfcreate_command = 'vdfcreate -dimension ' + dim_str $
		+ ' -numts ' + strtrim (num_files, 2) $
		+ ' -comment "Created by pc_write_vapor"' $
		+ ' -gridtype regular' $
		+ ' -coordsys cartesian' $
		+ ' -varnames ' + arraytostring (variables, list=':', /noleader) $
		+ ' ' + vdf_file 

	print, "Running: " + vdfcreate_command
	spawn, vdfcreate_command, exit_status=exit_status
	if (exit_status ne 0) then message, "ERROR: Failed to create VDF-file."

	; Create the metadata for the VDF-file.
	print,'vdf_create "' + vdf_file + '"'
	mfd = vdf_create (vdf_file)

	; Iterate over varfiles
	undefine, varcontent
	for file_pos = 0, num_files - 1 do begin

		; Read one varfile
		pc_read_var, obj=data, proc=proc, varfile=varfile[file_pos], $
			variables=variables, /trimall, /quiet, scalar=scalar, $
			varcontent=varcontent, additional=additional, _extra=_extra

		data_size = size (data, /dimensions)
		if (n_elements (data_size) ge 3) then max_z = data_size[3] else max_z = 0

		; Iterate over varcontent
		num_variables = n_elements (varcontent)
		for pos = 0, num_variables - 1 do begin
			varname = varcontent[pos].idlvar
			if (varcontent[pos].skip eq 2) then varname = varname + [ '_x', '_y', '_z' ]
			

			print, 'file: ', varfile[file_pos], ' variable: ', varnames[pos]
			dfd = vdc_bufwritecreate (mfd)
			vdc_openvarwrite, dfd, file_pos, varnames[pos], -1

			; Write the full slab to the VDF, iterate over z-slices
			min_a = 0
			max_a = 0
			for z = 0, max_z do begin
				for a = min_a, max_a do begin
					vdc_bufwriteslice, dfd, data[*,*,z,pos+a]
				end
			end
			pos = pos + 2

			; Now close the file
			vdc_closevar, dfd
			vdc_bufwritedestroy, dfd
		end
	end
end

