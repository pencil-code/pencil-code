; $Id$
;
; Description:
;   Saves a given data array to a VDF2 file for Vapor.
;
; Parameters:
;   * vdf_file       Output VDF2 file name.
;   * quantity       Quantity name to write.
;   * data           Data array as output from 'pc_get_quantity'.
;   * x              x-coordinates.
;   * y              y-coordinates.
;   * z              z-coordinates.
;   * size           Array with the box dimensions.
;   * origin         Array with the box origin coordinates.
;   * periodic       Array indicating the periodicity.
;   * equidist       Array indicating the equidistancy.
;   * time           Time of the data.
;
; Optional parameters:
;   * timestep       Number of the timestep (Default: 0).
;   * max_timesteps  Number of timesteps (Default: 1).
;   * reset          Force to overwrite an existing file (Default: 0 = off).
;   * coarsening     Number of coarsening levels (Default: 0 = off).
;   * reduce         Factor for reduction of the data (Default: 1 = off).
;
; Examples:
; =========
;
;   Load parts of a varfile and save the magnetic flux density to a VDF2 file:
;   IDL> pc_read_var_raw, obj=var, tags=tags, varfile='var.dat', var_list=['aa'], dim=dim, grid=grid
;   IDL> B_abs = pc_get_quantity ('B_abs', var, tags, dim=dim, grid=grid)
;   IDL> pc_save_vdf, 'B_abs.vdf', 'B_abs', B_abs, x, y, z, size, origin, periodic, equidist, time

pro pc_save_vdf, vdf_file, quantity, data, x, y, z, size, origin, periodic, equidist, time, timestep=timestep, max_timesteps=max_timesteps, reset=reset, coarsening=coarsening, reduce=reduce

	; default settings
	default, timestep, 0
	default, max_timesteps, max (timestep) + 1
	default, reset, 0
	default, coarsening, 0
	default, reduce, 1
	if (n_elements (reduce) eq 1) then reduce = replicate (reduce, 3)
	nx = n_elements (x)
	ny = n_elements (y)
	nz = n_elements (z)

	; consistency checks
	if (size (data, /type) eq 8) then message, 'pc_write_vdf: need a data array'

	if (reset or not file_test (vdf_file)) then begin
		; create a new VDF metadata object
		vdf_dim = [ nx, ny, nz ]
		mdo = vdf_create (vdf_dim, coarsening)

		; set the maximum number of timesteps
		vdf_setnumtimesteps, mdo, max_timesteps

		; set periodicity
		vdf_setperiodic, mdo, periodic

		; set box size
		vdf_setextents, mdo, [ origin, origin+size ]

		; set grid type
		if (any (equidist ne 1)) then grid_type = 'stretched' else grid_type = 'regular'
		vdf_setgridtype, mdo, grid_type

		; set the names of the variables
		vdf_setvarnames, mdo, [ quantity ]

		; store and close the metadata object
		vdf_write, mdo, vdf_file
		vdf_destroy, mdo

		; no reset when writing additional timesteps later
		reset = 0
	end

	; open existing VDF metadata object
	mdo = vdf_create (vdf_file)

	; set time of snapshot
	vdf_settusertime, mdo, timestep, [ time ] ; possible bug in VAPOR-IDL: time is expected to be a 1-element array

	; reduce grid resolution and the data, if necessary
	if (any (reduce ne 1)) then begin
		num_layers = round (nz/reduce[2])
		data = congrid (data, round (nx/reduce[0]), round (ny/reduce[1]), round (num_layers), /cubic, /interpolate)
		x = congrid (x, round (nx/reduce[0]), 1, 1, /cubic, /interpolate)
		y = congrid (y, round (ny/reduce[1]), 1, 1, /cubic, /interpolate)
		z = congrid (z, round (nz/reduce[2]), 1, 1, /cubic, /interpolate)
	end

	; set grid coordinates
	vdf_settxcoords, mdo, timestep, x
	vdf_settycoords, mdo, timestep, y
	vdf_settzcoords, mdo, timestep, z

	; store and close the changed metadata object
	vdf_write, mdo, vdf_file
	vdf_destroy, mdo

	; create a "buffered write" handle
	bwh = vdc_bufwritecreate (vdf_file)

	; write the quantity
	num_layers = nz

	; write the xy-slices
	vdc_openvarwrite, bwh, timestep, quantity, -1
	for pz = 0, num_layers-1 do begin
		vdc_bufwriteslice, bwh, reform (float (data[*,*,pz]))
	end
	vdc_closevar, bwh

	; destroy the "buffered write" handle
	vdc_bufwritedestroy, bwh
end

