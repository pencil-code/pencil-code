; $Id$
;
; Description:
;   Saves a given data array to a VDF2 file for Vapor.
;
; Parameters:
;   * vdf_file       Output VDF2 file name.
;   * quantity       Quantity name to write.
;   * data           Data array as output from 'pc_get_quantity'.
;   * var            Variable array as given to 'pc_get_quantity'.
;   * tags           'tags' structure as given to 'pc_get_quantity'.
;
; Optional parameters:
;   * coarsening     Number of coarsening levels (Default: 0 = off).
;   * reduce         Factor for reduction of the data (Default: 1 = off).
;
; Examples:
; =========
;
;   Load parts of a varfile and save the magnetic flux density to a VDF2 file:
;   IDL> pc_read_var_raw, obj=var, tags=tags, varfile='var.dat', var_list=['aa'], dim=dim, grid=grid
;   IDL> B_abs = pc_get_quantity ('B_abs', var, tags, dim=dim, grid=grid)
;   IDL> pc_write_vdf, 'B_abs.vdf', 'B_abs', B_abs, var, tags, dim=dim, grid=grid

pro pc_save_vdf, vdf_file, quantity, data, var, tags, timestep=timestep, max_timesteps=max_timesteps, reset=reset, coarsening=coarsening, reduce=reduce, dim=dim, grid=grid

	; default settings
	default, timestep, 0
	default, max_timesteps, max (timestep) + 1
	default, reset, 0
	default, coarsening, 0
	default, reduce, 1
	if (n_elements (reduce) eq 1) then reduce = replicate (reduce, 3)

	; consistency checks
	if (size (grid, /type) eq 0) then message, 'pc_write_vdf: need a grid structure'
	if (size (dim, /type) eq 0) then message, 'pc_write_vdf: need a dim structure'
	if (size (data, /type) eq 8) then message, 'pc_write_vdf: need a data array'

	if (reset or not file_test (vdf_file)) then begin
		; create a new VDF metadata object
		vdf_dim = [ dim.nx, dim.ny, dim.nz ]
		mdo = vdf_create (vdf_dim, coarsening)

		; set the maximum number of timesteps
		vdf_setnumtimesteps, mdo, max_timesteps

		; set periodicity
		vdf_setperiodic, mdo, start_param.lperi

		; set box size
		size_x = pc_get_quantity ('size_x', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
		size_y = pc_get_quantity ('size_y', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
		size_z = pc_get_quantity ('size_z', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
		origin_x = pc_get_quantity ('origin_x', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
		origin_y = pc_get_quantity ('origin_y', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
		origin_z = pc_get_quantity ('origin_z', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
		vdf_setextents, mdo, [ origin_x, origin_y, origin_z, origin_x+size_x, origin_y+size_y, origin_z+size_z ]

		; set grid type
		if (any (start_param.lequidist ne 1)) then grid_type = 'stretched' else grid_type = 'regular'
		vdf_setgridtype, mdo, grid_type

		; set the names of the variables
		vdf_setvarnames, mdo, quantities

		; store and close the metadata object
		vdf_write, mdo, vdf_file
		vdf_destroy, mdo

		; no reset when writing additional timesteps later
		reset = 0
	end

	; open existing VDF metadata object
	mdo = vdf_create (vdf_file)

	; set time of snapshot
	time = pc_get_quantity ('time', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
	vdf_settusertime, mdo, timestep, [ time ] ; possible bug in VAPOR-IDL: time is expected to be a 1-element array

	; get grid coordinates
	x = pc_get_quantity ('x', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
	y = pc_get_quantity ('y', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache)
	z = pc_get_quantity ('z', var, tags, units=units, dim=dim, grid=grid, start_param=start_param, run_param=run_param, /cache, /cleanup)

	; reduce grid resolution and the data, if necessary
	if (any (reduce ne 1)) then begin
		num_layers = round (dim.nz/reduce[2])
		data = congrid (data, round (dim.nx/reduce[0]), round (dim.ny/reduce[1]), round (num_layers), /cubic, /interpolate)
		x = congrid (x, round (dim.nx/reduce[0]), 1, 1, /cubic, /interpolate)
		y = congrid (y, round (dim.ny/reduce[1]), 1, 1, /cubic, /interpolate)
		z = congrid (z, round (dim.nz/reduce[2]), 1, 1, /cubic, /interpolate)
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
	data = quantity
	num_layers = dim.nz

	; write the xy-slices
	vdc_openvarwrite, bwh, timestep, quantity, -1
	for z = 0, num_layers-1 do vdc_bufwriteslice, bwh, data[*,*,z]
	vdc_closevar, bwh

	; destroy the "buffered write" handle
	vdc_bufwritedestroy, bwh
end

