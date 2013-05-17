;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_seed_points.pro       ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   GUI for selecting a seed points within given coordinates.
;;;   The returned array contains a list of seed points.
;;;   Optional parameters are:
;;;   * start (the starting coordinate, default: whole range)
;;;   * description (description text for the list of seed points)
;;;   * precision (the sub-grid precision for the tracing, default: 0.2)
;;;   * select (select only every select-th point for extraction, default: 5)
;;;
;;;   Example:
;;;   IDL> seeds = pc_seed_points (grid)
;;;   IDL> seeds = pc_seed_points (grid, start=[ px, py, pz ])
;;;


; Update seed points dialog window
pro pc_seed_points_update

	common seed_points_gui_common, sub_xs, sub_xe, sub_nx, sub_ys, sub_ye, sub_ny, sub_zs, sub_ze, sub_nz, sel_dx, sel_dy, sel_dz, num, descr, preci, selec
	common seed_points_common, coord, center, xs, xe, ys, ye, zs, ze, nx, ny, nz, num_x, num_y, num_z, dist_x, dist_y, dist_z, description_str, precision_float, select_int

	if ((dist_x eq 0) and (nx gt 1L)) then begin
		divisors = lindgen (xe - xs + 1L) + 1L
		divisors = divisors[where (((xe - xs + 1L) mod divisors) eq 0)]
		nx = divisors[pc_find_index (double (nx < max (divisors)), divisors, /round)]
	end
	if ((dist_y eq 0) and (ny gt 1L)) then begin
		divisors = lindgen (ye - ys + 1L) + 1L
		divisors = divisors[where (((ye - ys + 1L) mod divisors) eq 0)]
		ny = divisors[pc_find_index (double (ny < max (divisors)), divisors, /round)]
	end
	if ((dist_z eq 0) and (nz gt 1L)) then begin
		divisors = lindgen (ze - zs + 1L) + 1L
		divisors = divisors[where (((ze - zs + 1L) mod divisors) eq 0)]
		nz = divisors[pc_find_index (double (nz < max (divisors)), divisors, /round)]
	end

	WIDGET_CONTROL, sub_xs, SET_VALUE = xs + coord.x_off
	WIDGET_CONTROL, sub_xe, SET_VALUE = xe + coord.x_off
	WIDGET_CONTROL, sub_ys, SET_VALUE = ys + coord.y_off
	WIDGET_CONTROL, sub_ye, SET_VALUE = ye + coord.y_off
	WIDGET_CONTROL, sub_zs, SET_VALUE = zs + coord.z_off
	WIDGET_CONTROL, sub_ze, SET_VALUE = ze + coord.z_off
	WIDGET_CONTROL, sub_nx, SET_VALUE = nx
	WIDGET_CONTROL, sub_ny, SET_VALUE = ny
	WIDGET_CONTROL, sub_nz, SET_VALUE = nz

	WIDGET_CONTROL, num, SET_VALUE = nx*ny*nz
end


; Event handling of seed points dialog window
pro seed_points_event, event

	common seed_points_gui_common, sub_xs, sub_xe, sub_nx, sub_ys, sub_ye, sub_ny, sub_zs, sub_ze, sub_nz, sel_dx, sel_dy, sel_dz, num, descr, preci, selec
	common seed_points_common, coord, center, xs, xe, ys, ye, zs, ze, nx, ny, nz, num_x, num_y, num_z, dist_x, dist_y, dist_z, description_str, precision_float, select_int

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1

	SWITCH eventval of
	'RESET': begin
		xs = center[0]
		ys = center[0]
		zs = center[0]
		xe = xs
		ye = ys
		ze = zs
		nx = 1L
		ny = 1L
		nz = 1L
		break
	end
	'SUB_XS': begin
		WIDGET_CONTROL, event.id, GET_VALUE = xs
		xs = ((xs - coord.x_off) > 0L) < (num_x-1L)
		xe = xe > xs
		break
	end
	'SUB_YS': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ys
		ys = ((ys - coord.y_off) > 0L) < (num_y-1L)
		ye = ye > ys
		break
	end
	'SUB_ZS': begin
		WIDGET_CONTROL, event.id, GET_VALUE = zs
		zs = ((zs - coord.z_off) > 0L) < (num_z-1L)
		ze = ze > zs
		break
	end
	'SUB_XE': begin
		WIDGET_CONTROL, event.id, GET_VALUE = xe
		xe = ((xe - coord.x_off) > 0L) < (num_x-1L)
		xs = xs < xe
		break
	end
	'SUB_YE': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ye
		ye = ((ye - coord.y_off) > 0L) < (num_y-1L)
		ys = ys < ye
		break
	end
	'SUB_ZE': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ze
		ze = ((ze - coord.z_off) > 0L) < (num_z-1L)
		zs = zs < ze
		break
	end
	'SUB_NX': begin
		WIDGET_CONTROL, event.id, GET_VALUE = nx
		nx = nx > 1L
		break
	end
	'SUB_NY': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ny
		ny = ny > 1L
		break
	end
	'SUB_NZ': begin
		WIDGET_CONTROL, event.id, GET_VALUE = nz
		nz = nz > 1L
		break
	end
	'ALL_X': begin
		xs = 0L
		xe = num_x - 1L
		nx = xe - xs + 1L
		break
	end
	'ALL_Y': begin
		ys = 0L
		ye = num_y - 1L
		ny = ye - ys + 1L
		break
	end
	'ALL_Z': begin
		zs = 0L
		ze = num_z - 1L
		nz = ze - zs + 1L
		break
	end
	'SEL_DX': begin
		dist_x = event.index
		break
	end
	'SEL_DY': begin
		dist_y = event.index
		break
	end
	'SEL_DZ': begin
		dist_z = event.index
		break
	end
	'OK': begin
		WIDGET_CONTROL, descr, GET_VALUE = description_str
		WIDGET_CONTROL, preci, GET_VALUE = precision_float
		WIDGET_CONTROL, selec, GET_VALUE = select_int
		quit = event.top
		break
	end
	'CANCEL': begin
		nx = 0L
		ny = 0L
		nz = 0L
		quit = event.top
		break
	end
	endswitch

	pc_seed_points_update
	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; File selection dialog GUI.
function pc_seed_points, grid, start=start, description=description, precision=precision, select=select

	common seed_points_gui_common, sub_xs, sub_xe, sub_nx, sub_ys, sub_ye, sub_ny, sub_zs, sub_ze, sub_nz, sel_dx, sel_dy, sel_dz, num, descr, preci, selec
	common seed_points_common, coord, center, xs, xe, ys, ye, zs, ze, nx, ny, nz, num_x, num_y, num_z, dist_x, dist_y, dist_z, description_str, precision_float, select_int

	; Default settings:
	default, description, "pc_seed_points"
	default, precision, 0.2
	default, select, 5

	description_str = description
	precision_float = precision
	select_int = select

	coord = grid
	if (not any (strcmp (tag_names (coord), 'x_off', /fold_case))) then coord = create_struct (coord, 'x_off', 0)
	if (not any (strcmp (tag_names (coord), 'y_off', /fold_case))) then coord = create_struct (coord, 'y_off', 0)
	if (not any (strcmp (tag_names (coord), 'z_off', /fold_case))) then coord = create_struct (coord, 'z_off', 0)

	num_x = n_elements (coord.x)
	num_y = n_elements (coord.y)
	num_z = n_elements (coord.z)

	xs = 0L
	xe = num_x - 1L
	ys = 0L
	ye = num_y - 1L
	zs = 0L
	ze = num_z - 1L

	if (keyword_set (start)) then begin
		xs = (round (start[0]) > xs) < xe
		ys = (round (start[1]) > ys) < ye
		zs = (round (start[2]) > zs) < ze
		xe = xs
		ye = ys
		ze = zs
	end
	center = [ xs, ys, zs ]

	nx = xe - xs + 1L
	ny = ye - ys + 1L
	nz = ze - zs + 1L


	; Build GUI
	MOTHER	= WIDGET_BASE (title='PC seed points selector', EVENT_PRO=seed_points_event)
	BASE	= WIDGET_BASE (MOTHER, /row)

	CTRL	= WIDGET_BASE (BASE, /col)

	dist_names = [ 'exact grid points', 'equidistant distribution', 'random distribution' ]
	dist_x = 0
	dist_y = 0
	dist_z = 0

	tmp	= WIDGET_LABEL (CTRL, value='Seed Point Selection (start, end, number):', frame=0)
	SEL	= WIDGET_BASE (CTRL, frame=1, /align_center, /col)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_xs	= CW_FIELD (SUB, title='X:', uvalue='SUB_XS', value=(xs + coord.x_off), /long, /return_events, xsize=5)
	tmp	= WIDGET_BUTTON (SUB, xsize=100, value='<= MIN MAX =>', uvalue='ALL_X')
	sub_xe	= CW_FIELD (SUB, title='', uvalue='SUB_XE', value=(xe + coord.x_off), /long, /return_events, xsize=5)
	sel_dx	= WIDGET_DROPLIST (SUB, value=dist_names, uvalue='SEL_DX')
	sub_nx	= CW_FIELD (SUB, title='=', uvalue='SUB_NX', value=nx, /long, /return_events, xsize=5)
	tmp	= WIDGET_LABEL (SUB, value='seed points', frame=0)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_ys	= CW_FIELD (SUB, title='Y:', uvalue='SUB_YS', value=(ys + coord.y_off), /long, /return_events, xsize=5)
	tmp	= WIDGET_BUTTON (SUB, xsize=100, value='<= MIN MAX =>', uvalue='ALL_Y')
	sub_ye	= CW_FIELD (SUB, title='', uvalue='SUB_YE', value=(ye + coord.y_off), /long, /return_events, xsize=5)
	sel_dy	= WIDGET_DROPLIST (SUB, value=dist_names, uvalue='SEL_DY')
	sub_ny	= CW_FIELD (SUB, title='=', uvalue='SUB_NY', value=ny, /long, /return_events, xsize=5)
	tmp	= WIDGET_LABEL (SUB, value='seed points', frame=0)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_zs	= CW_FIELD (SUB, title='Z:', uvalue='SUB_ZS', value=(zs + coord.z_off), /long, /return_events, xsize=5)
	tmp	= WIDGET_BUTTON (SUB, xsize=100, value='<= MIN MAX =>', uvalue='ALL_Z')
	sub_ze	= CW_FIELD (SUB, title='', uvalue='SUB_ZE', value=(ze + coord.z_off), /long, /return_events, xsize=5)
	sel_dz	= WIDGET_DROPLIST (SUB, value=dist_names, uvalue='SEL_DZ')
	sub_nz	= CW_FIELD (SUB, title='=', uvalue='SUB_NZ', value=nz, /long, /return_events, xsize=5)
	tmp	= WIDGET_LABEL (SUB, value='seed points', frame=0)
	WIDGET_CONTROL, sub_xs, SENSITIVE = (num_x gt 1L)
	WIDGET_CONTROL, sub_xe, SENSITIVE = (num_x gt 1L)
	WIDGET_CONTROL, sub_nx, SENSITIVE = (num_x gt 1L)
	WIDGET_CONTROL, sub_ys, SENSITIVE = (num_y gt 1L)
	WIDGET_CONTROL, sub_ye, SENSITIVE = (num_y gt 1L)
	WIDGET_CONTROL, sub_ny, SENSITIVE = (num_y gt 1L)
	WIDGET_CONTROL, sub_zs, SENSITIVE = (num_z gt 1L)
	WIDGET_CONTROL, sub_ze, SENSITIVE = (num_z gt 1L)
	WIDGET_CONTROL, sub_nz, SENSITIVE = (num_z gt 1L)

	tmp	= WIDGET_BASE (CTRL, frame=0, /align_center, /col)
	BUT	= WIDGET_BASE (tmp, frame=0, /align_center, /row)
	descr	= CW_FIELD (BUT, title='Description:', xsize=40)
	BUT	= WIDGET_BASE (tmp, frame=0, /align_center, /row)
	preci	= CW_FIELD (BUT, title='Precision:', /floating, xsize=10)
	selec	= CW_FIELD (BUT, title='Select every:', /long, xsize=10)
	WIDGET_CONTROL, descr, SET_VALUE = description_str
	WIDGET_CONTROL, preci, SET_VALUE = precision_float
	WIDGET_CONTROL, selec, SET_VALUE = select_int

	BUT	= WIDGET_BASE (CTRL, frame=0, /align_center, /row)
	num	= CW_FIELD (BUT, title='Number of Streamlines:', /long, /noedit, xsize=10)
	WIDGET_CONTROL, num, SET_VALUE = nx*ny*nz
	SUB	= WIDGET_BASE (BUT, frame=0, /align_center, /row)
	tmp	= WIDGET_BUTTON (SUB, xsize=80, value='RESET', uvalue='RESET')
	GRP	= WIDGET_BASE (BUT, frame=1, /align_center, /row)
	tmp	= WIDGET_BUTTON (GRP, xsize=80, value='CANCEL', uvalue='CANCEL')
	tmp	= WIDGET_BUTTON (GRP, xsize=80, value='OK', uvalue='OK')

	WIDGET_CONTROL, MOTHER, /REALIZE
	wimg = !d.window

	WIDGET_CONTROL, BASE

	XMANAGER, "seed_points", MOTHER


	; Check for abortion or impossible values
	if (nx * ny * nz le 0L) then return, -1L

	; Build list of selected seed points
	if (dist_x eq 2) then rx = 1 else rx = 0
	if (dist_y eq 2) then ry = 1 else ry = 0
	if (dist_z eq 2) then rz = 1 else rz = 0
	x_start = coord.x[xs] - rx * 0.5 * coord.dx[xs]
	y_start = coord.y[ys] - ry * 0.5 * coord.dy[ys]
	z_start = coord.z[zs] - rz * 0.5 * coord.dz[zs]
	x_end = coord.x[xe] + rx * 0.5 * coord.dx[xe]
	y_end = coord.y[ye] + ry * 0.5 * coord.dy[ye]
	z_end = coord.z[ze] + rz * 0.5 * coord.dz[ze]
	if (grid.lperi[0]) then x_start = x_start > coord.x[0] & x_end = x_end < coord.x[num_x-1]
	if (grid.lperi[1]) then y_start = y_start > coord.y[0] & y_end = y_end < coord.y[num_y-1]
	if (grid.lperi[2]) then z_start = z_start > coord.z[0] & z_end = z_end < coord.z[num_z-1]
	dx = x_end - x_start
	dy = y_end - y_start
	dz = z_end - z_start
	seeds = dblarr (3, nx * ny * nz)
	seed = long (systime (/seconds))
	pos = 0L
	if ((dist_x eq 0) and (num_x ne nx)) then step_x = (xe - xs + 1L) / nx else step_x = 1L
	if ((dist_y eq 0) and (num_y ne ny)) then step_y = (ye - ys + 1L) / ny else step_y = 1L
	if ((dist_z eq 0) and (num_z ne nz)) then step_z = (ze - zs + 1L) / nz else step_z = 1L
	for pos_z = 0L, nz - 1L do begin
		for pos_y = 0L, ny - 1L do begin
			for pos_x = 0L, nx - 1L do begin
				if (rx ne 0) then px = randomu (seed, /double) else if (nx le 1L) then px = pos_x / (nx - 1.0d0) else px = 0.5d0
				if (ry ne 0) then py = randomu (seed, /double) else if (ny le 1L) then py = pos_y / (ny - 1.0d0) else py = 0.5d0
				if (rz ne 0) then pz = randomu (seed, /double) else if (nz le 1L) then pz = pos_z / (nz - 1.0d0) else pz = 0.5d0
				if (dist_x eq 0) then x = coord.x[xs + pos_x * step_x] else x = x_start + px * dx
				if (dist_y eq 0) then y = coord.y[ys + pos_y * step_y] else y = y_start + py * dy
				if (dist_z eq 0) then z = coord.z[zs + pos_z * step_z] else z = z_start + pz * dz
				seeds[*,pos] = [ x, y, z ]
				pos += 1L
			end
		end
	end

	description = description_str
	if (size (description, /n_dimensions) eq 1) then description = description[0]
	precision = precision_float
	select = select_int

	return, reform (seeds)
end

