;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_seed_points.pro      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   GUI for selecting a seed points within given coordinates.
;;;   The returned array contains a list of seed points.
;;;   Optional parameters are:
;;;   * start (the starting coordinate, default: whole range)
;;;
;;;   Example:
;;;   IDL> seeds = pc_seed_points (grid)
;;;   IDL> seeds = pc_seed_points (grid, start=[ px, py, pz ])
;;;


; Update seed points dialog window
pro pc_seed_points_update

	common seed_points_gui_common, sub_xs, sub_xe, sub_nx, sub_ys, sub_ye, sub_ny, sub_zs, sub_ze, sub_nz, sel_dx, sel_dy, sel_dz, num
	common seed_points_common, coord, center, xs, xe, ys, ye, zs, ze, nx, ny, nz, num_x, num_y, num_z, dist_x, dist_y, dist_z

	if (dist_x eq 0) then nx = xe - xs + 1L
	if (dist_y eq 0) then ny = ye - ys + 1L
	if (dist_z eq 0) then nz = ze - zs + 1L

	WIDGET_CONTROL, sub_xs, SET_VALUE = xs
	WIDGET_CONTROL, sub_xe, SET_VALUE = xe
	WIDGET_CONTROL, sub_ys, SET_VALUE = ys
	WIDGET_CONTROL, sub_ye, SET_VALUE = ye
	WIDGET_CONTROL, sub_zs, SET_VALUE = zs
	WIDGET_CONTROL, sub_ze, SET_VALUE = ze
	WIDGET_CONTROL, sub_nx, SET_VALUE = nx
	WIDGET_CONTROL, sub_ny, SET_VALUE = ny
	WIDGET_CONTROL, sub_nz, SET_VALUE = nz

	WIDGET_CONTROL, sub_nx, SENSITIVE = (num_x gt 1) and (dist_x ne 0)
	WIDGET_CONTROL, sub_ny, SENSITIVE = (num_y gt 1) and (dist_y ne 0)
	WIDGET_CONTROL, sub_nz, SENSITIVE = (num_z gt 1) and (dist_z ne 0)

	WIDGET_CONTROL, num, SET_VALUE = nx*ny*nz
end


; Event handling of seed points dialog window
pro seed_points_event, event

	common seed_points_gui_common, sub_xs, sub_xe, sub_nx, sub_ys, sub_ye, sub_ny, sub_zs, sub_ze, sub_nz, sel_dx, sel_dy, sel_dz, num
	common seed_points_common, coord, center, xs, xe, ys, ye, zs, ze, nx, ny, nz, num_x, num_y, num_z, dist_x, dist_y, dist_z

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
		pc_seed_points_update
		break
	end
	'SUB_XS': begin
		WIDGET_CONTROL, event.id, GET_VALUE = xs
		xs = (xs > 0) < (num_x-1)
		xe = xe > xs
		pc_seed_points_update
		break
	end
	'SUB_YS': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ys
		ys = (ys > 0) < (num_y-1)
		ye = ye > ys
		pc_seed_points_update
		break
	end
	'SUB_ZS': begin
		nz = ze - zs + 1
		WIDGET_CONTROL, event.id, GET_VALUE = zs
		zs = (zs > 0) < (num_z-1)
		ze = ze > zs
		pc_seed_points_update
		break
	end
	'SUB_XE': begin
		WIDGET_CONTROL, event.id, GET_VALUE = xe
		xe = (xe > 0) < (num_x-1)
		xs = xs < xe
		pc_seed_points_update
		break
	end
	'SUB_YE': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ye
		ye = (ye > 0) < (num_y-1)
		ys = ys < ye
		pc_seed_points_update
		break
	end
	'SUB_ZE': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ze
		ze = (ze > 0) < (num_z-1)
		zs = zs < ze
		pc_seed_points_update
		break
	end
	'SUB_NX': begin
		WIDGET_CONTROL, event.id, GET_VALUE = nx
		nx = nx > 1L
		pc_seed_points_update
		break
	end
	'SUB_NY': begin
		WIDGET_CONTROL, event.id, GET_VALUE = ny
		ny = ny > 1L
		pc_seed_points_update
		break
	end
	'SUB_NZ': begin
		WIDGET_CONTROL, event.id, GET_VALUE = nz
		nz = nz > 1L
		pc_seed_points_update
		break
	end
	'ALL_X': begin
		xs = 0
		xe = num_x - 1
		nx = xe - xs + 1L
		pc_seed_points_update
		break
	end
	'ALL_Y': begin
		ys = 0
		ye = num_y - 1
		ny = ye - ys + 1L
		pc_seed_points_update
		break
	end
	'ALL_Z': begin
		zs = 0
		ze = num_z - 1
		nz = ze - zs + 1L
		pc_seed_points_update
		break
	end
	'SEL_DX': begin
		dist_x = event.index
		pc_seed_points_update
		break
	end
	'SEL_DY': begin
		dist_y = event.index
		pc_seed_points_update
		break
	end
	'SEL_DZ': begin
		dist_z = event.index
		pc_seed_points_update
		break
	end
	'OK': begin
		quit = event.top
		break
	end
	'CANCEL': begin
		nx = 0
		ny = 0
		nz = 0
		quit = event.top
		break
	end
	endswitch

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; File selection dialog GUI.
function pc_seed_points, grid, start=start

	common seed_points_gui_common, sub_xs, sub_xe, sub_nx, sub_ys, sub_ye, sub_ny, sub_zs, sub_ze, sub_nz, sel_dx, sel_dy, sel_dz, num
	common seed_points_common, coord, center, xs, xe, ys, ye, zs, ze, nx, ny, nz, num_x, num_y, num_z, dist_x, dist_y, dist_z

	coord = grid

	num_x = n_elements (coord.x)
	num_y = n_elements (coord.y)
	num_z = n_elements (coord.z)

	xs = 0
	xe = num_x - 1
	ys = 0
	ye = num_y - 1
	zs = 0
	ze = num_z - 1

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
	sub_xs	= CW_FIELD (SUB, title='X:', uvalue='SUB_XS', value=xs, /integer, /return_events, xsize=5)
	tmp	= WIDGET_BUTTON (SUB, xsize=100, value='<= MIN MAX =>', uvalue='ALL_X')
	sub_xe	= CW_FIELD (SUB, title='', uvalue='SUB_XE', value=xe, /integer, /return_events, xsize=5)
	sel_dx	= WIDGET_DROPLIST (SUB, value=dist_names, uvalue='SEL_DX')
	sub_nx	= CW_FIELD (SUB, title='=', uvalue='SUB_NX', value=nx, /integer, /return_events, xsize=5)
	tmp	= WIDGET_LABEL (SUB, value='seed points', frame=0)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_ys	= CW_FIELD (SUB, title='Y:', uvalue='SUB_YS', value=ys, /integer, /return_events, xsize=5)
	tmp	= WIDGET_BUTTON (SUB, xsize=100, value='<= MIN MAX =>', uvalue='ALL_Y')
	sub_ye	= CW_FIELD (SUB, title='', uvalue='SUB_YE', value=ye, /integer, /return_events, xsize=5)
	sel_dy	= WIDGET_DROPLIST (SUB, value=dist_names, uvalue='SEL_DY')
	sub_ny	= CW_FIELD (SUB, title='=', uvalue='SUB_NY', value=ny, /integer, /return_events, xsize=5)
	tmp	= WIDGET_LABEL (SUB, value='seed points', frame=0)
	SUB	= WIDGET_BASE (SEL, frame=0, /align_center, /row)
	sub_zs	= CW_FIELD (SUB, title='Z:', uvalue='SUB_ZS', value=zs, /integer, /return_events, xsize=5)
	tmp	= WIDGET_BUTTON (SUB, xsize=100, value='<= MIN MAX =>', uvalue='ALL_Z')
	sub_ze	= CW_FIELD (SUB, title='', uvalue='SUB_ZE', value=ze, /integer, /return_events, xsize=5)
	sel_dz	= WIDGET_DROPLIST (SUB, value=dist_names, uvalue='SEL_DZ')
	sub_nz	= CW_FIELD (SUB, title='=', uvalue='SUB_NZ', value=nz, /integer, /return_events, xsize=5)
	tmp	= WIDGET_LABEL (SUB, value='seed points', frame=0)
	WIDGET_CONTROL, sub_xs, SENSITIVE = (num_x gt 1)
	WIDGET_CONTROL, sub_xe, SENSITIVE = (num_x gt 1)
	WIDGET_CONTROL, sub_nx, SENSITIVE = (num_x gt 1) and (dist_x ne 0)
	WIDGET_CONTROL, sub_ys, SENSITIVE = (num_y gt 1)
	WIDGET_CONTROL, sub_ye, SENSITIVE = (num_y gt 1)
	WIDGET_CONTROL, sub_ny, SENSITIVE = (num_y gt 1) and (dist_y ne 0)
	WIDGET_CONTROL, sub_zs, SENSITIVE = (num_z gt 1)
	WIDGET_CONTROL, sub_ze, SENSITIVE = (num_z gt 1)
	WIDGET_CONTROL, sub_nz, SENSITIVE = (num_z gt 1) and (dist_z ne 0)

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
	dx = x_end - x_start
	dy = y_end - y_start
	dz = z_end - z_start
	seeds = dblarr (nx * ny * nz, 3)
	seed = long (systime (/seconds))
	pos = 0L
	for pos_z = 0, nz - 1 do begin
		for pos_y = 0, ny - 1 do begin
			for pos_x = 0, nx - 1 do begin
				if (rx ne 0) then px = randomu (seed, /double) else if (nx gt 1) then px = pos_x / (nx - 1.0d0) else px = 0.5d0
				if (ry ne 0) then py = randomu (seed, /double) else if (ny gt 1) then py = pos_y / (ny - 1.0d0) else py = 0.5d0
				if (rz ne 0) then pz = randomu (seed, /double) else if (nz gt 1) then pz = pos_z / (nz - 1.0d0) else pz = 0.5d0
				if (dist_x eq 0) then x = coord.x[xs + pos_x] else x = x_start + px * dx
				if (dist_y eq 0) then y = coord.y[ys + pos_y] else y = y_start + py * dy
				if (dist_z eq 0) then z = coord.z[zs + pos_z] else z = z_start + pz * dz
				seeds[pos,*] = [ x, y, z ]
				pos += 1L
			end
		end
	end

	return, seeds
end

