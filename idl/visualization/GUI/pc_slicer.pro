;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_slicer.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;  $Id$
;;;
;;;  Description:
;;;    Tool for creating a configurable 2D slice of the given 3D quantity.
;;;


; Event handling of slicer window
pro pc_slicer_event, event

	common pc_slicer_common, x_anchor, y_anchor, z_anchor, theta, phi, cube, slice, coords, dims, slice_grid
	common pc_slicer_GUI_common, win, f_theta, s_theta, f_phi, s_phi, slice_sx, slice_sy

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)
	WIDGET_CONTROL, event.id, GET_UVALUE = eventval

	quit = -1
	DRAW_SLICE = 0

	SWITCH eventval of
	'THETA': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		theta = (pos > 0.0) < 360.0
		WIDGET_CONTROL, s_theta, SET_VALUE = theta
		DRAW_SLICE = 1
		break
	end
	'PHI': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		phi = (pos > (-180.0)) < 180.0
		WIDGET_CONTROL, s_phi, SET_VALUE = phi
		DRAW_SLICE = 1
		break
	end
	'SL_THETA': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		theta = (pos > 0.0) < 360.0
		WIDGET_CONTROL, f_theta, SET_VALUE = theta
		DRAW_SLICE = 1
		break
	end
	'SL_PHI': begin
		WIDGET_CONTROL, event.id, GET_VALUE = pos
		phi = (pos > (-180.0)) < 180.0
		WIDGET_CONTROL, f_phi, SET_VALUE = phi
		DRAW_SLICE = 1
		break
	end
	'RESET': begin
		pc_slicer_reset
		DRAW_SLICE = 1
		break
	end
	'IMAGE': begin
		WIDGET_CONTROL, event.id, SENSITIVE = 0
; WORK HERE:
;		pc_save_image, file_name+".png", window=win
;		save, z, t, prof_name, prof_mean, prof_min, prof_max, x_range, x_label, z_range, z_label, file=file_name+".xdr"
		WIDGET_CONTROL, event.id, SENSITIVE = 1
		break
	end
	'QUIT': begin
		quit = event.top
		win = -1
		break
	end
	endswitch

	if (DRAW_SLICE) then pc_slicer_draw

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	if (quit ge 0) then WIDGET_CONTROL, quit, /DESTROY

	return
end


; Draw the timeseries plots
pro pc_slicer_draw

	common pc_slicer_common, x_anchor, y_anchor, z_anchor, theta, phi, cube, slice, coords, dims, slice_grid
	common pc_slicer_GUI_common, win, f_theta, s_theta, f_phi, s_phi, slice_sx, slice_sy

	if (win lt 0) then return
	wset, win

;	slice[*,*] = !Values.D_NaN
;	slice_2D = pc_slice_2D (cube, [ x_anchor, y_anchor, z_anchor ], theta, phi, slice_grid=slice_grid, grid=coords, dim=dims)
;	if ((size (slice_2D, /n_dimensions) eq 1) and (slice_2D[0] eq -1)) then begin
;		print, "WARNING: Failed to create the slice for anchor:", x_anchor, y_anchor, z_anchor, " and theta,phi=", theta, phi
;		return
;	end
;	if ((slice_sx le 0.0) or (slice_sy le 0.0)) then begin
;		slice[0:slice_sx-1,0:slice_sy-1] = slice_2D
;	end else begin
;		slice[0:slice_sx-1,0:slice_sy-1] = congrid (slice_2D, slice_sx, slice_sy, 1, /center)
;	end
	slice = pc_slice_2D (cube, [ x_anchor, y_anchor, z_anchor ], theta, phi, slice_grid=slice_grid, grid=coords, dim=dims)
	if ((slice_sx gt 0.0) or (slice_sy gt 0.0)) then slice = congrid (slice, slice_sx, slice_sy, 1, /center)
	tvscl, slice, /NaN
end


; Update slice
pro pc_slicer_update, anchor

	common pc_slicer_common, x_anchor, y_anchor, z_anchor, theta, phi, cube, slice, coords, dims, slice_grid
	common pc_slicer_GUI_common, win, f_theta, s_theta, f_phi, s_phi, slice_sx, slice_sy

	x_anchor = anchor[0]
	y_anchor = anchor[1]
	z_anchor = anchor[2]

	if (n_elements (win) eq 0) then win = -1
	if (win ge 0) then pc_slicer_draw
end


; Reset to defaults
pro pc_slicer_reset

	common pc_slicer_common, x_anchor, y_anchor, z_anchor, theta, phi, cube, slice, coords, dims, slice_grid
	common pc_slicer_GUI_common, win, f_theta, s_theta, f_phi, s_phi, slice_sx, slice_sy

	theta = 0.0
	phi = 0.0

	if (win ge 0) then begin
		WIDGET_CONTROL, f_theta, SET_VALUE = theta
		WIDGET_CONTROL, f_phi, SET_VALUE = phi
		WIDGET_CONTROL, s_theta, SET_VALUE = theta
		WIDGET_CONTROL, s_phi, SET_VALUE = phi
		pc_slicer_draw
	end
end


; Calculate and draw a 2D slice of the given 3D data
;
; data:        3D data cube (can be including ghost cells)
; grid:        grid coordinate structure of the given 3D data cube,
;              asumed to be in the center of the data cube (eg. without ghost cells).
; dim:         dim coordinate structure of the given 3D data cube
; anchor:      anchor point for the rotation
; zoom:        magnification factor
;
pro pc_slicer, data, grid=grid, dim=dim, anchor=anchor, zoom=zoom

	common pc_slicer_common, x_anchor, y_anchor, z_anchor, theta, phi, cube, slice, coords, dims, slice_grid
	common pc_slicer_GUI_common, win, f_theta, s_theta, f_phi, s_phi, slice_sx, slice_sy

	if (n_elements (data) eq 0) then message, "You need to give some data to get slices."
	cube = data

	if (n_elements (anchor) ne 3) then message, "You need to give an anchor point."
	x_anchor = anchor[0]
	y_anchor = anchor[1]
	z_anchor = anchor[2]

	if (n_elements (dim) eq 0) then pc_read_dim, obj=dim
	if (n_elements (grid) eq 0) then pc_read_grid, obj=grid, dim=dim
	dims = dim
	coords = grid

	if (n_elements (file_name) eq 0) then file_name = ""
	if (file_name eq "") then file_name = "slicer"
	pos = stregex (file_name, '[ \/]+', length=len)
	while (pos gt 0) do begin
		file_name = strmid (file_name, 0, pos) + "_" + strmid (file_name, pos+len)
		pos = stregex (file_name, '[ \/]+', length=len)
	end
	file_label = file_name

	; GUI default settings
	win = -1
	pc_slicer_reset
	sl_width = 360

	; Determine size of the produced slice
	slice = pc_slice_2D (cube, [ x_anchor, y_anchor, z_anchor ], theta, phi, slice_grid=slice_grid, grid=coords, dim=dims)
	s = size (slice)
	slice_sx = s[1]
	slice_sy = s[2]
	if (n_elements (zoom) eq 1) then begin
		if ((zoom ne 1.0) and (zoom gt 0.0)) then begin
			slice_sx = round (slice_sx * zoom)
			slice_sy = round (slice_sy * zoom)
		end
	end

	MOTHER	= WIDGET_BASE (title="PC slicer")
	APP	= WIDGET_BASE (MOTHER, /col)

	BASE	= WIDGET_BASE (APP, /row)

	CTRL	= WIDGET_BASE (BASE, /row)

	BUT	= WIDGET_BASE (CTRL, /col)
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='RESET', uvalue='RESET')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='SAVE SLICE', uvalue='IMAGE')
	tmp	= WIDGET_BUTTON (BUT, xsize=100, value='CLOSE', uvalue='QUIT')

	SLIDE	= WIDGET_BASE (CTRL, /col)
	BUT	= WIDGET_BASE (SLIDE, /row, /base_align_center)
	f_theta	= CW_FIELD (BUT, title='THETA:', uvalue='THETA', value=theta, /floating, /return_events, xsize=12)
	s_theta	= CW_FSLIDER (BUT, uvalue='SL_THETA', /double, /suppress_value, min=0.0, max=360.0, /drag, value=0.0, xsize=sl_width)
	BUT	= WIDGET_BASE (SLIDE, /row, /base_align_center)
	f_phi	= CW_FIELD (BUT, title='  PHI:', uvalue='PHI', value=phi, /floating, /return_events, xsize=12)
	s_phi	= CW_FSLIDER (BUT, uvalue='SL_PHI', /double, /suppress_value, min=-180.0, max=180.0, /drag, value=0.0, xsize=sl_width)

	tmp	= WIDGET_BASE (BASE, /row)
	BUT	= WIDGET_BASE (tmp, /col)

	BASE	= WIDGET_BASE (APP, /row)

	PLOTS	= WIDGET_BASE (BASE, /col)
	tmp	= WIDGET_BASE (PLOTS, /row)
	d_prof	= WIDGET_DRAW (tmp, xsize=slice_sx, ysize=slice_sy, retain=2)

	BASE	= WIDGET_BASE (APP, /row)


	WIDGET_CONTROL, MOTHER, /REALIZE
	WIDGET_CONTROL, d_prof, GET_VALUE = win

	WIDGET_CONTROL, BASE

	pc_slicer_reset

	XMANAGER, "pc_slicer", MOTHER, /no_block

	pc_slicer_draw

end

