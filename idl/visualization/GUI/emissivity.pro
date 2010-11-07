;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   emissivity.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Fast and simple to use tool to view and compare emissivities of different ions.
;;;   This tool expects, that '.r analyse' has already been executed.
;;;
;;;  To do:
;;;   Add more comments


;;; Settings:

; Quantities to be used for emissivity (calculated in 'precalc_data'):
emissivity_quantities = { temperature:'Temp', logarithmic_density:'ln_rho' }


; Event handling of emissivity visualisation window
pro emissivity_event, event

	common emissive_common, parameter, selected_emissivity, em, em_x, em_y, em_z, cut_z, sub_horiz, aver_z, emin, emax
	common emigui_common, wem_x, wem_y, wem_z, val_t, val_b, sl_min, sl_max

	WIDGET_CONTROL, WIDGET_INFO(event.top, /CHILD)

	quit = -1
	DRAW_IMAGES = 0

	WIDGET_CONTROL, event.id, GET_UVALUE = eventval


	CASE eventval of
	'CUT_Z':  begin
		WIDGET_CONTROL, event.id, GET_VALUE = idx
		cut_z = idx
		DRAW_IMAGES = 1
	end
	'HORIZ':  begin
		sub_horiz = event.select
		precalc_emissivity
		DRAW_IMAGES = 1
	end
	'VAL_B': begin
		WIDGET_CONTROL, val_b, GET_VALUE = sl_min
		if (sl_min gt sl_max) then begin
			sl_min = sl_max
			WIDGET_CONTROL, val_b, SET_VALUE = sl_min
		end
		emin = 10^sl_min
		DRAW_IMAGES = 1
	end
	'VAL_T': begin
		WIDGET_CONTROL, val_t, GET_VALUE = sl_max
		if (sl_max lt sl_min) then begin
			sl_max = sl_min
			WIDGET_CONTROL, val_t, SET_VALUE = sl_max
		end
		emax = 10^sl_max
		DRAW_IMAGES = 1
	end
	'EMIS': begin
		last = selected_emissivity
		selected_emissivity = event.index
		if (last ne selected_emissivity) then precalc_emissivity
		DRAW_IMAGES = 1
	end
	'QUIT': begin
		quit = event.top
	end
	endcase

	if (DRAW_IMAGES) then plot_emissivity

	WIDGET_CONTROL, WIDGET_INFO (event.top, /CHILD)

	IF quit GE 0 THEN  WIDGET_CONTROL, quit, /DESTROY

	return
end


; Calculates emissivities
pro precalc_emissivity

	common emissive_common, parameter, selected_emissivity, em, em_x, em_y, em_z, cut_z, sub_horiz, aver_z, emin, emax
	common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, sources
	common emigui_common, wem_x, wem_y, wem_z, val_t, val_b, sl_min, sl_max
	common settings_common, idx1, idx2, idx3, cut, abs_scale, show_cross, show_cuts, sub_aver, selected_cube, selected_overplot, selected_snapshot, af_1, af_2, af_3
	common slider_common, bin1, bin2, bin3, num_1, num_2, num_3, pos_b, pos_t, csmin, csmax, dimensionality

	T_0 = parameter[selected_emissivity].T_ex
	dT = parameter[selected_emissivity].delta_T

;	for iy = 0, num_1-1 do em_x = em_x * exp (rho_0 - ((varsets[selected_snapshot].ln_rho[cut])[*,iy,cut_z:num_3-1] > rho_0)) + em[*,iy,cut_z:num_3-1]
;	for ix = 0, num_2-1 do em_y = em_y * exp (rho_0 - ((varsets[selected_snapshot].ln_rho[cut])[ix,*,cut_z:num_3-1] > rho_0)) + em[ix,*,cut_z:num_3-1]
;	for iz = cut_z, num_3-1 do em_z = em_z * exp (rho_0 - ((varsets[selected_snapshot].ln_rho[cut])[*,*,iz] > rho_0)) + em[*,*,iz]

;	em = (1 - cos (((1 - ((alog10 (varsets[selected_snapshot].temp[cut]) - T_0) / dT)^2) > 0) * !PI)) * exp ((varsets[selected_snapshot].ln_rho[cut]) * 2)
	em = ((1 - ((alog10 (varsets[selected_snapshot].temp[cut]) - T_0) / dT)^2) > 0) * exp ((varsets[selected_snapshot].ln_rho[cut]) * 2)

	em_x = total (em, 1)
	em_y = total (em, 2)
	em_z = total (em[*,*,cut_z:*], 3)

	; normalise to averages of maximum emissivities in horizontal layers
	if (sub_horiz) then begin
		m = 0
		std = 10^mean ([sl_min, sl_max])
		for z=cut_z, num_3-1 do begin
			zs = min ([num_3/2, num_3-aver_z-1])
			if (z lt zs) then begin
				m = 0
				for iz=zs, zs+aver_z do m += max (em_x[*,iz]) + max (em_y[*,iz])
				m /= aver_z * 2
			end
			if (m gt 0) then begin
				em_x[*,z] *= std / m
				em_y[*,z] *= std / m
			end
		end
	end

	if (bin1 ne 1 or bin3 ne 1) then em_x = congrid (em_x, fix (num_1*bin1), fix ((num_3-cut_z)*bin3), cubic = 0)
	if (bin2 ne 1 or bin3 ne 1) then em_y = congrid (em_y, fix (num_2*bin2), fix ((num_3-cut_z)*bin3), cubic = 0)
	if (bin1 ne 1 or bin2 ne 1) then em_z = congrid (em_z, fix (num_1*bin1), fix (num_2*bin2), cubic = 0)
end


; Plots integrated emissivities in x-, y- and z-direction
pro plot_emissivity

	common emissive_common, parameter, selected_emissivity, em, em_x, em_y, em_z, cut_z, sub_horiz, aver_z, emin, emax
	common emigui_common, wem_x, wem_y, wem_z, val_t, val_b, sl_min, sl_max
	common slider_common, bin1, bin2, bin3, num_1, num_2, num_3, pos_b, pos_t, csmin, csmax, dimensionality

	wset, wem_x
	tvscl, (em_x[*,cut_z:*] > emin) < emax, 0, cut_z

	wset, wem_y
	tvscl, (em_y[*,cut_z:*] > emin) < emax, 0, cut_z

	wset, wem_z
	tvscl, (em_z > emin) < emax
end


; Emissivitiy plotting GUI
pro emissivity, sets, limits, scaling=scaling

	common emissive_common, parameter, selected_emissivity, em, em_x, em_y, em_z, cut_z, sub_horiz, aver_z, emin, emax
	common emigui_common, wem_x, wem_y, wem_z, val_t, val_b, sl_min, sl_max
	common slider_common, bin1, bin2, bin3, num_1, num_2, num_3, pos_b, pos_t, csmin, csmax, dimensionality

	; Emissivities for different ions (values with only one decimal digit are UNVERIFIED)
	; Temperatures are logarithmic to the base of 10
	parameter = [ $
		{ title:'Si II',   lambda:1533, T_ion:4.60, T_ex:4.36, T_DEM:4.08, T_MHD:4.25, delta_T:0.3  }, $
		{ title:'Si IV',   lambda:1394, T_ion:4.81, T_ex:4.85, T_DEM:4.82, T_MHD:4.90, delta_T:0.29 }, $
		{ title:'C II',    lambda:1335, T_ion:4.67, T_ex:4.57, T_DEM:4.23, T_MHD:4.64, delta_T:0.22 }, $
		{ title:'C III',   lambda:977,  T_ion:4.78, T_ex:4.8,  T_DEM:4.71, T_MHD:4.84, delta_T:0.29 }, $
		{ title:'C IV',    lambda:1548, T_ion:5.00, T_ex:5.00, T_DEM:5.02, T_MHD:5.11, delta_T:0.25 }, $
		{ title:'O IV',    lambda:1401, T_ion:5.27, T_ex:5.14, T_DEM:5.15, T_MHD:5.18, delta_T:0.32 }, $
		{ title:'O V',     lambda:630,  T_ion:5.38, T_ex:5.35, T_DEM:5.40, T_MHD:5.44, delta_T:0.28 }, $
		{ title:'O VI',    lambda:1032, T_ion:5.45, T_ex:5.44, T_DEM:5.50, T_MHD:5.60, delta_T:0.23 }, $
		{ title:'Ne VIII', lambda:770,  T_ion:5.81, T_ex:5.76, T_DEM:5.82, T_MHD:5.89, delta_T:0.16 }, $
		{ title:'Mg X',    lambda:625,  T_ion:6.04, T_ex:6.01, T_DEM:6.01, T_MHD:6.06, delta_T:0.17 } $
	      ]
	n_emissivities = n_elements (parameter)

	; SETTINGS/DEFAULTS:
	selected_emissivity = 0
	cut_z = 0
	emin = -30.0
	emax =  10.0
	sl_min = -22.0
	sl_max = -18.0
	sub_horiz = 0
	aver_z = 10

	em = fltarr (num_1,num_2,num_3)
	em_x = fltarr (num_2,num_3)
	em_y = fltarr (num_1,num_3)
	em_z = fltarr (num_1,num_2)

	MOTHER	= WIDGET_BASE (title='emissivitiy')
	BASE    = WIDGET_BASE (MOTHER, /col)
	TOP     = WIDGET_BASE (base, /row)
	col     = WIDGET_BASE (top, /col)
	tmp     = WIDGET_SLIDER (col, uvalue='CUT_Z', value=cut_z, min=0, max=num_1-1, xsize=num_3*bin3, /drag)
	col     = WIDGET_BASE (top, /col)
	emis    = WIDGET_DROPLIST (col, value=(parameter[*].title), uvalue='EMIS', EVENT_PRO=emissivity_event, title='ion')
	col     = WIDGET_BASE (top, /col)
	b_sub   = CW_BGROUP (col, 'normalise averages', /nonexcl, uvalue='HORIZ', set_value=sub_horiz)
	col     = WIDGET_BASE (top, /col)
	tmp	= WIDGET_BUTTON (col, value='QUIT', UVALUE='QUIT', xsize=100)
	drow    = WIDGET_BASE (BASE, /row)
	tmp     = WIDGET_DRAW (drow, UVALUE='EM_X', xsize=num_2*bin2, ysize=num_3*bin3)
	WIDGET_CONTROL, tmp, /REALIZE
	wem_x   = !d.window
	tmp     = WIDGET_DRAW (drow, UVALUE='EM_Y', xsize=num_1*bin1, ysize=num_3*bin3)
	WIDGET_CONTROL, tmp, /REALIZE
	wem_y   = !d.window
	tmp     = WIDGET_DRAW (drow, UVALUE='EM_Z', xsize=num_1*bin1, ysize=num_2*bin2)
	WIDGET_CONTROL, tmp, /REALIZE
	wem_z   = !d.window
	TOP2    = WIDGET_BASE (base, /col)

	val_b   = CW_FSLIDER (TOP2, uvalue='VAL_B', /edit, min=emin, max=emax, drag=1, value=sl_min, xsize=(2*num_1*bin1+num_2*bin2)>max([num_1,num_2,num_3]) )
	val_t   = CW_FSLIDER (TOP2, uvalue='VAL_T', /edit, min=emin, max=emax, drag=1, value=sl_max, xsize=(2*num_1*bin1+num_2*bin2)>max([num_1,num_2,num_3]) )

	WIDGET_CONTROL, MOTHER, /REALIZE
	wimg = !d.window

	WIDGET_CONTROL, BASE

	XMANAGER, "emissivity", MOTHER, /no_block

	emin = 10^sl_min
	emax = 10^sl_max

	precalc_emissivity
	plot_emissivity

	return
end


emissivity, emissivity_quantities, lmn12, scaling=scaling

window, 0, xsize=8, ysize=8, retain=2
!P.MULTI = [0, 1, 1]
wdelete

end

