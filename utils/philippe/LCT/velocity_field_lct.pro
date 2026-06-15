; Get velocity field [m/s] by Local-Correlation-Tracking
;
; step = 4 correlates images 1 with 5, 2 with 6, 3 with 7, etc. (and averages afterwards)
; times are the observation times of the frames in SOT fits files format
; m_pixel is the number of meters per pixel

; Author: Philippe-A. Bourdin, 27 September 2012.

pro velocity_field_lct, in, vx, vy, step, times, m_pixel, delta_t

	NaN = !VALUES.D_NAN

	size_x = (size (in))[1]
	size_y = (size (in))[2]
	num_frames = (size (in))[3]

	; Time dimemsion:
	dim_t = 3
	; Standard deviation for the windowing function: std_dev = FWHM/pixelsize/2.35482 [pixel]:
	std_dev = 6
	; Rebinning [pixel]:
	rebin = 2
	; Window radius [timesteps]:
	win_radius = step

	vx = dblarr (size_x, size_y, num_frames - 2 * win_radius)
	vy = dblarr (size_x, size_y, num_frames - 2 * win_radius)

	delta_t = dblarr (num_frames)
	delta_t[*] = NaN
	timestamps = get_timestamp (times)
	for pos = win_radius, num_frames - 1 - win_radius do begin
		delta_t[pos] = 0.0
		for offset = -win_radius, win_radius-step do begin
			delta_t[pos] += timestamps[pos+offset+step] - timestamps[pos+offset]
		endfor
	endfor

	; Average timestep between two steps
	delta_t /= 2 * win_radius - step + 1

	for pos = step, num_frames - 1 - step do begin
		print, "POS: ", pos
		flow_reb_cube, in[*,*,pos-step:pos+step], dim_t, step, std_dev, rebin, tmp_vx, tmp_vy
		vx[*,*,pos - step] = tmp_vx / delta_t[pos] * m_pixel
		vy[*,*,pos - step] = tmp_vy / delta_t[pos] * m_pixel
	endfor

	; Average timestep between two subsequent frames
	delta_t /= step

end

