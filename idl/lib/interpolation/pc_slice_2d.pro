;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_slice_2d.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Extracts any 2D slice from a given 3D data array by interpolation.
;   This method implements the idea of Hermann Amandus Schwarz et al. (1865),
;   and extends it from a 2D minimal bilinear surface solution to an
;   interpolation in 3D between source data points given on a cartesian grid.
;   Even though we cheat here a bit for simplicity, a real minimal bilinear
;   surface would be more complicated than this edge-interpolation method,
;   still this tri-linear cuboid interpolation is a good approximation.
;
;  Parameters:
;   * in             1D or multi-dimensional data array
;   * source         original grid structure
;   * anchor         anchor point grid coordinates [array of three components]
;   * theta          rotation angle around Z-axis (0°-360°)
;   * phi            rotation angle around the horizontal axis (0°-360°)
;  Optional parameters:
;   * dim            original dim structure
;   * slice_gird     (optional) returns the grid coordinates of the slice
;   * crop           crop parts of the slice that contain no data
;
;   The horizontal axis for phi-rotation is defined by the cut of the XZ-plane
;   after theta-rotation with the XY-plane going though the anchor point.
;   For phi={0°|180°} the cut is {parallel|anti-parallel} to the Z-axis.
;   
;
;  Example:
;  ========
;
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid, dim=dim
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> anchor = [ 1.2, 3.4, 5.6 ]
;   IDL> theta = 75.6
;   IDL> phi = 12.3
;   IDL> Temp_slice = pc_slice_2d (Temp, grid, anchor, theta, phi, dim=dim, slice_grid=slice_grid)
;


; Extract a 2D slice from 3D data cube.
function pc_slice_2d, in, source, anchor, theta, phi, zoom=zoom, dh=dh, dv=dv, nh=nh, nv=nv, slice_grid=target, dim=dim, datadir=datadir, crop=crop

	if ((n_elements (in) eq 0) or (n_elements (source) eq 0)) then begin
		; Print usage
		print, "USAGE:"
		print, "======"
		print, "slice = pc_slice_2d (data, source, anchor, theta, phi, slice_grid=target_grid)"
		print, ""
		return, -1
	end

	NaN = !Values.D_NaN
	pi_180 = double (!Pi) / 180

	in_size = size (in)
	n_dim = in_size[0]
	if ((n_dim lt 3) or (n_dim gt 4)) then message, "ERROR: the input data must be 3- or 4-dimensional."

	if (n_elements (dim) eq 0) then pc_read_dim, obj=dim, datadir=datadir
	if (in_size[1] eq dim.nxgrid) then nghostx = 0 else nghostx = dim.nghostx
	if (in_size[2] eq dim.nygrid) then nghosty = 0 else nghosty = dim.nghosty
	if (in_size[3] eq dim.nzgrid) then nghostz = 0 else nghostz = dim.nghostz

	x_size = in_size[1]
	y_size = in_size[2]
	z_size = in_size[3]
	if (n_dim eq 4) then num_comp = in_size[4] else num_comp = 1

	if (n_elements (source) eq 0) then begin
		if (not keyword_set (dh)) then dh = 1.0
		if (not keyword_set (dv)) then dv = 1.0
		L_diagonal = sqrt ((x_size/dh)^2 + (y_size/dh)^2 + (z_size/dv)^2)
	end else begin
		if (not keyword_set (dh)) then dh = mean ([ source.dx, source.dy ])
		if (not keyword_set (dv)) then dv = mean (source.dz)
		L_diagonal = sqrt (source.Lx^2 + source.Ly^2 + source.Lz^2)
	end
	if (not keyword_set (nh)) then begin
		nh = round (L_diagonal / dh)
		dh = L_diagonal / nh
	end
	if (not keyword_set (nv)) then begin
		nv = round (L_diagonal / dv)
		dv = L_diagonal / nv
	end
	if (n_elements (zoom) ne 0) then begin
		dh /= zoom
		dv /= zoom
	end
	x = dh * (dindgen (nh) - 0.5*(nh-1))
	y = x
	z = dv * (dindgen (nv) - 0.5*(nv-1))
	if (n_elements (source) eq 0) then begin
		Lx = nh * dh
		Lz = nv * dv
		source = { x:x, y:y, z:z, lequidist:[1,1,1], lperi:[0,0,0], Lx:Lx, Ly:Lx, Lz:Lz, dx:dh, dy:dh, dz:dv, dx_1:1.0/dh, dy_1:1.0/dh, dz_1:1.0/dv }
	end

	; Construct output slice
	out = make_array ([ nh, nv, num_comp ], /double, value=NaN)

	if (n_elements (target) eq 0) then begin
		; Construct grid coordinates of output slice
		target = make_array ([ nh, nv, 3 ], /double, value=NaN)

		; Rotate XZ-plane around X-axis (phi)
		phi_par = any (phi eq [-360,-180,0,180,360])
		phi_ort = any (phi eq [-270,-90,90,270])
		if (phi_par or phi_ort) then begin
			cos_phi = 1.0 - phi_ort - 2*any (phi eq [-180,180])
			sin_phi = 0.0 + phi_ort - 2*any (phi eq [-90,270])
		end else begin
			cos_phi = cos (phi * pi_180)
			sin_phi = sin (phi * pi_180)
		end

		; Rotate phi-rotated XZ-plane around Z-axis (theta)
		theta_par = any (theta eq [-360,-180,0,180,360])
		theta_ort = any (theta eq [-270,-90,90,270])
		if (theta_par or theta_ort) then begin
			cos_theta = 1.0 - theta_ort - 2*any (theta eq [-180,180])
			sin_theta = 0.0 + theta_ort - 2*any (theta eq [-90,270])
		end else begin
			cos_theta = cos (theta * pi_180)
			sin_theta = sin (theta * pi_180)
		end

		; Construct rotated slice grid
		tx = x * cos_theta
		ty = x * sin_theta
		for pv = 0, nv-1 do begin
			; X and Y coordinates
			target[*,pv,0] = tx
			target[*,pv,1] = ty
		end
		add_x = z * sin_theta * sin_phi
		add_y = z * cos_theta * sin_phi
		tz = z * cos_phi
		for ph = 0, nh-1 do begin
			; Y and Z coordinates
			target[ph,*,0] += add_x
			target[ph,*,1] -= add_y
			target[ph,*,2] = tz
		end

		; Shift to the anchor position
		target[*,*,0] += anchor[0]
		target[*,*,1] += anchor[1]
		target[*,*,2] += anchor[2]
	end

	; Interpolation to the target slice grid coordinates
	for ph = 0, nh-1 do begin
		for pv = 0, nv-1 do begin

			; Find grid coordinates closest to the target coordinate
			px = pc_find_index (target[ph,pv,0], source.x, num=x_size, /round)
			py = pc_find_index (target[ph,pv,1], source.y, num=y_size, /round)
			pz = pc_find_index (target[ph,pv,2], source.z, num=z_size, /round)

			; Target position must still be inside the source data array
			if ((px lt 0) or (py lt 0) or (pz lt 0)) then continue
			if ((px ge x_size-1) or (py ge y_size-1) or (pz ge z_size-1)) then continue

			; Find lower neighbouring grid coordinates of the target coordinate
			if ((px ge 1) and (target[ph,pv,0] lt source.x[px])) then px -= 1
			if ((py ge 1) and (target[ph,pv,1] lt source.y[py])) then py -= 1
			if ((pz ge 1) and (target[ph,pv,2] lt source.z[pz])) then pz -= 1

			; Upper plane can be interpolated by using an inner grid cell
			if (target[ph,pv,0] eq source.x[x_size-1]) then px -= 1
			if (target[ph,pv,1] eq source.y[y_size-1]) then py -= 1
			if (target[ph,pv,2] eq source.z[z_size-1]) then pz -= 1

			; Target position must still be inside the source data array
			if ((px lt 0) or (py lt 0) or (pz lt 0)) then continue

			; Fraction of the target position to the side-length of the grid cell
			fx = (target[ph,pv,0] - source.x[px]) / (source.x[px+1] - source.x[px])
			fy = (target[ph,pv,1] - source.y[py]) / (source.y[py+1] - source.y[py])
			fz = (target[ph,pv,2] - source.z[pz]) / (source.z[pz+1] - source.z[pz])

			for pa = 0, num_comp-1 do begin
				; Find interpolated values of the target projection to the 12 edges
				; (Commented out duplicate definitions of the tri-linear interpolation)
;				xl_yl = (1.0 - fz) * in[px  ,py  ,pz  ,pa] + fz * in[px  ,py  ,pz+1,pa]
;				xl_yu = (1.0 - fz) * in[px  ,py+1,pz  ,pa] + fz * in[px  ,py+1,pz+1,pa]
;				xu_yl = (1.0 - fz) * in[px+1,py  ,pz  ,pa] + fz * in[px+1,py  ,pz+1,pa]
;				xu_yu = (1.0 - fz) * in[px+1,py+1,pz  ,pa] + fz * in[px+1,py+1,pz+1,pa]
				xl_zl = (1.0 - fy) * in[px  ,py  ,pz  ,pa] + fy * in[px  ,py+1,pz  ,pa]
				xl_zu = (1.0 - fy) * in[px  ,py  ,pz+1,pa] + fy * in[px  ,py+1,pz+1,pa]
				xu_zl = (1.0 - fy) * in[px+1,py  ,pz  ,pa] + fy * in[px+1,py+1,pz  ,pa]
				xu_zu = (1.0 - fy) * in[px+1,py  ,pz+1,pa] + fy * in[px+1,py+1,pz+1,pa]
;				yl_zl = (1.0 - fx) * in[px  ,py  ,pz  ,pa] + fx * in[px+1,py  ,pz  ,pa]
;				yl_zu = (1.0 - fx) * in[px  ,py  ,pz+1,pa] + fx * in[px+1,py  ,pz+1,pa]
;				yu_zl = (1.0 - fx) * in[px  ,py+1,pz  ,pa] + fx * in[px+1,py+1,pz  ,pa]
;				yu_zu = (1.0 - fx) * in[px  ,py+1,pz+1,pa] + fx * in[px+1,py+1,pz+1,pa]

				; Find interpolated values of the target projection to the 6 sides
				xy_l = (1.0 - fx) * xl_zl + fx * xu_zl
				xy_u = (1.0 - fx) * xl_zu + fx * xu_zu
;				xy_l = (1.0 - fy) * yl_zl + fy * yu_zl
;				xy_u = (1.0 - fy) * yl_zu + fy * yu_zu
;				xz_l = (1.0 - fx) * xl_yl + fx * xu_yl
;				xz_u = (1.0 - fx) * xl_yu + fx * xu_yu
;				xz_l = (1.0 - fz) * yl_zl + fz * yl_zu
;				xz_u = (1.0 - fz) * yu_zl + fz * yu_zu
;				yz_l = (1.0 - fy) * xl_yl + fy * xu_yl
;				yz_u = (1.0 - fy) * xl_yu + fy * xu_yu
;				yz_l = (1.0 - fz) * xl_zl + fz * xl_zu
;				yz_u = (1.0 - fz) * xl_zl + fz * xl_zu

				; Find interpolated value of the target point
				out[ph,pv,pa] = (1.0 - fz) * xy_l + fz * xy_u
;				out[ph,pv,pa] = (1.0 - fy) * xz_l + fy * xz_u
;				out[ph,pv,pa] = (1.0 - fx) * yz_l + fx * yz_u
			end
		end
	end

	if (keyword_set (crop)) then begin
		for ph = 0L, nh-1 do begin
			if (any (finite (out[ph,*,*]))) then begin
				chs = ph
				break
			end
		end
		for ph = nh-1L, 0, -1 do begin
			if (any (finite (out[ph,*,*]))) then begin
				che = ph
				break
			end
		end
		if ((n_elements (chs) eq 0) or (n_elements (che) eq 0)) then return, !Values.D_NaN
		for pv = 0L, nv-1 do begin
			if (any (finite (out[*,pv,*]))) then begin
				cvs = pv
				break
			end
		end
		for pv = nv-1L, 0, -1 do begin
			if (any (finite (out[*,pv,*]))) then begin
				cve = pv
				break
			end
		end
		if ((n_elements (cvs) eq 0) or (n_elements (cve) eq 0)) then return, !Values.D_NaN
		nh = che - chs + 1
		nv = cve - cvs + 1
		out = reform (out[chs:che,cvs:cve,*], nh, nv, num_comp)
		target = reform (target[chs:che,cvs:cve,*], nh, nv, 3)
	end

	return, out[*,*,0:num_comp-1]

end

