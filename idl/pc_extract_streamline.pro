;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_extract_streamline.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Extraction of any quantiy along a previously traced streamline.
;
;  Parameters:
;   * data           Data cube of a scalar field or vector field (3- or 4-dimensional).
;   * indices        Data indices of the traced streamline returned by 'pc_get_streamline'.
;
;  Returns:
;   * quantity       Extracted (interpolated) data values along the traced streamline.
;
;  Example:
;  ========
;
;   Load varfile and extract Temperature along a magnetic filedline:
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid
;   IDL> B = pc_get_quantity ('B', var, tags)
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> indices = pc_get_streamline (B, anchor=[2.0, 3.5, 1.2], grid=grid, distance=distance, length=length)
;   IDL> Temp_streamline = pc_extract_streamline (Temp, indices)
;   IDL> B_streamline = pc_extract_streamline (B, indices)
;   IDL> B_abs = B_streamline[0,*]^2 + B_streamline[1,*]^2 + B_streamline[2,*]^2
;   IDL> plot, distance, Temp_streamline, xtitle="coordinate along streamline", ytitle="temperature"
;   IDL> plot, distance, B_abs, xtitle="coordinate along streamline", ytitle="magnetic field", /ylog
;


; Calculation of streamline coordinates.
function pc_extract_streamline, data, indices

	if (n_elements (data) eq 0) then message, "ERROR: no data array given."
	if (n_elements (indices) eq 0) then message, "ERROR: no streamline indices given."

	n_dim = size (data, /n_dim)
	if ((n_dim lt 3) or (n_dim gt 4)) then message, "ERROR: data array dimension is invalid."

	s = size (data)
	nx = s[1]
	ny = s[2]
	nz = s[3]
	if (n_dim eq 3) then num = 1 else num = s[4]

	s = size (indices)
	if (s[0] le 1) then length = 1 else length = s[2]

	quantity = dblarr (num, length)

	; Follow the streamline
	for pos = 0, length - 1 do begin
		point = indices[*,pos]
		; Iterate over the data components
		for comp = 0, num - 1 do begin
			quantity[comp,pos] = interpolate (data[*,*,*,comp], point[0], point[1], point[2])
		end
	end

	return, quantity

end

