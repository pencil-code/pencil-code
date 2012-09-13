;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_destretch.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  $Id$
;
;  Description:
;   Destretches a given data array to an equidistant grid by interpolation.
;
;  Parameters:
;   * in             1D or multi-dimensional data array
;   * source         original grid
;   * target         optional: array=destination grid, scalar=destination grid spacing
;                    returns the equidistant target grid
;   * dim            Dimension to be destretched
;
;  Example:
;  ========
;
;   IDL> pc_read_var_raw, obj=var, tags=tags, grid=grid, dim=dim
;   IDL> Temp = pc_get_quantity ('Temp', var, tags)
;   IDL> z = grid.z[dim.n1:dim.n2]
;   IDL> Temp_equidist = pc_destretch (Temp, z, target=z_equidist, dim=3)
;
;  Credits:
;   Thanks to the cabin crew of LH962 for snacks and drinks...


; Destretch data from non-uniform to uniform grid.
function pc_destretch, in, source, target=target, dim=dim

	if ((n_elements (in) eq 0) or (n_elements (source) eq 0)) then begin
		; Print usage
		print, "USAGE:"
		print, "======"
		print, "destretched = pc_destretch (data, source_grid, target=destination_grid)"
		print, ""
		return, -1
	end

	in_size = size (in)
	in_layers = in_size[dim]
	n_dim = in_size[0]
	if (not keyword_set (dim)) then dim = 1
	if ((dim lt 1) or (dim gt n_dim)) then message, "ERROR: given dimension is invalid, should be between 1 and "+strtrim (n_dim, 2)+"."

	if (n_elements (source) eq 0) then message, "ERROR: the source grid must be defined."
	source_size = size (source)
	if (source_size[0] ne 1) then message, "ERROR: the source grid must be 1-dimensional."
	if (n_elements (source) ne in_layers) then message, "ERROR: the source grid doesn't fit to the given data."

	out_size = in_size
	if (n_elements (target) eq 0) then target = dindgen (in_layers) / (in_layers - 1.0) * max (source[in_layers-1] - source[0]) + source[0]
	if (n_elements (target) eq 1) then target = dindgen (in_layers) * target + source[0]
	out_size[dim] = n_elements (target)
	out = make_array (out_size[1:n_dim], type=out_size[n_dim+1])

	dim_str = ""
	for pos = 1, n_dim do begin
		if (pos eq dim) then dim_str += ",layer_pos" else dim_str += ",*"
	end
	dim_str = strmid (dim_str, 1)

	num_layers = (size (out))[dim]
	for pos = 0, num_layers-1 do begin

		pos_below = find_array_index (source, target[pos], /lower)
		layer_pos = pos_below
		result = execute ("below = in["+dim_str+"]")
		if (not result) then message, "ERROR: can't get below layer "+strtrim (layer_pos, 2)+"."

		pos_above = (pos_below + 1) < (in_layers - 1)
		layer_pos = pos_above
		result = execute ("above = in["+dim_str+"]")
		if (not result) then message, "ERROR: can't get above layer "+strtrim (layer_pos, 2)+"."

		above_weight = (target[pos] - source[pos_below]) / (source[pos_above] - source[pos_below])
		above_weight = (above_weight > 0.0) < 1.0
		below_weight = 1.0 - above_weight
		interpolated = below * below_weight + above * above_weight
		layer_pos = pos
		result = execute ("out["+dim_str+"] = interpolated")
		if (not result) then message, "ERROR: can't write layer "+strtrim (layer_pos, 2)+"."
	end

	return, out

end

