; Description:
;   Spread a scalar quantity into several components (default: 3 components).
;
; Examples:
;   data_3D = spread_scalar_to_components (scalar quantity)
;   data_7D = spread_scalar_to_components (scalar quantity, components=7)

function spread_scalar_to_components, scalar, components=components

	if (not keyword_set (components)) then components = 3

	data = scalar
	if (size (data, /n_dimensions) eq 1) then data = spread (data, 1, 1)
	if (size (data, /n_dimensions) eq 2) then data = spread (data, 2, 1)
	if (size (data, /n_dimensions) eq 3) then data = spread (data, 3, components)

	return, data
end

