; This function returns a given tag from a structure, if it exists.
;
; Parameters:
; structure     data structure
; tag           name of the tag to search for
;
; Returns:
; content       if tag exists
; NaN           otherwise
;
; History:
; 26.08.2016, PABourdin: coded

function get_tag, structure, tag

	COMPILE_OPT IDL2, HIDDEN

	index = find_tag (structure, tag)
	if (index lt 0) then return, !Values.D_NaN
	return, structure.(index)
end

