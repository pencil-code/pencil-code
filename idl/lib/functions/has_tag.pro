; This function checks for existence of a given tag with a structure.
;
; Parameters:
; structure     data structure
; tag           name of the tag to search for
;
; Returns:
; true          if tag exists
; flase         otherwise
;
; History:
; 24.10.2015, PABourdin: coded

function has_tag, structure, tag

	COMPILE_OPT IDL2, HIDDEN

	exists = strupcase (tag_names (structure))
	return, any (exists eq strupcase (tag))
end

