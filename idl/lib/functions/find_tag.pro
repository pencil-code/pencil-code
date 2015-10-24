; This function finds a tag within a structure.
;
; Parameters:
; structure     data structure
; tag           name of the tag to search for
;
; Returns:
; >= 0          position of tag, if tag exists
; -1            otherwise
;
; History:
; 24.10.2015, PABourdin: coded

function find_tag, structure, tag

	COMPILE_OPT IDL2, HIDDEN

	found = where (strupcase (tag_names (structure)) eq strupcase (tag))
	return, min (found)
end

