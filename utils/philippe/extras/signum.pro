; Get sign of a scalar or an array
;
; Parameters:
; data          scalar or array
;
; Returns:
; -1            where data is negative
; 0             where data is zero
; +1            where data is positive
;
; History:
; 15.06.2026, PABourdin: coded

function signum, data

	compile_opt idl2, hidden
	on_error, 2

	return, -fix (data lt 0.0d0) + fix (data gt 0.0d0)
end

