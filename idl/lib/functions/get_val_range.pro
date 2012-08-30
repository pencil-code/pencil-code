;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   get_val_range.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Get value range and extend it, if necessary.
;;;   IDL sliders won't accept min/max-settings that are equal.
;;;   Also for automated plotting, this might help to set good value ranges.
;;;

function get_val_range, data

	tmp = minmax (data)
	min = tmp[0]
	max = tmp[1]

	if (min eq max) then begin
		; extend value range a little, if necessary (must have min < max)
		; uniform data should stay in the middle of min and max
		if (min eq 0.0) then begin
			min = -1d-42
			max = 1d-42
		end else begin
			max *= 1.00001
			min -= max - min
		end
	end

	return, [min, max]
end

