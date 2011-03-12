;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   find_index.pro         ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Returns the index of a value in a sorted data array.
;;;   'num' is the maximum index upto that the value is searched.

function find_index, value, data, num

	pixel = 0
	if (num gt 1) then begin
		for i = 0, num-2 do begin
			if (value lt 0.5*(data[i]+data[i+1])) then begin
				pixel = i
				break
			endif
		end
		if (value ge 0.5*(data[num-2]+data[num-1])) then pixel = num - 1
	endif

	return, pixel
end
