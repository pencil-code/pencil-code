;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_find_index.pro         ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Returns the index of a 'needle' in a sorted (ascending values) 'haystack'.
;;;   'num' is the maximum index upto that the 'needle' is searched.
;;;   'sort' sorts the 'haystack' in the given order (default: ascending).
;;;   If the 'reverse' keyword is set, the 'haystack' is in descending order.
;;;   With the 'round' keyword, the index of the array element nearest to 'needle' is returned.
;;;   'upper' and 'lower' can be given as an estimate of the index range.
;;;   By default the interpolated position of 'needle' within 'haystack' is returned.

function pc_find_index, needle, haystack, num=num, sort=sort, reverse=reverse, round=round, upper=upper, lower=lower

	if (n_elements (num) eq 0) then num = n_elements (haystack)
	if (num le 0) then num = n_elements (haystack)
	default, lower, 0
	default, upper, num - 1

	if (keyword_set (sort)) then haystack = haystack[sort (haystack)]

        ; Search for lower limit, starting from estimated position
	inc = 2
	while ((lower gt 0) and (needle lt haystack[lower])) do begin
		upper = lower
		lower = (lower - inc) > 0
		inc *= 2
	end

        ; Search for upper limit, starting from estimated position
	inc = 2
	while ((upper lt num-1) and (needle gt haystack[upper])) do begin
		lower = upper
		upper = (upper + inc) < (num - 1)
		inc *= 2
	end

	if (needle gt haystack[upper]) then begin
		; needle is above range
		lower = num - 2
	end else if (needle lt haystack[lower]) then begin
		; needle is below range
		lower = 0
	end else begin
		; Find needle position
		while (lower+1 lt upper) do begin
			mid = lower + (upper - lower) / 2
			if (needle lt haystack[mid]) then begin
				upper = mid
			end else begin
				lower = mid
			end
		end
	end
	upper = lower + 1

	; Interpolate needle position
	index = lower + (needle - haystack[lower]) / (haystack[upper] - haystack[lower])

	if (keyword_set (round)) then index = round (index)
	if (keyword_set (reverse)) then begin
		if (keyword_set (sort)) then haystack = reverse (haystack)
		index = num - 1 - index
	end

	return, index
end

