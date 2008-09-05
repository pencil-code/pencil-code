;;
;;  $Id$
;;
;;  monotone_array - take array and return index of monotonous growth.
;;
function monotone_array, array, index

n=n_elements(array)

array2=dblarr(n)
index=lonarr(n)

array2[0]=array[0]
j=0L

for i=1L,n-1 do begin

  if (array[i] gt array[i-1]) then begin
;  Everything okay - value was higher than the previous.
    j=j+1
  endif else begin 
    if (array[i] eq array[i-1]) then begin
;  Sometimes, due to too low precision in time, there are consecutive
;  points with the same time value but different everything else.
;  in this case, print a warning, but continue.
      print,'WARNING: Equal vaules of time for consecutive points!'
      j=j+1
    endif else begin
;  Value lower than the previous - find out where the data point belongs.
      for k=1L,j do begin
        if ( (array2[k-1] lt array[i]) and (array2[k] ge array[i]) ) then break
;  If nowhere found, then the point belongs at the first place of the array.
        if (k eq j) then begin
          k=0
          break
        endif
      endfor
      j=k
    endelse
  endelse
  array2[j]=array[i]
  index[j]=i

endfor

return, index[0:j]

end
