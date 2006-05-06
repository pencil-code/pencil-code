;;
;;  $Id: monotone_array.pro,v 1.3 2006-05-06 10:03:04 ajohan Exp $
;;
;;  monotone_array - take array and return index of monotonous growth.
;;
function monotone_array, array, index

n=n_elements(array)

array2=fltarr(n)
index=intarr(n)

array2[0]=array[0]
j=0

for i=1L,n-1 do begin

  if (array[i] gt array[i-1]) then begin
;  Everything okay - value was higher than the previous.
    j=j+1
  endif else begin
;  Value lower than the previous - find out where the data point belongs.
    for k=1,j do begin
      if ( (array2[k-1] lt array[i]) and (array2[k] ge array[i]) ) then break
;  If nowhere found, then the point belongs at the first place of the array.
      if (k eq j) then begin
        k=0
        break
      endif
    endfor
    j=k
  endelse

  array2[j]=array[i]
  index[j]=i

endfor

return, index[0:j]

end
