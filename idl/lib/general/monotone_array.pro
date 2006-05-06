;;
;;  $Id: monotone_array.pro,v 1.2 2006-05-06 09:44:01 ajohan Exp $
;;
;;  monotone_array - take array and return index of monotonous growth.
;;
function monotone_array, array, index

n=n_elements(array)

array2=intarr(n)
index=intarr(n)

array2[0]=array[0]
j=0

for i=1L,n-1 do begin

  if (array[i] gt array[i-1]) then begin
    j=j+1
  endif else begin
    for k=1,j do begin
      if ( (array2[k-1] lt array[i]) and (array2[k] ge array[i]) ) then break
    endfor
    j=k
  endelse

  array2[j]=array[i]
  index[j]=i

endfor

return, index[0:j]

end
