;;
;;  $Id: monotone_array.pro,v 1.1 2006-03-15 15:12:46 ajohan Exp $
;;
;;  monotone_array - take array and return index of monotonous growth.
;;
function monotone_array, array, index

n=n_elements(array)

index=[0]
ilastmax=0

for i=1L,n-1 do begin

  if (array[i] gt array[ilastmax]) then begin
    index=[index,i]
    ilastmax=i
  endif

endfor

return, index

end
