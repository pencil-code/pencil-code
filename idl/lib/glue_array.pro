function glue_array, array1, array2, array3

size1=size(array1)
size2=size(array1)
size3=size(array1)

i=where( (size1 ne size2) or (size2 ne size3) )

if (max(i) gt -1) then begin
  print, 'ERROR: arrays must have same size'
  stop
endif

nx=n_elements(array1[*,0,0])
ny=n_elements(array1[0,*,0])
nz=n_elements(array1[0,0,*])

array=dblarr(nx,ny,nz,3)
array[*,*,*,0]=array1
array[*,*,*,1]=array2
array[*,*,*,2]=array3

return, array

end
