function count_polygons,polys
  npolys=0L

  for i=0L,n_elements(polys)-1 do begin 
    if polys[i] lt 0L then break
    if polys[i] eq 0L then continue
    npolys=npolys+1
    i=i+polys[i]
  end
  
  return, npolys
end
