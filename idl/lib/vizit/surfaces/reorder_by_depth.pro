pro free_polygon_array,poly_array

  npolys=n_elements(poly_array)

  for i=0L,npolys-1 do begin 
    if ptr_valid(poly_array[i]) then ptr_free,poly_array[i]
  end
end

function polygon_array_to_list_in_order,poly_array,order

  npolys=n_elements(order)
  if (npolys eq 0) then return, [0]

  polys=*(poly_array[order[0]])
  ptr_free,poly_array[order[0]]

  tenpercent=long(npolys/10.)

  for i=1L,npolys-1 do begin 
    ipoly=order[i]
    polys=[polys,*(poly_array[ipoly])]
    ptr_free,poly_array[ipoly]
  end

  return,polys
end

function polygon_array_to_list,poly_array

  npolys=n_elements(poly_array)
  if (npolys eq 0) then return, [0]

  polys=*(poly_array[0])
  ptr_free,poly_array[0]

  for i=1L,npolys-1 do begin 
    polys=[polys,*(poly_array[i])]
    ptr_free,poly_array[i]
    print,i
  end

  return,polys
end

function polygon_list_to_array,polys

  npolys=count_polygons(polys)
  poly_array=ptrarr(npolys) 

  ipoly=0L
  for i=0L,n_elements(polys)-1 do begin 
    if polys[i] lt 0L then break
    if polys[i] eq 0L then continue
    poly_array[ipoly]=ptr_new(polys[i:i+polys[i]])
  
    ipoly=ipoly+1
    i=i+polys[i]
  end
  return,poly_array
end

function get_polygon,polys,n

  ipoly=0L
  for i=0L,n_elements(polys)-1 do begin 
    if polys[i] lt 0L then break
    if polys[i] eq 0L then continue
  
    if ipoly eq n then begin
      return,polys[i:i+polys[i]]
    endif
  
    ipoly=ipoly+1
    i=i+polys[i]
  end

  return,[0]
end

function get_polygon,polys,n

  ipoly=0L
  for i=0L,n_elements(polys)-1 do begin 
    if polys[i] lt 0L then break
    if polys[i] eq 0L then continue
  
    if ipoly eq n then begin
      return,polys[i:i+polys[i]]
    endif
  
    ipoly=ipoly+1
    i=i+polys[i]
  end

  return,[0]
end

pro reorder_by_depth,object,tmat=tmat
  
  object->GetProperty, DATA=verts, POLYGONS=polys, $
                       XCOORD_CONV=xcc, $
                       YCOORD_CONV=ycc, $
                       ZCOORD_CONV=zcc

  if (keyword_set(tmat)) then tmat=reform(tmat)
  
  npolys=count_polygons(polys)
  distance=fltarr(npolys) 
 
  normverts=verts
  normverts[*,0]=xcc[0]+verts[*,0]*xcc[1]
  normverts[*,1]=ycc[0]+verts[*,1]*ycc[1]
  normverts[*,2]=zcc[0]+verts[*,2]*zcc[1]
  
  
  if keyword_set(tmat) then begin
    verts=vert_t3d(normverts,matrix=tmat,/no_copy)
  endif else begin
    verts=normverts
    normverts=0
  endelse
  
  vert_dist=reform(verts[2,*])
  
  ipoly=0L
  for i=0L,n_elements(polys)-1 do begin 
    if polys[i] lt 0L then break
    if polys[i] eq 0L then continue
  
    meandist=0.
    for j=i+1,i+polys[i] do begin
      meandist=meandist+vert_dist[polys[j]]  
    end
    meandist=meandist/polys[i]
    distance[ipoly]=meandist
  
    ipoly=ipoly+1
    i=i+polys[i]
  end
  
  order=sort(distance)

  poly_array=polygon_list_to_array(polys)
  newpolys=polygon_array_to_list_in_order(poly_array,order)

  free_polygon_array,poly_array
;  newpolys=get_polygon(polys,order[0L])
;  for ipoly=1L,npolys-1 do begin
;    newpolys=[newpolys,get_polygon(polys,order[ipoly])]
;  end

  object->SetProperty, POLYGONS=newpolys
  
end


