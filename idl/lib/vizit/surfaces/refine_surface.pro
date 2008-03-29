pro refine_surface,surf

  has_normals=PTR_VALID(surf.normals)
  newvertices=*(surf.vertices)
  polys=*(surf.triangles)
  if (has_normals) then newnormals=*(surf.normals)

  npolys=0L
  nvertices=(size(newvertices))[2]

  newpolys=0L

  for i=0L,n_elements(polys)-1 do begin 
    if polys[i] lt 0L then break
    if polys[i] eq 0L then continue
    npolys=npolys+1

; Refine only if it's a triangle!
    if (polys[i] eq 3) then begin    
     
      v1=polys[i+1]
      v2=polys[i+2]
      v3=polys[i+3]
      v4=nvertices
      v5=nvertices+1
      v6=nvertices+2
      refinedpolys=[3,v1,v4,v6, 3,v4,v2,v5, 3,v4,v5,v6, 3,v6,v5,v3]*1L

      verts=fltarr(3,3)
      verts[0:2,0]=(newvertices[0:2,v1]+newvertices[0:2,v2]) / 2.
      verts[0:2,1]=(newvertices[0:2,v2]+newvertices[0:2,v3]) / 2.
      verts[0:2,2]=(newvertices[0:2,v3]+newvertices[0:2,v1]) / 2.
      newvertices=transpose([transpose(newvertices),transpose(verts)])
      nvertices=nvertices+3
      if has_normals then begin
        norms=fltarr(3,3)
        norms[0:2,0]=(newnormals[0:2,v1]+newnormals[0:2,v2]) / 2.
        norms[0:2,1]=(newnormals[0:2,v2]+newnormals[0:2,v3]) / 2.
        norms[0:2,2]=(newnormals[0:2,v3]+newnormals[0:2,v1]) / 2.
        newnormals=transpose([transpose(newnormals),transpose(norms)])
      endif

      if n_elements(newpolys) gt 1 then begin
        newpolys=[newpolys,refinedpolys]
      endif else begin
        newpolys=refinedpolys
      endelse
   
    endif else begin
      if n_elements(newpolys) gt 1 then begin
        newpolys=[newpolys,polys[i:i+polys[i]]]
      endif else begin
        newpolys=polys[i:i+polys[i]]
      endelse
    endelse

    i=i+polys[i]
  end
 
  if has_normals then begin
    ptr_free,surf.normals
    ptr_free,surf.vertices
    ptr_free,surf.triangles
    surf=CREATE_STRUCT(['vertices','triangles','normals'],ptr_new(newvertices),ptr_new(newpolys),ptr_new(newnormals))
  endif else begin
    ptr_free,surf.vertices
    ptr_free,surf.triangles
    surf=CREATE_STRUCT(['vertices','triangles','normals'],ptr_new(newvertices),ptr_new(newpolys),ptr_new())
  endelse

end

