pro free_surfaces, data
  if (size (data, /type) eq 0) then return

  if (size (data, /type) eq 8) then begin
    if (has_tag (data, 'vertices')) then $
        if (ptr_valid (data.vertices)) then ptr_free, data.vertices
    if (has_tag (data, 'tangents')) then $
        if (ptr_valid (data.tangents)) then ptr_free, data.tangents
    if (has_tag (data, 'normals')) then $
        if (ptr_valid (data.normals)) then ptr_free, data.normals
  endif else begin
    for i=0, n_elements(data)-1 do begin
      free_surfaces, data[i]
    endfor
    data=0.
  endelse
end

