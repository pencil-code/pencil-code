pro free_surfaces,data
  if n_elements(data) lt 1 then return

  if n_elements(data) eq 1 then begin
    if n_tags(data) eq 0 then return
    tags=TAG_NAMES(data)
    for i=0,n_tags(data)-1 do begin
      if (tags[i] eq 'vertices') then begin
        if ptr_valid(data.vertices) then ptr_free,data.vertices
      endif else if (tags[i] eq 'tangents') then begin
        if ptr_valid(data.tangents) then ptr_free,data.tangents
      endif else if (tags[i] eq 'normals') then begin
        if ptr_valid(data.normals) then ptr_free,data.normals
      endif
    endfor
  endif else begin
    for i=0,n_elements(data)-1 do begin
      free_surfaces,data[i] 
    endfor
    data=0.
  endelse
end

