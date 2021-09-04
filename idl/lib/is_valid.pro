  function is_valid,  data, type
    type = strupcase(strtrim(type,2))
    if size(data,/type) ne 8 then $
      return, 0 $
    else begin
      name=strupcase(tag_names(data,/structure_name))
      if strmatch(name,'PC_'+type+'*') then $
        return, (getenv('PC_VALID_'+type) eq 'V') $
      else $
        return, 0
    endelse
  end
