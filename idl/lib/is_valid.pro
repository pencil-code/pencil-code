;+
; Checks validity of object "data" against "type" and "filename":
; Object has to be defined and of type "STRUCT", further its structure_name must be of the form
; "PC_<type>:<filename> to be valid.
; If object is defined, but invalid, it is undefined.
; Relevant vor data read by pc_read_[pqb]dim, pc_read_dim, pc_read_grid, pc_read_param (so far).
;-
;
; 01-09-21/MR: coded
;
  function is_valid,  data, type, filename

    if is_defined(data) then begin
      if size(data,/tname) eq 'STRUCT' then begin
        objname=tag_names(data,/structure_name) & colpos=strpos(objname,':')
        if colpos gt 0 then $
          if (strmid(objname,0,colpos) eq 'PC_'+strupcase(strtrim(type,2))) and $
             (strmid(objname,colpos+1) eq strupcase(filename)) then return, 1   ; object is valid
      endif
      undefine, data
    endif
    return, 0

  end
