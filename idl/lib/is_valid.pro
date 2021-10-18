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
  function is_valid,  data, type, filename, single=single

    common pc_precision, zero, one, precision, data_type, data_bytes, type_idl

    if is_defined(data) then begin
      if size(data,/tname) eq 'STRUCT' then begin
        objname=tag_names(data,/structure_name) & colpos=strpos(objname,':')
        col2pos = strpos(objname,':',/reverse_search) 
        lename = colpos eq col2pos ? strlen(objname)-colpos : col2pos-colpos-1
        if colpos gt 0 and col2pos gt 0 and lename gt 0 then begin
          if (strmid(objname,0,colpos) eq 'PC_'+strupcase(strtrim(type,2))) and $
             (strmid(objname,colpos+1,lename) eq strupcase(filename)) then begin
            if col2pos gt colpos and is_defined(single) then begin
              if is_defined(precision) then prec=precision else prec='N' 
              if (strmid(objname,col2pos+1) eq precision and not keyword_set(single)) or $
                 (strmid(objname,col2pos+1) eq 'S' and keyword_set(single)) then return, 1   ; object is valid
            endif else $
              return, 1
          endif
        endif
      endif
      undefine, data
    endif
    return, 0

  end
