 function sizeoftag, struct,tag
;
;  Returns size of member tag of structure struct without making any copy
;          -2 if type of struct or tag incorrect
;          -1 if tag not member of struct
;          scalar of value 1 if tag is scalar
;          vector of dimensions if tag is array
;
   if not is_struct(struct) then return, -2
   if not is_str(tag) then return, -2

   help, struct, out=res
  
   inds=where(stregex(res,'[^0-9a-zA-Z_]'+strupcase(tag)+'[^0-9a-zA-Z_]') ge 0)

   if inds[0] ge 0 then begin
     substr=stregex(res[inds[0]],'^.*Array(.+) *$',/SUBEXPR,/EXTRACT)
     if substr[1] eq '' then $
       return, 1 $
     else begin
       ret=execute('sz='+substr[1])
       return, sz
     endelse
   endif else $
     return, -1  

 end
