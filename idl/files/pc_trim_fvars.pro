; $Id: pc_trim_fvars.pro,v 1.1 2004-05-11 09:14:48 ajohan Exp $
;
;  Trim f array variables in structure of ghost zones and empty dimensions
;
;  Author: Anders Johansen (ajohan@astro.ku.dk)
;  $Date: 2004-05-11 09:14:48 $
;  $Revision: 1.1 $
;
;  11-may-04/anders: coded
;
  default, object, 'mydata'
;
; Extract variable names from structure
;
  string = 'tags=tag_names(' + object + ')'
  res = execute(string)
  ntags = n_elements(tags)

  for i=0,ntags-1 do begin
;
; Avoid non-f variables
;    
    if (tags(i) ne 'T' and $
        tags(i) ne 'X' and tags(i) ne 'Y' and tags(i) ne 'Z' and $
        tags(i) ne 'DX' and tags(i) ne 'DY' and tags(i) ne 'DZ') then begin
      string = strmid(tags(i),0,1)+tags(i) + ' = pc_trim_var(' + $
          object + '.' + tags(i) + ')'
      print, string
      res = execute(string)
    endif
  endfor 

end
