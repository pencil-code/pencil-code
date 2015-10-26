; $Id$
;
;  Trim f-array variables in structure of ghost zones and empty dimensions.
;  Copies multi-dimensional elements of the structure into the IDL environment.
;
;  Author: Anders Johansen (ajohan@astro.ku.dk)
;  $Date: 2004-05-25 09:19:32 $
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
; Avoid non-grid variables
;
    string = 'size_var = size(' + object + '.' + tags[i] + ')'
    res = execute(string)

    if (size_var[0] ne 0) then begin
;
; x     -> xxx
; y     -> yyy
; z     -> zzz
; uu    -> uuu
; rho   -> rrho
; lnrho -> llnrho
; TT    -> TTT
; lnTT  -> llnTT
;      
      if ((tags[i] eq 'X') or (tags[i] eq 'Y') or (tags[i] eq 'Z')) then $
          tag_extra = tags[i]
      string = strmid(tags[i],0,1) + tags[i] + tag_extra + ' = pc_trim_var(' + $
          object + '.' + tags[i] + ')'
      print, string
      res = execute(string)
    endif

    tag_extra = ''
    
  endfor 

end
