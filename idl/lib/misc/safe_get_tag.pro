; $Id$
;
;  Safely attempt to access a particular tag of a structure by ensuring 
;  that the tag exists first. 
;  
;  Returns !VALUES.F_NAN if the tag is unspecified, or the value of the element
;  if the tag is defined. 
;
;  If DEFAULT is set to some value then this value will be returned instead.
; 
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2005-10-20 09:18:32 $
;  $Revision: 1.2 $
;
;  03-jun-02/tony: coded 
;
;
function safe_get_tag,object,tag,DEFAULT=DEFAULT
COMPILE_OPT IDL2,HIDDEN

  found = where(tag_names(object) eq strupcase(tag))

  if (found[0] eq -1) then begin
    if (n_elements(DEFAULT) ne 0) then return, DEFAULT
    return, !VALUES.F_NAN
  endif

  return, object.(found[0])
end


