; $Id$
;
;  Concatenate an array of strings into a single string 
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2005-10-13 10:30:29 $
;  $Revision: 1.1 $
;
;  03-jun-02/tony: coded 
;
;
function arraytostring,array,QUOTE=QUOTE,LIST=LIST,NOLEADER=NOLEADER
COMPILE_OPT IDL2,HIDDEN


  default, QUOTE, ""
  default, LIST, ","
  result=""
  count=n_elements(array)

  for i=0,count-1 do begin
    if ((i eq 0) and (keyword_set(NOLEADER))) then begin 
      result=result+QUOTE+array[i]+QUOTE
    endif else begin
      result=result+LIST+QUOTE+array[i]+QUOTE
    endelse
  endfor

  return,result
end


