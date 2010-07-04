;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   running_gdl.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   10-Jun-2005
;;;
;;;  Description:
;;;    Returns -1 if we are running GDL, 0 otherwise (running expensive
;;;    IDL).

function running_gdl
  ;; Old test (see Marc Schellens' comment in June 2005, regarding
  ;; https://sourceforge.net/tracker/?func=detail&atid=618686&aid=1216778&group_id=97659)
  if (tag_names(/STRUCT, {a:0}) eq '$truct') then return, -1

  ;; Newer test, works with GDL-0.9 (as of Feb-2006):
  ;; Test for existence of !gdl system variable
  defsysv, '!GDL', EXISTS=exists
  return, -exists

end
; End of file running_gdl.pro
