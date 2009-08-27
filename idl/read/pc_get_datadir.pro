;
; $Id$
;
; pc_get_datadir - gets the name of the data directory from datadir .in
;                  or sets the default value './data'
; 
; Date: 03.08.2007 (Anders & Chao-Chin)
;
function pc_get_datadir
;
; Default value of datadir.
;
datadir='./data'
;
; Open './datadir.in' for reading.
;
close, 1
openr, 1, './datadir.in', error=ierr
;
; If './datadir.in', go on and read contents.
;
if (ierr eq 0) then begin
  readf, 1, datadir
endif
;
close, 1
;
return, datadir
;
end
