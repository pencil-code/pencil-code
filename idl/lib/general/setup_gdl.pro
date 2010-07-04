;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   setup_gdl.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   4-Jul-2010
;;;
;;;  Description:
;;;    Set up GDL, so we can use basic Pencil Code functionality with
;;;    gnudatalanguage version ~rc1-1.1ubuntu2 .
;;;    Currently, we just add the examples/pro directory to IDL_PATH
;;;  Usage:
;;;    .run setup_gdl
;;;  Limitations:
;;;    This was only ever tested on Ubuntu Jaunty.

@running_gdl ;; For some reason, this seems to be necessary with pc_read_var

dir = '/usr/share/doc/gnudatalanguage/examples/pro'
if (file_test(dir, /DIRECTORY)) then begin
  !path = !path + ':' + dir
endif

end
