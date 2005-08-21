;;;;;;;;;;;;;;;;;;;;;;
;;;   wdwset.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   06-Jun-2002
;;;
;;;  Description:
;;;   Like wset (switch to window with index WID), but creates the
;;;   window if it doesn't exist.
;;;   All extra parameters are handed on to WINDOW.

pro wdwset, wid, _EXTRA=extra

  okerrs = [-459, $ ; Window is closed and unavailable (Linux IDL5.5a)
            -480, $ ; Window is closed and unavailable (Linux IDL5.3)
            -386, $ ; Window is closed and unavailable (OSF1 IDL5.2)
            -531  $ ; Window is closed and unavailable (Linux IDL6.0)
           ]                    ; these error codes are OK

  ; Establish error handler
  catch, errstat
  if (errstat ne 0) then begin
    case (errstat) of
      -458: begin               ; can't happen
        message, "You can't use WDWSET for a PostScript device", /INFO
      end
      ;
      else: begin
        match = where(errstat eq okerrs)
        if (match[0] ge 0) then begin ; known error code
          message, 'Creating window with index ' + strtrim(wid,2), /INFO
        endif else begin
          message, 'Catching unexpected error', /INFO
          print, 'Error index: ', errstat
          print, 'Error message: ', !ERROR_STATE.MSG
          message, 'Nevertheless creating window with index ' $
              + strtrim(wid,2), /INFO
        endelse
        window, wid, _extra=EXTRA
      end
    endcase
    catch, /CANCEL              ; remove error handler
  endif

  if (!d.name eq 'X') then begin
    if (n_elements(wid) gt 0) then begin
      wset, wid
    endif else begin
      wset
    endelse
  endif else begin
    message, "makes only sense with device='X'", /INFO
  endelse

end
; End of file wdwset.pro
