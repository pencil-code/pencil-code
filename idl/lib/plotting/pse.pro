;;;;;;;;;;;;;;;;;
;;;  pse.pro  ;;;
;;;;;;;;;;;;;;;;;

;;; Author:  wd (Wolfgang.Dobler@Newcastle.ac.uk)
;;; Date:    9-May-2000
;;; Version: 1.1

;;;   Description: Finish output to PostScript device after psa, psl or epsf.
;;;     Switch back to device used before.
;;;     Insert sequence of IDL commands into PostScript file for
;;;     reproducing the plot.
;;;   Usage:
;;;     IDL> epsf    or     psa    or    psl
;;;     IDL> plot, ...
;;;     IDL> pse
;;;   Key words:
;;;     FIXBB  - run `psfixbb' on the postscript file to get
;;;              shrink-wrapped BoundingBox
;;;   Problems:
;;;     Extracting the history only works if psa/psl/epsf was called
;;;     interactively at the prompt -- the content of a sourced file
;;;     never  appears in the history.

pro pse, $
         FIXBB=fixbb, DEBUG=debug

  COMMON _ps_fonts, _old_P_FONT
  COMMON _ps1,_oldthick,_fname,_olddev
  ON_ERROR, 2

  default, fixbb, 0
  default, debug, 0

  if (!d.name eq 'PS') then begin
    ;; Reset PostScript device to default settings before closing
    device, XOFFSET=1.87, XSIZE=17.15, YOFFSET=12.70, YSIZE=11.99, $
        SCALE_FACTOR=1.
    device, /CLOSE
  endif
  dev = _olddev
  if (dev eq '') then dev=strtrim(GETENV("IDL_DEVICE"))
  if (dev eq '') then dev = 'X' 
  set_plot, dev
;; Reset !P.FONT
  IF (n_elements(_old_P_FONT) gt 1) THEN BEGIN
    !P.FONT = _old_P_FONT[1]
    _old_P_FONT = 0
  ENDIF
;; Reset drawing thicknesses
  IF (n_elements(_oldthick) gt 1) THEN BEGIN
    !P.CHARTHICK = _oldthick[1]
    !P.THICK = _oldthick[2]
    !X.THICK = _oldthick[3]
    !Y.THICK = _oldthick[4]
    !Z.THICK = _oldthick[5]
    _oldthick = 0
  ENDIF

  ;;;
  ;;; Insert some additional info into the just written (E)PS file:
  ;;;
  help, OUTPUT=hist, /RECALL_COMMANDS ; get history list
  ; drop uninteresting first line `Recall buffer length'; reverse to
  ; chronological order
  if (n_elements(hist) ge 2) then hist = hist[1:*]
  nhist = n_elements(hist)
  psa_class = ['psa','psl','epsf']
  cmdlist = '' & i=0 & last = 0 & no_psa = 0
  while (not last) do begin
    cmd = hist[i]
    ;; drop leading enumeration label
    ws = strpos(cmd,'	')        ; find first tab character
    cmd = strtrim(strmid(cmd,ws+1),2)
    cmdlist = cmd + '\n' + cmdlist
    ;; check for psa, psl and epsf and stop us there
    cmdup = strlowcase(cmd)
    for j=0,n_elements(psa_class)-1 do begin
      if (strpos(cmd,psa_class[j]) eq 0) then begin
        ; check for end-of-word (tedious whithout character classes)
        nextchar = strmid(cmd,strlen(psa_class[j]),1)
        if (nextchar eq '' or nextchar eq ',' or nextchar eq '&' $
	    or nextchar eq ';' or nextchar eq ' ') then last=1
      endif
    endfor
    i = i+1
    if ((not last) and (i ge nhist)) then begin
      last = 1
      no_psa = 1
      message, $
          "Couldn't find matching psa, psl or epsf in history", /INFORMATIONAL
    endif
  endwhile

  ;; document problems
  if (no_psa) then begin
    comment = " 'No matching psa, psl or epsf found -- showing last Cmd only'"
    cmd = hist[0]
    ws = strpos(cmd,'	')      ; find first tab character
    cmd = strtrim(strmid(cmd,ws+1),2)
    cmdlist = cmd
  endif else begin
    comment = ''
  endelse

  ;; quote double quote characters
  pos = 0
  iquot = strpos(cmdlist,'"',pos)
  if (debug) then print,'1.  pos, iquot = ', pos, iquot
  while (iquot ge pos) do begin
    cmdlist = strmid(cmdlist,0,iquot)+'\'+strmid(cmdlist,iquot)
    pos = iquot+2
    iquot = strpos(cmdlist,'"',pos)
    if (debug) then print, 'Loop: cmdlist=<', cmdlist,'>'
    if (debug) then print, 'Loop:   pos, iquot = ', pos, iquot
  endwhile

  ;; remove trailing \n
  l = strlen(cmdlist)
  if (strmid(cmdlist,l-2) eq '\n') then cmdlist = strmid(cmdlist,0,l-2)

  ;; modify PS file
  if (debug) then begin
    print, $
        '<' + 'ps-annotate' + ' --cmd "' + cmdlist + '" ' + _fname+comment + '>'
  endif
  spawn, /SH, 'ps-annotate' + ' --cmd "' + cmdlist + '" ' + _fname + comment
  if (fixbb) then begin
    if (debug) then begin
      print, 'Running psfixbb..'
      spawn, 'psfixbb ' + _fname
    endif else begin
      spawn, 'psfixbb -x 2 ' + _fname, dev_null, /stderr ; intercept output
    endelse
  endif

end
