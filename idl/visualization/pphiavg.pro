;;;;;;;;;;;;;;;;;;;;;;;
;;;   pphiavg.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   07-Apr-2004
;;;
;;;  Description:
;;;    Plot all known phi averages
;;;  Usage:
;;;    pphiavg                            ; read data from default location
;;;    pphiavg, avg                       ; use data from struct AVG 
;;;    pphiavg, 'data/averages/PHIAVG5'   ; read data/averages/PHIAVG5
;;;    pphiavg, 'PHIAVG5'                 ; read data/averages/PHIAVG5
;;;    pphiavg, PMULTI=[0,3,4], CHARSIZE=1.4

pro pphiavg, arg, QUIET=quiet, PMULTI=pmulti, CHARSIZE=charsize, $
             HELP=help

  if (keyword_set(help)) then extract_help, 'pphiavg'

  datatopdir ='data'
  avgdir = datatopdir+'/averages'
  phiavgfile = 'PHIAVG1'

  default, arg, avgdir+'/'+phiavgfile
  default, quiet, 0
  default, pmulti, [0,5,3]
  ; Crude heuristics for acceptable charsize -- could be improved:
  nmulti = ((max(pmulti[1:2])>2)-1)/2
  default, charsize, 1. + nmulti/(2.+nmulti)*2.

  s = size(arg)
  if (s[s[0]+1] eq 8) then begin ; ARG is a struct
    avg = arg
  endif else if(s[s[0]+1] eq 7) then begin ; ARG is a string
    ;; Try ARG, data/ARG, data/averages/ARG in that order:
    fnames = ['', datatopdir+'/', avgdir+'/']+arg
    found = 0
    filelist  =''
    for i=0,n_elements(fnames)-1 do begin
      if (not found) then begin
        file = fnames[i]
        filelist = filelist + file + ', '
        if (any(file_test(file))) then found = 1
      endif
    endfor
    if (not found) then message, "Couldn't open any file of "+filelist

    ;; Get data from file:
    if (not quiet) then print, 'Reading '+file

    avg = pc_read_phiavg(file)

    ;
    ; increase !p.multi if necessary
    ;
    if (avg.nvars gt pmulti[1]*pmulti[2]) then begin
      npy = sqrt(avg.nvars/1.6)
      pmulti = floor([0,1.6*npy+1,npy+1])
    endif

  endif else begin
    message, 'ARG must be a struct or a string'
  endelse

  save_state

  !p.multi = pmulti
  !p.charsize = charsize

  for i=0,n_elements(avg.labels)-1 do begin

    label = strtrim(avg.labels[i],2)
    cmd = 'var=avg.'+label
    if (execute(cmd) ne 1) then $
        message, /INFO, 'There was a problem executing "'+cmd+'"'
    pos = aspect_pos(2.,MARGIN=0.08)
    !x.style = 1.
    !y.style = 1.
    !p.title=label

    if (max(var)-min(var) gt 0) then begin
      contourfill, var, avg.rcyl, avg.z, POS=pos
    endif else begin
      contour, var, avg.rcyl, avg.z, POS=pos
    endelse

  endfor

  restore_state

end
; End of file pphiavg.pro
