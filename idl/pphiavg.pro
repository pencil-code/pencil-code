;;;;;;;;;;;;;;;;;;;;;;;
;;;   pphiavg.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   07-Apr-2004
;;;
;;;  Description:
;;;   Plot a number of phi averages for debugging purposes

default, datatopdir, 'data'

default, phiavgfile, 'PHIAVG1'

avg = read_phiavg(datatopdir+'/averages/'+phiavgfile)

save_state

!p.multi = [0,5,3]
!p.charsize = 1.8

for i=0,n_elements(avg.labels)-1 do begin

  label = strtrim(avg.labels[i],2)
  cmd = 'var=avg.'+label
  if (execute(cmd) ne 1) then $
      message, 'There was a problem with index.pro', /INFO

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
