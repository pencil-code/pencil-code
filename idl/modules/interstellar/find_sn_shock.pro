PRO find_SN_shock,var,x,pnts=pnts,leftpnt=leftpnt,rightpnt=rightpnt,tol=tol
;
; $Id$
;
; Routine to find shock front by finding points within a given
; tolerance of the maximum
;
;
default,tol,0.000000001

pnts = x[where(abs(var - max(var)) le tol)]
rightpnts = where(pnts gt 0.,nrightpnts)
leftpnts = where(pnts lt 0.,nleftpnts)
print,abs(var - max(var))
if (nrightpnts gt 1 ) then begin
  rightpnt = mean(pnts[rightpnts])
endif else begin
  if rightpnts eq -1 then rightpnt=-9999999. else rightpnt=pnts[rightpnts]
endelse


if (nleftpnts gt 1) then begin
  leftpnt = mean(pnts[leftpnts])
endif else begin
  if leftpnts eq -1 then leftpnt=-9999999. else leftpnt=pnts[leftpnts]
endelse

;leftpnts=where((var-var[0]) gt 10. and (var - var[0]) lt tol)
;rightpnts=where((var-var[n_elements(var)-1]) gt 10. and (var - var[n_elements(var)-1]) lt tol)
;
;if n_elements(rightpnts) gt 1 then begin
;  rightpnt = max(x[rightpnts])
;endif else begin
;  if rightpnts eq -1 then rightpnt=-9999999. else rightpnt=x[rightpnts]
;endelse
;
;
;if n_elements(leftpnts) gt 1 then begin
;  leftpnt = min(x[leftpnts])
;endif else begin
;  if leftpnts eq -1 then leftpnt=-9999999. else leftpnt=x[leftpnts]
;endelse


;;leftpnt = mean(pnts[where(pnts lt 0.)])

;print,'found points:',pnts
;print,'positive x shock at',mean(rightpnts)
;print,'negative x shock at',mean(leftpnts)

END
