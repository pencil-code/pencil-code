;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   read-ser-par.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   24-Jun-2003
;;;  $Id$
;;;
;;;  Description:
;;;   Read data from a run that was done in serial and parallel and
;;;   write them )in identical format) to vars-ser.sav and vars-par.sav.
;;;  Assumes:
;;;   A run with lnrho, uu, ss (possibly aa).
;;;  Usage:
;;;   Source with
;;;     @run-ser-par

spawn, 'ls -d data-ser data-par', check
if (n_elements(check) ne 2) then STOP, "Can't find data-{ser,par} directories"

default, quiet, 0
oldquiet = { q1: quiet, q2: !quiet }
quiet = 5
!quiet = 1

par2 = {nt: 0, t: 0}

;; Parallel run
datatopdir = 'data-par'
.r start
.r rall
;; Write to save file, so we can postprocess
lmn = {l1: l1, l2: l2, m1: m1, m2: m2, n1: n1, n2: n2}
if (lmagnetic) then $
  varp = { nt: par2.nt, t: t, uu: uu, lnrho: lnrho, ss: ss, aa: aa, lmn: lmn } $
else $
  varp = { nt: par2.nt, t: t, uu: uu, lnrho: lnrho, ss: ss, lmn: lmn }
save, FILE='vars-par.sav', varp

;; Serial run
datatopdir = 'data-ser'
.r start
.r rall
;; Write to save file, so we can postprocess
lmn = {l1: l1, l2: l2, m1: m1, m2: m2, n1: n1, n2: n2}
if (lmagnetic) then $
  vars = { nt: par2.nt, t: t, uu: uu, lnrho: lnrho, ss: ss, aa: aa, lmn: lmn } $
else $
  vars = { nt: par2.nt, t: t, uu: uu, lnrho: lnrho, ss: ss, lmn: lmn }
save, FILE='vars-ser.sav', vars

;; Calculate some differences
duu = vars.uu   -varp.uu
dll = vars.lnrho-varp.lnrho
dss = vars.ss   -varp.ss
if (lmagnetic) then daa = vars.aa - varp.aa

quiet = oldquiet.q1
!quiet = oldquiet.q2

;; Summarize differences found:
if (vars.t ne varp.t) then $
    print, 'WARNING: times differ by ', vars.t-varp.t
fmt = '(A4,A,5G14.5)'
print, '                        inner                   |          full'
print, ' var       min           max           rms      |    min           max'
; duu
for i=0,2 do begin $
    print, FORMAT=fmt, 'duu_', xyz[i], $
        minmax(duu[l1:l2,m1:m2,n1:n2,i]), rms(duu[l1:l2,m1:m2,n1:n2,i]), $
        minmax(duu[*,*,*,i])
; dll
    print, FORMAT=fmt, 'dll ', ' ', $
        minmax(dll[l1:l2,m1:m2,n1:n2]), rms(dll[l1:l2,m1:m2,n1:n2]), $
        minmax(dll)
; dss
    print, FORMAT=fmt, 'dss ', ' ', $
        minmax(dss[l1:l2,m1:m2,n1:n2] ), rms(dss[l1:l2,m1:m2,n1:n2]), $
        minmax(dss)
; daa
if (lmagnetic) then $
    for i=0,2 do begin $
        print, FORMAT=fmt, 'daa_', xyz[i], $
            minmax(daa[l1:l2,m1:m2,n1:n2,i]), rms(daa[l1:l2,m1:m2,n1:n2,i]), $
            minmax(daa[*,*,*,i])

; End of file read-ser-par.pro
