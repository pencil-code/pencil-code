;;;;;;;;;;;;;;;;;;;;;;;;
;;;   loadhist.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   20-Aug-2004
;;;
;;;  Description:
;;;   Read load info from data/proc*/top.log and plot as function of node
;;;   number and as histogram.
;;;   Do same thing for memory percentage (which is probably not interesting).

save_state
!p.multi = [0,2,2]


cmd_load = "grep -E 'run.x|lam-flatd1' data/proc{[0-9],[0-9][0-9],[0-9][0-9][0-9]}/top.log | sed 's/^[^:]*://' | awk '{print $9}' | head -115"
cmd_mem  = "grep -E 'run.x|lam-flatd1' data/proc{[0-9],[0-9][0-9],[0-9][0-9][0-9]}/top.log | sed 's/^[^:]*://' | awk '{print $10}' | head -115"

spawn, cmd_load, loadstr
load = float(loadstr)
;
spawn, cmd_mem, memstr
mem = float(memstr)

;; Plot as function of CPU #
!x.title = '!6CPU #!X'
!y.title = '!6%!X'
!x.style = 1
!y.style = 1
!y.range = [0,105]
!p.psym = 10
;
plot, load, TITLE='!6CPU Load!X'
ophline,[100]
;
plot, mem, TITLE='!6Memory percentage!X'
ophline,[100]


;; Plot histograms
!x.title = '!6%!X'
!y.title = '!6Relative frequency!X'
;
histo = histogram(load,MIN=0.,MAX=100.,BINSIZE=4)
N_hist = n_elements(histo)
perc = indgen(N_hist)*100./(N_hist-1)
plot, perc, histo, TITLE='!6CPU Load!X'
;
histo = histogram(mem,MIN=0.,MAX=100.,BINSIZE=4)
N_hist = n_elements(histo)
perc = indgen(N_hist)*100./(N_hist-1)
plot, perc, histo, TITLE='!6Memory percentage!X'


restore_state


end
; End of file loadhist.pro
