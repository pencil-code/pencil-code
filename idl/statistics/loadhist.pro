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


cmd_load = "grep -E 'run.x|lam-flatd1' data/proc{[0-9],[0-9][0-9],[0-9][0-9][0-9]}/top.log | sed 's/^[^:]*://' | awk '{print $9}'"
cmd_mem  = "grep -E 'run.x|lam-flatd1' data/proc{[0-9],[0-9][0-9],[0-9][0-9][0-9]}/top.log | sed 's/^[^:]*://' | awk '{print $10}'"

spawn, cmd_load, loadstr
load = float(loadstr)
;
spawn, cmd_mem, memstr
mem = float(memstr)

if (n_elements(load)*n_elements(mem) lt 4) then begin
  print, 'Not enough data'
endif else begin

  save_state
  !p.multi = [0,2,2]

  ;; Plot as function of CPU #
  !x.title = '!6CPU #!X'
  !y.title = '!6%!X'
  !x.style = 1
  !y.style = 1
  !y.range = [0,105]
  !p.psym = 10
  ;
  plot, load, TITLE='!6CPU Load!X'
  ophline,100
  ophline,mean(load),LINE=2
  ;
  plot, mem, TITLE='!6Memory percentage!X'
  ophline,[100]
  ophline,mean(mem),LINE=2


  ;; Plot histograms
  !x.title = '!6%!X'
  !y.title = '!6Absolute frequency!X'
  !y.range = 0
  ;
  histo = histogram(load,MIN=0.,MAX=100.,NBIN=25)
  N_hist = n_elements(histo)-1
  perc = indgen(N_hist)*100./(N_hist-1)
  plot, perc, [histo,0], $
      TITLE='!6CPU Load!X', YRANGE=[0,1.05*max(histo)]
  opvline,mean(load),LINE=2

  ;
  histo = histogram(mem,MIN=0.,MAX=100.,NBIN=25)
  N_hist = n_elements(histo)-1
  perc = indgen(N_hist)*100./(N_hist-1)
  plot, perc, histo, $
      TITLE='!6Memory percentage!X', YRANGE=[0,1.05*max(histo)]
  opvline,mean(mem),LINE=2

  restore_state

endelse


end
; End of file loadhist.pro
