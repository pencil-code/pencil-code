;;;;;;;;;;;;;;;;;;;
;;;   ts.pro   ;;;
;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  $Date: 2007-08-31 01:44:25 $
;;;  $Revision: 1.10 $
;;;  Description:
;;;   Read time series data from data/time_series.dat into the
;;;   structure `ts' and plot urms(t) and brms(t) (if available).
;;;  Usage:
;;;   .r ts
;;;   plot, ts.t, ts.brms, YRANGE=minmax([ts.brms,ts.bmax]), /YLOG
;;;   oplot, ts.t, ts.bmax, LINE=2

function parse_tsheader, hline
;
;  Split header line into the individual variable names.
;
  line  = strmid(hline,strpos(hline,'#')+1)
  line2 = strmid(hline,strpos(hline,'%')+1)
  if strlen(line2) lt strlen(line) then line=line2 
  labels = ['']
  ;
  ; strsplit() is not available in IDL prior to 5.3 (and str_sep() was
  ; obsoleted after 5.2..), so we do this manually:
  ;
  while (line ne '') do begin
    repeat begin
      line = strmid(line,1)     ; remove first character
    endrep until (strmid(line,0,1) ne '-')
    endlb = strpos(line,'-')
    if (endlb lt 0) then endlb=strlen(line) ; if no '-' at end of line
    labels = [labels, strmid(line,0,endlb)]
    line = strmid(line,endlb)
  endwhile
  ;
  ;  eliminate empty labels
  ;
  good = where(labels ne '')
  return, labels[good]
end
; ---------------------------------------------------------------------- ;
function list_idx, label, list
;
;  Return index if label is contained in list, else -1
;
  return, (where(list eq label))[0]
end
; ---------------------------------------------------------------------- ;
function in_list, label, list
;
;  Return 1 if label is contained in list, else 0
;
  return, (list_idx(label,list)+1) ne 0
end
; ---------------------------------------------------------------------- ;

;
;  Set path and file names if necessary
;
datatopdir = pc_get_datadir(datatopdir)
cd, cur=wd
if strtrim(wd,2)+'/data/' ne strtrim(datatopdir,2) then datatopdir=strtrim(wd,2)+'/data/'
tsfile=strtrim(datatopdir,2)+'/time_series.dat'
;
;  read header
;
close, 1                         ; just in case
openr, 1, tsfile
line = ''
repeat begin
  readf, 1, line
  line = strtrim(line)
  hashpos  = strpos(line,'#')
  hashpos2 = strpos(line,'%')
  if hashpos2 gt hashpos then hashpos=hashpos2 
  ;; identify header line as / *#--/
endrep until ((hashpos ge 0) and (strmid(line,hashpos+1,2) eq '--'))
labels = parse_tsheader(line)
ncols = n_elements(labels)
;
;  read table
;
data = input_table(tsfile,/DOUBLE)
if ((size(data))[1] ne ncols) then begin
  message, /INFO, 'Inconsistency: label number different from column number'
endif
;
;  assemble the data
;
cmd = 'ts = {'           ; build a command to execute
for i=0,ncols-1 do begin
  if (i eq 0) then pref=' ' else pref=', '
  cmd = cmd + pref + labels[i] + ': reform(data[' + strtrim(i,2) + ',*])'
endfor
cmd = cmd + ' }'

if (execute(cmd) ne 1) then $
    message, 'There was a problem executing <' + cmd + '>', /INFO

;
;  Do two plots
;  Try to plot urms(t) and brms(t), if any of these two is not
;  available or zero, fill the list with the first two variables other
;  than `it' and `dt'
;
if (in_list('t',labels)) then begin
  idxlist = [-1]                ; list of indices
  save_state                    ; save !p.multi, etc.
  for i=0,ncols-1 do begin      ; collect all non-special variables
    if (not in_list(labels[i], ['it','t','dt','urms','brms'])) then $
        idxlist = [idxlist, list_idx(labels[i],labels)]
;    print, strtrim(i,2), ' ', labels[i], $
;        in_list(labels[i], ['it','t','dt','urms','brms'])
  endfor
;  idxlist = [list_idx('brms',labels),idxlist]
;  idxlist = [list_idx('urms',labels),idxlist]
  if (in_list('brms',labels)) then begin
    if (max(ts.brms gt 0)) then idxlist = [list_idx('brms',labels),idxlist]
  endif
  if (in_list('urms',labels)) then begin
    if (max(ts.urms gt 0)) then idxlist = [list_idx('urms',labels),idxlist]
  endif
  idxlist = idxlist[where(idxlist ge 0)] ; clean list
  nplots = min([n_elements(idxlist),2])
  if (nplots eq 1) then !p.multi=[0,1,1] else !p.multi=[0,1,2]
  !x.title='!8t!X'
  for i=0,nplots-1 do begin
    lab = labels[idxlist[i]]
    if ((lab eq 'brms') or (lab eq 'bmax') or (lab eq 'urms') or (lab eq 'umax')) then ylog=',/YLOG' else ylog=''
    yrange=',YRANGE=minmax(ts.'+lab+'[1:*])>1.e-4*max(ts.'+lab+'), YSTYLE=3'
    !y.title = '!3'+lab+'!X'
    pcmd = 'plot, ts.t, ts.'+lab+ylog+yrange
    if (execute(pcmd) ne 1) then $
        message, 'There was a problem executing <' + pcmd + '>', /INFO
  endfor
  restore_state
  print, "Type   help,ts,/STRUCT   for a full list of available slots"
endif else begin
  message, "Odd data: cannot find time `t' in time series"
endelse


close,1
end
; End of file ts.pro
