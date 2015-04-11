function rtable,file,ncolumn,form=form,head=head,lines=lines,debug=debug,lun=lun,noclose=noclose
;
; IN: noclose - file is not closed
; IN: lines - number of lines to be read - 1 !!!
;
if n_params(0) eq 0 then begin
  print,'function table,file,ncolumn,form=form,head=head,lines=lines,debug=debug,lun=lun,noclose=noclose'
  return,0
end
;
;  read header
; 
if not keyword_set(lun) then lun=-1
if lun lt 0 then $
  
openr,lun,/get_lun,file

if (n_elements(head) ne 0) then begin
  card=''
  for i=0,abs(head)-1 do begin
    readf,lun,card
    if (head lt 0) then print,card
  endfor
endif
;
;  read: free format
;
a=dblarr(ncolumn) & aa=a & n=0L
if n_elements(form) eq 0 then begin
  readf,lun,a
  n=n+1L
  if n_elements(lines) eq 0 then begin
    while not eof(lun) do begin
      readf,lun,aa
      if (keyword_set(debug)) then print,n,aa
      a=[a,aa]
      n=n+1L
    endwhile
  endif else begin
    for ii=0,lines-1 do begin
      readf,lun,aa
      a=[a,aa]
      n=n+1L
    endfor
  endelse
endif
;
;  read: format given
;
if n_elements(form) ne 0 then begin
  readf,lun,a,form=form
  n=n+1L
  if n_elements(lines) eq 0 then begin
    while not eof(lun) do begin
      readf,lun,aa,form=form
      a=[a,aa]
      n=n+1L
    endwhile
  endif else begin
    for ii=0,lines-1 do begin
      readf,lun,aa,form=form
      a=[a,aa]
      n=n+1L
    endfor
  endelse
endif
;
if not keyword_set(noclose) then begin
  close,lun
  free_lun,lun
endif

return,reform(a,ncolumn,n)
end
