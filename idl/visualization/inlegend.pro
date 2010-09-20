
pro inlegend, linestyles, symbols, texts, xrange, ystart, offset, shift=shift, sep=sep, xlog=xlog, ylog=ylog, circsymsiz=circsymsiz, _extra=extra, color=color

size = n_elements(linestyles)<n_elements(texts)

if n_elements(symbols) ne 0 then begin
  size = size < n_elements(symbols)
  withsymbs=1
  symbols = abs(symbols)
endif else $
  withsymbs=0
 
if keyword_set(circsymsiz) then $
  circ = n_elements(circsymsiz) $
else $
  circ = 0

if not keyword_set(sep) then $
  sep = 0 $
else $
  sep = abs(sep)

if not keyword_set(shift) then $
  shift = 0 

if keyword_set(xlog) then begin 

  if sep ne 0 then $
    sep=sign(sep)*alog10(abs(sep))

  xcenter = 10^(0.5*(alog10(xrange(0))+alog10(xrange(1))))

endif else begin
  xcenter = 0.5*(xrange(0)+xrange(1))
  xlog = 0
endelse

if keyword_set(ylog) then begin 

  offset=sign(offset)*alog10(abs(offset))
  if shift ne 0 then $
    shift=sign(shift)*alog10(abs(shift))

endif else $
  ylog=0

xcenter = [xcenter,xcenter]
yrange=[ystart,ystart]

;stop

for i=0,size-1 do begin

  oplot, xrange, yrange, lines=linestyles(i), color=color
 
  if withsymbs then begin
 
    if symbols(i) eq 8 and i lt circ then $
      if circsymsiz(i) ne 0 then $
        circ_sym, abs(circsymsiz(i)), fill=circsymsiz(i) lt 0

    if symbols(i) ne 0 then $
      oplot, xcenter, yrange, ps=-symbols(i), lines=linestyles(i), color=color

  endif
  if xlog then $
    xpos=10^(alog10(xrange(1))+sep) $
  else $
    xpos=xrange(1)+sep

  if ylog then $
    ypos=10^(alog10(yrange)+shift) $
  else $
    ypos=yrange+shift

  xyouts, xpos, ypos, texts(i), _extra=extra, color=color

  if ylog then $
    yrange=10^(alog10(yrange)+offset) $
  else $
    yrange = yrange+offset

endfor

end
