;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vertex-cell-centered.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   17-Aug-2007
;;;
;;;  Description:
;;;    Sketch of vertex-centered (Pencil) vs cell-centered (textbooks)
;;;    multigrids.

N2 = 16
N1 = 16+1

Lx = 2.2
dx = 1./N2
dy = 0.3
ytop = 1.3

;scenario = 'simple'
scenario = 'smart'

plot, [0,Lx], [-dx/2,Lx], $
      XSTYLE=4, YSTYLE=4, /NODATA, $
      XMARGIN=[4,2], YMARGIN=[2,2]

y = 0
;
if (scenario eq 'simple') then begin
  for i=0,4 do begin
    x1 = linspace(0., 1., N1)
    plots, x1, 0*x1+y, PS=-4
    N1 = N1/2+1
    y = y + dy
  endfor
  xyouts, 0.2, y+dy/2, 'Vertex-centered grid'
endif else if (scenario eq 'smart') then begin
  N1 = N1-1
  for i=0,4 do begin
    x1 = linspace(0., 1., N1, ghost=-0.5)
    plots, x1, 0*x1+y, PS=-4
    plots, [0,1], 0*x1+y     ; make sure line is complete
    N1 = N1/2
    y = y + dy
  endfor
  xyouts, 0.2, y+dy/2, 'Modified "vertex-centered" grid'
endif
plots, [0, 0], [0, 1]*ytop
plots, [1, 1], [0, 1]*ytop

y = 0
;
if (scenario eq 'simple') then begin
  for i=0,3 do begin
    x2 = linspace(0.,1., N2)
    plots, Lx-1+x2, 0*x2+y, PS=-4
    N2 = N2/2
    y = y + dy
  endfor
  xyouts, 1.4, y+dy*3/2, 'Cell-centered grid'
endif else if (scenario eq 'smart') then begin
  for i=0,3 do begin
    x2 = linspace(0.,1., N2)
    plots, Lx-1+x2, 0*x2+y, PS=-4
    plots, [Lx-1,Lx], 0*x2+y     ; make sure line is complete
    N2 = N2/2 + 0.5
    y = y + dy
  endfor
  xyouts, 1.4, y+dy*3/2, 'Modified cell-centered grid'
endif else begin
  message, 'Scenario ' + scenario + ' unknown'
endelse

plots, Lx-[0, 0], [0, 1]*ytop
plots, Lx-[1, 1], [0, 1]*ytop


end
; End of file vertex-cell-centered.pro
