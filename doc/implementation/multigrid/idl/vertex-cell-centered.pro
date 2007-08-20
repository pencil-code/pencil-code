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

plot, [0,Lx], [-dx/2,Lx], $
      XSTYLE=4, YSTYLE=4, /NODATA, $
      XMARGIN=[4,2], YMARGIN=[2,2]

y = 0
;
for i=0,4 do begin
  x1 = linspace(0., 1., N1)
  plots, x1, 0*x1+y, PS=-4
  N1 = N1/2+1
  y = y + dy
endfor
xyouts, 0.2, y+dy/2, 'Vertex-centered grid'

y = 0
;
for i=0,3 do begin
  x2 = linspace(0.,1., N2)
  plots, Lx-x2, 0*x2+y, PS=-4
  N2 = N2/2
  y = y + dy
endfor
xyouts, 1.4, y+dy*3/2, 'Cell-centered grid'



end
; End of file vertex-cell-centered.pro
