;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   hsections.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Nov-2001
;;;
;;;  Description:
;;;   Plot velocity, density and entropy field in three horizontal
;;;   sections.

;;; Unfinished..

default, absolute, 0            ; flag four absolute colour scaling (i.e.
                                ; relative to absolute min and max of
                                ;  colour-represented data

nrholevs = 15                   ; No of isolines
nuulevs = 60                    ; No of colours
nentlevs = 60                   ; No of colours

nz1 = 0.2*nz > 4
nz2 = 0.5*nz
nz3 = 0.8*nz < nz-5

sz1 = '!8z!6='+strtrim(z[nz1],2)
sz2 = '!8z!6='+strtrim(z[nz2],2)
sz3 = '!8z=!6'+strtrim(z[nz3],2)

save_state

wput

!p.multi = [0,3,2]
!p.charsize = 1.8
!x.title = '!8x!X'
!y.title = '!8y!X'

tit = '!17u!6 at '

if (absolute) then begin
  zruu = minmax(uu[*,*,*,2])
endif else begin
  undefine, zruu                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

plot_3d_vect, uu[*,*,nz1,*],x,y, /KEEP, TITLE=tit+sz1+'!X', ZRANGE=zruu
opcircles, 1., LINE=2, THICK=2
plot_3d_vect, uu[*,*,nz2,*],x,y, /KEEP, TITLE=tit+sz2+'!X', ZRANGE=zruu
opcircles, 1., LINE=2, THICK=2
plot_3d_vect, uu[*,*,nz3,*],x,y, /KEEP, TITLE=tit+sz3+'!X', ZRANGE=zruu
opcircles, 1., LINE=2, THICK=2

tit = '!8s!6 and !7r!6 at '
!x.style = 1
!y.style = 1
!x.range = [x[3], x[nx-4]]      ; No ghost zones
!y.range = [y[3], y[ny-4]]

;
if (absolute) then begin
  levent = linspace(minmax(ent),nentlevs)
endif else begin
  undefine, levent                 ; LEVELS=<undef> is like no LEVELS kw at all
endelse

contourfill, ent[*,*,nz1],x,y, TITLE=tit+sz1+'!X', LEVELS=levent
var = reform(lam[*,*,nz1])
contour, var,x,y, /OVER, LEV=linspace(minmax(var),nrholevs)
opcircles, 1., LINE=2, THICK=2
;
contourfill, ent[*,*,nz2],x,y, TITLE=tit+sz2+'!X', LEVELS=levent
var = reform(lam[*,*,nz2])
contour, var,x,y, /OVER, LEV=linspace(minmax(var),nrholevs)
opcircles, 1., LINE=2, THICK=2
;
contourfill, ent[*,*,nz3],x,y, TITLE=tit+sz3+'!X', LEVELS=levent
var = reform(lam[*,*,nz3])
contour, var,x,y, /OVER, LEV=linspace(minmax(var),nrholevs)
opcircles, 1., LINE=2, THICK=2

wget

restore_state

end
; End of file hsections.pro
