;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vsections.pro   ;;;
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
                                ; colour-represented data
default, show_ghosts, 0

nrholevs = 15                   ; No of isolines
nuulevs = 60                    ; No of colours
nentlevs = 60                   ; No of colours

ny1 = 0.2*ny > 4
ny2 = 0.5*ny
ny3 = 0.8*ny < ny-5

sy1 = '!8y!6='+strtrim(y[ny1],2)
sy2 = '!8y!6='+strtrim(y[ny2],2)
sy3 = '!8y=!6'+strtrim(y[ny3],2)

save_state

wput

!p.multi = [0,3,2]
!p.charsize = 1.8
!x.title = '!8x!X'
!y.title = '!8z!X'

tit = '!17u!6 at '
!x.style = 1
!y.style = 1
if (show_ghosts) then begin
  !x.range = [x[0], x[nx-1]]    ; No ghost zones
  !y.range = [z[0], z[nz-1]]
endif else begin
  !x.range = [x[3], x[nx-4]]    ; No ghost zones
  !y.range = [z[3], z[nz-4]]
endelse

if (absolute) then begin
  zruu = minmax(uu[*,*,*,1])
endif else begin
  undefine, zruu                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

plot_3d_vect, uu[*,ny1,*,*],x,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sy1+'!X', ZRANGE=zruu
opcircles, 1., LINE=2, THICK=2
plot_3d_vect, uu[*,ny2,*,*],x,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sy2+'!X', ZRANGE=zruu
opcircles, 1., LINE=2, THICK=2
plot_3d_vect, uu[*,ny3,*,*],x,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sy3+'!X', ZRANGE=zruu
opcircles, 1., LINE=2, THICK=2

tit = '!8s!6 and !7r!6 at '

;
if (absolute) then begin
  levent = linspace(minmax(ent),nentlevs)
endif else begin
  undefine, levent                 ; LEVELS=<undef> is like no LEVELS kw at all
endelse

contourfill, ent[*,ny1,*],x,z, TITLE=tit+sy1+'!X', LEVELS=levent
var = reform(lam[*,ny1,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
opcircles, 1., LINE=2, THICK=2
;
contourfill, ent[*,ny2,*],x,z, TITLE=tit+sy2+'!X', LEVELS=levent
var = reform(lam[*,ny2,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
opcircles, 1., LINE=2, THICK=2
;
contourfill, ent[*,ny3,*],x,z, TITLE=tit+sy3+'!X', LEVELS=levent
var = reform(lam[*,ny3,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
opcircles, 1., LINE=2, THICK=2

wget

restore_state

end
; End of file vsections.pro
