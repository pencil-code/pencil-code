;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vsections.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Nov-2001
;;;
;;;  Description:
;;;   Plot velocity, density and entropy field in three horizontal
;;;   sections. Same as vsections, but cuts x=const are shown

@symbols

default, absolute, 0            ; flag four absolute colour scaling (i.e.
                                ; relative to absolute min and max of
                                ; colour-represented data
default, show_ghosts, 0

nrholevs = 15                   ; No of isolines
nuulevs = 60                    ; No of colours
nentlevs = 60                   ; No of colours

nx1 = 0.2*nx > 4
nx2 = 0.5*nx
nx3 = 0.8*nx < nx-5

sx1 = '!8x!6='+strtrim(x[nx1],2)
sx2 = '!8x!6='+strtrim(x[nx2],2)
sx3 = '!8x=!6'+strtrim(x[nx3],2)

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
  !x.range = [y[0], y[ny-1]]    ; No ghost zones
  !y.range = [z[0], z[nz-1]]
endif else begin
  !x.range = [y[3], y[ny-4]]    ; No ghost zones
  !y.range = [z[3], z[nz-4]]
endelse

if (absolute) then begin
  zruu = minmax(uu[*,*,*,1])
endif else begin
  undefine, zruu                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

hline = 2                       ; line style for vertical boundaries
hcol1 = !p.color                ; color for " "
hcol2 = !p.background           ; color for " "

plot_3d_vect, uu[nx1,*,*,*],y,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sx1+'!X', ZRANGE=zruu
;opcircles, 1., LINE=2, THICK=2
ophline,[z0,z1,z2,z3], LINE=hline, THICK=3, COLOR=hcol1
ophline,[z0,z1,z2,z3], LINE=hline, COLOR=hcol2
;
plot_3d_vect, uu[nx2,*,*,*],y,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sx2+'!X', ZRANGE=zruu
;opcircles, 1., LINE=2, THICK=2
ophline,[z0,z1,z2,z3], LINE=hline, THICK=3, COLOR=hcol1
ophline,[z0,z1,z2,z3], LINE=hline, COLOR=hcol2
;
plot_3d_vect, uu[nx3,*,*,*],y,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sx3+'!X', ZRANGE=zruu
;opcircles, 1., LINE=2, THICK=2
ophline,[z0,z1,z2,z3], LINE=hline, THICK=3, COLOR=hcol1
ophline,[z0,z1,z2,z3], LINE=hline, COLOR=hcol2

tit = '!8s!6 and '+s_varrho+'!6 at '

;
if (absolute) then begin
  levent = linspace(minmax(ent),nentlevs)
endif else begin
  undefine, levent                 ; LEVELS=<undef> is like no LEVELS kw at all
endelse

contourfill, ent[nx1,*,*],y,z, TITLE=tit+sx1+'!X', LEVELS=levent
var = reform(lam[nx1,*,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
;opcircles, 1., LINE=2, THICK=2
ophline,[z0,z1,z2,z3], LINE=hline, THICK=3, COLOR=hcol1
ophline,[z0,z1,z2,z3], LINE=hline, COLOR=hcol2
;
contourfill, ent[nx2,*,*],y,z, TITLE=tit+sx2+'!X', LEVELS=levent
var = reform(lam[nx2,*,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
;opcircles, 1., LINE=2, THICK=2
ophline,[z0,z1,z2,z3], LINE=hline, THICK=3, COLOR=hcol1
ophline,[z0,z1,z2,z3], LINE=hline, COLOR=hcol2
;
contourfill, ent[nx3,*,*],y,z, TITLE=tit+sx3+'!X', LEVELS=levent
var = reform(lam[nx3,*,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
;opcircles, 1., LINE=2, THICK=2
ophline,[z0,z1,z2,z3], LINE=hline, THICK=3, COLOR=hcol1
ophline,[z0,z1,z2,z3], LINE=hline, COLOR=hcol2

wget

restore_state

end
; End of file vsections.pro
