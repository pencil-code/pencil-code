;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vsections2.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Nov-2001
;;;
;;;  Description:
;;;   Plot velocity, density and entropy field in three horizontal
;;;   sections. Same as vsections, but cuts x=const are shown

; ---------------------------------------------------------------------- ;
pro _opstuff, z, r, LGRAVz=lgravz, lgravr=lgravr
;
;  Overplot some lines or circles (needed repeatedly). Plot white on
;  black, so lines can be distinguished on any background
;
  hline = 2                     ; line style
  hcol1 = !p.background         ; color1
  hcol2 = !p.color              ; color2
;
  default, lgravz, 0L
  default, lgravr, 0L
  if (lgravz)  then begin
    ophline,[z0,z1,z2,z3], LINE=hline, THICK=3, COLOR=hcol1
    ophline,[z0,z1,z2,z3], LINE=hline,          COLOR=hcol2
  endif
  if (lgravr) then begin
    opcircle, r, LINE=hline, THICK=3, COLOR=hcol1
    opcircle, r, LINE=hline,          COLOR=hcol2
  endif
end
; ---------------------------------------------------------------------- ;

s = texsyms()

default, absolute, 0            ; flag four absolute colour scaling (i.e.
                                ; relative to absolute min and max of
                                ; colour-represented data
default, show_ghosts, 0

nrholevs = 15                   ; No of isolines
nuulevs = 60                    ; No of colours
nentlevs = 60                   ; No of colours

nx1 = 0.25*nx > 4
nx2 = 0.5*nx
nx3 = 0.75*nx < (nx-5)

sx1 = '!8x!6='+strtrim(x[nx1],2)
sx2 = '!8x!6='+strtrim(x[nx2],2)
sx3 = '!8x=!6'+strtrim(x[nx3],2)

save_state

wput

!p.multi = [0,3,2]
!p.charsize = 2
!x.title = '!8y!X'
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

plot_3d_vect, uu[nx1,*,*,*],y,z, PERM=[1,2,0], $
    /KEEP, TITLE=tit+sx1+'!X', ZRANGE=zruu
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[nx2,*,*,*],y,z, PERM=[1,2,0], $
    /KEEP, TITLE=tit+sx2+'!X', ZRANGE=zruu
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[nx3,*,*,*],y,z, PERM=[1,2,0], $
    /KEEP, TITLE=tit+sx3+'!X', ZRANGE=zruu
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

tit = '!8s!6 and '+s.varrho+'!6 at '

;
if (absolute) then begin
  levent = linspace(minmax(ent),nentlevs)
endif else begin
  undefine, levent                 ; LEVELS=<undef> is like no LEVELS kw at all
endelse

contourfill, ent[nx1,*,*],y,z, TITLE=tit+sx1+'!X', LEVELS=levent
var = reform(lnrho[nx1,*,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ent[nx2,*,*],y,z, TITLE=tit+sx2+'!X', LEVELS=levent
var = reform(lnrho[nx2,*,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ent[nx3,*,*],y,z, TITLE=tit+sx3+'!X', LEVELS=levent
var = reform(lnrho[nx3,*,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

wget

restore_state

end
; End of file vsections2.pro
