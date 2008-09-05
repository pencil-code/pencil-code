;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vsections2.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Nov-2001
;;;  $Id$
;;;
;;;  Description:
;;;   Plot velocity, density and entropy field in three vertical
;;;   sections (x=const).

; ---------------------------------------------------------------------- ;
pro _opstuff, z, r, LGRAVZ=lgravz, LGRAVR=lgravr
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
  default, zlevels, [0]
  if (lgravz)  then begin
    ophline,z, LINE=hline, THICK=3, COLOR=hcol1
    ophline,z, LINE=hline,          COLOR=hcol2
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

default, z1,0                   ; neede later
default, z2,0
default, z3,0

nrholevs = 15                   ; No of isolines
nuulevs = 60                    ; No of colours
nsslevs = 60                    ; No of colours

nx1 = 0.25*mx > 4
nx2 = 0.5*mx
nx3 = 0.75*mx < (mx-5)

sx1 = '!8x!6='+strtrim(x[nx1],2)
sx2 = '!8x!6='+strtrim(x[nx2],2)
sx3 = '!8x=!6'+strtrim(x[nx3],2)

save_state

wput

!p.multi = [0,3,2]
if (!d.name eq 'X') then !p.charsize = 2
!x.title = '!8y!X'
!y.title = '!8z!X'

tit = '!17u!6 at '
!x.style = 1
!y.style = 1
if (show_ghosts) then begin
  !x.range = [y[0], y[my-1]]    ; Include ghost zones
  !y.range = [z[0], z[mz-1]]
endif else begin
  !x.range = [y[3], y[my-4]]    ; No ghost zones
  !y.range = [z[3], z[mz-4]]
endelse

if (absolute) then begin
  zruu = minmax(uu[*,*,*,0])
endif else begin
  undefine, zruu                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

ratio = (!y.range[1]-!y.range[0])/(!x.range[1]-!x.range[0])
plot_3d_vect, uu[nx1,*,*,*],y,z, PERM=[1,2,0], $
    /KEEP, TITLE=tit+sx1+'!X', ZRANGE=zruu, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[nx2,*,*,*],y,z, PERM=[1,2,0], $
    /KEEP, TITLE=tit+sx2+'!X', ZRANGE=zruu, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[nx3,*,*,*],y,z, PERM=[1,2,0], $
    /KEEP, TITLE=tit+sx3+'!X', ZRANGE=zruu, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

tit = '!8s!6 and '+s.varrho+'!6 at '

;
if (absolute) then begin
  levss = linspace(minmax(ss),nsslevs,/UNIQ)
endif else begin
  undefine, levss                 ; LEVELS=<undef> is like no LEVELS kw at all
endelse

contourfill, ss[nx1,*,*],y,z, TITLE=tit+sx1+'!X', LEVELS=levss, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
var = reform(lnrho[nx1,*,*])
contour, var,y,z, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ss[nx2,*,*],y,z, TITLE=tit+sx2+'!X', LEVELS=levss, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
var = reform(lnrho[nx2,*,*])
contour, var,y,z, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ss[nx3,*,*],y,z, TITLE=tit+sx3+'!X', LEVELS=levss, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
var = reform(lnrho[nx3,*,*])
contour, var,y,z, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, [z0,z1,z2,z3], sqrt(1-x[nx3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

wget

restore_state

end
; End of file vsections2.pro
