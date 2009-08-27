;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   hsections.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Nov-2001
;;;  $Id$
;;;
;;;  Description:
;;;   Plot velocity, density and entropy field in three horizontal
;;;   sections.

; ---------------------------------------------------------------------- ;
pro _opstuff, r, LGRAVZ=lgravz, LGRAVR=lgravr
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

nz1 = 0.25*mz > 4
nz2 = 0.5 *mz
nz3 = 0.75*mz < (mz-5)

sz1 = '!8z!6='+strtrim(z[nz1],2)
sz2 = '!8z!6='+strtrim(z[nz2],2)
sz3 = '!8z=!6'+strtrim(z[nz3],2)

save_state

wput

!p.multi = [0,3,2]
if (!d.name eq 'X') then !p.charsize = 2
!x.title = '!8x!X'
!y.title = '!8y!X'

tit = '!17u!6 at '
!x.style = 1
!y.style = 1
if (show_ghosts) then begin
  !x.range = [x[0], x[mx-1]]    ; Include ghost zones
  !y.range = [y[0], y[my-1]]
endif else begin
  !x.range = [x[3], x[mx-4]]    ; No ghost zones
  !y.range = [y[3], y[my-4]]
endelse

if (absolute) then begin
  zruu = minmax(uu[*,*,*,2])
endif else begin
  undefine, zruu                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

ratio = (!y.range[1]-!y.range[0])/(!x.range[1]-!x.range[0])
plot_3d_vect, uu[*,*,nz1,*],x,y, $
    /KEEP, TITLE=tit+sz1+'!X', ZRANGE=zruu, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
_opstuff, sqrt(1-z[nz1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[*,*,nz2,*],x,y, $
    /KEEP, TITLE=tit+sz2+'!X', ZRANGE=zruu, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
_opstuff, sqrt(1-z[nz2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[*,*,nz3,*],x,y, $
    /KEEP, TITLE=tit+sz3+'!X', ZRANGE=zruu, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
_opstuff, sqrt(1-z[nz3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

tit = '!8s!6 and '+s.varrho+'!6 at '

;
if (absolute) then begin
  levss = linspace(minmax(ss),nsslevs,/UNIQ)
endif else begin
  undefine, levss                 ; LEVELS=<undef> is like no LEVELS kw at all
endelse

contourfill, ss[*,*,nz1],x,y, TITLE=tit+sz1+'!X', LEVELS=levss, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
var = reform(lnrho[*,*,nz1])
contour, var, x, y, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, sqrt(1-z[nz1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ss[*,*,nz2],x,y, TITLE=tit+sz2+'!X', LEVELS=levss, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
var = reform(lnrho[*,*,nz2])
contour, var, x, y, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, sqrt(1-z[nz2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ss[*,*,nz3],x,y, TITLE=tit+sz3+'!X', LEVELS=levss, $
    POSITION=aspect_pos(ratio,MARGIN=0.1)
var = reform(lnrho[*,*,nz3])
contour, var, x, y, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, sqrt(1-z[nz3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

wget

restore_state

end
; End of file hsections.pro
