;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vsections.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Nov-2001
;;;  $Id: vsections.pro,v 1.18 2003-08-08 20:09:59 dobler Exp $
;;;
;;;  Description:
;;;   Plot velocity, density and entropy field in three horizontal
;;;   sections.

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
    ophline, z, LINE=hline, THICK=3, COLOR=hcol1
    ophline, z, LINE=hline,          COLOR=hcol2
  endif
  if (lgravr) then begin
    opcircle, r, LINE=hline, THICK=3, COLOR=hcol1
    opcircle, r, LINE=hline,          COLOR=hcol2
  endif
end
; ---------------------------------------------------------------------- ;
function _aspect_3x2, ratio, npos, MARGIN=margin
;
;  A hack of D. Fanning's aspect.pro that works for !p.multi=[0,3,2].
;  Should rather fix aspect.pro to handle all kinds of !p.multi.
;  NPOS is the index of the sub-image, i.e. !p.multi is in fact [npos,3,2]
;
  if (n_params() eq 0) then ratio = 1.0

  if (ratio eq 0) then begin
    message, 'Aspect Ratio of 0. Changing to 1...', /INFO
    ratio = 1.0
  endif
  default, npos, 0
  nposx = npos mod 3
  nposy = npos / 3

  s = size(ratio)
  type = s(s(0)+1)
  if (type ne 4 and type ne 5) then $
      message, 'Aspect Ratio is not a FLOAT. Take care...', /INFO
  
  ; Check for margins:
  if (n_elements(margin) eq 0) then margin = 0.15

  ; Error checking:
  if ((margin lt 0) or (margin ge 0.5)) then $
      message, 'The MARGIN keyword value must be between 0.0 and 0.5.'
 
  ; Calculate the aspect ratio of the current window.  
  wratio = float(!d.y_vsize*3) / (!d.x_vsize*2)

  ; Calculate normalized positions in window.
  if (ratio le wratio) then begin
    xstart = 0.333*nposx + margin
    ystart = 0.75 - 0.5*nposy - (0.25 - margin) * (ratio / wratio)
    xend   = 0.333*(nposx+1) - margin
    yend   = 0.75 - 0.5*nposy + (0.25 - margin) * (ratio / wratio)
  endif else begin
    xstart = 0.1667 + 0.333*nposx - (0.1667 - margin) * (wratio / ratio)
    ystart = 0.5*(1-nposy) + margin
    xend   = 0.1667 + 0.333*nposx + (0.1667 - margin) * (wratio / ratio)
    yend   = 0.5*(2-nposy) - margin
  endelse

  position = [xstart, ystart, xend, yend]
  return, position

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
nsslevs = 60                   ; No of colours

ny1 = 0.25*my > 4
ny2 = 0.5*my
ny3 = 0.75*my < (my-5)

sy1 = '!8y!6='+strtrim(y[ny1],2)
sy2 = '!8y!6='+strtrim(y[ny2],2)
sy3 = '!8y=!6'+strtrim(y[ny3],2)

save_state

wput

!p.multi = [0,3,2]
if (!d.name eq 'X') then !p.charsize = 2
!x.title = '!8x!X'
!y.title = '!8z!X'

tit = '!17u!6 at '
!x.style = 1
!y.style = 1
if (show_ghosts) then begin
  !x.range = [x[0], x[mx-1]]    ; Include ghost zones
  !y.range = [z[0], z[mz-1]]
endif else begin
  !x.range = [x[3], x[mx-4]]    ; No ghost zones
  !y.range = [z[3], z[mz-4]]
endelse

if (absolute) then begin
  zruu = minmax(uu[*,*,*,1])
endif else begin
  undefine, zruu                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

ratio = (!y.range[1]-!y.range[0])/(!x.range[1]-!x.range[0])
plot_3d_vect, uu[*,ny1,*,*],x,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sy1+'!X', ZRANGE=zruu, $
    POSITION=_aspect_3x2(ratio, 0, MARGIN=0.05)
_opstuff, [z0,z1,z2,z3], sqrt(1-y[ny1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[*,ny2,*,*],x,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sy2+'!X', ZRANGE=zruu, $
    POSITION=_aspect_3x2(ratio, 1, MARGIN=0.05)
_opstuff, [z0,z1,z2,z3], sqrt(1-y[ny2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
plot_3d_vect, uu[*,ny3,*,*],x,z, PERM=[0,2,1], $
    /KEEP, TITLE=tit+sy3+'!X', ZRANGE=zruu, $
    POSITION=_aspect_3x2(ratio, 2, MARGIN=0.05)
_opstuff, [z0,z1,z2,z3], sqrt(1-y[ny3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

tit = '!8s!6 and '+s.varrho+'!6 at '

;
if (absolute) then begin
  levss = linspace(minmax(ss),nsslevs,/UNIQ)
endif else begin
  undefine, levss                 ; LEVELS=<undef> is like no LEVELS kw at all
endelse

contourfill, ss[*,ny1,*],x,z, TITLE=tit+sy1+'!X', LEVELS=levss, $
    POSITION=_aspect_3x2(ratio, 3, MARGIN=0.05)
var = reform(lnrho[*,ny1,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, [z0,z1,z2,z3], sqrt(1-y[ny1]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ss[*,ny2,*],x,z, TITLE=tit+sy2+'!X', LEVELS=levss, $
    POSITION=_aspect_3x2(ratio, 4, MARGIN=0.05)
var = reform(lnrho[*,ny2,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, [z0,z1,z2,z3], sqrt(1-y[ny2]^2), LGRAVZ=lgravz, LGRAVR=lgravr
;
contourfill, ss[*,ny3,*],x,z, TITLE=tit+sy3+'!X', LEVELS=levss, $
    POSITION=_aspect_3x2(ratio, 5, MARGIN=0.05)
var = reform(lnrho[*,ny3,*])
contour, var,x,z, /OVER, LEVELS=linspace(minmax(var),nrholevs,/UNIQ)
_opstuff, [z0,z1,z2,z3], sqrt(1-y[ny3]^2), LGRAVZ=lgravz, LGRAVR=lgravr

wget

restore_state

end
; End of file vsections.pro
