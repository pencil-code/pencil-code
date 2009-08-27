;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   bbsections.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   25-Jul-2002
;;;
;;;  Description:
;;;   Plot magnetic field in various sections

; ---------------------------------------------------------------------- ;
pro _opstuff, r, LGRAVz=lgravz, lgravr=lgravr
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

default, t_bb, -1e20            ; time of current bb

if ((t_bb ne t) or (n_elements(bb) lt n_elements(aa))) then begin
  print, 'bb outdated.. recalculating'
  bb = curl(aa)
  t_bb = t
endif

default, absolute, 0            ; flag four absolute colour scaling (i.e.
                                ; relative to absolute min and max of
                                ; colour-represented data
default, show_ghosts, 0

default, z1,0                   ; possibly needed later
default, z2,0
default, z3,0

if (par.r_ext gt 0.) then begin ; choose planes based on geometry
  nx0 = (where(abs(x) eq min(abs(x))))[0]
  ny0 = (where(abs(y) eq min(abs(y))))[0]
  nz0 = (where(abs(z) eq min(abs(z))))[0]
  nxi = nx0 + [-1, 0, 1]*0.7*par.r_ext/dx + 1
  nyi = ny0 + [-1, 0, 1]*0.7*par.r_ext/dy + 1
  nzi = nz0 + [-1, 0, 1]*0.7*par.r_ext/dz + 1
endif else begin                ; choose plains based on grid
  nxi = [0.25,0.5,0.75]*mx
  nyi = [0.25,0.5,0.75]*my
  nzi = [0.25,0.5,0.75]*mz
endelse

;; Sanitize
nxi = nxi > 3 < (mx-4)
nyi = nyi > 3 < (my-4)
nzi = nzi > 3 < (mz-4)

sxi = '!8x!6='+strtrim(x[nxi],2)+'!X'
syi = '!8y!6='+strtrim(y[nyi],2)+'!X'
szi = '!8z!6='+strtrim(z[nzi],2)+'!X'

save_state

wput

!p.multi = [0,3,3]
!p.charsize = 2

tit = '!17B!6 at '
!x.style = 1
!y.style = 1
if (show_ghosts) then begin
  xrange = [x[0], x[mx-1]]    ; Include ghost zones
  yrange = [y[0], y[my-1]]
  zrange = [z[0], z[mz-1]]
endif else begin
  xrange = [x[3], x[mx-4]]    ; No ghost zones
  yrange = [y[3], y[my-4]]
  zrange = [z[3], z[mz-4]]
endelse

zrb = [-1,1]*max(sqrt(dot2(bb)))

if (absolute) then begin
  zrbb = zrb
endif else begin
  undefine, zrbb                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse


!x.title = '!8y!X'
!y.title = '!8z!X'
for i=0,2 do begin
  pos = aspect_pos((zrange[1]-zrange[0])/(yrange[1]-yrange[0]),MARGIN=0.09)
  plot_3d_vect, bb[nxi[i],*,*,*], y, z, $
      PERM=[1,2,0], /KEEP_CT, $
      TITLE=tit+sxi[i]+'!X', $
      XRANGE=yrange, XSTYLE=1, $
      YRANGE=zrange, YSTYLE=1, $
      ZRANGE=zrbb, $
      POS=pos
  _opstuff, sqrt(1-x[nxi[i]]^2), $
      LGRAVZ=lgravz, LGRAVR=lgravr
endfor

if (absolute) then begin
  zrbb = zrb
endif else begin
  undefine, zrbb                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

!x.title = '!8x!X'
!y.title = '!8z!X'
for i=0,2 do begin
  pos = aspect_pos((zrange[1]-zrange[0])/(xrange[1]-xrange[0]),MARGIN=0.09)
  plot_3d_vect, bb[*,nyi[i],*,*], x, z, $
      PERM=[0,2,1], /KEEP_CT, $
      TITLE=tit+syi[i]+'!X', $
      XRANGE=xrange, XSTYLE=1, $
      YRANGE=zrange, YSTYLE=1, $
      ZRANGE=zrbb, $
      POS=pos
  _opstuff, sqrt(1-y[nyi[i]]^2), $
      LGRAVZ=lgravz, LGRAVR=lgravr
endfor

if (absolute) then begin
  zrbb = zrb
endif else begin
  undefine, zrbb                ; ZRANGE=<undef> is like no ZRANGE kw at all
endelse

!x.title = '!8x!X'
!y.title = '!8y!X'
for i=0,2 do begin
  pos = aspect_pos((yrange[1]-yrange[0])/(xrange[1]-xrange[0]),MARGIN=0.09)
  plot_3d_vect, bb[*,*,nzi[i],*], x, y, $
      PERM=[0,1,2], /KEEP_CT, $
      TITLE=tit+szi[i]+'!X', $
      XRANGE=xrange, XSTYLE=1, $
      YRANGE=yrange, YSTYLE=1, $
      ZRANGE=zrbb, $
      POS=pos
  _opstuff, sqrt(1-z[nzi[i]]^2), $
      LGRAVZ=lgravz, LGRAVR=lgravr
endfor

wget

restore_state

end
; End of file bbsections.pro
