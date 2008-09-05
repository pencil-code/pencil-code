;;;;;;;;;;;;;;;;;;;;;;
;;;   urms-z.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   01-Mar-2002
;;;  $Id$
;;;
;;;  Description:
;;;   Plot the rms velocity as function of z

urms1 = (urms = fltarr(nz))
umeanz = fltarr(nz,3)

for i=0,nz-1 do begin
  for k=0,2 do begin
    umeanz[i,k] = mean(uu[nghx:nx-nghx-1,nghy:ny-nghy-1,i,k])
  endfor
  urms[i] = rms(uu[nghx:nx-nghx-1,nghy:ny-nghy-1,i,*],/VECT)
  urms1[i] = rms( uu[nghx:nx-nghx-1,nghy:ny-nghy-1,i,*] $
                  - spread(umeanz[i,*],[0,1],[nx-2*nghx,ny-2*nghy]), /VECT )
endfor
umeanabs = sqrt(umeanz[*,0]^2+umeanz[*,1]^2+umeanz[*,2]^2)

save_state

!p.multi = [0,2,1]
!y.style = 3

xr = minmax([urms1,urms,umeanabs])
plot, urms, z, /XLOG, $
    XTITLE='!8u!D!6rms!N!X', XRANGE=xr, YTITLE='!8z!X'
oplot, urms1, z, line=2
oplot, umeanabs, z, LINE=1
ophline, [z0,z1,z2,ztop]

xr = minmax([0,xr])
plot, urms, z, $
    XTITLE='!8u!D!6rms!N!X', XRANGE=xr, YTITLE='!8z!X'
oplot, urms1, z, line=2
oplot, umeanabs, z, LINE=1
ophline, [z0,z1,z2,ztop]

esrg_legend, $
    ['!6rms(!17u!6)!X', '!6rms(!17u!6-<!17u!6>)', '!3|<!17u!3>|!X'], $
    LINE=[0,2,1], SPOS='br', /BOX


restore_state

end
; End of file urms-z.pro
