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

ny1 = 0.2*ny > 4
ny2 = 0.5*ny
ny3 = 0.8*ny < ny-5

sy1 = '!8y!6='+strtrim(y[ny1],2)
sy2 = '!8y!6='+strtrim(y[ny2],2)
sy3 = '!8y=!6'+strtrim(y[ny3],2)

save_state
!p.multi = [0,3,2]
!p.charsize = 1.8
!x.title = '!8x!X'
!y.title = '!8z!X'

tit = '!17u!6 at '
plot_3d_vect, uu[*,ny1,*,*],x,z, PERM=[0,2,1], /KEEP, TITLE=tit+sy1+'!X'
plot_3d_vect, uu[*,ny2,*,*],x,z, PERM=[0,2,1], /KEEP, TITLE=tit+sy2+'!X'
plot_3d_vect, uu[*,ny3,*,*],x,z, PERM=[0,2,1], /KEEP, TITLE=tit+sy3+'!X'

tit = '!8s!6 and !7r!6 at '
!x.style = 1
!y.style = 1
!x.range = [x[3], x[nx-4]]      ; No ghost zones
!y.range = [z[3], z[nz-4]]

nrholevs = 15
;
contourfill, ent[*,ny1,*],x,z, TITLE=tit+sy1+'!X'
var = reform(lam[*,ny1,*])
contour, var,x,z, /OVER, LEV=linspace(minmax(var),nrholevs)
;
contourfill, ent[*,ny2,*],x,z, TITLE=tit+sy2+'!X'
var = reform(lam[*,ny2,*])
contour, var,x,z, /OVER, LEV=linspace(minmax(var),nrholevs)
;
contourfill, ent[*,ny3,*],x,z, TITLE=tit+sy3+'!X'
var = reform(lam[*,ny3,*])
contour, var,x,z, /OVER, LEV=linspace(minmax(var),nrholevs)


restore_state

end
; End of file hsections.pro
