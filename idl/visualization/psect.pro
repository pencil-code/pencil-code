;;;;;;;;;;;;;;;;;;;;;
;;;   psect.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   09-Sep-2001
;;;  $Id$
;;;
;;;  Description:
;;;   Plot a section of the box, indicating processor boundaries and
;;;   ghost zones.

default, cpu_grid, 1            ; flag for plotting CPU layout
default, iz, nz/2+nghx          ; position of vertical section
default, screen_bounds, 0       ; flag to screen out bb in boundary points

;!p.multi=0
;plot_3d_vect, uu[*,*,iz,0],uu[*,*,iz,1],lam[*,*,iz], x,y, $
;    XSTYLE=1,YSTYLE=1

save_state                      ; save current graphics state

!p.multi=[0,3,2]
!p.charsize=2
!x.style=1 & !y.style=1
!x.title='!8x!X'
!y.title='!8y!X'

;; 1. velocity
wdvelovect,uu[*,*,iz,0],uu[*,*,iz,1], x,y, $
    TITLE='!6Velocity!X'

;; 2. density
contourfill, lam[*,*,iz], x,y, $
    TITLE='!6Density!X'
if (cpu_grid) then begin
  for i=0,nprocx do begin
    opvline, x[i*nxloc+nghx]-dx/2, LINE=0
    opvline, x[i*nxloc+nghx]+2.5*dx, LINE=2
    opvline, x[i*nxloc+nghx]+(nxloc-3.4)*dx, LINE=2
  endfor
  ;
  for i=0,nprocy do begin
    ophline, y[i*nyloc+nghy]-dy/2, LINE=0
    ophline, y[i*nyloc+nghy]+2.5*dy, LINE=2
    ophline, y[i*nyloc+nghy]+(nyloc-3.4)*dy, LINE=2
  endfor
endif

;; 3. Magnetic field
aa_full=aa
aa_full[*,*,*,1] = aa[*,*,*,1]+xx*bz_ext
aa_full[*,*,*,2] = aa[*,*,*,2]-xx*by_ext+yy*bx_ext

;contour, aa_full[*,*,iz,2], x, y, $
;    LEVELS=linspace(minmax(aa_full[*,*,iz,2]),30),/OVER
;contour, aa[*,*,iz,2], x, y, $
;    LEVELS=linspace(minmax(aa[*,*,iz,2]),30),/OVER
bb_full = curl(aa_full)
; remove bb_full in ghost zones (which is irrelevant and inaccurate)
if (screen_bounds) then begin
  bb_full[0:nghx-1,*,*,*]=0 & bb_full[nxtot-nghx:nxtot-1,*,*,*]=0
  bb_full[*,0:nghx-1,*,*]=0 & bb_full[*,nxtot-nghx:nxtot-1,*,*]=0
  bb_full[*,*,0:nghx-1,*]=0 & bb_full[*,*,nxtot-nghx:nxtot-1,*]=0
endif
;
wdvelovect, bb_full[*,*,iz,0],bb_full[*,*,iz,1], x,y, $
    TITLE='!17B!6-field!X'

;; 4. Magnetic field minus background field
bb = curl(aa)
; remove bb_full in ghost zones (which is irrelevant and inaccurate)
if (screen_bounds) then begin
  bb[0:nghx-1,*,*,*]=0 & bb[nxtot-nghx:nxtot-1,*,*,*]=0
  bb[*,0:nghx-1,*,*]=0 & bb[*,nxtot-nghx:nxtot-1,*,*]=0
  bb[*,*,0:nghx-1,*]=0 & bb[*,*,nxtot-nghx:nxtot-1,*]=0
endif
;
wdvelovect, bb[*,*,iz,0],bb[*,*,iz,1], x,y, $
    TITLE='!17B!6-field!X'

;; 5. Profile of ux
plot, x, uu[10,*,10,0], $
    YSTYLE=3, $
    XTITLE='!8x!X', YTITLE='!8u!Dx!N!X', TITLE='!6Velocity!X'
oplot, x, uu[10,*,10,1], LINE=2

;; 6. Profile of bx
plot, x, bb[10,*,10,0], $
    YSTYLE=3, $
    XTITLE='!8x!X', YTITLE='!8B!Dx!N!X', TITLE='!17B!6-field!X'



restore_state                   ; restore graphics state

end
; End of file psect.pro
