;;;;;;;;;;;;;;;;;;;;;
;;;   pvert.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   11-Nov-2001
;;;
;;;  Description:
;;;   Plot vertical profiles of uz, lnrho and entropy.

nign = 3                        ; Number of close-to-bdry points to ignore

save_state

!p.multi = [0,3,1]
!p.charsize = 1.6
!y.title = '!8z!X'

for ivar = 0,2 do begin
  case ivar of
    0: begin
      var = lam
      title = 'log density'
    end
    1: begin
      var = uu[*,*,*,2]
      title = '!8u!Dz!N!X'
    end
    2: begin
      var = ent
      title = 'Entropy'
    end
  endcase
  plot, z, z, /NODATA,XRANGE=minmax(var),YRANGE=minmax(z),TITLE=title
  for ix=nign,nx-nign-1 do begin
    for iy=nign,ny-nign-1 do begin
      oplot, var[ix,iy,*], z
    endfor
    ophline, [z[3], z[nz-4]]
    if (ivar eq 1) then opvline
  endfor
endfor

restore_state

end
; End of file pvert.pro
