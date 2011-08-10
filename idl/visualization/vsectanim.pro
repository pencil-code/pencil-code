;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   vsectanim.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Nov-2001
;;;  $Id$
;;;
;;;  Description:
;;;   Scan through snapshots and plot velocity, density and entropy
;;;   field in three horizontal sections.


default, absolute, 0            ; flag four absolute colour scaling (i.e.
                                ; relative to absolute min and max of
                                ; colour-represented data)

nrholevs = 15                   ; No of isolines
nuulevs = 60                    ; No of colours
nentlevs = 60                   ; No of colours

ny1 = 0.2*ny > 4
ny2 = 0.5*ny
ny3 = 0.8*ny < ny-5

sy1 = '!8y!6='+strtrim(y[ny1],2)
sy2 = '!8y!6='+strtrim(y[ny2],2)
sy3 = '!8y=!6'+strtrim(y[ny3],2)

snapshots = file_search(datadir+"/VAR[0-9]*")
flist = strarr (1)
for len=min(strlen(snapshots)), max(strlen(snapshots)) do begin
	indices = where(strlen(snapshots) eq len)
	if (n_elements(indices) eq 1) then if (indices eq -1) then continue
	sub = snapshots[indices]
	reorder = sort(sub)
	flist = [ flist, sub[reorder] ]
end
flist = flist[1:*]

uu    = fltarr(nx,ny,nz,3)*one
lnrho = fltarr(nx,ny,nz)*one
ent   = fltarr(nx,ny,nz)*one

save_state

!p.multi = [0,3,2]
!p.charsize = 1.8
!x.title = '!8x!X'
!y.title = '!8z!X'

for i=0,n_elements(flist)-1 do begin
  file = flist[i]
  parts = strsplit(file,"/",/EXTRACT)
  file = parts[n_elements(parts)-1]

  if (file) then begin

    print, file

    close,1
    openr,1, datadir+'/'+file, /F77
    readu,1, uu, lnrho, ent
    if (i eq 0) then begin
      readu,1, t, x, y, z       ; Make sure corrdinates are right
    endif else begin
      readu,1, t
    endelse
    close,1

    wput

    tit = '!17u!6 at '

    if (absolute) then begin
      zruu = minmax(uu[*,*,*,1])
    endif else begin
      undefine, zruu            ; ZRANGE=<undef> is like no ZRANGE kw at all
    endelse

    plot_3d_vect, uu[*,ny1,*,*],x,z, PERM=[0,2,1], $
        /KEEP, TITLE=tit+sy1+' !8t!6='+strtrim(t,2)+'!X', ZRANGE=zruu
    opcircles, 1., LINE=2, THICK=2
    plot_3d_vect, uu[*,ny2,*,*],x,z, PERM=[0,2,1], $
        /KEEP, TITLE=tit+sy2+'!X', ZRANGE=zruu
    opcircles, 1., LINE=2, THICK=2
    plot_3d_vect, uu[*,ny3,*,*],x,z, PERM=[0,2,1], $
        /KEEP, TITLE=tit+sy3+'!X', ZRANGE=zruu
    opcircles, 1., LINE=2, THICK=2

    tit = '!8s!6 and !7r!6 at '
    !x.style = 1
    !y.style = 1
    !x.range = [x[3], x[nx-4]]  ; No ghost zones
    !y.range = [z[3], z[nz-4]]

    ;
    if (absolute) then begin
      levent = linspace(minmax(ent)+[-1,1]*1e-27,nentlevs,/UNIQ)
    endif else begin
      undefine, levent          ; LEVELS=<undef> is like no LEVELS kw at all
    endelse

    contourfill, ent[*,ny1,*],x,z, TITLE=tit+sy1+'!X', LEVELS=levent
    var = reform(lnrho[*,ny1,*])
    contour, var,x,z, /OVER, LEV=linspace(minmax(var),nrholevs,/UNIQ)
    opcircles, 1., LINE=2, THICK=2
    ;
    contourfill, ent[*,ny2,*],x,z, TITLE=tit+sy2+'!X', LEVELS=levent
    var = reform(lnrho[*,ny2,*])
    contour, var,x,z, /OVER, LEV=linspace(minmax(var),nrholevs,/UNIQ)
    opcircles, 1., LINE=2, THICK=2
    ;
    contourfill, ent[*,ny3,*],x,z, TITLE=tit+sy3+'!X', LEVELS=levent
    var = reform(lnrho[*,ny3,*])
    contour, var,x,z, /OVER, LEV=linspace(minmax(var),nrholevs,/UNIQ)
    opcircles, 1., LINE=2, THICK=2

    wget

  endif else print, "Empty file No. ", strtrim(i,2)

endfor

restore_state

end
; End of file vsectanim.pro
