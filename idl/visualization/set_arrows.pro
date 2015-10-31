pro set_arrows, info, xy, hsize=hsize, indmin=indmin, refpoint=refpoint, tol=tol, color=color, ref2point=ref2point
;
; sets arrows on isolines found by contour described by its output parameters path_info (info) and path_xy (xy)
;
; hsize     -  = hsize of routine arrow
; indmin    - number of first line to be considered
; refpoint  - arrows are set at minimum distance to refpoint on each line
; ref2point - if provided minimumdistance is also calculated w.r.t. refpoint2; arrow is set at point with the smaller of these two distances
; tol       - smallest relative distance to box boundaries at which arrow is still set
;
; apr-2014/MR: coded
;  
if not keyword_set(indmin) then indmin=0

if not keyword_set(refpoint) then refpoint=[0.,0.]

if not keyword_set(hsize) then $
  if !d.name eq 'PS' then $
    hsize=220 $
  else $
    hsize=12

if not keyword_set(tol) then tol=0.01

range=[min(xy(0,*)),min(xy(1,*)),max(xy(0,*)),max(xy(1,*))]
widths=[range(2)-range(0), range(3)-range(1)]
center=[mean([range(0),range(2)]), refpoint(1)]
;ref2point = 2.*center - refpoint
refpoint_dev=reform(convert_coord(refpoint,/data,/to_device))
if keyword_set(ref2point) then ref2point_dev=reform(convert_coord(ref2point,/data,/to_device))

for i=indmin,n_elements(info)-1 do begin

  np = info(i).n
  ia = info(i).offset

  dists=fltarr(np) & dists2=dists & mchord2=0. & firstp_dev=reform(convert_coord(xy(*,ia),/data,/to_device))

istop=21
;if i eq istop then print, 'ia, imin, np, ja, je, js=', ia, imin, np, ja, je, js
  for j=ia,ia+np-1 do begin
    ;if i eq istop then print, j
    runp_dev=reform(convert_coord(xy(*,j),/data,/to_device))
    dists(j-ia) = total((runp_dev-refpoint_dev)^2)
    if keyword_set(ref2point) then dists2(j-ia) = total((runp_dev-ref2point_dev)^2)
    mchord2 = mchord2 > total((firstp_dev-runp_dev)^2)
  endfor

  if hsize gt .7*sqrt(mchord2) then continue
  val=min(dists,imin)

  if keyword_set(ref2point) then begin
    val2=min(dists2,imin2)
    if val2 lt val then imin=imin2
    ;print, 'imin2=', imin2
  endif

  pt1=xy(*,ia+imin)
  pt1_dev=reform(convert_coord(pt1,/data,/to_device))

  if imin-ia lt np/2. then begin
    ja=-1 & je=-np+imin+1 & js=-1
  endif else begin
    ja=1 & je=imin & js=1
  endelse

;if i eq istop  then print, 'ia, imin, np, ja, je, js=', ia, imin, np, ja, je, js
  for j=ja,je,js do begin
    pt2=xy(*,ia+imin-j)
    pt2_dev=reform(convert_coord(pt2,/data,/to_device))
    if total((pt1_dev-pt2_dev)^2) gt hsize^2 then begin
      pt2=xy(*,ia+imin-j+js)
      break
    endif
  endfor

  if js eq -1 then begin
   temp=pt1 & pt1=pt2 & pt2=temp
  endif

  if abs(pt1(0)-range(0)) gt tol*widths(0) and abs(pt1(1)-range(1)) gt tol*widths(1) and $
     abs(pt1(0)-range(2)) gt tol*widths(0) and abs(pt1(1)-range(3)) gt tol*widths(1) and $
     abs(pt2(0)-range(0)) gt tol*widths(0) and abs(pt2(1)-range(1)) gt tol*widths(1) and $
     abs(pt2(0)-range(2)) gt tol*widths(0) and abs(pt2(1)-range(3)) gt tol*widths(1) then begin

    arrow, pt2(0), pt2(1), pt1(0), pt1(1), /data, /solid, hsize=hsize, color=color
  endif else begin
    ;arrow, pt1(0), pt1(1), pt2(0), pt2(1), /data, /solid, hsize=hsize, color=0
    ;print, 'black, i=', i
  endelse
    ;print, 'i=', i
    ;if i eq istop then stop
endfor

end

