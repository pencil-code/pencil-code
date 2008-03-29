pro fieldintegrator_set_field,newfield,x,y,z,xyz0=xyz0_,dx=dx_,dy=dy_,dz=dz_,nx=nx_,ny=ny_,nz=nz_
  common vizit_fieldintegrator,field,xx,yy,zz,dx,dy,dz,xyz0,nx,ny,nz,dx1,dy1,dz1

  field=newfield ;/sqrt(max(dot2(transpose(newfield))))

  nx=n_elements(x)
  ny=n_elements(y)
  nz=n_elements(z)

  dx=0.
  dy=0.
  dz=0.
  if (nx gt 1) then dx=x[1]-x[0]
  if (ny gt 1) then dy=y[1]-y[0]
  if (nz gt 1) then dz=z[1]-z[0]

  dx1=1./dx
  dy1=1./dy
  dz1=1./dz

  xyz0=[x[0],y[0],z[0]]

  xyz0_=xyz0
  dx_=dx
  dy_=dy
  dz_=dz
  nx_=nx
  ny_=ny
  nz_=nz

end

function fieldintegrator_get_field,t,x
  common vizit_fieldintegrator,field,xx,yy,zz,dx,dy,dz,xyz0,nx,ny,nz,dx1,dy1,dz1

;  if (nx gt 1) then px = (x[0] - xyz0[0])*dx1
;  if (ny gt 1) then py = (x[1] - xyz0[1])*dy1
;  if (nz gt 1) then pz = (x[2] - xyz0[2])*dz1

  result=fltarr(3)
  result[0] = interpolate(reform(field[0,*,*,*]),x[0],x[1],x[2])
  result[1] = interpolate(reform(field[1,*,*,*]),x[0],x[1],x[2])
  result[2] = interpolate(reform(field[2,*,*,*]),x[0],x[1],x[2])

  result=result/sqrt(dot2(result))

  return, result
end


function fieldintegrator,field,x=x,y=y,z=z,seeds=seeds,reverse=reverse,maxittr=maxittr

  s=size(field)

  default,seeds,[ [0.,0.,0.] ]
 
 
  sseeds=size(seeds)
  if sseeds[0] eq 1 then begin
    nseeds=1L
  endif else if sseeds[0] eq 2 then begin
    nseeds=sseeds[2]
  endif else begin
    message,"Seeds array must be of size [3] or [3,n]" 
  endelse

  npnts=100L
  default,maxittr,10000L

  path=fltarr(3,npnts,nseeds)
  tangents=fltarr(3,npnts,nseeds)
  lineel=fltarr(npnts,nseeds)
  h=fltarr(nseeds)
  pathlen=lonarr(nseeds)

  fieldintegrator_set_field,field,x,y,z,xyz0=xyz0,dx=dx,dy=dy,dz=dz,nx=nx,ny=ny,nz=nz

  dxmax=(x[1]-x[0])/2
  dxmin=(min([x[1]-x[0],y[1]-y[0],z[1]-z[0]]))

  h[*]=0.9
  h[*]=0.5
 ; h[*]=0.5*dxmin
;  h[*]=10.*dxmin

  if max(abs(h)) eq 0. then begin
    print,'Zero stepsize!! Exiting'
    return,0.
  endif 

  if keyword_set(reverse) then h = -h
 
  lineel[0,0:nseeds-1]=0.0
  path[0:2,0,0:nseeds-1]=(seeds-spread(xyz0,1,nseeds))/spread(spread(dxmin,0,3),1,nseeds)
  for iseed=0,nseeds-1 do begin
    tangents[0:2,0,iseed]=fieldintegrator_get_field(lineel[0,iseed],path[0:2,0,iseed])
  endfor 

  finished=0
  ipnt=0L
  while not finished do begin
    if (ipnt+10) ge npnts then begin
      old_path=path
      path=fltarr(3,npnts+100,nseeds)
      path[*,0:npnts-1,*]=old_path
      old_path=0

      old_tangents=tangents
      tangents=fltarr(3,npnts+100,nseeds)
      tangents[*,0:npnts-1,*]=old_tangents
      old_tangents=0

      old_lineel=lineel
      lineel=fltarr(npnts+100,nseeds)
      lineel[0:npnts-1,*]=old_lineel
      old_lineel=0

      npnts=npnts+100L
    endif

    for iseed=0L,nseeds-1 do begin
      if h[iseed] eq 0. then continue
      tangents[0:2,ipnt,iseed]=fieldintegrator_get_field(lineel[ipnt,iseed],path[0:2,ipnt,iseed])
      XX=path[0:2,ipnt,iseed]
      DXXDTT=tangents[0:2,ipnt,iseed]
      TT=lineel[ipnt,iseed]

      result=RK4(XX,DXXDTT,TT,h[iseed],'fieldintegrator_get_field') 

if iseed eq 0 then print,result
      ds2=dot2(result-XX)

      if ( ds2 eq 0. ) then begin
        h[iseed]=0.
        pathlen[iseed]=ipnt-1
        continue
      endif

      iresult=long(result)
      if (iresult[0] lt 0) or (iresult[0] gt nx) $
        or (iresult[1] lt 0) or (iresult[1] gt ny) $
        or (iresult[2] lt 0) or (iresult[2] gt nz) then begin
        h[iseed]=0.
        pathlen[iseed]=ipnt-1
        continue
      endif

;help,t
;print,ipnt,iseed
      lineel[ipnt+1,iseed]=lineel[ipnt,iseed]+h[iseed]
      path[0:2,ipnt+1,iseed]=result
      pathlen[iseed]=ipnt+1
;print,result

;      print,reform([lineel[ipnt,iseed],result])
      ;print,reform([lineel[ipnt,iseed],tangents[0:2,ipnt,iseed]])
      ;print,reform([lineel[ipnt,iseed],dx/dxmax])

    endfor 
    ipnt=ipnt+1
    if ipnt ge maxittr then finished=1
    if max(abs(h)) eq 0. then finished=1
    if (ipnt mod 50) eq 0 then print,ipnt
  endwhile

  for iseed=0,nseeds-1 do begin
    if (pathlen[iseed] lt 0) then continue
    tangents[0:2,pathlen[iseed]-1,iseed] = $
           fieldintegrator_get_field( $
                 lineel[pathlen[iseed]-1,iseed], $
                 path[0:2,pathlen[iseed]-1,iseed] )
  endfor

  result=objarr(nseeds)
  for iseed=0,nseeds-1 do begin
    if (pathlen[iseed] lt 0) then continue
    path[0:2,0:pathlen[iseed]-1,iseed]=path[0:2,0:pathlen[iseed]-1,iseed]*spread(spread(spread(dxmin,0,3),1,pathlen[iseed]),2,1)+spread(spread(xyz0,1,pathlen[iseed]),2,1)
    result[iseed]=obj_new('tony_polyline',npoints=pathlen[iseed],points=transpose(path[0:2,0:pathlen[iseed]-1,iseed]))
  endfor
  return,result
end
