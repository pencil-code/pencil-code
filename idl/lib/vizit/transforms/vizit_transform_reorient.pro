pro vizit_transform_reorient,surface=surface,upVector=upVector,direction=direction,scale=scale,position=position,recenter=recenter,matrix=tmat_in

 ; default,upVector,[0.,1.,0.]
  if keyword_set(direction) then begin
    v=direction/sqrt(dot2(direction))
    w=v 
    w[1]=0.

    ;theta = ATAN(v[1], v[0])
    theta = ATAN(v[1], v[0])
    ;phi   = ATAN(sqrt(dot2(v-w)),sqrt(v[2]^2+v[0]^2))
    phi   = ACOS(v[2]/sqrt(v[2]^2+v[1]^2+v[0]^2))
  endif else begin
    theta=0.
    phi=0.
  endelse

  if keyword_set(upVector) then begin
    unitup=upVector/dot2(upVector)
    T3D,/reset,MATRIX=t3d_twist,rotate=[0.,-theta,0.]
    T3D,t3d_twist,MATRIX=t3d_twist,rotate=[0.,0.,-phi]
    twist=VERT_T3D(unitup,MATRIX=t3d_twist)
    omega=-ATAN(twist[1], twist[2])
  endif else begin
    omega=0.
  endelse

  T3D,/reset,MATRIX=tmat
  if keyword_set(tmat_in) then begin
    tmat=tmat_in
  endif


omega=180.*omega/!pi
theta=180.*theta/!pi
phi=(180.*phi/!pi)
;print,omega,theta,phi
  if keyword_set(recenter) then begin
    T3D,tmat,MATRIX=tmat_out,translate=-recenter
    tmat=tmat_out
  endif

  if keyword_set(scale) then begin
    T3D,tmat,MATRIX=tmat_out,scale=scale
    tmat=tmat_out
  endif

  T3D,tmat,MATRIX=tmat_out,rotate=[omega,phi-90.,theta]
  tmat=tmat_out

;  T3D,tmat,MATRIX=tmat_out,rotate=[0.,theta,0.]
;  tmat=tmat_out

  if keyword_set(position) then begin
    T3D,tmat,MATRIX=tmat_out,translate=position
    tmat=tmat_out
  endif

  if keyword_set(surface) then begin
    new_vertices=VERT_T3D(*(surface.vertices),MATRIX=tmat,/NO_COPY)  
    ptr_free,surface.vertices
    surface.vertices=ptr_new(new_vertices)
  endif
end
