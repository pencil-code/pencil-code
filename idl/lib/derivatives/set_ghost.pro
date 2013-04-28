;***********************************************************************
pro set_ghost,f,debug=debug,add=add
;
;  set ghost zones
;  currently, only periodic boundary conditions are prepared
;  This procedure is probably obsolete now that I've checked in pc_setghost.
;  I keep this in case it proves advantageous with big data sets.
;
  if keyword_set(add) then begin
    s = size(f)
    omx = s[1]+2*nghostx
    omy = s[2]+2*nghosty
    omz = s[3]+2*nghostz
;
    if (s[0] eq 4) then begin
      w = make_array(size=[3,omx,omy,omz,s[4],s[5],omx*omy*omz*s[4]])
    endif else if (s[0] eq 4) then begin
      w = make_array(size=[4,omx,omy,omz,s[4],omx*omy*omz])
    endif else begin
      message, "set_ghost is only implemented for 3D or 4D arrays."
    endelse
;
    w[nghostx:omx-nghostx-1,nghosty:omy-nghosty-1,nghostz:omz-nghostz-1,*] = f
    f = w
    w = 0
  endif
;
  s = size(f)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
;  debug output
;
  if keyword_set(debug) then begin
    print,'$Id$'
    help,f
    print,'fmx,fmy,fmz=',fmx,fmy,fmz
  endif
;
;  set ghost zones on the left end
;
  f(0:l1-1,*,*,*) = f(fmx-2*nghostx:l2,*,*,*)
  f(*,0:m1-1,*,*) = f(*,fmy-2*nghosty:m2,*,*)
  f(*,*,0:n1-1,*) = f(*,*,fmz-2*nghostz:n2,*)
;
;  set ghost zones on the right end
;
  f(l2+1:fmx-1,*,*,*) = f(l1:2*nghostx-1,*,*,*)
  f(*,m2+1:fmy-1,*,*) = f(*,m1:2*nghosty-1,*,*)
  f(*,*,n2+1:fmz-1,*) = f(*,*,n1:2*nghostz-1,*)
;
end
