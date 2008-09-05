;;
;;  $Id$
;;
;;  Calculate derivative matrix.
;;
function gij,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  if (coord_system ne 'cartesian') then $
      message, $
        "gij not yet implemented for coord_system='" + coord_system + "'"
;
  s=size(f)
;
  if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,3)
    w[*,*,*,0,0]=xder(f[*,*,*,0])
    w[*,*,*,0,1]=yder(f[*,*,*,0])
    w[*,*,*,0,2]=zder(f[*,*,*,0])
    w[*,*,*,1,0]=xder(f[*,*,*,1])
    w[*,*,*,1,1]=yder(f[*,*,*,1])
    w[*,*,*,1,2]=zder(f[*,*,*,1])
    w[*,*,*,2,0]=xder(f[*,*,*,2])
    w[*,*,*,2,1]=yder(f[*,*,*,2])
    w[*,*,*,2,2]=zder(f[*,*,*,2])
;
  endif else begin
    print, 'error: gij only implemented for 4-D arrays'
  endelse
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
