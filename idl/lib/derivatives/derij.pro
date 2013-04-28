;;
;;  $Id$
;;
;;  Calculate second derivative matrix f_l,ij.
;;
function derij,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'error: derij not implemented for '+strtrim(s[0],2)+'-D arrays'
  fmx = s[1] & fmy = s[2] & fmz = s[3]
;
  if (s[0] eq 3) then begin
;
    w = make_array(size=[5,fmx,fmy,fmz,3,3,s[4],fmx*fmy*fmz*3*3])
    w[*,*,*,0,0] = xder2(f)
    w[*,*,*,1,1] = yder2(f)
    w[*,*,*,2,2] = zder2(f)
    w[*,*,*,0,1] = xderyder(f)
    w[*,*,*,0,2] = xderzder(f)
    w[*,*,*,1,2] = yderzder(f)
    w[*,*,*,1,0] = w[*,*,*,0,1]
    w[*,*,*,2,0] = w[*,*,*,0,2]
    w[*,*,*,2,1] = w[*,*,*,1,2]
;
  endif else begin
;
    w = make_array(size=[6,fmx,fmy,fmz,s[4],3,3,s[5],fmx*fmy*fmz*s[4]*3*3])
    w[*,*,*,*,0,0] = xder2(f)
    w[*,*,*,*,1,1] = yder2(f)
    w[*,*,*,*,2,2] = zder2(f)
    w[*,*,*,*,0,1] = xderyder(f)
    w[*,*,*,*,0,2] = xderzder(f)
    w[*,*,*,*,1,2] = yderzder(f)
    w[*,*,*,*,1,0] = w[*,*,*,*,0,1]
    w[*,*,*,*,2,0] = w[*,*,*,*,0,2]
    w[*,*,*,*,2,1] = w[*,*,*,*,1,2]
;
  endelse
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
