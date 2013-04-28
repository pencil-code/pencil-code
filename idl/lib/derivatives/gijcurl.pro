;;
;;  $Id$
;;
;;  Derivative matrix from curl in cartesian coordinates.
;;  Used to calculate B_{i,j} as bij=gijcurl(aa).
;;
;;  20-mar-04/axel: adapted from f90 routine bij_etc
;;
function gijcurl,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  if (coord_system ne 'cartesian') then message, $
      "gijcurl not yet implemented for coord_system='" + coord_system + "'"
;
  s = size(f)
  if (s[0] ne 4) then message, "gijcurl is only implemented for 4D arrays."
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  w = make_array(size=[5,fmx,fmy,fmz,3,3,s[5],fmx*fmy*fmz*3*3])
;
;  Calculate B_i,j = eps_ikl A_l,jk
;
  for i=0,2 do for j=0,2 do for k=0,2 do for l=0,2 do begin
    eps=levi_civita(i,k,l)
print, i, j, k, l, eps
    if (eps ne 0.0) then $
        w[*,*,*,i,j] += eps * derij_single(f[*,*,*,l],j,k)
  endfor
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
