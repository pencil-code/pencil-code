;;
;; Scalar product of two vector fields
;; Works for 1d, 2d or 3d fields
;;

function dot, f, g
  ;; Size and consistency checks
  sf = size(f)
  sg = size(g)
  if (sf[0] ne sg[0]) then message, 'Incompatible shapes'
  if ((sf[sf[0]] ne 3) or (sg[sg[0]] ne 3)) then message, 'Need two 3-vectors'

  case (sf[0]) of
    1: return, g[      0]*f[      0]+g[      1]*f[      1]+g[      2]*f[      2]
    2: return, g[*,    0]*f[*,    0]+g[*,    1]*f[*,    1]+g[*,    2]*f[*,    2]
    3: return, g[*,*,  0]*f[*,*,  0]+g[*,*,  1]*f[*,*,  1]+g[*,*,  2]*f[*,*,  2]
    4: return, g[*,*,*,0]*f[*,*,*,0]+g[*,*,*,1]*f[*,*,*,1]+g[*,*,*,2]*f[*,*,*,2]
    else: message, "Don't know how to treat fields of size " + strim(sf,2)
  endcase

end
