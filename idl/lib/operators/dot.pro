;;
;; Scalar product of two vector fields
;; Works for 1d, 2d or 3d fields
;;

function dot, f, g
  ;; Size and consistency checks
  s = size (f)
  if (any (s[1:s[0]] ne size (g, /dimensions))) then $
      message, 'dot: vectors must be of the same size'
  if ((s[s[0]] ne 3) and (s[0] ge 2) and (s[1] eq 3)) then $
      return, dot (transpose (f), transpose (g))
  if (s[s[0]] ne 3) then $
      message, 'cross: vectors must be of size 3 in last dimension'

  case (s[0]) of
    1: return, g[      0]*f[      0] + g[      1]*f[      1] + g[      2]*f[      2]
    2: return, g[*,    0]*f[*,    0] + g[*,    1]*f[*,    1] + g[*,    2]*f[*,    2]
    3: return, g[*,*,  0]*f[*,*,  0] + g[*,*,  1]*f[*,*,  1] + g[*,*,  2]*f[*,*,  2]
    4: return, g[*,*,*,0]*f[*,*,*,0] + g[*,*,*,1]*f[*,*,*,1] + g[*,*,*,2]*f[*,*,*,2]
    else: message, "dot: can't handle fields with " + strtrim (s, 2) + " dimensions."
  endcase

end
