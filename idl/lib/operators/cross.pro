;
;  $Id$
;
function cross, g, f
;
  s = size (f)
  if (any (s[1:s[0]] ne size (g, /dimensions))) then $
      message, 'cross: vectors must be of the same size'
  if ((s[s[0]] ne 3) and (s[0] ge 2) and (s[1] eq 3)) then $
      return, cross (transpose (g), transpose (f))
  if (s[s[0]] ne 3) then $
      message, 'cross: vectors must be of size 3 in last dimension'
;
  w=make_array(size=s)
;
  if (s[0] eq 2) then begin
    w[*,0] = g[*,1]*f[*,2] - g[*,2]*f[*,1]
    w[*,1] = g[*,2]*f[*,0] - g[*,0]*f[*,2]
    w[*,2] = g[*,0]*f[*,1] - g[*,1]*f[*,0]
  endif else if (s[0] eq 3) then begin
    w[*,*,0] = g[*,*,1]*f[*,*,2] - g[*,*,2]*f[*,*,1]
    w[*,*,1] = g[*,*,2]*f[*,*,0] - g[*,*,0]*f[*,*,2]
    w[*,*,2] = g[*,*,0]*f[*,*,1] - g[*,*,1]*f[*,*,0]
  endif else if (s[0] eq 4) then begin
    w[*,*,*,0] = g[*,*,*,1]*f[*,*,*,2] - g[*,*,*,2]*f[*,*,*,1]
    w[*,*,*,1] = g[*,*,*,2]*f[*,*,*,0] - g[*,*,*,0]*f[*,*,*,2]
    w[*,*,*,2] = g[*,*,*,0]*f[*,*,*,1] - g[*,*,*,1]*f[*,*,*,0]
  endif else begin
    message, "cross: can't cross vectors with "+strtim (s[0], 2)+" dimensions."
  endelse
;
  return, w
;
end
