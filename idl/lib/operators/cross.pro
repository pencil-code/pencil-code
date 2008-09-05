;
;  $Id$
;
function cross, g, f
;
s=size(f)
if (max(abs(size(f)-size(g))) ne 0) then begin
  print, 'cross: vectors must be of the same size'
  stop
endif
if (s[s[0]] ne 3) then begin
  print, 'cross: vectors must be of size 3 in last dimension'
  stop
endif
;
w=make_array(size=s)
;
if (s[0] eq 2) then begin
  w[*,0]=g[*,1]*f[*,2]-g[*,2]*f[*,1]
  w[*,1]=g[*,2]*f[*,0]-g[*,0]*f[*,2]
  w[*,2]=g[*,0]*f[*,1]-g[*,1]*f[*,0]
endif else if (s[0] eq 4) then begin
  w[*,*,*,0]=g[*,*,*,1]*f[*,*,*,2]-g[*,*,*,2]*f[*,*,*,1]
  w[*,*,*,1]=g[*,*,*,2]*f[*,*,*,0]-g[*,*,*,0]*f[*,*,*,2]
  w[*,*,*,2]=g[*,*,*,0]*f[*,*,*,1]-g[*,*,*,1]*f[*,*,*,0]
endif else begin
  print, 'cross: do not know how to cross vectors of size=', s
endelse
;
return,w
;
end
