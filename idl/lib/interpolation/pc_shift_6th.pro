;$Id$
function pc_shift_6th, f, x, deltax
;
;  20-jun-2008/anders: coded
;
;  Shift 1-D array by 6th order interpolation, similar to what is done in the
;  Pencil Code for shearing boundary conditions.
;
;  Syntax (example):
;     plot,zzz,pc_shift_6th(xyaver.bxmz(*,10),zzz,-c*t(it))
;
g=f

nx=n_elements(g)

dx=x[1]-x[0]
deltax_int=fix(deltax/dx)

g=shift(g,deltax_int)

deltax_fra=deltax/dx-fix(deltax/dx)

g_ghost= make_array(nx+6, type=size(g,/type))
l1=3 & l2=l1+nx-1
g_ghost[l1:l2]=g
g_ghost[0:2]=g_ghost[l2-2:l2]
g_ghost[l2+1:l2+3]=g_ghost[l1:l1+2]

frak=abs(deltax_fra)
c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
if (deltax_fra lt 0.0) then begin
  for ix=l1,l2 do begin
    g[ix-3]=c1*g_ghost[ix-2] + $
            c2*g_ghost[ix-1] + $
            c3*g_ghost[ix]   + $
            c4*g_ghost[ix+1] + $
            c5*g_ghost[ix+2] + $
            c6*g_ghost[ix+3]
  endfor
endif else begin
  for ix=l1,l2 do begin
    g[ix-3]=c1*g_ghost[ix+2] + $
            c2*g_ghost[ix+1] + $
            c3*g_ghost[ix]   + $
            c4*g_ghost[ix-1] + $
            c5*g_ghost[ix-2] + $
            c6*g_ghost[ix-3]
  endfor
endelse

return, g

end
