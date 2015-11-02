;
; IDL routine 
; Copyright (c) 2012 Oliver Gressel
;
; adopted from Python source:
;
; Streamline plotting like Mathematica.
; Copyright (c) 2011 Tom Flannaghan.

;
; call structure (same as IDL's "velovect")
; optional keyword 'density' specifies sampling density in [x,y]
;
; IDL> streamplot, u, v, x, y, density=density, _extra=extr
;

; Permission is hereby granted, free of charge, to any person obtaining a copy
; of this software and associated documentation files (the "Software"), to deal
; in the Software without restriction, including without limitation the rights
; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
; copies of the Software, and to permit persons to whom the Software is
; furnished to do so, subject to the following conditions:

; The above copyright notice and this permission notice shall be included in
; all copies or substantial portions of the Software.

; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
; THE SOFTWARE.

; ------------------------------------------------------------------------------

function blank_pos, xi, yi

  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby

  ; --- takes grid space coords and 
  ;     returns nearest space in the blank array
  return, [ long((xi / dbx) + 0.5), $
            long((yi / dby) + 0.5)  ]
end

; ------------------------------------------------------------------------------

function value_at, a, xi, yi, oob=oob, tstep=tstep

  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby

  ; --- linear interpolation - nice and quick because we are
  ;     working in grid-index coordinates.
  ;     oob: out-of-bounds flag

  if (xi lt 0) or (xi gt ngx-1) or $
     (yi lt 0) or (yi gt ngy-1) then oob=1 else oob=0

  x = long(xi > 0) &  xp = long(xi+1 < (ngx-1))
  y = long(yi > 0) &  yp = long(yi+1 < (ngy-1))

  a00 = a[y,  x ]
  a01 = a[y,  xp]
  a10 = a[yp, x ]
  a11 = a[yp, xp]

  xt = xi - x
  yt = yi - y
  a0 = a00*(1-xt) + a01*xt
  a1 = a10*(1-xt) + a11*xt

  res = a0*(1-yt) + a1*yt

  ; --- avoid zero return value if provided for setting of stepsize

  if keyword_set(tstep) and (res eq 0.) then begin
    if a01 ne 0. then begin
      xi=xp & res=a01 
    endif else if a10 ne 0. then begin
      yi=yp & res=a10 
    endif else if a11 ne 0. then begin
      xi=xp & yi=yp & res=a11
    endif
  endif

  return, res
end

; ------------------------------------------------------------------------------
        
function f, xi, yi, bw=bw, oob=oob

  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby
  common fields, u,v,x,y, speed, blank

  dt_ds = value_at(speed, xi, yi, oob=oob, /tstep)
  if dt_ds  eq 0 then stop,'streamplot: Interpolated value of speed is zero'
  dt_ds = 1./dt_ds

  ui = value_at(u, xi, yi)
  vi = value_at(v, xi, yi)

  dir = keyword_set(bw) ? -1. : +1. ; backward / forward
  return, [ui*dt_ds, vi*dt_ds] *dir
end

; ------------------------------------------------------------------------------

function check_bounds, xi, yi

  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby

  return, (xi ge 0) and (xi lt ngx-1) $
      and (yi ge 0) and (yi lt ngy-1)
end

; ------------------------------------------------------------------------------

function my_rk4, x0, y0, ytr=ytr, stotal=stotal, backward=bw

  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby
  common fields, u,v,x,y, speed, blank
  common changes, bx_changes, by_changes

  ds = min([1./ngx, 1./ngy, 0.005])
  stotal = 0
  xi = x0 &  xl = x0+ds
  yi = y0 &  yl = y0+ds

  bl = blank_pos(xi, yi) &  xb=bl[0] &  yb=bl[1]
  xf_traj = []
  yf_traj = []

  ; --- time step loop

  while check_bounds(xi, yi) do begin

     ; --- save point if incremental step 
     ;     is above round-off error 

     dl2 = (xi-xl)^2 + (yi-yl)^2
     if (dl2 gt 0.) then begin
        xf_traj = [ xf_traj, xi ] &  xl=xi
        yf_traj = [ yf_traj, yi ] &  yl=yi
     endif

     ; --- next, advance one using RK4

     k1 = f(xi,               yi              , bw=bw, oob=oob)
     k2 = f(xi + .5*ds*k1[0], yi + .5*ds*k1[1], bw=bw, oob=oob)
     k3 = f(xi + .5*ds*k2[0], yi + .5*ds*k2[1], bw=bw, oob=oob)
     k4 = f(xi +    ds*k3[0], yi +    ds*k3[1], bw=bw, oob=oob)

     ; --- check for out-of-bounds
     if (oob eq 1) then break

     ; --- update position
     xi += ds*(k1[0] + 2*k2[0] + 2*k3[0]+k4[0]) / 6.
     yi += ds*(k1[1] + 2*k2[1] + 2*k3[1]+k4[1]) / 6.

     ; --- final position might be out of the domain
     if not check_bounds(xi, yi) then break
     stotal += ds

     ; --- next, if s gets to thres, check blank.
     bl = blank_pos(xi, yi) &  new_xb=bl[0] &  new_yb=bl[1]

     if (new_xb ne xb) or $
        (new_yb ne yb) then begin

        ; --- new square, so check and colour. quit if required.
        if blank[new_yb,new_xb] eq 0 then begin
           blank[new_yb,new_xb] = 1

           bx_changes = [ bx_changes, new_xb ]
           by_changes = [ by_changes, new_yb ]
           xb = new_xb
           yb = new_yb
        endif else break
     endif

     ; --- avoid too long lines
     if (stotal gt 2.) then break
  end
  
  ytr = yf_traj
  return, xf_traj
end

; ------------------------------------------------------------------------------

function flow_integrate, x0, y0

  common fields, u,v,x,y, speed, blank
  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby
  common changes, bx_changes, by_changes

  ; --- this function does RK4 forward and back trajectories from
  ;     the initial conditions, with the odd 'blank array' termination

  bx_changes = []
  by_changes = []

  xf_traj = my_rk4(x0, y0, ytr=yf_traj, stotal=sf)
  xb_traj = my_rk4(x0, y0, ytr=yb_traj, stotal=sb, /backward)

  stotal = sf + sb
  x_traj = [ reverse(xb_traj), xf_traj ]
  y_traj = [ reverse(yb_traj), yf_traj ]

  ; --- tests to check length of traj. 
  ;     remember: s in units of axes.

  if (n_elements(x_traj) lt 1) then return, -1
  
  nbl = n_elements(bx_changes)

  ; --- only return long-enough line segments
  ;     which cover more than two pixels (avoids small doodle)

  if (stotal gt 0.1) and (nbl gt 2) then begin  

     bl = blank_pos(x0, y0) &  blank[bl[1],bl[0]] = 1
     return, [[x_traj], [y_traj]]

  endif else begin  ; --- otherwise erase the line segment's tracks

     for i=0,nbl-1 do $
        blank[ by_changes[i], bx_changes[i] ] = 0
     return, -1

  endelse

end

; ------------------------------------------------------------------------------

pro traj, xb,yb, _extra=extr, arrow=arrow

  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby
  common fields, u,v,x,y, speed, blank

  ; --- a quick function for integrating trajectories if blank==0.

  if (xb lt 0) or (xb ge nbx) or $
     (yb lt 0) or (yb ge nby) then return

  if (blank[yb, xb] eq 0) then begin

     ; --- do the integration
     tr = flow_integrate(xb*dbx, yb*dby)

     if ((size(tr))[0] gt 0) then begin

        tx = x[0] + tr[*,0]*dx
        ty = y[0] + tr[*,1]*dy

        ; plot line segment
        oplot, tx, ty, _extra=extr
        
        ; add arrow
        if keyword_set(arrow) then begin
          nn = n_elements(tx)/2
          if (nn gt 8) then $ ; jitter arrow position
             nn += long(0.16*nn*(randomn(system_seed)-0.5))
          arrow_pc, tx[nn-2], ty[nn-2], tx[nn], ty[nn], $
                    /data, _extra=extr
        endif
     endif
  endif

  return
end

; ------------------------------------------------------------------------------

pro streamplot, uu, vv, xx, yy, density=density, _extra=extr, hsize=hsize

  common geomet, ngx,ngy, dx,dy, nbx,nby, dbx,dby
  common fields, u,v,x,y, speed, blank

  ;  - x and y are 1d arrays defining an *evenly spaced* grid.
  ;  - u and v are 2d arrays (shape [y,x]) giving velocities.
  ;  - density controls the closeness of the streamlines. For different
  ;    densities in each direction, use a tuple [density_x, density_y].

  if not keyword_set(density) then density = 1.
  if ((size(density))[0] eq 0) then density=[density,density]

  ; --- set defaults for arrow lengths

  if not keyword_set(hsize) then hsize = (!d.name eq 'X') ? 8 : 160

  ; --- grid dimensions

  ngx = n_elements(xx)
  ngy = n_elements(yy)

  ; --- in polar coordinates: r d\phi/dt = v_phi !

  if keyword_set(extr.polar) then $
    for i=0,ngy-1 do vv(*,i) /= xx  

  ; --- copy input variables to common block

  u = transpose(uu)
  v = transpose(vv)
  x = xx
  y = yy

  ; --- grid properties

  dx = x[1] - x[0]
  dy = y[1] - y[0]

  ; --- sanity checks

  if( (size(u))[0] ne 2 or $
      (size(v))[0] ne 2 ) then message, 'fields must be 2d arrays'
  if( (size(x))[0] ne 1 or $
      (size(y))[0] ne 1 ) then message, 'coordinates must be 1d arrays'

  if( (size(u))[1] ne ngy or $
      (size(u))[2] ne ngx ) then message, "incompatible dimensions for 'u'"

  if( (size(v))[1] ne ngy or $
      (size(v))[2] ne ngx ) then message, "incompatible dimensions for 'v'"

  ; --- rescale velocity onto axes-coordinates
  ;     path length, s, will now be in axes-coordinate

  u /= ( x[ngx-1] - x[0] )
  v /= ( y[ngy-1] - y[0] )
  speed = sqrt( u^2 + v^2 )

  ; --- convert u and v into grid-coordinates

  u *= ngx
  v *= ngy

  ; --- Blank array: This is the heart of the algorithm. It begins life
  ;     zeroed, but is set to one when a streamline passes through each
  ;     box. Then streamlines are only allowed to pass through zeroed
  ;     boxes. The lower resolution of this grid determines the
  ;     approximate spacing between trajectories.

  nbx = byte(30*density[0])
  nby = byte(30*density[1])
  blank = bytarr([nby,nbx])

  ; --- Constants for conversion between grid-index space 
  ;     and blank-index space

  dbx = ngx/float(nbx-1)
  dby = ngy/float(nby-1)

  ; --- Now we build up the trajectory set. I've found it best to look
  ;     for blank==0 along the edges first, and work inwards.

  for indent = 0,max([nbx,nby])/2-1 do $
     for xi = 0,max([nbx,nby])-2*indent-1 do begin

        traj, xi + indent,  indent,       _extra=extr, hsize=hsize
        traj, xi + indent,  nby-1-indent, _extra=extr, hsize=hsize
        traj, indent,       xi + indent,  _extra=extr, hsize=hsize
        traj, nbx-1-indent, xi + indent,  _extra=extr, hsize=hsize

     endfor

  return
end

; ------------------------------------------------------------------------------

pro streamtest, dens=dens, vec=vec, _ref_extra=extr

  ; --- sample test program

  nx = 60
  ny = 40

  x = 2*findgen(nx)/(nx-1) - 1
  y = 2*findgen(ny)/(ny-1) - 1

  u = -abs(x) # abs(y)
  v =  abs(x) # y^2

  plot, [0], [0], xr=[-1,1], yr=[-1,1], xsty=1, ysty=1, /nodata, /iso

  if keyword_set(vec) then begin

     ; --- normalise vector length
     vel = sqrt(u^2+v^2)
     uu = u/vel
     vv = v/vel

     nn=4
     velovect, uu[0:nx-1:nn,0:ny-1:nn], $
               vv[0:nx-1:nn,0:ny-1:nn], /overplot, $
               x[0:nx-1:nn], y[0:ny-1:nn], len=1.5, th=1
  endif
  
  ; --- plot stream lines
  streamplot, u,v,x,y, dens=dens, _extra=extr

  return
end

; ------------------------------------------------------------------------------
