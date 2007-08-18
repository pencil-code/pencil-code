;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   multigrid.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   09-Nov-2005
;;;
;;;  Description:
;;;    Multigrid method for 3-d stationary heat-conduction problem
;;;      Δ T = -q(x) , x ∈ V
;;;      T(∂V) = 0
;;;  See also:
;;;    ${PENCIL_HOME}/doc/implementation/multigrid/idl/multigrid_1d/multigrid.pro
;;;    for a 1-d implementation that comes with notes and discussion of
;;;    different refinement schemes.

function q, xx, yy, zz, $
            CENTRE=cent
  ;; Heating profile -- an elliptical 3-d Gaussian
  w = [0.05, 0.1, 0.05] ;/ 10.
  w = [1,1,1]*0.05
  cent = [2./3., 0.5, 0.3]
  rad2 = ((xx-cent[0])/w[0])^2 + ((yy-cent[1])/w[1])^2 + ((zz-cent[2])/w[2])^2
  return, exp(-rad2)/(2*!pi)^1.5/w[0]/w[1]/w[2]

  ; (not really the same as our 1-d example, because of the boundary
  ; conditions:)
  ; return,  exp(-(xx-cent[0])^2/(2.*w[0]^2))/sqrt(2*!pi)/w[0]

end
; ---------------------------------------------------------------------- ;
function laplacian, f, dxyz
;
; 3-d Laplacian
;
  forward_function laplacian_x, laplacian_y, laplacian_z

  return, laplacian_x(f,dxyz) + laplacian_y(f,dxyz) + laplacian_z(f,dxyz)
;
end
; ---------------------------------------------------------------------- ;
function laplacian_x, f, dxyz
;
; Component of 3-d Laplacian
;
  s = size(f)
  Nx=s[1]

  d2f = 0.*f
  if (Nx gt 2) then begin
    d2f[1:Nx-2,*,*] = (     f[0:Nx-3,*,*] $
                        - 2*f[1:Nx-2,*,*] $
                        +   f[2:Nx-1,*,*] ) / dxyz[0]^2
  endif

  return, d2f
;
end
; ---------------------------------------------------------------------- ;
function laplacian_y, f, dxyz
;
; Component of 3-d Laplacian
;
  s = size(f)
  Ny=s[2]

  d2f = 0.*f
  if (Ny gt 2) then begin
    d2f[*,1:Ny-2,*] = (     f[*,0:Ny-3,*] $
                        - 2*f[*,1:Ny-2,*] $
                        +   f[*,2:Ny-1,*] ) / dxyz[1]^2
  endif

  return, d2f
;
end
; ---------------------------------------------------------------------- ;
function laplacian_z, f, dxyz
;
; Component of 3-d Laplacian
;
  s = size(f)
  Nz=s[3]

  d2f = 0.*f
  if (Nz gt 2) then begin
    d2f[*,*,1:Nz-2] = (     f[*,*,0:Nz-3] $
                        - 2*f[*,*,1:Nz-2] $
                        +   f[*,*,2:Nz-1] ) / dxyz[2]^2
  endif

  return, d2f
;
end
; ---------------------------------------------------------------------- ;
pro zero_boundaries, f
;
; Set boundary values to 0
;
  s = size(f)
  Nx=s[1] & Ny=s[2] & Nz=s[3]
  ;
  f[0,*,*]=0. & f[Nx-1,*   ,*   ]=0.
  f[*,0,*]=0. & f[*   ,Ny-1,*   ]=0.
  f[*,*,0]=0. & f[*   ,*   ,Nz-1]=0.
;
end
; ---------------------------------------------------------------------- ;
function residual, f, g, dxyz
;
; Return 
;   r = L f - g
; of linear system
;   L f_exact = g
; we try to solve
;
  s = size(f)
  res = make_array(SIZE=s)
  Nx=s[1] & Ny=s[2] & Nz=s[3]
  dx = dxyz[0]

  ;; do not update bdry values -> remain 0
  la = laplacian(f,dxyz)
  res = la - g
  zero_boundaries, res

  return, res
end
; ---------------------------------------------------------------------- ;
function gauss_seidel, f, g, dxyz
  s = size(f)
  Nx=s[1] & Ny=s[2] & Nz=s[3]
  ;
  dx=dxyz[0] & dy=dxyz[1] & dz=dxyz[2]
  dx_2 = 1./dx^2
  dy_2 = 1./dy^2
  dz_2 = 1./dz^2

  ;; do not update boundary values -> remain 0

;;   for ix=1L,Nx-2 do begin
;;     for iy=1L,Ny-2 do begin
;;       for iz=1L,Nz-2 do begin
;;         f[ix,iy,iz] $
;;             = (   (f[ix-1,iy  ,iz  ] + f[ix+1,iy  ,iz  ]) * dx_2 $
;;                 + (f[ix  ,iy-1,iz  ] + f[ix  ,iy+1,iz  ]) * dy_2 $
;;                 + (f[ix  ,iy  ,iz-1] + f[ix  ,iy  ,iz+1]) * dz_2 $
;;                 - g[ix,iy,iz] $
;;               ) / (2*(dx_2+dy_2+dz_2))
;;       endfor
;;     endfor
;;   endfor
;;   ;
;;   return, f


  ; Be clever about one (or more) of Nx, Ny, Nz = 2. If so, we just
  ; discard the corresponding second derivative from the Laplacian.
  if (Nx le 2) then dx_2 = 0.
  if (Ny le 2) then dy_2 = 0.
  if (Nz le 2) then dz_2 = 0.

  if (any([Nx,Ny,Nz] gt 2)) then begin
    for ix=1L,(Nx-2)>1 do begin
      for iy=1L,(Ny-2)>1 do begin
        for iz=1L,(Nz-2)>1 do begin
          sum = -g[ix,iy,iz]
          if (Nx gt 2) then $
              sum = sum + (f[ix-1,iy  ,iz  ] + f[ix+1,iy  ,iz  ]) * dx_2
          if (Ny gt 2) then $
              sum = sum + (f[ix  ,iy-1,iz  ] + f[ix  ,iy+1,iz  ]) * dy_2
          if (Nz gt 2) then $
              sum = sum + (f[ix  ,iy  ,iz-1] + f[ix  ,iy  ,iz+1]) * dz_2
                                ;
          f[ix,iy,iz] = sum / (2*(dx_2+dy_2+dz_2))
        endfor
      endfor
    endfor
  endif
  ;
  return, f
end
; ---------------------------------------------------------------------- ;
function jacobi, f, g, dxyz
  s = size(f)
  Nx=s[1] & Ny=s[2] & Nz=s[3]
  ;
  dx=dxyz[0] & dy=dxyz[1] & dz=dxyz[2]
  dx_2 = 1./dx^2
  dy_2 = 1./dy^2
  dz_2 = 1./dz^2

  if ( (min([Nx,Ny,Nz]) le 2) and (min([Nx,Ny,Nz]) gt 2)) then begin
    message, 'JACOBI does not handle this case correctly'
    ;; To fix this, a number of nested ifs would be required to get
    ;; all possible orderings of Nx, Ny, Nz right.
  endif

  if ((Nx gt 2) and (Ny gt 2) and (Nz gt 2)) then begin
  ;; do not update boundary values -> remain 0
    f[1:Nx-2,1:Ny-2,1:Nz-2] $
        = (   (f[0:Nx-3,1:Ny-2,1:Nz-2] + f[2:Nx-1,1:Ny-2,1:Nz-2]) * dx_2 $
            + (f[1:Nx-2,0:Ny-3,1:Nz-2] + f[1:Nx-2,2:Ny-1,1:Nz-2]) * dy_2 $
            + (f[1:Nx-2,1:Ny-2,0:Nz-3] + f[1:Nx-2,1:Ny-2,2:Nz-1]) * dz_2 $
            - g[1:Nx-2,1:Ny-2,1:Nz-2] $
          ) / (2*(dx_2+dy_2+dz_2))
  endif
  ;
  return, f
end
; ---------------------------------------------------------------------- ;
function interpolate_between_grids, f1, xyz1, xyz2
;
; Interpolate f1 trilinearly from grid xyz1 to xyz2.
; Used for interpolating to coarser (in restrict) and to finer grid (in refine)
;
  x0=xyz1[0,0,0,0] & dx=xyz1[1,0,0,0]-x0
  y0=xyz1[0,0,0,1] & dy=xyz1[0,1,0,1]-y0
  z0=xyz1[0,0,0,2] & dz=xyz1[0,0,1,2]-z0
  ;
  ixx = reform(xyz2[*,0,0,0] - x0) / dx
  iyy = reform(xyz2[0,*,0,1] - y0) / dy
  izz = reform(xyz2[0,0,*,2] - z0) / dz
  ;
  f2 = interpolate(f1, ixx, iyy, izz, /GRID)
  ;
  return, f2
end
; ---------------------------------------------------------------------- ;
function restrict, f, xyz, xyz_coarse, $
                   QUIET=quiet
;
; Restrict (= coarse-grain) f (the residual) from the finer to the next
; coarser grid
;
  s = size(f)
  Nx=s[1] & Ny=s[2] & Nz=s[3]
  ;
  s_coarse = size(xyz_coarse)
  Nx_coarse=s_coarse[1] & Ny_coarse=s_coarse[2] & Nz_coarse=s_coarse[3]

  ;; 1. Calculate `full-weighting' (i.e. [1 2 1]-weighted) restricted
  ;;    values on fine grid
  fr1 = f
  smooth_full_weight, fr1

  ;; 2. Interpolate those values trilinearly onto coarse grid points
  fr = interpolate_between_grids(fr1, xyz, xyz_coarse)
  ;; Dirichlet boundary conditions should automatically be preserved

  if (not keyword_set(quiet)) then $
      print_depth, [Nx,Ny,Nz], [Nx_coarse,Ny_coarse,Nz_coarse], /COARSER

  return, fr
end
; ---------------------------------------------------------------------- ;
pro smooth_full_weight, f
;
; Apply `fully weighted' smoothing to F in all three directions
;
  s = size(f)
  Nx=s[1] & Ny=s[2] & Nz=s[3]

  if (Nx gt 2) then f[1:Nx-2,*     ,*     ] = 0.25 * ( f[0:Nx-3,*     ,*     ] + 2*f[1:Nx-2,*     ,*     ] + f[2:Nx-1,*     ,*     ] ) 
  if (Ny gt 2) then f[*     ,1:Ny-2,*     ] = 0.25 * ( f[*     ,0:Ny-3,*     ] + 2*f[*     ,1:Ny-2,*     ] + f[*     ,2:Ny-1,*     ] ) 
  if (Nz gt 2) then f[*     ,*     ,1:Nz-2] = 0.25 * ( f[*     ,*     ,0:Nz-3] + 2*f[*     ,*     ,1:Nz-2] + f[*     ,*     ,2:Nz-1] ) 

end
; ---------------------------------------------------------------------- ;
function refine, f, xyz, xyz_fine
;
; Refine (= fine-grain by interpolating) f from coarser to the next finer
; grid
;
  s = size(f)
  Nx=s[1] & Ny=s[2] & Nz=s[3]
  ;
  s_fine = size(xyz_fine)
  Nx_fine=s_fine[1] & Ny_fine=s_fine[2] & Nz_fine=s_fine[3]

  ;
  ; trilinear interpolation
  ;
  fr = interpolate_between_grids(f, xyz, xyz_fine)

  print_depth, [Nx_fine,Ny_fine,Nz_fine], [Nx,Ny,Nz], /FINER

  return, fr
end
; ---------------------------------------------------------------------- ;
function pack_3d_vectors, xx, yy, zz
;
; Combine 3 arrays xx, yy, zz into one 4-d array
;
  s = size(xx)
  Nx=s[1] & Ny=s[2] & Nz=s[3]

  packed = fltarr(Nx,Ny,Nz,3)
  packed[*,*,*,0] = xx
  packed[*,*,*,1] = yy
  packed[*,*,*,2] = zz

  return, packed
  ;
end
; ---------------------------------------------------------------------- ;
function level_up, f, g, dxyz, n
;
; Move n levels up [towards coarser grid] (if possible) and do one
; Gauss-Seidel step at each coarser level.
; Not sure what n < ∞ is good for (other than testing), as `W' cycles
; would limit the number of levels to go _down_ from the coarsest, not
; _up_ from the finest level.
;
  s = size(f)
  Nx=s[1] & Ny=s[2] & Nz=s[3]

  Nx_2 = (Nx+1) / 2             ; number of points on the coarser grid
  Ny_2 = (Ny+1) / 2             ; ditto
  Nz_2 = (Nz+1) / 2             ; ditto
  ;
  Nx_2 = Nx_2 > 2               ; only >= 2 makes sense
  Ny_2 = Ny_2 > 2
  Nz_2 = Nz_2 > 2

  if ((n gt 0) and any([Nx,Ny,Nz] gt 2)) then begin

    fine   = mkgrid( Nx,      Ny,      Nz, $
                     [0.,1.], [0.,1.], [0.,1.] )
    coarse = mkgrid( Nx_2,    Ny_2,    Nz_2, $
                     [0.,1.], [0.,1.], [0.,1.] )

    ; Bloody hell: concatenation works only up to 3-d arrays
    xyz        = pack_3d_vectors(fine.xx, fine.yy, fine.zz)
    dxyz       = [fine.dx, fine.dy, fine.dz]
    xyz_coarse = pack_3d_vectors(coarse.xx, coarse.yy, coarse.zz)

    r = residual(f, g, dxyz)

    f1 = restrict(f, xyz, xyz_coarse)
    r1 = restrict(r, xyz, xyz_coarse, /QUIET)

    dx1 = fine.dx * (Nx-1.) / (Nx_2-1.) 
    dy1 = fine.dy * (Ny-1.) / (Ny_2-1.) 
    dz1 = fine.dz * (Nz-1.) / (Nz_2-1.) 
    dxyz1 = [dx1, dy1, dz1]

    df1 = level_up(0.*r1, r1, dxyz1, n-1)
    ;
    f = f - refine(df1, xyz_coarse, xyz)

  endif else begin

    ; print, 'LEVEL_UP:     Nothing to do for Nxyz<3: [Nx,Ny,Nz] = ', $
    ;        [Nx,Ny,Nz], $
    ;        FORMAT='(A,3I4)'

  endelse

  f = gauss_seidel(f, g, dxyz)
  ; f = jacobi(f,g,dxyz)
  return, f
end
; ---------------------------------------------------------------------- ;
pro print_depth, nn, mm, COARSER=coarser, FINER=finer
;
; Some diagnostic output that visually indicates the level we are at
;
  default, coarser, 0
  default, finer, 0

  indent = '                              '
  indent = strmid(indent,1,30-4*alog(max(nn))/alog(2))

  fmt = '(A,3I4,A,3I4)'

  if (coarser) then begin
    print, indent, nn, ' --> ', mm, FORMAT=fmt
  endif else if (finer) then begin
    print, indent, nn, ' <-- ', mm, FORMAT=fmt
  endif else begin
    print, indent, nn, ' --- ', mm, FORMAT=fmt
  endelse
;
end
; ---------------------------------------------------------------------- ;

niter = 20
Nx = 20
Ny = 40
Nz = 60

Nx = 64
Ny = 64
Nz = 64

grid = mkgrid( Nx,     Ny,     Nz, $
               [0.,1], [0.,1], [0.,1] )

dxyz = [grid.dx, grid.dy, grid.dz]


;f = 0.01*cos(!pi*x/dx) & f[0]=0. & f[Nx-1]=0.
f = 0.*grid.xx
g = -q(grid.xx, grid.yy, grid.zz, CENTRE=cent) ; heating 

ix0 = (cent[0]-grid.x[0]) / grid.dx + 0.5
iy0 = (cent[1]-grid.y[0]) / grid.dy + 0.5
iz0 = (cent[2]-grid.z[0]) / grid.dz + 0.5


save_state

plot, grid.x, f[*,iy0,iz0], $
      YRANGE=[0,1]*1, YSTYLE=3, $
      TITLE='!8N!6='+strtrim(Nx,2)+'!X'

colors = [-1]                   ; initialize, so we can append below
labels = ['ignore']

for i=0,niter-1 do begin
;   f = jacobi(f, g, dxyz)
;  f = gauss_seidel(f, g, dxyz)
  f = level_up(f, g, dxyz, 1000)
  color = i*150./niter
  oplot, grid.x, f[*,iy0,iz0], $
         COLOR=color, $
         PSYM=psym
  colors = [colors, color]
  labels = [labels, '!8i!6='+strtrim(i,2)+'!X']
  ;; Store f values
  if (i eq 0) then ff = [[f]] else ff = [[ff],[f]]

  print, strtrim(i,2), '/',strtrim(niter,2), ' : ', max(f[*,iy0,iz0])

endfor


ophline, [0.315802, 0.246309, 0.4223]
  esrg_legend, labels[1:*], $
               COLOR=colors[1:*], $
               /BOX, POS=[0.45, 0.1, 0.7, 0.4], $
               XSTRETCH=0.5


restore_state

la=laplacian(f,dxyz)
print, 'RMS error (abs, rel):', rms(la-g), rms(la-g)/rms(g)

end
; End of file multigrid.pro
