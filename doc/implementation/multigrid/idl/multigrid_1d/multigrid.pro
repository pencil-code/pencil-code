;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   multigrid.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   09-Nov-2005
;;;
;;;  Description:
;;;    Multigrid method for simple stationary heat-conduction problem
;;;      T'' = -q(x) , T(0)=T(1)=0
;;;  Note:
;;;   If we initialize f with a Nyquist signal, Jacobi iteration will not
;;;   get rid of it (it will oscillate from iteration to iteration), while
;;;   Gauss-Seidel will do so slowly.
;;;     This probably implies that one should use a low-pass filter (like
;;;   hyperdiffusion) after each iteration. 
;;;  Optimal refinement and restriction methods:
;;;    - For cell-centred multigrid, restrict1 (q*[1 2 1 0] + p*[0 1 2 1])
;;;      leads to better convergence, and a smoother error profile, than
;;;      restrict2 (4-point, Nyquist filtering Bessel interpolation),
;;;      although both are usable.
;;;    - For vertex-centred multigrid, restrict1 (which in this case is
;;;      identical to ``full weighting'', i.e. [1 2 1]) leads to
;;;      dramatically faster convergence (~1/2 the number of iterations
;;;      needed) than for the cell-centred case with one point less.
;;;    - For cell-centred multigrid, spline interpolation leads to
;;;      slightly improved convergence compared to linear interpolation.
;;;      Quadratic or `least-squares quadratic' interpolation don't
;;;      improve things at all.

@restrict

function q, x
  ;; Heating profile
  w = 0.05
  x0 = 2./3.
  return, exp(-(x-x0)^2/(2.*w^2))/sqrt(2*!pi)/w
end
; ---------------------------------------------------------------------- ;
function residual, f, g, dx
  ;; Return residual of linear system we try to solve
  res = 0.*f
  Nx = n_elements(f)             ; number of points
  ; for i=1L,Nx-2 do begin         ; do not update bdry values -> remain 0
  ;   res[i] = (f[i-1]-2*f[i]+f[i+1])/dx^2 - g[i]
  ; endfor
  res[1:Nx-2] = ( f[0:Nx-3] - 2.*f[1:Nx-2] + f[2:Nx-1] ) / dx^2 $
                - g[1:Nx-2]
  ;
  return, res
end
; ---------------------------------------------------------------------- ;
function gauss_seidel, f, g, dx
  Nx = n_elements(f)
  for i=1L,Nx-2 do begin         ; do not update bdry values -> remain 0
    f[i] = (f[i-1]+f[i+1] - dx^2*g[i])/2.
  endfor
  return, f
end
; ---------------------------------------------------------------------- ;
function jacobi, f, g, dx
  Nx = n_elements(f)
  f = (shift(f,1) + shift(f,-1) - dx^2*g)/2.
  f[0]    = 0.
  f[Nx-1] = 0.
  return, f
end
; ---------------------------------------------------------------------- ;
function restrict, f, x, x_coarse, $
                   QUIET=quiet
;
; Restrict (= coarse-grain) f (the residual) from the finer to the next
; coarser grid
;
;  fr = fltarr(N_coarse) * f[0] * 0. ; inherit float type from f
;  for i=1L,N_coarse-2 do begin      ; implies bundary conditions f=0
;    i_fine = 2*i                    ; integer in this case
;    fr[i] = (f[i_fine-1] + 2.*f[i_fine] + f[i_fine+1])/4.
;  endfor

  N_fine   = n_elements(x)
  N_coarse = n_elements(x_coarse)

print, '==> RESTRICT:  N_fine, N_coarse = ', N_fine, N_coarse

  fr = fltarr(N_coarse) * f[0] * 0. ; inherit float type from f
  for i=1L,N_coarse-2 do begin      ; implies bundary conditions f=0
    i_fine = i*1.*(N_fine-1)/(N_coarse-1) ; generally a non-integer
    ;; Use  [1, 2, 1]/4 restriction stencil at the two closest indices...
    i_l  = floor(i_fine)
    i_r  = i_l + 1
    fr_l = (f[i_l-1] + 2.*f[i_l] + f[i_l+1])/4.
    fr_r = (f[i_r-1] + 2.*f[i_r] + f[i_r+1])/4.
    ;; ...and interpolate linearly
    p = i_fine - i_l
    q = 1. - p
    fr[i] = p*fr_r + q*fr_l

print, '    ==> RESTRICT:  i_fine, i_l, i_r = ', i_fine, i_l, i_r
print, '    ==> RESTRICT:  p, q = ', p, q

  endfor

  if (not keyword_set(quiet)) then $
      print_depth, N_fine, N_coarse

  return, fr
end
; ---------------------------------------------------------------------- ;
function refine, f, x, x_fine
;
; Refine (= fine-grain by interpolating) f from coarser to the next finer
; grid
;
  N_coarse = n_elements(x)
  N_fine   = n_elements(x_fine)

  ;
  ; linear interpolation
  ;
  fr = interpol(f,N_fine)

;   ;
;   ; spline interpolation
;   ;
;    if (N_coarse gt 2) then begin
;      fr = spline(x, f, x_fine)
;    endif else begin
;      fr = interpol(f,N_fine)
;    endelse

  ;
  ; quadratic, or least-squares quadratic interpolation
  ;
;    if (N_coarse gt 4) then begin
;      fr = interpol(f,x,x_fine,/QUADRATIC)
; ;     fr = interpol(f,x,x_fine,/LSQUADRATIC)
;    endif else begin
;      fr = interpol(f,N_fine)
;    endelse

  print_depth, N_coarse, N_fine

  return, fr
end
; ---------------------------------------------------------------------- ;
function level_up, f, g, dx, n
  ;; Move n levels up [towards coarser grid] (if possible) and do one
  ;; Gauss-Seidel step there

  Nx   = n_elements(f)
  Nx_2 = (Nx+1) / 2             ; number of points on the coarser grid

  if ((n gt 0) and (n_elements(f) ge 3)) then begin

    if ((Nx_2)*2 eq Nx) then begin
      ;; Even number of points (= odd number of intervals -- cell-centred
      ;; multigrid).
      ;; New grid points don't (generally) coincide with old ones:
      ;; +---+---+---+---+---+---+---+
      ;; |        |         |        |
      ;; +--------+---------+--------+
      ;;
      odd = 0
    endif else begin
      ;; Odd number of points (= even number of intervals --
      ;; vertex-centred multigrid).
      ;; New grid is old grid with every other point removed:
      ;; +---+---+---+---+---+---+---+---+
      ;; |       |       |       |       |
      ;; +-------+-------+-------+-------+
      ;;
      odd = -1
    endelse

    x        = linspace([0.,1.],Nx)
    x_coarse = linspace([0.,1.],Nx_2)

    r = residual(f,g,dx)

    f1 = restrict2(f, x, x_coarse)
    r1 = restrict2(r, x, x_coarse, /QUIET)

    dx1 = dx * (Nx-1.) / (Nx_2 - 1.) 
    df1 = level_up(0.*r1,r1,dx1,n-1)
    ;
    f = f - refine(df1, x_coarse, x)

  endif
  f = gauss_seidel(f,g,dx)
  ; f = jacobi(f,g,dx)
  return, f
end
; ---------------------------------------------------------------------- ;
pro print_depth, n, m, ODD=odd
  ;; Some diagnostic output that visually indicates the level we are at
  default, odd, ''
  indent = '                 '
  indent = strmid(indent,1,17-n)
  if (odd eq '') then begin
    print, indent, 'n, m      = ', n, m
  endif else begin
    print, indent, 'n, m, odd = ', n, m, odd
  endelse
  ;
end
; ---------------------------------------------------------------------- ;

niter = 15
Nx = 34                          ; typically 2^k+1 or 2^k

x = linspace(0.D,1.,Nx)
dx = x[1]-x[0]
;f = 0.01*cos(!pi*x/dx) & f[0]=0. & f[Nx-1]=0.
f = 0.*x
g = -q(x)                        ; heating 

if (Nx le 33) then $
    psym = -1 $
else $
    psym = 0

plot, x, f, $
      YRANGE=[0,1]*0.25, YSTYLE=3, $
      TITLE='!8N!6='+strtrim(Nx,2)+'!X'

; for k=0,niter-1 do begin
;   f = gauss_seidel(f, g, dx)
;   oplot, x, f, COLOR=k*150./niter
; endfor

colors = [-1]                   ; initialize, so we can append below
labels = ['ignore']
for i=0,niter-1 do begin
;  f = gauss_seidel(f,g,dx)
  f = level_up(f,g,dx,1000)
  color = i*150./niter
  oplot, x, f, $
         COLOR=color, $
         PSYM=psym
  colors = [colors, color]
  labels = [labels, '!8i!6='+strtrim(i,2)+'!X']
  ;; Store f values
  if (i eq 0) then ff = [[f]] else ff = [[ff],[f]]
endfor

ophline, 0.204050
  esrg_legend, labels[1:*], $
               COLOR=colors[1:*], $
               /BOX, POS=[0.45, 0.1, 0.7, 0.4], $
               XSTRETCH=0.5


end
; End of file multigrid.pro
