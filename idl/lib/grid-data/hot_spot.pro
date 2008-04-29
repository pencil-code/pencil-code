;;;;;;;;;;;;;;;;;;;;;;;;
;;;   hot_spot.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   28-Jun-2001
;;;
;;;  Description:
;;;   Locate the maximal modulus (or the NaN/Inf points) of gridded
;;;   data. Returns the indices of the `hot spot', which is either the
;;;   location of max|f|, or the center of mass of degenerate
;;;   numerical values NaN, Inf.
;;;  Arguments:
;;;   F         -- scalar or vector field
;;;   X,Y,Z     -- coordinate vector
;;;  Keywords:
;;;   VECTORIAL -- flag indicating that argument is a vector and the
;;;                vector modulus is to be used
;;;                  If the last dimension of F is 2 or 3, vectorial
;;;                character is automatically assumed, but a warning
;;;                is issued (unless suppressed with /QUIET)
;;;   GHOST     -- data contain ghost cells to be excluded from the
;;;                maximum; can be a scalar or a vector (if different
;;;                directions have different number of ghost cells)
;;;   QUIET     -- set this to suppress warnings and any kind of
;;;                diagnostic messages
;;;   DEBUG     -- set this to get many diagnostic messages indicating
;;;                what is going on
;;;  Bugs:
;;;   - 2d vector arguments are unlikely to work correctly

function hot_spot, f,  $
                   xcoord, ycoord, zcoord, $
                   VECTORIAL=vector, $
                   GHOST=ghost, $
                   QUIET=quiet, DEBUG=debug

  default, vector, 0
  default, ghost, 0
  default, debug, 0
  default, quiet, 0
  if (debug) then quiet=0

  s = size(f)
  if ((s[s[0]] eq 2) or (s[s[0]] eq 3)) then begin
    vector=1
    if (not quiet) then message, 'assuming vector data', /INFO
  endif
  if (n_elements(ghost) eq 1) then ghost = [ghost,ghost,ghost]
  l1 = ghost[0]   & l2 = s[1]-ghost[0]-1
  if (s[0] ge 2) then begin
    m1 = ghost[0] & m2 = s[2]-ghost[0]-1
  endif else begin
    m1 = (m2 = 0)
  endelse
  if (s[0] ge 3) then begin
    n1 = ghost[0] & n2 = s[3]-ghost[0]-1
  endif else begin
    n1 = (n2 = 0)
  endelse
  
  g = f[l1:l2,m1:m2,n1:n2,*]

  if (vector) then begin
    g = sqrt(total(float(abs(1.D*g)^2), s[0]))
  endif else begin
    g = abs(g)
  endelse

  if (debug) then help, g

  ;; Create grid of index values
  s = size(g)                   ; may be smaller now than f
  nx = s[1] & ny = 0 & nz = 0   ; ny=0 means no y-direction
  xx = findgen(nx)
  yy = 0.*xx                    ; gets overwritten if necessary
  zz = 0.*xx
  if (s[0] gt 1) then begin
    ny = s[2]
    xx = spread(xx, 1, ny)
    yy = spread(findgen(ny), 0, nx)
    zz = 0.*xx
  endif
  if (s[0] gt 2) then begin
    nz = s[3]
    xx = spread(xx, 2, nz)
    yy = spread(yy, 2, nz)
    zz = spread(findgen(nz), [0,1], [nx,ny])
  endif

  if (debug) then help, xx, yy, zz

  nan = (finite(g) eq 0)
  idx_nan = where(nan)
  
  if (idx_nan[0] ge 0) then finite=0 else finite=1

  if (finite) then begin        ; all data are finite
    if (not quiet) then $
        print, '% HOT_SPOT: Showing location of maximum value'
    ;; Normalise g to interval [0,1]
    gmax = max(g,MIN=gmin)
    g = (g-gmin)/(gmax-gmin)
    mxidx = where(g eq max(g))
    if ((n_elements(mxidx) gt 1) and (not quiet)) then $
        message, 'More than one maximum', /INFO
    ipos = mxidx[0]
    if (ipos lt 0) then message, 'There is an inconsistency'
    pos = [ipos mod nx]
    if (ny gt 0) then pos = [pos, ipos/nx mod ny]
    if (nz gt 0) then pos = [pos, ipos/nx/ny mod nz]    
    if (n_elements(idx_g) le 1) then pos = long(pos)
    ;; Determine first and second momenta of x,y,z
    if (not quiet) then begin
      norm = 1./total(g)
      mom1 = [total(xx*abs(g))*norm]
      mom2 = [total((xx-mom1[0])^2*abs(g))*norm]
      if (ny gt 0) then begin
        mom1 = [mom1, total(yy*abs(g))*norm]
        mom2 = [mom2, total((yy-mom1[1])^2*abs(g))*norm]
      endif
      if (nz gt 0) then begin
        mom1 = [mom1, total(zz*abs(g))*norm]
        mom2 = [mom2, total((zz-mom1[2])^2*abs(g))*norm]
      endif
      print, 'center of mass:', mom1+ghost
      print, 'width         :', sqrt(mom2)
    endif
  endif else begin              ; there were NaNs
    if (not quiet) then $
        print, '% HOT_SPOT: encountered ',  $
        strtrim(n_elements(idx_nan),2), ' NaN/Inf values'
    gmax = g[min(where(not finite(g)))]
    norm = 1./total(nan)
    pos = [total(xx*abs(nan))*norm]
    if (ny gt 0) then pos = [pos, total(yy*abs(nan))*norm]
    if (nz gt 0) then pos = [pos, total(zz*abs(nan))*norm]
    if (n_elements(idx_nan) le 1) then pos = long(pos)
  endelse

  ;; Print coordinates if possible:
  if (not quiet) then begin
    if (nx gt 0) then begin
      x=xcoord[l1:l2]
      print, FORMAT='(A,F0,$)', 'X = ', x[pos[0]]
    endif
    if (ny gt 0) then begin
      ;; is it the full 2d or 3d coordinate array?
      if ((size(ycoord))[0] gt 1) then y=ycoord[0,m1:m2] else y=ycoord[m1:m2]
      print, FORMAT='(A,F0,$)', ', Y = ', y[pos[1]]
    endif
    if (nz gt 0) then begin
      if ((size(zcoord))[0] gt 2) then z=zcoord[0,0,n1:n2] else z=zcoord[n1:n2]
      print, FORMAT='(A,F0,$)', ', Z = ', z[pos[2]]
    endif
    if (nx gt 0) then print ; trailing carriage return
    print, 'value: ', gmax
  endif

  return, pos+ghost

end
; End of file hot_spot.pro
