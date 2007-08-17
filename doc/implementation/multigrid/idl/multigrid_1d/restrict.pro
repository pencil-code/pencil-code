;;;;;;;;;;;;;;;;;;;;;;;;
;;;   restrict.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   01-Mar-2007
;;;
;;;  Description:
;;;    Investigate different restriction (coarse-graining) schemes.
;;;    Run `.run restrict_test' to test these schemes, or use them in
;;;    multigrid.pro.
;;;  Discussion:
;;;    In what I call `classical multigrid' (vertex-centred multigrid),
;;;    the coarser grid's points represent every other point of the finer
;;;    grid, I.e. coarse-graining goes from 2^N+1 points to 2^{N-1}+1
;;;    points, with identical boundary points. Restriction is then
;;;    typically done using a (1,2,1)/4 stencil, which completely
;;;    eliminates any Nyquist signal on the finer grid.
;;;      For the Pencil Code (and other codes), however, it is way more
;;;    natural to go from 2^N to 2^{N-1} points, with identical boundary
;;;    points (cell-centred multigrid). This implies that the coarse-grid
;;;    points do not normally coincide with fine-grid points, and we need
;;;    to find a good way for doing the restriction (and probably later
;;;    also the refinement) step.

; ---------------------------------------------------------------------- ;
function mg_find_index, x, x1
;;
;; Return array of indices representing x_coarse's positions on the x
;; grid. Indices are floating point numbers, not (necessarily) integers.
;;
  N1   = n_elements(x1)
  idxx = 0.*x1

  for i=0, N1-1 do begin

    x1i = x1[i]
    idx = min(where(x1i le x))

    if (idx lt 0) then begin
      message, 'argument out of range for x1 = ' + wdtrim(x1i,2)
    endif

    if (idx eq 0.) then begin
      idxx[i] = 0.
    endif else begin
      idx_l = idx - 1
      idx_r = idx
      ;
      dx = x[idx_r] - x[idx_l]
      p = (x1i-x[idx_l]) / dx
      q = 1. - p
      if (abs(p-0.5) gt 0.5) then begin
        message, 'p = ' + wdtrim(p) + ' is off limits'
      endif
      idxx[i] = q*idx_l + p*idx_r
    endelse

  endfor

  return, idxx

end
; ---------------------------------------------------------------------- ;
function restrict1, f, x, x_coarse, $
                    QUIET=quiet
;;
;; Restrict by linearly interpolating the left and right (1,2,1)/4
;; restriction (which is the one normally used in classical multigrid).
;;
  N_coarse = n_elements(x_coarse)
  N_fine   = n_elements(x)

  idxx = mg_find_index(x, x_coarse)
  idxx_l = floor(idxx)
  f_l = 0.25 * (f[idxx_l-1] + 2*f[idxx_l] + f[idxx_l+1])
  idxx_r = idxx_l+1
  f_r = 0.25 * (f[idxx_r-1] + 2*f[idxx_r] + f[idxx_r+1])
  pp = idxx - idxx_l
  qq = 1. - pp
  f_coarse = qq*f_l + pp*f_r

;  f_coarse[0]          = 0.25 * (3*f[0] + 2*f[1] - f[2])
;  f_coarse[N_coarse-1] = 0.25 * (3*f[N_fine-1] + 2*f[N_fine-2] - f[N_fine-3])

  f_coarse[0]          = 0.
  f_coarse[N_coarse-1] = 0.

  return, f_coarse

end
; ---------------------------------------------------------------------- ;
function restrict2, f, x, x_coarse, $
                    QUIET=quiet
;;
;; Restriction operator using four nearest neighbours and modified Bessel
;; interpolation to cut off Nyquist signal.
;; Indices:
;;               + point
;;               |
;; ----+-------+-------+-------+----
;;     |       |       |       |
;;   idx_l-1  idx_l   idx_r   idx_r+1
;;
;; ----+-------+-------+-------+------> t
;;    -3/2    -1/2     1/2    3/2
;;
;; We interpolate this for (-0.5 <= t <= 0.5)
;
  N_coarse = n_elements(x_coarse)
  N_fine   = n_elements(x)

  idxx = mg_find_index(x, x_coarse)
  idxx_l = floor(idxx)
  idxx_r = idxx_l + 1

  y_mean = 0.5 * (f[idxx_l] + f[idxx_r])
  a      = 0.25 * (-f[idxx_l-1] - f[idxx_l] + f[idxx_r] + f[idxx_r+1])
;  a      = f[idxx_r] - f[idxx_l]
  b      = 0.25 * (f[idxx_l-1] - f[idxx_l] - f[idxx_r] + f[idxx_r+1])

  tt = -0.5 + idxx - idxx_l
  f_coarse = y_mean + a*tt + b*(tt^2-0.25)

;   f_coarse[0]          = 0.125 * (  7*f[0] $
;                                   + 3*f[1] $
;                                   - 3*f[2] $
;                                   +   f[3])
;   f_coarse[N_coarse-1] = 0.125 * (  7*f[N_fine-1] $
;                                   + 3*f[N_fine-2] $
;                                   - 3*f[N_fine-3] $
;                                   +   f[N_fine-4])
   f_coarse[0]          = 0.
   f_coarse[N_coarse-1] = 0.

  return, f_coarse

end
; ---------------------------------------------------------------------- ;

; End of file restrict.pro
