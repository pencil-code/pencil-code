;; This piece of code uses the IDL astro-lib routine READCOL,
;; which must be available in the IDL !PATH.
;;
pro pc_d2_dimension_plot, file, psfilename=psfilename
  compile_opt idl2

  if ~n_elements(file) then begin
    file = 'nnd_d2_beskow.dat'
    if ~file_test(file, /read, /regular) then file = 'nnd_d2.dat'
  endif

  if ~file_test(file, /read, /regular) then begin
    message, 'Cannot read the file "' + file + '".'
  endif
  if ~n_elements(psfilename) then psfilename = 'd2.eps'
  
  format = '(l, f, d, d, d)'
  readcol, file, model, ap0, d2, xmin_mean, xmin_var, format=format

  midx = model[uniq(model)]
  n_model = n_elements(midx)

  L = 6.28d0 * 2
  rho_gr = 1d3
  rho_mean = 1d0

  device, decomposed=0
  loadct, 0

  alpha = rho_gr / rho_mean * ap0[midx] / L
  alpha = ap0[midx]
  xrange = [0.9 * min(ap0), 1.1 * max(ap0)]
  yrange = [2.1d0, 3.2d0]
  if max(d2) gt yrange[1L] then yrange[1L] = max(d2)
  linestyle = [0, 3, 2, 4, 1, 5]
  thick = 2.0
  
  entrydevice = !d.name & entryfont = !p.font
  entrymulti = !p.multi & entryx = !x & entryy = !y
  !p.multi = [0L, 1L, 1L, 0, 0]

  for jj = 0L, 1L do begin

    if jj then begin
      set_plot, 'ps'
      device, filename=psfilename, /isolatin, /color, bits_per_pixel=8, $
              /landscape, xsize=psxsize, ysize=psysize, $
              xoffset=psxoffset, yoffset=psyoffset
      !p.font = 0L
    endif else begin
      window, 0
    endelse

    ;;!x.omargin = [0.0, 0.0]
    ;;!y.omargin = [1.0, 4.0]
    charsize = 2.0
    xmargin = [8.0, 5.0]
    ymargin = [3.5, 0.2]
    symsize = 1.0

    plot, [0], background=255, color=0, /nodata, $
          charsize=charsize, charthick=thick, /xlog, $
          xrange=xrange, /xstyle, xthick=thick, xtitle='!4a!x', $
          yrange=yrange, /ystyle, ythick=thick, ytitle='d!d2!n'

    oplot, 1d1 ^!x.crange, [3d0, 3d0], $
           linestyle=1, color=150, thick=thick

    for i = 0L, n_model - 1L do begin
      idx = where(model eq midx[i], count)

      alpha = rho_gr / rho_mean * ap0[idx] / L
      alpha = ap0[idx]

      linestyle_ = linestyle[i mod n_elements(linestyle)]
      oplot, alpha, d2[idx], linestyle=linestyle_, thick=thick, color=0
    endfor

    if jj then begin
      device, /close_file
      set_plot, entrydevice
    endif

  endfor

  !x = entryx
  !y = entryy
  !p.font = entryfont
  !p.multi = entrymulti

  return
end ;;; procedure: pc_d2_dimension_plot
