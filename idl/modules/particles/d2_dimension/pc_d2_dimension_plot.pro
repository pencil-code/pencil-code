;; This piece of code uses the IDL astro-lib routine READCOL,
;; which must be available in the IDL !PATH.
;;
pro pc_d2_dimension_plot, file, charsize=charsize, thick=thick, $
                          psxsize=psxsize, psysize=psysize, $
                          font_size=font_size, psfilename=psfilename, $
                          yrange=yrange, legendtext=legendtext, $
                          yright=yright, use_last_only=use_last_only
  compile_opt idl2

  if ~n_elements(file) then begin
    file = 'nnd_d2_beskow.dat'
    if ~file_test(file, /read, /regular) then file = 'nnd_d2.dat'
  endif

  if ~file_test(file, /read, /regular) then begin
    message, 'Cannot read the file "' + file + '".'
  endif
  if ~n_elements(psfilename) then psfilename = 'd2.eps'

  if ~n_elements(font_size) then font_size = 10.0
  if ~n_elements(charsize) then charsize = 1.0
  if ~n_elements(psxsize) then psxsize = 6.50
  if ~n_elements(psysize) then psysize = 4.00
  if ~n_elements(thick) then thick = 1.0
  use_last_only = keyword_set(use_last_only)

  format = '(l, f, d, d, d)'
  readcol, file, model, ap0, d2, xmin_mean, xmin_var, format=format

  midx = model[uniq(model)]
  n_model = n_elements(midx)

  ;; Retrieve model-specific values (instead of hard-coding them):
  pc_read_param, object=ps
  pc_read_param, object=ps2, /param2

  device, decomposed=0
  loadct, 0

  alpha = ps.rhopmat / ps.rho0 * ps.ap0 / ps.Lx0
  xrange = [0.9 * min(alpha), 1.1 * max(alpha)]

  if ~n_elements(yrange) then yrange = [2.1d0, 3.09d0]
  if max(d2) gt yrange[1L] then yrange[1L] = max(d2)
  linestyle = [0, 1, 2, 3, 5]

  entrydevice = !d.name & entryfont = !p.font
  entrymulti = !p.multi & entryx = !x & entryy = !y
  !p.multi = [0L, 1L, 1L, 0, 0]

  for jj = 0L, 1L do begin

    if jj then begin
      set_plot, 'ps'
      device, filename=psfilename, /isolatin, /color, bits_per_pixel=8, $
              /portrait, xsize=psxsize, ysize=psysize, /encapsulated, $
              font_size=font_size, inches=0
      !p.font = 0L
      xtitle = '!9a!x = !9r!x!dgr!na ' + string(215b) + $
               ' (!9' + string(225b) + 'r' + string(241b) + '!xL)!u-1!n'
      ytitle_pre = 'd!d2!n'
    endif else begin
      window, 0
      xtitle = '!4a!x = !4q!x!dgr!na ' + string(215b) + $
               ' (!13' + string(60b) + '!x!4q!x!13' + $
               string(62b) + '!xL)!u-1!n'
      ytitle_pre = 'd!d2!n'
    endelse

    xmargin = keyword_set(yright) ? [0.3, 6.5] : [6.5, 0.3]
    ymargin = [3.3, 0.3]
    symsize = 1.0

    ystyle = 1L + (keyword_set(yright) ? 8L : 0L)
    ytitle = keyword_set(yright) ? '' : ytitle_pre
    ytickformat = keyword_set(yright) ? '(a1)' : ''

    plot, [0], background=255, charsize=charsize, color=0, /nodata, $
          charthick=thick, /xlog, /xstyle, ystyle=ystyle, $
          xmargin=xmargin, xrange=xrange, xthick=thick, xtitle=xtitle, $
          ymargin=ymargin, yrange=yrange, ythick=thick, ytitle=ytitle, $
          ytickformat=ytickformat

    if keyword_set(yright) then axis, color=0, /yaxis, yrange=yrange, $
         charsize=charsize, charthick=charthick, $
         /ystyle, ythick=thick, ytitle=ytitle_pre

    oplot, 1d1 ^!x.crange, [3d0, 3d0], $
           linestyle=1, color=150, thick=thick

    i_start = use_last_only ? n_model - 1L : 0L
    for i = i_start, n_model - 1L do begin
      idx = where(model eq midx[i], count)

      alpha = ps.rhopmat / ps.rho0 * ps.ap0 / ps.Lx0

      linestyle_ = linestyle[i mod n_elements(linestyle)]
      oplot, alpha, d2[idx], linestyle=linestyle_, psym=-4, $
             thick=thick, color=0
    endfor

    if n_elements(legendtext) eq 1L then begin
      al_legend, legendtext, /bottom, /left, charsize=charsize, box=0
    endif

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
