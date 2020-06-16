pro pc_rhonp__height, height, charsize=charsize, $
        npx=npx, npy=npy, width=width, xmargin=xmargin, ymargin=ymargin
  compile_opt hidden, IDL2

  panel_wh = long((width - total(!x.omargin) * charsize * !d.x_ch_size) / $
                  double(npx) - total(xmargin) * charsize * !d.x_ch_size)

  height = long(npy * (panel_wh + total(ymargin) * charsize * !d.y_ch_size) + $
                total(!y.omargin) * charsize * !d.y_ch_size)


  return
end ;;; procedure: pc_rhonp__height


pro pc_rhonp__colorbar, charsize=charsize, ch2=ch2, charthick=charthick, font=font, $
        npx=npx, npy=npy, v_minmax=v_minmax, col_minmax=col_minmax, $
        colorbg=colorbg, colorfg=colorfg, thick=thick, title=title, top=top
  compile_opt hidden, IDL2

  x0 = !x.window[0L]
  x1 = !x.window[1L]
  if keyword_set(top) then begin
    y0 = !y.window[1L] + 2.7d0 * !d.y_ch_size * ch2 / !d.y_size
    y1 = !y.window[1L] + 3.7d0 * !d.y_ch_size * ch2 / !d.y_size
  endif else begin
    y0 = !y.window[0L] - 4.5d0 * !d.y_ch_size * ch2 / !d.y_size
    y1 = !y.window[0L] - 3.5d0 * !d.y_ch_size * ch2 / !d.y_size
  endelse

  bar_x = !d.x_size * (x1 - x0)
  bar_y = ceil(!d.y_size * (y1 - y0))

  xrange = v_minmax
  yrange = [0d0, 1d0]
  bar = byte(lindgen(bar_x) / (bar_x - 1d0) * (col_minmax[1L] - col_minmax[0L]) + col_minmax[0L])

  bar = rebin(bar, bar_x, bar_y)

  px0 = x0 * !d.x_size
  py0 = ceil(y0 * !d.y_size)

  xstyle = 1 + 4
  ystyle = 1 + 4
  ytickformat = '(a1)'
  plot, [0], /nodata, background=colorbg, charsize=charsize, charthick=charthick, $
        color=colorfg, font=font, /noerase, position=[x0, y0, x1, y1], $
        xrange=xrange, xstyle=xstyle, xthick=thick, $
        yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat

  tv, bar, px0, py0

  tickformat = '(a1)'
  ticklen = 0.2
  axis, yaxis=0, charsize=charsize, charthick=charthick, color=colorfg, $
        yminor=0L, yrange=yrange, /ystyle, ythick=thick, ytickformat=tickformat, ticklen=0
  axis, yaxis=1, charsize=charsize, charthick=charthick, color=colorfg, $
        yminor=0L, yrange=yrange, /ystyle, ythick=thick, ytickformat=tickformat, ticklen=0
  axis, xaxis=0, charsize=charsize, charthick=charthick, color=colorfg, $
        xrange=xrange, /xstyle, xthick=thick, xtitle=title, ticklen=ticklen, font=font
  axis, xaxis=1, charsize=charsize, charthick=charthick, color=colorfg, $
        xrange=xrange, /xstyle, xthick=thick, ticklen=ticklen, xtickformat=tickformat


  return
end ;;; procedure: pc_rhonp__colorbar


pro pc_rhonp, varfile, pvarfile, histogram=histogram, log=log, ylog=ylog, $
        nbinrho=nbinrho_, nbinnp=nbinnp_, $
        minlgrho=minlgrho, maxlgrho=maxlgrho, minnp=minnp_, maxnp=maxnp_, $
        charsize=charsize_in_, charthick=charthick, colortable=colortable, $
        psym=psym, symsize=symsize, set_font=set_font, t_start=t_start, $
        xsize=xsize, cmxsize=cmxsize, dpi=dpi, thick=thick_, font_size=font_size, $
        legend_log_aps=legend_log_aps, legendtext=legendtext, xrange=xrange, $
        ofilename=ofilename, common_range=common_range, colorbar_top=colorbar_top, $
        ret_val=ret_val, use_zbuffer=use_zbuffer, normalize_counts=normalize_counts
  compile_opt IDL2

  if ~n_elements(histogram) then histogram = 1L
  if ~n_elements(log) then log = 1L
  if ~n_elements(ylog) then ylog = 1L
  if ~n_elements(charsize_in_) then charsize_in_ = 1d0
  if ~n_elements(colortable) then colortable = 74
  if ~n_elements(psym) then psym = 3
  if ~n_elements(symsize) then symsize = 1.0
  if ~n_elements(xsize) then xsize = 800L
  if ~n_elements(cmxsize) then cmxsize = 13.5
  if n_elements(dpi) eq 1L then xsize = ceil(cmxsize * dpi / 2.54d0)
  if ~n_elements(thick_) then thick_ = 1.0
  if ~n_elements(t_start) then t_start = 0d0
  if ~n_elements(legend_log_aps) then legend_log_aps = 1L
  font = ~n_elements(set_font) ? - 1L : 1L
  tt_font = font eq 1L ? 1L : 0L
  common_range = keyword_set(common_range)
  colorbar_top = keyword_set(colorbar_top)
  if ~n_elements(use_zbuffer) then use_zbuffer = 1L
  if ~n_elements(normalize_counts) then normalize_counts = 1L

  pc_read_param, object=ps, /quiet
  pc_read_dim, object=pd, /quiet
  pc_read_var, object=v_object, varfile=varfile, /trimall
  pc_read_pvar, object=p_object, varfile=pvarfile
  pc_read_ts, object=ts

  np_ap = dblarr(pd.nx, pd.ny, pd.nz, n_elements(ps.ap0))

  x0 = ps.xyz0[0L] & dx = p_object.dx
  y0 = ps.xyz0[1L] & dy = p_object.dy
  z0 = ps.xyz0[2L] & dz = p_object.dz

  alpha = ps.rhopmat / ps.rho0 * ps.ap0 / ps.Lx0
  if legend_log_aps then alpha = alog10(alpha)

  ;; Calculate particle densities:
  message, 'Loop over all ' + strtrim(p_object.npar_found, 2L) + $
           ' particles.', /informational

  for i = 0UL, p_object.npar_found - 1UL do begin
    dx_x = (p_object.xp[i] - x0) / dx
    dx_y = (p_object.yp[i] - y0) / dy
    dx_z = (p_object.zp[i] - z0) / dz

    idx_xn = long(dx_x) & idx_xp = idx_xn + 1L
    idx_yn = long(dx_y) & idx_yp = idx_yn + 1L
    idx_zn = long(dx_z) & idx_zp = idx_zn + 1L

    value_xp = dx_x - idx_xn & value_xn = 1d0 - value_xp
    value_yp = dx_y - idx_yn & value_yn = 1d0 - value_yp
    value_zp = dx_z - idx_zn & value_zn = 1d0 - value_zp

    idx_a = where(ps.ap0 eq p_object.ap[i])

       np_ap[idx_xn, idx_yn, idx_zn, idx_a] += value_xn * value_yn * value_zn
    if idx_xp lt pd.nx then $
       np_ap[idx_xp, idx_yn, idx_zn, idx_a] += value_xp * value_yn * value_zn
    if idx_yp lt pd.ny then $
       np_ap[idx_xn, idx_yp, idx_zn, idx_a] += value_xn * value_yp * value_zn
    if idx_xp lt pd.nx && idx_yp lt pd.ny then $
       np_ap[idx_xp, idx_yp, idx_zn, idx_a] += value_xp * value_yp * value_zn
    if idx_zp lt pd.nz then $
       np_ap[idx_xn, idx_yn, idx_zp, idx_a] += value_xn * value_yn * value_zp
    if idx_xp lt pd.nx && idx_zp lt pd.nz then $
       np_ap[idx_xp, idx_yn, idx_zp, idx_a] += value_xp * value_yn * value_zp
    if idx_yp lt pd.ny && idx_zp lt pd.nz then $
       np_ap[idx_xn, idx_yp, idx_zp, idx_a] += value_xn * value_yp * value_zp
    if idx_xp lt pd.nx && idx_yp lt pd.ny && idx_zp lt pd.nz then $
       np_ap[idx_xp, idx_yp, idx_zp, idx_a] += value_xp * value_yp * value_zp
  endfor

  ;np_ap /= dx * dy * dz

  p_multi = !p.multi
  x_omargin = !x.omargin
  y_omargin = !y.omargin

  npx = ceil(sqrt(n_elements(ps.ap0)))
  npy = floor(sqrt(n_elements(ps.ap0)))
  !p.multi = [0L, npx, npy, 0L, 0L]


  ;;===========================================================================
  ;;===========================================================================
  ;; Time scales:

  ;; Drag force: stopping time:
  tau_stop = sqrt(!dpi / 8d0) * ps.rhopmat * ps.ap0 / (ps.rho0 * ps.cs0)

  ;; Mean turnover time:
  min = min(abs(ts.t - t_start), t_start_index)
  tau_turnover = 1d0 / mean(ts.urms[t_start_index : *])

  ;; Stokes; grain-size dependent:
  vStokes = strarr(n_elements(ps.ap0))
  for i = 0L, n_elements(ps.ap0) - 1L do vStokes[i] = $
     strtrim(string(tau_stop[i] / tau_turnover, format='(e9.2)'), 2L)


  ;;===========================================================================
  ;;===========================================================================
  ;; Set up color indices:

  colorbg = 1b
  colorfg = 0b
  colorrd = 2b
  colorgn = 4b
  colorbl = 5b
  coloror = 6b
  colorpr = 7b
  colorlb = 8b
  coloryl = 9b

  if use_zbuffer then begin
    scale_factor = 4L

    charsize_in = charsize_in_ * scale_factor
    charsize_label = charsize_in
    charsize = charsize_in * (npx * npy gt 2L ? 2d0 : 1d0)
    thick = thick_ * scale_factor
  endif else begin
    scale_factor = 1L

    charsize_in = charsize_in_
    charsize_label = charsize_in
    charsize = charsize_in * (npx * npy gt 2L ? 2d0 : 1d0)
    thick = thick_
  endelse

  device, get_decomposed=decomposed

  ;; Calculate the window height to get a 1.0 aspect ratio:

  !x.omargin = [2.0, 0.0]
  if normalize_counts then !x.omargin[1L] += 1.0
  !y.omargin = [0.0, 0.0]
  xmargin = [2.8, 0.3]
  ymargin = [histogram ? 4.3 : 0.3, 0.3]
  ymargin[colorbar_top ? 1L : 0L] += 3.5

  pc_rhonp__height, ysize__scale_factor, charsize=charsize_in, npx=npx, npy=npy, $
      width=xsize * scale_factor, xmargin=xmargin, ymargin=ymargin
  ysize = ysize__scale_factor / scale_factor

  if use_zbuffer then begin
    device_plot = !d.name
    set_plot, 'z'
    erase
    device, set_pixel_depth=24
    device, set_resolution=[xsize, ysize] * scale_factor
  endif else begin
    window, /free, xsize=xsize, ysize=ysize, /pixmap
    winid = !d.window
  endelse

  device, decomposed=0

  bottom = 10L
  ncolors = !d.table_size - bottom
  loadct, colortable, bottom=bottom, ncolors=ncolors, /silent
  tvlct,   0b,   0b ,  0b, colorfg   ;; Black
  tvlct, 255b, 255b, 255b, colorbg   ;; White
  tvlct, 213b,  94b,   0b, colorrd   ;; Vermillion
  tvlct,   0b, 158b, 115b, colorgn   ;; Bluish green
  tvlct,   0b, 114b, 178b, colorbl   ;; Blue
  tvlct, 230b, 159b,   0b, coloror   ;; Orange
  tvlct, 204b, 121b, 167b, colorpr   ;; Purple (Mulberry)
  tvlct,  86b, 180b, 233b, colorlb   ;; Sky Blue
  tvlct, 240b, 228b,  66b, coloryl   ;; Yellow

  if tt_font then device, set_font=set_font, /tt_font

  salpha = font eq 1 ? '!9a!x' : '!4a!x'
  srho = font eq 1 ? '!9r!x' : '!4q!x'

  lg_rho = v_object.lnrho / alog(1d1)

  xrange = minmax(lg_rho)
  if ~n_elements(minlgrho) then minlgrho = xrange[0L]
  if ~n_elements(maxlgrho) then maxlgrho = xrange[1L]
  nbinrho = ~n_elements(nbinrho) ? (xrange[1L] - xrange[0L]) / 4d2 : nbinrho_

  if common_range then begin
    yrange = minmax(np_ap)
    if ylog then yrange = alog10(yrange)
    minnp = ~n_elements(minnp_) ? yrange[0L] : minnp_
    maxnp = ~n_elements(maxnp_) ? yrange[1L] : maxnp_
    nbinnp = ~n_elements(nbinnp_) ? (yrange[1L] - yrange[0L]) / 4d2 : nbinnp_
  endif

  xstyle = histogram ? 5 : 1
  ystyle = histogram ? 5 : 1
  xtitle = 'log ' + srho
  ytitle_ = 'n!dp!ii!n'
  if ylog then ytitle_ = 'log ' + ytitle_
  symsize = 1.0

  if histogram && ~ylog then background_color = colorbg

  format = legend_log_aps ? '(f6.3)' : '(e10.3)'
  text = string(bindgen(26) + 97b)
  for j = 0UL, n_elements(ps.ap0) - 1UL do begin

    y_value = np_ap[*, *, *, j]

    if ~common_range then begin
      yrange = minmax(y_value)
      if ylog then begin
        tidx = where(y_value gt 0d0)
        yrange[0L] = min(y_value[temporary(tidx)])
        yrange = alog10(yrange)
      endif
      minnp = ~n_elements(minnp_) ? yrange[0L] : minnp_
      maxnp = ~n_elements(maxnp_) ? yrange[1L] : maxnp_
      nbinnp = ~n_elements(nbinnp_) ? (yrange[1L] - yrange[0L]) / 2d2 : nbinnp_
    endif

    null_value = - 10d0

    if histogram then begin
      hist_2d = hist_2d(lg_rho, ylog ? alog10(y_value) : y_value, $
                        bin1=nbinrho, bin2=nbinnp, $
                        min1=minlgrho, max1=maxlgrho, min2=minnp, max2=maxnp)
      if ~common_range then begin
        zidx = where(hist_2d ge 1L)
        zrange = minmax(hist_2d[zidx])
        if normalize_counts then zrange /= max(zrange)
      endif

      if log then begin
        idx = where(hist_2d eq 0d0, count, complement=cidx)
        max_h = double(max(hist_2d))
        hist_2d = alog10(normalize_counts ? hist_2d / max_h : hist_2d)
        min = min(hist_2d, /nan, max=max)
        if ~common_range then zrange = [min, max]
        color = byte((hist_2d - min) / (max - min) * (ncolors - 1d0) + bottom)
        if count ne 0L then color[idx] = colorbg
        color_range = minmax(color[cidx])
      endif else begin
        color = byte((hist_2d - min(hist_2d)) / $
                     (max(hist_2d) - min(hist_2d)) * (ncolors - 1d0) + bottom)
        color_range = minmax(color)
      endelse
    endif

    xtickformat = j ge n_elements(ps.ap0) - npx ? '' : '(a1)'
    ;;ytickformat = ~(j mod npx) ? '' : '(a1)'

    ;xtitle = j ge n_elements(ps.ap0) - npx ? xtitle_ : ''
    ytitle = ~(j mod npx) ? ytitle_ : ''

    plot, [0], /nodata, background=colorbg, charsize=charsize, charthick=charthick, $
          color=colorfg, font=font, xtitle=xtitle, ytitle=ytitle, $
          xmargin=xmargin, xrange=xrange, xstyle=xstyle, xthick=thick, xtickformat=xtickformat, $
          ymargin=ymargin, yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat

    if histogram then begin

      i_size_x = !d.x_size * (!x.window[1L] - !x.window[0L])
      i_size_y = !d.y_size * (!y.window[1L] - !y.window[0L])

      image = congrid(color, i_size_x, i_size_y, /center)

      tv, image, !x.window[0L] * !d.x_size, !y.window[0L] * !d.y_size

      axis, xaxis=0L, xrange=xrange, /xstyle, xthick=thick, xtitle=xtitle, $
            charsize=charsize, charthick=charthick, color=colorfg, font=font
      axis, xaxis=1L, xrange=xrange, /xstyle, xthick=thick, xtickformat='(a1)', $
            charsize=charsize, charthick=charthick, color=colorfg, font=font
      axis, yaxis=0L, yrange=yrange, /ystyle, ythick=thick, ytitle=ytitle, $
            charsize=charsize, charthick=charthick, color=colorfg, font=font
      axis, yaxis=1L, yrange=yrange, /ystyle, ythick=thick, ytickformat='(a1)', $
            charsize=charsize, charthick=charthick, color=colorfg, font=font

    endif else begin

      ;; Scatter plot, one color:
      plots, lg_rho, ylog ? alog10(np_ap[*, *, *, j] > 1d1 ^ null_value) : np_ap[*, *, *, j], $
             color=colorfg, psym=psym, symsize=symsize

    endelse


    item = legend_log_aps ? 'log(' + salpha + ') = ' : salpha + ' = '
    item += strtrim(string(alpha[j], format=format), 2L)
    item += legend_log_aps ? ', log(a) = ' : ', a = '
    item += strtrim(string(legend_log_aps ? alog10(ps.ap0[j]) : ps.ap0[j], format=format), 2L)
    if n_elements(ps.ap0) gt 1UL then item = strmid(text, j, 1L) + ')  ' + item
    item = [item, '    St = ' + vStokes[j]]

    xpos = 0.05 * (!x.window[1L] - !x.window[0L]) + !x.window[0L]
    ypos = 0.05 * (!y.window[1L] - !y.window[0L]) + !y.window[0L]
    xyouts, xpos, ypos, item[1L], alignment=0.0, charsize=charsize_label, charthick=charthick, $
            color=colorfg, font=font, /normal
    ypos += 1.4 * charsize_label * !d.y_ch_size / !d.y_size
    xyouts, xpos, ypos, item[0L], alignment=0.0, charsize=charsize_label, charthick=charthick, $
            color=colorfg, font=font, /normal

    ;al_legend, item, charsize=charsize_label, charthick=charthick, textcolors=colorfg, $
    ;           box=histogram && ~ylog, /left, top_legend=~ylog, font=font, $
    ;           background_color=background_color

    title = (log ? 'log' : '') + (normalize_counts ? ' normalized' : '') + ' counts'
    pc_rhonp__colorbar, charsize=charsize, ch2=charsize_in, charthick=charthick, $
        font=font, col_minmax=color_range, npx=npx, npy=npy, top=colorbar_top, $
        v_minmax=zrange, colorbg=colorbg, colorfg=colorfg, thick=thick, title=title

    ; [bottom, !d.table_size]

  endfor


  if use_zbuffer then begin

    ;; Allows for antialiasing, by plotting to a figure 4 times the
    ;; size and then scaling down the read out image:

    image = tvrd(/true)

    set_plot, device_plot
    window, xsize=xsize, ysize=ysize
    image = rebin(image, 3L, xsize, ysize)

    tv, image, /true

  endif else begin

    device, decomposed=1L
    image = tvrd(/true)

    copy = [0L, 0L, xsize, ysize, 0L, 0L, winid]
    window, xsize=xsize, ysize=ysize
    device, copy=copy
    wdelete, winid

  endelse

  device, decomposed=decomposed

  !p.multi = p_multi
  !x.omargin = x_omargin
  !y.omargin = y_omargin

  ;; Save the file:
  ofilename = 'pc_rhonp_' + varfile + '.png'
  if strpos(ofilename, '.png', /reverse_search) eq strlen(ofilename) - 4L then begin
    tvlct, red, green, blue, /get
    write_png, ofilename, image, red, green, blue
  endif else if strpos(ofilename, '.jpg', /reverse_search) eq strlen(ofilename) - 4L then begin
    write_jpeg, ofilename, image, /true, quality=95
  endif else if strpos(ofilename, '.bmp', /reverse_search) eq strlen(ofilename) - 4L then begin
    write_bmp, ofilename, image, /rgb
  endif else if strpos(ofilename, '.tif', /reverse_search) eq strlen(ofilename) - 4L then begin
    write_tiff, ofilename, reverse(image, 3), 1, /compression
  endif

end ;;; procedure: pc_rhonp
