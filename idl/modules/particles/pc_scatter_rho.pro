pro pc_scatter_rho__height, height, charsize=charsize, $
        npx=npx, npy=npy, width=width, xmargin=xmargin, ymargin=ymargin
  compile_opt hidden, IDL2

  panel_wh = long((width - total(!x.omargin) * charsize * !d.x_ch_size) / double(npx) - $
                  total(xmargin) * charsize * !d.x_ch_size)

  height = long(npy * (panel_wh + total(ymargin) * charsize * !d.y_ch_size) + $
                total(!y.omargin) * charsize * !d.y_ch_size)


  return
end ;;; procedure: pc_scatter_rho__height


pro pc_scatter_rho, varfile, pvarfile, $
        charsize=charsize_in, charthick=charthick, colortable=colortable, $
        psym=psym, symsize=symsize, set_font=set_font, t_start=t_start, $
        xsize=xsize, cmxsize=cmxsize, dpi=dpi, thick=thick, font_size=font_size, $
        legend_log_aps=legend_log_aps, legendtext=legendtext, xrange=xrange, $
        ofilename=ofilename, common_range=common_range

  compile_opt IDL2

  if ~n_elements(charsize_in) then charsize_in = 1d0
  if ~n_elements(colortable) then colortable = 74
  if ~n_elements(psym) then psym = 3
  if ~n_elements(symsize) then symsize = 1.0
  if ~n_elements(xsize) then xsize = 800L
  if ~n_elements(cmxsize) then cmxsize = 13.5
  if n_elements(dpi) eq 1L then xsize = cmxsize * dpi / 2.54d0
  if ~n_elements(t_start) then t_start = 0d0
  if ~n_elements(legend_log_aps) then legend_log_aps = 1L
  font = ~n_elements(set_font) ? - 1L : 1L
  tt_font = font eq 1L ? 1L : 0L

  pc_read_param, object=ps, /quiet
  pc_read_dim, object=pd, /quiet
  pc_read_var, object=v_object, varfile=varfile, /trimall
  pc_read_pvar, object=p_object, varfile=pvarfile
  pc_read_ts, object=ts

  np_ap = dblarr(pd.nx, pd.ny, pd.nz, n_elements(ps.ap0))

  x0 = ps.xyz0[0L] & dx = p_object.dx
  y0 = ps.xyz0[1L] & dy = p_object.dy
  z0 = ps.xyz0[2L] & dz = p_object.dz
  ;x0 = p_object.x[pd.nghostx] & dx = p_object.dx
  ;y0 = p_object.y[pd.nghosty] & dy = p_object.dy
  ;z0 = p_object.z[pd.nghostz] & dz = p_object.dz

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

    if value_xn lt 0d0 || value_yn lt 0d0 || value_zn lt 0d0 then begin
       print,"n ",i, value_xn, value_yn, value_zn
       stop
    endif
    if value_xp lt 0d0 || value_yp lt 0d0 || value_zp lt 0d0 then begin
       print,"p ",i, value_xp, value_yp, value_zp
       stop
    endif
    
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

  device, get_decomposed=decomposed
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

  ;; Calculate the window height to get a 1.0 aspect ratio:

  !x.omargin = [2.0, 0.0]
  !y.omargin = [3.0, 0.0]
  xmargin = [4.3, 0.3]
  ymargin = [0.3, 0.3]

  pc_scatter_rho__height, ysize, charsize=charsize_in, npx=npx, npy=npy, $
      width=xsize, xmargin=xmargin, ymargin=ymargin

  window, /free, xsize=xsize, ysize=ysize, /pixmap

  salpha = font eq 1 ? '!9a!x' : '!4a!x'
  srho = font eq 1 ? '!9r!x' : '!4q!x'

  charsize_label = charsize_in
  charsize = charsize_in * (npx * npy gt 2L ? 2d0 : 1d0)

  lg_rho = v_object.lnrho / alog(1d1)

  xrange = minmax(lg_rho)
  xstyle = 1
  ystyle = 1
  xtitle_ = 'log ' + srho
  ytitle_ = 'n!dp!ii!n'
  symsize = 1.0

  format = legend_log_aps ? '(f6.3)' : '(e10.3)'
  text = string(bindgen(26) + 97b)
  for j = 0UL, n_elements(ps.ap0) - 1UL do begin

    yrange = minmax(np_ap[*, *, *, j])

    xtickformat = j ge n_elements(ps.ap0) - npx ? '' : '(a1)'
    ;;ytickformat = ~(j mod npx) ? '' : '(a1)'

    xtitle = j ge n_elements(ps.ap0) - npx ? xtitle_ : ''
    ytitle = ~(j mod npx) ? ytitle_ : ''

    plot, [0], /nodata, background=colorbg, charsize=charsize, charthick=charthick, $
          color=colorfg, font=font, xtitle=xtitle, ytitle=ytitle, $
          xmargin=xmargin, xrange=xrange, xstyle=xstyle, xthick=thick, xtickformat=xtickformat, $
          ymargin=ymargin, yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat

    plots, lg_rho, np_ap[*, *, *, j], color=colorfg, psym=psym, symsize=symsize

    ;; Stokes number label:
    ;item = 'St = ' + vStokes[j]
    ;al_legend, item, charsize=charsize_label, charthick=charthick, box=0, textcolors=colorfg, $
    ;           /right, background_color=colorbg, outline_color=colorfg, font=font

    item = legend_log_aps ? 'log(' + salpha + ') = ' : salpha + ' = '
    item += strtrim(string(alpha[j], format=format), 2L)
    item += legend_log_aps ? ', log(a) = ' : ', a = '
    item += strtrim(string(legend_log_aps ? alog10(ps.ap0[j]) : ps.ap0[j], format=format), 2L)
    item = [item, '   St = ' + vStokes[j]]
    if n_elements(ps.ap0) gt 1UL then item = strmid(text, j, 1L) + ')  ' + item
    al_legend, item, charsize=charsize_label, charthick=charthick, textcolors=colorfg, box=0, $
                     /left, /top_legend, background_color=colorbg, outline_color=colorfg, font=font

  endfor

  !p.multi = p_multi
  !x.omargin = x_omargin
  !y.omargin = y_omargin

  device, decomposed=1L
  image = tvrd(/true)
  device, decomposed=decomposed

  ;; Save the file:
  ofilename = 'pc_scatter_rho_' + varfile + '.png'
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

end ;;; procedure: pc_scatter_rho
