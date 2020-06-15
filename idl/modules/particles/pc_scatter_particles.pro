pro pc_scatter_particles__height, height, charsize=charsize, $
        npx=npx, npy=npy, width=width, xmargin=xmargin, ymargin=ymargin
  compile_opt hidden, IDL2

  panel_wh = long((width - total(!x.omargin) * charsize * !d.x_ch_size) / double(npx) - $
                  total(xmargin) * charsize * !d.x_ch_size)

  height = long(npy * (panel_wh + total(ymargin) * charsize * !d.y_ch_size) + $
                total(!y.omargin) * charsize * !d.y_ch_size)


  return
end ;;; procedure: pc_scatter_particles__position


pro pc_scatter_particles__colorbar, charsize=charsize, ch2=ch2, charthick=charthick, $
        npx=npx, npy=npy, v_minmax=v_minmax, col_minmax=col_minmax, $
        colorbg=colorbg, colorfg=colorfg, thick=thick, font=font, top=top
  compile_opt hidden, IDL2

  x0 = !x.window[0L]
  x1 = !x.window[1L]
  if keyword_set(top) then begin
    y0 = !y.window[1L] + 2.7d0 * !d.y_ch_size * ch2 / !d.y_size
    y1 = !y.window[1L] + 3.7d0 * !d.y_ch_size * ch2 / !d.y_size
  endif else begin
    y0 = !y.window[0L] - 1.3d0 * !d.y_ch_size * ch2 / !d.y_size
    y1 = !y.window[0L] - 0.3d0 * !d.y_ch_size * ch2 / !d.y_size
  endelse

  bar_x = !d.x_size * (x1 - x0)
  bar_y = ceil(!d.y_size * (y1 - y0))

  xrange = v_minmax
  yrange = [0d0, 1d0]
  bar = byte(lindgen(bar_x) / (bar_x - 1d0) * (col_minmax[1L] - col_minmax[0L]) + col_minmax[0L])

  bar = rebin(bar, bar_x, bar_y)

  px0 = x0 * !d.x_size
  py0 = y0 * !d.y_size

  xstyle = 1 + 4
  ystyle = 1 + 4
  ytickformat = '(a1)'
  plot, [0], /nodata, background=colorbg, charsize=charsize, charthick=charthick, $
        color=colorfg, font=font, /noerase, position=[x0, y0, x1, y1], $
        xrange=xrange, xstyle=xstyle, xthick=thick, $
        yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat

  tv, bar, px0, py0

  tickformat = '(a1)'
  ticklen = 0.3 * !d.y_ch_size * ch2 / !d.y_size
  axis, yaxis=0, charsize=charsize, charthick=charthick, color=colorfg, $
        yminor=0L, yrange=yrange, /ystyle, ythick=thick, ytickformat=tickformat, ticklen=0
  axis, yaxis=1, charsize=charsize, charthick=charthick, color=colorfg, $
        yminor=0L, yrange=yrange, /ystyle, ythick=thick, ytickformat=tickformat, ticklen=0
  axis, xaxis=0, charsize=charsize, charthick=charthick, color=colorfg, $
        xrange=xrange, /xstyle, xthick=thick, ticklen=xticklen, font=font
  axis, xaxis=1, charsize=charsize, charthick=charthick, color=colorfg, $
        xrange=xrange, /xstyle, xthick=thick, ticklen=xticklen, xtickformat=tickformat


  return
end ;;; procedure: pc_scatter_particles__colorbar


pro pc_scatter_particles, pvar, ax=ax, az=az, fraction=fraction_, $
       center=center_, Lxw=Lxw_, Lyw=Lyw_, Lzw=Lzw_, $
       allagr=allagr, agr_index=agr_index, colorbar_top=colorbar_top, $
       charsize=charsize_in, charthick=charthick, colortable=colortable, $
       psym=psym, symsize=symsize, set_font=set_font, t_start=t_start, $
       xsize=xsize, cmxsize=cmxsize, dpi=dpi, thick=thick, font_size=font_size, $
       legend_log_aps=legend_log_aps, legendtext=legendtext, xrange=xrange, $
       ofilename=ofilename, common_range=common_range
  compile_opt IDL2

  if ~n_elements(Lxw_) then Lxw_ = 0L
  if ~n_elements(Lyw_) then Lyw_ = 0L
  if ~n_elements(Lzw_) then Lzw_ = 0L
  allagr = keyword_set(allagr)
  if ~n_elements(charsize_in) then charsize_in = 1d0
  if ~n_elements(colortable) then colortable = 74
  if ~n_elements(psym) then psym = 3
  if ~n_elements(symsize) then symsize = 1.0
  if ~n_elements(xsize) then xsize = 800L
  if ~n_elements(cmxsize) then cmxsize = 13.5
  if n_elements(dpi) eq 1L then xsize = cmxsize * dpi / 2.54d0
  if ~n_elements(t_start) then t_start = 0d0
  if ~n_elements(fraction_) then fraction_ = 1d0
  fraction = fraction_ < 1.0 > 0.0
  pos = strpos(pvar, '.h5', /reverse_search)
  hdf5 = pos eq strlen(pvar) - 3L ? 1L : 0L
  spvar = hdf5 ? file_basename(strmid(pvar, 0L, pos)) : pvar
  if ~n_elements(legend_log_aps) then legend_log_aps = 1L
  if ~keyword_set(common_range) then common_range = 0L
  font = ~n_elements(set_font) ? - 1L : 1L
  tt_font = font eq 1L ? 1L : 0L

  if Lxw_ gt 1L then begin
    str = '_YZ' + strtrim(Lxw_, 2L)
    if ~n_elements(ofilename) then ofilename = 'pc_scatter_particles_' + spvar + str + '.png'
  endif else if Lyw_ gt 1L then begin
    str = '_XZ' + strtrim(Lyw_, 2L)
    if ~n_elements(ofilename) then ofilename = 'pc_scatter_particles_' + spvar + str + '.png'
  endif else if Lzw_ gt 1L then begin
    str = '_XY' + strtrim(Lzw_, 2L)
    if ~n_elements(ofilename) then ofilename = 'pc_scatter_particles_' + spvar + str + '.png'
  endif else begin
    if ~n_elements(ofilename) then ofilename = 'pc_scatter_particles_' + spvar + '_3D.png'
  endelse

  ;; Determine if the input file is regular binary or HDF5:

  pc_read_param, object=ps, /quiet
  pc_read_dim, object=pd, /quiet
  pc_read_ts, object=ts, /quiet

  alpha = ps.rhopmat / ps.rho0 * ps.ap0 / ps.Lx0
  if legend_log_aps then alpha = alog10(alpha)

  pc_read_pvar, object=object, varfile=pvar
  n = object.npar_found

  if ~allagr then begin
    uniq = uniq(object.ap, sort(object.ap))
    agr = object.ap[uniq]
    nagr = ulong(n_elements(temporary(uniq)))
    if n_elements(agr_index) eq 1L then begin
      agr = agr[agr_index]
      nagr = 1UL
    endif else agr_index = 0L
  endif else begin
    agr = 0d0
    nagr = 1UL
    agr_index = 0L
  endelse

  nx = pd.nxgrid  & ny = pd.nygrid & nz = pd.nzgrid
  x_0 = ps.xyz0[0L] & y_0 = ps.xyz0[1L] & z_0 = ps.xyz0[2L]
  x_1 = ps.xyz1[0L] & y_1 = ps.xyz1[1L] & z_1 = ps.xyz1[2L]

  xrange = [x_0, x_1]
  yrange = [y_0, y_1]
  zrange = [z_0, z_1]


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
  for i = 0L, nagr - 1L do vStokes[i] = $
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


  ;;===========================================================================
  ;;===========================================================================
  ;; Setup the plot:

  oldway = 1L
  p_multi = !p.multi
  npx = 1L
  npy = 1L

  if ~allagr then begin
    if nagr ne 1UL then begin
      npx = ceil(sqrt(n_elements(ps.ap0)))
      npy = floor(sqrt(n_elements(ps.ap0)))
      if oldway then !p.multi = [0L, npx, npy, 0L, 0L]
    endif
  endif else begin
    if oldway then !p.multi = 0
  endelse

  charsize_label = charsize_in
  charsize = charsize_in * (npx * npy gt 2L ? 2d0 : 1d0)

  xmargin = [0.3, 0.3]
  ymargin = [4.3, 0.3]
  zmargin = [0.3, 0.3]
  if keyword_set(colorbar_top) then ymargin[1L] += 10.0

  ;; Calculate the window height to get a 1.0 aspect ratio:
  pc_scatter_particles__height, ysize, charsize=charsize_in, npx=npx, npy=npy, $
      width=xsize, xmargin=xmargin, ymargin=ymargin

  device, get_decomposed=decomposed

  plot_3d = 0L
  if ~(Lxw_ gt 1L || Lyw_ gt 1L || Lzw_ gt 1L) then begin
    plot_3d = 1L
    device_name = !d.name
    set_plot, 'z'
    erase
    device, set_pixel_depth=24
    device, set_resolution=[xsize, ysize]
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

  tvlct, red, green, blue, /get
  rgb_table = reform([red, green, blue], 3L, !d.table_size)

  if tt_font then device, set_font=set_font, /tt_font

  salpha = font ? '!9a!x' : '!4a!x'

  ;; Kinetic energy of each particle, one for all, and size-dependent:
  if common_range then begin
    v2 = object.vpx ^ 2 + object.vpy ^ 2 + object.vpz ^ 2
    mima = minmax(v2)
    mimi = lonarr(2L)
    v2 = v2 / max(v2, midx) * (ncolors - 1L)
    color = byte(temporary(v2) + bottom)
    mimi = [bottom, max(color)]
  endif else begin
    color = bytarr(n)
    mima = dblarr(2L, nagr)
    mimi = lonarr(2L, nagr)
    for i = 0UL, nagr - 1UL do begin
      idx = where(object.ap eq agr[i], count)
      if count gt 0L then begin
        v2 = object.vpx[idx] ^ 2 + object.vpy[idx] ^ 2 + object.vpz[idx] ^ 2
        mima[*, i] = minmax(v2)
        v2 = v2 / max(v2, midx) * (ncolors - 1L)
        color[idx] = byte(temporary(v2) + bottom)
        mimi[*, i] = [bottom, max(color[idx])]
      endif
    endfor
  endelse

  if Lxw_ gt 1L || Lyw_ gt 1L || Lzw_ gt 1L then begin

    ;;=================================================================================================
    ;;=================================================================================================
    ;;=================================================================================================
    ;; Slice:
    ;;=================================================================================================
    ;;=================================================================================================
    ;;=================================================================================================

    if Lxw_ gt 1L then begin
      xrange = [y_0, y_1]
      yrange = [z_0, z_1]

      object_c  = object.x
      object_cp = object.xp
      object_xp = object.yp
      object_yp = object.zp

      center = ~n_elements(center_) ? nx / 2L : center_
      center = center < (nx - 1L) > 0L
      Lxw = Lxw_ < (nx - 1L) > 0L
      vi_low = object_c[center - Lxw / 2L]
      vi_hig = object_c[center + Lxw / 2L]

      position = [0.95 * (y_1 - y_0) + y_0, 0.1 * (z_1 - z_0) + z_0]

    endif else if Lyw_ gt 1L then begin

      xrange = [x_0, x_1]
      yrange = [z_0, z_1]

      object_c  = object.y
      object_cp = object.yp
      object_xp = object.xp
      object_yp = object.zp

      center = ~n_elements(center_) ? ny / 2L : center_
      center = center < (ny - 1L) > 0L
      Lyw = Lyw_ < (ny - 1L) > 0L
      vi_low = object_c[center - Lyw / 2L]
      vi_hig = object_c[center + Lyw / 2L]

      position = [0.95 * (x_1 - x_0) + x_0, 0.1 * (z_1 - z_0) + z_0]

    endif else begin

      xrange = [x_0, x_1]
      yrange = [y_0, y_1]

      object_c  = object.z
      object_cp = object.zp
      object_xp = object.xp
      object_yp = object.yp

      center = ~n_elements(center_) ? nz / 2L : center_
      center = center < (nz - 1L) > 0L
      Lzw = Lzw_ < (nz - 1L) > 0L
      vi_low = object_c[center - Lzw / 2L]
      vi_hig = object_c[center + Lzw / 2L]

      position = [0.95 * (x_1 - x_0) + x_0, 0.1 * (y_1 - y_0) + y_0]

    endelse

    xstyle = 1 + 4 & xtickformat = '(a1)'
    ystyle = 1 + 4 & ytickformat = '(a1)'

    if allagr then begin
      plot, [0], /nodata, background=colorbg, charsize=charsize, charthick=charthick, $
            color=colorfg, font=font, $
            xmargin=xmargin, xrange=xrange, xstyle=xstyle, xthick=thick, xtickformat=xtickformat, $
            ymargin=ymargin, yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat

      i = - 1L
      while 1 do begin
        i += long(1d0 / fraction)
        if i ge n then break
        if object_cp[i] lt vi_low || object_cp[i] gt vi_hig then continue
        plots, object_xp[i], object_yp[i], psym=psym, symsize=symsize, color=color[i]
      endwhile

    endif else begin

      current = 0L
      margin = !d.x_ch_size / float(!d.x_size)
      xstyle = 1 + 4
      ystyle = 1 + 4

      format = legend_log_aps ? '(f6.3)' : '(e10.3)'
      object_ap = object.ap
      text = string(bindgen(26) + 97b)
      for j = 0UL, nagr - 1UL do begin

        if oldway then begin
          plot, [0], /nodata, background=colorbg, charsize=charsize, charthick=charthick, $
                color=colorfg, font=font, $
                xmargin=xmargin, xrange=xrange, xstyle=xstyle, xthick=thick, xtickformat=xtickformat, $
                ymargin=ymargin, yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat

          i = - 1L
          while 1 do begin
            i += long(1d0 / fraction)
            if i ge n then break
            if object_cp[i] lt vi_low || object_cp[i] gt vi_hig || object_ap[i] ne agr[j] then continue
            plots, object_xp[i], object_yp[i], psym=psym, symsize=symsize, color=color[i]
          endwhile

          ai = nagr eq 1UL ? agr_index : j

          ;; Stokes number label:
          item = 'St = ' + vStokes[ai]
          al_legend, item, charsize=charsize_label, charthick=charthick, /box, textcolors=colorfg, $
                 /right, background_color=colorbg, outline_color=colorfg, font=font

          item = legend_log_aps ? 'log(' + salpha + ') = ' : salpha + ' = '
          item += strtrim(string(alpha[ai], format=format), 2L)
          item += legend_log_aps ? ', log(a) = ' : ', a = '
          item += strtrim(string(legend_log_aps ? alog10(ps.ap0[ai]) : ps.ap0[ai], format=format), 2L)
          if nagr gt 1UL then item = strmid(text, j, 1L) + ')  ' + item
          al_legend, item, charsize=charsize_label, charthick=charthick, /box, textcolors=colorfg, $
                     /right, /bottom, background_color=colorbg, outline_color=colorfg, font=font

          pc_scatter_particles__colorbar, charsize=charsize, ch2=charsize_in, charthick=charthick, $
              font=font, col_minmax=mimi[*, j], npx=npx, npy=npy, top=colorbar_top, $
              v_minmax=mima[*, j], colorbg=colorbg, colorfg=colorfg, thick=thick

        endif else begin

          ;; Set up the value array at first:
          use_i = lonarr(n)
          i = - 1L & ii = 0L
          while 1 do begin
            i += long(1d0 / fraction)
            if i ge n then break
            if object_cp[i] lt vi_low || object_cp[i] gt vi_hig || object_ap[i] ne agr[j] then continue
            use_i[ii ++] = i
          endwhile
          use_i = use_i[0L : ii - 1L]

          layout = [npx, npy, j + 1L]
          x = object_xp[use_i]
          y = object_yp[use_i]
          magnitude = color[use_i]
          background_color = [255b, 255b, 255b]
          text_color = [0b, 0b, 0b]
          sym_size = 0.04

          rgb_table = 74
          myPlot = scatterplot(x, y, aspect_ratio=1.0, background_color=background_color, $
               current=current, dimensions=[xsize, ysize], layout=layout, magnitude=magnitude, $
               margin=margin, $
               overplot=overplot, rgb_table=rgb_table, symbol='o', /sym_filled, sym_size=sym_size, $
               xrange=xrange, xstyle=xstyle, xthick=thick, xtickformat=xtickformat, $
               yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat)
          current = 1L ;; Only the first plot is not overplotted.

          font_size = 9
          font_style = 'Bold'

          ;; Have to create a color bar to know its position...thereafter,
          ;; it is removed and created anew.
          colorBar = colorbar(target=myPlot, rgb_table=rgb_table)
          position_colorbar = colorBar.position
          colorBar.delete
          colorBar = colorbar(position=position_colorbar, rgb_table=rgb_table, range=mima[*, j], $
                             font_size=font_size, font_style=font_style)

          ai = nagr eq 1UL ? agr_index : j
          item = legend_log_aps ? 'log($\alpha$)=' : '$\alpha$='
          item += strtrim(string(alpha[ai], format=format), 2L)
          item += legend_log_aps ? ', log(a)=' : ', a='
          item += strtrim(string(legend_log_aps ? alog10(ps.ap0[ai]) : ps.ap0[ai], format=format), 2L)
          if nagr gt 1UL then item = strmid(text, j, 1L) + ')  ' + item

          legend = legend(target=myPlot, label=item, font_size=font_size, font_style=font_style, $
                          position=position, /data, text_color=text_color)
          ;                horizontal_alignment='RIGHT', /data, text_color=text_color, $
          ;                vertical_alignment='BOTTOM')

          ;al_legend, item, charsize=charsize, charthick=charthick, /box, textcolors=colorfg, $
          ;           /right, /bottom, background_color=colorbg, outline_color=colorfg
        endelse

      endfor

    endelse

  endif else if nx gt 1L && ny gt 1L && nz gt 1L then begin

    ;;=================================================================================================
    ;;=================================================================================================
    ;;=================================================================================================
    ;; 3D scatter plot:
    ;;=================================================================================================
    ;;=================================================================================================
    ;;=================================================================================================

    object_cp = object.xp
    object_xp = object.xp
    object_yp = object.yp
    object_zp = object.zp

    xtickformat = '(a1)'
    ytickformat = '(a1)'
    ztickformat = '(a1)'

    format = legend_log_aps ? '(f6.3)' : '(e10.3)'
    object_ap = object.ap
    text = string(bindgen(26) + 97b)
    for j = 0UL, nagr - 1UL do begin

      surface, [[0.0,0.0,0.0], [0.0,0.0,0.0]], /nodata, ax=ax, az=az, $
               background=colorbg, charsize=charsize, charthick=charthick, color=colorfg, $
               font=font, /save, $
               xmargin=xmargin, xrange=xrange, /xstyle, xthick=thick, xtickformat=xtickformat, $
               ymargin=ymargin, yrange=yrange, /ystyle, ythick=thick, ytickformat=ytickformat, $
               zmargin=zmargin, zrange=zrange, /zstyle, zthick=thick, ztickformat=ztickformat

      axis, xaxis=1, x_0, y_1, z_1, color=colorfg, charsize=charsize, charthick=charthick, $
            /t3d, xthick=thick, xtickformat=xtickformat
      axis, xaxis=1, x_0, y_1, z_0, color=colorfg, charsize=charsize, charthick=charthick, $
            /t3d, xthick=thick, xtickformat=xtickformat
      axis, yaxis=1, x_1, y_0, z_0, color=colorfg, charsize=charsize, charthick=charthick, $
            /t3d, ythick=thick, ytickformat=ytickformat
      axis, yaxis=1, x_1, y_0, z_1, color=colorfg, charsize=charsize, charthick=charthick, $
            /t3d, ythick=thick, ytickformat=ytickformat
      axis, zaxis=0, x_1, y_1, z_0, color=colorfg, charsize=charsize, charthick=charthick, $
            /t3d, zthick=thick, ztickformat=ztickformat
      axis, zaxis=0, x_1, y_0, z_0, color=colorfg, charsize=charsize, charthick=charthick, $
            /t3d, zthick=thick, ztickformat=ztickformat

      i = - 1L
      while 1 do begin
        i += long(1d0 / fraction)
        if i ge n then break
        if object_ap[i] ne agr[j] then continue
        plots, object_xp[i], object_yp[i], object_zp[i], /t3d, $
               psym=psym, symsize=symsize, color=color[i]
      endwhile

      axis, zaxis=1, x_0, y_0, z_0, /t3d, color=colorfg, charsize=charsize, charthick=charthick, $
            zthick=thick, ztickformat=ztickformat
      axis, yaxis=0, x_0, y_0, z_1, /t3d, color=colorfg, charsize=charsize, charthick=charthick, $
            ythick=thick, ytickformat=ytickformat
      axis, xaxis=0, x_0, y_0, z_1, /t3d, color=colorfg, charsize=charsize, charthick=charthick, $
            xthick=thick, xtickformat=xtickformat

      ai = nagr eq 1UL ? agr_index : j

      ;; Stokes number label:
      item = 'St = ' + vStokes[ai]
      al_legend, item, charsize=charsize_label, charthick=charthick, /box, textcolors=colorfg, $
                 /right, background_color=colorbg, outline_color=colorfg, font=font

      ;; Alpha and a values:
      item = legend_log_aps ? 'log(' + salpha + ') = ' : salpha + ' = '
      item += strtrim(string(alpha[ai], format=format), 2L)
      item += legend_log_aps ? ', log(a) = ' : ', a = '
      item += strtrim(string(legend_log_aps ? alog10(ps.ap0[ai]) : ps.ap0[ai], format=format), 2L)
      if nagr gt 1UL then item = strmid(text, j, 1L) + ')  ' + item
      al_legend, item, charsize=charsize_label, charthick=charthick, /box, textcolors=colorfg, $
                 /right, /bottom, background_color=colorbg, outline_color=colorfg, font=font

      ;; Color bar:
      pc_scatter_particles__colorbar, charsize=charsize, ch2=charsize_in, charthick=charthick, $
          font=font, col_minmax=mimi[*, j], npx=npx, npy=npy, top=colorbar_top, $
          v_minmax=mima[*, j], colorbg=colorbg, colorfg=colorfg, thick=thick

    endfor ;; j = 0UL, nagr - 1UL

  endif

  if plot_3d then begin
    image = tvrd(/true)
    size = size(image, /dimensions)
    xsize_2 = size[1L]
    ysize_2 = size[2L]
    set_plot, device_name
    device, decomposed=decomposed
  endif else begin
    copy = [0L, 0L, xsize, ysize, 0L, 0L, winid]
    xsize_2 = xsize
    ysize_2 = ysize
  endelse

  window, xsize=xsize_2, ysize=ysize_2
  if plot_3d then begin
    tv, image, /true
  endif else if oldway then begin
    wset, !d.window
    device, copy=copy
    wdelete, winid
    device, decomposed=1L
    image = tvrd(/true)
    device, decomposed=decomposed
  endif

  ;; Save the file:
  if oldway then begin
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
  endif else begin
    resolution = 300L
    MyPlot.save, ofilename, /bitmap, resolution=resolution, /compression
  endelse

  if ~allagr && oldway then !p.multi = p_multi
  device, decomposed=decomposed


  return
end ;;; procedure: pc_scatter_particles
