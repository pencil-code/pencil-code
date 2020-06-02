function pc_power_particles_range, i, n=n, k=k, specptr=specptr, xrange=xrange
  compile_opt hidden, IDL2

  for ii = 0UL, n - 1UL do begin
    idx = where(k gt xrange[0L] and k le xrange[1L], count)
    min = min((*specptr[ii])[idx, i], max=max)
    yrange = ~ii ? [min, max] : $
             [min([yrange[0L], min]), max([yrange[1L], max])]
  endfor


  return, yrange
end ;;; function: pc_power_particles_range


function pc_power_particles_range_er, i, n=n, k=k, spec=spec, xrange=xrange, div=div
  compile_opt hidden, IDL2

  idx = where(k gt xrange[0L] and k le xrange[1L], count)
  if ~n_elements(div) then begin
    min = min(spec[idx, i]           , max=max)
  endif else begin
    min = min(spec[idx, i] * div[idx], max=max)
  endelse
  yrange = [min, max]


  return, yrange
end ;;; function: pc_power_particles_range_er


pro pc_power_dimension_annotate, k, data, text, left=left, $
                                 charsize=charsize, charthick=charthick
  compile_opt hidden, IDL2

  left = keyword_set(left)

  dpx = charsize * !d.x_ch_size / !d.x_size
  dpy = charsize * !d.y_ch_size / !d.y_size

  ;; Place text just above lower axis:
  value = 0.03 * (!y.crange[1L] - !y.crange[0L]) + !y.crange[0L]

  ;; Locate data point near this position:
  min = min(abs(alog10(data) - value), idx, /nan)

  pos = convert_coord(k[idx], 1d1 ^ value, /data, /to_normal)
  ypos = pos[1L]

  ;; Put the text to the left or the right of the point:
  if keyword_set(left) then begin
    xpos = pos[0L] - 1.0 * dpx
    alignment = 1.0
  endif else begin
    xpos = pos[0L] + 0.2 * dpx
    alignment = 0.0
  endelse
  xyouts, xpos, ypos, text, alignment=alignment, /normal, $
          charsize=charsize, charthick=charthick, color=colorfg

  return
end ;;; procedure: pc_power_dimension_annotate

pro pc_power_particles, noekin=noekin, norho=norho, noomega=noomega, $
        t_start=t_start, t_relaxed=t_relaxed, $
        stride=stride, noloop=noloop, snapshot_i=snapshot_i, $
        urms_relaxed=urms_relaxed, sleep=sleep, charsize=charsize, $
        charthick=charthick, thick=thick, postscript=postscript, $
        psfilename=psfilename, shift_yrange=shift_yrange, $
        scalefactor=scalefactor, separate_yrange=separate_yrange_, $
        xrange=xrange, xsize=xsize, xtitle=xtitle, psxsize=psxsize, $
        yrange=yrange, ysize=ysize, ytitle=ytitle, psysize=psysize, $
        legend_linsize=legend_linsize_, legend_spacing=legend_spacing_, $
        legend_top=legend_top, legend_log_aps=legend_log_aps, $
        legend_short=legend_short, legend_text=legend_text, $
        leftann_ekin=left_ekin, leftann_rho=left_lnrho, leftann_omega=left_omega
  compile_opt IDL2

  ;on_error, 2

  if ~n_elements(t_start) then t_start = 0d0
  if ~n_elements(t_relaxed) then t_relaxed = t_start
  if ~n_elements(stride) then stride = 1L
  noloop = keyword_set(noloop)
  if ~n_elements(sleep) then sleep = 0.1
  if ~n_elements(postscript) then postscript = 1L
  if ~n_elements(psfilename) then psfilename = 'power.eps'
  if ~n_elements(charsize) then charsize = 1.0
  if ~n_elements(charthick) then charthick = 1.0
  if ~n_elements(thick) then thick = 1.0
  if ~n_elements(psxsize) then psxsize = 20.0
  if ~n_elements(psysize) then psysize = 17.0
  if ~n_elements(shift_yrange) then shift_yrange = 0d0
  if ~n_elements(scalefactor) then scalefactor = 1d0
  separate_yrange = keyword_set(separate_yrange_)
  separate_yrange_shift = shift_yrange ne 0d0 ? 1L : 0L
  separate_yrange = keyword_set(separate_yrange_) or  shift_yrange ne 0d0
  scalefactor_use = separate_yrange ? 1d0 : scalefactor
  if ~n_elements(xtitle) then xtitle = 'k/k!d1!n'
  legend_spacing = charsize * (~n_elements(legend_spacing_) ? 1.5 : legend_spacing_)
  legend_linsize = ~n_elements(legend_linsize_) ? 0.4 : legend_linsize_
  legend_top = keyword_set(legend_top)
  if ~n_elements(legend_log_aps) then legend_log_aps = 1L
  legend_short = ~n_elements(legend_short) ? 1L : 0L

  ;; Color indices:
  colorbg = 1b
  colorfg = 0b
  colorrd = 2b
  colorgn = 4b
  colorbl = 5b
  coloror = 6b
  colorpr = 7b
  colorlb = 8b
  coloryl = 9b

  na_colors = [colorbl, colorlb, colorgn, colorpr, coloror]
  
  plot_kin = 0L
  plot_lnrho = 0L
  plot_omega = 0L
  plot_na = 0L

  filename_kin = '_kin'
  filename_kin_full = 'data' + path_sep() + 'power' + filename_kin + '.dat'
  if file_test(filename_kin_full, /read, /regular) && ~keyword_set(noekin) then begin
    power, v1=filename_kin, k=k, spec1=spec_kin, i=n, tt=t, /noplot
    ytitle_s = 'E!dRHO!n(k)'
    plot_kin = 1L
  endif

  filename_lnrho = '_lr'
  filename_lnrho_full = 'data' + path_sep() + 'power' + filename_lnrho + '.dat'
  if file_test(filename_lnrho_full, /read, /regular) && ~keyword_set(norho) then begin
    power, v1=filename_lnrho, k=k, spec1=spec_lr, i=n, tt=t, /noplot
    ytitle_s_tmp = 'P!dRHO!n'
    if ~n_elements(ytitle_s) then begin
      ytitle_s = temporary(ytitle_s_tmp)
    endif else begin
      ytitle_s += ' :: ' + temporary(ytitle_s_tmp)
    endelse
    plot_lnrho = 1L
  endif

  filename_omega = 'o'
  filename_omega_full = 'data' + path_sep() + 'power' + filename_omega + '.dat'
  if file_test(filename_omega_full, /read, /regular) && ~keyword_set(noomega) then begin
    power, v1=filename_omega, k=k, spec1=spec_o, i=n, tt=t, /noplot
    ytitle_s_tmp = 'P!dOMEGA!n'
    if ~n_elements(ytitle_s) then begin
      ytitle_s = temporary(ytitle_s_tmp)
    endif else begin
      ytitle_s += ' :: ' + temporary(ytitle_s_tmp)
    endelse
    plot_omega = 1L
  endif

  filename_na1 = file_search('data' + path_sep() + 'power_na-?.dat')
  if filename_na1[0L] ne '' then begin
    filename_na = temporary(filename_na1)

    filename_na2 = file_search('data' + path_sep() + 'power_na-??.dat')
    if filename_na2[0L] ne '' then $
       filename_na = [filename_na, temporary(filename_na2)]

    ;; Only keep the suffix:
    filename_na = file_basename(filename_na)
    pos = strpos(filename_na, '.')
    for i = 0UL, n_elements(filename_na) - 1UL do $
       filename_na[i] = strmid(filename_na[i], 5L, pos[i] - 5L)

    spec_na = ptrarr(n_elements(filename_na))
    for i = 0UL, n_elements(filename_na) - 1UL do begin
      power, v1=filename_na[i], k=k, spec1=spec_na_tmp,i=n,tt=t,/noplot
      spec_na[i] = ptr_new(spec_na_tmp)
    endfor

    if separate_yrange then begin
      ytitle_na = 'P!dn!ia,i!n'
    endif else begin
      if scalefactor_use ne 1d0 then begin
        mantissa = floor(alog10(scalefactor_use))
        remainder = scalefactor_use / (1d1 * mantissa)
        str = string(remainder, format=format) + string(215b) + $
              '10!u' + strtrim(string(mantissa, format='(i3)'), 2L) + '!n '
        ytitle_s += ' :: ' + str + string(215b) + ' n!dp,i!n'
      endif
    endelse

    ;; For the k^2 guide:
    k_l1_na = 1d1 & k_l2_na = 4d1

    plot_na = 1L
  endif

  k_l1_nn = 2.5d0 & k_l2_nn = 9d0

  if ~plot_lnrho && ~plot_kin && ~plot_omega && ~plot_na then $
     message, 'There are no power spectra to plot.'

  
  ;; Read model parameters:
  pc_read_param, object=ps
  pc_read_param, object=ps2, /param2
  pc_read_ts, object=ts

  alpha = ps.rhopmat / ps.rho0 * ps.ap0 / ps.Lx0
  if legend_log_aps then alpha = alog10(alpha)

  digits = legend_log_aps ? 1UL : 1UL
  f_done = 0UL
  repeat begin
    if legend_log_aps then begin
      format = '(f' + strtrim(digits + 4UL, 2L) + '.' + strtrim(digits, 2L) + ')'
    endif else begin
      format = '(e' + strtrim(digits + 7UL, 2L) + '.' + strtrim(digits, 2L) + ')'
    endelse
    aps = strtrim(string(alpha, format=format), 2L)
    idx = uniq(aps, sort(aps))
    if n_elements(idx) eq n_elements(alpha) or digits ge (legend_log_aps ? 6UL : 4UL) then begin
      f_done = 1UL
    endif else begin
      digits ++
    endelse
  endrep until f_done


  ;; Check for the existence of the "epsk" variable:
  epsk_available = 0L
  tag_names = strlowcase(tag_names(ts))
  idx = where(~strpos(tag_names, 'epsk'), count)
  if count eq 1L then epsk_available = 1L

  ;;spawn,'if (! -f param.pro) touch parameters.pro'
  ;;@parameters

  if ~n_elements(xrange) then xrange = [1d0, max(k)]


  k1 = 1d0 / k


  ;; Finding the starting index of the time interval to use:
  min = min(abs(t - t_start), i_0)

  ;; Final index:
  i_1 = n - 2L

  ;; Determine which snapshot to use:
  i_start = ~noloop ? i_0 : i_1
  i_end = i_1
  if n_elements(snapshot_i) eq 1L then begin
    i_start = snapshot_i
    i_end = snapshot_i
  endif

  i_b = 0L
  if n_elements(urms_relaxed) ne 0L then begin
    min = min(abs(ts.urms - urms_relaxed), idx)
    if idx ne 0L then i_b = idx
  endif

  if n_elements(t_relaxed) ne 0L then begin
    min = min(abs(t - t_relaxed), idx)
    if idx ne 0L then i_b = idx
  endif

  min_urms = min(ts.urms[i_b : *], max=max_urms)


  ;;===========================================================================
  ;;===========================================================================
  ;; Time scales:
  
  ;; Irrotational forcing: stopping time:
  tau_stop = sqrt(!dpi / 8d0) * ps.rhopmat * ps.ap0 / (ps.rho0 * ps.cs0)
  str_tau_stop = strtrim(string(tau_stop, format='(e9.2)'), 2L)

  ;; Vortical forcing: mean turnover time:
  tau_turnover = 1d0 / mean(ts.urms[i_b : *])
  str_tau_turnover = strtrim(string(tau_turnover, format='(e9.2)'), 2L)

  ;;tau_stop = ps.rhopmat * 4.0 * ps.ap0 ^ 2 / (18d0 * ps.rho0 * ps2.nu)

  ;; Kolmogorov:
  if epsk_available then begin
    tau_Kolmogorov = (ps.rho0 * ps2.nu ^ 3 / mean(ts.epsK[i_b : *])) ^ 0.25d0
    str_tau_Kolmogorov = strtrim(string(tau_Kolmogorov, format='(e9.2)'), 2L)
  endif

  ;; Stokes; grain-size dependent:
  vStokes = strarr(n_elements(ps.ap0))
  for i = 0L, n_elements(ps.ap0) - 1L do vStokes[i] = $
     strtrim(string(tau_stop[i] / tau_turnover, format='(e9.2)'), 2L)


  ;;===========================================================================
  ;;===========================================================================
  ;; Set up colors and device-dependent titles:

  if ~(noloop && postscript) then begin
    bottom = 10L
    ncolors = !d.table_size - bottom
    loadct, 74, bottom=bottom, ncolors=ncolors
    tvlct,   0b,   0b ,  0b, colorfg ;; Black
    tvlct, 255b, 255b, 255b, colorbg ;; White
    tvlct, 213b,  94b,   0b, colorrd ;; Vermillion
    tvlct,   0b, 158b, 115b, colorgn ;; Bluish green
    tvlct,   0b, 114b, 178b, colorbl ;; Blue
    tvlct, 230b, 159b,   0b, coloror ;; Orange
    tvlct, 204b, 121b, 167b, colorpr ;; Purple (Mulberry)
    tvlct,  86b, 180b, 233b, colorlb ;; Sky Blue
    tvlct, 240b, 228b,  66b, coloryl ;; Yellow

    ;; Y title:
    pos = strpos(ytitle_s, 'OMEGA')
    repeat begin
      if pos ge 0L then ytitle_s = strmid(ytitle_s, 0L, pos) + $
                                   '!4X!x' + strmid(ytitle_s, pos + 5L)
      pos = strpos(ytitle_s, 'OMEGA')
    endrep until pos eq - 1L

    pos = strpos(ytitle_s, 'RHO')
    repeat begin
      if pos ge 0L then ytitle_s = strmid(ytitle_s, 0L, pos) + $
                                   '!4q!x' + strmid(ytitle_s, pos + 3L)
      pos = strpos(ytitle_s, 'RHO')
    endrep until pos eq - 1L

  endif else begin

    ;; Y title:
    pos = strpos(ytitle_s, 'OMEGA')
    repeat begin
      if pos ge 0L then ytitle_s = strmid(ytitle_s, 0L, pos) + $
                                   '!9w!x' + strmid(ytitle_s, pos + 5L)
      pos = strpos(ytitle_s, 'OMEGA')
    endrep until pos eq - 1L

    pos = strpos(ytitle_s, 'RHO')
    repeat begin
      if pos ge 0L then ytitle_s = strmid(ytitle_s, 0L, pos) + $
                                   '!9r!x' + strmid(ytitle_s, pos + 3L)
      pos = strpos(ytitle_s, 'RHO')
    endrep until pos eq - 1L

  endelse

  ytitle = temporary(ytitle_s)


  ;;===========================================================================
  ;;===========================================================================
  ;; Loop over the spectra:

  scaling = 2.0

  if ~(noloop && postscript) then begin
    window, /free, xsize=xsize, ysize=ysize, /pixmap
    win_id_pixmap = !d.window
    window, /free, xsize=xsize, ysize=ysize
    win_id = !d.window

    if ~n_elements(xsize) then xsize = !d.x_size
    if ~n_elements(ysize) then ysize = !d.y_size
    copy = [0L, 0L, xsize, ysize, 0L, 0L, win_id_pixmap]
  endif

  for i = i_start, i_end, stride do begin
 
    if noloop && postscript then begin
      pd = !d.name
      set_plot, 'ps'

      pos = strpos(psfilename, '.eps', /reverse_search)
      fils = strmid(psfilename, 0L, pos)
      file = strmid(psfilename, pos)
      psfilename_use = fils + '_' + strtrim(i, 2L) + file

      device, filename=psfilename_use, /isolatin1, /color, bits_per_pixel=8, $
              /encapsulated, /portrait, xsize=psxsize, ysize=psysize, font_size=10
      pf = !p.font
      !p.font = 0L

      bottom = 10L
      ncolors = !d.table_size - bottom
      loadct, 74, bottom=bottom, ncolors=ncolors
      tvlct,   0b,   0b ,  0b, colorfg ;; Black
      tvlct, 255b, 255b, 255b, colorbg ;; White
      tvlct, 213b,  94b,   0b, colorrd ;; Vermillion
      tvlct,   0b, 158b, 115b, colorgn ;; Bluish green
      tvlct,   0b, 114b, 178b, colorbl ;; Blue
      tvlct, 230b, 159b,   0b, coloror ;; Orange
      tvlct, 204b, 121b, 167b, colorpr ;; Purple (Mulberry)
      tvlct,  86b, 180b, 233b, colorlb ;; Sky Blue
      tvlct, 240b, 228b,  66b, coloryl ;; Yellow

      staus = '!9t!x!ds!n='
      staut = '!9t!x!dt!n='
      stKol = '!9t!x!dK!n='
      sStokes = ', St='
      sdens = '!9r!x'
      sdens = '!9r!x'
      somega = '!9w!x'
      salpha = '!9a!x'

    endif else begin

      sdens = '!4q!x'
      somega = '!4X!x'
      salpha = '!4a!x'

      charthick_use = charthick
      charsize_use = charsize

      staus = '!4s!x!ds!n='
      staut = '!4s!x!dt!n='
      stKol = '!4s!x!dK!n='
      sStokes = ', St='

      wset, win_id_pixmap
    endelse

    dpx = charsize * !d.x_ch_size / !d.x_size
    dpy = charsize * !d.y_ch_size / !d.y_size

    n_spec_na = n_elements(spec_na)
    icolors = intarr(n_spec_na + 1L)
    ilinestyle = intarr(n_spec_na + 1L) - 1L
    iitems = strarr(n_spec_na + 1L)
    linestyle = 2
    for ii = 0UL, n_spec_na - 1UL do begin
      icolors[ii] = na_colors[ii / 5L]
      iitems[ii] = (legend_log_aps ? 'log(' : '') + salpha + ')=' + aps[ii]
      if ~legend_short then $
         iitems[i] += ', ' + staus + str_tau_stop[ii] + sStokes + vStokes[ii]
      ilinestyle[ii] = linestyle

      linestyle ++
      if linestyle eq 6 then linestyle = 1
    endfor

    title = 't=' + strtrim(string(t[i], format='(f10.3)'), 2L) + $
            ' :: ' + staut + str_tau_turnover
    if epsk_available then title += ' :: ' + stKol + str_tau_Kolmogorov

    yrange_nn = [1d50, -1d50]
    ki2r = 1d0 / k ^ 2
    if plot_kin then begin
      xrange_ = [k_l1_nn, k_l2_nn]
      yrange_ne  = pc_power_particles_range_er(i, k=k, spec=spec_kin, xrange=xrange_, div=ki2r)
      yrange_nn = [min([yrange_nn[0L], yrange_ne[0L]]), max([yrange_nn[1L], yrange_ne[1L]])]
    endif
    if plot_lnrho then begin
      xrange_ = [k_l1_nn, k_l2_nn]
      yrange_nr  = pc_power_particles_range_er(i, k=k, spec=spec_lr  , xrange=xrange_)
      yrange_nn = [min([yrange_nn[0L], yrange_nr[0L]]), max([yrange_nn[1L], yrange_nr[1L]])]
    endif
    if plot_omega then begin
      xrange_ = [k_l1_nn, k_l2_nn]
      yrange_no  = pc_power_particles_range_er(i, k=k, spec=spec_o  , xrange=xrange_, div=ki2r)
      yrange_nn = [min([yrange_nn[0L], yrange_no[0L]]), max([yrange_nn[1L], yrange_no[1L]])]
    endif

    if plot_na then begin
      yrange_na  = pc_power_particles_range(i, n=n_spec_na, $
                                            k=k, specptr=spec_na, xrange=xrange)
      xrange_ = [k_l1_na, k_l2_na]
      yrange_nax = pc_power_particles_range(i, n=n_spec_na, $
                                            k=k, specptr=spec_na, xrange=xrange_)
      if ~n_elements(yrange) then $
         yrange = [min([3d-12, yrange_na[0L]]), $
                   max([3d-2 , yrange_na[1L]])]
    endif else begin
      if ~n_elements(yrange) then yrange = [3d-12, 3d-2]
    endelse

    yminor = 9
    xstyle = 1
    ystyle = plot_na && separate_yrange ? 9 : 1
    xmargin = plot_na && separate_yrange ? [8.0, 8.0] : [8.0, 0.3]
    ymargin = [3.3, 0.3]
    plot, k, /nodata, background=colorbg, color=colorfg, $
          charsize=charsize, charthick=charthick, $
          /xlog, xmargin=xmargin, xrange=xrange, xstyle=xstyle, $
          xthick=thick, xtitle=xtitle, $
          /ylog, ymargin=ymargin, yminor=yminor, yrange=yrange, ystyle=ystyle, $
          ythick=thick, ytitle=ytitle

    if ~legend_short then begin
      xpos = 0.50d0 * (!x.window[1L] - !x.window[0L]) + !x.window[0L]
      ypos = 0.95d0 * (!y.window[1L] - !y.window[0L]) + !y.window[0L]
      xyouts, xpos, ypos, title, charsize=charsize, charthick=charthick, $
              color=colorfg, alignment=0.5, /normal
    endif


    ;;=========================================================================
    ;;=========================================================================
    ;; Guides:

    if plot_kin || plot_lnrho || plot_omega then begin
      ;; Guide: 1 / k^2
      k_arr = [k_l1_nn, k_l2_nn]
      y_arr = k_arr ^ (- 2d0) / k_l2_nn ^ (- 2d0) * yrange_nn[0L] * 0.5
      oplot, k_arr, y_arr, linestyle=0, color=coloror, thick=thick * 2
  
      ;; Annotation:
      xpos = k_arr[0L]
      ypos = y_arr[0L]
      pos = convert_coord(xpos, ypos, /data, /to_normal)
      xpos = pos[0L] - 0.25 * dpx
      ypos = pos[1L] - 0.25 * dpy
      alignment = 1.0
      text = 'k!u-2!n'
      xyouts, xpos, ypos, text, alignment=alignment, /normal, $
              charsize=charsize, charthick=charthick, color=colorfg

      ;; Guide: 1 / k^(5/3)
      y_arr = k_arr ^ (- 5d0 / 3d0) / k_l2_nn ^(- 5d0 / 3d0) * yrange_nn[0L] * 0.5
      oplot, k_arr, y_arr, linestyle=0, color=colorgn, thick=thick * 2

      ;; Annotation:
      xpos = 1d1 ^ (0.25 * (alog10(k_arr[1L]) - alog10(k_arr[0L])) + alog10(k_arr[0L]))
      ypos = y_arr[0L]
      pos = convert_coord(xpos, ypos, /data, /to_normal)
      xpos = pos[0L]
      ypos = pos[1L] - 2.5 * dpy
      alignment = 0.5
      text = 'k!u-5/3!n '
      xyouts, xpos, ypos, text, alignment=alignment, /normal, $
              charsize=charsize, charthick=charthick, color=colorfg

    endif

    if plot_kin then begin
      ekin = spec_kin[*, i] / k ^ 2
      oplot, k, ekin, color=colorfg, thick=thick, linestyle=3

      pc_power_dimension_annotate, k, ekin, 'E(k)', $
          left=keyword_set(left_ekin), charsize=charsize, charthick=charthick
    endif

    if plot_omega then begin
      omega = spec_o[*, i] / k ^ 2
      oplot, k, omega, color=colorpr, thick=thick

      pc_power_dimension_annotate, k, omega, 'P!d' + somega + '!n', $
          left=keyword_set(left_omega), charsize=charsize, charthick=charthick
    endif

    if plot_lnrho then begin
      oplot, k, spec_lr[*, i], color=colorlb, thick=thick

      pc_power_dimension_annotate, k, spec_lr[*, i], 'P!d' + sdens + '!n', $
          left=keyword_set(left_lnrho), charsize=charsize, charthick=charthick
    endif

    ;if plot_kin then begin
    ;  l_items = ['E!d' + sdens + '!n(k)']
    ;  l_colors = [colorfg]
    ;  l_linestyle = [2]
    ;  l_thick = thick * [1.0]
    ;endif

    ;if plot_lnrho then begin
    ;  l_items = ~n_elements(l_items) ? ['P!d' + sdens + '!n'] : $
    ;                                   [l_items, 'P!d' + sdens + '!n']
    ;  l_colors = ~n_elements(l_colors) ? [colorlb] : [l_colors, colorlb]
    ;  l_linestyle = ~n_elements(l_linestyle) ? [0] : [l_linestyle, 0]
    ;  l_thick = ~n_elements(l_thick) ? [1.0] : [l_thick, 1.0]
    ;endif
    ;
    ;if plot_omega then begin
    ;  l_items = ~n_elements(l_items) ? ['P!d' + somega + '!n'] : $
    ;                                   [l_items, 'P!d' + somega + '!n']
    ;  l_colors = ~n_elements(l_colors) ? [colorpr] : [l_colors, colorpr]
    ;  l_linestyle = ~n_elements(l_linestyle) ? [0] : [l_linestyle, 0]
    ;  l_thick = ~n_elements(l_thick) ? [1.0] : [l_thick, 1.0]
    ;endif

    ;if plot_kin || plot_lnrho then begin
    ;  l_items = [l_items, 'k!u-2!n', 'k!u5/3!n']
    ;  l_colors = [l_colors, coloror, colorgn]
    ;  l_linestyle = [l_linestyle, 0, 0]
    ;  l_thick = [l_thick, 2.0, 2.0]
    ;endif

    ;if plot_na then begin
    ;  l_items = ~n_elements(l_items) ? 'k!u2!n' : [l_items, 'k!u2!n']
    ;  l_colors = [l_colors, colorrd]
    ;  l_linestyle = [l_linestyle, 0]
    ;  l_thick = [l_thick, 2.0]
    ;endif

    ;if legend_short then begin
    ;  l_items = [l_items, iitems]
    ;  l_colors = [l_colors, icolors]
    ;  l_linestyle = [l_linestyle, ilinestyle]
    ;  l_thick = [l_thick, fltarr(n_elements(iitems)) + thick]
    ;endif
    al_legend, iitems, linestyle=ilinestyle, top=legend_top, /right, box=0, $
               charsize=charsize, charthick=charthick, color=icolors, $
               thick=thick, linsize=legend_linsize, $
               textcolors=colorfg, spacing=legend_spacing

    if n_elements(legend_text) eq 1L then begin
      al_legend, legend_text, /bottom, /left, textcolors=colorfg, $
                 charsize=charsize, charthick=charthick, box=0
    endif


    if plot_na then begin

      if separate_yrange then begin

        yrange_yaxis = separate_yrange_shift ? yrange * shift_yrange[0L] : yrange_na
        axis, charsize=charsize, charthick=charthick, color=colorfg, /yaxis, /ylog, yminor=yminor, $
              /ystyle, yrange=yrange_yaxis, ythick=thick, ytitle=ytitle_na, /save

      endif


      ;;=======================================================================
      ;;=======================================================================
      ;; Guide: k^2

      k_arr = [k_l1_na, k_l2_na]
      y_arr = k_arr ^ 2 / k_l1_na ^ 2 * yrange_nax[0L] * 0.5
      oplot, k_arr, y_arr, linestyle=0, color=colorrd, thick=thick * 2

      ;; Annotation:
      xpos = 0.2 * (k_arr[1L] - k_arr[0L]) + k_arr[0L]
      ypos = y_arr[0L]
      pos = convert_coord(xpos, ypos, /data, /to_normal)
      xpos = pos[0L]
      ypos = pos[1L] - 0.5 * dpy
      alignment = 1.0
      text = 'k!u2!n '
      xyouts, xpos, ypos, text, alignment=alignment, /normal, $
              charsize=charsize, charthick=charthick, color=colorfg


      ;;=======================================================================
      ;;=======================================================================
      ;; Plot particle density spectra:

      linestyle = 2
      for ii = 0UL, n_spec_na - 1UL do begin
        color = na_colors[ii / 5L]
        oplot, k, scalefactor_use * (*spec_na[ii])[*, i], $
               linestyle=linestyle, color=color, thick=thick

        linestyle ++
        if linestyle eq 6 then linestyle = 1
      endfor

      if ~legend_short then $
         al_legend, iitems, linestyle=ilinestyle, /bottom, box=0, thick=thick, $
                    charsize=charsize, charthick=charthick, color=icolors, $
                    linsize=legend_linsize, textcolors=colorfg, spacing=legend_spacing

    endif ;; plot_na


    if noloop && postscript then begin

      !p.font = pf
      device, /close_file
      set_plot, pd

    endif else begin

      ;; Copy plotted image to shown window:
      wset, win_id
      device, copy=copy

      print, i, t[i]
      wait, sleep

    endelse

  endfor ;; i = i_start, i_end, stride

  if plot_na then begin
    for i = 0UL, n_elements(spec_na) - 1UL do begin
      if ptr_valid(spec_na[i]) then ptr_free, spec_na[i]
    endfor
  endif
  
  if ~(noloop && postscript) then wdelete, win_id_pixmap

end ;; procedure: pc_power
