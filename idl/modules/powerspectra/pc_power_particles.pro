pro pc_power_particles, t_start=t_start, t_relaxed=t_relaxed, $
        stride=stride, noloop=noloop, snapshot_i=snapshot_i, $
        urms_relaxed=urms_relaxed, sleep=sleep, charsize=charsize, $
        charthick=charthick, thick=thick, postscript=postscript, $
        psfilename=psfilename, shift_yrange=shift_yrange, $
        scalefactor=scalefactor, separate_yrange=separate_yrange_, $
        xrange=xrange, xsize=xsize, xtitle=xtitle, psxsize=psxsize, $
        yrange=yrange, ysize=ysize, ytitle=ytitle, psysize=psysize
  compile_opt IDL2

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
  if shift_yrange ne 0d0 then separate_yrange_shift = 1L
  scalefactor_use = separate_yrange ? 1d0 : scalefactor
  if ~n_elements(xtitle) then xtitle = 'k/k!d1!n'

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
  if file_test(filename_kin_full, /read, /regular) then begin
    power, v1=filename_kin, k=k, spec1=spec_kin, i=n, tt=t, /noplot
    ytitle_s = 'E!dRHO!n(k)'
    plot_kin = 1L
  endif

  filename_lnrho = '_lr'
  filename_lnrho_full = 'data' + path_sep() + 'power' + filename_lnrho + '.dat'
  if file_test(filename_lnrho_full, /read, /regular) then begin
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
  if file_test(filename_omega_full, /read, /regular) then begin
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

    plot_na = 1L
  endif

  if ~plot_lnrho && ~plot_kin && ~plot_omega && ~plot_na then $
     message, 'There are no power spectra to plot.'

  
  ;; Read model parameters:
  pc_read_param, object=ps
  pc_read_param, object=ps2, /param2
  pc_read_ts, object=ts

  aps = strtrim(string(ps.ap0, format='(e10.3)'), 2L)

  ;; Check for the existence of the "epsk" variable:
  epsk_available = 0L
  tag_names = strlowcase(tag_names(ts))
  idx = where(~strpos(tag_names, 'epsk'), count)
  if count eq 1L then epsk_available = 1L

  ;;spawn,'if (! -f param.pro) touch parameters.pro'
  ;;@parameters

  if ~n_elements(xrange) then xrange = [1d0, max(k)]
  if ~n_elements(yrange) then yrange = [3d-12, 3d-2]

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
                                   '!9W!x' + strmid(ytitle_s, pos + 5L)
      pos = strpos(ytitle_s, 'OMEGA')
    endrep until pos eq - 1L

    pos = strpos(ytitle_s, 'RHO')
    repeat begin
      if pos ge 0L then ytitle_s = strmid(ytitle_s, 0L, pos) + $
                                   '!9r!x' + strmid(ytitle_s, pos + 3L)
    endrep unitl pos eq - 1L

  endelse

  ytitle = temporary(ytitle_s)


  ;;===========================================================================
  ;;===========================================================================
  ;; Loop over the spectra:

  n_spec_na = n_elements(spec_na)
  icolors = intarr(n_spec_na)
  ilinestyle = intarr(n_spec_na)
  iitems = strarr(n_spec_na)
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
              /encapsulated, /portrait, xsize=psxsize, ysize=psysize
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
      sdens = '!9r!x!dg!n'
      somega = '!9W!x!drms!n'

    endif else begin

      sdens = '!4q!x!dg!n'
      somega = '!4X!x!drms!n'

      charthick_use = charthick
      charsize_use = charsize

      staus = '!4s!x!ds!n='
      staut = '!4s!x!dt!n='
      stKol = '!4s!x!dK!n='
      sStokes = ', St='

      wset, win_id_pixmap
    endelse

    title = 't=' + strtrim(string(t[i], format='(f10.3)'), 2L) + $
            ' :: ' + staut + str_tau_turnover
    if epsk_available then title += ' :: ' + stKol + str_tau_Kolmogorov

    xstyle = 1
    ystyle = plot_na && separate_yrange ? 9 : 1
    xmargin = plot_na && separate_yrange ? [10.0, 10.0] : [10.0, 0.3]
    ymargin = [3.3, 0.3]
    plot, k, /nodata, background=colorbg, color=colorfg, $
          charsize=charsize, charthick=charthick, $
          /xlog, xmargin=xmargin, xrange=xrange, xstyle=xstyle, $
          xthick=thick, xtitle=xtitle, $
          /ylog, ymargin=ymargin, yrange=yrange, ystyle=ystyle, $
          ythick=thick, ytitle=ytitle

    xpos = 0.50d0 * (!x.window[1L] - !x.window[0L]) + !x.window[0L]
    ypos = 0.95d0 * (!y.window[1L] - !y.window[0L]) + !y.window[0L]
    xyouts, xpos, ypos, title, charsize=charsize, charthick=charthick, $
            color=colorfg, alignment=0.5, /normal

    if plot_kin || plot_lnrho then begin
      ;; guide: 1 / k^2
      k_l1 = 2.5d0 & k_l2 = 9d0
      k_arr = [k_l1, k_l2]
      y_arr = 1d0 / k_arr ^ 2 * 1d-2
      oplot, k_arr, y_arr, linestyle=0, color=coloror, thick=thick * 2
  
      ;; guide: 1 / k^(5/3)
      k_l1 = 2.5d0 & k_l2 = 9d0
      k_arr = [k_l1, k_l2]
      y_arr = 1d0 / k_arr ^ (5d0 / 3d0) * 1d-2
      oplot, k_arr, y_arr, linestyle=0, color=colorgn, thick=thick * 2
    endif

    if plot_kin then begin
      oplot, k, spec_kin[*, i] / k ^ 2, color=colorfg, thick=thick, linestyle=5
    endif

    if plot_omega then begin
      oplot, k, spec_o[*, i] / k ^ 2, color=colorpr, thick=thick
    endif

    if plot_lnrho then begin
      oplot, k, spec_lr[*, i], color=colorlb, thick=thick
    endif

    if plot_kin then begin
      l_items = ['E!dk!n']
      l_colors = [colorfg]
      l_linestyle = [5]
      l_thick = thick * [1.0]
    endif

    if plot_lnrho then begin
      l_items = ~n_elements(l_items) ? [sdens] : [l_items, sdens]
      l_colors = ~n_elements(l_colors) ? [colorlb] : [l_colors, colorlb]
      l_linestyle = ~n_elements(l_linestyle) ? [0] : [l_linestyle, 0]
      l_thick = ~n_elements(l_thick) ? [1.0] : [l_thick, 1.0]
    endif
    
    if plot_omega then begin
      l_items = ~n_elements(l_items) ? [somega] : [l_items, somega]
      l_colors = ~n_elements(l_colors) ? [colorpr] : [l_colors, colorpr]
      l_linestyle = ~n_elements(l_linestyle) ? [0] : [l_linestyle, 0]
      l_thick = ~n_elements(l_thick) ? [1.0] : [l_thick, 1.0]
    endif

    if plot_kin || plot_lnrho then begin
      l_items = [l_items, 'guide: k!u-2!n', 'guide: k!u5/3!n']
      l_colors = [l_colors, coloror, colorgn]
      l_linestyle = [l_linestyle, 0, 0]
      l_thick = [l_thick, 2.0, 2.0]
    endif

    if plot_na then begin
      l_items = ~n_elements(l_items) ? 'guide: k!u2!n' : $
                [l_items, 'guide: k!u2!n']
      l_colors = [l_colors, colorrd]
      l_linestyle = [l_linestyle, 0]
      l_thick = [l_thick, 2.0]
    endif

    al_legend, l_items, linestyle=l_linestyle, /top, /right, box=0, $
               charsize=charsize, charthick=charthick, color=l_colors, $
               thick=l_thick * thick, linsize=0.5, $
               textcolors=colorfg, spacing=spacing


    if plot_na then begin

      if separate_yrange then begin

        if separate_yrange_shift then begin
          yrange_na = yrange + shift_yrange[0L]
        endif else begin
          for ii = 0UL, n_spec_na - 1UL do begin
            idx = where(k gt xrange[0L] and k le xrange[1L], count)
            min = min((*spec_na[ii])[idx, i], max=max)
            yrange_na = ~ii ? [min, max] : $
                        [min([yrange_na[0L], min]), max([yrange_na[1L], max])]
          endfor
        endelse

        axis, charsize=charsize, charthick=charthick, color=colorfg, /yaxis, $
              /ystyle, yrange=yrange_na, ythick=thick, ytitle=ytitle_na, /save

      endif


      ;; guide: k^2
      k_l1 = 1d1 & k_l2 = 4d1
      k_arr = [k_l1, k_l2]
      y_arr = k_arr ^ 2 / k_l1 ^ 2 * 0.5 * 1d1 ^ !y.crange[0L]
      oplot, k_arr, y_arr, linestyle=0, color=colorrd, thick=thick * 2

      linestyle = 2

      for ii = 0UL, n_spec_na - 1UL do begin
        color = na_colors[ii / 5L]
        oplot, k, scalefactor_use * (*spec_na[ii])[*, i], $
               linestyle=linestyle, color=color, thick=thick

        icolors[ii] = color
        iitems[ii] = 'a=' + aps[ii] + ', ' + $
                     staus + str_tau_stop[ii] + sStokes + vStokes[ii]
        ilinestyle[ii] = linestyle

        linestyle ++
        if linestyle eq 6 then linestyle = 1
      endfor

      al_legend, iitems, linestyle=ilinestyle, /bottom, box=0, thick=thick, $
                 charsize=charsize, charthick=charthick, color=icolors, $
                 linsize=0.5, textcolors=colorfg, spacing=spacing

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
