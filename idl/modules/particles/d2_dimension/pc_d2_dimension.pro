;; This piece of code calculates nearest-neighbor statistics, and the
;; D2 value. This is accomplished by calling the cKDTree Python /
;; Scipy module and feeding it with the particle position data that
;; are available in the PVAR files...and then create a histogram of
;; the distances.
;;
;; The calculation of the nearest neighbor statistics is very time
;; consuming (in the way implemented here), which is why the
;; resulting histograms are saved to (compressed) text files.
;;
;; The final part of this tool creates a plot of the histogram
;; and fits a line to the slope of the lowermost side. The slope
;; is the D2 value.
;;
;;
;; The tool handles both regular pencil-code binary files as well
;; as HDF5 files. Using regular binary particle files, set the
;; paramater file argument PVARS as, for example:
;;
;;  IDL> pvars = ['PVAR' + ['20', '24', '28'], 'pvar.dat']
;;
;; where as with HDF5 files, you could use:
;;
;;  IDL> pvars = 'allprocs/' + ['PVAR' + ['20', '24', '28'], 'pvar'] + '.h5'
;;
;;
;; To use this code with the full Python routine, use:
;;  pc_d2_dimension, pvars, /allpython
;;
;; To use this code with the IDL, plotting and saving figures,
;; ensure that mpfit.pro and mpfitfun.pro are in the IDL path.
;; Call the tool using:
;;  pc_d2_dimension, pvars
;;
;;
;; Output:
;;   If called will ALLPYTHON=0, all histogram data are saved to
;;   compressed text files, which are read automatically at the
;;   next call to this routine, if they are present in the current
;;   directory.
;;
;;   The D2 statistics are saved to a plain-text file named
;;   'nnd_d2.dat', which contents can then be plotted using, for
;;   example, pc_d2_dimension_plot.
;;
;;   Image files are created showing the histogram and the fit of
;;   the slope; the format is PNG. The plots are LIN-LIN, unless
;;   the keyword LOG is set, and in this case they are LOG-LOG plots.
;;
;;
;; Requirements:
;;   IDL version 8.5 (for the Python bridge), or higher; especially
;;   if you want to use Python version >= 3.5, then you need IDL version
;;   8.6 or newer. Also, IDL needs to be launched from an environment
;;   where "python" invokes a version of Python that is supported by the
;;   used version of IDL.
;;
;;   In addition to IDL, it is necessary to dowmload mpfit.pro and
;;   mpfitfun.pro from the Markwardt library (and put them in the IDL
;;   path):
;;     http://cow.physics.wisc.edu/~craigm/idl/fitting.html
;;
;;   Also, the IDL astronomy library is required:
;;     https://idlastro.gsfc.nasa.gov/
;;
;;   It is perhaps the easiest to install a recent version of
;;   Anaconda and then inside of Anaconda install a recent version
;;   of Python that is supported by the  used version of IDL. Before
;;   starting IDL, one needs to activate the correct Anaconda Python
;;   environment. For example,
;;
;;   > conda deactivate # To deactivate any currently used Anaconda
;;                      # Python version that is incompatible with
;;                      # IDL.
;;   > conda activate py36 # Activate Python version 3.6, which would
;;                      # work with, e.g., IDL 8.7.3. Before trying
;;                      # this, it is necessary to install Python 3.6
;;                      # inside of Anaconda (see the documentation;
;;                      # this is easy and straightforward). Also,
;;                      # install the two packages "scipy" and "numpy"
;;                      # inside of this package.
;;   > idl
;;   IDL> pvars = 'allprocs/' + ['PVAR' + ['20', '24', '28'], 'pvar'] + '.h5'
;;   IDL> pc_d2_dimension, pvars
;;
;;   Note!
;;     On some machines, such as certain clusters (SNIC/PDC/Beskow and Tegner),
;;     it is necessary to adjust environment variables in the correct order
;;     to use the Python bridge from within IDL. Here is an approach that I
;;     found is working for me (CS):
;;
;;     1) Assuming Anaconda is installed and that it uses a version of
;;       Python that isn't supported by IDL as yet (currently, 2020-06-03,
;;       IDL can only use Python 3.6, while Anaconda comes with Python 3.7),
;;       install an environment using Python 3.6:
;;
;;       $ conda create -n py36 python=3.6
;;       $ conda deactivate
;;       $ conda activate py36
;;       $ conda install scipy
;;
;;
;;     2) Create a script that sets up the environment variables that need
;;       to be correctly set. Create a file that looks something like the
;;       following; I call this script "py36_runme.sh":
;;
;;       === Contents of py36_runme.sh:
;;       #!/bin/bash
;;       # Adjust the following line no point at the location of your own Anaconda installation:
;;       export PYTHONHOME=/cfs/klemming/nobackup/c/csandin/anaconda3/envs/py36
;;       export PYTHONPATH=$PYTHONPATH:${IDL_DIR}/lib/bridges:${IDL_DIR}/bin/bin.linux.x86_64
;;       export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$PYTHONHOME/lib:/usr/lib64:${IDL_DIR}/bin/bin.linux.x86_64
;;       ===
;;
;;
;;     3) Whenever you log in to the machine in question and wish to use
;;       the Python bridge from IDL, change environment and adjust the
;;       environment variables. Finally, launch IDL:
;;
;;       $ conda deactivate
;;       $ conda activate py36
;;       $ source ./py36_runme.sh
;;       $ idl
;;
;;
;; Written by: Christer Sandin, 2020
;;
pro pc_d2_dimension_write, str, append=append, file=file
  compile_opt hidden, IDL2

  on_ioerror, err
  openw, unit, file, /get_lun, append=append
  printf, unit, str
  free_lun, unit
  on_ioerror, NULL

  if ~append then append = 1L


  return

err:
  if n_elements(unit) eq 1L then begin
    fstat = fstat(unit)
    if fstat.open then begin
      free_lun, unit
      on_ioerror, NULL
    endif
  endif

  return
end


pro pc_d2_dimension_function, x, a, f, pder
  compile_opt hidden, IDL2

  z = - a[0L] * x ^ a[1L]
  w = x ^ (a[1L] - 1d0) * exp(z)
  q = a[1L] * w

  f = a[0L] * q

  if n_params() ge 4 then $
     pder = [[q * (1d0 + z)], [a[0L] * (w + q * alog(x) * (1d0 + z))]]
  
  return
end


function pc_d2_dimension_mpfit_f, x, p
  compile_opt hidden, IDL2

  return, p[0L] * p[1L] * x ^ (p[1L] - 1d0) * exp(- p[0L] * x ^ p[1L])
end


pro pc_d2_dimension_hplot, locations, obs, i, j, fitx=fitx, fity=fity, $
        filename=filename, log=log, charsize=charsize, font_size=font_size, $
        noshow=noshow, psfilename=psfilename, psxsize=psxsize, $
        psysize=psysize, xsize=xsize, ysize=ysize, legendtext=legendtext, $
        comparison_fit=comparison_fit, thick=thick, yright=yright, xrange=xrange, $
        animate=animate, data=data, pstruc=pstruc, j_max=j_max, $
        nplot_x=nplot_x, nplot_y=nplot_y, do_postscript=do_postscript, $
        ap=ap, file_string=file_string, warr=warr, do_grab=do_grab
  compile_opt hidden, IDL2

  do_postscript = keyword_set(do_postscript)
  if ~n_elements(charsize) then charsize = 1.0
  if ~n_elements(thick) then thick = 1.0
  charthick = thick
  if ~n_elements(log) then log = 1L
  if ~n_elements(noshow) then noshow = 1L
  if ~n_elements(psxsize) then psxsize = 6.5
  if ~n_elements(psysize) then psysize = 4.0
  if ~n_elements(font_size) then font_size = 10.0
  do_grab = keyword_set(do_grab)

  if ~j then begin
    pf = !p.font
    device_name = !d.name

    if do_postscript then begin

      set_plot, 'ps'
      device, xsize=psxsize, ysize=psysize
      device, filename=psfilename, /isolatin, /color, bits_per_pixel=8, $
              /encapsulated, inches=0, font_size=font_size
      !p.font = 0L

    endif else begin

      if noshow || animate then begin
        window, /free, xsize=xsize, ysize=ysize, /pixmap
        if animate && ~do_grab then warr[i] = !d.window
      endif else begin
        window, 0L, xsize=xsize, ysize=ysize
      endelse

    endelse

    xo = !x.omargin
    yo = !y.omargin
    !x.omargin = keyword_set(yright) ? [0.5, 3.0] : [3.0, 0.5]
    !y.omargin = [3.0, 0.0]

    device, get_decomposed=decomposed
    device, decomposed=0L
    bottom = 7L & ncolors = !d.table_size - bottom
    loadct, 74, bottom=bottom, ncolors=ncolors, /silent
    tvlct,   0b,   0b ,  0b, 0L         ;; Black
    tvlct, 255b, 255b, 255b, 1L         ;; White
    tvlct, 213b,  94b,   0b, 2L         ;; Vermillion
    tvlct,   0b, 158b, 115b, 3L         ;; Bluish green
    tvlct,   0b, 114b, 178b, 4L         ;; Blue
    tvlct, 230b, 159b,   0b, 5L         ;; Orange
    tvlct, 204b, 121b, 167b, 6L         ;; Purple (Mulberry)

    pm = !p.multi
    !p.multi = [0L, nplot_x, nplot_y, 0L, 0L]

    pstruc = {pm:pm, pf:pf, xo:xo, yo:yo, decomposed:decomposed, $
              device_name:device_name, window:!d.window}

  endif else begin
    !p.multi = [nplot_x * nplot_y - j, nplot_x, nplot_y, 0L, 0L]
  endelse

  background = 1L
  color = 0L

  jx_pos = j mod nplot_x
  jy_pos = j / nplot_x
  jx_pos_max = (j_max - 1L) mod nplot_x
  show_x = (jy_pos eq (nplot_y - 1L)) or $
           ((jy_pos eq (nplot_y - 2L)) and jx_pos gt jx_pos_max)

  if keyword_set(yright) then begin
    show_y = (jx_pos eq nplot_x - 1L) or $
             ((jy_pos eq (nplot_y - 2L)) and jx_pos gt jx_pos_max)
  endif else begin
    show_y = ~jx_pos
  endelse

  charsize_leg = charsize * (nplot_x gt 2 or nplot_y gt 2 ? 0.5 : 1.0)
  xmargin = keyword_set(yright) ? [0.5, 5.0] : [5.0, 0.5]
  ymargin = [0.5, 0.5]
  xtitle_ = 'Minimum distance bins'
  if keyword_set(log) then xtitle_ = 'log ' + xtitle_
  xtitle = show_x ? xtitle_ : ''

  locations_ = keyword_set(log) ? alog10(locations) : locations

  yrange = keyword_set(log) ? alog10([min(obs) > 1d-3, max(obs)]) : [min(obs), max(obs)]
  if keyword_set(log) then begin
    obs_ = obs * 0
    idx = where(obs gt 1d-10, count, complement=cidx, ncomplement=ncount)
    if  count ne 0L then obs_[ idx] = alog10(obs[idx])
    if ncount ne 0L then obs_[cidx] = -10
  endif else obs_ = obs
  obs_[0L] = - 10

  ystyle = 1L + (keyword_set(yright) ? 8L : 0L)

  ytitle_ = 'Normalized counts'
  if keyword_set(log) then ytitle_ = 'log ' + ytitle_
  ytitle = show_y && ~keyword_set(yright) ? ytitle_ : ''
  ytickformat = keyword_set(yright) ? '(a1)' : ''

  xtickformat = show_x ? '' : '(a1)'

  plot, locations_, obs_, charthick=charthick, /nodata, $
        background=background, charsize=charsize, color=color, $
        xmargin=xmargin, xrange=xrange, /xstyle, xthick=thick, $
        xtickformat=xtickformat, xtitle=xtitle, $
        ymargin=ymargin, yrange=yrange, ystyle=ystyle, ythick=thick, $
        ytitle=ytitle, ytickformat=ytickformat

  if keyword_set(yright) then begin
     ytitle = show_y ? ytitle_ : ''
     axis, color=color, /yaxis, yrange=yrange, charsize=charsize, $
           charthick=charthick, /ystyle, ythick=thick, ytitle=ytitle
  endif

  oplot, locations_, obs_, color=4, thick=thick

  if n_elements(fitx) gt 0L then begin
    fitx_ = keyword_set(log) ? alog10(fitx) : fitx
    fity_ = keyword_set(log) ? alog10(fity) : fity
    oplot, fitx_, fity_, linestyle=3, color=2, thick=thick * 1.5

    fity_ = keyword_set(log) ? alog10(comparison_fit) : comparison_fit
    oplot, fitx_, fity_, linestyle=1, color=6, thick=thick

    items = ['d!d2!n fit', 'd!d2!n=3.0']

    al_legend, items, charsize=charsize_leg, charthick=charthick, box=0, $
               linestyle=[3, 1], colors=[2, 6], $
               thick=thick, linsize=0.33, textcolors=color

  endif

  if animate then begin
    xpos = 0.05 * (!x.crange[1L] - !x.crange[0L]) + !x.crange[0L]
    ypos = 0.08 * (!y.crange[1L] - !y.crange[0L]) + !y.crange[0L]
    a_string = 'a=' + strtrim(string(ap, format='(e8.1)'), 2L)
    xyouts, xpos, ypos, a_string, alignment=0.0, color=color, $
            charsize=charsize_leg, charthick=charthick
    if ~j then begin
      ypos += 1.5 * !d.y_ch_size * charsize / ((!y.window[1L] - !y.window[0L]) * !d.y_size)
      xyouts, xpos, ypos, file_string, alignment=0.0, color=color, $
              charsize=charsize_leg, charthick=charthick
    endif
  endif

  if n_elements(legendtext) eq 1L then begin
    al_legend, legendtext, /bottom, /left, charsize=charsize, box=0, textcolor=color
  endif

  if j eq j_max - 1L then begin
    if do_postscript then begin

      device, /close_file
      set_plot, pstruc.device_name
      !p.font = pstruc.pf

    endif else begin

      ;; Read image buffer and save to PNG file:
      device, decomposed=1L
      grab = tvrd(true=1)
      if animate then begin
        if do_grab then data[i] = ptr_new(grab)
      endif else begin
        tvlct, red, green, blue, /get
        write_png, filename, grab, red, green, blue
      endelse
      device, decomposed=pstruc.decomposed

      if noshow || (animate && do_grab) then wdelete, pstruc.window

    endelse
  endif


  return
end ;;; pc_d2_dimension_hplot


pro pc_d2_dimension_animate, data, warr=warr, xsize=xsize, ysize=ysize, sleep=sleep
  compile_opt hidden, IDL2

  size = size(data, /dimensions)
  if ~n_elements(sleep) then sleep = 0.1
  do_grab = min(warr) eq - 1L ? 1L : 0L

  if ~n_elements(xsize) then xsize = !d.x_size
  if ~n_elements(ysize) then ysize = !d.y_size

  window, /free, xsize=xsize, ysize=ysize
  wid = !d.window

  print, 'Loop: Interrupt with Ctrl-C.'
  print, '   Enter "loop=0" followed by ".cont" to finish animation.'
  print, '   Enter "pingpong=1" followed by ".cont" to have animation go back and forward.'
  print, '   Enter "sleep=x.x" followed by ".cont" to set time between draws [s] (x >= 0.01).'

  loop = 1L
  i_plus = 1L
  pingpong = 0L

  device, get_decomposed=decomposed
  device, decomposed=1L
  while loop do begin

    if i_plus then begin
      i_0 = 0L
      i_1 = size[0L] - 1L
      i_stride = 1L
    endif else begin
      i_0 = size[0L] - 2L
      i_1 = 1L
      i_stride = - 1L
    endelse

    for i = i_0, i_1, i_stride do begin
      if do_grab then begin
        tv, *data[i], /true
      endif else begin
        copy = [0L, 0L, xsize, ysize, 0L, 0L, warr[i]]
        device, copy=copy
      endelse
      wait, sleep

    endfor

    if pingpong then i_plus = ~i_plus else i_plus = 1L

  endwhile
  device, decomposed=decomposed

  return
end ;;; pc_d2_dimension_animate


pro pc_d2_dimension, pvars, allpython=allpython, recalculate=recalculate, $
                     cumulative=cumulative, fixed_bins=fixed_bins, charsize=charsize, $
                     thick=thick, font_size=font_size, psxsize=psxsize, psysize=psysize, $
                     xsize=xsize, ysize=ysize, noshow=noshow, yright=yright, log=log, $
                     legendtext=legendtext, xrange=xrange, sleep=sleep, $
                     postscript=postscript, animate=animate, grab=do_grab
  compile_opt IDL2

  python_all = keyword_set(allpython)
  animate = keyword_set(animate) && ~python_all
  cumulative = keyword_set(cumulative)
  cumul_str = ''
  if ~n_elements(fixed_bins) then fixed_bins = 1L
  do_postscript = ~n_elements(postscript) ? 1L : 0L

  ;; Determine if the input file is regular binary or HDF5:

  hdf5 = 0L

  qpvars_1 = file_search('data/proc0/PVAR?')
  qpvars_1 = file_basename(qpvars_1)
  qpvars_2 = file_search('data/proc0/PVAR??')
  qpvars_2 = file_basename(qpvars_2)
  qpvars_3 = file_search('data/proc0/PVAR???')
  qpvars_3 = file_basename(qpvars_3)
  if qpvars_1[0L] eq '' && $
     qpvars_2[0L] eq '' && $
     qpvars_3[0L] eq '' then hdf5 = 1L


  pc_read_param, object=ps

  ap0 = ps.ap0 & nap0 = n_elements(ap0)
  j_str = ' / ' + strtrim(nap0, 2L)

  nplot_x = ceil(sqrt(nap0))
  nplot_y = nplot_x - 1L

  ;; Compile (Import) Python program (change directory to find the routine):
  if python_all then begin
    cd, current=pwd
    cd, file_dirname(routine_filepath())
    python_aq = Python.Import('pc_d2_dimension__import_CKD')
    cd, pwd
  endif


  ;; Cumulative sum, create a pointer array:
  if cumulative then obs_sum = ptrarr(nap0)


  ;; Using ANIMATE, one doesn't set PVARS, but the tool finds out which ones
  ;; exist by itself:

  if animate and ~n_elements(pvars) then begin

    pstr = hdf5 ? '.' : '_'

    for ij = 0UL, 2UL do begin
      case ij of
        0UL: files = file_search('histogram_PVAR?' + pstr + '*_hist.dat.gz')
        1UL: files = file_search('histogram_PVAR??' + pstr + '*_hist.dat.gz')
        2UL: files = file_search('histogram_PVAR???' + pstr + '*_hist.dat.gz')
      endcase
      if files[0L] eq '' then continue

      nfiles = n_elements(files)
      for i = 0UL, nfiles - 1UL do begin
        if hdf5 then begin
          pos = strpos(files[i], '.')
          files[i] = strmid(files[i], 10L, pos - 10L)
        endif else begin
          pos = strpos(strmid(files[i], 1L + strpos(files[i], '_')), '_')
          files[i] = strmid(files[i], 10L, pos)
        endelse
      endfor
      if ~n_elements(pvars) then begin
        pvars = files[uniq(files, sort(files))]
      endif else begin
        pvars = [pvars, files[uniq(files, sort(files))]]
      endelse

    endfor

    if hdf5 then pvars = 'allprocs' + path_sep() + pvars + '.h5'

    ;; Create pointer array that stores the data:
    data = ptrarr(n_elements(pvars))
    warr = lonarr(n_elements(pvars)) - 1L
  endif

  ;; Open the output file:
  file = 'nnd_d2.dat'
  append = 0L
  stb = string(9b)

  ;; Loop over the parameter files:

  npvars = n_elements(pvars)
  for i = 0UL, npvars - 1UL do begin

    print, "PVAR = ", pvars[i]

    ;; Histogram filenames are setup here:

    filenames = 'histogram_' + file_basename(pvars[i]) + '_' + $
                strtrim(lindgen(nap0), 2L) + '_hist.dat.gz'
    load_data = 0UL
    for j = 0UL, nap0 - 1UL do $
       if ~file_test(filenames[j], /read, /regular) then load_data = 1UL
    if python_all then load_data = 1UL

    if load_data then begin
      if hdf5 then begin
        xp = pc_read('part/xp', file=pvars[i])
        yp = pc_read('part/yp', file=pvars[i])
        zp = pc_read('part/zp', file=pvars[i])
        ap = pc_read('part/ap', file=pvars[i])
        nap = n_elements(ap)

        lap = alog10(ap)
      endif else begin
        pc_read_pvar, varfile=pvars[i], rmv=rmv, array=array
        s_array = size(array, /dimensions)
        nap = s_array[0L]

        lap = alog10(array[*, 6L])
      endelse

      ;; Calculate bin indices for all particles - fast:
      diff = lap[1L : nap0 - 1L] - lap[0L : nap0 - 2L]
      fast_indices = stdev(diff) eq 0d0 ? 1L : 0L
      if fast_indices then $
         h = histogram(lap, min=min(lap), nbins=nap0, reverse_indices=ri)
      print, 'Fast indices: ' + (fast_indices ? 'yes' : 'no')

      ;; Statistics values:
      N_gr = nap / nap0
      xdes = 0.55396 * ps.lx0 / N_gr ^ (1d0 / 3d0)

    endif ;; load_data


    ;; Loop over the grain sizes:

    for j = 0UL, nap0 - 1UL do begin
      print, "    Particle size " + strtrim(j + 1L, 2L) + j_str + $
             ' [' + string(ap0[j], format='(e10.3)') + ']'

      if load_data then begin

        if fast_indices then begin

          ;; Retrieving the size-specific indices in a fast way:
          idx = ri[ri[j] : ri[j + 1L] - 1L]
          count = ri[j + 1L] - ri[j]

        endif else begin

          if hdf5 then begin
            idx = where(ap gt 0.999d0 * ap0[j] and $
                        ap lt 1.001d0 * ap0[j], count)
          endif else begin
            idx = where(array[*, 6L] gt 0.999d0 * ap0[j] and $
                        array[*, 6L] lt 1.001d0 * ap0[j], count)
          endelse

        endelse

        nbins = long(ceil(sqrt(count)))

        if ~count then begin
          message, 'No particles of this size, skip to the next size.', $
                   /informational
          continue
        endif

      endif ;; load_data

      if python_all then begin

        ;; Approach A: Use the Python module:

        str = python_aq.pc_d2_dimension_import_CKD(i, hdf5 ? xp[idx] : array[idx, 0L], $
                                          hdf5 ? yp[idx] : array[idx, 1L], $
                                          hdf5 ? zp[idx] : array[idx, 2L], $
                                          nap0, xdes, ap0[j])
        print, str

      endif else begin

        ;;========================================
        ;; Retrieve or calculate histogram:

        if ~file_test(filenames[j], /read, /regular) || $
           keyword_set(recalculate) then begin

          ;; Approach B: Call Python commands, but stay in IDL:

          Python.Run, 'import numpy as np'
          Python.Run, 'from scipy.spatial import cKDTree'
          Python.Run, 'from scipy.optimize import curve_fit'

          ;; Transfer to Python as Numpy arrays:
          Python.pxp = Python.Wrap(hdf5 ? xp[idx] : array[idx, 0L])
          Python.pyp = Python.Wrap(hdf5 ? yp[idx] : array[idx, 1L])
          Python.pzp = Python.Wrap(hdf5 ? zp[idx] : array[idx, 2L])
          Python.Run, 'p = np.stack((pxp,pyp,pzp)).T'
          print, ' Tree: create'
          Python.Run, 'tree = cKDTree(p)'
          print, ' Tree: created'

          Python.xdes = xdes

          Python.Run, 'xmil = np.zeros(' + strtrim(count, 2L) + ')'
          Python.Run, 'for rps in range(0, ' + strtrim(count, 2L) + '):\n' + $
                      '  rps3d = np.stack((pxp[rps], pyp[rps], pzp[rps]))\n' + $
                      '  xmin_kdt, index = tree.query(rps3d, k=2)\n' + $
                      '  xmil[rps] = xmin_kdt[1] / xdes'

          xmil = Python.xmil


          ;;========================================

          if fixed_bins then begin

            binsize = 3.3333d-3
            min = 0d0
            nbins = 1100L & nbins_fit = nbins
            obs = histogram(xmil, binsize=binsize, $
                            min=min, nbins=nbins, locations=locations)

          endif else begin

            ;;========================================
            ;; Calculate the histogram: round 1:

            obs = histogram(xmil, nbins=nbins, locations=locations)

            norm = dblarr(n_elements(obs))
            for iii = 0L, n_elements(obs) - 2L do $
               norm[iii] = (locations[iii + 1L] - locations[iii]) * obs[iii]
            obs /= total(norm) ;; Normalize data

            for jj = 0UL, nbins - 2UL do locations[jj] += locations[jj + 1UL]
            locations /= 2d0

            ;; Locate maximum, and step back to a third:
            obs_max = max(obs[1L : *], obs_max_idx)
            obs_max_idx_13 = long(mean(obs_max_idx) / 3d0)

            ;;filename = 'nnd_d2_' + file_basename(pvars[i]) + '_' + $
            ;;           strtrim(j, 2L) + '_h1.png'
            ;;pc_d2_dimension_plot, locations, obs, filename=filename


            ;; Calculate new number of bins:
            nbins_fit = long(1d2 / (obs_max_idx_13 + 1L) * nbins)


            ;;========================================
            ;; Histogram: round 2:

            obs = histogram(xmil, nbins=nbins_fit, locations=locations)

          endelse

          norm = dblarr(n_elements(obs))
          for iii = 0L, n_elements(obs) - 2L do $
             norm[iii] = (locations[iii + 1L] - locations[iii]) * obs[iii]
          obs /= total(norm) ;; Normalize data

          ;; Save the histogram data for later use:
          openw, unit, filenames[j], /get_lun, /compress
          printf, unit, ';; Format: location, normalized histogram value'
          printf, unit, '; normalization factor: ' + $
                  string(format='(e15.8)', total(norm))
          format = '(e15.8, 2x, e15.8)'
          for jj = 0UL, n_elements(obs) - 1UL do $
             printf, unit, format=format, locations[jj], obs[jj]
          free_lun, unit

        endif else begin

          ;;========================================
          ;; Retrieve the histogram from the saved file:

          readcol, filenames[j], /compress, format='(f,f)', locations, obs, /silent
          nbins_fit = n_elements(obs)

          if keyword_set(cumulative) then begin

            if ~ptr_valid(obs_sum[j]) then begin
              obs_sum[j] = ptr_new(temporary(obs))
            endif else begin
              *obs_sum[j] += obs
            endelse

            ;; Normalize the total sum once all data were loaded:
            if i eq npvars - 1UL then begin

              norm = dblarr(n_elements(*obs_sum[j]))
              for iii = 0L, n_elements(*obs_sum[j]) - 2L do $
                 norm[iii] = (locations[iii + 1L] - locations[iii]) * (*obs_sum[j])[iii]
              obs = *obs_sum[j] / total(norm)

              ptr_free, obs_sum[j]

              cumul_str = '_cumul'

            endif else begin

              ;; Skip to the next file:
              continue

            endelse

          endif ;; cumulative

        endelse

        for jj = 0UL, nbins_fit - 2UL do locations[jj] += locations[jj + 1UL]
        locations /= 2d0

        ;; Locate maximum in histogram, and step back to a third:
        obs_max = max(obs[1L : *], obs_max_idx)
        obs_max_idx_13 = long(mean(obs_max_idx) / 3d0)


        ;;========================================
        ;; Fit a curve to the lower end of values:

        x =  locations[1L : obs_max_idx_13]
        y = double(obs[1L : obs_max_idx_13])

        ftol = 1d-6
        maxiter = 2000L
        dy = dblarr(n_elements(x)) + 1d0
        pi_s = {fixed:0L, limited:[0L, 0L], limits:[0d0, 0d0], value:0d0}
        pi = replicate(temporary(pi_s), 2L)
        pi[1L].limited = [1L, 1L]
        pi[1L].limits  = [0.5d0, 3.5d0]
        pi.value = [1d0, 1d0] ;; starting values
        p = mpfitfun('pc_d2_dimension_mpfit_f', x, y, dy, $
                     ftol=ftol, maxiter=maxiter, niter=niter, $
                     parinfo=pi, status=status, /quiet)
        d2 = p[1L]
        yfit = pc_d2_dimension_mpfit_f(x, p)

        yfit_3 = pc_d2_dimension_mpfit_f(x, [p[0L], 3d0])

        ;; Another way to do it using the IDL-to-Python bridge (*very* slow):
        ;;
        ;;Python.xp_fit  = Python.Wrap(x)
        ;;Python.obs_fit = Python.Wrap(y)
        ;;Python.Run, 'f_D = lambda xp_fit, a, D: a*D*xp_fit**(D - 1)*' + $
        ;;   'np.exp(-1.0*a*xp_fit**D)'
        ;;Python.Run, 'popt, pcov = curve_fit(f_D, xp_fit, obs_fit)'
        ;;Python.Run, 'yfit = popt[0] * popt[1] * xp_fit**(popt[1]-1) * ' + $
        ;;   'np.exp(-1.0 * popt[0] * xp_fit**popt[1])'
        ;;popt = Python.popt
        ;;yfit = Python.yfit
        ;;D2 = popt[1L]

        filename = 'nnd_d2_' + file_basename(pvars[i]) + cumul_str + '_hist.png'
        psfilename = 'nnd_d2_' + file_basename(pvars[i]) + cumul_str + '_hist.eps'
        file_string = file_basename(pvars[i])
        pos = strpos(file_string, '.')
        if pos gt 0L then file_string = strmid(file_string, 0L, pos)
        pc_d2_dimension_hplot, locations, obs, i, j, fitx=x, fity=yfit, log=log, $
            filename=filename, psfilename=psfilename, psxsize=psxsize, $
            psysize=psysize, xsize=xsize, ysize=ysize, noshow=noshow, $
            charsize=charsize, comparison_fit=yfit_3, legendtext=legendtext, $
            thick=thick, font_size=font_size, yright=yright, xrange=xrange, $
            animate=animate, data=data, pstruc=pstruc, $
            j_max=nap0, nplot_x=nplot_x, nplot_y=nplot_y, do_postscript=do_postscript, $
            ap=ap0[j], file_string=file_string, warr=warr, do_grab=do_grab


        if load_data then begin
          eformat = '(e10.3)'
          fformat = '(f6.3)'
          str = strtrim(i, 2L) + $
                stb + string(ap0[j], format=format) + $
                stb + strtrim(string(d2, format=fformat), 2L) + $
                stb + string(mean(xmil), format=eformat) + $
                stb + string(stdev(xmil), format=eformat)
        endif

      endelse

      if load_data then pc_d2_dimension_write, str, append=append, file=file
      
    endfor ;; j = 0UL, nap0 - 1UL

    if load_data then pc_d2_dimension_write, '', append=append, file=file

  endfor

  if animate then $
     pc_d2_dimension_animate, data, xsize=xsize, ysize=ysize, sleep=sleep, warr=warr

  ;; Release allocated memory, if any:
  if animate then begin

    if keyword_set(do_grab) then begin
       for i = 0UL, n_elements(data) - 1UL do $
          if ptr_valid(data[i]) then ptr_free, data[i]
    endif else begin
       for i = 0UL, n_elements(pvars) - 1UL do wdelete, warr[i]
    endelse

  endif


  return
end ;;; procedure: pc_d2_dimension
