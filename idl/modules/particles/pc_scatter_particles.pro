pro pc_scatter_particles, pvar, ax=ax, az=az, fraction=fraction_, $
       center=center_, Lxw=Lxw_, Lyw=Lyw_, Lzw=Lzw_, $
       charsize=charsize, charthick=charthick, colortable=colortable, $
       psym=psym, symsize=symsize, $
       xsize=xsize, ysize=ysize, thick=thick, font_size=font_size, $
       legendtext=legendtext, xrange=xrange, $
       ofilename=ofilename
  compile_opt IDL2

  if ~n_elements(Lxw_) then Lxw_ = 0L
  if ~n_elements(Lyw_) then Lyw_ = 0L
  if ~n_elements(Lzw_) then Lzw_ = 0L
  if ~n_elements(colortable) then colortable = 74
  if ~n_elements(psym) then psym = 3
  if ~n_elements(symsize) then symsize = 1.0
  if ~n_elements(xsize) then xsize = 800L
  if ~n_elements(ysize) then ysize = xsize
  if ~n_elements(fraction_) then fraction_ = 1d0
  fraction = fraction_ < 1.0 > 0.0
  if ~n_elements(ofilename) then ofilename = 'pc_scatter_particles.png'

  ;; Determine if the input file is regular binary or HDF5:

  pos = strpos(pvar, '.h5', /reverse_search)
  hdf5 = pos eq strlen(pvar) - 3L ? 1L : 0L

  pc_read_param, object=ps, /quiet
  pc_read_dim, object=pd, /quiet

  pc_read_pvar, object=object, varfile=pvar
  n = object.npar_found

  nx = pd.nxgrid  & ny = pd.nygrid & nz = pd.nzgrid
  x_0 = ps.xyz0[0L] & y_0 = ps.xyz0[1L] & z_0 = ps.xyz0[2L]
  x_1 = ps.xyz1[0L] & y_1 = ps.xyz1[1L] & z_1 = ps.xyz1[2L]

  xrange = [x_0, x_1]
  yrange = [y_0, y_1]
  zrange = [z_0, z_1]


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

  plot_3d = 0L
  device, get_decomposed=decomposed
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


  ;; Kinetic energy of each particle:
  v2 = object.vpx ^ 2 + object.vpy ^ 2 + object.vpz ^ 2
  v2 = v2 / max(v2) * (ncolors - 1L)
  color = byte(temporary(v2) + bottom)

  if Lxw_ gt 1L || Lyw_ gt 1L || Lzw_ gt 1L then begin

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
    endelse

    xmargin = [0.3, 0.3] & xstyle = 1 + 4 & xtickformat = '(a1)'
    ymargin = [0.3, 0.3] & ystyle = 1 + 4 & ytickformat = '(a1)'

    plot, [0], /nodata, background=colorbg, color=colorfg, $
          xmargin=xmargin, xrange=xrange, xstyle=xstyle, xthick=thick, xtickformat=xtickformat, $
          ymargin=ymargin, yrange=yrange, ystyle=ystyle, ythick=thick, ytickformat=ytickformat

    i = - 1L
    while 1 do begin
      i += long(1d0 / fraction)
      if i ge n then break
      if object_cp[i] lt vi_low || object_cp[i] gt vi_hig then continue
      plots, object_xp[i], object_yp[i], psym=psym, symsize=symsize, color=color[i]
    endwhile
     
  endif else if nx gt 1L && ny gt 1L && nz gt 1L then begin

    object_xp = object.xp
    object_yp = object.yp
    object_zp = object.zp

    xtickformat = '(a1)'
    ytickformat = '(a1)'
    ztickformat = '(a1)'

    surface, [[0.0,0.0,0.0], [0.0,0.0,0.0]], /nodata, ax=ax, az=az, $
             background=colorbg, color=colorfg, /save, $
             xrange=xrange, /xstyle, xthick=thick, xtickformat=xtickformat, $
             yrange=yrange, /ystyle, ythick=thick, ytickformat=ytickformat, $
             zrange=zrange, /zstyle, zthick=thick, ztickformat=ztickformat

    axis, xaxis=1, x_0, y_1, z_1, /t3d, xthick=thick, xtickformat=xtickformat, color=colorfg
    axis, xaxis=1, x_0, y_1, z_0, /t3d, xthick=thick, xtickformat=xtickformat, color=colorfg
    axis, yaxis=1, x_1, y_0, z_0, /t3d, ythick=thick, ytickformat=ytickformat, color=colorfg
    axis, yaxis=1, x_1, y_0, z_1, /t3d, ythick=thick, ytickformat=ytickformat, color=colorfg
    axis, zaxis=0, x_1, y_1, z_0, /t3d, zthick=thick, ztickformat=ztickformat, color=colorfg
    axis, zaxis=0, x_1, y_0, z_0, /t3d, zthick=thick, ztickformat=ztickformat, color=colorfg

    i = - 1L
    while 1 do begin
      i += 1d0 / fraction
      if i ge n then break
      plots, object_xp[i], object_yp[i], object_zp[i], /t3d, $
             psym=psym, symsize=symsize, color=color[i]
    endwhile

    axis, zaxis=1, x_0, y_0, z_0, /t3d, zthick=thick, ztickformat=ztickformat,col=colorfg
    axis, yaxis=0, x_0, y_0, z_1, /t3d, ythick=thick, ytickformat=ytickformat,col=colorfg
    axis, xaxis=0, x_0, y_0, z_1, /t3d, xthick=thick, xtickformat=xtickformat,col=colorfg

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
  endif else begin
    wset, !d.window
    device, copy=copy
    wdelete, winid
    device, decomposed=1L
    image = tvrd(/true)
    device, decomposed=decomposed
  endelse

  ;; Save the file:
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

  device, decomposed=decomposed


  return
end ;;; procedure: pc_scatter_particles
