;
;  $Id$
;+
;  Reads in HDF5 slice files as they are generated by the pencil code.
;  The variable "field" can be changed. Default is 'lnrho'.
;
;  If the keyword '/mpeg' is given, the file movie.mpg is written.
;  'tmin' is the time after which data are written
;  'nrepeat' is the number of repeated images (to slow down movie)
;  An alternative is to set the '/png_truecolor' flag and postprocess the
;  PNG images with ${PENCIL_HOME}/utils/makemovie (requires imagemagick
;  and mencoder to be installed).
;
;  Typical calling sequence
;    rvid_box, 'bz', tmin=190, tmax=200, min=-.35, max=.35, /mpeg
;    rvid_box, 'ss', max=6.7, fo='(1e9.3)'
;
;  For 'gyroscope' look:
;    rvid_box, 'oz', tmin=-1, tmax=1001, /shell, /centred, r_int=0.5, r_ext=1.0
;
;  For centred slices, but without masking outside shell
;    rvid_box, 'oz', tmin=-1, tmax=1010, /shell, /centred, r_int=0.0, r_ext=5.0
;
;  For slice position m
;    rvid_box, field, tmin=0, min=-amax, max=amax, /centred, /shell, $
;       r_int=0., r_ext=8., /z_topbot_swap
;
;  For masking of shell, but leaving projections on edges of box
;    rvid_box, 'oz', tmin=-1, tmax=1001, /shell, r_int=0.5, r_ext=1.0
;
;  If using '/centred', the optional keywords '/z_bot_twice' and '/z_top_twice'
;  plot the bottom or top xy-planes (i.e. xy and xy2, respectively) twice.
;  (Once in the centred position, once at the bottom of the plot; cf. the
;  default, which plots the top slice in the centred position.)
;  (This can be useful for clarifying hidden features in gyroscope plots.)
;-
  pro rvid_box_hdf5, field, $
    mpeg=mpeg, png=png, truepng=png_truecolor, tmin=tmin, tmax=tmax, $
    max=amax, min=amin, noborder=noborder, imgprefix=imgprefix, imgdir=imgdir, $
    dev=dev, nrepeat=nrepeat, wait=wait, stride=stride, datadir=datadir, $
    noplot=noplot, fo=fo, swapz=swapz, xsize=xsize, ysize=ysize, $
    title=title, itpng=itpng, global_scaling=global_scaling, proc=proc, $
    exponential=exponential, sqroot=sqroot, logarithmic=logarithmic, $
    shell=shell, centred=centred, r_int=r_int, r_ext=r_ext, colmpeg=colmpeg, $
    z_bot_twice=z_bot_twice, z_top_twice=z_top_twice, $
    z_topbot_swap=z_topbot_swap, xrot=xrot, zrot=zrot, zof=zof, $
    magnify=magnify, ymagnify=ymagnify, zmagnify=zmagnify, xpos=xpos, zpos=zpos, $
    xmax=xmax, ymax=ymax, sample=sample, $
    xlabel=xlabel, ylabel=ylabel, tlabel=tlabel, label=label, $
    size_label=size_label, $
    monotonous_scaling=monotonous_scaling, symmetric_scaling=symmetric_scaling, $
    automatic_scaling=automatic_scaling, roundup=roundup, $
    nobottom=nobottom, oversaturate=oversaturate, cylinder=cylinder, $
    tunit=tunit, qswap=qswap, bar=bar, nolabel=nolabel, norm=norm, $
    divbar=divbar, blabel=blabel, bsize=bsize, bformat=bformat, thlabel=thlabel, $
    bnorm=bnorm, newwindow=newwindow, axes=axes, ct=ct, neg=neg, scale=scale, $
    colorbarpos=colorbarpos, quiet=quiet, single=single, help=help

  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
  if (keyword_set(help)) then begin
    doc_library, 'rvid_box_hdf5'
    return
  endif
;
  datadir = pc_get_datadir (datadir)
  default, amax, 0.05
  default, amin, -amax
  default, field, 'lnrho'
  default, dimfile, 'dim.dat'
  default, varfile, 'var.dat'
  default, nrepeat, 0
  default, stride, 0
  default, tmin, 0.0
  default, tmax, 1e38
  default, wait, 0.0
  default, fo, "(f6.1)"
  default, xsize, 512
  default, ysize, 448
  default, title, ''
  default, noborder, [ 0, 0, 0, 0, 0, 0 ]
  default, r_int, 0.5
  default, r_ext, 1.0
  default, imgprefix, 'img_'
  default, imgdir, '.'
  default, dev, 'x'
  default, magnify, 1.0
  default, ymagnify, 1.0
  default, zmagnify, 1.0
  default, xpos, 0.0
  default, zpos, 0.34
  default, xrot, 30.0
  default, zrot, 30.0
  default, zof, 0.7
  default, xmax, 1.0
  default, ymax, 1.0
  default, xlabel, 0.08
  default, ylabel, 1.18
  default, tlabel, '!8t'
  default, label, ''
  default, size_label, 1.4
  default, thlabel, 1.0
  default, nobottom, 0.0
  default, monotonous_scaling, 0.0
  default, oversaturate, 1.0
  default, tunit, 1.0
  default, scale, 1.0
  default, colorbarpos, [ 0.80, 0.15, 0.82, 0.85 ]
  default, norm, 1.0
  default, yinyang, 0
  default, single, 0

  if (keyword_set (png_truecolor)) then png = 1

  ; if PNGs are requested do not open a window
  if (not keyword_set (png)) then begin
    if (keyword_set (newwindow)) then window, /free, xsize=xsize, ysize=ysize, title=title
  end

  file_slice1 = field+'_xy2'+'.h5'
  file_slice2 = field+'_xy'+'.h5'
  file_slice3 = field+'_xz'+'.h5'
  file_slice4 = field+'_yz'+'.h5'
  if (keyword_set (swapz)) then begin
    ; swap z slices
    tmp = file_slice1
    file_slice1 = file_slice2
    file_slice2 = tmp
  end

  if (not file_test (datadir+'/slices/'+file_slice1)) then begin
    print, 'Slice file "'+datadir+'/slices/'+file_slice1+'" does not exist!'
    pos = strpos (file_slice3, '_xz')
    compfile = strmid (file_slice3, 0, pos)+'1_xz.h5'
    if (file_test (datadir+'/slices/'+compfile)) then print, 'Field name "'+field+'" refers to a vectorial quantity -> select component!'
    return
  end

  if (not keyword_set (quiet)) then print, 'Reading "'+datadir+'/slices/'+file_slice1+'", ...'
  last = pc_read ('last', filename=file_slice1, datadir=datadir+'/slices', /close)

  pc_read_dim, obj=dim, proc=proc, datadir=datadir, /quiet
  pc_set_precision, dim=dim, /quiet
  mx = dim.mx
  my = dim.my
  mz = dim.mz
  nx = dim.nx
  ny = dim.ny
  nz = dim.nz
  nghostx = dim.nghostx
  nghosty = dim.nghosty
  nghostz = dim.nghostz
  ncpus = dim.nprocx * dim.nprocy * dim.nprocz

  pc_read_grid, obj=grid, proc=proc, dim=dim, datadir=datadir, /quiet, single=single
  x = grid.x
  y = grid.y
  z = grid.z

  pc_read_param, obj=par, dim=dim, datadir=datadir, /quiet, single=single
  if (not all (par.lequidist)) then begin
    ; consider non-equidistant grid
    destretch = 1
    iix = spline (x[dim.l1:dim.l2], findgen (nx), par.xyz0[0] + (findgen (nx)+0.5) * (par.lxyz[0] / nx))
    iiy = spline (y[dim.m1:dim.m2], findgen (ny), par.xyz0[1] + (findgen (ny)+0.5) * (par.lxyz[1] / ny))
    iiz = spline (z[dim.n1:dim.n2], findgen (nz), par.xyz0[2] + (findgen (nz)+0.5) * (par.lxyz[2] / nz))
  end
  yinyang = par.lyinyang

  if (keyword_set (shell)) then begin
    xx = spread (x, [1,2], [my,mz])
    yy = spread (y, [0,2], [mx,mz])
    if (keyword_set (cylinder)) then begin
      zz = 0
    end else begin
      zz = spread (z, [0,1], [mx,my])
    end

    ; assume slices are all central for now -- perhaps generalize later
    ; nb: need pass these into boxbotex_scl for use after scaling of image;
    ;     otherwise pixelisation can be severe...
    ; nb: at present using the same z-value for both horizontal slices.
    ix = mx / 2
    iy = my / 2
    iz = mz / 2
    dist_xz = reform (sqrt (xx[*,iy,*]^2 + yy[*,iy,*]^2 + zz[*,iy,*]^2), nx, nz)
    dist_yz = reform (sqrt (xx[ix,*,*]^2 + yy[ix,*,*]^2 + zz[ix,*,*]^2), ny, nz)
    dist_xy = reform (sqrt (xx[*,*,iz]^2 + yy[*,*,iz]^2 + zz[*,*,iz]^2), nx, ny)
  end

  if single then begin
    t = 0.
    xy2 = fltarr (nx, ny)
    xy = fltarr (nx, ny)
    xz = fltarr (nx, nz)
    if (not yinyang) then yz = fltarr (ny, nz)
  endif else begin
    t = zero
    xy2 = make_array(nx, ny, type=type_idl)
    xy = make_array(nx, ny, type=type_idl)
    xz = make_array(nx, nz, type=type_idl)
    if (not yinyang) then yz = make_array(ny, nz, type=type_idl)
  endelse

  if (keyword_set (png)) then begin
    ; switch to Z buffer
    set_plot, 'z'
    ; set window size
    device, set_resolution=[ xsize, ysize ]
    itpng = 0
    dev = 'z'
  end else if (keyword_set (mpeg)) then begin
    ; open MPEG file
    if (!d.name eq 'X') then wdwset, 2, xs=xsize, ys=ysize
    mpeg_name = 'movie.mpg'
    print, 'write mpeg movie: ', mpeg_name
    mpegID = mpeg_open([ xsize, ysize ], filename=mpeg_name)
    itmpeg = 0
  end else begin
    ; default: screen output
    if (!d.name eq 'X') then wdwset, xs=xsize, ys=ysize
  end

  if (keyword_set (ct)) then loadct, ct

  ; set processor dimensions
  if (size (proc, /type) ne 0) then begin
    ipx = proc mod dim.nprocx
    ipy = (proc / dim.nprocx) mod dim.nprocy
    ipz = proc / (dim.nprocx * dim.nprocy)
    indices_xy = [ 0, 1 ]
    indices_xz = [ 0, 2 ]
    indices_yz = [ 1, 2 ]
    start = [ ipx*nx, ipy*ny, ipz*nz ]
    count = [ nx, ny, nz ]
    start_xy = start[indices_xy]
    count_xy = count[indices_xy]
    start_xz = start[indices_xz]
    count_xz = count[indices_xz]
    start_yz = start[indices_yz]
    count_yz = count[indices_yz]
  end

  if (keyword_set (global_scaling)) then begin
    ; find global min and max
    amax = single ? !Values.F_NaN : !Values.D_NaN*zero & amin=amax
    for pos = 1, last, stride+1 do begin
      frame = str (pos)
      xy2 = pc_read (frame+'/data', start=start_xy, count=count_xy, filename=file_slice1, datadir=datadir+'/slices', /close, single=single)
      xy = pc_read (frame+'/data', start=start_xy, count=count_xy, filename=file_slice2, datadir=datadir+'/slices', /close, single=single)
      xz = pc_read (frame+'/data', start=start_xz, count=count_xz, filename=file_slice3, datadir=datadir+'/slices', /close, single=single)
      yz = pc_read (frame+'/data', start=start_yz, count=count_yz, filename=file_slice4, datadir=datadir+'/slices', /close, single=single)
      amax = max ([ amax, max(xy2), max(xy), max(xz), max(yz) ], /NaN)
      amin = min ([ amin, min(xy2), min(xy), min(xz), min(yz) ], /NaN)
    end
    if (keyword_set (exponential)) then begin
      amax = exp (amax)
      amin = exp (amin)
    end else if (keyword_set (log)) then begin
      tiny = 1e-30
      amax = alog10 (amax)
      amin = alog10 (amin > tiny)
    end else if (keyword_set (sqroot)) then begin
      amax = sqrt (amax)
      amin = sqrt (amin)
    end
    if (not keyword_set (quiet)) then print, 'Scale using global min, max: ', amin, amax
  end

  header = 1
  for pos = 1, last, stride+1 do begin
    frame = str (pos)
    t = pc_read (frame+'/time', filename=file_slice1, datadir=datadir+'/slices', /close, single=single)
    ; skip if outside of time interval
    if (t lt tmin) then continue
    if (t gt tmax) then break

    xy2 = pc_read (frame+'/data', start=start_xy, count=count_xy, filename=file_slice1, datadir=datadir+'/slices', /close, single=single)
    xy = pc_read (frame+'/data', start=start_xy, count=count_xy, filename=file_slice2, datadir=datadir+'/slices', /close, single=single)
    xz = pc_read (frame+'/data', start=start_xz, count=count_xz, filename=file_slice3, datadir=datadir+'/slices', /close, single=single)
    yz = pc_read (frame+'/data', start=start_yz, count=count_yz, filename=file_slice4, datadir=datadir+'/slices', /close, single=single)

    if (keyword_set (neg)) then begin
      xy2 = -xy2
      xy = -xy
      xz = -xz
      yz = -yz
    end

    if (keyword_set (destretch)) then begin
      xy = interpolate (xy, iix, iiy, /grid)
      xy2 = interpolate (xy2, iix, iiy, /grid)
      xz = interpolate (xz, iix, iiz, /grid)
      yz = interpolate (yz, iiy, iiz, /grid)
    end

    ; perform preset mathematical operation on data before plotting
    if (keyword_set (sqroot)) then begin
      xy2 = sqrt (xy2)
      xy = sqrt (xy)
      xz = sqrt (xz)
      yz = sqrt (yz)
    end

    if (keyword_set (exponential)) then begin
      xy2 = exp (xy2)
      xy = exp (xy)
      xz = exp (xz)
      yz = exp (yz)
    end

    if (keyword_set (logarithmic)) then begin
      xy2 = alog (xy2)
      xy = alog (xy)
      xz = alog (xz)
      yz = alog (yz)
    end

    if (keyword_set (monotonous_scaling)) then begin
      ; increase the range
      amax1 = max ([ amax, max(xy2), max(xy), max(xz), max(yz) ])
      amin1 = min ([ amin, min(xy2), min(xy), min(xz), min(yz) ])
      amax = (4 * amax + amax1) * 0.2
      amin = (4 * amin + amin1) * 0.2
    end else if (keyword_set (automatic_scaling)) then begin
      amax = max ([ max(xy2), max(xy), max(xz), max(yz) ])
      amin = min ([ min(xy2), min(xy), min(xz), min(yz) ])
    end

    if (keyword_set (symmetric_scaling)) then begin
      ; symmetric scaling about zero
      amax = amax > abs (amin)
      amin = -amax
    end

    if (keyword_set (roundup)) then begin
      ; round up the amax and amin value
      amax = pc_round (amax)
      amin = pc_round (amin)
    end

    ; if noborder is set
    s = size (xy)
    l1 = noborder[0]
    l2 = s[1] - 1 - noborder[1]
    s = size (yz)
    m1 = noborder[2]
    m2 = s[1] - 1 - noborder[3]
    n1 = noborder[4]
    n2 = s[2] - 1 - noborder[5]

    if (keyword_set (qswap)) then begin
      ; swap xy2 and xy
      xy2s = rotate (xy[l1:l2, m1:m2], 2)
      xys = rotate (xy2[l1:l2, m1:m2], 2)
      xzs = rotate (xz[l1:l2, n1:n2], 5)
      yzs = rotate (yz[m1:m2, n1:n2], 5)
    end else begin
      xy2s = xy2[l1:l2, m1:m2]
      xys = xy[l1:l2, m1:m2]
      xzs = xz[l1:l2, n1:n2]
      yzs = yz[m1:m2, n1:n2]
    end

    if (keyword_set (noplot)) then begin
      ; output min and max without plotting
      if (keyword_set (header)) then print, '       it        t           min          max'
      header = 0
      print, pos, t, $
          min ([ min(xy2), min(xy), min(xz), min(yz) ]), $
          max ([ max(xy2), max(xy), max(xz), max(yz) ]), format='(i9, e12.4, 2f13.7)'
    end else begin
      ; plot normal box
      if (not keyword_set (shell)) then begin
        boxbotex_scl, xy2s, xys, xzs, yzs, xmax, ymax, zof=zof, zpos=zpos, ip=3, $
            amin=amin/oversaturate, amax=amax/oversaturate, dev=dev, $
            xpos=xpos, magnify=magnify, ymagnify=ymagnify, zmagnify=zmagnify, scale=scale, $
            nobottom=nobottom, norm=norm, xrot=xrot, zrot=zrot, sample=sample
        if (keyword_set (nolabel)) then begin
          if (label ne '') then xyouts, xlabel, ylabel, label, col=1, siz=size_label, charthick=thlabel
        end else begin
          if (label eq '') then begin
            xyouts, xlabel, ylabel, tlabel+'!3='+string(t/tunit, fo=fo)+'!c!6'+title, col=1, siz=size_label, charthick=thlabel
          end else begin
            xyouts, xlabel, ylabel, label+'!c!8t!3='+string(t, fo=fo)+'!c!6'+title, col=1, siz=size_label, charthick=thlabel
          end
        end
      end else begin
        ; draw axes box
        if (keyword_set (centred)) then begin
          if (keyword_set (z_bot_twice)) then begin
            xy2s = xys
          end else if (keyword_set (z_top_twice)) then begin
            xys = xy2s
          end else if (keyword_set (z_topbot_swap)) then begin
            xys = xy2s
            xy2s = xys
          end
          boxbotex_scl, xy2s, xys, xzs, yzs, 1.0, 1.0, zof=0.36, zpos=0.25, ip=3, $
              amin=amin, amax=amax, dev=dev, $
              shell=shell, centred=centred, scale=scale, sample=sample, $
              r_int=r_int, r_ext=r_ext, zrr1=dist_xy, zrr2=dist_xy, yrr=dist_xz, xrr=dist_yz, $
              nobottom=nobottom, norm=norm, xrot=xrot, zrot=zrot
          xyouts, 0.08, 0.81, '!8t!6='+string(t/tunit, fo=fo)+'!c'+title, col=1, siz=1.6
        end else begin
          boxbotex_scl, xy2s, xys, xzs, yzs, xmax, ymax, zof=0.65, zpos=0.34, ip=3, $
              amin=amin, amax=amax, dev=dev, sample=sample, shell=shell, $
              r_int=r_int, r_ext=r_ext, zrr1=dist_xy, zrr2=dist_xy, yrr=dist_xz, xrr=dist_yz, $
              nobottom=nobottom, norm=norm, xrot=xrot, zrot=zrot
          xyouts, 0.08, 1.08, '!8t!6='+string(t/tunit, fo=fo)+'!c'+title, col=1, siz=1.6
        end
      end

      if (keyword_set (bar)) then begin
        ; draw color bar
        default, bsize, 1.5
        default, bformat, '(f5.2)'
        default, bnorm, 1.
        default, divbar, 2
        default, blabel, ''
        !p.title = blabel
        colorbar, pos=colorbarpos, color=1, div=divbar, range=[ amin, amax ]/bnorm, /right, /vertical, format=bformat, charsize=bsize, title=title
        !p.title=''
      end

      if (keyword_set (axes)) then begin
        ; draw axes
        xx = !d.x_size
        yy = !d.y_size
        aspect_ratio = (1.0 * yy) / xx
        ; length of the arrow
        length = 0.1
        xlength = length
        ylength = xlength / aspect_ratio
        ; rotation angles. This 0.7 is an ugly hack that looks good for most angles.
        gamma = 0.7 * xrot * !pi/180.0
        alpha = zrot * !pi/180.0
        ; position of the origin
        x0 = 0.12
        y0 = 0.25

        ; x arrow
        x1 = x0 + xlength * (cos(gamma) * cos(alpha))
        y1 = y0 + ylength * (sin(gamma) * sin(alpha))
        angle = atan ((y1-y0) / (x1-x0))
        if ((x1-x0 le 0) and (y1-y0 ge 0)) then angle = angle + !pi
        if ((x1-x0 le 0) and (y1-y0 le 0)) then angle = angle - !pi
        x2 = x0 + length * cos(angle)
        y2 = y0 + length * sin(angle)
        arrow, x0, y0, x2, y2, col=1, /normal, thick=thlabel, hthick=thlabel
        xyouts, x2-0.01, y2-0.045, '!8x!x', col=1, /normal, siz=size_label, charthick=thlabel

        ; y arrow
        x1 = x0 + xlength * (-cos(gamma) * sin(alpha))
        y1 = y0 + ylength * ( sin(gamma) * cos(alpha))
        angle = atan ((y1-y0) / (x1-x0))
        if ((x1-x0 le 0) and (y1-y0 ge 0)) then angle = angle + !pi
        if ((x1-x0 le 0) and (y1-y0 le 0)) then angle = angle - !pi
        x2 = x0 + length * cos(angle)
        y2 = y0 + length * sin(angle)
        arrow, x0, y0, x2, y2, col=1, /normal, thick=thlabel, hthick=thlabel
        xyouts, x2-0.03, y2-0.01, '!8y!x', col=1, /normal, siz=size_label, charthick=thlabel

        ; z arrow
        x1 = x0
        y1 = y0 + ylength
        arrow, x0, y0, x1, y1, col=1, /normal, thick=thlabel, hthick=thlabel
        xyouts, x1-0.015, y1+0.01, '!8z!x', col=1, /normal, siz=size_label, charthick=thlabel
      end

      if (keyword_set (png)) then begin
        ; save as png file
        istr2 = string (itpng, '(I4.4)') ; show maximum 9999 frames
        image = tvrd ()
        ; make background white, and write png file
        bad = where (image eq 0, num)
        if (num gt 0) then image[bad] = 255
        tvlct, red, green, blue, /GET
        imgname = imgprefix+istr2+'.png'
        write_png, imgdir+'/'+imgname, image, red, green, blue
        if (keyword_set (png_truecolor)) then spawn, 'mogrify -type TrueColor ' + imgdir+'/'+imgname
        itpng++
      end else if (keyword_set (mpeg)) then begin
        ; write mpeg file directly
        ; NOTE: for idl_5.5 and later this requires the mpeg license
        image = tvrd (true=1)
        if (keyword_set (colmpeg)) then begin
          ; it seems we need 24-bit to get color MPEGs on some machines
          image24 = bytarr (3, xsize, ysize)
          tvlct, red, green, blue, /GET
        end
        for irepeat=0, nrepeat do begin
          if (keyword_set (colmpeg)) then begin
            image24[0, *, *] = red (image[0, *, *])
            image24[1, *, *] = green (image[0, *, *])
            image24[2, *, *] = blue (image[0, *, *])
            mpeg_put, mpegID, image=image24, frame=itmpeg, /order
          end else begin
            mpeg_put, mpegID, window=2, frame=itmpeg, /order
          end
        end

        if (not keyword_set (quiet)) then begin
          if (keyword_set (header)) then print, '       it        t        min/norm     max/norm         amin         amax'
          header = 0
          print, pos, t, $
              min ([ min(xy2), min(xy), min(xz), min(yz) ]) / norm, $
              max ([ max(xy2), max(xy), max(xz), max(yz) ]) / norm, $
              amin, amax
        end
        itmpeg++
      end else begin
        ; default: output on the screen
        if (not keyword_set (quiet)) then begin
          if (keyword_set (header)) then   print, '       it        t        min/norm     max/norm         amin         amax'
          header = 0
          print, pos, t, $
              min ([ min(xy2), min(xy), min(xz), min(yz) ]) / norm, $
              max ([ max(xy2), max(xy), max(xz), max(yz) ]) / norm, $
              amin, amax, format='(i9, e12.4, 4f13.7)'
        end
      end

      ; wait in case movie runs to fast
      wait, wait

      ; check whether file has been written
      if (not keyword_set (quiet) and keyword_set (png)) then spawn, 'ls -l '+imgdir+'/'+imgname

    end
  end

  if (not keyword_set (quiet)) then begin
    ; inform the user of why program stopped
    if (t gt tmax) then begin
      print, 'Stopping since t>', str (tmax)
    end else begin
      print, 'Read last slice at t=', str (t)
    end
  end

  if (keyword_set (mpeg)) then begin
    ; write and close MPEG file
    print, 'Writing MPEG file...'
    mpeg_save, mpegID, filename=mpeg_name
    mpeg_close, mpegID
  end

end
