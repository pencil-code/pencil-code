;
; $Id$
;+
;  Reads and displays data in a plane (currently with tvscl) and plots a
;  curve as well (cross-section through iy).
;
;  If the keyword /mpeg is given, the file movie.mpg is written.
;
;  tmin is the time after which data are written.
;
;  nrepeat is the number of repeated images (to slow down movie)
;
;  An alternative is to set the /png_truecolor flag and postprocess the
;  PNG images with ${PENCIL_HOME}/utils/makemovie (requires imagemagick
;  and mencoder to be installed).
;
;  The /polar option is for sphere/cylinder-in-a-box simulations only.
;
;  Typical calling sequence:
;    rvid_plane_hdf5, 'uz', min=-1e-1, max=1e-1
;
;  ... and for spherical slices
;    rvid_plane_hdf5, 'bb1', min=-.5, max=.5, /sph
;-
pro rvid_plane_hdf5, field, mpeg=mpeg, png=png, truepng=png_truecolor, tmin=tmin, $
    tmax=tmax, max=amax, swap_endian=swap_endian, quiet=quiet, $
    min=amin, extension=extension, nrepeat=nrepeat, wait=wait, $
    stride=stride, datadir=datadir, oldfile=oldfile, debug=debug, $
    proc=proc, ix=ix, iy=iy, ps=ps, iplane=iplane, imgdir=imgdir, $
    global_scaling=global_scaling, automatic_scaling=automatic_scaling, $
    shell=shell, r_int=r_int, r_ext=r_ext, zoom=zoom, colmpeg=colmpeg, exponential=exponential, $
    contourplot=contourplot, color=color, sqroot=sqroot, tunit=tunit, $
    nsmooth=nsmooth, cubic=cubic, textsize=textsize, help=help, $
    polar=polar, anglecoord=anglecoord, style_polar=style_polar, $
    spherical_surface=spherical_surface, nlevels=nlevels, $
    doublebuffer=doublebuffer, wsx=wsx, wsy=wsy, title=title, log=log, $
    interp=interp, savefile=savefile, rotate=rotate, phi_shift=phi_shift, $
    Beq=Beq, taud=taud, grid=gridplot, maptype=maptype, single=single, _extra=_extra
;
  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
;
  if (keyword_set(help)) then begin
    doc_library, 'rvid_plane_hdf5'
    return
  endif

  if (keyword_set (swap_endian)) then print, "rvid_plane_hdf5: WARNING: the 'swap_endian' parameter is ignored for HDF5 files."
;
  default, ix, -1
  default, iy, -1
  default, ps, 0
  default, quiet, 0
  default, single, 0
  default, gridplot, 0
;
;  default extension
;
  if (keyword_set (spherical_surface)) then begin
    default, extension, 'yz'
  end else begin
    default, extension, 'xz'
  end
;
  datadir = pc_get_datadir(datadir)
;
  default, amax, 0.05
  default, amin, -amax
  default, field, 'lnrho'
  default, nrepeat, 0
  default, stride, 0
  default, tmin, 0.0
  default, tmax, 1e38
  default, tunit, 1
  default, iplane, 0
  default, wait, 0.03
  default, r_int, 0.5
  default, r_ext, 1.0
  default, zoom, 1.0
  default, imgdir, '.'
  default, color, 1
  default, pixelsize, 1
  default, ximg, 1
  default, yimg, 1
  default, textsize, 1.0
  default, anglecoord, 'z'
  default, style_polar, 'fill'
  default, wsx, 640
  default, wsy, 480
  default, title, 'rvid_plane'
  default, nlevels, 30
  default, phi_shift, 0.0
  default, Beq, 1.0
  default, yinyang, 0
  default, triangles, 0
  default, maptype, 'orthographic'
;
  sample = ~keyword_set (interp)
;
  if (keyword_set (doublebuffer)) then begin
    ; set up a window for double buffering
    base = widget_base (title=title)
    draw = widget_draw (base, xsize=wsx, ysize=wsy)
    widget_control, /realize,base
    widget_control, draw, get_value=windex
  end else begin
    if (keyword_set (png_truecolor)) then png=1
  end
;
  pc_read_dim, obj=dim, proc=proc, datadir=datadir, /quiet
  pc_set_precision, dim=dim, /quiet

  nx = dim.nx
  ny = dim.ny
  nz = dim.nz
  mx = dim.mx
  my = dim.my
  mz = dim.mz
  nxgrid = dim.nxgrid
  nygrid = dim.nygrid
  nzgrid = dim.nzgrid
  mxgrid = dim.mxgrid
  mygrid = dim.mygrid
  mzgrid = dim.mzgrid
  nghostx = dim.nghostx
  nghosty = dim.nghosty
  nghostz = dim.nghostz
  nprocx = dim.nprocx
  nprocy = dim.nprocy
  nprocz = dim.nprocz
  ncpus = nprocx*nprocy*nprocz
;
  pc_read_grid, obj=grid, proc=proc, dim=dim, datadir=datadir, /quiet, single=single
  x = grid.x(dim.l1:dim.l2)
  y = grid.y(dim.m1:dim.m2)
  z = grid.z(dim.n1:dim.n2)
;
  ; adjust extension for 2D runs
  if ((nx ne 1) and (ny ne 1) and (nz eq 1)) then extension = 'xy'
  if ((nx ne 1) and (ny eq 1) and (nz ne 1)) then extension = 'xz'
  if ((nx eq 1) and (ny ne 1) and (nz ne 1)) then extension = 'yz'
;
  pc_read_param, obj=par, dim=dim, datadir=datadir, /quiet, single=single
  if (not all (par.lequidist)) then begin
    ; consider non-equidistant grid
    destretch = 1
;
    if ((nx gt 1) and not par.lequidist[0]) then begin
      x0 = 0.5 * (grid.x[dim.l1-1] + grid.x[dim.l1])
      x1 = 0.5 * (grid.x[dim.l2] + grid.x[dim.l2+1])
      dx = (x1 - x0) / nx
      iix = spline (grid.x, findgen(mx) - nghostx, x0 + (findgen(nx) + 0.5) * dx)
    end else begin
      iix = findgen (nx)
    end
;
    if ((ny gt 1) and not par.lequidist[1]) then begin
      y0 = 0.5 * (grid.y[dim.m1-1] + grid.y[dim.m1])
      y1 = 0.5 * (grid.y[dim.m2] + grid.y[dim.m2+1])
      dy = (y1 - y0) / ny
      iiy = spline (grid.y, findgen(my) - nghosty, y0 + (findgen(ny) + 0.5) * dy)
    end else begin
      iiy = findgen (ny)
    end
;
    if ((nz gt 1) and not par.lequidist[2]) then begin
      z0 = 0.5 * (grid.z[dim.n1-1] + grid.z[dim.n1])
      z1 = 0.5 * (grid.z[dim.n2] + grid.z[dim.n2+1])
      dz = (z1 - z0) / nz
      iiz = spline (grid.z, findgen(mz) - nghostz, z0 + (findgen(nz) + 0.5) * dz)
    end else begin
      iiz = findgen (nz)
    end
;
    if (extension eq 'xy') then begin
      ii1 = iix
      ii2 = iiy
    end else if (extension eq 'xz') then begin
      ii1 = iix
      ii2 = iiz
    end else begin
      ii1 = iiy
      ii2 = iiz
    end
  end
;
  yinyang = par.lyinyang
;
  file_slice = field+'_'+extension+'.h5'
;
  if (not file_test (datadir+'/slices/'+file_slice)) then begin
    print, 'Slice file "'+datadir+'/slices/'+file_slice+'" does not exist!'
    pos = strpos (file_slice, '_'+extension)
    compfile = strmid (file_slice, 0, pos)+'1_'+extension+'.h5'
    if (file_test (datadir+'/slices/'+compfile)) then print, 'Field name "'+field+'" refers to a vectorial quantity -> select component!'
    return
  end
;
  if (not quiet) then print, 'Reading "'+datadir+'/slices/'+file_slice+'"...'
  last = pc_read ('last', filename=file_slice, datadir=datadir+'/slices')
;
  case field of
    'uu1': quan = '!8u!dx!n!6'
    'uu2': quan = '!8u!dy!n!6'
    'uu3': quan = '!8u!dz!n!6'
    'bb1': quan = '!8B!dx!n!6'
    'bb2': quan = '!8B!dy!n!6'
    'bb3': quan = '!8B!dz!n!6'
    'rho': quan = '!7q!6'
    'lnrho': quan = 'ln!7q!6'
    'ss': quan = '!8s!6'
    else: quan = ''
  end
;
  if (keyword_set  (polar)) then begin
    if (yinyang) then begin
      print, 'Polar plot not meaningful for Yin-Yang data. Use contourplot or spherical_surface.'
      return
    end
    if (anglecoord eq 'y') then begin
      theta = y
    end else if (anglecoord eq 'z') then begin
      theta = z
    end
    xx = fltarr(nx,nz)
    yy = fltarr(nx,nz)
    for i = 0, nx-1 do begin
      for j = 0, nz-1 do begin
        xx[i,j] = x[i]*cos(theta[j])
        yy[i,j] = x[i]*sin(theta[j])
      end
    end
  end
;
  if (keyword_set (shell)) then begin
    ; mask the outside shell
    xx = spread (x, [1,2], [ny,nz])
    yy = spread (y, [0,2], [nx,nz])
    zz = spread (z, [0,1], [nx,ny])
    ; assume slices are all central for now - perhaps generalize later.
    ix = nx / 2
    iy = ny / 2
    iz = nz / 2
    ; nb: need pass these into boxbotex_scl for use after scaling of image; otherwise pixelisation can be severe...
    ; nb: at present using the same z-value for both horizontal slices; hardwired into boxbotex_scl, also.
    if (extension eq 'xz') then begin
      rr = reform (sqrt (xx[*,iy,*]^2 + yy[*,iy,*]^2 + zz[*,iy,*]^2), nx, nz)
    end else if (extension eq 'yz') then begin
      rr = reform (sqrt (xx[ix,*,*]^2 + yy[ix,*,*]^2 + zz[ix,*,*]^2), ny, nz)
    end else begin
      rr = reform (sqrt (xx[*,*,iz]^2 + yy[*,*,iz]^2 + zz[*,*,iz]^2), nx, ny)
    end
    dist = rebinbox (rr, zoom)
  end
;
  size_plane = [ nx, ny ]
  if (extension eq 'xz') then size_plane = [ nx, nz ]
  if (extension eq 'yz') then size_plane = [ ny, nz ]
  if (yinyang) then size_plane = [ 2, nz, nz ]
;
  dev = 'x'
  if (keyword_set (png)) then begin
    Nwx = zoom * size_plane[0]
    Nwy = zoom * size_plane[1]
    Nwy = (Nwx * 15) / 20
    resolution = [ Nwx, Nwy ]
    if (not quiet) then print, 'z-buffer resolution (in pixels, set with zoom='+str (zoom)+') = '+str (resolution[0])+' * '+str (resolution[1])
    ; switch to Z buffer
    set_plot, 'z'
    ; set window size
    device, set_resolution=resolution
    itpng = 0
    dev = 'z'
  end else if (keyword_set (mpeg)) then begin
    Nwx = zoom * size_plane[0]
    Nwy = zoom * size_plane[1]
    resolution = [ Nwx, Nwy ]
    if (not quiet) then print, 'z-buffer resolution (in pixels) = '+str (resolution[0])+' * '+str (resolution[1])
    ; switch to Z buffer
    set_plot, 'z'
    ; set window size
    device, set_resolution=resolution
    dev = 'z'
    if (!d.name eq 'X') then window, 2, xs=Nwx, ys=Nwy
    mpeg_name = 'movie.mpg'
    if (not quiet) then print, 'write mpeg movie: ', mpeg_name
    ; open MPEG file
    mpegID = mpeg_open (resolution, filename=mpeg_name)
    itmpeg = 0
  end else if (not keyword_set (doublebuffer)) then begin
    if (keyword_set (spherical_surface) and not keyword_set (contourplot)) then begin
      q = 1.0
    end else begin
      Nwx = zoom * size_plane[0]
      Nwy = zoom * size_plane[1]
      q = Nwy / Nwx
      if (yinyang) then q = 0.5
    end
    window, xsize=700, ysize=700, title=title
  end

  ; read auxiliary data for Yin-Yang grid: number of points in merged (irregular) grid ngrid and merged grid itself
  if (yinyang and (extension eq 'yz')) then begin
    ;ngrid = 0L
    ;readu, lun, ngrid
    ;yz_yy = fltarr (2, ngrid) * one
    ;readu, lun, yz_yy
    ;plane = fltarr (ngrid) * one
    ;triangulate, yz_yy[0,*], yz_yy[1,*], triangles
    ;size_plane = size (plane, /dimensions)
  end

  ; set processor dimensions
  if (is_defined(proc)) then begin
    ipx = proc mod dim.nprocx
    ipy = (proc / dim.nprocx) mod dim.nprocy
    ipz = proc / (dim.nprocx * dim.nprocy)
    indices = [ 0, 1 ]
    if (extension eq 'xz') then indices = [ 0, 2 ]
    if (extension eq 'yz') then indices = [ 1, 2 ]
    start = ([ ipx*nx, ipy*ny, ipz*nz ])[indices]
    count = ([ nx, ny, nz ])[indices]
  end

  if (not quiet) then print, 'Array size: ', size_plane, yinyang ? '(Yin-Yang grid)' : ''

  if (keyword_set (global_scaling)) then begin
    amax = !Values.F_NaN & amin=amax
    for pos = 1, last, stride+1 do begin
      frame = str (pos)
      plane = pc_read (frame+'/data', start=start, count=count, single=single)
      if (keyword_set (nsmooth)) then plane = smooth (plane, nsmooth)
      amax = max ([ amax, max (plane) ], /NaN)
      amin = min ([ amin, min (plane) ], /NaN)
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
    if (not quiet) then print, 'Scale using global min, max: ', amin, amax
  end

  for pos = 1, last, stride+1 do begin
    frame = str (pos)
    index = (pos - 1) / (stride + 1)
    plane = pc_read (frame+'/data', start=start, count=count, single=single)
    t = pc_read (frame+'/time', single=single)
;
    if (keyword_set (destretch)) then plane = interpolate(plane, ii1, ii2, /grid)
;
;  Rescale data with optional parameter zoom.
;  WARNING: the scaling can produce artifacts at shearing boundaries. Contour
;  plots give better results in that case (/contour).
;  In future, we might want to choose better names than x2 and y2,
;    especially if they are later theta (theta2) and phi.
;
    if (not (yinyang and (extension eq 'yz'))) then begin
      nx_plane = size_plane[0]
      ny_plane = size_plane[1]
    end
    if (zoom eq 1.) then begin
      if (yinyang and (extension eq 'yz')) then begin
        x2 = reform (yz_yy[0,*])
        y2 = reform (yz_yy[1,*])
      end else begin
        x2 = x
        y2 = y
      end
    end else if (extension eq 'xy') then begin
      x2 = rebin (x, zoom*nx_plane, sample=sample)
      y2 = rebin (y, zoom*ny_plane, sample=sample)
    end else if (extension eq 'yz') then begin
      if (yinyang) then begin
        ; yet incorrect
        x2 = rebin (yz_yy[0,*], zoom*ngrid, sample=sample)
        y2 = rebin (yz_yy[1,*], zoom*ngrid, sample=sample)
      end else begin
        x2 = rebin (y, zoom*nx_plane, sample=sample)
        y2 = rebin (z, zoom*ny_plane, sample=sample)
      end
    end
;
    ; our y2 is not right for extension 'xz'
    if (extension eq 'xz') then y2 = rebin (z, zoom*ny_plane, sample=sample)
;
    ; other options
    if (keyword_set (exponential)) then begin
      plane2 = rebin (exp (plane), zoom*nx_plane, zoom*ny_plane, sample=sample)
    end else if (keyword_set (nsmooth)) then begin
      plane2 = rebin (smooth (plane, nsmooth), zoom*nx_plane, zoom*ny_plane, sample=sample)
    end else if (keyword_set (sqroot)) then begin
      plane2 = rebin (sqrt (plane), zoom*nx_plane, zoom*ny_plane, sample=sample)
    end else if (keyword_set (log)) then begin
       tiny = 1e-30
       plane2 = rebin (alog10 (plane > tiny), zoom*nx_plane, zoom*ny_plane, sample=sample)
    end else if (keyword_set(cubic)) then begin
       if (cubic gt 0.0) then cubic = -0.5
       plane2 = congrid (plane, zoom*nx_plane, zoom*ny_plane, /center, cubic=cubic, interp=interp)
    end else if (zoom ne 1.) then begin
       plane2 = congrid (plane, zoom*nx_plane, zoom*ny_plane, /center, interp=interp)
    end else begin
       plane2 = plane
    end
;
    if (keyword_set (shell)) then begin
      ; do masking
      white = 255
      indices = where ((dist lt r_int) or (dist gt r_ext), num)
      if (num gt 0) then plane2[indices] = white
    end
;
    if (keyword_set (debug)) then begin
      print, t, min ([ plane2, xy, xz, yz ]), max ([ plane2, xy, xz, yz ])
    end else begin
      if ((t ge tmin) and (t le tmax)) then begin
        if (keyword_set (automatic_scaling)) then begin
          amax = max (plane2)
          amin = min (plane2)
        end
;
        ; show image scaled between amin and amax and filling whole screen
        if (keyword_set (doublebuffer)) then begin
          ; paint into buffer
          window, xsize=wsx, ysize=wsy, /pixmap, /free
          pixID = !D.Window
        end
;
        if (keyword_set (contourplot)) then begin
          lev = grange (amin, amax, 60)
          xmargin = !x.margin - [ 4, -6 ]
          title = '!8t!6 = '+strtrim (string (t/tunit, fo="(f7.1)"), 2)
          if (keyword_set (spherical_surface)) then begin
            xtitle = '!7u!6'
            ytitle = '!7h!6'
            if (yinyang) then begin
              contourfill, plane2, y2, x2, lev=lev, _extra=_extra, tri=triangles, title=title, xtitle=xtitle, ytitle=ytitle, xmar=xmargin, grid=gridplot
            end else begin
              contourfill, transpose(plane2), y2, x2, lev=lev, _extra=_extra, tri=triangles, xtitle=xtitle, ytitle=ytitle, title=title, xmar=xmargin, grid=gridplot 
            end
          end else begin
            xtitle = '!8y!6'
            ytitle = '!8z!6'
            contourfill, plane2, x2, y2, lev=lev, title=title, _extra=_extra, tri=triangles, xtitle=xtitle, ytitle=ytitle, xmar=xmargin, grid=gridplot
          end
          colorbar_co,range=[min(lev),max(lev)],pos=[0.975,0.18,0.99,0.93],/vert, yticks=4, yminor=1, char=1.5, col=0, xtit=quan, xchars=1. ;, format='(f6.4)', $;ytickv=[min(lev),0.,max(lev)]
;         end else if (keyword_set (polar)) then begin
          if (style_polar eq 'fill') then begin
            contourfill, plane2, x2, y2, levels=grange(amin,amax,60), $
                tit='!8t!6 ='+string(t/tunit,fo="(f7.1)"), _extra=_extra, grid=gridplot
          end else if (style_polar eq 'lines') then begin
            contour, plane2, x2, y2, nlevels=nlevels, $
                tit='!8t!6 ='+string(t/tunit,fo="(f7.1)"), _extra=_extra
            if (keyword_set (gridplot)) then oplot, x2, y2, psym=3
          end
;
        end else if (keyword_set (spherical_surface)) then begin
          ; spherical surface plot in a good projection still need to check whether /rotate is correct (see below)
          theta2 = x2 / !dtor
          phi = y2 / !dtor
          if (keyword_set (phi_shift)) then phi -= phi_shift
          !p.background = 255
          map_set, 15, 60, /noborder, /isotropic, latdel=15, londel=15, limit=[-80,-30,90,160], xmargin=[0.5,7], ymargin=1.5, color=0, name=maptype
;
          lev = grange (.8*amin, amax, nlevels)
          if (yinyang) then begin
            tmp = plane2
          end else begin
            if (keyword_set (rotate)) then tmp = rotate (plane2,3) else tmp = transpose (plane2)
          end
          if (keyword_set (Beq)) then tmp = tmp / Beq
          if (keyword_set (taud)) then t = t / taud
;
          mima = minmax (lev)
          contour, (tmp > mima(0)) < mima(1), phi, 90.-theta2, lev=lev, /fill, /overplot, col=0, _extra=_extra, tri=triangles
          if (keyword_set (gridplot)) then begin
            oplot, phi, 90.-theta2, psym=3
          end else begin
            map_grid, latdel=15, londel=15
          end
          colorbar_co, range=[min(lev),max(lev)], pos=[0.97,0.15,0.99,0.85], /vert, format='(F6.3)', yticks=4, yminor=1, charsize=1.5, col=0, xtit=quan, xchars=1.5 ;, ytickv=[min(lev),0.,max(lev)]
          xyouts, 0.06, 0.06, '!8t!6 = '+str(t, fo='(f5.1)')+'', col=0, /normal, charsize=2
          wait, wait
        end else begin
          tv, bytscl (plane2, min=amin, max=amax), iplane
        end
        if (keyword_set (doublebuffer)) then begin
          wset, windex
          device, copy=[ 0, 0, !D.X_Size, !D.Y_Size, 0, 0, pixID ]
          wdelete, pixID
        end
        if (keyword_set (savefile)) then begin
          if (not is_defined(slices)) then begin
            slices = plane
            times = t
          end else begin
            slices = [ [[slices]], [[plane]] ]
            times = [ times, t ]
          end
        end
        if (keyword_set (png)) then begin
          istr2 = string (itpng, '(I4.4)') ; show maximum 9999 frames
          image = tvrd ()
          ; make background white, and write png file
          tvlct, red, green, blue, /get
          imgname = imgdir+'/img_'+istr2+'.png'
          write_png, imgname, image, red, green, blue
          if (keyword_set (png_truecolor)) then spawn, 'mogrify -type TrueColor ' + imgname
          itpng += 1
        end else if (keyword_set (mpeg)) then begin
          ; write MPEG file directly
          ; for IDL v5.5+ this requires a MPEG license
          image = tvrd (true=1)
          if (keyword_set (colmpeg)) then begin
            ; ngrs seem to need to work explictly with 24-bit color to generate color MPEGs
            image24 = bytarr (3, Nwx, Nwy)
            tvlct, red, green, blue, /get
          end
          for irepeat = 0, nrepeat do begin
            if (keyword_set (colmpeg)) then begin
              image24[0,*,*] = red (image[0,*,*])
              image24[1,*,*] = green (image[0,*,*])
              image24[2,*,*] = blue (image[0,*,*])
              mpeg_put, mpegID, image=image24, frame=itmpeg, /order
            end else begin
              mpeg_put, mpegID, window=2, frame=itmpeg, /order
            end
            itmpeg += 1
          end
          if (not quiet) then print, pos, itmpeg, t, min([plane2]), max([plane2])
        end else begin
          ; default: output on the screen
          if (not quiet) then begin
            if (index eq 0) then print, '------it----------t----------min------------max--------'
            print, pos, t, min ([plane2]), max ([plane2])
          end
        end
        wait, wait
        ; check whether file has been written
        if (keyword_set (png)) then spawn, 'ls -l '+imgname
      end
    end
  end
  h5_close_file
;
  if (keyword_set (mpeg)) then begin
    ; write and close mpeg file
    if (not quiet) then print, 'Writing MPEG file...'
    mpeg_save, mpegID, filename=mpeg_name
    mpeg_close, mpegID
    set_plot, 'X'
  end
  if (keyword_set (png)) then set_plot, 'X'
  if (keyword_set (savefile)) then begin
    num_slices = 1 + (last - first) / (stride + 1)
    save, file=savefile, slices, num_slices, times, x, y, z
  end
;
END
