;
;  Plot plane of 2D-averages.
;
pro pc_read_aver_plot_plane, array_plot=array_plot, nxg=nxg, nyg=nyg, $
    min=min, max=max, zoom=zoom, xax=xax, yax=yax, $
    xtitle=xtitle, ytitle=ytitle, title=title, $
    subbox=subbox, subcen=subcen, $
    subpos=subpos, rsubbox=rsubbox, subcolor=subcolor, tsubbox=tsubbox,$
    submin=submin, submax=submax, sublog=sublog, $
    noaxes=noaxes, thick=thick, charsize=charsize, $
    loge=loge, log10=log10, $
    t_title=t_title, t_scale=t_scale, t_zero=t_zero, interp=interp, $
    ceiling=ceiling, position=position, fillwindow=fillwindow, $
    tformat=tformat, xxstalk=xxstalk, yystalk=yystalk, $
    tmin=tmin, ps=ps, png=png, imgdir=imgdir, noerase=noerase, $
    xsize=xsize, ysize=ysize, lwindow_opened=lwindow_opened, $
    colorbar=colorbar, bartitle=bartitle, quiet=quiet, $
    x0=x0, x1=x1, y0=y0, y1=y1, Lx=Lx, Ly=Ly, time=time, itimg=itimg, wait=wait
;
;  Define line and character thickness (will revert to old settings later).
;
  oldthick=thick
  !p.charthick=thick
  !p.thick=thick
  !x.thick=thick
  !y.thick=thick
;
  if (fillwindow) then position=[0.1,0.1,0.9,0.9]
;
;  Possible to zoom.
;
  if (zoom ne 1.0) then begin
    array_plot=congrid(array_plot,zoom*nxg,zoom*nyg,interp=interp,/center)
    xax=congrid(xax,zoom*nxg,interp=interp,/center)
    yax=congrid(yax,zoom*nyg,interp=interp,/center)
  endif
;
;  Take logarithm.
;
  if (loge)  then array_plot=alog(array_plot)
  if (log10) then array_plot=alog10(array_plot)
;
;  Set the maximum value of the array.
;
  if (ceiling) then begin
    ii=where(array_plot gt ceiling)
    if (max(ii) ne -1) then array_plot[ii]=ceiling
  endif
;  Plot to post script (eps).
  if (ps) then begin
    set_plot, 'ps'
    imgname='img_'+strtrim(string(itimg,'(i20.4)'),2)+'.eps'
    device, filename=imgdir+'/'+imgname, xsize=xsize, ysize=ysize, $
        color=1, /encapsulated, bits_per_pixel=8
    ps_fonts
    !p.font=0
  endif else if (png) then begin
;  Plot to png.
  endif else begin
;  Plot to X.
    if (not noerase) then begin
      window, retain=2, xsize=xsize > zoom*nxg, ysize=ysize > zoom*nyg
    endif else begin
      if (not lwindow_opened) then $
          window, retain=2, xsize=xsize > zoom*nxg, ysize=ysize > zoom*nyg
      lwindow_opened=1
    endelse
  endelse
  sym=texsyms()
;  Put current time in title if requested.
  if (t_title) then begin
    timestr='t='+strtrim(string(time/t_scale-t_zero,format=tformat),2)
    title = strtrim(title,2) eq '' ? timestr : title+'('+timestr+')'
  endif
;  tvscl-type plot with axes.
  if min gt max then begin
    min=min(array_plot) & max=max(array_plot)
  endif
  plotimage, array_plot, $
      range=[min, max], imgxrange=[x0,x1], imgyrange=[y0,y1], $
      xtitle=xtitle, ytitle=ytitle, title=title, $
      position=position, noerase=noerase, noaxes=noaxes, $
      interp=interp, charsize=charsize, thick=thick
;
;  Enlargement of ``densest'' point.
;
  if ( subbox and (time ge tsubbox) ) then begin
    if (subcen[0] eq -1) then begin
      isub=where(array_plot eq max(array_plot))
      isub=array_indices(array_plot,isub)
    endif else begin
      isub=subcen
    endelse
;  Plot box indicating enlarged region.
    oplot, [xax[isub[0]]-rsubbox,xax[isub[0]]+rsubbox, $
            xax[isub[0]]+rsubbox,xax[isub[0]]-rsubbox, $
            xax[isub[0]]-rsubbox], $
           [yax[isub[1]]-rsubbox,yax[isub[1]]-rsubbox, $
            yax[isub[1]]+rsubbox,yax[isub[1]]+rsubbox, $
            yax[isub[1]]-rsubbox], color=subcolor, thick=thick
;  Box crosses lower boundary.
    if (yax[isub[1]]-rsubbox lt y0) then begin
      oplot, [xax[isub[0]]-rsubbox,xax[isub[0]]+rsubbox, $
              xax[isub[0]]+rsubbox,xax[isub[0]]-rsubbox, $
              xax[isub[0]]-rsubbox], $
             [yax[isub[1]]+Ly-rsubbox,yax[isub[1]]+Ly-rsubbox, $
              yax[isub[1]]+Ly+rsubbox,yax[isub[1]]+Ly+rsubbox, $
              yax[isub[1]]+Ly-rsubbox], color=subcolor, thick=thick
    endif
;  Box crosses upper boundary.
    if (yax[isub[1]]+rsubbox gt y1) then begin
      oplot, [xax[isub[0]]-rsubbox,xax[isub[0]]+rsubbox, $
              xax[isub[0]]+rsubbox,xax[isub[0]]-rsubbox, $
              xax[isub[0]]-rsubbox], $
             [yax[isub[1]]-Ly-rsubbox,yax[isub[1]]-Ly-rsubbox, $
              yax[isub[1]]-Ly+rsubbox,yax[isub[1]]-Ly+rsubbox, $
              yax[isub[1]]-Ly-rsubbox], thick=thick
    endif
;  Subplot and box.
    if ( (xax[isub[0]]-rsubbox lt x0) or $
         (xax[isub[0]]+rsubbox gt x1) ) then begin
      array_plot=shift(array_plot,[nxg/2,0])
      if (xax[isub[0]]-rsubbox lt x0) then isub[0]=(isub[0]+nxg/2) mod nxg
      if (xax[isub[0]]+rsubbox gt x1) then isub[0]=(isub[0]-nxg/2) mod nxg
    endif
    if ( (yax[isub[1]]-rsubbox lt y0) or $
         (yax[isub[1]]+rsubbox gt y1) ) then begin
      array_plot=shift(array_plot,[0,nyg/2])
      if (yax[isub[1]]-rsubbox lt y0) then isub[1]=(isub[1]+nyg/2) mod nyg
      if (yax[isub[1]]+rsubbox gt y1) then isub[1]=(isub[1]-nyg/2) mod nyg
    endif
    if (sublog) then array_plot=alog10(array_plot)

    xrange=xax[isub[0]]+[-rsubbox,rsubbox]
    yrange=yax[isub[1]]+[-rsubbox,rsubbox]

    reset=0
    if submin gt submax then begin
      sumin=min & submax=max
      indsx=where(xax[isub[0]] ge xrange[0] and xax[isub[0]] le xrange[1], countx)
      if countx gt 0 then begin
        indsy=where(yax[isub[1]] ge yrange[0] and yax[isub[1]] le yrange[1], county)
        if county gt 0 then $
          submin=min(array_plot[indsx[0]:indsx[countx-1],indsy[0]:indsy[county-1]],max=submax) 
      endif
      reset=1
    endif
    plotimage, array_plot, $
        xrange=xrange, yrange=yrange, $
        range=[submin,submax], imgxrange=[x0,x1], imgyrange=[y0,y1], $
        position=subpos, /noerase, /noaxes, $
        interp=interp, charsize=charsize, thick=thick
    if reset then begin
      submin=1. & submax=0.
    endif
    plots, [subpos[0],subpos[2],subpos[2],subpos[0],subpos[0]], $
           [subpos[1],subpos[1],subpos[3],subpos[3],subpos[1]], /normal, $
           thick=thick, color=subcolor
  endif
;  Overplot stalked particles.
  if (n_elements(xxstalk) ne 0) then oplot, [xxstalk], [yystalk], ps=1, color=255, thick=2.0
;  Colorbar indicating range.
  if (colorbar) then begin
    colorbar_co, range=[min,max], pos=[0.89,0.15,0.91,0.35], divisions=1, $
        title=bartitle, /normal, /vertical
  endif
;  For png output, take image from z-buffer.
  if (png) then begin
    image = tvrd()
    tvlct, red, green, blue, /get
    imgname='img_'+strtrim(string(itimg,'(i20.4)'),2)+'.png'
    write_png, imgdir+'/'+imgname, image, red, green, blue
  endif else if (ps) then $
;  Close postscript device.
    device, /close $
  else $
    wait, wait

  itimg=itimg+1
  if (ps or png and not quiet) then print, 'Written image '+imgdir+'/'+imgname
;
;  Revert to old settings for line and character thickness.
;
  thick=oldthick
  !p.charthick=thick
  !p.thick=thick
  !x.thick=thick
  !y.thick=thick
;
end

