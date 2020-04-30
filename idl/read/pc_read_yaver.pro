;
; $Id$
;
;   Read y-averages from file.
;
;   Default is to only plot the data (with tvscl), not to save it in memory.
;   The user can get the data returned in an object by specifying nit, the
;   number of snapshots to save.
;
;  We start with a script for plotting data plane - the main program follows
;  below.
;
pro pc_read_yaver, object=object, varfile=varfile, datadir=datadir, dim=dim, grid=grid, $
    nit=nit, iplot=iplot, min=min, max=max, zoom=zoom, xax=xax, zax=zax, $
    ipxread=ipxread, ipzread=ipzread, $
    xtitle=xtitle, ztitle=ztitle, title=title, subbox=subbox, subcen=subcen, $
    subpos=subpos, rsubbox=rsubbox, subcolor=subcolor, tsubbox=tsubbox,$
    submin=submin, submax=submax, sublog=sublog, $
    noaxes=noaxes, thick=thick, charsize=charsize, loge=loge, log10=log10, $
    t_title=t_title, t_scale=t_scale, t_zero=t_zero, interp=interp, $
    ceiling=ceiling, position=position, fillwindow=fillwindow, $
    tformat=tformat, stalk=stalk, nstalk=nstalk, swap_endian=swap_endian, $
    tmin=tmin, njump=njump, ps=ps, png=png, imgdir=imgdir, noerase=noerase, $
    xsize=xsize, zsize=zsize, it1=it1, variables=variables, $
    colorbar=colorbar, bartitle=bartitle, xshift=xshift, timefix=timefix, $
    readpar=readpar, readgrid=readgrid, debug=debug, quiet=quiet, wait=wait, write=write
;
COMPILE_OPT IDL2,HIDDEN
;
  pc_read_2d_aver, 'y', object=object, varfile=varfile, datadir=datadir, dim=dim, grid=grid, $
      nit=nit, iplot=iplot, min=min, max=max, zoom=zoom, xax=xax, yax=zax, $
      ipxread=ipxread, ipyread=ipzread, $
      xtitle=xtitle, ytitle=ztitle, title=title, subbox=subbox, subcen=subcen, $
      subpos=subpos, rsubbox=rsubbox, subcolor=subcolor, tsubbox=tsubbox,$
      submin=submin, submax=submax, sublog=sublog, $
      noaxes=noaxes, thick=thick, charsize=charsize, loge=loge, log10=log10, $
      t_title=t_title, t_scale=t_scale, t_zero=t_zero, interp=interp, $
      ceiling=ceiling, position=position, fillwindow=fillwindow, $
      tformat=tformat, stalk=stalk, nstalk=nstalk, swap_endian=swap_endian, $
      tmin=tmin, njump=njump, ps=ps, png=png, imgdir=imgdir, noerase=noerase, $
      xsize=xsize, ysize=zsize, it1=it1, variables=variables, $
      colorbar=colorbar, bartitle=bartitle, xshift=xshift, timefix=timefix, $
      readpar=readpar, readgrid=readgrid, debug=debug, quiet=quiet, wait=wait, write=write
;
end
