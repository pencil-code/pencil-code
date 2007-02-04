;; $Id: pc_read_yaver.pro,v 1.2 2007-02-04 21:19:19 brandenb Exp $
;;
;;   Read y-averages from file.
;;   Default is to only plot the data (with tvscl), not to save it in memory.
;;   The user can get the data returned in an object by specifying nit, the
;;   number of snapshots to save.
;;
;;   27-01-07/axel: adapted from pc_read_zaver, but doesn't seem to work yet
;;
pro pc_read_yaver, object=object, varfile=varfile, datadir=datadir, $
    nit=nit, lplot=lplot, iplot=iplot, min=min, max=max, zoom=zoom, $
    xax=xax, zax=zax, xtitle=xtitle, ytitle=ytitle, title=title, $
    lsubbox=lsubbox, rsubbox=rsubbox, noaxes=noaxes, $
    t_title=t_title, t_scale=t_scale, t_zero=t_zero, $
    position=position, fillwindow=fillwindow, tformat=tformat, $
    tmin=tmin, njump=njump, ps=ps, png=png, imgdir=imgdir, noerase=noerase, $
    xsize=xsize, ysize=ysize, it1=it1, nt=nt, quiet=quiet
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;;
;;  Default values.
;;
default, datadir, './data'
default, varfile, 'yaverages.dat'
default, nit, 0
default, lplot, 0 ;(changed default to 0/axel)
default, iplot, 0
default, zoom, 1
default, min, 0.0
default, max, 1.0
default, tmin, 0.0
default, njump, 0
default, ps, 0
default, png, 0
default, noerase, 0
default, imgdir, '.'
default, xsize, 10.0
default, ysize, 10.0
default, title, ''
default, t_title, 0
default, t_scale, 1.0
default, t_zero, 0.0
default, lsubbox, 0
default, rsubbox, 5
default, fillwindow, 0
default, tformat, '(f5.1)'
default, it1, 10
default, quiet, 0
;
if (fillwindow) then position=[0.1,0.1,0.9,0.9]
;;
;;  Get necessary dimensions.
;;
pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
pc_set_precision, dim=dim, quiet=quiet
;;
;;  Derived dimensions.
;;
nx=dim.nx
nz=dim.nz
;;
;;  Define axes (default to indices if no axes are supplied).
;;
if (n_elements(xax) eq 0) then xax=findgen(nx)
if (n_elements(zax) eq 0) then zax=findgen(nz)
x0=xax[0] & x1=xax[nx-1] & z0=zax[0] & z1=zax[nz-1]
Lx=xax[nx-1]-xax[0] & Lz=zax[nz-1]-zax[0]
;;
;;  Read variables from yaver.in
;;  Note however that nvar should be equal to nnamexz in data/index.pro
;;  This will no be the case when more than mnamexz variables are read.
;;
spawn, "echo "+datadir+" | sed -e 's/data\/*$//g'", datatopdir
spawn, 'cat '+datatopdir+'/yaver.in', varnames
if (not quiet) then print, 'Preparing to read y-averages ', $
    arraytostring(varnames,quote="'",/noleader)
nvar=n_elements(varnames)
;;  Die if attempt to plot variable that does not exist.
if (iplot gt nvar-1) then message, 'iplot must not be greater than nvar-1!'
;;
;;  Define arrays to put data in.
;;  The user must supply the length of the time dimension (nit).
;;
if (nit gt 0) then begin

  if (not quiet) then print, 'Returning averages at ', strtrim(nit,2), ' times'

  tt=fltarr(nit)*one
  for i=0,nvar-1 do begin
    cmd=varnames[i]+'=fltarr(nx,nz,nit)*one'
    if (execute(cmd,0) ne 1) then message, 'Error defining data arrays'
  endfor

endif
;;
;;  Variables to put single time snapshot in.
;;
array=fltarr(nx,nz,nvar)*one
help,array
t  =0.0*one
;;
;;  Prepare for read
;;
GET_LUN, file
filename=datadir+'/'+varfile 
if (not quiet) then print, 'Reading ', filename
dummy=findfile(filename, COUNT=countfile)
if (not countfile gt 0) then begin
  print, 'ERROR: cannot find file '+ filename
  stop
endif
close, file
openr, file, filename, /f77
;;
;; For png output, open z buffer already now.
;;
if (png) then begin
  set_plot, 'z'
  device, set_resolution=[zoom*nx,zoom*nz]
  print, 'Deleting old png files in the directory ', imgdir
  spawn, '\rm -f '+imgdir+'/img_*.png'
endif
;;
;;  Read y-averages and put in arrays if requested.
;;
it=0 & itimg=0
lwindow_opened=0
while (not eof(file)) do begin

  readu, file, t
  if ( (t ge tmin) and (it mod njump eq 0) ) then begin
    readu, file, array
;;
;;  Plot requested variable, variable number 0 by default.
;;
    if (lplot) then begin
      array_plot=array[*,*,iplot]
;;  Plot to post script (eps).      
      if (ps) then begin
        set_plot, 'ps'
        imgname='img_'+strtrim(string(itimg,'(i20.4)'),2)+'.eps'
        device, filename=imgdir+'/'+imgname, xsize=xsize, ysize=ysize, $
            color=1, /encapsulated
      endif else if (png) then begin
;;  Plot to png.
      endif else begin
;;  Plot to X.
        if (not noerase) then begin
          window, retain=2, xsize=zoom*nx, ysize=zoom*nz
        endif else begin
          if (not lwindow_opened) then $
              window, retain=2, xsize=zoom*nx, ysize=zoom*nz
          lwindow_opened=1
        endelse
      endelse
;;  Put current time in title if requested.      
      if (t_title) then $
          title='t='+strtrim(string(t/t_scale-t_zero,format=tformat),2)
;;  tvscl-type plot with axes.        
      plotimage, array_plot, $
          range=[min, max], imgxrange=[x0,x1], imgyrange=[z0,z1], $
          xtitle=xtitle, ytitle=ytitle, title=title, $
          position=position, noerase=noerase, noaxes=noaxes
;;
;;  Enlargement of ``densest'' point.          
;;
      if (lsubbox) then begin
        imax=where(array_plot eq max(array_plot))
        imax=array_indices(array_plot,imax)
        subpos=[0.7,0.7,0.9,0.9]
;;  Plot box indicating enlarged region.
        oplot, [xax[imax[0]]-rsubbox,xax[imax[0]]+rsubbox, $
                xax[imax[0]]+rsubbox,xax[imax[0]]-rsubbox, $
                xax[imax[0]]-rsubbox], $
               [zax[imax[1]]-rsubbox,zax[imax[1]]-rsubbox, $
                zax[imax[1]]+rsubbox,zax[imax[1]]+rsubbox, $
                zax[imax[1]]-rsubbox]
;;  Box crosses lower boundary.
        if (zax[imax[1]]-rsubbox lt z0) then begin
          oplot, [xax[imax[0]]-rsubbox,xax[imax[0]]+rsubbox, $
                  xax[imax[0]]+rsubbox,xax[imax[0]]-rsubbox, $
                  xax[imax[0]]-rsubbox], $
                 [zax[imax[1]]+Lz-rsubbox,zax[imax[1]]+Lz-rsubbox, $
                  zax[imax[1]]+Lz+rsubbox,zax[imax[1]]+Lz+rsubbox, $
                  zax[imax[1]]+Lz-rsubbox]
        endif
;;  Box crosses upper boundary.
        if (zax[imax[1]]+rsubbox gt z1) then begin
          oplot, [xax[imax[0]]-rsubbox,xax[imax[0]]+rsubbox, $
                  xax[imax[0]]+rsubbox,xax[imax[0]]-rsubbox, $
                  xax[imax[0]]-rsubbox], $
                 [zax[imax[1]]-Lz-rsubbox,zax[imax[1]]-Lz-rsubbox, $
                  zax[imax[1]]-Lz+rsubbox,zax[imax[1]]-Lz+rsubbox, $
                  zax[imax[1]]-Lz-rsubbox]
        endif
;;  Subplot and box.
        if ( (xax[imax[0]]-rsubbox lt x0) or $
             (xax[imax[0]]+rsubbox gt x1) ) then begin
          array_plot=shift(array_plot,[nx/2,0])
          if (xax[imax[0]]-rsubbox lt x0) then imax[0]=imax[0]+nx/2
          if (xax[imax[0]]+rsubbox gt x1) then imax[0]=imax[0]-nx/2
        endif
        if ( (zax[imax[1]]-rsubbox lt z0) or $
             (zax[imax[1]]+rsubbox gt z1) ) then begin
          array_plot=shift(array_plot,[0,nz/2])
          if (zax[imax[1]]-rsubbox lt z0) then imax[1]=imax[1]+nz/2
          if (zax[imax[1]]+rsubbox gt z1) then imax[1]=imax[1]-nz/2
        endif
        plotimage, array_plot, $
            xrange=xax[imax[0]]+[-rsubbox,rsubbox], $
            yrange=zax[imax[1]]+[-rsubbox,rsubbox], $
            range=[min,max], imgxrange=[x0,x1], imgyrange=[z0,z1], $
            position=subpos, /noerase, /noaxes
        plots, [subpos[0],subpos[2],subpos[2],subpos[0],subpos[0]], $
               [subpos[1],subpos[1],subpos[3],subpos[3],subpos[1]], /normal
      endif

;;  For png output, take image from z-buffer.          
      if (png) then begin
        image = tvrd()
        tvlct, red, green, blue, /get
        imgname='img_'+strtrim(string(itimg,'(i20.4)'),2)+'.png'
        write_png, imgdir+'/'+imgname, image, red, green, blue
      endif
;;  Close postscript device.      
      if (ps) then begin
        device, /close
      endif
      itimg=itimg+1
      if (ps or png and not quiet) $
          then print, 'Written image '+imgdir+'/'+imgname
    endif
;;
;;  Diagnostics.
;;
    if (it mod it1 eq 0) then begin
      if (it eq 0 ) then $
          print, '  ------- it ------- ivar -------- t --------- min(var) ------- max(var) -----'
      for ivar=0,nvar-1 do begin
          print, it, ivar, t, min(array[*,*,ivar]), max(array[*,*,ivar])
      endfor
    endif
;;
;;  Split read data into named arrays.
;;
    if ( it le nit-1 ) then begin
      tt[it]=t
      for ivar=0,nvar-1 do begin
        cmd=varnames[ivar]+'[*,*,it]=array[*,*,ivar]'
        if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
      endfor
    endif
  endif else begin
    readu, file, dummy
  endelse
;;
    it=it+1
;;
endwhile
nt=it
;;
;;  Put data in structure.
;;
if (nit ne 0) then begin
  makeobject="object = CREATE_STRUCT(name=objectname,['t'," + $
      arraytostring(varnames,QUOTE="'",/noleader) + "]," + $
      "tt,"+arraytostring(varnames,/noleader) + ")"
;
  if (execute(makeobject) ne 1) then begin
    message, 'ERROR Evaluating variables: ' + makeobject, /INFO
    undefine,object
  endif
endif
;
end
