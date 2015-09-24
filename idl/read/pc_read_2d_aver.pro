;
; $Id$
;
;  Read 2D-averages from file.
;
pro pc_read_2d_aver, dir, object=object, varfile=varfile, datadir=datadir, $
    nit=nit, iplot=iplot, min=min, max=max, zoom=zoom, xax=xax, yax=yax, $
    ipxread=ipxread, ipyread=ipyread, $
    xtitle=xtitle, ytitle=ytitle, title=title, subbox=subbox, subcen=subcen, $
    subpos=subpos, rsubbox=rsubbox, subcolor=subcolor, tsubbox=tsubbox,$
    submin=submin, submax=submax, sublog=sublog, $
    noaxes=noaxes, thick=thick, charsize=charsize, loge=loge, log10=log10, $
    t_title=t_title, t_scale=t_scale, t_zero=t_zero, interp=interp, $
    ceiling=ceiling, position=position, fillwindow=fillwindow, $
    tformat=tformat, stalk=stalk, nstalk=nstalk, swap_endian=swap_endian, $
    tmin=tmin, njump=njump, ps=ps, png=png, imgdir=imgdir, noerase=noerase, $
    xsize=xsize, ysize=ysize, it1=it1, variables=variables, $
    colorbar=colorbar, bartitle=bartitle, xshift=xshift, timefix=timefix, $
    readpar=readpar, readgrid=readgrid, debug=debug, quiet=quiet
;
COMPILE_OPT IDL2,HIDDEN
;
  common pc_precision, zero, one
;
  if ((dir ne 'y') and (dir ne 'z')) then message, 'ERROR: direction "'+dir+'" unknown.'
;
;  Default values.
;
  if (not keyword_set(datadir)) then datadir=pc_get_datadir()
  default, varfile, dir+'averages.dat'
  default, nit, 0
  default, ipxread, -1
  default, ipyread, -1
  default, iplot, -1
  default, zoom, 1
  default, min, 0.0
  default, max, 1.0
  default, tmin, 0.0
  default, njump, 1
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
  default, stalk, 0
  default, interp, 0
  default, ceiling, 0.0
  default, subbox, 0
  default, subcen, -1
  default, subpos, [0.7,0.7,0.9,0.9]
  default, submin, min
  default, submax, max
  default, sublog, 0
  default, rsubbox, 5.0
  default, subcolor, 255
  default, tsubbox, 0.0
  default, thick, 1.0
  default, charsize, 1.0
  default, loge, 0
  default, log10, 0
  default, fillwindow, 0
  default, tformat, '(f5.1)'
  default, swap_endian, 0
  default, it1, 10
  default, variables, ''
  default, colorbar, 0
  default, bartitle, ''
  default, xshift, 0
  default, timefix, 0
  default, readpar, 0
  default, readgrid, 0
  default, debug, 0
  default, quiet, 0
  default, in_file, dir+'aver.in'
;
;  Read additional information.
;
  pc_read_dim, obj=dim, datadir=datadir, /quiet
  pc_read_dim, obj=locdim, datadir=datadir+'/proc0', /quiet
  if (stalk) then begin
    pc_read_pstalk, obj=pst, datadir=datadir, /quiet
    default, nstalk, n_elements(pst.ipar)
  endif
  pc_set_precision, dim=dim, /quiet
;
;  Need to know box size for proper axes and for shifting the shearing
;  frame.
;
  if (readpar or (xshift ne 0)) then begin
    pc_read_param, obj=par, datadir=datadir, /quiet
  endif
  if (readgrid or (xshift ne 0)) then begin
    pc_read_grid, obj=grid, /trim, datadir=datadir, /quiet
    xax=grid.x
    yax=grid.y
    if (dir eq 'y') then yax=grid.z
  endif
;
;  Some z-averages are erroneously calculated together with time series
;  diagnostics. The proper time is thus found in time_series.dat.
;
  if (timefix) then pc_read_ts, obj=ts, datadir=datadir, /quiet
;
;  Derived dimensions.
;
  nx=locdim.nx
  ny=locdim.ny
  nxgrid=dim.nx
  nygrid=dim.ny
  nprocx=dim.nprocx
  nprocy=dim.nprocy
  if (dir eq 'y') then begin
    ny=locdim.nz
    nygrid=dim.nz
    nprocy=dim.nprocz
  endif
;
;  Read variables from '*aver.in' file
;
  run_dir = stregex ('./'+datadir, '^(.*)data\/', /extract)
  variables_all = strarr(file_lines(run_dir+in_file))
  openr, lun, run_dir+in_file, /get_lun
  readf, lun, variables_all
  close, lun
  free_lun, lun
;
; Remove commented and empty elements from variables_all
;
  variables_all = strtrim (variables_all, 2)
  inds = where (stregex (variables_all, '^[a-zA-Z]', /boolean), nvar_all)
  if (nvar_all le 0) then message, "ERROR: there are no variables found."
  variables_all = variables_all[inds]
;
  if (variables[0] eq '') then begin
    variables = variables_all
    nvar = nvar_all
    ivarpos = indgen(nvar)
  endif else begin
    nvar=n_elements(variables)
    ivarpos=intarr(nvar)
;
;  Find the position of the requested variables in the list of all
;  variables.
;
    for ivar=0,nvar-1 do begin
      ivarpos_est=where(variables[ivar] eq variables_all)
      if (ivarpos_est[0] eq -1) then $
          message, 'ERROR: can not find the variable '''+variables[ivar]+'''' + $
                   ' in '+arraytostring(variables_all,/noleader)
      ivarpos[ivar]=ivarpos_est[0]
    endfor
  endelse
;  Die if attempt to plot a variable that does not exist.
  if (iplot gt nvar-1) then message, 'iplot must not be greater than nvar-1!'
;
;  Define filenames to read data from - either from a single processor or
;  from all of them.
;
  if ((ipxread eq -1) and (ipyread eq -1)) then begin
    ipxarray = indgen(nprocx*nprocy) mod nprocx
    ipyarray = (indgen(nprocx*nprocy)/nprocx) mod nprocy
    iproc = ipxarray+ipyarray*nprocx
    if (dir eq 'y') then iproc = ipxarray+ipyarray*nprocx*nprocy
    filename = datadir+'/proc'+strtrim(iproc,2)+'/'+varfile
    nxg = nxgrid
    nyg = nygrid
  endif else begin
    if ((ipxread lt 0) or (ipxread gt nprocx-1)) then begin
      print, 'ERROR: ipx is not within the proper bounds'
      print, '       ipx, nprocx', ipxread, nprocx
      stop
    endif
    if ((ipyread lt 0) or (ipyread gt nprocy-1)) then begin
      print, 'ERROR: ip'+dir+' is not within the proper bounds'
      print, '       ip'+dir+', nproc'+dir, ipyread, nprocy
      stop
    endif
    iproc = ipxread+ipyread*nprocx
    if (dir eq 'y') then iproc = ipxarray+ipyarray*nprocx*nprocy
    filename = datadir+'/proc'+strtrim(iproc,2)+'/'+varfile
    ipxarray = intarr(1)
    ipyarray = intarr(1)
    nxg = nx
    nyg = ny
  endelse
;
;  Define arrays to put data in. The user may supply the length of
;  the time dimension (nit). Otherwise nit will be read from the
;  file data/t2davg.dat.
;
  if (not keyword_set(nit)) then begin
;
;  Test if the file that stores the time and number of 2d-averages, t2davg,
;  exists.
;
    file_t2davg=datadir+'/t2davg.dat'
    if (not file_test(file_t2davg)) then begin
      print, 'ERROR: cannot find file '+ file_t2davg
      stop
    endif
    openr, lun, file_t2davg, /get_lun
    dummy_real=0.0d0
    dummy_int=0L
    readf, lun, dummy_real, dummy_int
    close, lun
    free_lun, lun
    nit=dummy_int-1
  endif
;
  if (nit gt 0) then begin
;
    if (not quiet) then begin
      if (njump eq 1) then begin
        print, 'Returning averages at '+strtrim(nit,2)+' times'
      endif else begin
        print, 'Returning averages at '+strtrim(nit,2)+' times'+ $
            ' at an interval of '+strtrim(njump,2)+' steps'
      endelse
    endif
;
    tt=fltarr(ceil(nit/double(njump)))*one
    for i=0,nvar-1 do begin
      cmd=variables[i]+'=fltarr(nxg,nyg,ceil(nit/double(njump)))*one'
      if (execute(cmd,0) ne 1) then message, 'Error defining data arrays'
    endfor
;
  endif
;
;  Define axes (default to indices if no axes are supplied).
;
  if (n_elements(xax) eq 0) then xax=findgen(nxg)
  if (n_elements(yax) eq 0) then yax=findgen(nyg)
  if (n_elements(par) ne 0) then begin
    x0=par.xyz0[0]
    x1=par.xyz1[0]
    y0=par.xyz0[1]
    y1=par.xyz1[1]
    Lx=par.Lxyz[0]
    Ly=par.Lxyz[1]
    if (dir eq 'y') then begin
      y0=par.xyz0[2]
      y1=par.xyz1[2]
      Ly=par.Lxyz[2]
    endif
  endif else begin
    x0=xax[0]
    x1=xax[nxg-1]
    y0=yax[0]
    y1=yax[nyg-1]
    Lx=xax[nxg-1]-xax[0]
    Ly=yax[nyg-1]-yax[0]
  endelse
;
;  Prepare for read.
;
  if (not quiet) then print, 'Preparing to read '+dir+'-averages ', $
      arraytostring(variables,quote="'",/noleader)
;
  num_files = n_elements(filename)
  for ip=0, num_files-1 do begin
    if (not file_test(filename[ip])) then begin
      print, 'ERROR: cannot find file '+ filename[ip]
      stop
    endif
  endfor
;
;  Method 1: Read in full data processor by processor. Does not allow plotting.
;
  if (iplot eq -1) then begin
    array_local=fltarr(nx,ny,nvar_all)*one
    array_global=fltarr(nxg,nyg,ceil(nit/double(njump)),nvar_all)*one
    t=zero
;
    for ip=0, num_files-1 do begin
      if (not quiet) then print, 'Loading chunk ', strtrim(ip+1,2), ' of ', $
          strtrim(num_files,2), ' (', filename[ip], ')'
      ipx=ipxarray[ip]
      ipy=ipyarray[ip]
      openr, lun, filename[ip], /f77, swap_endian=swap_endian, /get_lun
      it=0
      while ( not eof(lun) and ((nit eq 0) or (it lt nit)) ) do begin
;
;  Read time.
;
        readu, lun, t
        if (it eq 0) then t0=t
;
;  Read full data, close and free lun after endwhile.
;
        if ( (t ge tmin) and (it mod njump eq 0) ) then begin
          readu, lun, array_local
          array_global[ipx*nx:(ipx+1)*nx-1,ipy*ny:(ipy+1)*ny-1,it/njump,*]=array_local
          tt[it/njump]=t
        endif else begin
          dummy=zero
          readu, lun, dummy
        endelse
        it=it+1
      endwhile
      close, lun
      free_lun, lun
;
    endfor
;
;  Diagnostics.
;
    if (not quiet) then begin
      for it=0,ceil(nit/double(njump))-1 do begin
        if (it mod it1 eq 0) then begin
          if (it eq 0 ) then begin
            print, '  ------ it -------- t ---------- var ----- min(var) ------- max(var) ------'
          endif
          for ivar=0,nvar-1 do begin
            print, it, tt[it], variables[ivar], $
                min(array_global[*,*,it,ivarpos[ivar]]), $
                max(array_global[*,*,it,ivarpos[ivar]]), $
                format='(i11,e17.7,A12,2e17.7)'
          endfor
        endif
      endfor
    endif
;
;  Possible to shift data in the x-direction.
;
    if (xshift ne 0) then begin
      for it=0, nit/njump-1 do begin
        pc_read_aver_shift_plane, nvar_all, array_global[*,*,it,*], xshift, par, timefix, ts, t, t0
      endfor
    endif
;
;  Split read data into named arrays.
;
    for ivar=0,nvar-1 do begin
      cmd=variables[ivar]+'=array_global[*,*,*,ivarpos[ivar]]'
      if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
    endfor
;
  endif else begin
;
;  Method 2: Read in data slice by slice and plot the results.
;
    if (png) then begin
      set_plot, 'z'
      device, set_resolution=[zoom*nx,zoom*nyg]
      print, 'Deleting old png files in the directory ', imgdir
      spawn, '\rm -f '+imgdir+'/img_*.png'
    endif
;
;  Open all ipz=0 processor directories.
;
    luns = intarr(num_files)
    for ip=0, num_files-1 do begin
      openr, lun, filename[ip], /f77, swap_endian=swap_endian, /get_lun
      luns[ip] = lun
    endfor
;
;  Variables to put single time snapshot in.
;
    array=fltarr(nxg,nyg,nvar_all)*one
    array_local=fltarr(nx,ny,nvar_all)*one
    tt_local=fltarr(num_files)*one
;
;  Read 2D-averages and put in arrays if requested.
;
    it=0
    itimg=0
    lwindow_opened=0
    while ( not eof(luns[0]) and ((nit eq 0) or (it lt nit)) ) do begin
;
;  Read time.
;
      for ip=0, num_files-1 do begin
        tt_tmp=zero
        readu, luns[ip], tt_tmp
        tt_local[ip]=tt_tmp
        if (ip ge 1) then begin
          if (max(tt_local[ip]-tt_local[0:ip-1])) then begin
            print, 'Inconsistent times between processors!'
            print, tt_local
            stop
          endif
        endif
      endfor
      t=tt_local[0]
      if (it eq 0) then t0=t
;
;  Read data one slice at a time.
;
      if ( (t ge tmin) and (it mod njump eq 0) ) then begin
        for ip=0, num_files-1 do begin
          readu, luns[ip], array_local
          ipx=ipxarray[ip]
          ipy=ipyarray[ip]
          array[ipx*nx:(ipx+1)*nx-1,ipy*ny:(ipy+1)*ny-1,*]=array_local
        endfor
;
;  Shift plane in the radial direction.
;
        if (xshift ne 0) then $
            pc_read_aver_shift_plane, nvar_all, array, xshift, par, timefix, ts, t
;
;  Interpolate position of stalked particles to time of snapshot.
;
        if (stalk) then begin
          ist=where(abs(pst.t-t) eq min(abs(pst.t-t)))
          ist=ist[0]
          if (max(tag_names(pst) eq 'APS') eq 1) then begin
            if (n_elements(pst.ipar) eq 1) then begin
              if (pst.aps[ist] eq 0.0) then nstalk=0 else nstalk=-1
            endif else begin
              ipar=where(pst.aps[*,ist] ne 0.0)
              if (ipar[0] eq -1) then nstalk=0 else nstalk=n_elements(ipar)
            endelse
          endif else begin
            ipar=indgen(nstalk)
          endelse
          if (nstalk eq -1) then begin
            xxstalk=pst.xp[ist]
            if (dir eq 'z') then yystalk=pst.yp[ist]
            if (dir eq 'y') then yystalk=pst.zp[ist]
          endif else if (nstalk ne 0) then begin
            xxstalk=pst.xp[ipar,ist]
            if (dir eq 'z') then yystalk=pst.yp[ipar,ist]
            if (dir eq 'y') then yystalk=pst.zp[ipar,ist]
          endif
        endif
;
;  Plot requested variable (plotting is turned off by default).
;
        pc_read_aver_plot_plane, array_plot=array[*,*,ivarpos[iplot]], nxg=nxg, nyg=nyg, $
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
            x0=x0, x1=x1, y0=y0, y1=y1, Lx=Lx, Ly=Ly, time=t, itimg=itimg
;
;  Diagnostics.
;
        if ( (not quiet) and (it mod it1 eq 0) ) then begin
          if (it eq 0 ) then $
              print, '  ------ it -------- t ---------- var ----- min(var) ------- max(var) ------'
          for ivar=0,nvar-1 do begin
              print, it, t, variables[ivar], $
                  min(array[*,*,ivarpos[ivar]]), max(array[*,*,ivarpos[ivar]]), $
                  format='(i11,e17.7,A12,2e17.7)'
          endfor
        endif
;
;  Split read data into named arrays.
;
        if (it le nit-1) then begin
          tt[it/njump]=t
          for ivar=0,nvar-1 do begin
            cmd=variables[ivar]+'[*,*,it/njump]=array[*,*,ivarpos[ivar]]'
            if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
          endfor
        endif
      endif else begin
        for ip=0, num_files-1 do begin
          dummy=zero
          readu, luns[ip], dummy
        endfor
      endelse
;
      it=it+1
;
    endwhile
;
;  Close files and free luns.
;
    for ip=0, num_files-1 do begin
      close, luns[ip]
      free_lun, luns[ip]
    endfor
  endelse
;
;  Put data in structure. One must set nit>0 to specify how many lines
;  of data that should be saved.
;
  if (nit ne 0) then begin
    makeobject="object = create_struct(name=objectname,['t'," + $
        arraytostring(variables,quote="'",/noleader) + "],"+"tt,"+$
        arraytostring(variables,/noleader) + ")"
;
    if (execute(makeobject) ne 1) then begin
      message, 'ERROR evaluating variables: ' + makeobject, /info
      undefine,object
    endif
  endif
;
end
