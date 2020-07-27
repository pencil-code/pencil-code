;
; $Id$
;
;  Read 2D-averages from file.
;
;  1-mar-17/MR: added new keyword parameter write: allows to write the averages,
;               just read in, to destination write; if write eq '.' or write eq './'
;               the averages in the cirrent working directory will be overwritten - 
;               meaningful to reduce their size, when reading was for restricted time
;               interval/coarser grained/with less variables. In the latter case the *aver.in
;               has to be adjusted by hand.
;               If write is another working directory, the individual processor subdirectories
;               of data must exist.
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
    readpar=readpar, readgrid=readgrid, debug=debug, quiet=quiet, wait=wait, $
    write=write, single=single, dim=dim, grid=grid
;
COMPILE_OPT IDL2,HIDDEN
;
  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
  if ((dir ne 'y') and (dir ne 'z')) then message, 'ERROR: averaging direction "'+strtrim(dir,2)+'" unknown.'
;
;  Default values.
;
  datadir = pc_get_datadir(datadir)
  default, varfile, dir+'averages.dat'
  default, nit, 0
  default, ipxread, -1
  default, ipyread, -1
  default, iplot, -1
  default, zoom, 1
  default, min, 1.0
  default, max, 0.0
  default, tmin, 0.0
  default, njump, 1
  default, ps, 0
  default, png, 0
  default, noerase, 0
  default, imgdir, '.'
  default, xsize, 10.0
  default, ysize, 10.0
  default, title, 'notitle'
  default, t_title, 0
  default, t_scale, 1.0
  default, t_zero, 0.0
  default, stalk, 0
  default, interp, 0
  default, ceiling, 0.0
  default, subbox, 0
  default, subcen, -1
  default, subpos, [0.7,0.7,0.9,0.9]
  default, submin, 1.
  default, submax, 0.
  default, sublog, 0
  default, rsubbox, 5.0
  default, subcolor, 255
  default, tsubbox, 0.0
  default, thick, 1.0
  default, charsize, 1.0
  default, loge, 0
  default, log10, 0
  default, fillwindow, 0
  default, tformat, '(g8.2)'
  default, swap_endian, 0
  default, it1, 10
  default, variables, ''
  default, colorbar, 0
  default, bartitle, ''
  default, xshift, 0
  default, timefix, 0
  default, readpar, 1
  default, readgrid, 0
  default, debug, 0
  default, quiet, 0
  default, wait, 0
  default, yinyang, 0
  default, xtitle, 'x'
  default, ytitle, dir eq 'y' ? 'z' : 'y'
;
  if (size (grid, /type) eq 0) then pc_read_grid, obj=grid, dim=dim, datadir=datadir, /quiet
;
  ; load HDF5 averages, if available
  if (file_test (datadir+'/averages/'+dir+'.h5')) then begin
    message, "pc_read_2d_aver: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    last = pc_read ('last', filename=dir+'.h5', datadir=datadir+'/averages/')
    groups = str (lindgen (last + 1))
    times = pc_read (groups[0]+'/time')
    for it=1, last do times=[times,pc_read(groups[it]+'/time')]
    in = (where (times ge tmin))[0:*:njump]
    num = n_elements (in)
    if (num le 0) then message, 'pc_read_2d_aver: ERROR: "'+h5_file+'" no data available after given "tmin"!'
    groups = groups[in]
    times = times[in]
    if (size (vars, /type) ne 7) then vars = h5_content (groups[0])
    found = where (strlowcase (vars) ne 'time', numb)
    vars = vars[found]
    object = { t:times, last:last, pos:1+long (groups), num_quantities:numb, labels:vars }
    object = create_struct (object, 'x', grid.x[dim.l1:dim.l2])
;
; Note: dir is the direction of averaging, hence for dir='z', the second variable is y,
;       for dir='y', the second variable is z.
;
    if (dir eq 'z') then begin
      object = create_struct (object, 'y', grid.y[dim.m1:dim.m2])
    end else if (dir eq 'y') then begin
      object = create_struct (object, 'z', grid.z[dim.n1:dim.n2])
    end
    for pos = 0, numb-1 do begin
       variable=pc_read (groups[0]+'/'+vars[pos])
       variable=spread(variable,2,num)
       for it=1, num-1 do variable[*,*,it]=pc_read (groups[it]+'/'+vars[pos])
      object = create_struct (object, vars[pos], variable)
    end
    h5_close_file
    return
  end

;
  if strpos(varfile,'0.dat') ge 0 then begin
    default, in_file, dir+'aver0.in'
    file_t2davg=datadir+'/t2davg0.dat'
  endif else begin
    default, in_file, strtrim(dir,2)+'aver.in'
    file_t2davg=datadir+'/t2davg.dat'
  endelse
;
;  Read additional information.
;
  pc_read_dim, obj=dim, datadir=datadir, /quiet
  pc_read_dim, obj=locdim, datadir=datadir+'/proc0', /quiet
  if (stalk) then begin
    pc_read_pstalk, obj=pst, datadir=datadir, /quiet
    default, nstalk, n_elements(pst.ipar)
  endif
;
;  Need to know box size for proper axes and for shifting the shearing
;  frame.
;
  if (readpar or (xshift ne 0)) then begin
    pc_read_param, obj=par, datadir=datadir, /quiet
;
; We know from param whether we have a Yin-Yang grid.
;
    if tag_exists(par,'LYINYANG') then $
      if (par.lyinyang) then yinyang=1
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
  nxgrid=dim.nx
  nprocx=dim.nprocx
  nprocy=dim.nprocy
  nprocz=dim.nprocz
  if (dir eq 'y') then begin
    ny=locdim.nz
    nygrid=dim.nz
  endif else begin
    ny=locdim.ny
    nygrid=dim.ny
  endelse

  nycap=0

  if (readgrid or (xshift ne 0) or yinyang) then begin
    pc_read_grid, obj=grid, /trim, datadir=datadir, /quiet
    xax=grid.x
    if (dir eq 'y') then $
      yax=grid.z $
    else $
      yax=grid.y

    if yinyang then begin
      nycap=floor(grid.y[0]/grid.dy)-1
      yintl=reverse(yax[0]-grid.dy*(indgen(nycap)+1))       ; only valid fo requidistant grid
      yinto=yax[nygrid-1]+grid.dy*(indgen(nycap)+1)
      yax=[yintl,yax,yinto]
    endif

  endif
;
;  Read variables from '*aver.in' file
;
;  run_dir = stregex ('./'+datadir, '^(.*)data\/', /extract)
  if (datadir eq 'data') then begin
    rd='./'
   endif else begin
    rd=strsplit(datadir,'data',/extract,/regex)
  endelse
  run_dir=rd[0]
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
    nfound=0
    for ivar=0,nvar-1 do begin
      ivarpos_est=where(variables[ivar] eq variables_all)
      if (ivarpos_est[0] eq -1) then begin
        message, 'ERROR: cannot find the variable '''+variables[ivar]+'''' + $
                 ' in '+arraytostring(variables_all,quote="'",/noleader)+'.', /continue
        variables[ivar]=''
      endif else begin
        ivarpos[nfound]=ivarpos_est[0]
        nfound++
      endelse
    endfor
    if nfound eq 0 then $
      message, 'ERROR: No valid variable names in request.' $
    else begin
      nvar=nfound
      ivarpos=ivarpos[0:nfound-1]
      variables = variables[where(variables ne '')]
    endelse
  endelse
;
;  Die if attempt to plot a variable that does not exist.
;
  if (iplot gt nvar-1) then message, 'iplot must not be greater than nvar-1!'
;
;  Define filenames to read data from - either from a single processor or
;  from all of them.
;
  if ((ipxread eq -1) and (ipyread eq -1)) then begin
    nxg = nxgrid
    if (dir eq 'y') then begin
      ipxarray = indgen(nprocx*nprocz) mod nprocx
      ipyarray = (indgen(nprocx*nprocz)/nprocx) mod nprocz
      iproc = ipxarray+ipyarray*nprocx*nprocy
      nyg = nygrid
    end else begin
      ipxarray = indgen(nprocx*nprocy) mod nprocx
      ipyarray = (indgen(nprocx*nprocy)/nprocx) mod nprocy
      iproc = ipxarray+ipyarray*nprocx
      if yinyang then begin
        ncpus=nprocx*nprocy*nprocz
        iproc=[iproc,ncpus+indgen(nprocx),2*ncpus-nprocx+indgen(nprocx)]
        ipxarray=[ipxarray,indgen(nprocx),indgen(nprocx)]
        nyg=nygrid+2*nycap
      endif else $
        nyg=nygrid
    end
    filename = datadir+'/proc'+strtrim(iproc,2)+'/'+varfile
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
    if (dir eq 'y') then begin
      iproc = ipxread+ipyread*nprocx*nprocy
      nyg = nz
    endif else begin
      if yinyang then begin
        print, 'Reading from individual procs not implemented for Yin-Yang grid!'
        stop
      endif
      iproc = ipxread+ipyread*nprocx
      nyg = ny
    endelse
    filename = datadir+'/proc'+strtrim(iproc,2)+'/'+varfile
    ipxarray = intarr(1)
    ipyarray = intarr(1)
    nxg = nx
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
    if nit le 0 then message, 'Error: Supposedly no averages existent'
  endif
;
  if (nit gt 0) then begin
;
    nret=ceil(nit/double(njump))
    relchar = tmin gt 0 ? '<=' : ''
    if (not quiet) then $
      if (njump eq 1) then $
        print, 'Returning averages at '+relchar+strtrim(nit,2)+' times.' $
      else $
        print, 'Returning averages at '+relchar+strtrim(nret,2)+' times out of '+strtrim(nit,2)+ $
               ' at an interval of '+strtrim(njump,2)+' steps.'
;
    tt=fltarr(nret)*one
;
  endif
;
;  Define axes (default to indices if no axes are supplied).
;
  if (n_elements(par) gt 1) then begin
    x0=par.xyz0[0]
    x1=par.xyz1[0]
    Lx=par.Lxyz[0]
    if (dir eq 'y') then begin
      y0=par.xyz0[2]
      y1=par.xyz1[2]
      Ly=par.Lxyz[2]
    endif else begin
      y0=par.xyz0[1]
      y1=par.xyz1[1]
      Ly=par.Lxyz[1]
    endelse
    if (n_elements(xax) eq 0) then xax=findgen(nxg)/(nxg-1.)*Lx+x0
    if (n_elements(yax) eq 0) then yax=findgen(nyg)/(nyg-1.)*Ly+y0
  endif else begin
    if (n_elements(xax) eq 0) then xax=findgen(nxg)
    if (n_elements(yax) eq 0) then yax=findgen(nyg)
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
    if (not file_test(filename[ip])) then $
      message, 'ERROR: cannot find file '+ filename[ip]
  endfor
;
;  Method 1: Read in full data processor by processor. Does not allow plotting.
;
  if (iplot lt 0) then begin
;
    array_local=fltarr(nx,ny,nvar_all)*one
    if keyword_set(single) then $
      array_global=fltarr(nxg,nyg,nret,nvar) $
    else $
      array_global=fltarr(nxg,nyg,nret,nvar)*one
    type_as_on_disk=not keyword_set(single) or size(one,/type) eq 4

    t=zero
    if njump gt 1 then incr=njump*(double(n_elements(array_local))*data_bytes+8L + data_bytes+8L)
;
    for ip=0, num_files-1 do begin
;
      if (not quiet) then print, 'Loading chunk ', strtrim(ip+1,2), ' of ', $
          strtrim(num_files,2), ' (', filename[ip], ')'

      if yinyang and ip ge num_files-2*nprocx then begin
        nyl=nycap
        if ip eq num_files-2*nprocx then $
          array_local=fltarr(nx,nyl,nvar_all)*one
        if ip lt num_files-nprocx then $
          iya=nygrid+nycap $
        else $
          iya=0
      endif else begin
        nyl=ny
        iya=ipyarray[ip]*ny+nycap
      endelse

      ipx=ipxarray[ip]
      iye=iya+nyl-1
      openr, lun, filename[ip], /f77, swap_endian=swap_endian, /get_lun
      it=0 & nread=0 & filepos=-1.d0

      while (not eof(lun) and ((nit eq 0) or (it lt nit))) do begin
;
;  Read time.
;
        if filepos lt 0. then point_lun, -lun, filepos
        readu, lun, t

        if (it eq 0) then t0=t
;
;  Read full data, close and free lun after endwhile.
;
        if t ge tmin then begin
          readu, lun, array_local
          array_global[ipx*nx:(ipx+1)*nx-1,iya:iye,nread,*]=array_local[*,*,ivarpos]
          tt[nread]=t
          nread++
          if njump gt 1 then begin
            filepos+=incr
            point_lun, lun, filepos
          endif
          it=it+njump
        endif else begin
          dummy=zero
          readu, lun, dummy
          it=it+1
          filepos=-1.d0
        endelse
      endwhile
      close, lun
      free_lun, lun
;
    endfor
;
    if nread gt 0 then begin
      tt=tt[0:nread-1]
      array_global=array_global[*,*,0:nread-1,*]
      if (not quiet) then $
        print, 'Returned averages at '+strtrim(nread,2)+' times.'
    endif else $
      message, 'No averages found for t>='+strtrim(tmin,2)+'.'

    if nread gt 0 and is_defined(write) then begin

      if write ne '' then begin
; 
;  If a path for writing the averages is given.
;
        if yinyang then $
          print, 'Warning: writing of averages not implemented for Yin-Yang grid!' $
        else begin

          warn=0
          if write eq '.' or write eq './' then begin
            warn=1 & write='.'
          endif else begin
            cd, current=wdr
            len=strlen(strtrim(write,2))
            if strmid(strtrim(write,2),len,1) eq '/' then len-=1
            if wdr eq strmid(write,0,len) then warn=1
          endelse

          if warn then begin
            yn=''
            read, prompt='Do you really want to overwrite the averages (yes/no)?', yn 
            if yn eq 'yes' and not type_as_on_disk then $
              read, prompt='Data was read in in single precision, so you loose precision in overwriting. Continue (yes/no)?', yn
          endif else $
            yn='yes'

          if yn eq 'yes' then begin
              
            for ip=0, num_files-1 do begin
              ipx=ipxarray[ip]
              iya=ipyarray[ip]*ny
              iye=iya+ny-1
              openw, lun, write+'/'+filename[ip], /f77, swap_endian=swap_endian, /get_lun
              for itt=0,nread-1 do begin
                writeu, lun, tt[itt], tt[itt]
                if type_as_on_disk then $
                  writeu, lun, array_global[ipx*nx:(ipx+1)*nx-1,iya:iye,itt,*] $
                else $
                  writeu, lun, double(array_global[ipx*nx:(ipx+1)*nx-1,iya:iye,itt,*])
              endfor
              close, lun
              free_lun, lun
            endfor
          endif
        endelse
      endif
    endif
;
;  Diagnostics.
;
    if (not quiet) then begin
      varnamsiz=(max(strlen(strtrim(variables,2)))+7)/2
      filler=arraytostring(replicate('-',varnamsiz),list='')

      for it=0,nread-1 do begin
        if (it mod it1 eq 0) then begin
          if (it eq 0 ) then $
            print, '  ------ it -------- t ------'+filler+' var '+filler+'- min(var) ------- max(var) ------'
          for ivar=0,nvar-1 do $
            print, it, tt[it], variables[ivar], $
                min(array_global[*,*,it,ivar]), $
                max(array_global[*,*,it,ivar]), $
                format='(i11,e17.7,A12,2e17.7)'
        endif
      endfor
    endif
;
;  Possible to shift data in the x-direction.
;
    if (xshift ne 0) then $
      for it=0, nread-1 do $
        pc_read_aver_shift_plane, nvar, array_global[*,*,it,*], xshift, par, timefix, ts, t, t0
;
;  Split read data into named arrays.
;
    for ivar=0,nvar-1 do begin
      cmd=variables[ivar]+'=array_global[*,*,*,ivar]'
      if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
    endfor
;
  endif else begin
;
;  Method 2: Read in data slice by slice and plot the results.
;
    for i=0,nvar-1 do begin
      cmd=variables[i]+'=fltarr(nxg,nyg,ceil(nit/double(njump)))*one'
      if (execute(cmd,0) ne 1) then message, 'Error defining data arrays'
    endfor

    if (png) then begin
      set_plot, 'z'
      device, set_resolution=[zoom*nxg,zoom*nyg]
      print, 'Deleting old png files in the directory ', imgdir
      spawn, 'rm -f '+imgdir+'/img_*.png'
    endif
;
;  Open all ip[yz]=0 processor directories.
;
    luns = intarr(num_files)
    for ip=0, num_files-1 do begin
      openr, lun, filename[ip], /f77, swap_endian=swap_endian, /get_lun
      luns[ip] = lun
    endfor
;
;  Variables to put single time snapshot in.
;
    array=fltarr(nxg,nyg,nvar)*one
    array_local=fltarr(nx,ny,nvar_all)*one
    tt_local=fltarr(num_files)*one
;
;  Read 2D-averages and put in arrays if requested.
;
    it=0 & nread=0 & dummy=zero
    itimg=0
    lwindow_opened=0
    while (not eof(luns[0]) and ((nit eq 0) or (it lt nit))) do begin
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
          array[ipx*nx:(ipx+1)*nx-1,ipy*ny:(ipy+1)*ny-1,*]=array_local[*,*,ivarpos]
        endfor
;
;  Shift plane in the x direction.
;
        if (xshift ne 0) then $
            pc_read_aver_shift_plane, nvar, array, xshift, par, timefix, ts, t, t0
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
          endif else $
            ipar=indgen(nstalk)

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
        pc_read_aver_plot_plane, array_plot=array[*,*,iplot], nxg=nxg, nyg=nyg, $
            min=min, max=max, zoom=zoom, xax=xax, yax=yax, $
            xtitle=xtitle, ytitle=ytitle, title=title eq 'notitle' ? strtrim(variables[iplot],2) : title, $
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
            x0=x0, x1=x1, y0=y0, y1=y1, Lx=Lx, Ly=Ly, time=t, itimg=itimg, wait=wait
;
;  Diagnostics.
;
        if ( (not quiet) and (it mod it1 eq 0) ) then begin
          if (it eq 0 ) then $
              print, '  ------ it -------- t ---------- var ----- min(var) ------- max(var) ------'
          for ivar=0,nvar-1 do $
              print, it, t, variables[ivar], $
                  min(array[*,*,ivar]), max(array[*,*,ivar]), $
                  format='(i11,e17.7,A12,2e17.7)'
        endif
;
;  Split read data into named arrays.
;
        tt[nread]=t
        for ivar=0,nvar-1 do begin
          cmd=variables[ivar]+'[*,*,nread]=array[*,*,ivar]'
          if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
        endfor
        nread++
      endif else $
        for ip=0, num_files-1 do readu, luns[ip], dummy
;
      it=it+1
;
    endwhile

    if (nread gt 0) then begin
;
;  Trimming of time and variables arrays.
;
      tt=tt[0:nread-1]
      for ivar=0,nvar-1 do begin
        cmd=variables[ivar]+'='+variables[ivar]+'[*,*,0:nread-1]'
        if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
      endfor
      if (not quiet) then $
        print, 'Returned averages at '+strtrim(nread,2)+' times.'
    endif else $
      message, 'No averages found for t>='+strtrim(tmin,2)+'.'
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
    if (yinyang) then begin 
      coornames="'r','theta',"
      coors="xax,yax,"
    endif else begin
      coornames="" & coors=""
    endelse

    makeobject="object = create_struct(name=objectname,['t',"+coornames + $
        arraytostring(variables,quote="'",/noleader) + "],"+"tt,"+coors+$
        arraytostring(variables,/noleader) + ")"
;
    if (execute(makeobject) ne 1) then begin
      message, 'ERROR evaluating variables: ' + makeobject, /info
      undefine,object
    endif
  endif
;
end
