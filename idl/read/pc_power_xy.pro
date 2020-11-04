FUNCTION gen_exinds, vec

;  7-jan-11/MR: coded

; determines the indices of the local extrema in the vector vec
  
  if vec(1) ge vec(0) then $
    grows=-1 $
  else $
    grows=1
    
  imax=n_elements(vec)-1
  imin=0
  s0=1
  
  while imin lt imax do begin

    ind = (where(grows*(vec(imin:imax-1)-vec(imin+1:imax)) lt 0.))(0)
    
    if ind ne -1 then begin
    
      ind = ind+imin
      if s0 then begin
  	exinds = [ind] 
  	s0 = 0
      endif else $
  	exinds = [exinds, ind]
  	
      imin = ind
      grows = -grows
      
    endif else $
      imin = imax
    
  endwhile
	
  if n_elements(exinds) eq 0 then exinds = -1
  
  return, exinds
  
END
;******************************************************************************************
PRO plot_segs, x, y, seginds, styles=styles, colors=colors, overplot=overplot, _extra=extra

;  7-jan-11/MR: coded

;  plots y vs x with segemented coloring/linestyling according to the indices of the segment boundaries seginds

  if not keyword_set(colors) then $
    colors = [!p.color,!p.color] $
  else if n_elements(colors) eq 1 then begin
    if colors(0) eq 250 then $
      colors = [colors, !p.color] $
    else $
      colors = [colors, 250]
  endif

  if not keyword_set(styles) then begin
    if colors(0) eq colors(1) then $
      styles = [0,3] $
    else $
      styles = [0,0]
  endif else if n_elements(styles) eq 1 then begin
    if styles(0) eq 3 then $
      styles = [styles, 0] $
    else $
      styles = [styles, 3]
  endif
      
  np = n_elements(x) < n_elements(y)
  
  if np le 1 then return
  
  if not keyword_set(overplot) then $
    plot, x(0:np-1), y(0:np-1), /nodata, _extra=extra
  
  ia = 0
  toggle = 0
  
  for i=0,n_elements(seginds)-1 do begin
  
    ie = seginds(i) > 0
    oplot, x(ia:ie), y(ia:ie), lines=styles(toggle), color=colors(toggle)
    ia = ie
    toggle = 1-toggle
    
  endfor
  
  if ia lt np-1 then oplot, x(ia:np-1), y(ia:np-1), lines=styles(toggle), color=colors(toggle)

END
;******************************************************************************************
PRO alloc_specs, spec, lint_shell, lint_z, lcomplex, fullspec=fullspec

common pars,  nx, ny, nz, nk, ncomp, nt, Lx, Ly, Lz

  if lint_shell then begin
    if lint_z then $
      dims='nk' $
    else $
      dims='nk,nz'      
  endif else if lint_z then $
    dims='nx,ny' $
  else $
    dims='nx,ny,nz'

  if lcomplex then $
    cmd = 'complex' $
  else $
    cmd = 'flt'

  cmd = 'spec='+cmd+'arr('+dims

  ierr = execute(cmd+')') 
  if ~ierr then begin
    print, 'alloc_specs: Allocation '+cmd+') failed!!!'
    stop
  endif
  
  if keyword_set(fullspec) then begin
    ierr = execute('full'+cmd+',ncomp,nt)')
    if ~ierr then begin
      print, 'alloc_specs: Allocation '+'full'+cmd+',nt) failed!!!'
      stop
    endif
  endif
END
;******************************************************************************************
FUNCTION find_no_wavenos, headline, text, symbol, no

  opos = strpos(headline,'(') & cpos = strpos(headline,')')
  len = cpos-opos-1
  
  if len le 0 then begin
    print, 'Warning - header corrupt: Number of '+strtrim(text,2)+' wavenumbers missing -'
    print, '         will use '+strtrim(symbol,2)+'='+strtrim(string(no),2)+' which is perhaps by one too large.'
    print, '         Correct data file by hand (or '+strtrim(symbol,2)+' in code if necessary).'
  endif else begin
    nol=no
    on_ioerror, readerr
    reads, strmid(headline,opos+1,len), nol
    no=nol
    return, cpos+1

readerr:
    print, ' Error when reading '+strtrim(symbol,2)+': '+!ERROR_STATE.MSG
    on_ioerror, NULL
    stop
    return, -1
  endelse

END
;******************************************************************************************
FUNCTION read_firstpass, file, lint_shell, lint_z, lcomplex, extr, startpos, fmt=fmt, zpos=zpos

common pars,  nx, ny, nz, nk, ncomp, nt, Lx, Ly, Lz
common wavenrs, kxs, kys, kshell
;
;  Looping through all the data to get number of spectral snapshots
;
close,2

openr,2, file, error=err
if err ne 0 then return, 0

  headline = ''
  readf, 2, headline
  headline = strlowcase(headline)
  witheader=strpos(headline,'spectrum') ne -1

  if witheader then begin
   
    warn = 'Warning: File header of '+strtrim(file,2)
    if (strpos(headline,'shell-integrated') ne -1) ne lint_shell then begin
    
      if lint_shell then $
        print, warn+' says "not shell-integrated" - assume that' $
      else $
        print, warn+' says "shell-integrated" - assume that'
        
      lint_shell = ~lint_shell
      
    endif
    
    if (strpos(headline,'z-integrated') ne -1) ne lint_z then begin

      if lint_z then $
        print, warn+' says "not z-integrated" - assume that' $
      else $
        print, warn+' says "z-integrated" - assume that'
        
      lint_z = ~lint_z
      
    endif

    lcomplex = strpos(headline,'complex') ne -1 

    if strpos(headline,'componentwise') ne -1 then begin
      ia = strpos(headline,'(')+1 & ie = strpos(headline,')')-1
      ncomp = fix(strmid(headline,ia,ie-ia+1))      
    endif else $
      ncomp = 1
  
    readf, 2, headline

    if strpos(strlowcase(headline),'wavenumbers') eq -1 then begin
      print, warn+' corrupt (keyword "wavenumbers" missing!' 
      lcorr=1
    endif else $
      lcorr=0
   
    if (lint_shell) then begin

      if not lcorr then $
       ind=find_no_wavenos(headline, 'shell', 'nk', nk)
 
      kshell = fltarr(nk) 
      readf, 2, kshell 

    endif else begin

      if not lcorr then begin
        
        ind=find_no_wavenos(headline, 'x', 'nx', nx)
        if ind ge 0 then ind=find_no_wavenos(strmid(headline,ind), 'y', 'ny', ny)

      endif

      kxs = fltarr(nx) & kys = fltarr(ny)

      readf, 2, kxs, kys 
      kxs = shift(kxs,nx/2)
      kys = shift(kys,ny/2)  

    endelse

    if not lint_z then begin
      point_lun, -2, pos     
      readf, 2, headline
      if strpos(strlowcase(headline),'positions') eq -1 then begin
;
; no z-positions in data file -> spectra given for the whole z grid
;
        point_lun, 2, pos     
        pc_read_grid,obj=grid,/trim,/quiet
        zpos=grid.z
        nz = n_elements(zpos)
        lzpos_exist=0
;
      endif else begin
        ia = strpos(headline,'(')+1 & ie = strpos(headline,')')-1
        nz = fix(strmid(headline,ia,ie-ia+1))
          
        if nz le 0 then begin
          print, warn+' corrupt! -- no positive number of z positions given!'
          stop
        endif else begin
          zpos=fltarr(nz)
          readf, 2, zpos
        endelse  
        lzpos_exist=1
      endelse
    endif
  
  endif else begin

    kxs=(findgen(nx)-nx/2)*2*!pi/Lx                ;,(nx+1)/2)*2*!pi/Lx
    kys=(findgen(ny)-ny/2)*2*!pi/Ly                ;,(ny+1)/2)*2*!pi/Ly
     
    if lint_shell then begin                       ; valid for old style files

      kshell = fltarr(nk)
       
      for ix=0,nx-1 do $
        for iy=0,ny-1 do begin
          ik = round( sqrt( kxs(ix)^2+kys(iy)^2 )/(2*!pi/Lx) )
	  kshell(ik)=ik*2*!pi/Lx
        endfor

    endif

    ncomp=1 
    pos=0L
    point_lun, 2, pos        ; set file position to top of file (i.e., at first time stamp)
    
  endelse

  alloc_specs, spectrum1, lint_shell, lint_z, lcomplex
  alloc_specs, spectrum1y, lint_shell, lint_z, lcomplex
  alloc_specs, spectrum1z, lint_shell, lint_z, lcomplex
   
  if lcomplex then begin  
    globalmin=1e12 +complexarr(ncomp)
    globalmax=1e-30+complexarr(ncomp)
  endif else begin
    globalmin=1e12 +fltarr(ncomp)
    globalmax=1e-30+fltarr(ncomp)
  endelse

  nt=0L & time=0. & s0=1 & s0form=1 & nseg=0 
  
  while not eof(2) do begin

    if s0 then begin
      if witheader then begin
        
        point_lun, -2, pos
        readf,2,headline
      
        if strpos(strlowcase(headline),'spectrum') eq -1 then $
          point_lun, 2, pos $
        else begin
          readf,2,headline
          readf,2, kxs, kys
          if lint_shell then begin       
            readf, 2, headline
            readf, 2, kshell 
          endif
          if not lint_z and lzpos_exist then begin
            readf,2,headline
            readf,2, zpos
          endif
        endelse

      nz = abs(nz)
      endif

    endif

    if witheader then begin
      on_ioerror, newheader
      point_lun, -2, timepos
    endif
    readf,2,time 
    if witheader then begin
      on_ioerror,  NULL
      if s0 then begin
        nseg +=1 
        if nseg eq 1 then startpos = timepos else startpos = [startpos,timepos]
        s0=0
      endif
    endif

    if s0form then begin 
      if lcomplex then begin
;
;  derive format for reading from data line (only needed for complex data)
;
        dataline = ''
        point_lun, -2, pos 
        readf, 2, dataline
        point_lun, 2, pos 

        inds = strsplit(dataline,' ',length=lens) 
        num = n_elements(inds) & ends = inds+lens 
        lens = ends-[0,ends(0:num-2)]
        fmt = '(('
        for ii=0,num-1,2 do $
          fmt += '(e'+strtrim(string(lens(ii)),2)+''+'.1,e'+strtrim(string(lens(ii+1)),2)+'.1),'
        fmt = strmid(fmt,0,strlen(fmt)-1)+'))'

      endif else $
        fmt=''

      s0form = 0
    endif

    for i=0,ncomp-1 do begin

      if fmt eq '' then $
        readf,2,spectrum1 $
      else $
        readf,2,spectrum1,format=fmt 

      ;readf,2,spectrum1y,format=fmt 
      ;readf,2,spectrum1z,format=fmt 

      globalmax(i)=max(spectrum1) > globalmax(i)
      globalmin(i)=min(spectrum1) < globalmin(i)
    endfor

    nt=nt+1L
    continue
newheader:
    if eof(2) then exit $
    else begin
      s0=1
      point_lun, 2, timepos
    endelse
  
  endwhile

  close, 2

; end first reading

  extr = [globalmin,globalmax]
  
  return, nt

END
;******************************************************************************************
PRO pc_power_xy,var1,var2,last,w,v1=v1,v2=v2,all=all,wait=wait,k=k,spec1=spec1, $
          spec2=spec2,i=i,tt=tt,noplot=noplot,tmin=tmin,tmax=tmax, $
          tot=tot,lin=lin,png=png,yrange=yrange,norm=norm,helicity2=helicity2, $
          compensate1=compensate1,compensate2=compensate2,datatopdir=datatopdir, $ 
	  lint_shell=lint_shell, lint_z=lint_z, print=prnt, obj=obj
;
;  $Id$
;
;  This routine reads in the power spectra generated during the run
;  (provided dspec is set to a time interval small enough to produce
;  enough spectra.) 
;  By default, spec1 is kinetic energy and spec2 magnetic energy.
;  The routine plots the spectra for the last time saved in the file.
;  All times are stored in the arrays spec1 and spec2.
;  The index of the first time is 1, and of the last time it is i-2.
;
;  var1  : Use v1 instead (Kept to be backward compatible)
;  var2  : Use v2 instead (Kept to be backward compatible)
;  last  : Use all instead (Kept to be backward compatible)
;  w     : Use wait instead (Kept to be backward compatible)
;
;  v1    : First variable to be plotted (Ex: 'u')
;  v2    : Second variable to be plotted (Ex: 'b') 
;  all   : Plot all snapshots if /all is set, otherwise only last snapshot
;  wait  : Time to wait between each snapshot (if /all is set) 
;  k     : Returns the wavenumber vector
;  spec1 : Returns all the spectral snapshots of the first variable 
;  spec2 : Returns all the spectral snapshots of the second variable 
;  i     : The index of the last time is i-2
;  tt    : Returns the times for the different snapshots (vector)
;  noplot: Do not plot if set
;  tmin  : First time for plotting snapshots  (if /all is set) 
;  tmax  : Last time for plotting snapshots  (if /all is set) 
;  tot   : Plots total power spectrum if tot=1
;  lin   : Plots the line k^lin
;  png   : to write png file for making a movie
;  yrange: y-range for plot
;  compensate: exponent for compensating power spectrum (default=0)
;  lint_shell: shell-integrated spectrum (default=1)
;  lint_z    : z-integrated spectrum (default=0)
;  print     : flag for print into PS file (default=0)
;
;  24-sep-02/nils: coded
;   5-oct-02/axel: comments added
;  29-nov-10/MR  : adaptions to different types of spectra (with/without headers), 
;                  keyword parameters lint_shell, lint_z and print added
;
common pars,  nx, ny, nz, nk, ncomp, nt, Lx, Ly, Lz
common wavenrs, kxs, kys, kshell

default,var1,'u'
default,var2, ''      ;'b'
default,last,1
default,w,0.1
default,tmin,0
default,tmax,1e34
default,tot,0
default,lin,0
default,dir,''
default,compensate1,0
default,compensate2,compensate1
default,compensate,compensate1
default,datatopdir,'data'
default,lint_shell,1
default,lint_z,0
default, prnt, 0
;
;  This is done to make the code backward compatible.
;
if  keyword_set(all) then begin
    last=0
end
if  keyword_set(wait) then begin
    w=wait
end
if  keyword_set(v1) then begin
    file1='power'+v1
    if  keyword_set(v2) then begin
        file2='power'+v2
    end else begin
        file2=''
    end
end else begin
    if  keyword_set(v2) then begin
        print,'In order to set v2 you must also set v1!'
        print,'Exiting........'
        stop
    end
    file1='power'+var1
    if var2 ne '' then $
      file2='power'+var2 $
    else $
      file2=''
end 

;;file1=file1+'_xy'
;;if file2 ne '' then file2=file2+'_xy'

file1=file1+'.dat'
if file2 ne '' then file2=file2+'.dat'
;
;  plot only when iplot=1 (default)
;  can be turned off by using /noplot
;
if keyword_set(noplot) then iplot=0 else iplot=1
;
;!p.multi=[0,1,1]
;!p.charsize=2

;
;  Reading number of grid points from 'data/dim.dat'
;
  pc_read_param,o=param,/quiet
  
  Lx=param.Lxyz[0]
  Ly=param.Lxyz[1]
  Lz=param.Lxyz[2]

  pc_read_dim,obj=param,/quiet
  nx=param.nx
  ny=param.ny
  nz=param.nz
;
;  Calculating some variables
;
first='true'

kxs = fltarr(nx) & kys = fltarr(ny) 
  
k0=2.*!pi/Lz
nk0=round( sqrt( ((nx+1)*!pi/Lx)^2+((ny+1)*!pi/Ly)^2)/(2*!pi/Lx) )+1 
nk=nk0

;if  keyword_set(v1) then begin
;  if ((v1 EQ "_phiu") OR (v1 EQ "_phi_kin") $
;       OR (v1 EQ "hel_phi_kin") OR (v1 EQ "hel_phi_mag") $
;      OR (v1 EQ "_phib") OR (v1 EQ "_phi_mag") ) then begin
;    k = k*k0   
;  end
;end

lcomplex=0
nt1 = read_firstpass( datatopdir+'/'+file1, lint_shell, lint_z, lcomplex, global_ext1, startpos1, fmt=fmt1, zpos=zpos1 )

if nt1 eq 0 then begin
  print, 'Error when reading '+datatopdir+'/'+file1+'!'
  stop
endif else $
  nk1 = nk

default,yrange,[10.0^(floor(alog10(float(global_ext1(0))))),10.0^ceil(alog10(float(global_ext1(1))))]

;
;  Reading file 2 if it is defined
;
if (file2 ne '') then begin

  nk = nk0
  nt2 = read_firstpass( datatopdir+'/'+file2, lint_shell, lint_z, lcomplex, global_ext2, startpos2, fmt=fmt2, zpos=zpos2 )
 
  if nt2 eq 0 or nk ne nk1 then begin

    if nt2 eq 0 then $ 
      print, 'Error: No data readable from file '+datatopdir+'/'+file2+'!' $
    else $
      print, 'Error: Number of shell wavenumbers different in files '+datatopdir+'/'+file1+' and '+datatopdir+'/'+file2+'!'

    print, '       File2 ignored.' 
    file2 = ''
    close, 1
    nk = nk1
    nt = nt1

  endif else begin

    if nt1 ne nt2 then begin 
      print, 'Warning: Number of time steps different in files '+datatopdir+'/'+file1+' and '+datatopdir+'/'+file2+'!'
      print, '         Adopting minimum.' 
      nt = nt1<nt2
    endif else $
      nt = nt1
    spec2=1
    alloc_specs, spectrum2, lint_shell, lint_z, lcomplex, fullspec=spec2 
    openr, 2, datatopdir+'/'+file2
    point_lun, 2, startpos2(0)

  endelse

endif

tt=fltarr(nt)
lasti=nt-1
spec1=1
alloc_specs, spectrum1, lint_shell, lint_z, lcomplex, fullspec=spec1
;
;  Plotting the results for last time frame
;
;  check whether we want png files (for movies)
;
if lint_shell then begin

  if keyword_set(png) then begin
  
    set_plot, 'z'                   ; switch to Z buffer
    device, SET_RESOLUTION=[!d.x_size,!d.y_size] ; set window size
    itpng=0 ;(image counter)
    ;
    ;  set character size to bigger values
    ;
    !p.charsize=2
    !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
  
  endif else if prnt then begin
      set_plot, 'PS'
      device, file='power.ps'
  endif else begin
      set_plot, 'x'   
  endelse
  
  !x.title='!8k!3'
  fo='(f4.2)'
  if compensate eq 0. then !y.title='!8E!3(!8k!3)' else !y.title='!8k!6!u'+string(compensate,fo=fo)+' !8E!3(!8k!3)'
  !x.range=[1,nk*k0]

endif

  openr, 1, datatopdir+'/'+file1
  point_lun, 1, startpos1(0) & iseg=0

  it=0L
  while not eof(1) do begin 
    
    on_ioerror, newheader
    readf,1,time
    on_ioerror, NULL
    tt(it)=time

    for ic=0,ncomp-1 do begin
      if lcomplex then $
        readf,1,spectrum1, format=fmt1 $
      else $
        readf,1,spectrum1
    
      if lint_shell then $
        spec1(*,*,ic,it)=spectrum1 $
      else if lint_z then $
        spec1(*,*,ic,it)=shift(spectrum1,nx/2,ny/2) $     ; why transpose necessary?
      else $
        spec1(*,*,*,ic,it)=spectrum1
    
      maxy=max(spectrum1(1:*))
      miny=min(spectrum1(1:*))

    end

    if (file2 ne '') then begin

      readf,2,time
      if time ne tt(it) then $
        print, 'Warning: Times in '+datatopdir+'/'+file1+' and '+datatopdir+'/'+file2+' different!'

      if lcomplex then $
        readf,2,spectrum2,format=fmt2 $
      else $
        readf,2,spectrum2
       

      if lint_shell then $
        spec2(*,*,it)=spectrum2 $
      else if lint_z then $
        spec2(*,*,it)=shift(spectrum2,nx/2,ny/2) $     ; why transpose necessary?
      else $
        spec2(*,*,*,it)=spectrum2

    endif
    ;
    ;  normalize?
    ;
    if keyword_set(norm) then begin
      spectrum2=spectrum2/total(spectrum2)
      print,'divide spectrum2 by total(spectrum2)'
    endif
    
    if (last eq 0) then begin 

    ; creates movie with time dependent spectra (only implemented for lint_shell=T, lint_z=T!)
     
      if (time ge tmin and time le tmax) then begin

    	if lint_shell and lint_z and not lcomplex then begin

    	  xrr=[1,nk*k0]
    	  yrr=global_ext
    	  !p.title='t='+str(time)

    	  if iplot eq 1 then begin

    	    plot_oo,kshell,spectrum1*kshell^compensate1,back=255,col=0,yr=yrange
    	    oplot,kshell,spectrum2*kshell^compensate1,col=122
             ;xx=[1,5] & oplot,xx,2e-6*xx^1.5,col=55
             ;xyouts,2,1e-5,'!8k!6!u3/2!n',siz=1.8,col=55
            xx=[1,5] & oplot,xx,2e-8*xx^1.5,col=55
            xyouts,2,1e-7,'!8k!6!u3/2!n',siz=1.8,col=55
               
	    if (file2 ne '') then begin
              ;
              ; possibility of special settings for helicity plotting
              ; of second variable
              ;
              if keyword_set(helicity2) then begin
                oplot,kshell,.5*abs(spectrum2)*kshell^compensate2,col=122
                oplot,kshell,+.5*spectrum2*kshell^compensate2,col=122,ps=6
                oplot,kshell,-.5*spectrum2*kshell^compensate2,col=55,ps=6
              endif else begin
                oplot,kshell,abs(spectrum2)*kshell^compensate2,col=122
              endelse
	
              if (tot eq 1) then $
                oplot,kshell,(spectrum1+spectrum2)*kshell^compensate,col=47
        
            endif
	    
            if (lin ne 0) then begin
              fac=spectrum1(2)/kshell(2)^(lin)*1.5
              oplot,kshell(2:*),kshell(2:*)^(lin)*fac,lin=2,col=0
            endif
            wait,w
          endif 

        endif else begin
          if lint_z then $
            stop, 'Warning: No implementation for movie with lint_shell=F, lint_z=T!' $
          else $
            stop, 'Warning: No implementation for movie with lint_shell=T, lint_z=F!'
        endelse
      endif
    endif
    
    ;
    ;  check whether we want to write png file (for movies)
    ;
    if keyword_set(png) then begin
      istr2 = strtrim(string(itpng,'(I20.4)'),2) ;(only up to 9999 frames)
      image = tvrd()
      ;
      ;  write png file
      ;
      tvlct, red, green, blue, /GET
      imgname = dir+'img_'+istr2+'.png'
      write_png, imgname, image, red, green, blue
      print,'itpng=',itpng
      itpng=itpng+1 ;(counter)
    endif
    ;
    it=it+1L
    if it eq nt then break
    continue
newheader:
    if eof(1) then break
    iseg += 1
    point_lun, 1, startpos1(iseg)

  endwhile
  
;stop,'AXEL1'
  if arg_present(obj) then $
    if n_elements(zpos1) gt 0 then $
      obj = {FILE: file1, TT: tt, ZPOS: zpos1, SPEC1: reform(spec1)} $
    else $
      obj = {FILE: file1, TT: tt, SPEC1: reform(spec1)}
; MR: would prefer the name 'T' for the time as it is also use elsewhere
close, 1
return

  if not lint_shell then begin
    if lint_z then begin
  
      ;goto, loop1
      
      set_plot, 'X'
      
      for it=0,n_elements(tt)-1 do begin 
        contour, spec1(*,*,it),kxs,kys, nlevels=30, /fill, c_colors=8.*indgen(30), xrange=[-15.,15.], yr=[-5.,5.], xtitle='kx', ytitle='ky'
        xyouts, 15.1, 4., string(tt(it))
        wait, .1
      endfor
      ;goto, loop2;
loop1:
  
      if keyword_set(png) then begin
      
        set_plot, 'z'                   ; switch to Z buffer
        device, SET_RESOLUTION=[!d.x_size,!d.y_size] ; set window size
        itpng=0 ;(image counter)
        ;
        ;  set character size to bigger values
        ;
        !p.charsize=2
        !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
      
      endif else if prnt then begin
        set_plot, 'PS'
        device, file='power.ps'
      endif else begin
        set_plot, 'X'   
      endelse
  
      ikx0 = 0
      while ikx0 ge 0 do begin
      
        read, ikx0, iky0, prompt='Wellenzahlen:'
        
        if ikx0 ge 0 and iky0 ge 0 then begin
        
          ikx1 = 16-ikx0 & ikx2 = 16+ikx0
          iky1 = 16-iky0 & iky2 = 16+iky0
          
          for ikx=ikx1,ikx2,8 do begin
            
            if ikx eq ikx1 then $
              exinds1 = gen_exinds(spec1(ikx,iky2,*)) $
            else $
              exinds = gen_exinds(spec1(ikx,iky2,*))
          
          endfor 
loop2:       
          itmin=0  & itmax=n_elements(tt)-1
          while itmax ge 0 do begin
          
            read, itmin, prompt='itmin= (>0)'
            read, itmax, prompt='itmax= (<'+string(n_elements(tt), format='(i4)')+', Stop by negative itmax)'

            if itmin ge 0 and itmax ge 0 then begin
            
              plot_segs, tt(itmin:itmax), spec1(ikx1,iky2,itmin:itmax), exinds1-itmin, colors=[!p.color], /overplot
              ;oplot, tt(itmin:itmax), spec1(ikx1,iky2,itmin:itmax)      
              ;oplot, tt(itmin:itmax), spec1(ikx2,iky1,itmin:itmax), linest=3, color=250

              if ikx2 eq ikx1+8 then begin
              
                plot_segs, tt(itmin:itmax), spec1(ikx2,iky2,itmin:itmax), exinds-itmin, colors=[!p.color], /ylog, xtitle='!8t', ytitle='!8E!Dk!N!3'
                ;plot , tt(itmin:itmax), spec1(ikx2,iky2,itmin:itmax), /ylog, xtitle='!8t', ytitle='!8E!Dk!N!3'
                ;oplot, tt(itmin:itmax), spec1(ikx1,iky1,itmin:itmax), linest=3, color=250
              
                xyouts, .5*tt(itmax), 0.5*spec1(ikx2,iky2,0.5*itmax), '++'
                xyouts, .5*tt(itmax),  3.*spec1(ikx2,iky1,0.5*itmax), '+-'

                xyouts, .8*tt(itmax), 1.e-10, '!7k!D!8x!N!3='+string(2.*!pi/kxs(ikx2)/Lx, format='(f4.1)')+'!8L!Dx!N!3'
                xyouts, .8*tt(itmax), 1.e-11, '!7k!D!8y!N!3='+string(2.*!pi/kys(iky2)/Ly, format='(f4.1)')+'!8L!Dy!N!3'

                openw, 12, 'powjxb.dat'
                printf, 12, itmin, itmax
                printf, 12, spec1(ikx2,iky2,itmin:itmax)
                printf, 12, spec1(ikx1,iky2,itmin:itmax)
                close, 12
              endif
            endif  
          endwhile
        endif
      endwhile
      stop
cont1:    
    endif
    
  endif else $
  
    if (last eq 1 and iplot eq 1) then begin
    
      if lint_shell and lint_z then begin

        !p.title='t=' + string(time)
        !y.range=[miny,maxy]
        
        plot_oo,kshell,spectrum1*kshell^compensate,back=255,col=0,yr=yrange
        
        if (file2 ne '') then begin
        
          oplot,kshell,spectrum2*kshell^compensate,col=122
          if (tot eq 1) then $
            oplot,kshell,(spectrum1+spectrum2)*kshell^compensate,col=47
    
        endif
        
        if (lin ne 0) then begin
          fac=spectrum1(2)/kshell(2)^(lin)*1.5
          oplot,k(2:*),kshell(2:*)^(lin)*fac,lin=2,col=0
        endif

      endif else begin
        if lint_z then $
          stop, 'Warning: No implementation for last spectrum plot with lint_shell=F, lint_z=T!' $
        else $
          stop, 'Warning: No implementation for last spectrum plot with lint_shell=T, lint_z=F!'
      endelse
      
    endif
    
close,1
close,2

if !d.name eq 'PS' then device, /close

!x.title='' 
!y.title=''
!p.title=''
!x.range=''
!y.range=''

END

