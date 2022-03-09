PRO power,var1,var2,last,w,v1=v1,v2=v2,v3=v3,all=all,wait=wait,k=k, $
          spec1=spec1,spec2=spec2,spec3=spec3,scal2=scal2,scal3=scal3, $
          i=i,tt=tt,noplot=noplot,tmin=tmin,tmax=tmax, $
          tot=tot,lin=lin,png=png,yrange=yrange,norm=norm,helicity2=helicity2, $
          compensate1=compensate1,compensate2=compensate2, $
          compensate3=compensate3,datatopdir=datatopdir,double=double, $
          lkscale=lkscale,cyl=cyl,zwav=zwav
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
;  v3    : Third variable to be plotted (Ex: 'cc') 
;  all   : Plot all snapshots if /all is set, otherwise only last snapshot
;  wait  : Time to wait between each snapshot (if /all is set) 
;  k     : Returns the wavenumber vector
;  spec1 : Returns all the spectral snapshots of the first variable 
;  spec2 : Returns all the spectral snapshots of the second variable 
;  spec3 : Returns all the spectral snapshots of the third variable
;  scal2 : Scaling factor for the second spectrum
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
;  lkscale: apply scaling factor if box is not (2pi)^3
;
;  24-sep-02/nils: coded
;   5-oct-02/axel: comments added
;
default,file3,''
default,var1,'u'
default,var2,'b'
default,last,1
default,w,0.1
default,tmin,0
default,tmax,1e34
default,scal2,1.
default,scal3,1.
default,tot,0
default,lin,0
default,dir,''
default,compensate1,0
default,lkscale,0
default,compensate2,compensate1
default,compensate3,compensate1
default,compensate,compensate1
default,datatopdir,'data'
;
;  This is done to make the code backward compatible.
;
if  keyword_set(all) then begin
    last=0
end
if  keyword_set(wait) then begin
    w=wait
end
if  n_elements(v1) eq 1L then begin
    file1='power'+v1+'.dat'
;
;  second spectrum
;
    if n_elements(v2) eq 1L then begin
        file2='power'+v2+'.dat'
    end else begin
        file2=''
    end
;
;  third spectrum
;
    if n_elements(v3) eq 1L then begin
        file3='power'+v3+'.dat'
    end else begin
        file3=''
    end
;
end else begin
    if n_elements(v2) eq 1L then begin
        print,'In order to set v2 you must also set v1!'
        print,'Exiting........'
        stop
    end
    if keyword_set(cyl) then begin
      file1='cyl_power'+var1+'.dat'
      file2='cyl_power'+var2+'.dat'
    endif else begin
      file1='power'+var1+'.dat'
      file2= var2 eq '' ? '' : 'power'+var2+'.dat'
    endelse
end 
;
;  plot only when iplot=1 (default)
;  can be turned off by using /noplot
;
iplot = keyword_set(noplot) ? 0 : 1
;
;!p.multi=[0,1,1]
;!p.charsize=2

;  Reading number of grid points from 'data/dim.dat'
;  Need both mx and nghostx to work out nx.
;  Assume nx=ny=nz
;
pc_read_dim, obj=dim, datadir=datatopdir, /quiet
;
;  Calculating some variables
;
first='true'
nx=dim.mx-dim.nghostx*2
size=2*!pi
if  keyword_set(v1) then begin
  if ((v1 EQ "_phiu") OR (v1 EQ "_phi_kin") $
       OR (v1 EQ "hel_phi_kin") OR (v1 EQ "hel_phi_mag") $
       OR (v1 EQ "_phib") OR (v1 EQ "_phi_mag") ) then begin
     nx=dim.mz-dim.nghostz*2
     pc_read_grid,o=grid,dim=dim,/quiet,datadir=datatopdir
     size=grid.Lz
  end
end
imax=nx/2
if keyword_set(cyl) then begin
  if keyword_set(double) then begin
    spectrum1=dblarr(imax,dim.nz)
  endif else begin
    spectrum1=fltarr(imax,dim.nz)
  endelse
endif else begin
  if keyword_set(double) then begin
    spectrum1=dblarr(imax)
  endif else begin
    spectrum1=fltarr(imax)
  endelse
endelse
;
if n_elements(size) eq 0 then size=2.*!pi
k0=2.*!pi/size
wavenumbers=indgen(imax)*k0 
;k=findgen(imax)+1.
k=findgen(imax)
zwav=findgen(dim.nz)-dim.nz/2
if  keyword_set(v1) then begin
  if ((v1 EQ "_phiu") OR (v1 EQ "_phi_kin") $
       OR (v1 EQ "hel_phi_kin") OR (v1 EQ "hel_phi_mag") $
       OR (v1 EQ "_phib") OR (v1 EQ "_phi_mag") ) then begin
    k = k*k0   
  end
end
;
;  Looping through all the data to get number of spectral snapshots
;
globalmin=1e12
globalmax=1e-30
i=1L
openr, unit, datatopdir+'/'+file1, /get_lun
  while ~eof(unit) do begin
    readf,unit,time
    readf,unit,spectrum1
    ;if (max(spectrum1(1:*)) gt globalmax) then globalmax=max(spectrum1(1:*))
    ;if (min(spectrum1(1:*)) lt globalmin) then globalmin=min(spectrum1(1:*))
    ;AB: can include 0, I think.
    if (max(spectrum1) gt globalmax) then globalmax=max(spectrum1)
    if (min(spectrum1) lt globalmin) then globalmin=min(spectrum1)
    i++
  endwhile
free_lun, unit
if keyword_set(cyl) then begin
  if keyword_set(double) then begin
    spec1=dblarr(imax,dim.nz,i-1)
  endif else begin
    spec1=fltarr(imax,dim.nz,i-1)
  endelse
endif else begin
  if keyword_set(double) then begin
    spec1=dblarr(imax,i-1)
  endif else begin
    spec1=fltarr(imax,i-1)
  endelse
endelse
;
tt=fltarr(i-1)
lasti=i-2
default,yrange,[10.0^(floor(alog10(min(spectrum1(1:*))))),10.0^ceil(alog10(max(spectrum1(1:*))))]
;
;  Opening file 2 if it is defined
;
unit_2__open = 0L
if (file2 ne '') then begin
  if keyword_set(cyl) then begin
    if keyword_set(double) then begin
      spectrum2=dblarr(imax,dim.nz)
      spec2=dblarr(imax,dim.nz,i-1)
    endif else begin
      spectrum2=fltarr(imax,dim.nz)
      spec2=fltarr(imax,dim.nz,i-1)
    endelse
  endif else begin
    if keyword_set(double) then begin
      spectrum2=dblarr(imax)
      spec2=dblarr(imax,i-1)
    endif else begin
      spectrum2=fltarr(imax)
      spec2=fltarr(imax,i-1)
    endelse
  endelse
  openr, unit_2,datatopdir+'/'+file2, /get_lun
  unit_2__open = 1L
endif
;
;  Opening file 3 if it is defined
;
unit_3__open = 0L
if (file3 ne '') then begin
  if keyword_set(cyl) then begin
    if keyword_set(double) then begin
      spectrum3=dblarr(imax,dim.nz)
      spec3=dblarr(imax,dim.nz,i-1)
    endif else begin
      spectrum3=fltarr(imax,dim.nz)
      spec3=fltarr(imax,dim.nz,i-1)
    endelse
  endif else begin
    if keyword_set(double) then begin
      spectrum3=dblarr(imax)
      spec3=dblarr(imax,i-1)
    endif else begin
      spectrum3=fltarr(imax)
      spec3=fltarr(imax,i-1)
    endelse
  endelse
  openr, unit_3, datatopdir+'/'+file3, /get_lun
  unit_3__open = 1L
  ;spec3=fltarr(imax,i-1)
endif
;
;  Plotting the results for last time frame
;
!x.title='!8k!3'
fo='(f4.2)'
if compensate eq 0. then !y.title='!8P!3(!8k!3)' else !y.title='!8k!6!u'+string(compensate,fo=fo)+' !8P!3(!8k!3)'
!x.range=[1,imax*k0]
;
;  check whether we want png files (for movies)
;
if keyword_set(png) then begin
  set_plot, 'z'                   ; switch to Z buffer
  device, SET_RESOLUTION=[!d.x_size,!d.y_size] ; set window size
  itpng=0 ;(image counter)
  ;
  ;  set character size to bigger values
  ;
  !p.charsize=2
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
endif else if (!d.name eq 'PS') then begin
    set_plot, 'PS'
endif else begin
    set_plot, 'x'   
endelse
;
i=1L
openr, unit_1, datatopdir+'/'+file1, /get_lun
    while ~eof(unit_1) do begin
      	readf,unit_1,time
       	readf,unit_1,spectrum1
	tt(i-1)=time
  ;if keyword_set(cyl) then begin
  ;endif else begin
  ;endelse
        if keyword_set(cyl) then $
       	  spec1(*,*,i-1)=spectrum1 $
        else $
       	  spec1(*,i-1)=spectrum1
;
       	;maxy=max(spectrum1(1:*))
       	;miny=min(spectrum1(1:*))
       	maxy=max(spectrum1)
       	miny=min(spectrum1)
        ;
        ;  read second spectrum
        ;
       	if unit_2__open then begin

          on_ioerror, filend2
	  readf,unit_2,time
	  readf,unit_2,spectrum2
          if keyword_set(cyl) then $
            spec2(*,*,i-1)=spectrum2 $
          else $
            spec2(*,i-1)=spectrum2

          ;if (max(spectrum2(1:*)) gt maxy) then maxy=max(spectrum2(1:*))
          ;if (min(spectrum2(1:*)) lt miny) then miny=min(spectrum2(1:*))
          if (max(spectrum2) gt maxy) then maxy=max(spectrum2)
          if (min(spectrum2) lt miny) then miny=min(spectrum2)
          ;
          ;  normalize?
          ;
          if keyword_set(norm) then begin
            spectrum2=spectrum2/total(spectrum2)
            print,'divide spectrum by total(spectrum2)'
          endif
          goto, cont2
filend2:
          on_ioerror, NULL
          unit_2__open=0
          if i gt 1 then begin
            if keyword_set(cyl) then $
              spec2=spec2(*,*,0:i-2) $
            else $
              spec2=spec2(*,0:i-2)
          endif else $
            undefine, spec2
cont2:
       	endif
        ;
        ;  read third spectrum
        ;
       	if unit_3__open then begin
          on_ioerror, filend3
	  readf,unit_3,time
	  readf,unit_3,spectrum3
          if keyword_set(cyl) then $
            spec3(*,*,i-1)=spectrum3 $
          else $
            spec3(*,i-1)=spectrum3

          ;
          ;if (max(spectrum3(1:*)) gt maxy) then maxy=max(spectrum3(1:*))
          ;if (min(spectrum3(1:*)) lt miny) then miny=min(spectrum3(1:*))
          if (max(spectrum3) gt maxy) then maxy=max(spectrum3)
          if (min(spectrum3) lt miny) then miny=min(spectrum3)
          ;
          ;  normalize?
          ;
          if keyword_set(norm) then begin
            spectrum3=spectrum3/total(spectrum3)
            print,'divide spectrum by total(spectrum3)'
          endif
          goto, cont3
filend3:
          on_ioerror, NULL
          unit_3__open=0
          if i gt 1 then begin
            if keyword_set(cyl) then $
              spec3=spec3(*,*,0:i-2) $
            else $
              spec3=spec3(*,0:i-2)
          endif else $
            undefine, spec3
cont3:
       	endif
        ;
       	if (last eq 0) then begin
	  if (time ge tmin) then begin
	    if (time le tmax) then begin
	      xrr=[1,imax*k0]
	      yrr=[globalmin,globalmax]
	      !p.title='t='+str(time)
              if iplot eq 1 then begin
		plot_oo,k,spectrum1*k^compensate1,back=255,col=0,yr=yrange
                ;
                ;  second spectrum
                ;
         	if unit_2__open then begin
                  ;
		  ; possibility of special settings for helicity plotting
		  ; of second variable
                  ;
                  if keyword_set(helicity2) then begin
		    oplot,k,.5*scal2*abs(spectrum2)*k^compensate2,col=122
		    oplot,k,+.5*scal2*spectrum2*k^compensate2,col=122,ps=6
		    oplot,k,-.5*scal2*spectrum2*k^compensate2,col=55,ps=6
                  endif else begin
		    oplot,k,scal2*abs(spectrum2)*k^compensate2,col=122
                  endelse
		  if (tot eq 1) then begin
		    oplot,k,(spectrum1+spectrum2)*k^compensate,col=47
		  endif
		endif
                ;
                ;  third spectrum
                ;
         	if unit_3__open then begin
                  ;
		  ; possibility of special settings for helicity plotting
		  ; of second variable
                  ;
		  oplot,k,scal3*abs(spectrum3)*k^compensate3,col=55
		endif
                ;
	        if (lin ne 0) then begin
		  fac=spectrum1(2)/k(2)^(lin)*1.5
		  oplot,k(2:*),k(2:*)^(lin)*fac,lin=2,col=0
		endif
              	wait,w
	      endif
            endif
	  endif
       	endif
       	i++
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
          itpng++ ;(counter)
        endif
      ;
    endwhile
    if (last eq 1 and iplot eq 1) then begin

      !p.title='t=' + string(time)
      !y.range=[miny,maxy]
      plot_oo,k,spectrum1*k^compensate,back=255,col=0,yr=yrange
      if is_defined(spec2) then begin
        oplot,k,spectrum2*k^compensate,col=122
        if (tot eq 1) then oplot,k,(spectrum1+spectrum2)*k^compensate,col=47
      endif
      if (lin ne 0) then begin
        fac=spectrum1(2)/k(2)^(lin)*1.5
        oplot,k(2:*),k(2:*)^(lin)*fac,lin=2,col=0
      endif

    endif
if n_elements(unit_1) eq 1L then begin
  fstat = fstat(unit_1)
  if fstat.open then free_lun, unit_1
endif
if n_elements(unit_2) eq 1L then begin
  fstat = fstat(unit_2)
  if fstat.open then free_lun, unit_2
endif
if n_elements(unit_3) eq 1L then begin
  fstat = fstat(unit_3)
  if fstat.open then free_lun, unit_3
endif
;
;  scale k correctly
;
if lkscale then begin
  pc_read_param,obj=param,/quiet,datadir=datatopdir
  Lx=param.Lxyz[0]
  kscale_factor=2.*!pi/Lx
  k=k*kscale_factor
endif
;
!x.title='' 
!y.title=''
!p.title=''
!x.range=''
!y.range=''
END
