PRO power,var1,var2,last,w,v1=v1,v2=v2,all=all,wait=wait,k=k,spec1=spec1, $
          spec2=spec2,i=i,tt=tt,noplot=noplot,tmin=tmin,tmax=tmax, $
          tot=tot,lin=lin,png=png,yrange=yrange,norm=norm,helicity2=helicity2, $
          compensate1=compensate1,compensate2=compensate2,datatopdir=datatopdir
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
;
;  24-sep-02/nils: coded
;   5-oct-02/axel: comments added
;
default,var1,'u'
default,var2,'b'
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
    file1='power'+v1+'.dat'
    if  keyword_set(v2) then begin
        file2='power'+v2+'.dat'
    end else begin
        file2=''
    end
end else begin
    if  keyword_set(v2) then begin
        print,'In order to set v2 you must also set v1!'
        print,'Exiting........'
        stop
    end
    file1='power'+var1+'.dat'
    file2='power'+var2+'.dat'
end 
;
;  plot only when iplot=1 (default)
;  can be turned off by using /noplot
;
if keyword_set(noplot) then iplot=0 else iplot=1
;
;!p.multi=[0,1,1]
;!p.charsize=2

mx=0L & my=0L & mz=0L & nvar=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
;  Reading number of grid points from 'data/dim.dat'
;  Need both mx and nghostx to work out nx.
;  Assume nx=ny=nz
;
close,1
openr,1,datatopdir+'/'+'dim.dat'
readf,1,mx,my,mz,nvar
readf,1,prec
readf,1,nghostx,nghosty,nghostz
close,1
size=2*!pi
;
;  Calculating some variables
;
first='true'
nx=mx-nghostx*2
if  keyword_set(v1) then begin
  if ((v1 EQ "_phiu") OR (v1 EQ "_phi_kin") $
       OR (v1 EQ "hel_phi_kin") OR (v1 EQ "hel_phi_mag") $
       OR (v1 EQ "_phib") OR (v1 EQ "_phi_mag") ) then begin
     nx=mz-nghostz*2
     pc_read_grid,o=grid,/quiet
     size=grid.Lz
  end
end
;print,'nx=',nx
imax=nx/2
spectrum1=fltarr(imax)
if n_elements(size) eq 0 then size=2.*!pi
k0=2.*!pi/size
wavenumbers=indgen(imax)*k0 
k=findgen(imax)+1.
k=findgen(imax)
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
close,1
openr,1, datatopdir+'/'+file1
  while not eof(1) do begin
    readf,1,time
    readf,1,spectrum1
    if (max(spectrum1(1:*)) gt globalmax) then globalmax=max(spectrum1(1:*))
    if (min(spectrum1(1:*)) lt globalmin) then globalmin=min(spectrum1(1:*))
    i=i+1L
  endwhile
close,1
spec1=fltarr(imax,i-1)
tt=fltarr(i-1)
lasti=i-2
default,yrange,[10.0^(floor(alog10(min(spectrum1(1:*))))),10.0^ceil(alog10(max(spectrum1(1:*))))]
;
;  Opening file 2 if it is defined
;
if (file2 ne '') then begin
  close,2
  spectrum2=fltarr(imax)
  openr,2,datatopdir+'/'+file2
  spec2=fltarr(imax,i-1)
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
openr,1, datatopdir+'/'+file1
    while not eof(1) do begin 
      	readf,1,time
       	readf,1,spectrum1
	tt(i-1)=time
       	spec1(*,i-1)=spectrum1
       	maxy=max(spectrum1(1:*))
       	miny=min(spectrum1(1:*))
       	if (file2 ne '') then begin
	  readf,2,time
	  readf,2,spectrum2
          spec2(*,i-1)=spectrum2
          if (max(spectrum2(1:*)) gt maxy) then maxy=max(spectrum2(1:*))
          if (min(spectrum2(1:*)) lt miny) then miny=min(spectrum2(1:*))
       	endif
        ;
        ;  normalize?
        ;
        if keyword_set(norm) then begin
          spectrum2=spectrum2/total(spectrum2)
          print,'divide spectrum by total(spectrum2)'
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
         	if (file2 ne '') then begin
                  ;
		  ; possibility of special settings for helicity plotting
		  ; of second variable
                  ;
                  if keyword_set(helicity2) then begin
		    oplot,k,.5*abs(spectrum2)*k^compensate2,col=122
		    oplot,k,+.5*spectrum2*k^compensate2,col=122,ps=6
		    oplot,k,-.5*spectrum2*k^compensate2,col=55,ps=6
                  endif else begin
		    oplot,k,abs(spectrum2)*k^compensate2,col=122
                  endelse
		  if (tot eq 1) then begin
		    oplot,k,(spectrum1+spectrum2)*k^compensate,col=47
		  endif
		endif
	        if (lin ne 0) then begin
		  fac=spectrum1(2)/k(2)^(lin)*1.5
		  oplot,k(2:*),k(2:*)^(lin)*fac,lin=2,col=0
		endif
              	wait,w
	      endif
            endif
	  endif
       	endif
       	i=i+1L
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
    endwhile
    if (last eq 1 and iplot eq 1) then begin
	!p.title='t=' + string(time)
	!y.range=[miny,maxy]
	plot_oo,k,spectrum1*k^compensate,back=255,col=0,yr=yrange
      	if (file2 ne '') then begin
		oplot,k,spectrum2*k^compensate,col=122
		if (tot eq 1) then begin
			oplot,k,(spectrum1+spectrum2)*k^compensate,col=47
		endif
	endif
	if (lin ne 0) then begin
		fac=spectrum1(2)/k(2)^(lin)*1.5
		oplot,k(2:*),k(2:*)^(lin)*fac,lin=2,col=0
	endif
    endif
close,1
close,2

!x.title='' 
!y.title=''
!p.title=''
!x.range=''
!y.range=''
END





