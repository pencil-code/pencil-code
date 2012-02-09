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

PRO pc_power_xy,var1,var2,last,w,v1=v1,v2=v2,all=all,wait=wait,k=k,spec1=spec1, $
          spec2=spec2,i=i,tt=tt,noplot=noplot,tmin=tmin,tmax=tmax, $
          tot=tot,lin=lin,png=png,yrange=yrange,norm=norm,helicity2=helicity2, $
          compensate1=compensate1,compensate2=compensate2,datatopdir=datatopdir, $ 
	  lint_shell=lint_shell, lint_z=lint_z, print=prnt
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
    file2='power'+var2
end 

file1=file1+'_xy'
if file2 ne '' then file2=file2+'_xy'

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
;
;  Calculating some variables
;
first='true'
nx=mx-nghostx*2
ny=my-nghosty*2
nz=mz-nghostz*2
;print,'nx=',nx

;if  keyword_set(v1) then begin
  ;if ((v1 EQ "_phiu") OR (v1 EQ "_phi_kin") $
  ;     OR (v1 EQ "hel_phi_kin") OR (v1 EQ "hel_phi_mag") $
  ;     OR (v1 EQ "_phib") OR (v1 EQ "_phi_mag") ) then begin
     pc_read_grid,o=grid,/quiet
     
     Lz=grid.Lz
     Lx=grid.Lx
     Ly=grid.Ly
  ;end
;end

kxs = fltarr(nx) & kys = fltarr(ny) 
  
k0=2.*!pi/Lz

nk=round( sqrt( ((nx+1)*!pi/Lx)^2+((ny+1)*!pi/Ly)^2)/(2*!pi/Lx) )+1 

;if  keyword_set(v1) then begin
;  if ((v1 EQ "_phiu") OR (v1 EQ "_phi_kin") $
;       OR (v1 EQ "hel_phi_kin") OR (v1 EQ "hel_phi_mag") $
;      OR (v1 EQ "_phib") OR (v1 EQ "_phi_mag") ) then begin
;    k = k*k0   
;  end
;end
;
;  Looping through all the data to get number of spectral snapshots
;
close,1

openr,1, datatopdir+'/'+file1

  point_lun, -1, pos0
  headline = ''
  readf, 1, headline
  witheader=strpos(strlowcase(headline),'power spectrum') ne -1
  
  if witheader then begin
    
    if (strpos(strlowcase(headline),'shell-integrated') ne -1) ne lint_shell then begin
    
      if lint_shell then $
        print, 'Warning: File header says "not shell-integrated" - assume that' $
      else $
        print, 'Warning: File header says "shell-integrated" - assume that'
        
      lint_shell = ~lint_shell
      
    endif
    
    if (strpos(strlowcase(headline),'z-integrated') ne -1) ne lint_z then begin

      if lint_z then $
        print, 'Warning: File header says "not z-integrated" - assume that' $
      else $
        print, 'Warning: File header says "z-integrated" - assume that'
        
      lint_z = ~lint_z
      
    endif
  
    readf, 1, headline

    if lint_shell then kshell = fltarr(nk)
 
    if strpos(strlowcase(headline),'wavenumbers') eq -1 then $
      print, 'Warning: File header corrupt!'
     
    readf, 1, kxs, kys 
    if (lint_shell) then begin
    
      readf, 1, headline
      if strpos(strlowcase(headline),'wavenumbers') eq -1 then $
        print, 'Warning: File header corrupt!'
        
      readf, 1, kshell 
    endif
    
    point_lun, -1, pos 
   
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
 
    pos=0L
    point_lun, 1, pos
    
  endelse

  if lint_shell then begin
    
    if lint_z then $
      spectrum1=fltarr(nk) $
    else $
      spectrum1=fltarr(nk,nz)
      
  endif else if lint_z then $
    spectrum1=fltarr(nx,ny) $
  else $
    spectrum1=fltarr(nx,ny,nz)
    
  globalmin=1e12
  globalmax=1e-30
  nt=0L & time=0.
  
  while not eof(1) do begin
  
    if witheader then begin
      
      point_lun, -1, pos
      readf,1,headline
    
      if strpos(strlowcase(headline),'power spectrum') eq -1 then $
        point_lun, 1, pos $
      else begin
        readf,1,headline
        readf,1, kxs, kys
        if (lint_shell) then begin
        
          readf, 1, headline
          if strpos(strlowcase(headline),'wavenumbers') eq -1 then $
            print, 'Warning: File header corrupt!'
            
          readf, 1, kshell 
        endif
      endelse
 
      readf,1,time
     
    endif else $
      readf,1,time
      
    readf,1,spectrum1
    
    if (max(spectrum1(1:*,*)) gt globalmax) then globalmax=max(spectrum1(1:*,*))
    if (min(spectrum1(1:*,*)) lt globalmin) then globalmin=min(spectrum1(1:*,*))
    nt=nt+1L
  
  endwhile
  
; end first reading
 
if lint_shell then begin
  if lint_z then $
    spec1=fltarr(nk,nt) $
  else $
    spec1=fltarr(nk,nz,nt)
endif else if lint_z then $
  spec1=fltarr(nx,ny,nt) $
else $
  spec1=fltarr(nx,ny,nz,nt)
  
tt=fltarr(nt)
lasti=nt-1

default,yrange,[10.0^(floor(alog10(min(spectrum1(1:*))))),10.0^ceil(alog10(max(spectrum1(1:*))))]

;
;  Opening file 2 if it is defined
;
if (file2 ne '') then begin
  close,2
  if lint_shell then $
    spectrum2=fltarr(nk,nz)
  openr,2,datatopdir+'/'+file2, ERROR = err
  if err eq 0 then $
    spec2=fltarr(nk,nz,nt) $
  else $
    file2=''
endif
;
;  Plotting the results for last time frame
;;
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

  i=0L
  point_lun, 1, pos0
  while not eof(1) do begin 
 
    if witheader then begin
      
      point_lun, -1, pos
      readf,1,headline
    
      if strpos(strlowcase(headline),'power spectrum') eq -1 then $
        point_lun, 1, pos $
      else begin
        
        readf,1,headline
        readf,1, kxs, kys
        
        if (lint_shell) then begin
          readf, 1, headline
          readf, 1, kshell 
        endif

      endelse
      readf,1,time
      
    endif else $
      readf,1,time
  
    readf,1,spectrum1
    tt(i)=time
    
    if lint_shell then $
      spec1(*,*,i)=spectrum1 $
    else if lint_z then $
      spec1(*,*,i)=shift(spectrum1,nx/2,ny/2) $     ; why transpose necessary?
    else $
      spec1(*,*,*,i)=spectrum1
    
    maxy=max(spectrum1(1:*))
    miny=min(spectrum1(1:*))

    if (file2 ne '') then begin
      readf,2,time
      readf,2,spectrum2
      spec2(*,*,i)=spectrum2
      if (max(spectrum2(1:*,*)) gt maxy) then maxy=max(spectrum2(1:*,*))
      if (min(spectrum2(1:*,*)) lt miny) then miny=min(spectrum2(1:*,*))
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
    	  if lint_shell then begin
    	    xrr=[1,nk*k0]
    	    yrr=[globalmin,globalmax]
    	    !p.title='t='+str(time)
    	    if iplot eq 1 then begin
    	      plot_oo,k,spectrum1*k^compensate1,back=255,col=0,yr=yrange
    	      oplot,k,spectrum2*k^compensate1,col=122
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
        	  oplot,k,.5*abs(spectrum2)*k^compensate2,col=122
        	  oplot,k,+.5*spectrum2*k^compensate2,col=122,ps=6
        	  oplot,k,-.5*spectrum2*k^compensate2,col=55,ps=6
        	endif else begin
        	  oplot,k,abs(spectrum2)*k^compensate2,col=122
        	endelse
		
        	if (tot eq 1) then $
        	  oplot,k,(spectrum1+spectrum2)*k^compensate,col=47
        	
              endif
	      
              if (lin ne 0) then begin
        	fac=spectrum1(2)/k(2)^(lin)*1.5
                oplot,k(2:*),k(2:*)^(lin)*fac,lin=2,col=0
              endif
              wait,w
            endif 
          endif else if lint_z then begin
               stop
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
  
  if witheader then begin
    kxs = shift(kxs,nx/2)
    kys = shift(kys,ny/2)  
  endif
   
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
    
      !p.title='t=' + string(time)
      !y.range=[miny,maxy]
      
      plot_oo,kshell,spectrum1*kshell^compensate,back=255,col=0,yr=yrange
      
      if (file2 ne '') then begin
      
    	oplot,k,spectrum2*kshell^compensate,col=122
    	if (tot eq 1) then $
    	  oplot,kshell,(spectrum1+spectrum2)*kshell^compensate,col=47
    	
      endif
      
      if (lin ne 0) then begin
    	fac=spectrum1(2)/kshell(2)^(lin)*1.5
    	oplot,k(2:*),kshell(2:*)^(lin)*fac,lin=2,col=0
      endif
      
    endif
    
cont:
close,1
close,2

if !d.name eq 'PS' then device, /close

!x.title='' 
!y.title=''
!p.title=''
!x.range=''
!y.range=''

END

