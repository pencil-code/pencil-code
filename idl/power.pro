PRO power,var1,var2,last,w,k=k,spec1=spec1,spec2=spec2,i=i,tt=tt, $
  noplot=noplot,tmin=tmin,tmax=tmax,tot=tot,lin=lin
;
;  $Id: power.pro,v 1.17 2003-03-05 09:32:18 nilshau Exp $
;
;  This routine reads in the power spectra generated during the run
;  (provided dspec is set to a time interval small enough to produce
;  enough spectra.) 
;  By default, spec1 is kinetic energy and spec2 magnetic energy.
;  The routine plots the spectra for the last time saved in the file.
;  All times are stored in the arrays spec1 and spec2.
;  The index of the first time is 1, and of the last time it is i-2.
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
;
;  plot only when iplot=1 (default)
;  can be turned off by using /noplot
;
if keyword_set(noplot) then iplot=0 else iplot=1
;
file1='power'+var1+'.dat'
file2='power'+var2+'.dat'
;
;!p.multi=[0,1,1]
!p.charsize=2

mx=0L & my=0L & mz=0L & nvar=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
;  Reading number of grid points from 'data/dim.dat'
;  Need both mx and nghostx to work out nx.
;  Assume nx=ny=nz
;
datatopdir='data'
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
print,'nx=',nx
imax=nx/2
spectrum1=fltarr(imax)
if n_elements(size) eq 0 then size=2.*!pi
k0=2.*!pi/size
wavenumbers=indgen(imax)*k0 
k=findgen(imax)+1.
k=findgen(imax)
;
;  Looping through all the data to get number of spectral snapshots
;
globalmin=1e12
globalmax=1e-30
i=1
close,1
openr,1, datatopdir+'/'+file1
    while not eof(1) do begin
       readf,1,time 
       readf,1,spectrum1
       if (max(spectrum1(1:*)) gt globalmax) then globalmax=max(spectrum1(1:*))
       if (min(spectrum1(1:*)) lt globalmin) then globalmin=min(spectrum1(1:*))
       i=i+1
    endwhile
close,1
spec1=fltarr(imax,i-1)
tt=fltarr(i-1)
lasti=i-2
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
!x.title='k'
!y.title='P(k)'
!x.range=[1,imax]
;
i=1 
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
       	if (last eq 0) then begin
	  if (time ge tmin) then begin
	    if (time le tmax) then begin
	      xrr=[1,imax]
	      yrr=[globalmin,globalmax]
	      !p.title='t=' + string(time)
	      !y.range=[globalmin,globalmax]
              if iplot eq 1 then begin
		plot_oo,k,spectrum1
         	if (file2 ne '') then begin
		  oplot,k,spectrum2,col=122
		  if (tot eq 1) then begin
		    oplot,k,spectrum1+spectrum2,col=47
		  endif
		endif
	        if (lin ne 0) then begin
		  fac=spectrum1(2)/k(2)^(lin)*1.5
		  oplot,k(2:*),k(2:*)^(lin)*fac,lin=2
		endif
              	wait,w
	      endif
            endif
	  endif
       	endif
       	i=i+1
    endwhile
    if (last eq 1 and iplot eq 1) then begin
	!p.title='t=' + string(time)
	!y.range=[miny,maxy]
	plot_oo,k,spectrum1
      	if (file2 ne '') then begin
		oplot,k,spectrum2,col=122
		if (tot eq 1) then begin
			oplot,k,spectrum1+spectrum2,col=47
		endif
	endif
	if (lin ne 0) then begin
		fac=spectrum1(2)/k(2)^(lin)*1.5
		oplot,k(2:*),k(2:*)^(lin)*fac,lin=2
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





