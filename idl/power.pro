PRO power,file1,file2,k=k,spec1=spec1,spec2=spec2,i=i

default,file1,'poweru.dat'

!p.multi=[0,1,1]
!p.charsize=2

mx=0L & my=0L & mz=0L & nvar=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
;  Reading number of grid points from 'data/dim.dat'
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
i=1
close,1
openr,1, datadir+'/'+file1
    while not eof(1) do begin 
       readf,1,spectrum
       i=i+1
    endwhile
close,1
spec1=fltarr(imax,i-1)
;
;  Opening file 2 if it is defined
;
if (file2 ne '') then begin
  close,2
  spectrum2=fltarr(imax)
  openr,2,datadir+'/'+file2
  spec2=fltarr(imax,i-1)
endif
;
;  Plotting the results
;
i=1 
openr,1, datadir+'/'+file1
    while not eof(1) do begin 
       readf,1,spectrum1
       spec1(*,i-1)=spectrum1
       maxy=max(spectrum1)
       miny=min(spectrum1)
       if (file2 ne '') then begin
	   readf,2,spectrum2
           spec2(*,i-1)=spectrum2
           if (max(spectrum2) gt maxy) then maxy=max(spectrum2)
           if (min(spectrum2) lt miny) then miny=min(spectrum2)
       endif
       plot_oo,xrange=[1,imax],yrange=[miny,maxy],k,spectrum1
       if (file2 ne '') then oplot,k,spectrum2,col=122
       wait,.1
       i=i+1
    endwhile
close,1
close,2



END
