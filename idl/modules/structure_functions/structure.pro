;  $Idg: structure.pro,v 1.15 2003/09/29 08:07:19 nilshau Exp $
;  procedure to read structure functions
;
;  17-dec-02/nils: coded
;  12-may-03/wolf: fixed param.pro stuff 
;
;  variable  : Which variable to consider (sfu, sfb, sfz1, sfz2)
;  moment1   : The moment of the struct. func. (often denoted q or p)
;  moment2   : The moment to plot against if zeta is set (moment2 is
;              normally 3)
;  mplot     : Plotting mode:
;              -average : plot the average between t1_struct and t2_struct
;              -all     : plot all snapshots bet. t1_struct and t2_struct
;              -last    : plot last snapshot only
;  w         : Time to wait between each snapshot if mplot=all
;  derivate  : When derivate='T' is set the derivative of the function
;              instead of the function itself is considered
;  transverse: Longitudinal or transversal:
;              -transverse='1' => longitudinal 
;              -transverse='2' => transversal
;  directions: must tell if the SFs have been calculated in all dirs
;              (directions=3)
;              or only in x-dir (directions=1), which is most common
;  idir      : if idir=0 transvers=2 give the sum of the two
;              directions, while transverse=3 is empty
;  xd        : returns the x-values
;  yd        : returns the y-values
;  noplot    : do not plot the result if set
;  zeta      : calculates the structure function scaling exponents if
;              set (the best result is obtained if mplot='all')
;  cut1      : ignores cut1 number of points on the left side
;  cut2      : ignores cut2 number of points on the right side
;  noprint   : do not print results if set
;
; For example; to obtain the structure function scaling exponents of
; the longitudinal magnetic field you type:
; structure,derivate='T',transverse='1',/zeta,mplot='all',idir=0, $
;           w=0,directions=1,/noplot,/noprint,variable='sfb'
; and the zero'th to eight' scaling exponent is found in the file
; 'structb1.sav'. Setting transverse='2' in the call to structure 
; you'll find the sum of the 
; two transversal directions in the file 'structb2.sav'.
;
PRO structure,variable=variable,moment1=moment1,moment2=moment2,w=w, $
	derivate=derivate,mplot=mplot,transverse=transverse, $
	directions=directions,idir=idir,xd=xd,yd=yd,noplot=noplot, $
	zeta=zeta,cut1=cut1,cut2=cut2,noprint=noprint, $
        sf1=sf1,sf2=sf2,sf3=sf3,nsnap=nsnap,tt=tt,t1_struct=t1_struct, $
        t2_struct=t2_struct

;
;; Cannot work (for several reasons):
; spawn,'if (! -f param.pro) touch param.pro'
; @param
;; may work:
  dummy = findfile('param.pro',count=cpar)
if (cpar gt 0) then begin
  print, 'Reading param.pro..'
  spawn, 'cat param.pro', param_code
  for i=0,n_elements(param_code)-1 do begin
    if (execute(param_code[i]) ne 1) then $
        message, 'Cannot execute param.pro (line ' + strtrim(i+1,2) + ')'
  endfor
endif

datatopdir='data'
datadir='proc0'
close,1
;; openr,1,datatopdir+'/'+datadir+'/'+'dim.dat'
openr,1,datatopdir+'/'+'dim.dat'
readf,1,mx,my,mz,nvar
close,1
size=2*!pi
;
;  Calculating some variables
;
nghostx=3
nx=fix(mx-nghostx*2)
dx=size/nx
imax=nint(alog(nx)/alog(2))*2-2
print,'nx,imax,dx=',nx,imax,dx
dir=3
if keyword_set(noplot) then iplot=0 else iplot=1
;
;  Some default values
;
default,directions,1
default,variable,'sfu'
default,moment1,1
default,moment2,3
default,idir,0
default,w,0.1
default,derivate,'F'
default,mplot,'average'
default,transverse,'1'
default,t1_struct,0
default,t2_struct,10000
default,cut1,0
default,cut2,0
;
file=variable+'-'
;
;  Reading in qmax and time
;
file1=file+'1.dat'
openr,1, datatopdir+'/'+file1
   readf,1,time,qmax
close,1
;
;  Printouts
;
print,'qmax=',qmax
print,'time=',time
print,'directions=',directions
print,'variable=',variable
print,''
print,'Variables related to plotting:'
print,'moment1=',moment1
if (derivate eq 'ESS') then begin
  print,'moment2=',moment2
end
print,'transverse=',transverse
print,'idir=',idir
print,'w=',w
print,'derivate=',derivate
print,'mplot=',mplot
;
;  Defining temperary arrays
;
temp1=fltarr(imax,qmax,dir)
temp2=fltarr(imax,qmax,dir)
temp3=fltarr(imax,qmax,dir)
;
;  Determining number of snapshots
;
nsnap=0
file1=file+'1.dat'
openr,1, datatopdir+'/'+file1
  while not eof(1) do begin
    nsnap=nsnap+1
    readf,1,time,qmax
    readf,1,temp1
  endwhile
close,1
;
;  Defining final arrays
;
sf1=fltarr(imax,qmax+1,dir+1,nsnap)
sf2=fltarr(imax,qmax+1,dir+1,nsnap)
sf3=fltarr(imax,qmax+1,dir+1,nsnap)
sf1_aver=fltarr(imax,qmax+1)
sf2_aver=fltarr(imax,qmax+1)
sf3_aver=fltarr(imax,qmax+1)
sf_plot=fltarr(imax)
LL=fltarr(imax)
tt=fltarr(nsnap)
;
;  Finding LL (the separation)
;
for lb_ll=1,imax do begin
  exp2=lb_ll mod 2
  if (lb_ll eq 1) then exp2=0
  exp1=fix((lb_ll)/2)-exp2
  separation=(2^exp1)*(3^exp2)
  LL(lb_ll-1)=separation*dx
end
;
;  Reading data files for structure functions
;
file1=file+'1.dat'
isnap=0
openr,1, datatopdir+'/'+file1
  while not eof(1) do begin
    readf,1,time,qmax
    readf,1,temp1
    sf1(*,1:qmax,1:3,isnap)=temp1
    tt(isnap)=time
    isnap=isnap+1
  endwhile
close,1
;
file1=file+'2.dat'
isnap=0
openr,1, datatopdir+'/'+file1
  while not eof(1) do begin
    readf,1,time,qmax
    readf,1,temp2
    sf2(*,1:qmax,1:3,isnap)=temp2
    isnap=isnap+1
  endwhile
close,1
;
file1=file+'3.dat'
isnap=0
openr,1, datatopdir+'/'+file1
  while not eof(1) do begin
    readf,1,time,qmax
    readf,1,temp3
    sf3(*,1:qmax,1:3,isnap)=temp3
    isnap=isnap+1
  endwhile
close,1
;
print,'The data set contains data from t=',tt(0),' to t=',tt(nsnap-1)
if mplot eq 'average' then begin
     mintime=max([t1_struct,tt(0)])
     maxtime=min([t2_struct,tt(nsnap-1)])
     print,'Averaged between t='+string(mintime)+' and t= '+string(maxtime)
end
;
;  Averageing over different directions if nr_directions=3
;
if directions eq 1 then begin ; No averageing if only one direction is calc.
  sf1(*,1:qmax,0,*)=sf1(*,1:qmax,1,*)
  sf2(*,1:qmax,0,*)=(sf2(*,1:qmax,1,*)+sf3(*,1:qmax,1,*))/2
  sf3(*,1:qmax,0,*)=0
endif else begin
  ;
  ; Must average between sf1, sf2 and sf3 here because sf1(*,1:qmax,1,*),
  ; sf2(*,1:qmax,2,*) and sf3(*,1:qmax,3,*) are all longitudinal. The 
  ; longitudinal average is placed in sf1(*,1:qmax,0,*). 
  ; The transverse averages are similarily put in sf2(*,1:qmax,0,*).
  ; sf3(*,1:qmax,0,*) is left empty.
  ;
  sf1(*,1:qmax,0,*)=(sf1(*,1:qmax,1,*)+sf2(*,1:qmax,2,*)+sf3(*,1:qmax,3,*))/3
  sf2(*,1:qmax,0,*)=(sf1(*,1:qmax,2,*)+sf1(*,1:qmax,3,*)+ $
  	sf2(*,1:qmax,1,*)+sf2(*,1:qmax,3,*)+ $
	sf3(*,1:qmax,1,*)+sf3(*,1:qmax,2,*))/6
end
;
;  Averaging in time (t1_struct and t2_struct is set in param.pro) 
;
if mplot eq 'average' then begin
  good=where(tt ge t1_struct and tt le t2_struct)
  for qq=1,qmax do begin
    for ii=0,imax-1 do begin
      sf1_aver(ii,qq)=mean(sf1(ii,qq,idir,good))
      sf2_aver(ii,qq)=mean(sf2(ii,qq,idir,good))
      sf3_aver(ii,qq)=mean(sf3(ii,qq,idir,good))
    end
  end
end
;
;  Plotting
;
if (mplot eq 'last') then begin
  istart=nsnap-1
  istop =nsnap-1
end
if (mplot eq 'all') then begin
  print,'Printing data between t=',t1_struct,' and t=',t2_struct 
  good=where(tt ge t1_struct and tt le t2_struct)
  istart=good(0)
  gooddim=size(good)
  istop=gooddim(1)+istart-1
end
if (mplot eq 'average') then begin
  istart=nsnap-1
  istop =nsnap-1
end
;
;  Calculate zeta_p if /zeta is set
;
if keyword_set(zeta) then begin
  zeta=fltarr(qmax+1)
  zeta_all_aver=fltarr(qmax+1)
  izeta_min=1
  izeta_max=qmax
  if (qmax eq 9) then izeta_max=8
  derivate='ESS'
endif else begin
  izeta_min=moment1
  izeta_max=moment1
end
;
for mom=izeta_min,izeta_max do begin
  ;
  ;  Starting loop for plotting
  ;
  for i=istart,istop do begin
    ;
    ;  Plot average or not?
    ;
    if (mplot eq 'average') then begin
      ;
      ;  Transverse or longitudinal sf's
      ;
      ;print,'Plotting average between t=',t1_struct,' and t=',t2_struct
      if transverse eq '1' then begin
        sf_plot=sf1_aver(*,mom)
        if (derivate eq 'ESS') then sf_plot2=sf1_aver(*,moment2)
	if (mom eq 9) then begin
	  pos=where(sf1_aver(*,mom) gt 0)
	  sf_plot=sf1_aver(pos,mom)
	end
      end
      if transverse eq '2' then begin
        sf_plot=sf2_aver(*,mom)
        if (derivate eq 'ESS') then sf_plot2=sf2_aver(*,moment2)
	if (mom eq 9) then begin
	  pos=where(sf2_aver(*,mom) gt 0)
	  sf_plot=sf2_aver(pos,mom)
	end
      end
      if transverse eq '3' then begin
        sf_plot=sf3_aver(*,mom)
        if (derivate eq 'ESS') then sf_plot2=sf3_aver(*,moment2)
	if (mom eq 9) then begin
	  pos=where(sf3_aver(*,mom) gt 0)
	  sf_plot=sf3_aver(pos,mom)
	end
      end
    endif else begin
      ;
      ;  Transverse or longitudinal sf's
      ;
      if (not keyword_set(noprint)) then print,'time=',tt(i)
      if transverse eq '1' then begin
        sf_plot=sf1(*,mom,idir,i)
        if (derivate eq 'ESS') then sf_plot2=sf1(*,moment2,idir,i)
	if (mom eq 9) then begin
	  pos=where(sf1(*,mom,idir,i) gt 0)
	  sf_plot=sf1(pos,mom,idir,i)
	end
      end
      if transverse eq '2' then begin
        sf_plot=sf2(*,mom,idir,i)
        if (derivate eq 'ESS') then sf_plot2=sf2(*,moment2,idir,i)
	if (mom eq 9) then begin
	  pos=where(sf2(*,mom,idir,i) gt 0)
	  sf_plot=sf2(pos,mom,idir,i)
	end
      end
      if transverse eq '3' then begin
        sf_plot=sf3(*,mom,idir,i)
        if (derivate eq 'ESS') then sf_plot2=sf3(*,moment2,idir,i)
	if (mom eq 9) then begin
	  pos=where(sf3(*,mom,idir,i) gt 0)
	  sf_plot=sf3(pos,mom,idir,i)
	end
      end
    end
    ;
    ;  Plot derivate, plain or extended self similarity (ESS) 
    ;
    xd=alog10(LL)
    if (mom eq 9) then xd=alog10(LL(pos))
    if (derivate eq 'T') then begin
      yd=alog10(sf_plot)
      yd=deriv(xd,yd)
      if (iplot eq 1) then plot,xd,yd
      if (iplot eq 1) then oplot,xd,xd^0,lin=2
    end
    if (derivate eq 'F') then begin
      if (min(sf_plot) le 0.) then begin
        print,'WARNING: min(sf_plot) le 0'      
      endif else begin
        yd=alog10(sf_plot)
        if (iplot eq 1) then plot,xd,yd
        if (iplot eq 1) then oplot,xd,yd,ps=1,col=122
      end
    end
    if (derivate eq 'ESS') then begin
      if (min(sf_plot) le 0.) then begin
        print,'WARNING: min(sf_plot) le 0'      
      endif else begin
        xd=alog10(sf_plot2)
        yd=alog10(sf_plot)
        if (iplot eq 1) then plot,xd,yd
      end
    end
    wait,w
    ;
    ;  Calculate zeta_p if /zeta is set
    ;
    if (keyword_set(zeta) and mom le 8) then begin
      xd1=cut1
      xd2=imax-1-cut2
      xy=poly_fit(xd(xd1:xd2),yd(xd1:xd2),1)
      if (iplot eq 1) then oplot,xd(xd1:xd2),yd(xd1:xd2),ps=2,col=122
      if (iplot eq 1) then oplot,xd(xd1:xd2),xy(0)+xy(1)*xd(xd1:xd2),lin=2,col=122
      zeta(mom)=xy(1)
      if (not keyword_set(noprint)) then print,'p=',mom,' zeta=',zeta(mom)
      if (mplot eq 'all')  then zeta_all_aver(mom)=zeta_all_aver(mom)+xy(1)
      wait,w
    end
  end
  if (mplot eq 'all') and (keyword_set(zeta)) then begin
    zeta_all_aver(mom)=zeta_all_aver(mom)/(istop-istart+1)
  end
end
if (mplot eq 'all') and (keyword_set(zeta)) then begin
  zeta=zeta_all_aver
  print,'zeta=',zeta
end
;
;  Save zeta_p to file
;
if (variable eq 'sfz1') or (variable eq 'sfz2') then begin
  if keyword_set(zeta) then begin
    if (transverse eq '1') then save,file='struct1.sav',zeta
    if (transverse eq '2') then save,file='struct2.sav',zeta
    if (transverse eq '3') then save,file='struct3.sav',zeta
  end
end
if (variable eq 'sfu') then begin
  if keyword_set(zeta) then begin
    if (transverse eq '1') then save,file='structu1.sav',zeta
    if (transverse eq '2') then save,file='structu2.sav',zeta
    if (transverse eq '3') then save,file='structu3.sav',zeta
  end
end
if (variable eq 'sfb') then begin
  if keyword_set(zeta) then begin
    if (transverse eq '1') then save,file='structb1.sav',zeta
    if (transverse eq '2') then save,file='structb2.sav',zeta
    if (transverse eq '3') then save,file='structb3.sav',zeta
  end
end
;
;
;
END

