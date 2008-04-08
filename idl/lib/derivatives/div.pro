function div,f,coord_system
COMPILE_OPT IDL2,HIDDEN
common cdat,x,y
;
default,coord_system,'cartesian'
;
 tmp=size(f) &  my=tmp[2] &  mz=tmp[3]
 lsystem=-1
 if (coord_system eq 'cartesian') then lsystem=0
 if (coord_system eq 'cylindric') then lsystem=1
 if (coord_system eq 'spherical') then lsystem=2
 if (lsystem eq -1) then $
   print,'coord_system= ',coord_system,' is not valid'
;
 if (lsystem eq 0) then corr = 0.
 if (lsystem eq 1) then corr = f[*,*,*,0]/xx
 if (lsystem eq 2) then begin
    yy=spread(y,[1,2],[my,mz])
    cotth=cos(yy)/sin(yy)      
    i_sin=where(abs(sin(yy)) lt 1e-5) ;sinth_min=1e-5
    if (i_sin ne -1) then cotth[i_sin]=0.
    corr = (2.*f[,*,*,0]+cotth*f[*,*,*,1])/xx 
 endif
;
return,xder(f[*,*,*,0])+$
       yder(f[*,*,*,1],lsystem)+ $
       zder(f[*,*,*,2],lsystem)+ corr
end
