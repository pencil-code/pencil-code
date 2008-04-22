function div,f
COMPILE_OPT IDL2,HIDDEN
common cdat,x,y
common cdat_coords,coord_system
;
 tmp=size(f) &  my=tmp[2] &  mz=tmp[3]
;
 if (coord_system eq 'cartesian') then corr = 0.
;
 if (coord_system eq 'cylindric') then begin
     xx=spread(x,[1,2],[my,mz])
     corr = f[*,*,*,0]/xx
 endif
;
 if (coord_system eq 'spherical') then begin
    yy=spread(y,[1,2],[my,mz])
    cotth=cos(yy)/sin(yy)      
    i_sin=where(abs(sin(yy)) lt 1e-5) ;sinth_min=1e-5
    if (i_sin ne -1) then cotth[i_sin]=0.
    corr = (2.*f[*,*,*,0]+cotth*f[*,*,*,1])/xx 
 endif
;
return,xder(f[*,*,*,0])+ $
       yder(f[*,*,*,1])+ $
       zder(f[*,*,*,2])+ corr
end
