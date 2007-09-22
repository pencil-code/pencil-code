function div,f
COMPILE_OPT IDL2,HIDDEN
common cdat,x
pc_read_param,obj=par

if (par.coord_system eq 'cylindric') then begin
    tmp=size(f) &  my=tmp[2] &  mz=tmp[3]
    xx=spread(x,[1,2],[my,mz])
    return,xder(f[*,*,*,0])+yder(f[*,*,*,1])/xx+zder(f[*,*,*,2])+f[*,*,*,0]/xx 
endif else if (par.coord_system eq 'cartesian') then begin
    return,xder(f[*,*,*,0])+yder(f[*,*,*,1])+zder(f[*,*,*,2])
endif else begin
    print, 'error: div not implemented for spherical polars'
    return,0
endelse
end
