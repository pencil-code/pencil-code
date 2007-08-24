function div,f
COMPILE_OPT IDL2,HIDDEN
common cdat,x
pc_read_param,obj=par


endif else if (par.coord_system eq 'cartesian') then begin
    return,xder(f[*,*,*,0])+yder(f[*,*,*,1])+zder(f[*,*,*,2])
endif else begin
    print, 'error: div not implemented for spherical polars'
endelse


end
