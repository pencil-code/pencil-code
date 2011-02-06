;$Id: pc_read_xyaver_coarsegrain.pro 14056 2010-06-06 09:42:25Z AxelBrandenburg $
;
;  return z array for interior
;
;  30-dec-10/axel
;
pro pc_read_zzz,zzz
;
pc_read_dim,o=dim
pc_read_grid,o=grid
;
zzz=grid.z(dim.n1:dim.n2)
END
