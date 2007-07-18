function pc_unshear,f,t
;
;  Transform a quantity back from the sheared frame
;  in order to make it fully periodic.
;
;  16-Jul-2007/tobi: coded
;

  pc_read_param,obj=param,/quiet
  pc_read_grid,obj=grid,/quiet
  pc_read_dim,obj=dim,/quiet
  
  tshear = 1./abs(param.Sshear)
  deltay = -param.Sshear*grid.x*(t - floor(t/tshear)*tshear)
  
  deltay_dy = deltay/grid.dy
  displs = floor(deltay_dy)
  frak = deltay_dy - displs
  
  c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
  c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
  c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
  c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
  c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
  c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
  
  g = f
  for l = 0,dim.mx-1 do begin
    g[l,dim.m1:dim.m2,*,*] = $
        c1[l]*shift(reform(f[l,dim.m1:dim.m2,*,*]),displs[l]-2) $
      + c2[l]*shift(reform(f[l,dim.m1:dim.m2,*,*]),displs[l]-1) $
      + c3[l]*shift(reform(f[l,dim.m1:dim.m2,*,*]),displs[l]  ) $
      + c4[l]*shift(reform(f[l,dim.m1:dim.m2,*,*]),displs[l]+1) $
      + c5[l]*shift(reform(f[l,dim.m1:dim.m2,*,*]),displs[l]+2) $
      + c6[l]*shift(reform(f[l,dim.m1:dim.m2,*,*]),displs[l]+3)
  endfor
  
  g[*,dim.m1-3:dim.m1-1,*,*] = g[*,dim.m2-2:dim.m2,*,*]
  g[*,dim.m2+1:dim.m2+3,*,*] = g[*,dim.m1:dim.m1+2,*,*]
  
  return,g

end
