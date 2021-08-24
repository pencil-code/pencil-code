;
;  Shift plane in the x-direction.
;
pro pc_read_aver_shift_plane, nvar_all, array, xshift, par, timefix, ts, t, t0, grid
;
  for ivar=0, nvar_all-1 do $
    array[*,*,ivar] = shift(array[*,*,ivar],xshift,0)
;
;  With shear we need to displace part of the plane.
;
  if (par.Sshear ne 0.0) then begin
;
;  Some z-averages are erroneously calculated together with time series
;  diagnostics. The proper time is thus found in time_series.dat.
;
    if (timefix) then begin
      ii = where( abs(ts.t-t) eq min(abs(ts.t-t)) )
      ii = ii[0]
      if (ts.t[ii] ge (t-ts.dt[ii])) then ii = ii-1
      if (debug) then print, 'it, t, t_ts, dt_ts=', it, t, ts.t[ii], ts.dt[ii]
      deltay = (-par.Sshear*ts.t[ii]*par.Lxyz[0]) mod par.Lxyz[1]
    endif else $
      deltay = (-par.Sshear*t*par.Lxyz[0]) mod par.Lxyz[1]
    
    deltay_int = fix(deltay/grid.dy)
    if (debug) then print, 'it, t, deltay, deltay_int, deltay_frac=', it, t, deltay, deltay_int, deltay/grid.dy-deltay_int
    for ivar=0, nvar_all-1 do begin
      array2 = array[*,*,ivar]
      for ix=0, xshift-1 do begin
        array2[ix,*] = shift(reform(array2[ix,*]),deltay_int)
        array2[ix,*] = pc_shift_6th(reform(array2[ix,*]),grid.y,deltay-deltay_int*grid.dy)
      endfor
;
;  Comove with central grid point.
;
      print, 'AJ/2010-05-26: the following appears to be wrong'
      print, '               xax[0] should be xshift*dx'   ; xax,yax anyway not defined
      stop
      for ix=0, nxg-1 do $
        array2[ix,*] = pc_shift_6th(reform(array2[ix,*]),yax,-par.Sshear*(t-t0)*xax[0])
      
      array[*,*,ivar] = array2
    endfor
  endif
;
end
