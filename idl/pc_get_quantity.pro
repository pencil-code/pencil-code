;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_get_quantity.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Calculation of physical quantities. If the unit-structure is given,
;;;   this unit system is chosen, otherwise the result is in Pencil-units.
;;;
;;;  Available physical quantities are:
;;;
;;;   Label            Description
;;;  =============================================================
;;;   u_abs            absolute velocity
;;;   uu               velocity
;;;   Temp             temperature
;;;   Spitzer_q        absolute value of Spitzer heat flux vector
;;;   ln_rho           natural logarithm of density
;;;   log_rho          decatic logarithm of density
;;;   rho              density
;;;   P                pressure
;;;   HR_viscous       volumetric viscous heating rate
;;;   B                magnetic field
;;;   HR_ohm           volumetric Ohmic heating rate
;;;   j                current density
;;;


; Calculation of physical quantities.
function pc_get_quantity, vars, index, quantity, unit=unit, dim=dim, grid=grid, param=param, run_param=run_param

  ; First and last physical value, excluding ghost cells
  l1 = dim.l1
  l2 = dim.l2
  m1 = dim.m1
  m2 = dim.m2
  n1 = dim.n1
  n2 = dim.n2
  
  sources = tag_names (index)

  if (n_elements(grid) eq 0) then begin
    print, "WARNING: using unity as unit."
    unit = { velocity:1, temperature:1, density:1, magnetic_field:1, time:1, length:1, mass:1 }
  endif

  ; Absolute velocity
  if (strcmp (quantity, 'u_abs', /fold_case)) then begin
    if (not any (strcmp (sources, 'ux', /fold_case)) or not any (strcmp (sources, 'uy', /fold_case)) or not any (strcmp (sources, 'uz', /fold_case))) then $
        message, "Can't compute '"+quantity+"' without 'ux', 'uy', or 'uz'"
    return, sqrt (dot2 (vars[l1:l2,m1:m2,n1:n2,index.ux:index.uz])) * unit.velocity
  end

  ; Velocity
  if (strcmp (quantity, 'uu', /fold_case)) then begin
    if (not any (strcmp (sources, 'ux', /fold_case)) or not any (strcmp (sources, 'uy', /fold_case)) or not any (strcmp (sources, 'uz', /fold_case))) then $
        message, "Can't compute '"+quantity+"' without 'ux', 'uy', or 'uz'"
    return, vars[l1:l2,m1:m2,n1:n2,index.ux:index.uz] * unit.velocity
  end

  ; Temperature
  if (strcmp (quantity, 'Temp', /fold_case)) then begin
    if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
      return, exp (vars[l1:l2,m1:m2,n1:n2,index.lnTT]) * unit.temperature
    end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
      return, vars[l1:l2,m1:m2,n1:n2,index.TT] * unit.temperature
    end else begin
      message, "Can't compute '"+quantity+"' without 'lnTT' or 'TT'"
    end
  end

  ; Absolute value of the Spitzer heat flux density vector q [W/m^2] = [kg/s^3]
  if (strcmp (quantity, 'Spitzer_q', /fold_case)) then begin
    if (not any (tag_names (run_param) eq "K_SPITZER")) then message, "Can't compute '"+quantity+"' without parameter 'K_SPITZER'"
    Temp = pc_get_quantity (vars, index, 'Temp', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param)
    return, run_param.K_spitzer * Temp^2.5 * sqrt (dot2 ((grad (vars[*,*,*,index.TT]))[l1:l2,m1:m2,n1:n2,*])) * unit.density * unit.velocity^3
  end

  ; Density
  if (strcmp (quantity, 'ln_rho', /fold_case)) then begin
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      return, alog (exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * unit.density)
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      return, alog (vars[l1:l2,m1:m2,n1:n2,index.rho] * unit.density)
    end else begin
      message, "Can't compute '"+quantity+"' without 'lnrho' or 'rho'"
    end
  end else if (strcmp (quantity, 'log_rho', /fold_case)) then begin
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      return, alog10 (exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * unit.density)
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      return, alog10 (vars[l1:l2,m1:m2,n1:n2,index.rho] * unit.density)
    end else begin
      message, "Can't compute '"+quantity+"' without 'lnrho' or 'rho'"
    end
  end else if (strcmp (quantity, 'rho', /fold_case)) then begin
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      return, exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * unit.density
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      return, vars[l1:l2,m1:m2,n1:n2,index.rho] * unit.density
    end else begin
      message, "Can't compute '"+quantity+"' without 'lnrho' or 'rho'"
    end
  end

  ; Pressure
  if (strcmp (quantity, 'P', /fold_case)) then begin
    if (not any (tag_names (run_param) eq "CP") or not any (tag_names (run_param) eq "GAMMA")) then message, "Can't compute '"+quantity+"' without parameter 'CP' or 'GAMMA'"
    rho = pc_get_quantity (vars, index, 'rho', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param)
    Temp = pc_get_quantity (vars, index, 'Temp', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param)
    return, param.cp * (param.gamma - 1.0) / param.gamma * rho * Temp * unit.density * unit.velocity^2
  end

  ; Viscous heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
  if (strcmp (quantity, 'HR_viscous', /fold_case)) then begin
    if (not any (strcmp (sources, 'ux', /fold_case)) or not any (strcmp (sources, 'uy', /fold_case)) or not any (strcmp (sources, 'uz', /fold_case))) then $
        message, "Can't compute '"+quantity+"' without 'ux', 'uy', or 'uz'"
    if (not any (tag_names (run_param) eq "NU")) then message, "Can't compute '"+quantity+"' without parameter 'NU'"
    u_xx = (xder (vars[*,*,*,index.ux]))[l1:l2,m1:m2,n1:n2]
    u_xy = (yder (vars[*,*,*,index.ux]))[l1:l2,m1:m2,n1:n2]
    u_xz = (zder (vars[*,*,*,index.ux]))[l1:l2,m1:m2,n1:n2]
    u_yx = (xder (vars[*,*,*,index.uy]))[l1:l2,m1:m2,n1:n2]
    u_yy = (yder (vars[*,*,*,index.uy]))[l1:l2,m1:m2,n1:n2]
    u_yz = (zder (vars[*,*,*,index.uy]))[l1:l2,m1:m2,n1:n2]
    u_zx = (xder (vars[*,*,*,index.uz]))[l1:l2,m1:m2,n1:n2]
    u_zy = (yder (vars[*,*,*,index.uz]))[l1:l2,m1:m2,n1:n2]
    u_zz = (zder (vars[*,*,*,index.uz]))[l1:l2,m1:m2,n1:n2]
    div_u3 = (u_xx + u_yy + u_zz) / 3.0
    rho = pc_get_quantity (vars, index, 'rho', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param)
    return, run_param.nu * rho * ( 2*((u_xx - div_u3)^2 + (u_yy - div_u3)^2 + (u_zz - div_u3)^2) + (u_xy + u_yx)^2 + (u_xz + u_zx)^2 + (u_yz + u_zy)^2 ) * unit.density * unit.velocity^3 / unit.length
  end

  mu0_SI = 4.0 * !Pi * 1.e-7

  ; Magnetic field
  if (strcmp (quantity, 'B', /fold_case)) then begin
    if (not any (strcmp (sources, 'ax', /fold_case)) or not any (strcmp (sources, 'ay', /fold_case)) or not any (strcmp (sources, 'az', /fold_case))) then $
        message, "Can't compute '"+quantity+"' without 'ax', 'ay', or 'az'"
    return, curl (vars[*,*,*,index.ax:index.az]) * unit.magnetic_field
  end

  ; Magnetic energy density [WORK HERE: unfinished, currently only computes B^2]
  if (strcmp (quantity, 'rho_mag', /fold_case)) then begin
    if (not any (strcmp (sources, 'ax', /fold_case)) or not any (strcmp (sources, 'ay', /fold_case)) or not any (strcmp (sources, 'az', /fold_case))) then $
        message, "Can't compute '"+quantity+"' without 'ax', 'ay', or 'az'"
    return, dot2 (pc_get_quantity (vars, index, 'B', unit=unit, dim=dim, grid=grid, param=param, run_param=run_param))
  end

  ; Ohming heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
  if (strcmp (quantity, 'HR_ohm', /fold_case)) then begin
    if (not any (strcmp (sources, 'ax', /fold_case)) or not any (strcmp (sources, 'ay', /fold_case)) or not any (strcmp (sources, 'az', /fold_case))) then $
        message, "Can't compute '"+quantity+"' without 'ax', 'ay', or 'az'"
    if (not any (tag_names (run_param) eq "ETA")) then message, "Can't compute '"+quantity+"' without parameter 'ETA'"
    return, run_param.eta * param.mu0 * dot2 ((curlcurl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] / param.mu0) * unit.density * unit.velocity^3 / unit.length
  end

  ; Current density [A / m^2]
  if (strcmp (quantity, 'j', /fold_case)) then begin
    if (not any (strcmp (sources, 'ax', /fold_case)) or not any (strcmp (sources, 'ay', /fold_case)) or not any (strcmp (sources, 'az', /fold_case))) then $
        message, "Can't compute '"+quantity+"' without 'ax', 'ay', or 'az'"
    return, sqrt (dot2 ((curlcurl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] / param.mu0)) * unit.velocity * sqrt (param.mu0 / mu0_SI * unit.density) / unit.length
  end

  message, "Unknown quantity '"+quantity+"'"
  return, !Values.D_NaN

end

