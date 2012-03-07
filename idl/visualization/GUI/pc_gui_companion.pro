;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_gui_companion.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Framework for precalculation and comparision of output in pencil units.
;;;   Companion procedures needed by 'pc_gui.pro'.
;;;
;;;  To do:
;;;   Add more comments


; Prepares the varset
pro prepare_varset, num, units, coords, varset, overset, dir, params, run_params

  common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param
  
  datadir = dir
  
  unit = units
  coord = coords
  param = params
  run_param = run_params
  
  varfiles = { title:"-", time:0.0d0, loaded:0, number:-1, precalc_done:0 }
  varfiles = replicate (varfiles, num)
  
  varsets = replicate (varset, num)
  oversets = replicate (overset, num)
end


; Precalculates a data set and loads data, if necessary
pro precalc, i, number=number, varfile=varfile, datadir=dir, dim=dim, param=par, run_param=run_par, varcontent=varcontent, allprocs=allprocs, show_aver=show_aver, time=time
  
  common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param

  ; Default settings
  default, show_aver, 0
  default, number, i
  default, dir, pc_get_datadir()
  default, datadir, dir
  default, time, 0.0d0
  if (keyword_set (par)) then param = par
  if (keyword_set (run_par)) then run_param = run_par
  
  if (varfiles[i].number le 0) then varfiles[i].number = number
  
  if (varfiles[i].loaded eq 0) then begin
    default, varfile, "var.dat"
    if (n_elements (vars) eq 0) then begin
      print, 'Reading: ', varfile, ' ... please wait!'
      if (i eq 0) then begin
        pc_read_var_raw, varfile=varfile, object=vars, tags=tags, datadir=datadir, dim=dim, param=param, par2=run_param, varcontent=varcontent, allprocs=allprocs, time=time
      endif else begin
        pc_read_var_raw, varfile=varfile, object=vars, tags=tags, datadir=datadir, dim=dim, param=param, par2=run_param, varcontent=varcontent, allprocs=allprocs, time=time, /quiet
      endelse
      sources = varcontent.idlvar
      sources = sources[where (varcontent.idlvar ne 'dummy')]
      precalc_data, number, vars, tags
      print, 'Ready.'
    end
    varfiles[i].title = varfile
    varfiles[i].loaded = 1
    varfiles[i].precalc_done = 1
    varfiles[i].time = time * unit.time / unit.default_time
    vars = 0
  end
  
  if (show_aver) then draw_averages, number
  if (keyword_set (par)) then par = param
  if (keyword_set (run_par)) then run_par = run_param
end


; Precalculates a data set
pro precalc_data, i, vars, index

  common varset_common, set, overplot, oversets, unit, coord, varsets, varfiles, datadir, sources, param, run_param
  
  ; First and last physical value, excluding ghost cells
  l1 = coord.l1
  l2 = coord.l2
  m1 = coord.m1
  m2 = coord.m2
  n1 = coord.n1
  n2 = coord.n2
  
  tags = tag_names (varsets[i])

  ; Compute all desired quantities from available source data
  if (any (strcmp (sources, 'uu', /fold_case))) then begin
    if (any (strcmp (tags, 'u_abs', /fold_case))) then begin
      ; Absolute velocity
      varsets[i].u_abs = sqrt (dot2 (vars[l1:l2,m1:m2,n1:n2,index.ux:index.uz])) * unit.velocity / unit.default_velocity
    end
    if (any (strcmp (tags, 'u_x', /fold_case))) then begin
      ; Velocity x-component
      varsets[i].u_x = vars[l1:l2,m1:m2,n1:n2,index.ux] * unit.velocity / unit.default_velocity
    end
    if (any (strcmp (tags, 'u_y', /fold_case))) then begin
      ; Velocity y-component
      varsets[i].u_y = vars[l1:l2,m1:m2,n1:n2,index.uy] * unit.velocity / unit.default_velocity
    end
    if (any (strcmp (tags, 'u_z', /fold_case))) then begin
      ; Velocity z-component
      varsets[i].u_z = vars[l1:l2,m1:m2,n1:n2,index.uz] * unit.velocity / unit.default_velocity
    end
  end
  if (any (strcmp (tags, 'Temp', /fold_case))) then begin
    ; Temperature
    if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
      varsets[i].Temp = exp (vars[l1:l2,m1:m2,n1:n2,index.lnTT]) * unit.temperature
    end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
      varsets[i].Temp = vars[l1:l2,m1:m2,n1:n2,index.TT] * unit.temperature
    end
  end
  if (any (strcmp (tags, 'Spitzer_q', /fold_case)) and any (tag_names (run_param) eq "K_SPITZER")) then begin
    ; Absolute value of the Spitzer heat flux density vector q [W/m^2] = [kg/s^3]
    if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
      varsets[i].Spitzer_q = run_param.K_spitzer * exp (vars[l1:l2,m1:m2,n1:n2,index.lnTT]) ^ 2.5 * sqrt (dot2 ((grad (exp (vars[*,*,*,index.lnTT])))[l1:l2,m1:m2,n1:n2,*])) * unit.density * unit.velocity^3
    end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
      varsets[i].Spitzer_q = run_param.K_spitzer * (vars[l1:l2,m1:m2,n1:n2,index.TT]) ^ 2.5 * sqrt (dot2 ((grad (vars[*,*,*,index.TT]))[l1:l2,m1:m2,n1:n2,*])) * unit.density * unit.velocity^3
    end
  end
  if (any (strcmp (tags, 'ln_rho', /fold_case))) then begin
    ; Natural logarithmic density
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      varsets[i].ln_rho = alog (exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * unit.density / unit.default_density)
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      varsets[i].ln_rho = alog (vars[l1:l2,m1:m2,n1:n2,index.rho] * unit.density / unit.default_density)
    end
  end else if (any (strcmp (tags, 'log_rho', /fold_case))) then begin
    ; Logarithmic density
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      varsets[i].log_rho = alog10 (exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * unit.density / unit.default_density)
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      varsets[i].log_rho = alog10 (vars[l1:l2,m1:m2,n1:n2,index.rho] * unit.density / unit.default_density)
    end
  end else if (any (strcmp (tags, 'rho', /fold_case))) then begin
    ; Density
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      varsets[i].rho = exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * unit.density / unit.default_density
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      varsets[i].rho = vars[l1:l2,m1:m2,n1:n2,index.rho] * unit.density / unit.default_density
    end
  end
  if (any (strcmp (tags, 'P', /fold_case))) then begin
    ; Pressure [N/m^2]
    varsets[i].P = param.cp * (param.gamma - 1.0) / param.gamma * unit.density * unit.velocity^2
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      varsets[i].P *= exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho])
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      varsets[i].P *= vars[l1:l2,m1:m2,n1:n2,index.rho]
    endif
    if (any (strcmp (sources, 'lnTT', /fold_case))) then begin
      varsets[i].P *= exp (vars[l1:l2,m1:m2,n1:n2,index.lnTT])
    end else if (any (strcmp (sources, 'TT', /fold_case))) then begin
      varsets[i].P *= vars[l1:l2,m1:m2,n1:n2,index.TT]
    endif
  end
  if (any (strcmp (tags, 'rho_u_z', /fold_case)) and any (strcmp (sources, 'uu', /fold_case))) then begin
    ; Vertical component of impulse density
    if (any (strcmp (sources, 'lnrho', /fold_case))) then begin
      varsets[i].rho_u_z = exp (vars[l1:l2,m1:m2,n1:n2,index.lnrho]) * vars[l1:l2,m1:m2,n1:n2,index.uz] * unit.density*unit.velocity / (unit.default_density*unit.default_velocity)
    end else if (any (strcmp (sources, 'rho', /fold_case))) then begin
      varsets[i].rho_u_z = vars[l1:l2,m1:m2,n1:n2,index.rho] * vars[l1:l2,m1:m2,n1:n2,index.uz] * unit.density*unit.velocity / (unit.default_density*unit.default_velocity)
    endif
  end
  if (any (strcmp (sources, 'aa', /fold_case))) then begin
    if (any (strcmp (tags, 'Ax', /fold_case))) then begin
      ; Magnetic vector potential x-component
      varsets[i].Ax = vars[l1:l2,m1:m2,n1:n2,index.ax] * unit.magnetic_field
    end
    if (any (strcmp (tags, 'Ay', /fold_case))) then begin
      ; Magnetic vector potential y-component
      varsets[i].Ay = vars[l1:l2,m1:m2,n1:n2,index.ay] * unit.magnetic_field
    end
    if (any (strcmp (tags, 'Az', /fold_case))) then begin
      ; Magnetic vector potential z-component
      varsets[i].Az = vars[l1:l2,m1:m2,n1:n2,index.az] * unit.magnetic_field
    end
    ; Magnetic field
    bb = curl (vars[*,*,*,index.ax:index.az]) * unit.magnetic_field
    if (any (strcmp (tags, 'bx', /fold_case))) then begin
      ; Magnetic field x-component
      varsets[i].bx = reform (bb[l1:l2,m1:m2,n1:n2,0]) / unit.default_magnetic_field
    end
    if (any (strcmp (tags, 'by', /fold_case))) then begin
      ; Magnetic field y-component
      varsets[i].by = reform (bb[l1:l2,m1:m2,n1:n2,1]) / unit.default_magnetic_field
    end
    if (any (strcmp (tags, 'bz', /fold_case))) then begin
      ; Magnetic field z-component
      varsets[i].bz = reform (bb[l1:l2,m1:m2,n1:n2,2]) / unit.default_magnetic_field
    end
    if (any (strcmp (tags, 'rho_mag', /fold_case))) then begin
      ; Magnetic energy density
      varsets[i].rho_mag = dot2 (bb[l1:l2,m1:m2,n1:n2,*])
    end
    if (any (strcmp (tags, 'spitzer_ratio', /fold_case))) then begin
      ; Ratio of perpendicular and parallel Spitzer heat conduction coefficients
      m_p = 1.6726218e-27 ; [kg]
      if (any (strcmp (tags, 'ln_rho', /fold_case))) then begin
        n = exp (varsets[i].ln_rho) / (m_p * param.mu)
      end else if (any (strcmp (tags, 'log_rho', /fold_case))) then begin
        n = 10^varsets[i].log_rho / (m_p * param.mu)
      end else if (any (strcmp (tags, 'rho', /fold_case))) then begin
        n = varsets[i].rho / (m_p * param.mu)
      end
      varsets[i].spitzer_ratio = 2.e-31 * n^2 / (varsets[i].Temp^3 * dot2 (bb[l1:l2,m1:m2,n1:n2,*]))
      n = 0
    end
    mu0_SI = 4.0 * !Pi * 1.e-7
    if (any (strcmp (tags, 'HR_ohm', /fold_case)) and any (tag_names (run_param) eq "ETA")) then begin
      ; Ohming heating rate [W / m^3] = [kg/m^3] * [m/s]^3 / [m]
      varsets[i].HR_ohm = run_param.eta * param.mu0 * dot2 ((curlcurl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] / param.mu0) * unit.density * unit.velocity^3 / unit.length
    end
    if (any (strcmp (tags, 'j', /fold_case))) then begin
      ; Current density [A / m^2]
      if (any (strcmp (tags, 'HR_ohm', /fold_case)) and any (tag_names (run_param) eq "ETA")) then begin
        varsets[i].j = sqrt (varsets[i].HR_ohm / (run_param.eta * mu0_SI * unit.velocity * unit.length))
      end else begin
        varsets[i].j = sqrt (dot2 ((curlcurl (vars[*,*,*,index.ax:index.az]))[l1:l2,m1:m2,n1:n2,*] / param.mu0)) * unit.velocity * sqrt (param.mu0 / mu0_SI * unit.density) / unit.length
      end
    end
  end
  
  over_tags = tag_names (oversets[i])
  if (any (strcmp (sources, 'uu', /fold_case))) then begin
    if (any (strcmp (over_tags, 'u', /fold_case))) then begin
      ; Velocity overplot
      oversets[i].u[*,*,*,0] = float (vars[l1:l2,m1:m2,n1:n2,index.ux] * unit.velocity / unit.default_velocity)
      oversets[i].u[*,*,*,1] = float (vars[l1:l2,m1:m2,n1:n2,index.uy] * unit.velocity / unit.default_velocity)
      oversets[i].u[*,*,*,2] = float (vars[l1:l2,m1:m2,n1:n2,index.uz] * unit.velocity / unit.default_velocity)
    end
  end
  if (any (strcmp (tags, 'Rn_mag', /fold_case)) and any (strcmp (tags, 'bx', /fold_case)) and any (strcmp (tags, 'by', /fold_case)) and any (strcmp (tags, 'bz', /fold_case)) and any (strcmp (sources, 'uu', /fold_case)) and any (tag_names (run_param) eq "ETA")) then begin
    ; Magnetic mesh Reynolds number
    bb_abs_1 = 1.0 / (sqrt (dot2 (bb[l1:l2,m1:m2,n1:n2,*])) / unit.default_magnetic_field)
    Rx = reform (coord.dx, coord.nx, 1, 1) * abs (vars[l1:l2,m1:m2,n1:n2,index.ux]) * (1 - abs (varsets[i].bx * bb_abs_1))
    Ry = reform (coord.dy, 1, coord.ny, 1) * abs (vars[l1:l2,m1:m2,n1:n2,index.uy]) * (1 - abs (varsets[i].by * bb_abs_1))
    Rz = reform (coord.dz, 1, 1, coord.nz) * abs (vars[l1:l2,m1:m2,n1:n2,index.uz]) * (1 - abs (varsets[i].bz * bb_abs_1))
    varsets[i].Rn_mag = ((Rx > Ry) > Rz) / run_param.eta / unit.length
    Rx = 0
    Ry = 0
    Rz = 0
  end
  if (any (strcmp (sources, 'aa', /fold_case))) then begin
    if (any (strcmp (over_tags, 'b', /fold_case))) then begin
      ; Magnetic field overplot
      oversets[i].b[*,*,*,0] = float (bb[l1:l2,m1:m2,n1:n2,0] / unit.default_magnetic_field)
      oversets[i].b[*,*,*,1] = float (bb[l1:l2,m1:m2,n1:n2,1] / unit.default_magnetic_field)
      oversets[i].b[*,*,*,2] = float (bb[l1:l2,m1:m2,n1:n2,2] / unit.default_magnetic_field)
    end
    bb = 0
    if (any (strcmp (over_tags, 'a_contour', /fold_case))) then begin
      ; Magnetic field lines overplot
      oversets[i].a_contour[*,*,*,0] = float (vars[l1:l2,m1:m2,n1:n2,index.ax] * unit.magnetic_field)
      oversets[i].a_contour[*,*,*,1] = float (vars[l1:l2,m1:m2,n1:n2,index.ay] * unit.magnetic_field)
      oversets[i].a_contour[*,*,*,2] = float (vars[l1:l2,m1:m2,n1:n2,index.az] * unit.magnetic_field)
    end
  end
  if (any (strcmp (over_tags, 'grad_P', /fold_case))) then begin
    ; Gradient of pressure
    grad_P_fact = param.cp * (param.gamma - 1.0) / param.gamma * unit.density*unit.temperature/unit.length
    if (any (strcmp (sources, 'lnrho', /fold_case)) and any (strcmp (sources, 'lnTT', /fold_case))) then begin
      varsets[i].grad_P = float (grad_P_fact * (grad (exp (vars[*,*,*,index.lnrho])) * exp (vars[*,*,*,index.lnTT]) + exp (vars[*,*,*,index.lnrho]) * grad (exp (vars[*,*,*,index.lnTT])))[l1:l2,m1:m2,n1:n2])
    end else if (any (strcmp (sources, 'rho', /fold_case)) and any (strcmp (sources, 'lnTT', /fold_case))) then begin
      varsets[i].grad_P = float (grad_P_fact * (grad (vars[*,*,*,index.rho]) * exp (vars[*,*,*,index.lnTT]) + vars[*,*,*,index.rho] * grad (exp (vars[*,*,*,index.lnTT])))[l1:l2,m1:m2,n1:n2])
    end else if (any (strcmp (sources, 'lnrho', /fold_case)) and any (strcmp (sources, 'TT', /fold_case))) then begin
      varsets[i].grad_P = float (grad_P_fact * (grad (exp (vars[*,*,*,index.lnrho])) * vars[*,*,*,index.TT] + exp (vars[*,*,*,index.lnrho]) * grad (vars[*,*,*,index.TT]))[l1:l2,m1:m2,n1:n2])
    end else if (any (strcmp (sources, 'rho', /fold_case)) and any (strcmp (sources, 'TT', /fold_case))) then begin
      varsets[i].grad_P = float (grad_P_fact * (grad (vars[*,*,*,index.rho]) * vars[*,*,*,index.TT] + vars[*,*,*,index.rho] * grad (vars[*,*,*,index.TT]))[l1:l2,m1:m2,n1:n2])
    endif
  end
  if (any (strcmp (tags, 'rho_c', /fold_case))) then begin
    ; Minimum density for an Alfvén speed below the speed of light
    varsets[i].rho_c = varsets[i].log_rho - alog10 (dot2 (bb[l1:l2,m1:m2,n1:n2,*]) / (2 * mu0_SI * (299792458.0 * run_param.cdtv)^2))
  end
end


; Dummy routine
pro pc_gui_companion

  pc_gui_companion_loaded = 1
end

