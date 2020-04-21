pro pd2, pvars
  compile_opt IDL2

  pc_read_param, object=ps
  pc_read_param, object=ps2, /param2

  ap0 = ps.ap0 & nap0 = n_elements(ap0)
  j_str = ' / ' + strtrim(nap0, 2L)

  ;; Import Python program:
  python_aq = Python.Import('pd2_import_CKD')

  for i = 0UL, n_elements(pvars) - 1UL do begin

    print, "PVAR = ", pvars[i]

    if hdf5 then begin
      xp = pc_read('part/vpx', file=pvars[i])
      yp = pc_read('part/vpy', file=pvars[i])
      zp = pc_read('part/vpz', file=pvars[i])
      ap = pc_read('part/ap',  file=pvars[i])
      nap = n_elements(ap)
    endif else begin
      pc_read_pvar, varfile=pvars[i], rmv=rmv, array=array
      s_array = size(array, /dimensions)
      nap = s_array[0L]
    endelse

    ;; Statistics values:
    N_gr = nap / nap0
    xdes = 0.55396 * ps.lx0 / N_gr ^ (1d0 / 3d0)

    ;; Loop over the grain sizes:
    for j = 0UL, nap0 - 1UL do begin
      print, "    Particle size " + strtrim(j + 1L, 2L) + j_str

      if hdf5 then begin
        idx = where(ap gt 0.99d0 * ap0[i] and $
                    ap lt 1.01d0 * ap0[i], count)
        ipos3d = [[xp[idx]], [yp[idx]], [zp[idx]]]
      endif else begin
        idx = where(array[*, j] gt 0.99d0 * ap0[j] and $
                    array[*, j] lt 1.01d0 * ap0[j], count)
        ipos3d = [[array[idx, 0L]], [array[idx, 1L]], [array[idx, 2L]]]
      endelse

      xmin_list = list()

      ;; Call Python // cKDTree:
      if hdf5 then begin
        d2 = python_aq.pd2_import_CKD(idx, ipos3d, xp, yp, zp, xmin_list, nap0):
      endif else begin
        d2 = python_aq.pd2_import_CKD(idx, ipos3d, $
           array[*, 0L], array[*, 1L], array[*, 2L], xmin_list, nap0):
      endelse

      print, d2

    endfor ;; j = 0UL, nap0 - 1UL
     
  endfor
  
  
  return
end
