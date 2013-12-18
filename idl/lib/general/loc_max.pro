  function loc_max, series
;
; returns indices of local maxima of series
;
    ns = n_elements(series)

    s1 = series(1:ns-1)-series(0:ns-2) 
    s2 = series(1:ns-2)-series(2:ns-1)

    return, where(sign(s1) gt 0. and sign(s2) gt 0.)+1

  end
