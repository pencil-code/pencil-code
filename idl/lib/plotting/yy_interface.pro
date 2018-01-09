  pro yy_interface, yprocs=yprocs, zprocs=zprocs, procs=procs, _extra=extra

;  Plots interface between Yin and Yang grid and optionally processor boundaries in Yin grid
;  onto an existing longitude-latitude plot of a full spherical surface.
;  Graphics keywords are conveyed to the plots calls.

    plots, [1.,7.,7.,1.,1.]*!pi/4, [1.,1.,3.,3.,1.]*!pi/4, _extra=extra

    if keyword_set(yprocs) or keyword_set(zprocs) or keyword_set(procs) then begin

      pc_read_dim, obj=d, /quiet
      if keyword_set(procs) then begin yprocs=1 & zprocs=1 & endif
      if keyword_set(yprocs) then $
        for i=1,d.nprocy-1 do $
          plots, [1.,7.]*!pi/4, [1.,1.]*!pi/4 + !pi/2./d.nprocy*i, [1.,3.]*!pi/4, _extra=extra
      if keyword_set(zprocs) then $
        for i=1,d.nprocz-1 do $
          plots, [1.,1.]*!pi/4 + 3*!pi/2./d.nprocz*i, [1.,3.]*!pi/4, _extra=extra

    endif

  end
