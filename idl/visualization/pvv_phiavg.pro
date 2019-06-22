;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pvv_phiavg.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   16-Apr-2004
;;;
;;;  Description:
;;;    Plot azimuthal averages of a vector field. Currently, only 'uu'
;;;    and 'bb' are supported.
;;;  Usage:
;;;    pvv_phiavg [,filename]   [,/BB] [,/UU] [,NLINES=nlines]
;;;    pvv_phiavg [,avg_struct] [,/BB] [,/UU] [,NLINES=nlines]
;;;  IF FILENAME doesn't exist, try 'data/FILENAME' and
;;;  'data/averages/FILENAME'.
;;;  Keywords:
;;;    BB            -- if set, plot magnetic field bb (default is uu)
;;;    UU            -- if set, plot velocity field uu (also the default):
;;;                     poloidal uu as arrows, angular velocity omega as
;;;                     colour
;;;    NEGATE        -- if set, negate fields
;;;    PLOT_UPHI     -- colour-code u_phi, not omega
;;;    PLOT_LZ       -- colour-code l_z = r*u_phi, not omega
;;;    RADII         -- radii for overplotting circles;
;;;                     for statistics, assume [R_sphere] (whole sphere),
;;;                     or [R_in, R_out] (spherical shell)
;;;    MARGIN        -- margin in x and y (see ASPECT_POS)
;;;    AXIS          -- include the z axis in the plot
;;;    STAT          -- print some statistics
;;;    OMEGA0        -- angular velocity Omega0 for getting statistics right
;;;    MAXVEC        -- maximum number of vectors to plot (see WD_VELOVECT)
;;;    NLINES        -- if >0, plot NLINES field lines for the poloidal
;;;                     component instead of arrows
;;;    ABSZERO       -- always use the same color for value zero when
;;;                     NLINES is set
;;;    BOTTOMCOLFRAC -- don't use colors below bottomcolfrac*!d.table_size
;;;    QUIET         -- keep quiet
;;;  Examples:
;;;    pvv_phiavg                    ; plot uu data from PHIAVG1,
;;;    pvv_phiavg, /BB               ; plot bb data from PHIAVG1,
;;;    pvv_phiavg, 'PHIAVG3'         ; plot uu data from PHIAVG3,
;;;    pvv_phiavg, 'PHIAVG3', /BB    ; plot bb data from PHIAVG3,
;;;    pvv_phiavg, avg               ; use data from struct AVG
;;;    pvv_phiavg, avg, /STAT, OMEGA0=.5 ; print some rms values
;;;  Advanced example:
;;;    avg=tavg_phiavg([800.,1e10])
;;;    pvv_phiavg, avg, /BB, RADII=par2.r_ext, /STAT, NLINES=20
;;;    pvv_phiavg, avg, /BB, RADII=par2.r_ext, /STAT, NLINES=20, $
;;;                /ABSZERO, BOTTOMCOLFRAC=0.3, /AXIS

; pro _new_ct, ct
; ;; Build custom colour table, based on colour table 5
;   common _ct_data, orig_ct

;   loadct, 5
  
; end

; pro _restore_ct, ct
; ;; Restore previous colour table
;   common _ct_data, orig_ct

;   loadct, 
  
; end

; ---------------------------------------------------------------------- ;
pro pvv_phiavg, arg, $
                BB=bb, UU=uu, $
                NEGATE=negate, $
                PLOT_UPHI=plot_uphi, PLOT_LZ=plot_lz, $
                RADII=radii, MARGIN=margin, AXIS=axis, $
                STAT=stat, OMEGA0=omega, $
                MAXVEC=maxvec, $
                NLINES=nlines, $
                ABSZERO=abszero, BOTTOMCOLFRAC=bottomcolfrac, $
                QUIET=quiet, HELP=help, $
                CHARSIZE=charsize, $
                XRANGE=xrange, YRANGE=yrange, $
                _EXTRA=extra

  if (keyword_set(help)) then extract_help, 'pvv_phiavg'

  datatopdir ='data'
  avgdir = datatopdir+'/averages'
  phiavgfile = 'PHIAVG1'

  default, arg, avgdir+'/'+phiavgfile
  default, uu, 1
  default, bb, 0
  default, negate, 0
  default, plot_uphi, 0
  default, plot_lz, 0
  default, margin, 0.05
  default, axis, 0
  default, stat, 0
  default, maxvec, [20,40]
  default, nlines, 0
  default, abszero, 0
  default, bottomcolfrac, 0
  default, quiet, 0
  if (n_elements(omega) eq 0) then begin
    omega = 1.
    noomega = 1
  endif else begin
    noomega = 0
  endelse
  default, charsize, !p.charsize

  s = size(arg)
  if (s[s[0]+1] eq 8) then begin
    ;; ARG is struct

    avg = arg

  endif else if(s[s[0]+1] eq 7) then begin
    ;; ARG is string

    ;; Try ARG, data/ARG, data/averages/ARG in that order:
    paths = ['', datatopdir+'/', avgdir+'/']
    found = 0
    filelist = ''
    for i=0,n_elements(paths)-1 do begin
      if (not found) then begin
        file = paths[i]+arg[0]
        filelist = filelist + file + ', '
        if (any(file_test(file))) then found = 1
      endif
    endfor
    if (not found) then message, "Couldn't open any file of "+filelist

    ;; Get data from file:
    if (not quiet) then print, 'Reading '+file

    avg = pc_read_phiavg(file)

  endif else begin

    message, 'ARG must be a struct or a string'

  endelse

  default, xrange, minmax(avg.rcyl)
  default, yrange, minmax(avg.z)
  if (axis) then xrange[0] = 0.

  nr = n_elements(avg.rcyl)
  nz = n_elements(avg.z)

  if (bb) then begin
    var1 = avg.brmphi
    var2 = avg.bzmphi
    var3 = avg.bpmphi
    var3_plot = var3
    varname = 'B'
    Etype = 'mag'
    var3name = 'B_phi'
  endif else if (uu) then begin
    var1 = avg.urmphi
    var2 = avg.uzmphi
    var3 = avg.upmphi
    if (plot_uphi ne 0) then begin
      var3_plot = var3                              ; colour-code u_phi
    endif else if (plot_lz ne 0) then begin
      var3_plot = var3*spread(avg.rcyl,1,nz)        ; colour-code L_z=u_phi*r
    endif else begin
      var3_plot = var3/spread(avg.rcyl>1.e-20,1,nz) ; colour-code omega=u_phi/r
    endelse
    varname = 'v'
    Etype = 'kin'
    var3name = 'Omega'
  endif else begin
    message, 'Need one of UU or BB to be true'
  endelse

  if (negate) then begin
    var1 = -var1
    var2 = -var2
    var3 = -var3
    var3_plot = -var3_plot
  endif

  ;
  pos = aspect_pos(minmax(yrange,/RANGE)/minmax(xrange,/RANGE), $
                   MARGIN=margin)
  ;
  ; Set colour range:
  ;
  data_range = minmax(var3_plot)
  if (abszero) then data_range = [-1,1]*max(abs(data_range))
  if (bottomcolfrac ne 0) then begin
    level_bot = (data_range[0]-bottomcolfrac*data_range[1]) $
                / (1.-bottomcolfrac)
    level_top = data_range[1]
    levels = linspace([level_bot,level_top], 60/(1-bottomcolfrac), GHOST=0.5)
  endif else begin
    levels = linspace(data_range, 60, GHOST=0.5)
  endelse

  if (nlines ne 0) then begin
    ;; Plot fieldlines on color
    ;; 1. Construct vector potential
    pot = var1*0.
    if (quiet eq 0) then print, 'Calculating stream function..'
    for iiz=0,nz-1 do begin
      pot[*,iiz] = integral(avg.rcyl,avg.rcyl*var2[*,iiz],/accumulate)
    endfor
    ;; 2. Plot
    contourfill, var3_plot, $
        avg.rcyl, avg.z, $
        POS=pos, XRANGE=xrange, XSTYLE=1, $
        YRANGE=yrange, YSTYLE=1, $
        XTITLE='!8r!X', YTITLE='!8z!X', $
        LEVELS=levels, $
        CHARSIZE=charsize, $
        _EXTRA=extra
    if (nlines gt 0) then begin
      ; linearly spaced
      levs = linspace(minmax(pot),abs(nlines),GHOST=0.5)
    endif else begin
      ; sqrt-spaced such that uniform field would have equidistant lines
      levs = fllevels(pot,abs(nlines))
    endelse
    contour, pot, avg.rcyl, avg.z, $
        /OVERPLOT, $
        LEVELS=levs
  endif else begin
    ;; Plot arrows on color
    plot_3d_vect, var1, var2, var3_plot, $
        avg.rcyl, avg.z, $
        POS=pos, XRANGE=xrange, XSTYLE=1, $
        YRANGE=yrange, YSTYLE=1, $
        /KEEP, LEVELS=levels, $
        XTITLE='!8r!X', YTITLE='!8z!X', $
        MAXVEC=maxvec, CHARSIZE=charsize, $
        /QUIET, $
        _EXTRA=extra
  endelse

  ;; Overplot spheres (well, circles) if that makes sense
  if (n_elements(radii) gt 0) then opcircle, radii

  ;; Print some statistics. Should ideally be in another place (or
  ;; rather standalone), but it comes handy to do this here.
  if (stat) then begin
    ; ;; weights for use with rms(weight*var):
    ; weight = spread(sqrt(avg.rcyl),1,nz)
    ; weight = weight/rms(weight)
    ;; weights for use with sqrt(mean(weight*var^2)):
    weight = spread(avg.rcyl,1,nz)
    weight = weight/mean(weight)

    ;; weights for averaging over sphere
    ;;   r_spher < radii[0]
    ;; or
    ;;   radii[0] < r_spher < radii[1]:
    ;;
    if (n_elements(radii) gt 0) then begin
      r_spher = sqrt(spread(avg.rcyl,1,nz)^2 + spread(avg.z,0,nr)^2)
      if (n_elements(radii) eq 1) then begin ; just one radius -> outer
        rout = radii[0]
        weight2 = weight*(r_spher lt rout)
        spher = 'sphere  '
      endif else begin          ; at least two radii -> [inner, outer]
        rout = radii[1]
        weight2 = weight*((r_spher gt radii[0]) and (r_spher lt rout))
        spher = 'shell   '
      endelse
      ; weight2 = weight2/rms(weight2)
      weight2 = weight2/mean(weight2)
    endif else begin
      rout = 1.e37
      r_spher = 0.*var1
      weight2 = 0.
      spher = '<ignore>'
    endelse
    ; vrm=rms(weight*var1) & vrm2=rms(weight2*var1)
    ; vzm=rms(weight*var2) & vzm2=rms(weight2*var2)
    ; vpm=rms(weight*var3) & vpm2=rms(weight2*var3)
    vrm=sqrt(mean(weight*var1^2)) & vrm2=sqrt(mean(weight2*var1^2))
    vzm=sqrt(mean(weight*var2^2)) & vzm2=sqrt(mean(weight2*var2^2))
    vpm=sqrt(mean(weight*var3^2)) & vpm2=sqrt(mean(weight2*var3^2))
    vv = sqrt(var1^2+var2^2+var3^2)

    sph = where(r_spher le rout)
    if (sph[0] ge 0) then begin
      vrmin=min(abs(var1[sph])) & vrmax=max(abs(var1))
      vzmin=min(abs(var2[sph])) & vzmax=max(abs(var2))
      vpmin=min(abs(var3[sph])) & vpmax=max(abs(var3))
      vvmin=min(vv[sph])        & vvmax=max(vv)
    endif else begin
      message, /INFO, 'No points satisfying r_spher < rout'
      NaN = !values.f_nan
      vrmin = (vrmax = NaN)
      vzmin = (vzmax = NaN)
      vpmin = (vpmax = NaN)
      vvmin = (vvmax = NaN)
    endelse

    ;; Net differential rotation Omega_max-Omega_min
    om = var3 / spread(avg.rcyl>1.e-20, 1, nz)
    ; Exclude points closest to axis
    dom0 = minmax(om[0:*,*], /RANGE)
    dom1 = minmax(om[1:*,*], /RANGE)
    dom2 = minmax(om[2:*,*], /RANGE)
    dom3 = minmax(om[3:*,*], /RANGE)


    print, 'component:    min         rms_global    rms_'+spher+'   max'
    print, '------------------------------------------------------------------'
    print, varname+'_rcyl   :', vrmin, vrm, vrm2, vrmax
    print, varname+'_phi    :', vpmin, vpm, vpm2, vpmax
    print, varname+'_z      :', vzmin, vzm, vzm2, vzmax
    print, '------------------------------------------------------------------'
    print, varname+varname+'       :', vvmin, sqrt(vrm^2+vpm^2+vzm^2), $
        sqrt(vrm2^2+vpm2^2+vzm2^2), vvmax

    r_cyl = spread(avg.rcyl,1,nz)
    rho = avg.rhomphi
    E_pol = mean(weight2*0.5*rho*(var1^2+var2^2))
    E_tor = mean(weight2*0.5*rho*var3^2)
    E_rot = mean(weight2*0.5*rho*r_cyl^2*Omega^2)
    print
    print, '<E_'+Etype+'>_'+spher+' = '
    print, '  pol ( /Erot):               ', E_pol, E_pol/E_rot
    print, '  tor ( /Erot):               ', E_tor, E_tor/E_rot
    print, 'E_rot ='
    print, '  <1/2*rho r_cyl^2.>_'+spher+' = ', E_rot
    print
    print, 'Delta ' + var3name $
        + ', excluding   0            1            2            3        pts'
    print, '  closest to axis:  ', [dom0,dom1,dom2,dom3]
    print

    if (noomega) then print, 'WARNING: Assuming Omega=1.'

    print, mean(weight2*r_spher^2)
    print, mean(weight2*spread(avg.rcyl,1,nz)^2)
    print, 'total mass: ', mean(weight2*rho)*4./3.*!pi

  endif

end
; End of file pvv_phiavg.pro
