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
;;;    BB     -- if set, plot magnetic field bb (default is uu)
;;;    UU     -- if set, plot velocity field uu (also the default)
;;;    RADII  -- radii for overplotting circles
;;;    MARGIN -- margin in x and y (see ASPECT_POS)
;;;    STAT   -- print some statistics
;;;    OMEGA  -- angular velocity for getting statistics right
;;;    MAXVEC -- maximum number of vectors to plot (see WD_VELOVECT)
;;;    NLINES -- if >0, plot NLINES field lines for the poloidal
;;;              component instead of arrows
;;;    QUIET  -- keep quiet
;;;  Examples:
;;;    pvv_phiavg                    ; plot uu data from PHIAVG1,
;;;    pvv_phiavg, /BB               ; plot bb data from PHIAVG1,
;;;    pvv_phiavg, 'PHIAVG3'         ; plot uu data from PHIAVG3,
;;;    pvv_phiavg, 'PHIAVG3', /BB    ; plot bb data from PHIAVG3,
;;;    pvv_phiavg, avg               ; use data from struct AVG
;;;    pvv_phiavg, avg, /STAT, OMEGA=.5 ; print some rms values
;;;  Advanced example:
;;;    avg=tavg_phiavg([800.,1e10])
;;;    pvv_phiavg, avg, /BB, RADII=par2.r_ext, /STAT, NLINES=20

pro pvv_phiavg, arg, BB=bb, UU=uu, $
                MARGIN=margin, RADII=radii, OMEGA=omega, $
                CHARSIZE=charsize, MAXVEC=maxvec, $
                STAT=stat, NLINES=nlines, $
                QUIET=quiet

  datatopdir ='data'
  avgdir = datatopdir+'/averages'
  phiavgfile = 'PHIAVG1'

  default, arg, avgdir+'/'+phiavgfile
  default, uu, 1
  default, bb, 0
  default, margin, 0.05
  default, charsize, !p.charsize
  default, maxvec, [20,40]
  default, stat, 0
  default, quiet, 0
  default, nlines, 0
  if (n_elements(omega) eq 0) then begin
    omega = 1.
    noomega = 1
  endif else begin
    noomega = 0
  endelse


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
        if (any(findfile(file) ne '')) then found = 1
      endif
    endfor
    if (not found) then message, "Couldn't open any file of "+filelist

    ;; Get data from file:
    if (not quiet) then print, 'Reading '+file

    avg = read_phiavg(file)

  endif else begin

    message, 'ARG must be a struct or a string'

  endelse

  nr = n_elements(avg.rcyl)
  nz = n_elements(avg.z)

  if (bb) then begin
    var1 = avg.brmphi
    var2 = avg.bzmphi
    var3 = avg.bpmphi
    var3_plot = var3
  endif else if (uu) then begin
    var1 = avg.urmphi
    var2 = avg.uzmphi
    var3 = avg.upmphi
    var3_plot = var3/spread(avg.rcyl>1.e-20,1,nz) ; plot omega, not u_phi
  endif else begin
    message, 'Need one of UU or BB to be true'
  endelse

  pos = aspect_pos((max(avg.z)-min(avg.z))/(max(avg.rcyl)-min(avg.rcyl)), $
                   MARGIN=margin)
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
        POS=pos, XRANGE=minmax(avg.rcyl), XSTYLE=1, $
        YRANGE=minmax(avg.z), YSTYLE=1, $
        XTITLE='!8r!X', YTITLE='!8z!X', $
        CHARSIZE=charsize
    if (nlines gt 0) then begin
      ; linearly spaced
      levels = linspace(minmax(pot),abs(nlines),GHOST=0.5)
    endif else begin
      ; sqrt-spaced such that uniform field would have equidistant lines
      levels = fllevels(pot,abs(nlines))
    endelse
    contour, pot, avg.rcyl, avg.z, $
        /OVERPLOT, $
        LEVELS=levels
  endif else begin
    ;; Plot arrows on color
    plot_3d_vect, var1, var2, var3_plot, $
        avg.rcyl, avg.z, $
        POS=pos, XRANGE=minmax(avg.rcyl), XSTYLE=1, $
        YRANGE=minmax(avg.z), YSTYLE=1, $
        /KEEP, $
        XTITLE='!8r!X', YTITLE='!8z!X', $
        MAXVEC=maxvec, CHARSIZE=charsize, $
        /QUIET
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
      endif else begin          ; at least two radii -> [inner, outer]
        rout = radii[1]
        weight2 = weight*((r_spher gt radii[0]) and (r_spher lt rout))
      endelse
      ; weight2 = weight2/rms(weight2)
      weight2 = weight2/mean(weight2)
    endif else begin
      rout = 1.e37
      r_spher = 0.*var1
      weight2 = 0.
    endelse
    ; vrm=rms(weight*var1) & vrm2=rms(weight2*var1)
    ; vzm=rms(weight*var2) & vzm2=rms(weight2*var2)
    ; vpm=rms(weight*var3) & vpm2=rms(weight2*var3)
    vrm=sqrt(mean(weight*var1^2)) & vrm2=sqrt(mean(weight2*var1^2))
    vzm=sqrt(mean(weight*var2^2)) & vzm2=sqrt(mean(weight2*var2^2))
    vpm=sqrt(mean(weight*var3^2)) & vpm2=sqrt(mean(weight2*var3^2))

    sph = where(r_spher<rout)
    vrmin=min(abs(var1[sph])) & vrmax=max(abs(var1))
    vzmin=min(abs(var2[sph])) & vzmax=max(abs(var2))
    vpmin=min(abs(var3[sph])) & vpmax=max(abs(var3))
    vv = sqrt(var1^2+var2^2+var3^2)
    vvmin=min(vv[sph])        & vvmax=max(vv)

    print, 'component:    min         overall rms   rms over sphere   max'
    print, '------------------------------------------------------------------'
    print, 'v_rcyl   :', vrmin, vrm, vrm2, vrmax
    print, 'v_phi    :', vpmin, vpm, vpm2, vpmax
    print, 'v_z      :', vzmin, vzm, vzm2, vzmax
    print, '------------------------------------------------------------------'
    print, 'vv       :', vvmin, sqrt(vrm^2+vpm^2+vzm^2), $
        sqrt(vrm2^2+vpm2^2+vzm2^2), vvmax

    r_cyl = spread(avg.rcyl,1,nz)
    rho = avg.rhomphi
    E_pol = mean(weight2*0.5*rho*(var1^2+var2^2))
    E_tor = mean(weight2*0.5*rho*var3^2)
    E_rot = mean(weight2*0.5*rho*r_cyl^2*Omega^2)
    print
    print, '<E_kin>_sphere = '
    print, '  pol ( /Erot):               ', E_pol, E_pol/E_rot
    print, '  tor ( /Erot):               ', E_tor, E_tor/E_rot
    print, 'E_rot ='
    print, '  <1/2*rho r_cyl^2.>_sphere = ', E_rot

    if (noomega) then print, 'WARNING: Assuming Omega=1.'

    print, mean(weight2*r_spher^2)
    print, mean(weight2*spread(avg.rcyl,1,nz)^2)
    print, 'total mass: ', mean(weight2*rho)*4./3.*!pi

  endif

end
; End of file pvv_phiavg.pro
