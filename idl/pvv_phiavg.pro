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
;;;    pvv_phiavg, [filename], [/BB], [/UU]
;;;    pvv_phiavg, [avg_struct], [/BB], [/UU]
;;;  IF FILENAME doesn't exist, try 'data/FILENAME' and
;;;  'data/averages/FILENAME'.
;;;  Examples:
;;;    pvv_phiavg                   ; plot uu data from PHIAVG1,
;;;    pvv_phiavg, /BB              ; plot bb data from PHIAVG1,
;;;    pvv_phiavg, 'PHIAVG3'        ; plot uu data from PHIAVG3,
;;;    pvv_phiavg, 'PHIAVG3', /BB   ; plot bb data from PHIAVG3,
;;;    pvv_phiavg, avg              ; use data from struct AVG
;;;    pvv_phiavg, avg, /STAT       ; print some rms values
;;;  Advanced example:
;;;    avg=tavg_phiavg([800.,1e10])
;;;    pvv_phiavg, avg, /BB, RADII=par2.r_ext, /STAT

pro pvv_phiavg, arg, BB=bb, UU=uu, $
                MARGIN=margin, RADII=radii, $
                CHARSIZE=charsize, MAXVEC=maxvec, $
                STAT=stat, QUIET=quiet

  datatopdir ='data'
  avgdir = datatopdir+'/averages'
  phiavgfile = 'PHIAVG1'

  default, arg, avgdir+'/'+phiavgfile
  default, uu, 1
  default, bb, 0
  default, margin, 0.05
  default, charsize, !p.charsize
  default, maxvec, [20,40]
  default, STAT, stat
  default, quiet, 0


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

  if (bb) then begin
    var1 = avg.brmphi
    var2 = avg.bzmphi
    var3 = avg.bpmphi
  endif else if (uu) then begin
    var1 = avg.urmphi
    var2 = avg.uzmphi
    var3 = avg.upmphi
  endif else begin
    message, 'Need one of UU or BB to be true'
  endelse


  pos = aspect_pos((max(avg.z)-min(avg.z))/(max(avg.rcyl)-min(avg.rcyl)), $
                   MARGIN=margin)
  plot_3d_vect, var1, var2, var3, $
      avg.rcyl, avg.z, $
      POS=pos, XRANGE=minmax(avg.rcyl), XSTYLE=1, $
      YRANGE=minmax(avg.z), YSTYLE=1, $
      /KEEP, $
      XTITLE='!8r!X', YTITLE='!8z!X', $
      MAXVEC=maxvec, CHARSIZE=charsize, $
      /QUIET

  ;; Overplot spheres (well, circles) if that makes sense
  if (n_elements(radii) gt 0) then opcircle, radii

  ;; Print some statistics. Should ideally be in another place (or
  ;; rather standalone), but it comes handy to do this here.
  if (stat) then begin
    nr = n_elements(avg.rcyl)
    nz = n_elements(avg.z)
    ;; weights for use with rms(weight*var):
    weight = spread(sqrt(avg.rcyl),1,nz)
    weight = weight/rms(weight)
    ; ;; weights for use with sqrt(mean(weight*var^2)):
    ;  weight = spread(avg.rcyl,1,nz)
    ;  weight = weight/mean(weight)

    ;; weights for averaging over sphere
    ;;   r_sphere < radii[0]
    ;; or
    ;;   radii[0] < r_sphere < radii[1]:
    ;;
    if (n_elements(radii) gt 0) then begin
      r_spher = sqrt(spread(avg.rcyl,1,nz)^2 + spread(avg.z,0,nr))
      if (n_elements(radii) eq 1) then begin ; just one radius -> outer
        weight2 = weight*(r_spher lt radii[0])
      endif else begin          ; at least two radii -> [inner, outer]
        weight2 = weight*((r_spher gt radii[0]) and (r_spher lt radii[1]))
      endelse
      weight2 = weight2/rms(weight2)
      ; weight2 = weight2/mean(weight2)
    endif else begin
      weight2 = 0.
    endelse
    vrm=rms(weight*var1) & vrm2=rms(weight2*var1)
    vpm=rms(weight*var2) & vpm2=rms(weight2*var2)
    vzm=rms(weight*var3) & vzm2=rms(weight2*var3)
    ; vrm=sqrt(mean(weight*var1^2)) & vrm2=sqrt(mean(weight2*var1^2))
    ; vpm=sqrt(mean(weight*var2^2)) & vpm2=sqrt(mean(weight2*var2^2))
    ; vzm=sqrt(mean(weight*var3^2)) & vzm2=sqrt(mean(weight2*var3^2))

    print, 'component:  overall rms    rms over sphere'
    print, '--------------------------------------------'
    print, 'v_rcyl   :', vrm, vrm2
    print, 'v_phi    :', vpm, vpm2
    print, 'v_z      :', vzm, vzm2
    print, '--------------------------------------------'
    print, 'vv       :', sqrt(vrm^2+vpm^2+vzm^2), sqrt(vrm2^2+vpm2^2+vzm2^2)

  endif


end
; End of file vv_phiavg.pro
