;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   plot_3d_vect.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   18-Aug-2000
;;;  Version: 1.37
;;;  CVS: $Revision: 1.8 $
;;;
;;;  Description:
;;;   Given either a 3d array on a 2d grid (e.g. B(0:64,0:64,15,0:3)),
;;;   or three 1d arrays on a 2d grid (e.g. Bx(0:64,0:64,15),
;;;   By(..), Bz(..)), plot the vertical component (the third cone) as
;;;   coloured background and overlay the horizontal components with
;;;   wdvelovect.
;;;
;;;  Usage:
;;;    PLOT_3D_VECT, VEC, [X, Y,]
;;;                  SAMPLE=sample,
;;;                  PERMUTE=permute, AUTOXYTITLE=autoxy, $
;;;                  PS_RGB=ps_rgb, PS_COL_MIN=ps_col_min, $
;;;                  TV_RGB=tv_rgb, $
;;;                  PS_ENH=ps_enh, $
;;;                  BITS=bits, $
;;;                  ZRANGE=zrange, $
;;;                  /KEEP_CT, $
;;;                  /TRUE_COLOR, $
;;;                  /PS_FULL_RANGE
;;;
;;;  or
;;;    PLOT_3D_VECT, VEC_x, VEC_y, VEC_z, [X, Y,] $
;;;                 .........
;;;
;;;  Arguments:
;;;   VEC       ---  An array of size (NX, NY, 3) containing the vector
;;;                  field to be plotted
;;;   VEC_X/Y/Z ---  Three arrays of size (NX,NY) containing the vector
;;;                  field to be plotted
;;;   X, Y      ---  Arrays of size (NX) and (NY) containing the
;;;                  (equidistant) coordinates of the grid axis
;;;
;;;  Keywords:
;;;   SAMPLE        --- A scalar or 1-d array containing a factor by
;;;                     which the vec_x and vec_y should be downsampled.
;;;   PERMUTE       --- A 3-d array, specifying the order in which the
;;;                     data array components are to be used, e.g.
;;;                       plot_3d_vect, uu[*,0,*], PERMUTE=[0,2,1], x, z
;;;                     will plot an x-z section
;;;   AUTOXYTITLE   --- Automatically generate XTITLE and YTITLE based
;;;                     on PERMUTE (and assuming the canonical order x-y-z)
;;;   NEGATE3       --- Multiply 3rd component (the colour-coded one) with -1
;;;   PS_RGB        --- A (n_cols,3) array defining the red, green, and
;;;                     blue values of the PostScript colortable
;;;   PS_COL_MIN    --- Index of the first colour (usually the darkest
;;;                     one) to be used in PostScript.This is to make the
;;;                     overlayed arrows visible.
;;;   TV_RGB        --- Either: a (n_cols,3) array defining the red, green,
;;;                     and blue values of the computer display colour
;;;                     table;
;;;                     or: a number, indexing one of several predefined
;;;                     colour tables
;;;   PS_ENH        --- Factor by which the number of grid points should
;;;                     be enhanced (by interpolating data) for the
;;;                     vertical component
;;;   BITS          --- Colour depth for PostSCript device (default = 8)
;;;   ZRANGE        --- interval for z component (mapped to colours
;;;                     [0, ncolours])
;;;   KEEP_CT       --- Flag for keeping the colour table instead of
;;;                     creating a new one.
;;;   TRUE_COLOR    --- Flag for explicitly choosing TrueColor (if set to
;;;                     non-zero value) or PseudoColor (if set to zero).
;;;                     If this keyword is not set, plot_3d_vect
;;;                     determines the visual on its own.
;;;   PS_FULL_RANGE --- Flag telling plot_3d_vect to use the full
;;;                     available colour range for PostScript. This
;;;                     will in general shift the colour for value zero
;;;   XRANGE,YRANGE --- Set x- and y-range
;;;   DEBUG         --- Write diagnostic output
;;;
;;;  Any other keywords are handed on to WDVELOVECT.
;;;

PRO plot_3d_vect, arg1, arg2, arg3, $
                  arg4, arg5, $
                  SAMPLE=sample, $
                  PERMUTE=perm, AUTOXYTITLE=autoxy, NEGATE3=negate3, $
                  PS_RGB=ps_rgb, PS_COL_MIN=ps_col_min, $
                  TV_RGB=tv_rgb, $
                  PS_ENH=ps_enh, $
                  BITS=bits, $
                  ZRANGE=zrange, $
                  KEEP_CT=keep_ct, $
                  TRUE_COLOR=true_col, $
                  PS_FULL_RANGE=ps_full_range, $
                  XRANGE=xrange, YRANGE=yrange, $
                  DEBUG=debug, HELP=help, _EXTRA=extra

  COMMON display1, visual, called, tv_rgb_list ; Remember these variables

  if (keyword_set(help)) then extract_help, 'plot_3d_vect'

;  ON_ERROR, 2                   ;Return to caller if an error occurs

  pmulti0 = !p.multi

  default, ps_col_min, 0
  default, bits, 8
  default, perm, [0,1,2]
  default, negate3, 0
  default, debug, 0

  xytit = '!8' + ['x','y','z'] + '!X'

  if (debug) then $
      print, FORMAT='(A, 10(I3))', $
      'PLOT_3D_VECT: Initial !p.multi = ', !p.multi

  if (n_elements(xrange) gt 1) then begin
    xbounds = xrange
  endif else begin
    xbounds = !x.range
  endelse
  if (n_elements(yrange) gt 1) then begin
    ybounds = yrange
  endif else begin
    ybounds = !y.range
  endelse
  ;; Make sure these ranges make sense
  if ((xbounds[1]-xbounds[0]) eq 0) then xbounds = [-1e30,1e30]
  if ((ybounds[1]-ybounds[0]) eq 0) then ybounds = [-1e30,1e30]

  IF (N_ELEMENTS(called) LT 1) THEN BEGIN ; Do this just at the first call
    if (!p.multi[0] eq 0) then begin
      erase                     ; Sometimes necessary to get n_cols right
    endif
  ENDIF
  n_cols = !D.TABLE_SIZE        ; Number of colours

  IF (NOT KEYWORD_SET(ps_rgb)) THEN BEGIN
    ps_rgb = intarr(n_cols,3)
    ps_rgb(*,0) = ps_col_min + $
        findgen(n_cols)*(n_cols-ps_col_min-1.)/(n_cols-1.)
    ps_rgb(*,1) = ps_rgb(*,0)
    ps_rgb(*,2) = ps_rgb(*,1)
  ENDIF

;; Set up a few colour tables
  IF (N_ELEMENTS(called) LT 1) THEN BEGIN ; Do this just at the first call
    tv_rgb_list = intarr(n_cols,3,2)
    phi = (findgen(n_cols-2)+1)/(n_cols-1)*!PI/2
    ;; First colour table
    tv_rgb_list(*,0,0) = byte([0, 255*sin(phi), 255])
    tv_rgb_list(*,1,0) = byte([0, 255/sqrt(2.)*sin(2*phi)^2, 255])
    tv_rgb_list(*,2,0) = byte([0, 255*cos(phi), 255])
    ;; Second colour table
    tv_rgb_list(*,0,1) = [0, bindgen(n_cols-2)+1, 255]
    tv_rgb_list(*,1,1) = [0, bindgen(n_cols/2-1)+1, $
                  n_cols-bindgen((n_cols-1)/2)-1-n_cols/2, 255]
    tv_rgb_list(*,2,1) = [0, n_cols-bindgen(n_cols-2)-1, 255]
  ENDIF
  IF (N_ELEMENTS(tv_rgb) EQ 0) THEN BEGIN 
    tv_rgb = reform(tv_rgb_list(*,*,0))
  ENDIF ELSE BEGIN
    IF ((SIZE(tv_rgb))(0) EQ 0) THEN BEGIN ; If tv_rgb is scalar then..
      IF (tv_rgb GE (size(tv_rgb_list))(3)) THEN BEGIN
        MESSAGE, "Not enough colour tables in list."
      ENDIF
      tv_rgb = reform(tv_rgb_list(*,*,tv_rgb)) ; ..use it as index for the list
    ENDIF
  ENDELSE

  IF (KEYWORD_SET(ps_full_range)) THEN absolut = 0 ELSE absolut = 1

  ;; Get the name of the visual from IDL
  if (!d.name eq 'X') then begin
    if ((n_elements(called) lt 1) or (n_elements(visual) le 0)) then begin
      ;; Do this just at the first call or if the first call was for a
      ;; non-X device
      device, GET_VISUAL_NAME=visual
    endif
    if (visual eq 'TrueColor') then true_color = 1
  endif

  IF (N_ELEMENTS(true_col) NE 0) THEN BEGIN ; Overrides any previous choices
    IF (true_col NE 0) THEN true_color = 1 ELSE true_color = 0
  ENDIF

  IF (N_ELEMENTS(true_color) EQ 0) THEN true_color = 0

  IF (NOT KEYWORD_SET(keep_ct)) THEN keep_ct = 0 
  IF (keep_ct) THEN true_color = 0 ; Helps a lot

  ;; Process the data
  vec_x = reform(arg1)
  IF ((size(vec_x))[0] EQ 3) THEN BEGIN
    ;; If the first argument has dimensions (NX, NY, 3)
    IF (N_ELEMENTS(arg3) NE 0) THEN BEGIN
      ;; three arguments: data_matrix, x, y
      x = arg2
      y = arg3
    ENDIF
    vec = vec_x
    vec_x = vec(*,*,0)
    vec_y = vec(*,*,1)
    vec_z = vec(*,*,2)
  ENDIF ELSE BEGIN              ; First ARG1 has dimensions (NX, NY)
    vec_x = arg1
    vec_y = arg2
    vec_z = arg3
    IF (N_ELEMENTS(arg5) NE 0) THEN BEGIN
      ;; five arguments: data_x, data_y, data_z, x, y
      x = arg4
      y = arg5
    ENDIF
  ENDELSE
  vec_all = [[[reform(vec_x)]], [[reform(vec_y)]], [[reform(vec_z)]]] 
;  vec_x = reform(vec_x)
;  vec_y = reform(vec_y)
;  vec_z = reform(vec_z)
  ;; Apply permutation
  ;; 1. Determine parity of permutation.
  ;;    Even permutation (parity=+1) --> right-handed system
  ;;    Odd  permutation (parity=-1) --> left-handed system
  parity = sign(perm[1]-perm[0])+sign(perm[2]-perm[1])+sign(perm[0]-perm[2])
  if (negate3) then parity = -parity
  vec_x = vec_all[*,*,perm[0]]
  vec_y = vec_all[*,*,perm[1]]
  vec_z = vec_all[*,*,perm[2]]*parity
  v_size = size(vec_x)
  nx = v_size(1)
  ny = v_size(2)
  nx_s = nx
  ny_s = ny

  IF (N_ELEMENTS(x) EQ 0) THEN  x = findgen(nx)
  IF (N_ELEMENTS(y) EQ 0) THEN  y = findgen(ny)
  x_s = congrid(x, nx_s)
  y_s = congrid(y, ny_s)
  ;; Honour xrange, yrange:
  goodx = where((x_s ge xbounds[0]) and (x_s le xbounds[1]))
  if (goodx[0] lt 0) then message, 'No valid points within xrange'
  x_s = x_s[goodx]
  nx_s = n_elements(x_s)
  vec_x = vec_x[goodx,*]
  vec_y = vec_y[goodx,*]
  vec_z = vec_z[goodx,*]
  ;
  goody = where((y_s ge ybounds[0]) and (y_s le ybounds[1]))
  if (goody[0] lt 0) then message, 'No valid points within yrange'
  y_s = y_s[goody]
  ny_s = n_elements(y_s)
  vec_x = vec_x[*,goody]
  vec_y = vec_y[*,goody]
  vec_z = vec_z[*,goody]

  dx_s = x_s(1)-x_s(0)
  dy_s = y_s(1)-y_s(0)

  IF (NOT KEYWORD_SET(ps_enh)) THEN ps_enh = 1

  IF (!D.NAME EQ 'PS') THEN  DEVICE, /COLOR, BITS=bits

  called = 1

;;;;;;;

;; Establish coordinates and size of plot
  xr = minmax(x_s)+dx_s*[-0.5,0.5]
  yr = minmax(y_s)+dy_s*[-0.5,0.5]
  idxr = (xr-x[0])/(x[1]-x[0])  ; index range to be plotted
  idyr = (yr-y[0])/(y[1]-y[0])  ; index range to be plotted
  ;;
  default, xrange, !x.range
  default, yrange, !y.range
  ;; Do the sampling if needed
  if (n_elements(sample) ne 0) then begin
    samp_x = sample[0] > 1
    if (n_elements(sample) eq 1) then begin
      samp_y = sample[0] > 1
    endif else begin
      samp_y = sample[1] > 1
    endelse
  endif else begin
    samp_x = 1
    samp_y = 1
  endelse
  ixvel = indgen(nx_s/samp_x)*samp_x
  iyvel = indgen(ny_s/samp_y)*samp_y
  ;; Awkward, but IDL does not support start:stop:stride notation, nor
  ;; does it allow for the array subscripts a la vec[ixvel,iyvel]
  vec_x = vec_x[ixvel,*]
  vec_x = vec_x[*,iyvel]
  vec_y = vec_y[ixvel,*]
  vec_y = vec_y[*,iyvel]
  wdvelovect, vec_x, vec_y, x_s[ixvel], y_s[iyvel],  $
      XRANGE=xrange, YRANGE=yrange, /NOPLOT, _EXTRA=extra
  sx = (xr[1]-xr[0])*!D.X_VSIZE*!X.S[1] ; Size ...
  sy = (yr[1]-yr[0])*!D.Y_VSIZE*!Y.S[1] ; Size ...
  xpos = (!X.S[1]*xr[0] + !X.S[0]) * !D.X_VSIZE + 1 ; ..and position of image
  ypos = (!Y.S[1]*yr[0] + !Y.S[0]) * !D.Y_VSIZE + 1 ; ..and position of image
  IF (!D.NAME EQ 'PS') THEN BEGIN
;    TVLCT, ps_rgb(*,0), ps_rgb(*,1), ps_rgb(*,2)
    if (n_elements(zrange) gt 0) then message, $
        "Warning: keyword ZRANGE not yet implemented for PS output", $
        /INFORMATIONAL
    wdtvscl, $
        interpolate(vec_z, $
                    (findgen(nx_s*ps_enh)+0.5)/ps_enh-0.5, $
                    (findgen(ny_s*ps_enh)+0.5)/ps_enh-0.5, $
                    /grid), $
        xpos, ypos, $
        XSIZE=sx, YSIZE=sy, $
        BOTTOM=PS_COL_MIN, $
        ABSOLUT=absolut
  ENDIF ELSE BEGIN
    IF ((NOT true_color) AND (NOT keep_ct)) THEN  $
        TVLCT,  tv_rgb(*,0), tv_rgb(*,1), tv_rgb(*,2)
    if (n_elements(zrange) eq 0) then begin
      ;; need a copy, so PLOT_#D_VECT, ..,ZRANGE=undef
      ;; will not set undef to a value if it was previously undefined
      zrange1=[-1,1]*MAX([-MIN(vec_z), MAX(vec_z)])
    endif else begin
      zrange1 = zrange
    endelse
    x_enh = sx/(idxr[1]-idxr[0])
    y_enh = sy/(idyr[1]-idyr[0])
    nx1 = floor(nx_s*x_enh)
    ny1 = floor(ny_s*y_enh)

    image = 1+(interpolate(vec_z, $
                           (findgen(nx1)+0.5)/x_enh-0.5, $
                           (findgen(ny1)+0.5)/y_enh-0.5, $
                           /grid) $
               -zrange1[0])/(zrange1[1]-zrange1[0])*(n_cols-3.)
    IF (true_color) THEN BEGIN
      true_image = lonarr(nx1,ny1,3)
      FOR j = 0,2 DO BEGIN
        true_image[*,*,j] = reform(tv_rgb[image,j], nx1, ny1)
      ENDFOR
      TV, true_image, xpos, ypos, true=3
    ENDIF ELSE BEGIN
      TV, image, xpos, ypos
    ENDELSE
  ENDELSE

  pmulti1 = !p.multi
  !p.multi = pmulti0
  if (keyword_set(autoxy)) then begin
    wdvelovect, vec_x , vec_y, x_s[ixvel], y_s[iyvel], /NOERASE, $
        XRANGE=xrange, YRANGE=yrange, $
        XTITLE=xytit[perm[0]], YTITLE=xytit[perm[1]], $
        _EXTRA=extra
  endif else begin
    wdvelovect, vec_x , vec_y, x_s[ixvel], y_s[iyvel], /NOERASE, $
        XRANGE=xrange, YRANGE=yrange, _EXTRA=extra
  endelse
  !p.multi = pmulti1

  if (debug) then begin
    print, FORMAT='(A, 10(I3))', $
        'PLOT_3D_VECT: Final !p.multi   = ', !p.multi
    print, '------------------------------------------------'
  endif


END

; End of file plot_3d_vect.pro
