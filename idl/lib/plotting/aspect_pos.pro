;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   aspect_pos.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  $Date: 2004-05-14 11:59:24 $
;;;  $Revision: 1.7 $
;;;
;;;  Description:
;;;   Return position [X0, Y0, X1, Y1] (in normalized coordinates) for
;;;   placing a graph with given aspect ratio RATIO. Works with
;;;   !p.multi and for all devices.
;;;     Based on aspect.pro by David Fanning.
;;;   Usage:
;;;     plot,a,b,POS=aspect_pos(Lb/La)
;;;     plot,a,b,POS=aspect_pos(Lb/La,MARGIN=0.1)
;;;     plot,a,b,POS=aspect_pos(Lb/La,npos,nx,ny)
;;;   Arguments:
;;;     NPOS   -- Position number in the NX-NY grid as indicated by !p.multi[0]
;;;     NX     -- Number of columns in layout grid as indicated by !p.multi[1]
;;;     NY     -- Number of columns in layout grid as indicated by !p.multi[2]
;;;     All three default to the values obtained from !p.multi
;;;   Keywords:
;;;     MARGIN -- Leave (at least) this fraction as margin on either side.
;;;               Can be (a) one number, specifying all four margins
;;;                      (b) a vector [margx, margy], specifying
;;;                          different margins for x and y
;;;                      (c) a vector [mxleft, mxright, mybot, mytop],
;;;                          specifying all four margins explicitly
;;;               Defaults to 0.05

function aspect_pos, ratio, npos, nx, ny, MARGIN=margin, $
                     DEBUG=debug
  default, debug, 0

  ;; Ensure safe value of ratio
  if (n_params() eq 0) then ratio = 1.0
  if (ratio eq 0) then begin
    message, 'ratio=0 is meaningless. Defaulting to 1.', /INFO
    ratio = 1.0
  endif

  ;
  ;  [npos,nx,ny] = !p.multi
  ;
  default, npos, !p.multi[0]
  default, nx, !p.multi[1]>1
  default, ny, !p.multi[2]>1
  ntot = nx*ny
  ;; Transform npos into x and y grid positions. Not as trivial as
  ;; expected, since IDL uses a layout like
  ;;    0  5
  ;;    4  3
  ;;    2  1
  ;; for !p.multi=[*,2,3]
  nposx = nx-1 - ((npos+nx-1) mod nx)
  nposy = ((npos+ntot-1) mod ntot) / nx

  default, margin, 0.05 
  ;; Turn margin into 4-vector [x_left, x_right, y_bot, y_top]
  if (n_elements(margin) eq 1) then begin
    margin = [1,1,1,1]*margin[0]
  endif else if (n_elements(margin) eq 2) then begin
    margin = [[1,1]*margin[0], [1,1]*margin[1]]
  endif
  if (n_elements(margin) ne 4) then $
    message, 'margin must have the form MXY or [MX, MY] or [MX0, MX1, MY0, MY1]'
  margin = ((margin > 0) < 0.5) ; sanitize values

  s = size(ratio)
  type = s[s[0]+1]
  if (type ne 4 and type ne 5) then $
      message, 'Warning: RATIO is not a FLOAT -- may have been rounded.', /INFO

  ;; Aspect ratio of current window.  
  wratio = float(!d.y_vsize)*(1-margin[2]-margin[3]) $
           / (!d.x_vsize*(1-margin[0]-margin[1]))
  wratio = wratio

  ;; Normalized positions in window.
  xcent = (nposx+0.5)/nx
  ycent = (nposy+0.5)/ny
  if (ratio*ny/nx le wratio) then begin
    if (debug) then message, 'Bounded by xwidth', /INFO
    xwidth = 1./nx
    ywidth = xwidth*(ratio/wratio)
  endif else begin
    if (debug) then message, 'Bounded by ywidth', /INFO
    ywidth = 1./ny
    xwidth = ywidth*(wratio/ratio)
  endelse
  xstart = xcent - (0.5-margin[0])*xwidth
  xend   = xcent + (0.5-margin[1])*xwidth
  ystart = ycent - (0.5-margin[2])*ywidth
  yend   = ycent + (0.5-margin[3])*ywidth
  position = [xstart, ystart, xend, yend]
  if (debug) then message, string('POSITION = ',position), /INFO
  return, position

end
; ---------------------------------------------------------------------- ;

; End of file aspect_pos.pro
