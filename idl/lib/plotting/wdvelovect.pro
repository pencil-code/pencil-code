; $Id$
;
; Copyright (c) 1983-1998, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.
; Changes to original routine velovect by Martin Schultz and Franz Rohrer.
;
; Changes to this improved routine by wd (Wolfgang.Dobler@ncl.ac.uk)

PRO WDVELOVECT,U_in,V_in,X_in,Y_in, $
        Missing = Missing, Length = length, Dots = dots,  $
        Color=color, CLIP=clip, NOCLIP=noclip, OVERPLOT=overplot,  $
        NOZERO=nozero, CENTER=center, XRANGE=xrange, YRANGE=yrange,  $
        POSITION=position, POLAR=POLAR, LINEONLY=LINEONLY, $
        ScaLength = scalength, MAXVEC=maxvec, NOPLOT=noplot, QUIET=quiet, $
        _EXTRA=extra
;
;+ 
; NAME:
;	WDVELOVECT
;
; PURPOSE:
;	Produce a two-dimensional velocity field plot.
;
;	A directed arrow is drawn at each point showing the direction and 
;	magnitude of the field.
;
;       This routine was developed from the original VELOVECT routine
;       and now allows two dimensional arrays for X and Y for irregularily
;       spaced data.
;               
; CATEGORY:
;	2D Graphics
;
; CALLING SEQUENCE:
;	WDVELOVECT, U, V [, X, Y]
;
; INPUTS:
;	U:	The X component of the two-dimensional field.  
;		U must be a two-dimensional array.
;
;	V:	The Y component of the two dimensional field.  Y must have
;		the same dimensions as X.  The vector at point [i,j] has a 
;		magnitude of:
;
;			(U[i,j]^2 + V[i,j]^2)^0.5
;
;		and a direction of:
;
;			ATAN2(V[i,j],U[i,j]).
;
;  OR, if the keyword POLAR is set then:
;
;	U:	The magnitude of the two-dimensional field.  
;		U must be a two-dimensional array.
;
;	V:	The angular direction of the two dimensional field (in radians).
;               Y must have same dimensions as X.  
;               Angles are measured as per the computing standard i.e. 
;               anti-clockwise from the +ve x-axis.
;
; OPTIONAL INPUT PARAMETERS:
; 	X:	Optional abcissae values.  X must be a vector with a length 
;		equal to the first dimension of U and V *OR* a 2-dimensional
;               array with the same dimensions as U and V.
;
;	Y:	Optional ordinate values.  Y must be a vector with a length
;		equal to the first dimension of U and V *OR* a 2-dimensional
;               array with the same dimensions as U and V.
;
; KEYWORD INPUT PARAMETERS:
;	COLOR:	The color index used for the plot.
;
;	DOTS:	Set this keyword to 1 to place a dot at each missing point. 
;		Set this keyword to 0 or omit it to draw nothing for missing
;		points.  Has effect only if MISSING is specified.
;
;        POLAR: Use the U and V specified as the magnitude and direction,
;               respectively, of the two-dimensional vector field. 
;               (See U and V above)
;
;
;	LENGTH:	Length factor.  The default of 1.0 makes the longest (U,V)
;		vector the length of a cell.
;
;       SCALENGTH:  Franz Rohrer's modification of the LENGTH keyword:
;              SCALENGTH applies a scale factor relative to the data
;              values, not just a relative scaling.
;              [wd: renamed this from LENGTH. The new behaviour is
;              counterintuitive and certainly should never have
;              changed the behaviour of LENGTH.]
;
;       MISSING: Missing data value.  Vectors with a LENGTH greater
;		than MISSING are ignored.
;
;       OVERPLOT: Set this keyword to make VELOVECT "overplot".  That is, the
;               current graphics screen is not erased, no axes are drawn, and
;               the previously established scaling remains in effect.
;
;       NOZERO: Do not plot zero vectors as dots.
;
;       CENTER: Position of the grid point on the arrow (CENTER=0
;               corresponds to the default behavior of velovect,
;               CENTER=0.5 centers the arrow w.r.t. the grid point and
;               is the default here).
;
;       MAXVEC: Maximum number of arrows; either a scalar (limit total
;               number) or an array [maxvecx,maxvecy]. A value of zero
;               means no limit. The default tries to keep the arrows
;               visible, based on screen resolution.
;
;       NOPLOT: Do not plot any arrows. Useful for just establishing
;               the coordinate system like PLOT_3D_VECT does.
;
;       QUIET:  Don't print informational messages
;
;	Note:   All other keywords are passed directly to the PLOT procedure
;		and may be used to set option such as TITLE, POSITION, 
;		NOERASE, etc.
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Plotting on the selected device is performed.  System
;	variables concerning plotting are changed.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Straightforward.  Unrecognized keywords are passed to the PLOT
;	procedure.  
;
; MODIFICATION HISTORY:
;	DMS, RSI, Oct., 1983.
;	For Sun, DMS, RSI, April, 1989.
;	Added TITLE, Oct, 1990.
;	Added POSITION, NOERASE, COLOR, Feb 91, RES.
;	August, 1993.  Vince Patrick, Adv. Visualization Lab, U. of Maryland,
;		fixed errors in math.
;	August, 1993. DMS, Added _EXTRA keyword inheritance.
;	January, 1994, KDB. Fixed integer math which produced 0 and caused
;		            divide by zero errors.
;	December, 1994, MWR. Added _EXTRA inheritance for PLOTS and OPLOT.
;	June, 1995, MWR. Removed _EXTRA inheritance for PLOTS and changed
;			 OPLOT to PLOTS.
;       September, 1996, GGS. Changed denominator of x_step and y_step vars. 
;       February, 1998, DLD.  Add support for CLIP and NO_CLIP keywords.
;       June, 1998, DLD.  Add support for OVERPLOT keyword.
;       16 Sep 1999: Martin Schultz added support for 2D U and V arrays
;              cleaned up the routine some and added the NOZERO keyword.
;              (renamed as msvelovect.pro on Sep 23)
;              Also included Franz Rohrer's modification of the LENGTH keyword:
;              LENGTH now applies a scale factor relative to the data values,
;              not just a relative scaling. If you don't specify LENGTH it acts
;              as before.
;       29 Feb 2000: mgs. Bug fix: for loop needs long int!
;              Also sx and sy are now always computed even if no x or
;              y are passed into routine.
;       26 Jun 2001: wd. Center arrows and add CENTER keyword
;       27 Jun 2001: wd. Made angles of arrow heads device independent
;        4 Jul 2001: wd. Intercept keywords [XY]RANGE
;       19 Jul 2001: wd. Apply reform() to arguments U and V
;       22 Mar 2002: wd. Reverted F. Rohrer's silly LENGTH behaviour to original
;       22 Mar 2002: wd. Added MAXVEC, QUIET and NOPLOT keywords
;       07 Aug 2004: Antony Mee, U. of Newcastle upon Tyne
;              Added POLAR keyword to allow specification of vectors as
;              magnitude and direction.  Also added LINEONLY to prevent
;              drawing an arrow head.
;-
;
        on_error,2                      ;Return to caller if an error occurs
;
        if (n_elements(quiet) le 0) then quiet = 0
        if (n_elements(noplot) le 0) then noplot = 0
;
        u = reform(u_in)
        v = reform(v_in)
        s = size(u)
        t = size(v)
        if s[0] ne 2 then begin 
baduv:     message, 'U and V parameters must be 2D and same size.'
        endif
        if total(abs(s[0:2]-t[0:2])) ne 0 then goto,baduv
;
        if n_params() lt 3 then begin
           x = findgen(s[1]) 
           sx = size(x)
        endif else begin
           x = x_in
           sx = size(x)
           if (sx[0] eq 2) then begin
              if total(abs(sx[0:2]-s[0:2])) ne 0 then begin
badx:            message, 'X array has incorrect size.'
              endif
           endif else $
              if n_elements(x) ne s[1] then goto,badx
        endelse
;
        if n_params() lt 4 then begin
           y = findgen(s[2]) 
           sy = size(y)
        endif else begin
           y = y_in
           sy = size(y)
           if (sy[0] eq 2) then begin
              if (sx[0] ne 2) then goto,bady
              if total(abs(sy[0:2]-s[0:2])) ne 0 then begin
bady:            message, 'Y array has incorrect size.'
              endif
           endif else $
              if n_elements(y) ne s[2] then goto,bady
        endelse
;
; wd: regrid data if too many points
;
        if (not keyword_set(overplot)) then begin
          ;---------  pretend to plot to get plot window size right  ---------
          if (n_elements(position) lt 4) then begin
            plot, x, y, /NODATA, /NOERASE,XSTYLE=4, YSTYLE=4, TITLE=''
          endif else begin
            plot, x, y, /NODATA, /NOERASE,XSTYLE=4, YSTYLE=4, TITLE='', $
                POSITION=position
          endelse
        endif
        ; only now we can extract position
        position=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
        xwidth = !x.window[1] - !x.window[0]
        ywidth = !y.window[1] - !y.window[0]
        case n_elements(maxvec) of
          ;
          0: begin
            ;; default settings: maxvec depends on resolution
            if ((!d.flags and 1) eq 0) then begin ; fixed-size pixels
              npix = 9     ; make longest arrows that many pixels long
              maxvec = [!d.x_vsize*xwidth, !d.y_vsize*ywidth]/npix
            endif else begin    ; scalable pixels
              nmms = 2.0      ;  make longest arrows that many mm long
              maxvec = [!d.x_vsize*xwidth, !d.y_vsize*ywidth]/(nmms*100)
            endelse
          end
          ;
          1: begin
            ;; One arg given : replicate for both directions according
            ;; to number of points in x and y
            maxvec = floor([sqrt(maxvec*sx/sy), sqrt(maxvec*sy/sx)])
          end
          ;
          else:                 ; nothing to do
          ;
        endcase

        regrid = 0
        nx = s[1]
        ny = s[2]
        if ( (maxvec[0] gt 0) and (s[1] gt maxvec[0])) then begin
          nx = maxvec[0]
          regrid = 1
        endif
        if ( (maxvec[1] gt 0) and (s[2] gt maxvec[1])) then begin
          ny = maxvec[1]
          regrid = 1
        endif
        nx = floor(nx)
        ny = floor(ny)
        if (regrid) then begin
          if ((quiet eq 0) and (noplot eq 0)) then begin
            message, /INFO, $
                'Regridding from ' + strtrim(s[1],2) + 'x' + strtrim(s[2],2) $
                + ' to ' + strtrim(nx,2) + 'x' + strtrim(ny,2)
          endif
          x = linspace(minmax(x),nx,GHOST=0.5)
          y = linspace(minmax(y),ny,GHOST=0.5)
          u = congrid(u, nx, ny, /CUBIC)
          v = congrid(v, nx, ny, /CUBIC)
          ; Recalculate these:
          s = size(u)
          t = size(v)
          sx = size(x)
          sy = size(y)
        endif
;
        if (n_elements(center) le 0) then center = 0.5
;
        if n_elements(missing) le 0 then missing = 1.0e30
; ### FR: use LENGTH differently -- allows absolute scaling
;       if n_elements(length) le 0 then length = 1.0

        if keyword_set(POLAR) then begin
          mag = u
        endif else begin
          mag = sqrt(u^2.+v^2.)             ;magnitude.
        endelse
                ;Subscripts of good elements
        nbad = 0                        ;# of missing points
; ## mgs: because of defaulting 5 lines above, missing always has a value!!
;       if n_elements(missing) gt 0 then begin
           good = where(mag lt missing) 
           if keyword_set(dots) then bad = where(mag ge missing, nbad)
;       endif else begin
;               good = lindgen(n_elements(mag))
;       endelse
        ugood = u[good]
        vgood = v[good]
        x0 = min(x)                     ;get scaling
        x1 = max(x)
        y0 = min(y)
        y1 = max(y)
	x_step=(x1-x0)/(s[1]-1.0)   ; Convert to float. Integer math
	y_step=(y1-y0)/(s[2]-1.0)   ; could result in divide by 0

        if keyword_set(POLAR) then begin
	  maxmag=max(abs(u))/min([x_step,y_step])
;AJWM:   Strange but true!
;        sina and cosa are switched to use angles measured CCW from the +ve 
;        +ve x-axis switching them back gives angles measured CW from 
;        the +ve y-axis; the (x,y)-component based code considers the
;        angles of vectors such.
;
;        Could consider some keywords to allow switching between the two
;        methods.  But will stick with the IDL/computing standard for now.
;
          default,length,1.0
          cosa = length*ugood/maxmag*sin(vgood)
          sina = length*ugood/maxmag*cos(vgood)
        endif else begin
	  maxmag=max([max(abs(ugood/x_step)),max(abs(vgood/y_step))])
; ### FR:
          if n_elements(scalength) gt 0 then maxmag=scalength/x_step
          sina = (ugood/maxmag)
          cosa = (vgood/maxmag)
; ### wd: LENGTH overrides SCALENGTH
          if (n_elements(length) gt 0) then begin
            sina = length * (ugood/maxmag)
            cosa = length * (vgood/maxmag)
          endif
        endelse
;
        if n_elements(title) le 0 then title = ''
        ;--------------  plot to get axes  ---------------
        if n_elements(color) eq 0 then color = !p.color
        if n_elements(noclip) eq 0 then noclip = 0
        x_b0=x0-x_step
	x_b1=x1+x_step
	y_b0=y0-y_step
	y_b1=y1+y_step
        if (not keyword_set(overplot)) then begin
;         if n_elements(position) eq 0 then begin
          ;; The following would be clearer in the form
          ;;   if ((n_elts(xr) le 1) or (xr[0] eq xr[1])) then xr = [..] ,
          ;; but IDL doesn't support logical shortcircuiting
          if (n_elements(xrange) le 1) then xrange=[x_b0,x_b1]
          if (xrange[0] eq xrange[1])  then xrange=[x_b0,x_b1]
          if (n_elements(yrange) le 1) then yrange=[y_b0,y_b1]
          if (yrange[0] eq yrange[1])  then yrange=[y_b0,y_b1]
          plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
              color=color, xrange=xrange, yrange=yrange, POSITION=position, $
              _EXTRA = extra
;         endif else begin
;           plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
;             color=color, _EXTRA = extra
;         endelse
        endif
; Only plot if NOPLOT is not set
        if (noplot eq 0) then begin
            if n_elements(clip) eq 0 then $
                clip = [!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]]
;
            r = .3                ;len of arrow head
            angle = 22.5 * !dtor  ;Angle of arrowhead
            st = r * sin(angle)   ;sin 22.5 degs * length of head
            ct = r * cos(angle)
;
; Make angles of arrow head device independent:
            xd=(!x.crange[1]-!x.crange[0])/(!x.window[1]-!x.window[0])/!d.x_size
            yd=(!y.crange[1]-!y.crange[0])/(!y.window[1]-!y.window[0])/!d.y_size
;
            for i=0L,n_elements(good)-1 do begin ;Each point
                if (sx[0] eq 2) then begin
                    x0 = x[good[i]] ;get coords of start & end
                    y0 = y[good[i]]
                endif else begin
                    x0 = x[good[i] mod s[1]] ;get coords of start & end
                    y0 = y[good[i] / s[1]]
                endelse
                dx = sina[i]
                x0 = x0 - center*dx
                x1 = x0 + dx
                dy = cosa[i]
                y0 = y0 - center*dy
                y1 = y0 + dy
;	         xd=x_step
;	         yd=y_step
                ; plot zero vectors as dots
                if (mag[good[i]] eq 0.) then begin
                    if  (not keyword_set(NOZERO)) then $
                        plots,x0,y0,psym=3,color=color,clip=clip,  $
                        noclip=noclip  
                endif else $
                  if keyword_set(LINEONLY) then begin
                    plots,[x0,x1],[y0,y1], $
                           color=color,clip=clip,noclip=noclip
                  endif else begin
                    plots,[x0,x1,x1-(ct*dx/xd-st*dy/yd)*xd, $
                           x1,x1-(ct*dx/xd+st*dy/yd)*xd], $
                          [y0,y1,y1-(ct*dy/yd+st*dx/xd)*yd, $
                           y1,y1-(ct*dy/yd-st*dx/xd)*yd], $
                           color=color,clip=clip,noclip=noclip
                  endelse
            endfor
            if nbad gt 0 then begin
                if (sx[0] eq 2) then begin ;Dots for missing?
                    PLOTS, x[bad], y[bad], psym=3, $
                        color=color, clip=clip,noclip=noclip   
                endif else begin
                    PLOTS, x[bad mod s[1]], y[bad / s[1]], psym=3, $
                        color=color, clip=clip,noclip=noclip
                endelse
            endif
        endif
;
end
