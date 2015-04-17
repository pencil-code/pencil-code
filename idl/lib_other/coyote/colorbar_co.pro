;+
; NAME:
;   COLORBAR_CO
;
; PURPOSE:
;       The purpose of this routine is to add a color bar to the current
;       graphics window.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;       Graphics, Widgets.
;
; CALLING SEQUENCE:
;       COLORBAR_CO
;
; INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;
;       BOTTOM:   The lowest color index of the colors to be loaded in
;                 the bar.
;
;       CHARSIZE: The character size of the color bar annotations. Default is 1.0.
;
;       COLOR:    The color index of the bar outline and characters. Default
;                 is !P.Color..
;
;       DIVISIONS: The number of divisions to divide the bar into. There will
;                 be (divisions + 1) annotations. The default is 6.
;
;       FONT:     Sets the font of the annotation. Hershey: -1, Hardware:0, True-Type: 1.
;
;       FORMAT:   The format of the bar annotations. Default is '(I5)'.
;
;       INVERTCOLORS: Setting this keyword inverts the colors in the color bar.
;
;       MAXRANGE: The maximum data value for the bar annotation. Default is
;                 NCOLORS.
;
;       MINRANGE: The minimum data value for the bar annotation. Default is 0.
;
;       MINOR:    The number of minor tick divisions. Default is 2.
;
;       NCOLORS:  This is the number of colors in the color bar.
;
;       POSITION: A four-element array of normalized coordinates in the same
;                 form as the POSITION keyword on a plot. Default is
;                 [0.88, 0.15, 0.95, 0.95] for a vertical bar and
;                 [0.15, 0.88, 0.95, 0.95] for a horizontal bar.
;;
;       RANGE:    A two-element vector of the form [min, max]. Provides an
;                 alternative way of setting the MINRANGE and MAXRANGE keywords.
;
;       RIGHT:    This puts the labels on the right-hand side of a vertical
;                 color bar. It applies only to vertical color bars.
;
;       TITLE:    This is title for the color bar. The default is to have
;                 no title.
;
;       TOP:      This puts the labels on top of the bar rather than under it.
;                 The keyword only applies if a horizontal color bar is rendered.
;
;       VERTICAL: Setting this keyword give a vertical color bar. The default
;                 is a horizontal color bar.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       Color bar is drawn in the current graphics window.
;
; RESTRICTIONS:
;       The number of colors available on the display device (not the
;       PostScript device) is used unless the NCOLORS keyword is used.
;
; EXAMPLE:
;       To display a horizontal color bar above a contour plot, type:
;
;       LOADCT, 5, NCOLORS=100
;       CONTOUR, DIST(31,41), POSITION=[0.15, 0.15, 0.95, 0.75], $
;          C_COLORS=INDGEN(25)*4, NLEVELS=25
;       COLORBAR_CO, NCOLORS=100, POSITION=[0.15, 0.85, 0.95, 0.90]
;
; MODIFICATION HISTORY:
;       Written by: David Fanning, 10 JUNE 96.
;       10/27/96: Added the ability to send output to PostScript. DWF
;       11/4/96: Substantially rewritten to go to screen or PostScript
;           file without having to know much about the PostScript device
;           or even what the current graphics device is. DWF
;       1/27/97: Added the RIGHT and TOP keywords. Also modified the
;            way the TITLE keyword works. DWF
;       7/15/97: Fixed a problem some machines have with plots that have
;            no valid data range in them. DWF
;       12/5/98: Fixed a problem in how the colorbar image is created that
;            seemed to tickle a bug in some versions of IDL. DWF.
;       1/12/99: Fixed a problem caused by RSI fixing a bug in IDL 5.2. Sigh... DWF.
;       3/30/99: Modified a few of the defaults. DWF.
;       3/30/99: Used NORMAL rather than DEVICE coords for positioning bar. DWF.
;       3/30/99: Added the RANGE keyword. DWF.
;       3/30/99: Added FONT keyword. DWF
;       5/6/99: Many modifications to defaults. DWF.
;       5/6/99: Removed PSCOLOR keyword. DWF.
;       5/6/99: Improved error handling on position coordinates. DWF.
;       5/6/99. Added MINOR keyword. DWF.
;       5/6/99: Set Device, Decomposed=0 if necessary. DWF.
;       2/9/99: Fixed a problem caused by setting BOTTOM keyword, but not NCOLORS. DWF.
;       8/17/99. Fixed a problem with ambiguous MIN and MINOR keywords. DWF
;       8/25/99. I think I *finally* got the BOTTOM/NCOLORS thing sorted out. :-( DWF.
;       10/10/99. Modified the program so that current plot and map coordinates are
;            saved and restored after the colorbar is drawn. DWF.
;       3/18/00. Moved a block of code to prevent a problem with color decomposition. DWF.
;       4/28/00. Made !P.Font default value for FONT keyword. DWF.
;       9/26/00. Made the code more general for scalable pixel devices. DWF.
;       1/16/01. Added INVERTCOLORS keyword. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2000 Fanning Software Consulting.
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################


PRO COLORBAR_CO, BOTTOM=bottom, CHARSIZE=charsize, COLOR=color, DIVISIONS=divisions, $
   FORMAT=format, POSITION=position, MAXRANGE=maxrange, MINRANGE=minrange, NCOLORS=ncolors, $
   TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right, MINOR=minor, $
   RANGE=range, FONT=font, TICKLEN=ticklen, _EXTRA=extra, INVERTCOLORS=invertcolors

   ; Return to caller on error.

On_Error, 2

   ; Save the current plot state.

bang_p = !P
bang_x = !X
bang_Y = !Y
bang_Z = !Z
bang_Map = !Map

   ; Are scalable pixels available on the device?

IF (!D.Flags AND 1) NE 0 THEN scalablePixels = 1 ELSE scalablePixels = 0

   ; Which release of IDL is this?

thisRelease = Float(!Version.Release)

    ; Check and define keywords.

IF N_ELEMENTS(ncolors) EQ 0 THEN BEGIN

   ; Most display devices to not use the 256 colors available to
   ; the PostScript device. This presents a problem when writing
   ; general-purpose programs that can be output to the display or
   ; to the PostScript device. This problem is especially bothersome
   ; if you don't specify the number of colors you are using in the
   ; program. One way to work around this problem is to make the
   ; default number of colors the same for the display device and for
   ; the PostScript device. Then, the colors you see in PostScript are
   ; identical to the colors you see on your display. Here is one way to
   ; do it.

   IF scalablePixels THEN BEGIN
      oldDevice = !D.NAME

         ; What kind of computer are we using? SET_PLOT to appropriate
         ; display device.

      thisOS = !VERSION.OS_FAMILY
      thisOS = STRMID(thisOS, 0, 3)
      thisOS = STRUPCASE(thisOS)
      CASE thisOS of
         'MAC': SET_PLOT, thisOS
         'WIN': SET_PLOT, thisOS
         ELSE: SET_PLOT, 'X'
      ENDCASE

         ; Here is how many colors we should use.

      ncolors = !D.TABLE_SIZE
      SET_PLOT, oldDevice
    ENDIF ELSE ncolors = !D.TABLE_SIZE
ENDIF
IF N_ELEMENTS(bottom) EQ 0 THEN bottom = 0B
IF N_ELEMENTS(charsize) EQ 0 THEN charsize = 1.0
IF N_ELEMENTS(format) EQ 0 THEN format = '(I5)'
IF N_ELEMENTS(color) EQ 0 THEN color = !P.Color
IF N_ELEMENTS(minrange) EQ 0 THEN minrange = 0
IF N_ELEMENTS(maxrange) EQ 0 THEN maxrange = ncolors
IF N_ELEMENTS(ticklen) EQ 0 THEN ticklen = 0.2
IF N_ELEMENTS(minor) EQ 0 THEN minor = 2
IF N_ELEMENTS(range) NE 0 THEN BEGIN
   minrange = range[0]
   maxrange = range[1]
ENDIF
IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 6
IF N_ELEMENTS(font) EQ 0 THEN font = !P.Font
IF N_ELEMENTS(title) EQ 0 THEN title = ''

IF KEYWORD_SET(vertical) THEN BEGIN
   bar = REPLICATE(1B,20) # BINDGEN(ncolors)
   IF Keyword_Set(invertcolors) THEN bar = Reverse(bar, 2)
   IF N_ELEMENTS(position) EQ 0 THEN BEGIN
      position = [0.88, 0.1, 0.95, 0.9]
   ENDIF ELSE BEGIN
      IF position[2]-position[0] GT position[3]-position[1] THEN BEGIN
         position = [position[1], position[0], position[3], position[2]]
      ENDIF
      IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
      IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
   ENDELSE
ENDIF ELSE BEGIN
   bar = BINDGEN(ncolors) # REPLICATE(1B, 20)
   IF Keyword_Set(invertcolors) THEN bar = Reverse(bar, 1)
   IF N_ELEMENTS(position) EQ 0 THEN BEGIN
      position = [0.1, 0.88, 0.9, 0.95]
   ENDIF ELSE BEGIN
      IF position[3]-position[1] GT position[2]-position[0] THEN BEGIN
         position = [position[1], position[0], position[3], position[2]]
      ENDIF
      IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
      IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
   ENDELSE
ENDELSE

   ; Scale the color bar.

 bar = BYTSCL(bar, TOP=(ncolors-1 < (255-bottom))) + bottom

   ; Get starting locations in NORMAL coordinates.

xstart = position(0)
ystart = position(1)

   ; Get the size of the bar in NORMAL coordinates.

xsize = (position(2) - position(0))
ysize = (position(3) - position(1))

   ; Alternative using CONTOUR

;bar=grange(range(0),range(1),60) # REPLICATE(1., 2)
;contour, bar, position=[xstart, ystart, xstart+xsize, ystart+ysize], nlev=255, /noerase,$
;         /fill, xticks=1, yticks=1, xminor=1, yminor=1,xtickname=[' ',' '],ytickname=[' ',' ']

   ; Display the color bar in the window. Sizing is
   ; different for PostScript and regular display.

IF scalablePixels THEN BEGIN
  
   TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize, /Normal

ENDIF ELSE BEGIN

   bar = CONGRID(bar, CEIL(xsize*!D.X_VSize), CEIL(ysize*!D.Y_VSize), /INTERP)

        ; Decomposed color off if device supports it.

   CASE  StrUpCase(!D.NAME) OF
        'X': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        'WIN': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        'MAC': BEGIN
            IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
            Device, Decomposed=0
            ENDCASE
        ELSE:
   ENDCASE

   TV, bar, xstart, ystart, /Normal

      ; Restore Decomposed state if necessary.

   CASE StrUpCase(!D.NAME) OF
      'X': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      'WIN': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      'MAC': BEGIN
         IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
         ENDCASE
      ELSE:
   ENDCASE

ENDELSE

   ; Annotate the color bar.

IF KEYWORD_SET(vertical) THEN BEGIN

   IF KEYWORD_SET(right) THEN BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', YTICKLEN=ticklen , $
         YRANGE=[minrange, maxrange], FONT=font, _EXTRA=extra, YMINOR=minor

      AXIS, YAXIS=1, YRANGE=[minrange, maxrange], YTICKFORMAT=format, YTICKS=divisions, $
         YTICKLEN=ticklen, YSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         FONT=font, YTITLE=title, _EXTRA=extra, YMINOR=minor

   ENDIF ELSE BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, YMINOR=minor, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT=format, XTICKFORMAT='(A1)', YTICKLEN=ticklen , $
         YRANGE=[minrange, maxrange], FONT=font, YTITLE=title, _EXTRA=extra

      AXIS, YAXIS=1, YRANGE=[minrange, maxrange], YTICKFORMAT='(A1)', YTICKS=divisions, $
         YTICKLEN=ticklen, YSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         FONT=font, _EXTRA=extra, YMINOR=minor

   ENDELSE

ENDIF ELSE BEGIN

   IF KEYWORD_SET(top) THEN BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=9, YSTYLE=1, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', XTICKLEN=ticklen, $
         XRANGE=[minrange, maxrange], FONT=font, _EXTRA=extra, XMINOR=minor

      AXIS, XTICKS=divisions, XSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         XTICKFORMAT=format, XTICKLEN=ticklen, XRANGE=[minrange, maxrange], XAXIS=1, $
         FONT=font, XTITLE=title, _EXTRA=extra, XCHARSIZE=charsize, XMINOR=minor

   ENDIF ELSE BEGIN

      PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=1, YSTYLE=1, TITLE=title, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT=format, XTICKLEN=ticklen, $
         XRANGE=[minrange, maxrange], FONT=font, XMinor=minor, _EXTRA=extra

    ENDELSE

ENDELSE

   ; Restore the previous plot and map system variables.

!P = bang_p
!X = bang_x
!Y = bang_y
!Z = bang_z
!Map = bang_map

END
