; $Id$
;
; Copyright (c) 1997-1999, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.
;+
; NAME:
;	vizit_trackball
;
; PURPOSE:
;	This object translates widget events for draw widgets into
;       transformations that emulate a virtual trackball (for transforming
;       object graphics in three dimensions).
;
; CATEGORY:
;	Object Graphics.
;
; CALLING SEQUENCE:
;	To initially create:
;		oTrackball = OBJ_NEW('vizit_trackball', Center, Radius)
;
;	To update the trackball state based on a widget event:
;		oTrackball->Update, sEvent
;
;	To re-initialize the trackball state:
;		oTrackball->Reset, Center, Radius
;
;	To destroy:
;		OBJ_DESTROY, oTrackball
;
; INPUTS:
;   TRACKBALL::INIT:
;	Center:	A two-dimensional vector, [x,y], representing the requested
;		center (measured in device units) of the trackball.
;       Radius:	The requested radius (measured in device units) of the
; 		trackball.
;
;   TRACKBALL::UPDATE:
;        sEvent: The widget event structure.  The event type indicates
;		how the trackball state should be updated.
;
;   TRACKBALL::RESET:
;	Center:	A two-dimensional vector, [x,y], representing the requested
;		center (measured in device units) of the trackball.
;       Radius:	The requested radius (measured in device units) of the
; 		trackball.
;
; KEYWORD PARAMETERS:
;   TRACKBALL::INIT:
;	AXIS:		Set this keyword to indicate the axis about which
;			rotations are to be constrained if the CONSTRAIN
;			keyword is set to a nonzer value.  Valid values
;			include:
;				0 = X-Axis
;				1 = Y-Axis
;				2 = Z-Axis (default)
;	CONSTRAIN:	Set this keyword to a nonzero value to indicate that
;			the trackball transformations are to be constrained
;			about a given axis (as specified by the AXIS
;			keyword).  The default is zero (no constraints).
;	MOUSE		Silently ignored, for compatibilty with IDL's Trackball
;
;   TRACKBALL::UPDATE:
;	MOUSE:	Silently ignored, for compatibilty with IDL's Trackball
;
;	TRANSFORM:	Set this keyword to a named variable that upon
;			return will contain a floating point 4x4 array
;			if a transformations matrix is calculated as
;			a result of the widget event.
;	TRANSLATE:	Set this keyword to indicate that the trackball
;			movement should be constrained to x and y translation
;			rather than rotation about an axis.
;
;   TRACKBALL::RESET:
;	AXIS:		Set this keyword to indicate the axis about which
;			rotations are to be constrained if the CONSTRAIN
;			keyword is set to a nonzer value.  Valid values
;			include:
;				0 = X-Axis
;				1 = Y-Axis
;				2 = Z-Axis (default)
;	CONSTRAIN:	Set this keyword to a nonzero value to indicate that
;			the trackball transformations are to be constrained
;			about a given axis (as specified by the AXIS
;			keyword).  The default is zero (no constraints).
;	MOUSE:	Silently ignored, for compatibilty with IDL's Trackball
;
; OUTPUTS:
;   TRACKBALL::UPDATE:
;	This function returns a 1 if a transformation matrix is calculated
;	as a result of the widget event, or 0 otherwise.
;
; EXAMPLE:
;	Create a trackball centered on a 512x512 pixel drawable area, and
;	a view containing the model to be manipulated:
;		xdim = 512
;		ydim = 512
;		wBase = WIDGET_BASE()
;		wDraw = WIDGET_DRAW(wBase, XSIZE=xdim, YSIZE=ydim, $
;		                    GRAPHICS_LEVEL=2, /BUTTON_EVENTS, $
;		                    /MOTION_EVENTS, /EXPOSE_EVENTS, RETAIN=0 )
;		WIDGET_CONTROL, wBase, /REALIZE
;		WIDGET_CONTROL, wDraw, GET_VALUE=oWindow
;
;		oTrackball = OBJ_NEW('Trackball', [xdim/2.,ydim/2.], xdim/2.)
;		oView = OBJ_NEW('IDLgrView')
;		oModel = OBJ_NEW('IDLgrModel')
;		oView->Add, oModel
;
;		XMANAGER, 'TrackEx', wBase
;
;	In the widget event handler, handle trackball updates.
;	As the trackball transformation changes, update the transformation
;	for a model object (instance of IDLgrModel), and redraw the view:
;
;	PRO TrackEx_Event, sEvent
;		...
;		bHaveXform = oTrackball->Update( sEvent, TRANSFORM=TrackXform )
;		IF (bHaveXform) THEN BEGIN
;		    oModel->GetProperty, TRANSFORM=ModelXform
;		    oModel->SetProperty, TRANSFORM=ModelXform # TrackXform
;		    oWindow->Draw, oView
;		ENDIF
;		...
;	END
;
; MODIFICATION HISTORY:
; 	Written by:	DD, December 1996
; 	Modified by: AJWM, September 2006
;      Added the possibility to specify the behaviour of
;      the buttons independently, including new scale function
;
;-

;----------------------------------------------------------------------------
; vizit_trackball::TRACKBALL_CONSTRAIN
;
; Purpose:
;  Given a point and a constraint vector, map the point to its constrained
;  equivalent.
;
; Arguments:
;  pt - The unconstrained point.
;  vec - A three-element vector, [x,y,z], representing the unit vector about
;        which rotations are constrained.
;
FUNCTION vizit_trackball::TRACKBALL_CONSTRAIN, pt, vec
    ; Project the point.
    proj = pt - TOTAL(vec * pt) * vec

    ; Normalizing factor.
    norm = SQRT(TOTAL(proj^2))

    IF (norm GT 0.0) THEN BEGIN
        IF (proj[2] LT 0.0) THEN $
            cpt = -1.0 / norm * proj $
        ELSE $
            cpt = 1.0 / norm * proj
    ENDIF ELSE IF vec[2] EQ 1.0 THEN $
        cpt = [1.0, 0.0, 0.0] $
    ELSE $
        cpt = [-vec[1], vec[0], 0.0] / SQRT(TOTAL(vec[0:1]^2))

    RETURN, cpt
END

;----------------------------------------------------------------------------
; vizit_trackball::motion_rotate_grab
;
; Purpose:
;  Calculate a transformation matrix based upon a pair of (x,y) coordinates
;  This function handles panning of an object parallel the the screen
;
;  The return value is nonzero if a transformation matrix is calculated
;  as a result of the event, or zero otherwise.
;
; Arguments:
;  pt0 - Old coords
;  pt1 - New coords
;
FUNCTION vizit_trackball::motion_rotate_grab, pt0, pt1
  ; translation only
  ; Compute transformation.
    q = [CROSSP(pt0,pt1), TOTAL(pt0*pt1)]

    x = q[0]
    y = q[1]
    z = q[2]
    w = q[3]

    return, [[ w^2+x^2-y^2-z^2, 2*(x*y-w*z), 2*(x*z+w*y), 0], $
             [ 2*(x*y+w*z), w^2-x^2+y^2-z^2, 2*(y*z-w*x), 0], $
             [ 2*(x*z-w*y), 2*(y*z+w*x), w^2-x^2-y^2+z^2, 0], $
             [ 0          , 0          , 0              , 1]]
END

;----------------------------------------------------------------------------
; vizit_trackball::motion_scale_grab
;
; Purpose:
;  Calculate a transformation matrix based upon a pair of (x,y) coordinates
;  This function handles scaling of an object as if moving in and out of
;  the screen
;
; Arguments:
;  pt0 - Old coords
;  pt1 - New coords
;
FUNCTION vizit_trackball::motion_scale_grab, pt0, pt1
  ; translation only
    s = exp(-(pt1[1]-pt0[1]))
    return, [[ s, 0, 0,             0 ], $
             [ 0, s, 0,             0 ], $
             [ 0, 0, s,             0 ], $
             [ 0, 0, 0,             1 ]  ]
END

;----------------------------------------------------------------------------
; vizit_trackball::motion_translate_grab
;
; Purpose:
;  Calculate a transformation matrix based upon a pair of (x,y) coordinates
;  This function handles panning of an object parallel the the screen
;
;  The return value is nonzero if a transformation matrix is calculated
;  as a result of the event, or zero otherwise.
;
; Arguments:
;  pt0 - Old coords
;  pt1 - New coords
;
FUNCTION vizit_trackball::motion_translate_grab, pt0, pt1
  ; translation only
   return, [[ 1, 0, 0, pt1[0]-pt0[0] ], $
            [ 0, 1, 0, pt1[1]-pt0[1] ], $
            [ 0, 0, 1,             0 ], $
            [ 0, 0, 0,             1 ]  ]
END

;----------------------------------------------------------------------------
; vizit_trackball::motion_translatez_grab
;
; Purpose:
;  Calculate a transformation matrix based upon a pair of (x,y) coordinates
;  This function handles panning of an object parallel the the screen
;
;  The return value is nonzero if a transformation matrix is calculated
;  as a result of the event, or zero otherwise.
;
; Arguments:
;  pt0 - Old coords
;  pt1 - New coords
;
FUNCTION vizit_trackball::motion_translatez_grab, pt0, pt1
  ; translation only
   return, [[ 1, 0, 0,             0 ], $
            [ 0, 1, 0,             0 ], $
            [ 0, 0, 1, pt1[1]-pt0[1] ], $
            [ 0, 0, 0,             1 ]  ]
END

;----------------------------------------------------------------------------
; TRACKBALL::UPDATE
;
; Purpose:
;  Given a widget event structure, updates the trackball state.
;
;  The return value is nonzero if a transformation matrix is calculated
;  as a result of the event, or zero otherwise.
;
; Arguments:
;  sEvent - The widget event structure.
;
; Keywords:
;  TRANSFORM - If a transformation matrix is calculated, upon return
;              transform will contain a 4x4 matrix.
;
FUNCTION vizit_trackball::UPDATE, sEvent, TRANSFORM=transform, $
     MOUSE=mouse, TRANSLATE=translate
  ; Initialize return value.
  bHaveTransform=0

  IF (N_ELEMENTS(mouse) NE 0) THEN BEGIN
      if (mouse ne 1) and (mouse ne 2) and (mouse ne 4) then begin
          PRINT, 'Trackball: invalid value for MOUSE keyword.'
          RETURN, 0
      ENDIF ELSE $
          self.mouse = mouse
  ENDIF ; Don't set self.mouse if no argument, keep setting from INIT

  ; Ignore non-Draw-Widget events.
  IF (TAG_NAMES(sEvent, /STRUCTURE_NAME) NE 'WIDGET_DRAW') THEN $
    RETURN, bHaveTransform
  ; Determine event type.

  CASE sEvent.type OF
    0: BEGIN    ;Button press.
         ; Only handle event if appropriate mouse button.
         ;aa IF (sEvent.press EQ self.mouse) THEN BEGIN
           self.mouse = sEvent.press

           ; Calculate distance of mouse click from center of unit circle.
           xy = ([sEvent.x,sEvent.y] - self.center) / self.radius
           r = TOTAL(xy^2)
           IF (r GT 1.0) THEN $
             self.pt1 = [xy/SQRT(r) ,0.0] $
           ELSE $
             self.pt1 = [xy,SQRT(1.0-r)]

           ; Constrain if necessary.
           IF (self.constrain NE 0) THEN BEGIN
               vec = [0.,0.,0.]
               vec[self.axis] = 1.0
               self.pt1 = self->TRACKBALL_CONSTRAIN(self.pt1, vec)
           ENDIF
           self.pt0 = self.pt1
           self.btndown = 1b
         ;aa ENDIF
       END

    2: BEGIN    ;Button motion.
         IF (self.btndown EQ 1b) THEN BEGIN
           ;aa print, sEvent.type, sEvent.press, sEvent.release, sEvent.x
           ; Calculate distance of mouse click from center of unit circle.
           xy = ([sEvent.x,sEvent.y] - self.center) / $
                self.radius
           r = TOTAL(xy^2)
           IF (r GT 1.0) THEN $
             pt1 = [xy/SQRT(r) ,0.0] $
           ELSE $
             pt1 = [xy,SQRT(1.0-r)]

           ; Constrain if necessary.
           IF (self.constrain NE 0) THEN BEGIN
               vec = [0.,0.,0.]
               vec[self.axis] = 1.0
               pt1 = self->TRACKBALL_CONSTRAIN( pt1, vec)
           ENDIF

           ; Update the transform only if the mouse button has actually
           ; moved from its previous location.
           pt0 = self.pt0
           IF ((pt0[0] NE pt1[0]) OR $
               (pt0[1] NE pt1[1]) OR $
               (pt0[2] NE pt1[2])) THEN BEGIN

             bHaveTransform = 0b

             ; 1: begin		; left button
             ; 2: begin		; middle button
             ; 4: begin		; right button
             ; 5: begin		; both buttons
             if (self.mouse lt n_elements(self.button_functions)) then begin
                func=self.button_functions[self.mouse]
                if (func ne '') then begin
                  bHaveTransform = 1b
                  res = execute('transform=self->motion_'+func+'(pt0,pt1)')
                endif
             endif

             self.pt0 = pt1
           ENDIF

           self.pt1 = pt1
         ENDIF
       END

    1: BEGIN    ;Button Release.
         IF (self.btndown EQ 1b) THEN $
           self.btndown = 0b
       END

    ELSE: RETURN, bHaveTransform

   ENDCASE

  RETURN, bHaveTransform
END

;----------------------------------------------------------------------------
; vizit_trackball::INIT
;
; Purpose:
;   Initializes the trackball state.
;
; Arguments:
;  center - A two-dimensional vector, [x,y], representing the requested
;           center (measured in device units) of the trackball.
;  radius - The requested radius (measured in device units) of the trackball.
;
; Keywords:
;  AXIS -  Set this keyword to indicate the axis about which rotations
;          are to be constrained if the CONSTRAIN keyword is set to a
;          nonzero value).  Valid values include:
;               0 = X
;               1 = Y
;               2 = Z (default)
;  CONSTRAIN - Set this keyword to a nonzero value to indicate that the
;          trackball transformations are to be constrained about a given
;          axis (as specified by the AXIS keyword).  The default is zero
;          (no constraints).
;  MOUSE - Set this keyword to a bitmask to indicate which mouse button to
;          honor for trackball events.  The least significant bit represents
;          the leftmost button, the next highest bit represents the middle
;          button, and the next highest bit represents the right button.
;          The default is 1b, for the left moust button.
;
FUNCTION vizit_trackball::INIT, center, radius, AXIS=axis, CONSTRAIN=constrain, $
                          MOUSE=mouse,FUNCTIONS=functions
    IF (N_ELEMENTS(center) NE 2) THEN BEGIN
        PRINT, 'Trackball: center must be a two-dimensional array.'
        RETURN, 0
    ENDIF

    IF (N_ELEMENTS(radius) NE 1) THEN BEGIN
        PRINT, 'Trackball: invalid radius.'
        RETURN, 0
    ENDIF

    IF (N_ELEMENTS(axis) NE 0) THEN BEGIN
        IF (axis lt 0) OR (axis gt 2) THEN BEGIN
            PRINT, 'Trackball: invalid value for AXIS keyword.'
            RETURN, 0
        ENDIF ELSE $
            self.axis = axis
    ENDIF ELSE $
        self.axis = 2

    IF (N_ELEMENTS(constrain) NE 0) THEN $
        self.constrain = constrain $
    ELSE $
        self.constrain = 0

    IF (N_ELEMENTS(mouse) NE 0) THEN BEGIN
        if (mouse ne 1) and (mouse ne 2) and (mouse ne 4) then begin
            PRINT, 'Trackball: invalid value for MOUSE keyword.'
            RETURN, 0
        ENDIF ELSE $
            self.mouse = mouse
    ENDIF ELSE $
        self.mouse = 1 ; Default is left mouse button

    IF (N_ELEMENTS(functions) NE 0) THEN $
        self.button_functions=functions $
    ELSE $
        self.button_functions=['','rotate_grab','translate_grab','','scale_grab','translate_grab']

    self.center = center
    self.radius = radius
    self.btndown = 0b

    RETURN, 1
END

;----------------------------------------------------------------------------
; vizit_trackball::CLEANUP
;
; Purpose:
;   Cleanup the trackball when it is destroyed.
;
; Arguments:
;  <None>
;
; Keywords:
;  <None>
;
PRO vizit_trackball::CLEANUP
    ; No work needs to be done.  We provide this method to avoid any
    ; problems resolving the Cleanup method call on Windows 3.11
    ; which has short filename restrictions.
END

;----------------------------------------------------------------------------
; vizit_trackball::RESET
;
; Purpose:
;   Resets the trackball state.
;
; Arguments:
;  center - A two-dimensional vector, [x,y], representing the requested
;           center (measured in device units) of the trackball.
;  radius - The requested radius (measured in device units) of the trackball.
;
; Keywords:
;  MOUSE - Set this keyword to a bitmask to indicate which mouse button to
;          honor for trackball events.  The least significant bit represents
;          the leftmost button, the next highest bit represents the middle
;          button, and the next highest bit represents the right button.
;          The default is 1b, for the left moust button.
;
PRO vizit_trackball::RESET, center, radius, AXIS=axis, CONSTRAIN=constrain, $
                      MOUSE=mouse, FUNCTIONS=functions
    IF (N_ELEMENTS(center) NE 2) THEN BEGIN
        PRINT, 'vizit_trackball: center must be a two-dimensional array.'
        RETURN
    ENDIF

    IF (N_ELEMENTS(radius) NE 1) THEN BEGIN
        PRINT, 'vizit_trackball: invalid radius.'
        RETURN
    ENDIF

    IF (N_ELEMENTS(axis) NE 0) THEN BEGIN
        IF (axis lt 0) OR (axis gt 2) THEN BEGIN
            PRINT, 'vizit_trackball: invalid value for AXIS keyword.'
            RETURN
        ENDIF ELSE $
            self.axis = axis
    ENDIF ELSE $
        self.axis = 2

    IF (N_ELEMENTS(constrain) NE 0) THEN $
        self.constrain = constrain $
    ELSE $
        self.constrain = 0

    IF (N_ELEMENTS(functions) NE 0) THEN $
        self.button_functions=functions $
    ELSE $
        self.button_functions=['','rotate_grab','translate_grab','','scale_grab','translate_grab']

    IF (N_ELEMENTS(mouse) NE 0) THEN BEGIN
        if (mouse ne 1) and (mouse ne 2) and (mouse ne 4) then begin
            PRINT, 'vizit_trackball: invalid value for MOUSE keyword.'
            RETURN
        ENDIF ELSE $
            self.mouse = mouse
    ENDIF ELSE $
        self.mouse = 1 ; Default is left mouse button

    self.center = center
    self.radius = radius

    self.btndown = 0b

END

;----------------------------------------------------------------------------
; vizit_trackball__DEFINE
;
; Purpose:
;  Defines the object structure for a trackball object.
;
PRO vizit_trackball__define

  struct = {vizit_trackball, $
            btndown: 0b,    $
            axis: 0, $
            constrain: 0b, $
            mouse: 0b, $
            button_functions: ['','rotate_grab','translate_grab','','scale_grab','translate_grab'], $
            center: LONARR(2), $
            radius: 0.0, $
            pt0: FLTARR(3), $
            pt1: FLTARR(3) $
           }
END
