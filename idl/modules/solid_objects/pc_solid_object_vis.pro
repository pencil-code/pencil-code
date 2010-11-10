FUNCTION CIRCLE_, xcenter, ycenter, radius
points = (2 * !PI / 99.0) * FINDGEN(100)
x = xcenter + radius * COS(points )
y = ycenter + radius * SIN(points )
RETURN, TRANSPOSE([[x],[y]])
END

pro pc_solid_object_vis,field
;device,decompose=0
;loadct,5
  ;; for icyl=0,ncylinders-1 do begin
  ;;    POLYFILL, CIRCLE_(xpos[icyl], $
  ;;                      ypos[icyl], $
  ;;                      radius[icyl]),color=122
  ;; end
;pc_read_var,obj=obj
;contour,obj.uu(*,*,3,0),obj.x,obj.y,/fill,nlev=60,/iso
POLYFILL, CIRCLE_(0.2,0.0,0.04),color=255,/DATA

END
