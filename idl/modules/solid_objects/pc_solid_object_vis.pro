FUNCTION CIRCLE_, xcenter, ycenter, radius
points = (2 * !PI / 99.0) * FINDGEN(100)
x = xcenter + radius * COS(points )
y = ycenter + radius * SIN(points )
RETURN, TRANSPOSE([[x],[y]])
END

pro pc_solid_object_vis,start,color=color
;
default,color,255
;
; Check if cylinders and spheres are defined
;
print,start.ncylinders
;
; Insert cylinders
;
if (has_tag (start, 'ncylinders')) then begin
    if (start.ncylinders > 0) then begin
        for icyl=0,start.ncylinders-1 do begin
            r=start.cylinder_radius[icyl]
            x0=start.cylinder_xpos[icyl]
            y0=start.cylinder_ypos[icyl]
            POLYFILL, CIRCLE_(x0,y0,r),color=color,/DATA
        end
    endif
endif
;
; Insert spheres
;
if (has_tag (start, 'nspheres')) then begin
    if (start.nspheres > 0) then begin
        for icyl=0,start.nspheres-1 do begin
            r=start.sphere_radius[icyl]
            x0=start.sphere_xpos[icyl]
            y0=start.sphere_ypos[icyl]
            POLYFILL, CIRCLE_(x0,y0,r),color=color,/DATA
        end
    endif
endif


END
