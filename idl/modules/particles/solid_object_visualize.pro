FUNCTION CIRCLE_, xcenter, ycenter, radius
points = (2 * !PI / 99.0) * FINDGEN(100)
x = xcenter + radius * COS(points )
y = ycenter + radius * SIN(points )
RETURN, TRANSPOSE([[x],[y]])
END

pro solid_object_visualize,ncylinders,n_steps,obj,tmin,tmax,startparam,xpos,$
                           ypos,radius,oneradius,png,w,psym,Stokes=Stokes,$
                           r_i=r_i,i=i
for i=0,n_steps-1 do begin
;
; Check if we want to plot for this time
;
if ((obj.t[i] gt tmin) and (obj.t[i] lt tmax)) then begin
    titlestring='t='+str(obj.t[i])
    if (startparam.coord_system eq 'cylindric') then begin
        xp=obj.xp(*,i)*cos(obj.yp(*,i))
        yp=obj.xp(*,i)*sin(obj.yp(*,i))
        plot,xp,yp,psym=psym,symsize=1,/iso,title=titlestring 
        POLYFILL, CIRCLE_(xpos, $
                          ypos, $
                          radius),color=122        
    endif else begin
        if (oneradius eq 0) then begin
            plot,obj.xp(*,i),obj.yp(*,i),psym=psym,symsize=1,/iso,$
              title=titlestring
        endif else begin
            titlestring=titlestring+',St='+String(Stokes[r_i],$
                                                  format='(F4.2)')
            rad_indices = where(float(obj.ap[*,i]) eq $
                                float(startparam.ap0[r_i]))
            n_particles = n_elements(rad_indices)
            if (i eq 50) then print, 'n_particles = ', n_particles
            if (n_particles gt 1) then begin
                plot,obj.xp(rad_indices,i),obj.yp(rad_indices,i),$
                  psym=psym,symsize=1,/iso,title=titlestring
            endif
        endelse
        for icyl=0,ncylinders-1 do begin
            POLYFILL, CIRCLE_(xpos[icyl], $
                              ypos[icyl], $
                              radius[icyl]),color=122
        end
    endelse
    ;
    ; Do we want to write png files or to show results on screen
    ;
    if png eq 1 then begin
        store_png_frame,i
    endif else begin
        wait,w
    endelse
endif
end
END
