;
;  $Id:$
;
; This script will visualise the flow of particles around one or more
; solid objects.  
;
;  Author: Nils Erland L. Haugen
;
FUNCTION CIRCLE_, xcenter, ycenter, radius
points = (2 * !PI / 99.0) * FINDGEN(100)
x = xcenter + radius * COS(points )
y = ycenter + radius * SIN(points )
RETURN, TRANSPOSE([[x],[y]])
END
;
pro pc_visualize_particles,png=png,removed=removed
;
device,decompose=0
loadct,5
;
; Set defaults
;
default,writepng,0
default,removed,0
;
; Read dimensions and namelists
;
pc_read_dim,obj=procdim
pc_read_param, object=param
pc_read_pvar,obj=objpvar,/solid_object,irmv=irmv,theta_arr=theta_arr
pc_read_pstalk,obj=obj
;
; Set some auxillary variables
;
nx=procdim.nx
ny=procdim.ny
dims=size(obj.xp)
n_parts=dims(1)
n_steps=dims(2)
radius=param.cylinder_radius[0]
xpos=param.cylinder_xpos[0]
ypos=param.cylinder_ypos[0]
;
; Find positions of removed particles (to be used later for plotting them).
;
if (removed eq 1) then begin
    removed_pos=objpvar.xx(irmv,*)
    skin=nx
    collision_radius=sqrt((removed_pos(*,0)-xpos)^2+(removed_pos(*,1)-ypos)^2)
    solid_colls=where(collision_radius lt radius+1e-3)
endif
;
; Find where (in radians) the particles hit the surface of the cylinder as a
; function of time
;
theta_=theta_arr[*,0]
time_=theta_arr[*,1]
here=where(theta_ ne 0)
print,here
if (here[0] ne -1) then begin
    WINDOW,4,XSIZE=128*2,YSIZE=256*2
    theta=theta_[here]
    timereal=time_[here]
    time=timereal-min(timereal)
    dims=size(theta)
    ind=indgen(dims[1])
    !x.range=[0,max(ind)]
    !x.range=[min(time),max(time)]
    !y.range=[min(theta),max(theta)]
    plot,time,theta,ps=2,ytit='!4h!6',xtit='time'
    print,'The first particle hit the surface at t=',min(timereal)
    print,'The last particle hit the surface at t =',max(timereal)
    save,time,theta,filename='./data/theta.sav'
endif else begin
    print,'No particles has hit the cylinder surface!'
endelse
print,'The final time of the simulation is  t =',objpvar.t
;
; Set window size
;
xr=param.xyz1[0]-param.xyz0[0]
yr=param.xyz1[1]-param.xyz0[1]
WINDOW,3,XSIZE=1024*xr/yr*1.6,YSIZE=1024
!x.range=[param.xyz0[0],param.xyz1[0]]
!y.range=[param.xyz0[1],param.xyz1[1]]
;
; Show results
;
for i=0,n_steps-1 do begin
    plot,obj.xp(*,i),obj.yp(*,i),psym=sym(1),symsize=1,/iso
    POLYFILL, CIRCLE_(xpos, ypos, radius),color=122
    ;
    ; Do we want to write png files or to show results on screen
    ;
    if writepng eq 1 then begin
        istr2 = strtrim(string(i,'(I20.4)'),2) ;(only up to 9999 frames)
        file='img_'+istr2+'.png'
        write_png,file,tvrd()
    endif else begin
        wait,0.03
    endelse
end
;
; Plot the removed particles as blue dots
;
if (removed eq 1) then begin
    oplot,removed_pos(solid_colls,0),removed_pos(solid_colls,1),col=45,ps=sym(1)
endif
;
; Write png files if required
;
if writepng eq 1 then begin
    istr2 = strtrim(string(i+1,'(I20.4)'),2) ;(only up to 9999 frames)
    file='img_'+istr2+'.png'
    write_png,file,tvrd()
endif
;
END
