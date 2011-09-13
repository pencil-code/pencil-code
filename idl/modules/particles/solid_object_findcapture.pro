;
; Finds deposition and capture efficiency on solid_object (only
; cylinder so far).
;
; Procedures:
; <solid_object_findcapture>: main procedure, communicates with
;       pc_visualize_particles
; <cylinder_particle_deposition>: finds amount of particles deposited
;       on cylinder(s), as well as their position
; <cylinder_capture_efficiency>: calculates capture efficiency of
;       cylinder
; <cylinder_find_removed_positions>: finds removed positions to be
;       used later for plotting them
; <cylinder_plot_deposition>: Plots theta(t) and deposition
;       distribution between particle radii
;
; Authors: Nils E. Haugen and Anders Granskogen Bjørnstad
;
;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; main procedure, communicates with  pc_visualize_particles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro solid_object_findcapture,ncylinders=ncylinders,npart_radii=npart_radii,$
                             startparam=startparam,runparam=runparam,$
              linsert_particles_continuously=linsert_particles_continuously,$
                             objpvar=objpvar,irmv=irmv,trmv=trmv,xdir=xdir,$
                             radii_arr=radii_arr,savefile=savefile,$
                             removed=removed,oneradius=oneradius,$
                             radius=radius,xpos=xpos,ypos=ypos,rmv_pos=rmv_pos,$
                             on_cylinder_indices=on_cylinder_indices,r_i=r_i,$
                             material_density=material_density
;

print_remove_data=0
pc_read_ts,obj=ts
if (startparam.coord_system eq 'cylindric') then begin
    init_uu=startparam.ampluu
endif else begin
    init_uu=startparam.init_uu
endelse

; <cylinder_particle_deposition>: finds amount of particles deposited
;       on cylinder(s), as well as their position
cylinder_particle_deposition,ncylinders,npart_radii,startparam,objpvar,irmv,$
  trmv,xdir,xpos,ypos,radius,init_uu,print_remove_data,solid_colls=solid_colls,$
  front_colls=front_colls,back_colls=back_colls,theta_arr=theta_arr,$
  impact_vel_arr=impact_vel_arr

; <cylinder_capture_efficiency>: calculates capture efficiency of
;       cylinder
cylinder_capture_efficiency,ncylinders,linsert_particles_continuously,ts,$
  objpvar,startparam,runparam,init_uu,radius,radii_arr,xdir,npart_radii,$
  solid_colls,front_colls,back_colls,savefile,material_density

dims=size(irmv)
if ((removed eq 1) and (dims[0] ne 0) ) then begin
    ; <cylinder_find_removed_positions>: finds removed positions to be
    ;       used later for plotting them
    cylinder_find_removed_positions,ncylinders,objpvar,irmv,xpos,ypos,radius,$
      rmv_pos=rmv_pos, on_cylinder_indices=on_cylinder_indices
endif

; <cylinder_plot_deposition>: Plots theta(t) and deposition
;       distribution between particle radii
cylinder_plot_deposition,ts,oneradius,npart_radii,objpvar,theta_arr,savefile,$
  startparam,r_i=r_i,impact_vel_arr=impact_vel_arr

END
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; finds amount of particles deposited
;       on cylinder(s), as well as their position
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cylinder_particle_deposition,ncylinders,npart_radii,startparam,objpvar,irmv,$
                                 trmv,xdir,xpos,ypos,radius,init_uu,$
                                 print_remove_data,solid_colls=solid_colls,$
                                 front_colls=front_colls,back_colls=back_colls,$
                                 theta_arr=theta_arr,$
                                 impact_vel_arr=impact_vel_arr

dims=size(irmv)
solid_colls=fltarr(ncylinders,npart_radii)
solid_colls[*,*]=0
front_colls=solid_colls
back_colls=solid_colls
theta_arr=fltarr(ncylinders,npart_radii,1400000,2)
impact_vel_arr=fltarr(ncylinders,npart_radii,1400000,4)
;
;Loop over all cylinders
;
for icyl=0,ncylinders-1 do begin
    ;
    ; Loop over all removed particles
    ;
    for k=long(0),long(dims[1])-1 do begin        
        if (dims[0]>0) then begin              
            x0=objpvar.xx[irmv[k],0]-xpos[icyl]
            y0=objpvar.xx[irmv[k],1]-ypos[icyl]
            vx=objpvar.vv[irmv[k],0]
            vy=objpvar.vv[irmv[k],1]
            vz=objpvar.vv[irmv[k],2]
            deposition_radius2=x0^2+y0^2
            deposition_radius=sqrt(deposition_radius2)
            ;
            ; Check if the solid collision happened in 
            ; the vicinity of the
            ; cylinder. Here this is defined as being closer than 
            ; 0.2 radii away from the surface.
            ;
            if (deposition_radius lt radius[icyl]*1.2) then begin
                ;
                ; Find the radius of the particle
                ;
                ipart_radii=0
                small_number=1e-8
                while ((objpvar.a[irmv[k]] gt $
                        (startparam.ap0(ipart_radii)+small_number)) or $
                       (objpvar.a[irmv[k]] lt $
                        (startparam.ap0(ipart_radii)-small_number))) do begin
                    ipart_radii=ipart_radii+1
                end
                ;
                ; Find the angle of the impaction
                ;
                if (xdir) then begin
                    theta_tmp=acos(x0/deposition_radius)
                endif else begin
                    theta_tmp=acos(y0/deposition_radius)
                end
                theta=3.1415-theta_tmp
                if (print_remove_data) then begin
                    print,'time,k,r,x,y,theta=',$
                      trmv[k],irmv[k],deposition_radius,x0,y0,theta
                endif
                ;
                ; Store the angle of impaction in an array. 
                ;
                if (total(solid_colls[icyl]) lt 100000) then begin
                    theta_arr[icyl,ipart_radii,$
                              solid_colls[icyl,ipart_radii],0]=theta
                    theta_arr[icyl,ipart_radii,$
                              solid_colls[icyl,ipart_radii],1]=trmv[k]
                endif
                ;
                ; Store the velocity at impaction in an array
                ;
                if (total(solid_colls[icyl]) lt 100000) then begin
                    impact_vel_arr[icyl,ipart_radii,$
                              solid_colls[icyl,ipart_radii],0]=trmv[k]
                    impact_vel_arr[icyl,ipart_radii,$
                              solid_colls[icyl,ipart_radii],1]=vx
                    impact_vel_arr[icyl,ipart_radii,$
                              solid_colls[icyl,ipart_radii],2]=vy
                    impact_vel_arr[icyl,ipart_radii,$
                              solid_colls[icyl,ipart_radii],3]=vz
                endif
                ;
                ;Find the number of solid collisions
                ;for a given cylinder 
                ;and a given particle radius.
                ;
                solid_colls[icyl,ipart_radii]= $
                  solid_colls[icyl,ipart_radii]+1
                ;
                ; Find the number of back side and front side collisions.
                ;
                if (xdir) then begin
                    if (objpvar.xx[irmv[k],0] gt xpos[icyl]) then begin
                        back_colls[icyl,ipart_radii]=$
                          back_colls[icyl,ipart_radii]+1
                    endif else begin
                        front_colls[icyl,ipart_radii]=$
                          front_colls[icyl,ipart_radii]+1
                    endelse
                endif else begin
                    if (objpvar.xx[irmv[k],1] gt ypos[icyl]) then begin
                        back_colls[icyl,ipart_radii]=$
                          back_colls[icyl,ipart_radii]+1
                    endif else begin
                        front_colls[icyl,ipart_radii]=$
                          front_colls[icyl,ipart_radii]+1
                    endelse
                endelse
            endif
        endif
    endfor
endfor

END
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  calculates capture efficiency of
;       cylinder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cylinder_capture_efficiency,ncylinders,linsert_particles_continuously,ts,$
                                objpvar,startparam,runparam,init_uu,radius,$
                                radii_arr,xdir,npart_radii,solid_colls,$
                                front_colls,back_colls,savefile,$
                                material_density
; 
; Find how many particles have been inserted
;
if (linsert_particles_continuously) then begin
    initial_time=ts.t[0]
    final_time=min([objpvar.t,runparam.max_particle_insert_time])
    npar_inserted=(final_time-initial_time)*runparam.particles_insert_rate
    dimsp=size(where(objpvar.a gt 0))
    npar_inserted=dimsp[1]
endif else begin
   pc_read_pdim, npar=npar
   npar_inserted=npar
endelse
print,'Total number of inserted particles:',npar_inserted
lambda=67e-9
;                                ;
; Loop over all particle diameters
;
for i=0,npart_radii-1 do begin
    diameter=2*startparam.ap0[i]
    Stokes_Cunningham=1+2*lambda/diameter*$
      (1.257+0.4*exp(-1.1*diameter/(2*lambda)))
    tau_p=material_density*diameter^2/(18.0*runparam.nu)
    ;
    ; Assume that the radii of all cylinders are the same
    ;
    Stokes=tau_p*init_uu/radius[0]
    ;                            ;
    ; Check how large the box for the initial particle positions is
    ; compared to the radius of the cylinder.
    ;
    if (xdir) then begin
        fractional_area=-startparam.yp0/radius[0]
    endif else begin
        fractional_area=-startparam.xp0/radius[0]
    end
    ;
    ; Print header
    ;
    if (i eq 0) then begin
        print,'Part. dia.','icyl','Stokes','eta',$
          'eta_front','eta_back','n_colls',$
          'Cs',FORMAT='(A12,A6,5A12)'
    endif
    ;
    ; Find the capture efficiency
    ;
    eta=float(solid_colls[*,i])*fractional_area/radii_arr(i)
    front_eta=float(front_colls[*,i])*fractional_area/radii_arr(i)
    back_eta=float(back_colls[*,i])*fractional_area/radii_arr(i)
    for icyl=0,ncylinders-1 do begin
        print,$
          diameter,$
          icyl+1,$
          Stokes,$
          eta[icyl],$
          front_eta[icyl],$
          back_eta[icyl],$
          solid_colls[icyl,i],$
          Stokes_Cunningham,FORMAT='(E12.3,I6,F12.5,3E12.3,I12,E12.3)'
    endfor
    if (ncylinders gt 1) then begin
        print,$
          diameter,$
          'All',$
          Stokes,$
          total(front_eta),$
          total(back_eta),$
          total(solid_colls[*,i]),$
          Stokes_Cunningham,FORMAT='(E12.3,A6,F12.5,2E12.3,I12,E12.3)'
    endif
    if (savefile) then begin
        filename='./data/capture_eff.sav'
        if (i gt 0) then begin
            filename='./data/capture_eff'+str(i)+'.sav'           
        endif
        save,Stokes,eta,Stokes_Cunningham,front_eta,back_eta,ncylinders,$
          filename=filename
    endif
end                             ; Particle diameter loop
END
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; finds removed positions to be
;       used later for plotting them
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cylinder_find_removed_positions,ncylinders,objpvar,irmv,xpos,ypos,radius,$
                                    rmv_pos=rmv_pos,$
                                    on_cylinder_indices=on_cylinder_indices
;
; Find positions of removed particles 
; (to be used later for plotting them).
;
rmv_pos=objpvar.xx(irmv,*)
for icyl=0,ncylinders-1 do begin
    collision_radius=sqrt($
                           (rmv_pos(*,0)-xpos[icyl])^2+$
                           (rmv_pos(*,1)-ypos[icyl])^2)
    if ( icyl eq 0) then begin
        on_cylinder_indices=where(collision_radius lt radius[icyl]*1.1)
    endif else begin
        on_cylinder_indices=[on_cylinder_indices,where(collision_radius lt $
                                       radius[icyl]*1.1)]
    endelse
end
END
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots theta(t) and deposition
;       distribution between particle radii
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cylinder_plot_deposition,ts,oneradius,npart_radii,objpvar,theta_arr,$
                             savefile,startparam,r_i=r_i,$
                             impact_vel_arr=impact_vel_arr
;
; Find where (in radians) the particles hit 
; the surface of the cylinder as a
; function of time
;
particle_hit=0
tmin2=1e38
tmax2=-1e38
print,'The initial time of the simulation is  t =',min(ts.t)
if (oneradius) then begin
    print,'npart_radii=',npart_radii
    READ, r_i, PROMPT = $
      'Enter # of particle radius of interest [0,npart_radii-1]'
endif else begin
    r_i=-1
endelse
first=1
WINDOW,4,XSIZE=256*2,YSIZE=256*2
!p.multi=[0,1,2]
!x.range=[min(ts.t),objpvar.t]
!y.range=[0,!pi]
for i=0,npart_radii-1 do begin
    if (oneradius and i ne r_i) then continue ; jump to next iteration
    theta_=total(theta_arr[*,i,*,0],1)
    time_=total(theta_arr[*,i,*,1],1)
    here=where(theta_ ne 0)
    if (here[0] ne -1) then begin
        particle_hit=1
        theta=theta_[here]
        timereal=time_[here]
        dims=size(theta)
        ind=indgen(dims[1])
        if (min(timereal) lt tmin2) then begin
            tmin2=min(timereal)
        endif
        if (max(timereal) gt tmax2) then begin
            tmax2=max(timereal)
        endif
        if (first) then begin
            plot,timereal,theta,ps=1,ytit='!4h!6',xtit='time'
            if (savefile) then begin
                ap0=startparam.ap0[0:npart_radii-1]
                save,timereal,theta_arr,ap0,$
                  filename='./data/theta_arr.sav'
            endif
            first=0
        endif else begin
            oplot,timereal,theta,ps=1
        end
    endif
    if (oneradius and i eq r_i) then break
end
first=1
!x.range=[min(startparam.ap0[0:npart_radii-1])*0.5,$
          max(startparam.ap0[0:npart_radii-1])*2.0]
for i=0,npart_radii-1 do begin
    theta_=total(theta_arr[*,i,*,0],1)
    rad_=startparam.ap0[i]
    here=where(theta_ ne 0)
    if (here[0] ne -1) then begin
        particle_hit=1
        theta=theta_[here]
        rad=time_[here]
        rad[*]=rad_
        dims=size(theta)
        ind=indgen(dims[1])
        if (first) then begin
            plot_oi,rad,theta,ps=3,ytit='!4h!6',xtit='particle radius'
            first=0
        endif else begin
            oplot,rad,theta,ps=3
        end
    endif
end


WINDOW,5,XSIZE=256*2,YSIZE=256*2
!p.multi=[0,1,2]
!x.range=[min(startparam.ap0[0:npart_radii-1])*0.5,$
          max(startparam.ap0[0:npart_radii-1])*2.0]
!y.range=[min(impact_vel_arr[*,*,*,1:3]),max(impact_vel_arr[*,*,*,1:3])]
first=1
for i=0,npart_radii-1 do begin
    theta_=total(impact_vel_arr[*,i,*,1],1)
    time_=total(impact_vel_arr[*,i,*,0],1)
    here=where(theta_ ne 0)
    if (here[0] ne -1) then begin
        particle_hit=1
        theta=theta_[here]
        rad_=startparam.ap0[i]
        timereal=time_[here]
        rad=time_[here]
        rad[*]=rad_
        dims=size(theta)
        ind=indgen(dims[1])
        if (min(timereal) lt tmin2) then begin
            tmin2=min(timereal)
        endif
        if (max(timereal) gt tmax2) then begin
            tmax2=max(timereal)
        endif
        if (first) then begin
            plot_oi,rad,theta,ps=1,ytit='!8v!dx!n!6',xtit='particle radius'
            if (savefile) then begin
                ap0=startparam.ap0[0:npart_radii-1]
                save,timereal,impact_vel_arr,ap0,$
                  filename='./data/impact_vel_arr.sav'
            endif
            first=0
        endif else begin
            oplot,rad,theta,ps=1
        end
    endif
    if (oneradius and i eq r_i) then break
end

first=1
for i=0,npart_radii-1 do begin
    theta_=total(impact_vel_arr[*,i,*,2],1)
    time_=total(impact_vel_arr[*,i,*,0],1)
    here=where(theta_ ne 0)
    if (here[0] ne -1) then begin
        particle_hit=1
        theta=theta_[here]
        rad_=startparam.ap0[i]
        timereal=time_[here]
        rad=time_[here]
        rad[*]=rad_
        dims=size(theta)
        ind=indgen(dims[1])
        if (min(timereal) lt tmin2) then begin
            tmin2=min(timereal)
        endif
        if (max(timereal) gt tmax2) then begin
            tmax2=max(timereal)
        endif
        if (first) then begin
            plot_oi,rad,theta,ps=1,ytit='!8v!dy!n!6',xtit='particle radius'
            first=0
        endif else begin
            oplot,rad,theta,ps=1
        end
    endif
    if (oneradius and i eq r_i) then break
end


!p.multi=[0,1,1]
!x.range=''
!y.range=''
if (particle_hit) then begin
    print,'The first particle hit the surface at  t =',tmin2
    print,'The last particle hit the surface at   t =',tmax2
endif else begin
    print,'No particle has hit the surface!'
end
print,'The final time of the simulation is    t =',objpvar.t
END 
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
