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
pro pc_visualize_particles,png=png,removed=removed, savefile=savefile,xmin=xmin,$
                           xmax=xmax,ymin=ymin,ymax=ymax,tmin=tmin,tmax=tmax,$
                           w=w,trace=trace,velofield=velofield,dots=dots,$
                           finalpng=finalpng,xtrace=xtrace,ytrace=ytrace,$
                           store_vec=store_vec
;
device,decompose=0
loadct,5
;
; Set defaults
;
default,png,0
default,finalpng,0
default,removed,0
default,savefile,1
default,tmin,-1e37
default,tmax,1e37
default,w,0.03
default,trace,0
default,xtrace,0
default,ytrace,0
default,velofield,0
default,dots,0
;
; Read dimensions and namelists
;
pc_read_dim,obj=procdim
pc_read_param, object=param
pc_read_pvar,obj=objpvar,/solid_object,irmv=irmv,theta_arr=theta_arr,$
  savefile=savefile,trmv=trmv
pc_read_pstalk,obj=obj
dims=size(obj.xp)
npar=dims[1]
print,'npar=',npar
;
; Find number of cylinders
;
if (param.coord_system eq 'cylindric') then begin
    ncylinders=1
endif else begin
    ncylinders=param.ncylinders
endelse
;
; Set some auxillary variables
;
radius=fltarr(ncylinders,1)
xpos=fltarr(ncylinders,1)
ypos=fltarr(ncylinders,1)
nx=procdim.nx
ny=procdim.ny
dims=size(obj.xp)
n_parts=dims(1)
n_steps=dims(2)
print,'param.coord_system=',param.coord_system
if (param.coord_system eq 'cylindric') then begin
    radius[0]=param.xyz0[0]
    xpos[0]=0.0
    ypos[0]=0.0
    default,xmax,param.xyz1[0]
    default,xmin,-xmax
    default,ymax, xmax
    default,ymin,-xmax
endif else begin
    radius=param.cylinder_radius[0:ncylinders-1]
    xpos=param.cylinder_xpos[0:ncylinders-1]
    ypos=param.cylinder_ypos[0:ncylinders-1]
    default,xmax,param.xyz1[0]
    default,xmin,param.xyz0[0]
    default,ymax,param.xyz1[1]
    default,ymin,param.xyz0[1]
endelse
print,'xmin=',xmin
print,'xmax=',xmax
print,'ymin=',ymin
print,'ymax=',ymax
;
; Find how many different particle radii we have
;
npart_radii=0
dims=size(param.ap0)
if (dims[0] gt 0) then begin
    ninit=dims[1]
    for i=0,ninit-1 do begin
        if (param.ap0[i] ne 0) then begin
            npart_radii=npart_radii+1
        endif
    end
endif else begin
    npart_radii=1
end
print,'npart_radii=',npart_radii
;
; Check if we are inserting particles continuously
;
dummy=findfile('./data/param2.nml', COUNT=countfile)
if (countfile gt 0) then begin
    pc_read_param, object=param2,/param2
    linsert_particles_continuously=param2.linsert_particles_continuously
endif else begin
    pc_read_param, object=param
    linsert_particles_continuously=param.linsert_particles_continuously
endelse
;
;  Check how many particles have collided with a solid object
;
solid_object=1
print_remove_data=0
if (solid_object) then begin
  pc_read_ts,obj=ts
  dims=size(irmv)
  solid_colls=fltarr(ncylinders,npart_radii)
  solid_colls[*,*]=0
  front_colls=solid_colls
  back_colls=solid_colls
  maxy=param.xyz1[1]
  theta_arr=fltarr(npart_radii,10000,2)
  for icyl=0,ncylinders-1 do begin
      if (param.coord_system eq 'cylindric') then begin
          init_uu=param.ampluu
      endif else begin
          init_uu=param.init_uu
      endelse
      for k=long(0),long(dims[1])-1 do begin        
          if (dims[0]>0) then begin              
             x0=objpvar.xx[irmv[k],0]-xpos[icyl]
             y0=objpvar.xx[irmv[k],1]-ypos[icyl]
             deposition_radius2=x0^2+y0^2
             deposition_radius=sqrt(deposition_radius2)
             if (deposition_radius lt radius[icyl]*1.1) then begin
                ipart_radii=0
                while (objpvar.a[irmv[k]] ne param.ap0(ipart_radii)) do begin
                   ipart_radii=ipart_radii+1
                end
                theta_tmp=acos(y0/deposition_radius)
                theta=3.1415-theta_tmp
                if (print_remove_data) then begin
                   print,'time,k,r,x,y,theta=',$
                         trmv[k],irmv[k],deposition_radius,x0,y0,theta
                endif
                if (total(solid_colls[icyl]) lt 10000) then begin
                   theta_arr[ipart_radii,solid_colls[icyl,ipart_radii],0]=theta
                   theta_arr[ipart_radii,solid_colls[icyl,ipart_radii],1]=trmv[k]
                endif
                solid_colls[icyl,ipart_radii]=solid_colls[icyl,ipart_radii]+1
                if (objpvar.xx[irmv[k],1] gt ypos[icyl]) then begin
                   back_colls[icyl,ipart_radii]=back_colls[icyl,ipart_radii]+1
                endif else begin
                   front_colls[icyl,ipart_radii]=front_colls[icyl,ipart_radii]+1
                endelse
             endif
          endif
       endfor
   endfor
  ; 
  ; Find how many particles have been inserted
  ;
  if (linsert_particles_continuously) then begin
     initial_time=ts.t[0]
     final_time=min([objpvar.t,param2.max_particle_insert_time])
     npar_inserted=(final_time-initial_time)*param2.particles_insert_rate
  endif else begin
     npar_inserted=npar
  endelse
  print,'Total number of inserted particles:',npar_inserted
  lambda=67e-9
  ;
  ; Loop over all particle diameters
  ;
  for i=0,npart_radii-1 do begin
     diameter=2*param.ap0[i]
     Stokes_Cunningham=1+2*lambda/diameter*$
                       (1.257+0.4*exp(-1.1*diameter/(2*lambda)))
     tau_p=param.rhops*diameter^2/(18.0*param2.nu)
     ;
     ; Assume that the radii of all cylinders are the same
     ;
     Stokes=tau_p*init_uu/radius[0]
     ;
     ; Check how large the box for the initial particle positions is
     ; compared to the radius of the cylinder.
     ;
     fractional_area=-param.xp0/radius[0]
     ;
     ; Print header
     ;
     if (i eq 0) then begin
        print,'Part. dia.','icyl','Stokes','eta_front','eta_back','n_colls',$
        'Cs',FORMAT='(A12,A6,5A12)'
     endif
     ;
     ; Find the capture efficiency
     ;
     eta=float(solid_colls[*,i])*fractional_area/npar_inserted
     front_eta=float(front_colls[*,i])*fractional_area/npar_inserted
     back_eta=float(back_colls[*,i])*fractional_area/npar_inserted
     for icyl=0,ncylinders-1 do begin
        print,$
           diameter,$
           icyl+1,$
           Stokes,$
           front_eta[icyl],$
           back_eta[icyl],$
           solid_colls[icyl,i],$
           Stokes_Cunningham,FORMAT='(E12.3,I6,F12.5,2E12.3,I12,E12.3)'
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
  end ; Particle diameter loop
endif
;
; Find positions of removed particles (to be used later for plotting them).
;
dims=size(irmv)
if ((removed eq 1) and (dims[0] ne 0) ) then begin
    removed_pos=objpvar.xx(irmv,*)
    skin=nx
    for icyl=0,ncylinders-1 do begin
        collision_radius=sqrt($
                               (removed_pos(*,0)-xpos[icyl])^2+$
                               (removed_pos(*,1)-ypos[icyl])^2)
        if ( icyl eq 0) then begin
            solid_colls=where(collision_radius lt radius[icyl]*1.1)
        endif else begin
            solid_colls=[solid_colls,where(collision_radius lt radius[icyl]*1.1)]
        endelse
    end
endif
;
; Find where (in radians) the particles hit the surface of the cylinder as a
; function of time
;
print,'The initial time of the simulation is  t =',min(obj.t)
for i=0,npart_radii-1 do begin
   theta_=theta_arr[i,*,0]
   time_=theta_arr[i,*,1]
   here=where(theta_ ne 0)
   if (here[0] ne -1) then begin
      WINDOW,4,XSIZE=128*2,YSIZE=256*2
      theta=theta_[here]
      timereal=time_[here]
      dims=size(theta)
      ind=indgen(dims[1])
      !x.range=[0,max(ind)]
      !x.range=[min(obj.t),objpvar.t]
      !y.range=[min(theta),max(theta)]
      if (i eq 0) then begin
         plot,timereal,theta,ps=i+1,ytit='!4h!6',xtit='time'
         print,'The first particle hit the surface at  t =',min(timereal)
         print,'The last particle hit the surface at   t =',max(timereal)
         if (savefile) then begin
            save,timereal,theta,filename='./data/theta.sav'
         endif
      endif else begin
         oplot,timereal,theta,ps=i+1
      end
   endif else begin
      print,'No particles has hit the cylinder surface!'
   endelse
end
print,'The final time of the simulation is    t =',objpvar.t
;
; Set window size
;
xr=xmax-xmin
yr=ymax-ymin
WINDOW,3,XSIZE=1024*xr/yr*1.6,YSIZE=1024
!x.range=[xmin,xmax]
!y.range=[ymin,ymax]
;
; Choose symbol for representing particles
;
if (dots) then begin
    psym=3
endif else begin
    psym=sym(1)
endelse
;
; Show results
;
for i=0,n_steps-1 do begin
   ;
   ; Check if we want to plot for this time
   ;
   if ((obj.t[i] gt tmin) and (obj.t[i] lt tmax)) then begin
      titlestring='t='+str(obj.t[i])
      if (param.coord_system eq 'cylindric') then begin
         xp=obj.xp(*,i)*cos(obj.yp(*,i))
         yp=obj.xp(*,i)*sin(obj.yp(*,i))
         plot,xp,yp,psym=psym,symsize=1,/iso,title=titlestring 
         POLYFILL, CIRCLE_(xpos, $
                           ypos, $
                           radius),color=122        
      endif else begin
         plot,obj.xp(*,i),obj.yp(*,i),psym=psym,symsize=1,/iso,$
           title=titlestring
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
         istr2 = strtrim(string(i,'(I20.4)'),2) ;(only up to 9999 frames)
         file='img_'+istr2+'.png'
         print,'Writing file: ',file
         write_png,file,tvrd(/true)
      endif else begin
         wait,w
      endelse
   endif
end
;
; Do we want to give x and y values of the trace as output?
;
xytrace=0
if (arg_present(xtrace) or arg_present(ytrace)) then begin
    dims=size(obj.xp)
    timeiter=dims[2]
    xtrace=dblarr(npar,timeiter)
    ytrace=dblarr(npar,timeiter)
    xytrace=1
endif
;
; Do we want to show the trace of the particles?
;
ddt=2e-4
if (trace) then begin
   for ipar=0,npar-1 do begin
       if (param.coord_system eq 'cylindric') then begin
           xx0=obj.xp(ipar,*)*cos(obj.yp(ipar,*))
           yy0=obj.xp(ipar,*)*sin(obj.yp(ipar,*))
           ux_loc=cos(obj.yp(ipar,*))*obj.ux(ipar,*)$
                 -sin(obj.yp(ipar,*))*obj.uy(ipar,*)
           uy_loc=sin(obj.yp(ipar,*))*obj.ux(ipar,*)$
                 +cos(obj.yp(ipar,*))*obj.uy(ipar,*)
           ARROW, xx0, yy0,$
             xx0+ux_loc*ddt, $
             yy0+uy_loc*ddt, $
             /data,col=255,HSIZE=4           
       endif else begin
           xx0=obj.xp(ipar,*)
           yy0=obj.yp(ipar,*)
           ARROW, xx0, yy0,$
             xx0+obj.ux(ipar,*)*ddt, $
             yy0+obj.uy(ipar,*)*ddt, $
             /data,col=255,HSIZE=4
       endelse
       oplot,xx0,yy0,ps=3
       if (xytrace) then begin
           xtrace(ipar,*)=xx0
           ytrace(ipar,*)=yy0
       endif
   end 
end
;
; Do we want to overplot the velocity field?
;
if (velofield) then begin
    store_vec=fltarr(10000,4)
    store_count=0
   pc_read_var,obj=objvar
   pc_read_dim,obj=objdim
   l1=objdim.l1
   l2=objdim.l2
   m1=objdim.m1
   m2=objdim.m2
   n1=objdim.n1
   for iii=l1,l2 do begin
       for jjj=m1,m2 do begin
           if (param.coord_system eq 'cylindric') then begin
               xx0=objvar.x(iii)*cos(objvar.y(jjj))
               yy0=objvar.x(iii)*sin(objvar.y(jjj))
               ux=cos(objvar.y(jjj))*objvar.uu(iii,jjj,n1,0)$
                 -sin(objvar.y(jjj))*objvar.uu(iii,jjj,n1,1)
               uy=sin(objvar.y(jjj))*objvar.uu(iii,jjj,n1,0)$
                 +cos(objvar.y(jjj))*objvar.uu(iii,jjj,n1,1)
           endif else begin
               xx0=objvar.x(iii)
               yy0=objvar.y(jjj)
               ux=objvar.uu(iii,jjj,n1,0)
               uy=objvar.uu(iii,jjj,n1,1)
           endelse
           if ((xx0 gt xmin) and (xx0 lt xmax)) then begin
               if ((yy0 gt ymin) and (yy0 lt ymax)) then begin
                   ARROW, xx0, yy0,$
                     xx0+ux*ddt, $
                     yy0+uy*ddt, $
                     /data,col=122,HSIZE=4  
                   store_vec(store_count,0)=xx0
                   store_vec(store_count,1)=yy0
                   store_vec(store_count,2)=xx0+ux*ddt
                   store_vec(store_count,3)=yy0+uy*ddt
                   store_count=store_count+1
               endif
           end
       end
   end
endif
;
; Plot the removed particles as blue dots
;
if (removed eq 1) then begin
    oplot,removed_pos(solid_colls,0),removed_pos(solid_colls,1),col=45,ps=sym(1)
endif
;
; Write png files if required
;
if (png or finalpng) then begin
    istr2 = strtrim(string(i+1,'(I20.4)'),2) ;(only up to 9999 frames)
    file='img_'+istr2+'.png'
    print,'Writing file: ',file
    write_png,file,tvrd(/true)
endif
;
END
