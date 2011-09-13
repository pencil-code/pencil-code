
;
; This script will visualise the flow of particles around one or more
; solid objects or inside a channel, and will find capture/deposition 
; efficiencies
;
;  Author: Nils Erland L. Haugen, Anders G. Bjørnstad and Solveig S. Alnes
;

;
pro pc_visualize_particles,png=png,removed=removed,savefile=savefile,$
                                xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,$
                                tmin=tmin,tmax=tmax,w=w,trace=trace,$
                                velofield=velofield,dots=dots,finalpng=finalpng,$
                                xtrace=xtrace,ytrace=ytrace,store_vec=store_vec,$
                                theta_arr=theta_arr,npart=npart,$
                                noviz=noviz, oneradius=oneradius,xdir=xdir,$
                                nopvar=nopvar,channel=channel,plotvort=plotvort,$
                           solid_object=solid_object,plane=plane
;
; compile needed procedures
resolve_routine,'solid_object_findcapture',/compile_full_file
resolve_routine,'solid_object_visualize',/compile_full_file

device,decompose=0
loadct,5
!p.multi=[0,1,1]
;
; Set defaults
;
default,plane,'xy'
default,png,0
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
default,finalpng,0
default,xdir,0       ; flow in x-direction
default,noviz,0  ; don't visualize stalked particles
default,oneradius,0  ; show only one radius if more than one are present
default,nopvar,0 ; don't read pvar to save time
default,noread,0
default,channel,0 ; solid walls
default,solid_object,1
default,plotvort,0; plot vorticity plot
print,'noread=',noread
;
; Read dimensions and namelists
;
pc_read_dim,obj=procdim
pc_read_param, object=param
if (nopvar eq 0) then begin
    pc_read_pvar,obj=objpvar,solid_object=solid_object,irmv=irmv,$
      savefile=savefile,trmv=trmv
endif
pc_read_pstalk,obj=obj

;
; Find number of solid objects
;

if (param.coord_system eq 'cylindric') then begin
    ncylinders=1
endif else begin
    if(channel) then begin 
        solid_object=0
    endif else begin
        if (solid_object) then begin
            ncylinders=param.ncylinders
        endif
    endelse
endelse

if (solid_object) then begin
    ; Find position and radius of cylinder(s)
    radius=fltarr(ncylinders,1)
    xpos=fltarr(ncylinders,1)
    ypos=fltarr(ncylinders,1)
    if (param.coord_system eq 'cylindric') then begin
        radius[0]=param.xyz0[0]
        xpos[0]=0.0
        ypos[0]=0.0
    endif else begin
        radius=param.cylinder_radius[0:ncylinders-1]
        xpos=param.cylinder_xpos[0:ncylinders-1]
        ypos=param.cylinder_ypos[0:ncylinders-1]
    endelse
endif 
;
; Set some auxillary variables
;
nx=procdim.nx
ny=procdim.ny
dims=size(obj.xp)
n_parts=dims(1)
n_steps=dims(2)
print,'param.coord_system=',param.coord_system
if (param.coord_system eq 'cylindric') then begin
    default,xmax,param.xyz1[0]
    default,xmin,-xmax
    default,ymax, xmax
    default,ymin,-xmax
endif else begin
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
if (file_test('./data/param2.nml')) then begin
    pc_read_param, object=param2,/param2
    linsert_particles_continuously=param2.linsert_particles_continuously
endif else begin
    pc_read_param, object=param
    linsert_particles_continuously=param.linsert_particles_continuously
endelse
;
; Check if we hare using the old version, 'rhops', or the new one; 'rhopmat'.
;
res=tag_names(param)
dims=size(res)
ndims=dims[1]
found=0
target='RHOPS'
for i=0,ndims-1 do begin
    print,res[i],' ',target
    if res[i] eq target then begin
        found=1
    endif
end
print,'found=',found
if (found) then begin
    material_density=param.rhops
endif else begin
    material_density=param.rhopmat
end
;
; Some things should be done only if we have read pvar.dat
;
if (not nopvar) then begin
    ;
    ; Find total number of particles with the different radii
    ;
    npart=fltarr(npart_radii)
    for ipart_radii=0,npart_radii-1 do begin
        npart[ipart_radii]=$
          n_elements(where(float(objpvar.a) eq float(param.ap0[ipart_radii])))
    end
    print,'npart=    ',npart
    print,'param.ap0=',param.ap0
    ;
    ;  Check how many particles have collided with a solid object
    ;
     if (solid_object) then begin
         solid_object_findcapture,ncylinders=ncylinders,npart_radii=npart_radii,$
           startparam=param,runparam=param2,$
           linsert_particles_continuously=linsert_particles_continuously,$
           objpvar=objpvar,irmv=irmv,trmv=trmv,xdir=xdir,radii_arr=npart,$
           savefile=savefile,removed=removed,oneradius=oneradius,radius=radius,$
           xpos=xpos,ypos=ypos,rmv_pos=removed_pos,$
           on_cylinder_indices=on_cylinder_indices,r_i=r_i,$
           material_density=material_density
    endif ; solid_object
    ;
    ;  Check how many particles have collided with a the walls
    ;
    if (channel) then begin
        channel_deposition,npart_radii,irmv,trmv,$
          objpvar=objpvar,startparam=param,savefile=savefile
    endif ; channel

endif ;nopvar=0 endif



if (noviz eq 0) then begin
    if (nopvar and oneradius) then begin
        print,'npart_radii=',npart_radii
        READ, r_i, PROMPT = $
          'Enter # of particle radius of interest [0,npart_radii-1]'
    endif
    if (oneradius) then begin
        if(solid_object) then begin
            diameter=2*param.ap0
            tau_p=material_density*diameter^2/(18.0*param2.nu)
            tau_f=radius[0]/param.init_uu
            Stokes=tau_p/tau_f
            print,'Only visualizing particles with Stokes no St=',Stokes[r_i]
        endif else if (channel) then begin
            print,'Only visualizing particles with radius=',param.ap0(r_i)
        endif
    endif
    ;
    ; Set window size
    ;
    xr=xmax-xmin
    yr=ymax-ymin
    WINDOW,3;,XSIZE=1024*xr/yr*1.6,YSIZE=1024
    !x.range=[xmin,xmax]
    !y.range=[ymin,ymax]
    !x.style=1
    !y.style=1
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
;   if (solid_object) then begin
;    endif
   if (channel) then begin
       channel_visualize,irmv,trmv,n_steps,tmin,tmax,png,w,psym,$
         obj=obj,pvar=objpvar
       plot_vorticity,npart_radii=npart_radii,irmv=irmv,trmv=trmv,$
         totpar=n_parts,nts=n_steps,savefile=savefile,$
         prmv=prvm,$
         startparam=startparam,obj=obj
   endif else begin
       if (plane eq 'xy') then begin
           solid_object_visualize,ncylinders,n_steps,obj,tmin,tmax,param,$
             xpos,ypos,radius,oneradius,png,w,psym,obj.xp,obj.yp,$
             Stokes=Stokes,r_i=r_i,i=i,solid_object=solid_object
       endif else if (plane eq 'xz') then begin
           solid_object_visualize,ncylinders,n_steps,obj,tmin,tmax,param,$
             xpos,ypos,radius,oneradius,png,w,psym,obj.xp,obj.zp,$
             Stokes=Stokes,r_i=r_i,i=i,solid_object=solid_object
       endif else if (plane eq 'yz') then begin
           solid_object_visualize,ncylinders,n_steps,obj,tmin,tmax,param,$
             xpos,ypos,radius,oneradius,png,w,psym,obj.yp,obj.zp,$
             Stokes=Stokes,r_i=r_i,i=i,solid_object=solid_object
       end
   end

  ;
  ; Do we want to give x and y values of the trace as output?
  ;
  xytrace=0
  if (arg_present(xtrace) or arg_present(ytrace)) then begin
      dims=size(obj.xp)
      timeiter=dims[2]
      xtrace=dblarr(npart,timeiter)
      ytrace=dblarr(npart,timeiter)
      xytrace=1
  endif
  
  ;
  ; Do we want to show the trace of the particles?
  ;
  ddt=2e-4
  if (trace) then begin
      for ipar=0,npart-1 do begin
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
          endif else if (tracermv) then begin
              rmvdims=size(irmv)
              for k=0,rmvdims[1]-1 do begin        
                  xx0=obj.xp(irmv(k),*)
                  yy0=obj.yp(irmv(k),*)
          
                  ARROW, xx0, yy0,$
                    xx0+obj.ux(irmv(k),*)*ddt, $
                    yy0+obj.uy(irmv(k),*)*ddt, $
                    /data,col=255,HSIZE=4
              endfor
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
      endfor
  endif

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
              endif 
          endfor  
      endfor 
  endif 
;
; Plot the removed particles as blue dots
;
if (removed eq 1) then begin
    oplot,removed_pos(on_cylinder_indices,0),removed_pos(on_cylinder_indices,1),col=45,ps=sym(1)
endif
;
; Write png files if required
;
  if (png or finalpng) then begin
      store_png_frame,i+1
  endif
endif ; end noviz

!x.range=''
!y.range=''
;
END

pro store_png_frame,i
 istr2 = strtrim(string(i,'(I20.4)'),2) ;(only up to 9999 frames)
 file='img_'+istr2+'.png'
 print,'Writing file: ',file
 write_png,file,tvrd(/true)
END
