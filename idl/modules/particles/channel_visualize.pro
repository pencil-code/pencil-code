;
;Script visualizing particle deposition in a channel
;

pro channel_visualize,irmv,trmv,n_steps,$
                      tmin,tmax,png,w,psym,$
                      obj=obj,pvar=pvar

counter_up = 0 ; counter upper wall
counter_lo = 0 ; counter lower wall
rmvdims=size(irmv)

;
;Sort index and time in ascending order of time
;
trmv_sorted = trmv[sort(trmv)] 
irmv_sorted = irmv[sort(trmv)]

;
;Sort the removed particles by upper and lower wall
;
irmv_sorted_lo = irmv_sorted(where(pvar.xx(irmv_sorted(*),1) lt 0))
irmv_sorted_up = irmv_sorted(where(pvar.xx(irmv_sorted(*),1) gt 0))



for i=0,n_steps-1 do begin
    ;
    ; Check if we want to plot for this time
    ;
    if ((obj.t[i] gt tmin) and (obj.t[i] lt tmax)) then begin
        titlestring='t='+str(obj.t[i])       
          ;  plot,obj.xp(*,i),obj.yp(*,i),psym=psym,symsize=1,$
          ;    title=titlestring,position=[0.1,0.1,0.9,0.30],$
          ;    xtitle='x',ytitle='y',charsize=2.5
            
        ; only plot particles that will hit the wall
        plot,obj.xp(irmv,i),obj.yp(irmv,i),$
          psym=psym,symsize=1,title=titlestring,$
          position=[0.1,0.1,0.9,0.30],xtitle='x',ytitle='y',$
          charsize=2.5
        
        if i eq n_steps -1 then begin
            ; Plot the removed particles as blue dots
            oplot,pvar.xx[irmv,0],pvar.xx[irmv,1],col=45,ps=sym(1)
        endif
    endif

    ;plot frames of wallplot
    plot,indgen(10),indgen(10),$
      psym=psym,symsize=1, yr=[0.0,0.07],position=[0.1,0.5,0.36,0.9], $
      /iso,xtitle='x',ytitle='z',charsize=2.2,title='Upper wall', /noerase
    
    plot,indgen(10),indgen(10),$
          psym=psym,symsize=1, yr=[0.0,0.07],position=[0.42,0.5,0.68,0.9], $
          /iso,xtitle='x',ytitle='z',charsize=2.2,title='Lower wall', /noerase

    ; plot the position on the wall where the particle hit
    if (obj.t[i] gt trmv_sorted(counter_up)) $
      or (obj.t[i] gt trmv_sorted(counter_lo)) then begin
        ;upper wall
        plot,pvar.xx[irmv_sorted_up(0:counter_up),0],$
          pvar.xx[irmv_sorted_up(0:counter_up),2],$
          psym=psym,symsize=1, yr=[0.0,0.07],position=[0.1,0.5,0.36,0.9], $
          /iso,xtitle='x',ytitle='z',charsize=2.2,title='Upper wall', /noerase
        ;lowerwall
        plot,pvar.xx[irmv_sorted_lo(0:counter_lo),0],$
          pvar.xx[irmv_sorted_lo(0:counter_lo),2],$
          psym=psym,symsize=1, yr=[0.0,0.07],position=[0.42,0.5,0.68,0.9], $
          /iso,xtitle='x',ytitle='z',charsize=2.2,title='Lower wall', /noerase
        if(counter_up lt size(irmv_sorted_up,/N_ELEMENTS)-1) then begin
            counter_up += 1
        endif
        if(counter_lo lt size(irmv_sorted_lo,/N_ELEMENTS)-1) then begin
            counter_lo += 1
        endif
    endif
    ;
    ; Do we want to write png files or to show results on screen
    ;
    if png eq 1 then begin
        store_png_frame,i
    endif else begin
        wait,w
    endelse
endfor



plot,pvar.xx[irmv_sorted_up(*),0],$
  pvar.xx[irmv_sorted_up(*),2],$
  col=45,psym=psym,symsize=1, yr=[0.0,0.07],position=[0.1,0.5,0.36,0.9], $
  /iso,xtitle='x',ytitle='z',charsize=2.2,title='Upper wall', /noerase

plot,pvar.xx[irmv_sorted_lo(*),0],$
  pvar.xx[irmv_sorted_lo(*),2],$
  col=45,psym=psym,symsize=1, yr=[0.0,0.07],position=[0.42,0.5,0.68,0.9], $
  /iso,xtitle='x',ytitle='z',charsize=2.2,title='Lower wall', /noerase
    
;
; Make array of accumulated deposited particles
;
    
acc_particles = indgen(size(trmv_sorted,/N_ELEMENTS))+1

;tpluss
nu = 2.0e-5
u_tau = 0.27
tpluss = (trmv_sorted-trmv_sorted(0))*u_tau^2/nu
    
;
; Plot number of deposited particles versus time
; 

plot,tpluss,acc_particles,xr=[0,tpluss(size(trmv_sorted,/N_ELEMENTS)-1)],$
  yr=[0,size(trmv_sorted,/N_ELEMENTS)],xtitle='t!u+',$
  ytitle='Number of particles',position=[0.75,0.5,0.98,0.9],$
  charsize=2.2,/noerase

        
END
