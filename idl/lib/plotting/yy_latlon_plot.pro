   pro yy_latlon_plot, v, quan, ir, colat=colat, lat=lat, lon=lon, interface=interface, oplot=oplot, _extra=extra
;
;  For Yin-Yang grids, plots quantity quan (string of the form <name>_[xyz] for vectors or <name> for scalars, 
;  with <name> a structure element of the output of pc_read_var, v)
;  along a parallel (colat/lat present) or a meridian (lon present) on sphere with index ir. If both lat and lon are defined, lon is ignored.
;  'rho' can be used for quan even if v contains merely lnrho.
;   _extra parameters are handed over to plot.
;
;  title and xtitle of the plot are generated internally unless they are provided as keyword parameters in the call
;
;  6-dec-17/MR: coded
;  
      pc_read_dim, obj=d, /quiet
      if ir lt 0 or ir gt d.nx-1 then begin
        print, 'Invalid ir: must be in [0,'+strtrim(string(d.nx-1),2)+']!'
        return
      endif

      if n_elements(v.y) eq d.my then begin     ; data not trimmed
        m1=d.m1 & m2=d.m2 & n1=d.n1 & n2=d.n2
        ir_full=ir+d.l1
      endif else begin				; data trimmed
        m1=0 & m2=d.ny-1 & n1=0 & n2=d.nz-1
        ir_full=ir
      endelse

      dy=v.y[m1+1]-v.y[m1]
      dz=v.z[n1+1]-v.z[n1]

      comp=0 & pos=strpos(quan,'_') & len=strlen(quan)
      if pos ge 0 and pos lt len-1 then begin
        case (strmid(quan,pos+1)) of
        'x': comp=',0'
        'y': comp=',1'
        'z': comp=',2'
        endcase
        qcomp=strmid(quan,pos+1)
      endif else begin
        pos=len
        comp='' & qcomp=''
      endelse

      quan_=strtrim(strmid(quan,0,pos),2)
      quanv='v.'+quan_ & cbrack=''
      if quan_ eq 'rho' then $
        if where( tag_names(v) eq 'RHO' ) eq -1 then begin
          quanv='exp(v.lnrho' & cbrack=')'
        endif
      
      if strmid(quan_,0,1) eq strmid(quan_,1,1) then quan_=strmid(quan_,1)

      lcolat=is_defined(colat) or is_defined(lat) 
      if lcolat then begin
        if is_defined(colat) then begin
          if (colat lt 0) or (colat ge 180) then begin
            print, 'Error: colatitude < 0 or >= 180!!!'
            return
          endif
          colat *=!dtor
        endif else begin
          if is_defined(lat) then begin
            if (lat lt -90) or (lat ge 90) then begin
              print, 'Error: colatitude < -90 or >= +90!!!'
              return
            endif
            colat=!pi/2-lat*!dtor
          endif
        endelse

        zintu=reverse(v.z[n1]-dz*(indgen(d.nz/6-1)+1))
        zinto=v.z[n2]+dz*(indgen(d.nz/6-1)+1)
        xplt=[zintu,v.z[n1:n2],zinto]

        icolat=(where(v.y[m1:m2] gt colat))[0]

        if icolat le 0 then begin
          cmd1='yplt=reform(yindat)'
          yout=xplt
        endif else begin
          if abs(v.y[m1+1+icolat]-colat) > abs(v.y[m1+icolat]-colat) then icolat-=1
          cmd1='yplt=[(reform(yindat))[0:d.nz/6-2],reform('+quanv+'[ir_full,icolat,n1:n2'+comp+',0]'+cbrack+'),(reform(yindat))[d.nz/6-1:*]]'
          yout=[zintu,zinto]
        endelse
        cmd='yindat=griddata(reform(v.yz(0,*)),reform(v.yz(1,*)),reform('+quanv+'_merge[ir,*'+comp+'])'+cbrack+',xout=[colat],yout=yout,triangles=v.triangles,/linear,/grid)'

      endif else if is_defined(lon) then begin

        if lon lt 0 then lon+=2*!pi $
        else if lon gt 2*!pi then lon = lon mod 2*!pi

        yintu=reverse(v.y[m1]-dy*(indgen(d.ny/2-2)+1))
        yinto=v.y[m2]+dy*(indgen(d.ny/2-2)+1)
        xplt=[yintu,v.y[m1:m2],yinto]

        ilon=(where(v.z[n1:n2] gt lon))[0]

        if ilon le 0 then begin
          cmd1='yplt=reform(yindat)'
          xout=xplt
        endif else begin
          if abs(v.z[n1+1+ilon]-lon) > abs(v.z[n1+ilon]-lon) then ilon-=1
          cmd1='yplt=[(reform(yindat))[0:d.ny/2-2],reform('+quanv+'[ir_full,m1:m2,ilon'+comp+',0])' $
               +cbrack+',(reform(yindat))[d.ny/2-1:*]]'
          xout=[yintu,yinto]
        endelse
        cmd='yindat=griddata(reform(v.yz(0,*)),reform(v.yz(1,*)),reform('+quanv+'_merge[ir,*'+comp+'])' $
            +cbrack+',yout=[lon],xout=xout,triangles=v.triangles,/linear,/grid)'

      endif else begin
        print, 'No interpolation specified.'
        return
      endelse

      ok=execute(cmd)
      if not ok then begin
        print, 'Error when executing command "'+strtrim(cmd,2)+'"'
        stop
        return
      endif
      ok1=execute(cmd1)
      if not ok1 then begin
        print, 'Error when executing command "'+strtrim(cmd1,2)+'"'
        stop
        return
      endif

      if ok and ok1 then begin

        if keyword_set(oplot) then $
          oplot, xplt/!dtor, yplt, _extra=extra $
        else begin

          title=quan_+'!D'+qcomp+'!N at'+ $
                (lcolat ? (is_defined(lat) ? ' latitude='+strtrim(string(lat),2) : ' colatitude=' $
                +strtrim(string(colat/!dtor),2)) : ' longitude='+strtrim(string(lon/!dtor),2)) $
                +'!Uo!N'
          xtitle=(lcolat ? 'longitude' : 'colatitude') + ' [!Uo!N]'
          yrange=[min(yplt),max(yplt)]

          if n_elements(extra) gt 0 then begin
            if has_tag(extra,'title') then $
              title=extra.title
              
            if has_tag(extra,'xtitle') then $
              xtitle=extra.xtitle 
            
            if has_tag(extra,'yrange') then $
              yrange=extra.yrange 
          endif

          plot, xplt/!dtor, yplt, title=title, xtitle=xtitle, yrange=yrange, _extra=extra

        endelse

        if keyword_set(interface) then begin
          plots, [1./4,1./4]*!pi/!dtor, !y.crange, linest=1
          if lcolat then $
            plots, [7./4,7./4]*!pi/!dtor, !y.crange, linest=1 $
          else $
            plots, [3./4,3./4]*!pi, !y.crange, linest=1
        endif

      endif

   end
