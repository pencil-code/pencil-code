   pro yy_latlon_plot, v, quan, ir, lat=lat, lon=lon, _extra=extra
;
;  For Yin-Yang grids, plots quantity quan (string of the form <name>_[xyz] for vectors or <name> for scalars, 
;  with <name> a structure element of the output of pc_read_var, vv)
;  along a parallel (lat present) or a meridian (lon present) on sphere with index ir. _extra parameters are handed over to plot.
;
;  6-dec-17/MR: coded
;  
      pc_read_dim, obj=d

      if n_elements(v.y) eq d.my then begin
        m1=d.m1 & m2=d.m2 & n1=d.n1 & n2=d.n2
      endif else begin
        m1=0 & m2=d.ny-1 & n1=0 & n2=d.nz-1
      endelse

      dy=v.y[m1+1]-v.y[m1]
      dz=v.z[n1+1]-v.z[n1]

      comp=0 & pos=strpos(quan,'_') & len=strlen(quan)
      if pos ge 0 and pos lt len-1 then begin
        case (strmid(quan,pos+1)) of
        'x': comp=0
        'y': comp=1
        'z': comp=2
        endcase
        qcomp=strmid(quan,pos+1)
      endif else begin
        pos=len
        qcomp=''
      endelse

      quan_=strtrim(strmid(quan,0,pos),2)
      quanv='v.'+quan_ & cbrack=''
      if quan_ eq 'rho' then $
        if where( tag_names(v) eq 'RHO' ) eq -1 then begin
          quanv='exp(v.lnrho' & cbrack=')'
        endif
      
      if strmid(quan_,0,1) eq strmid(quan_,1,1) then quan_=strmid(quan_,1)
 
      if is_defined(lat) then begin

        if (lat le 0) or (lat ge !pi) then begin
          print, 'lat <= 0 or lat >= Pi'
          return
        endif

        zintu=reverse(v.z[n1]-dz*(indgen(d.nz/6-1)+1))
        zinto=v.z[n2]+dz*(indgen(d.nz/6-1)+1)
        xplt=[zintu,v.z[n1:n2],zinto]

        ilat=(where(v.y[m1:m2] gt lat))[0]

        if ilat le 0 then begin
          cmd1='yplt=reform(yindat)'
          yout=xplt
        endif else begin
          if abs(v.y[m1+1+ilat]-lat) > abs(v.y[m1+ilat]-lat) then ilat-=1
          cmd1='yplt=[(reform(yindat))[0:d.nz/6-2],reform('+quanv+'[ir,ilat,n1:n2,0]'+cbrack+'),(reform(yindat))[d.nz/6-1:*]]'
          yout=[zintu,zinto]
        endelse
        cmd='yindat=griddata(reform(v.yz(0,*)),reform(v.yz(1,*)),reform('+quanv+'_merge[ir,*,comp])'+cbrack+',xout=[lat],yout=yout,triangles=v.triangles,/linear,/grid)'

      endif else if is_defined(lon) then begin

        if lon lt 0 then lon+=2*!pi $
        else if lon gt 2*!pi then lon = lon mod 2*!pi

        yintu=reverse(v.y[m1]-dy*(indgen(d.ny/2-1)+1))
        yinto=v.y[m2]+dy*(indgen(d.ny/2-1)+1)
        xplt=[yintu,v.y[m1:m2],yinto]

        ilon=(where(v.z[n1:n2] gt lon))[0]

        if ilon le 0 then begin
          cmd1='yplt=reform(yindat)'
          xout=xplt
        endif else begin
          if abs(v.z[n1+1+ilon]-lon) > abs(v.z[n1+ilon]-lon) then ilon-=1
          cmd1='yplt=[(reform(yindat))[0:d.ny/3-2],reform('+quanv+'[ir,m1:m2,ilon,0])'+cbrack+',(reform(yindat))[d.ny/3-1:*]]'
          xout=[yintu,yinto]
        endelse
        cmd='yindat=griddata(reform(v.yz(0,*)),reform(v.yz(1,*)),reform('+quanv+'_merge[ir,*,comp])'+cbrack+',yout=[lon],xout=xout,triangles=v.triangles,/linear,/grid)'

      endif else begin
        print, 'No interpolation specified.'
        return
      endelse

      ok=execute(cmd)
      if not ok then begin
        print, 'Error when executing command "'+strtrim(cmd,2)+'"'
        stop
      endif
      ok1=execute(cmd1)
      if not ok1 then begin
        print, 'Error when executing command "'+strtrim(cmd1,2)+'"'
        stop
      endif

      if ok and ok1 then $
        plot, xplt, yplt, title=quan_+'!D'+qcomp+'!N at'+ $
             (is_defined(lat) ? ' latitude='+strtrim(string(lat),2) : ' longitude='+strtrim(string(lon),2)),xtitle=(is_defined(lat) ? 'longitude' : 'latitude'), $
              _extra=extra

   end
