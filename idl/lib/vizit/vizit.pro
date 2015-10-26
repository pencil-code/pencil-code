; $Id$

pro safe_obj_destroy,obj
  obj_size = n_elements(obj)
  if obj_size eq 1 then begin
    if obj_valid(obj) then obj_destroy,obj
  endif else begin
    for i = 0L, obj_size-1 do safe_obj_destroy,obj[i]
  endelse
end

PRO vizit_isoupd, data, state
    wBase = state.wBase
    widget_control, wBase, bad=bad, get_uvalue=state
    if bad ne 0 then return
    if nlevels eq 1 then begin
      shade_volume, data, state.level, vert ,poly
      *(state.data) = data
      ;if n_elements(vert) le 1 then return
      if n_elements(vert) le 1 then begin
        vert=fltarr(3,3)
        vert[*,*]=0.
        poly=[3,0,1,2]
      endif
      if state.nonequidistant eq 1 then begin
        vert[0,*]=interpol(state.xx,findgen(n_elements(state.xx)),reform(vert[0,*]))
        vert[1,*]=interpol(state.yy,findgen(n_elements(state.yy)),reform(vert[1,*]))
        vert[2,*]=interpol(state.zz,findgen(n_elements(state.zz)),reform(vert[2,*]))
      endif
      state.oPolygon->Setproperty, data=vert, poly=poly
      oModel->add, oPolygon
    endif else begin
      ;for eachlevel=0,nlevels-1 do oModel->add, oPolygon[eachlevel]
    endelse
    ;state.oModel->Rotate, [0,1,0], 5
    state.oMainWindow->draw, state.oMainView

    min = min(data, max=max)
    hist = alog10(histogram(data, min=min, max=max, binsize=(max-min)/255.)+0.1)
    xcc = [0.,1./n_elements(hist)] & ycc = [0.,1./max(hist)]
    state.min = min
    state.max = max
    state.oPlot->SetProperty, datay=hist, XCOORD=xcc, YCOORD=ycc
    index = byte(255.*(state.level-min)/(max-min))
    state.oLine->SetProperty, datay=[0,max(hist)], datax=[index,index], $
      color=255, XCOORD=xcc, YCOORD=ycc
    state.oPlotWindow->draw, state.oPlotView
    widget_control, wBase, set_uvalue=state
END

;----------------------------------------------------------------
PRO vizit, object, xsize=xsz, ysize=ysz, colortable=coltab, $
       dx=dx, dy=dy, dz=dz, state=state, bgcolor=bgcolor, $
       Lx=Lx,Ly=Ly,Lz=Lz,xyz0=xyz0, $
       plevel=plevel, level=level, index=index, $
       pVector=pVector, alpha=alpha, $
       bboxcolor=bboxcolor, bboxfgcolor=bboxfgcolor, bboxbgcolor=bboxbgcolor, $
       xtitle=xtitle,ytitle=ytitle,ztitle=ztitle,title=title, $
       xrange=xrange,yrange=yrange,zrange=zrange, $
       surfaces=surfaces, surfacealpha=surfacealpha, $
       surfaceindex=surfaceindex, noiso=noiso, $
       png=png,gif=gif,jpg=jpg,filename=filename, $
       xx=xx,yy=yy,zz=zz, nonequidistant=nonequidistant, $
       calcnormals=calcnormals, $
       hideXAxis, hideYaxis, hideZaxis, $
       drawbbox, drawcontent, $
       drawdatabox, $
       tmat=tmat, bboxthick=bboxthick, $
       frameno=frameno, format_frameno=format_frameno,$
       DEPTH_REORDER=DEPTH_REORDER, zclip=zclip, $
       clip_planes=clip_planes, gradvar=gradvar, $
       _extra=extra

  COMMON colors, ro,go, bo, r,g,b
  COMMON vizit, curr_coltab, curr_r, curr_g, curr_b


  IF N_PARAMS() EQ 0 THEN BEGIN
    print, "VIZIT, object, xsize=xsz, ysize=ysz, colortable=coltab, dx=dx, dy=dy, dz=dz"
    print, "    object: a 3D data array, or a structure {data:d, dx:dx, dy:dy, dz:dz}"
    print, "     xsize: window x-size"
    print, "     ysize: window y-size"
    print, "colortable: color table number (default=40), or an RGB triplet"
    print, " "
    print, "  left mousebutton: rotate"
    print, "middle mousebutton: translate"
    print, " right mousebutton: dolly"
    return
  ENDIF

  default, xsz, 600
  default, ysz, 600

  default, clip_planes, -1

  default, bgcolor, [50,50,50]
  default, bboxcolor, [255,255,255]
  default, coltab, 39
  default, curr_coltab, -1
  default, bboxthick, 1.

  default, bboxfgcolor, bboxcolor
  default, bboxbgcolor, uint(bgcolor*0.25+bboxfgcolor*0.75)

  default, xtitle, ''
  default, ytitle, ''
  default, ztitle, ''
  default, hideXAxis, 1
  default, hideYAxis, 1
  default, hideZAxis, 1
  default, drawdatabox, 0
  if keyword_set(PNG) then begin
    default, drawbbox, 1
    default, drawcontent, 1
  endif else begin
    default, drawbbox, 0
    default, drawcontent, 0
  endelse

  default, format_frameno,'(I06)'

  if keyword_set(nonequidistant) then nonequidistant=1 else nonequidistant=0
  if keyword_set(calcnormals) then calcnormals=1 else calcnormals=0

  s=size(object)
  IF s[s[0]+1] EQ 8 THEN BEGIN
    default,field,'data'
    if tag_exists(object, field) then begin
     res=execute(strjoin(['data=object.',field]),1)
    endif
    if tag_exists(object, 'vector') then pVector = object.vector
    if tag_exists(object, 'dx'    ) then default, dx, object.dx
    if tag_exists(object, 'dy'    ) then default, dy, object.dy
    if tag_exists(object, 'dz'    ) then default, dz, object.dz
    if tag_exists(object,'xyz0'   ) then default, xyz0, object.xyz0
  END ELSE BEGIN
    data=object
  END

  s = size(data)
  mx = s[1]  & my = s[2]    & mz = s[3]

  if not keyword_set(nonequidistant) then begin
    default, dx, 1.
    default, dy, 1.
    default, dz, 1.
    default, Lx, dx*mx
    default, Ly, dy*my
    default, Lz, dz*mz
    default, xyz0, [-Lx/2.,-Ly/2.,-Lz/2.]

    default,xrange,[xyz0[0],xyz0[0]+Lx]
    default,yrange,[xyz0[1],xyz0[1]+Ly]
    default,zrange,[xyz0[2],xyz0[2]+Lz]
  endif else begin
    default, xx, findgen(mx)
    default, yy, findgen(my)
    default, zz, findgen(mz)
    default, Lx, max(xx)-min(xx)
    default, Ly, max(yy)-min(yy)
    default, Lz, max(zz)-min(zz)
    default, xyz0, [min(xx),min(yy),min(zz)]
  endelse

  default,xrange,[xyz0[0],xyz0[0]+Lx]
  default,yrange,[xyz0[1],xyz0[1]+Ly]
  default,zrange,[xyz0[2],xyz0[2]+Lz]


  size = sqrt((Lx^2+Ly^2+Lz^2)/3.)

  box_size=[Lx,Ly,Lz]/size
  box_xyz0=-(box_size/2.)

  xcc_surfaces = [box_xyz0[0]-box_size[0]*(xyz0[0]/Lx),box_size[0]/Lx]
  ycc_surfaces = [box_xyz0[1]-box_size[1]*(xyz0[1]/Ly),box_size[1]/Ly]
  zcc_surfaces = [box_xyz0[2]-box_size[2]*(xyz0[2]/Lz),box_size[2]/Lz]

  xcc_polygon = [box_xyz0[0],box_size[0]/mx]
  ycc_polygon = [box_xyz0[1],box_size[1]/my]
  zcc_polygon = [box_xyz0[2],box_size[2]/mz]

  xcc_axis = [box_xyz0[0]-box_size[0]*(xyz0[0]/Lx),box_size[0]/Lx]
  ycc_axis = [box_xyz0[1]-box_size[1]*(xyz0[1]/Ly),box_size[1]/Ly]
  zcc_axis = [box_xyz0[2]-box_size[2]*(xyz0[2]/Lz),box_size[2]/Lz]

  default,axes_intersect,[min(xrange),min(yrange),min(zrange)]

  dev = !d.name
  set_plot,'z'

  IF n_elements(coltab) eq 1 THEN BEGIN
      if (coltab ne curr_coltab) then begin
        loadct, coltab, ncol=256
        tvlct, curr_r, curr_g, curr_b, /GET
        curr_coltab=coltab
      endif
  END ELSE IF n_elements(coltab) eq 3 THEN BEGIN
      curr_r = replicate(coltab[0], 256)
      curr_g = replicate(coltab[1], 256)
      curr_b = replicate(coltab[2], 256)
      curr_coltab=-1
  endif else begin
    if (curr_coltab ne 0) then begin
      loadct, 0, ncol=256
      tvlct, curr_r, curr_g, curr_b, /GET
      curr_coltab=0
      print, 'COLORTABLE should be either a colortable number or an RGB triple'
    endif
  endelse

  set_plot, dev

  oPalette = obj_new('IDLgrPalette', RED=curr_r, GREEN=curr_g, BLUE=curr_b)

  min = min(data,max=max)
  if min eq max then begin              ; create fake variation
    eps = 1e-6*abs(max)+1e-20
    data(mx/2,my/2,mz/2) = data(mx/2,my/2,mz/2)+eps
    min = min-eps
    default,level,min+0.5*eps
    level[*]=min+0.5*eps
    plevel=level
    plevel[*]=0.5
    nlevels=n_elements(level)
    default, index, 128
  endif else begin
    default, plevel, 0.9
    default, level, (1.-plevel)*min + plevel*max
    plevel=(level-min)/(max-min)
    nlevels=n_elements(level)
    default, index, byte(255.*plevel)
  endelse

;
; Generate extra surfaces
;
  nsurfaces=n_elements(surfaces)
  oSurfaces=0
  oTexture=0
  if (nsurfaces gt 0) then begin
    surfnumbers=findgen(nsurfaces)/nsurfaces
    default, surfaceindex, byte((127.*surfnumbers)+128)
    oSurfaces=objarr(nsurfaces)
    default,surfacealpha,make_array(nsurfaces,value=255,/int)
    for eachsurface=0,nsurfaces-1 do begin	
      texture=fltarr(2,2,2)
      texture[0,*,*]=surfaceindex[eachsurface]
      texture[1,*,*]=surfacealpha[eachsurface]
      oTexture=obj_new('IDLgrImage',texture,interleave=0)

      undefine,normals
      if (has_tag(surfaces[eachsurface], 'normals')) then begin
        if (ptr_valid(surfaces[eachsurface].normals)) then begin
          normals=*(surfaces[eachsurface].normals)
        endif
      endif

      oSurfaces[eachsurface] = obj_new('IDLgrPolygon', $
                       *(surfaces[eachsurface].vertices), $
              POLYGONS=*(surfaces[eachsurface].triangles), $
              NORMALS=normals, $
              COLOR=surfaceindex[eachsurface], $
              BOTTOM=surfaceindex[eachsurface], $
              HIDE=1-drawcontent,  $
              PALETTE=oPalette, /SHADING, $
              XCOORD_CONV=xcc_surfaces, $
              YCOORD_CONV=ycc_surfaces, $
              ZCOORD_CONV=zcc_surfaces, $
              TEXTURE_MAP=oTexture)
       if ((((size(tmat))[0] eq 2) and n_elements(tmat) eq 16)) then begin
         if keyword_set(DEPTH_REORDER) then reorder_by_depth,oSurfaces[eachsurface],tmat=tmat
       endif
    endfor
  endif

if keyword_set(noiso) then begin
  oPolygon=0
endif else begin
;
; Render Isosurface Objects (can be regenerated dynamically later)
;
  if (nlevels eq 1) then begin
    shade_volume, data, level, vert, poly
    if n_elements(vert) le 1 then begin
      vert=fltarr(3,3)
      vert[*,*]=0.
      poly=[3,0,1,2]
    endif
    default,alpha,255
    index_shape=n_elements(index)
    if (index_shape eq 1) then begin
      print,"Color from color table"
      texture=fltarr(2,2,2)
      texture[0,*,*]=index
      texture[1,*,*]=alpha
    endif else if (index_shape eq 3) then begin
      print,"Color from RGB values"
      texture=fltarr(4,2,2)
      texture[0,*,*]=255
      texture[1,*,*]=255
      texture[2,*,*]=255
      ;texture[0,*,*]=index[0]
      ;texture[1,*,*]=index[1]
      ;texture[2,*,*]=index[2]
      texture[3,*,*]=alpha
    endif else begin
       print,"Color index must be scalar colour index or [3] R,G,B value"
       exit
    endelse

    oTexture=obj_new('IDLgrImage',texture,interleave=0)

    if calcnormals eq 1 then begin
      default,gradvar,grad(data)
      normals=vert
      normals[0,*]=interpolate(gradvar[*,*,*,0],vert[0,*],vert[1,*],vert[2,*])
      normals[1,*]=interpolate(gradvar[*,*,*,1],vert[0,*],vert[1,*],vert[2,*])
      normals[2,*]=interpolate(gradvar[*,*,*,2],vert[0,*],vert[1,*],vert[2,*])
    endif

    if nonequidistant eq 1 then begin
      vert[0,*]=interpol(xx,findgen(n_elements(xx)),reform(vert[0,*]))
      vert[1,*]=interpol(yy,findgen(n_elements(yy)),reform(vert[1,*]))
      vert[2,*]=interpol(zz,findgen(n_elements(zz)),reform(vert[2,*]))
      oPolygon = obj_new('IDLgrPolygon', vert, POLYGONS=poly, $
        COLOR=index, BOTTOM=index, PALETTE=oPalette, $
        /SHADING, TEXTURE_MAP=oTexture, $
        HIDE=1-drawcontent,  $
        NORMALS=normals,  $
        XCOORD_CONV=xcc_surfaces, $
        YCOORD_CONV=ycc_surfaces, $
        ZCOORD_CONV=zcc_surfaces)
        if (((size(tmat))[0] eq 2) and n_elements(tmat) eq 16) then begin
          if keyword_set(DEPTH_REORDER) then reorder_by_depth,oPolygon,tmat=tmat
        endif
    endif else begin
;      print,index,curr_r[index],curr_g[index],curr_b[index]
      oPolygon = obj_new('IDLgrPolygon', vert, POLYGONS=poly, $
        COLOR=index, BOTTOM=index, PALETTE=oPalette, $
;        COLOR=255, BOTTOM=255, PALETTE=oPalette, $
        /SHADING, TEXTURE_MAP=oTexture, $
        HIDE=1-drawcontent,  $
        NORMALS=normals,  $
        XCOORD_CONV=xcc_polygon, $
        YCOORD_CONV=ycc_polygon, $
        ZCOORD_CONV=zcc_polygon)
        if (((size(tmat))[0] eq 2) and n_elements(tmat) eq 16) then begin
          if keyword_set(DEPTH_REORDER) then reorder_by_depth,oPolygon,tmat=tmat
        endif
    endelse
  endif else begin
    oPolygon=objarr(nlevels)
    for eachlevel=0,nlevels-1 do begin	
      default,alpha,make_array(nlevels,value=255,/int)
      shade_volume, data, level[eachlevel], vert, poly
      if n_elements(vert) le 1 then begin
        vert=fltarr(3,3)
        vert[*,*]=0.
        poly=[3,0,1,2]
      endif
    index_shape=size(index)
    if (index_shape[0] eq 1) then begin
      texture=fltarr(2,2,2)
      texture[0,*,*]=index[eachlevel]
      texture[1,*,*]=alpha[eachlevel]
    endif else if (index_shape[0] eq 2) then begin
      texture=fltarr(4,2,2)
;      print,index[*,eachlevel]
      ;texture[0,*,*]=255
      ;texture[1,*,*]=255
      ;texture[2,*,*]=255
      texture[0,*,*]=index[0,eachlevel]
      texture[1,*,*]=index[1,eachlevel]
      texture[2,*,*]=index[2,eachlevel]
      texture[3,*,*]=alpha
    endif else begin
       print,"Color index must be scalar colour index or [3] R,G,B value"
       exit
    endelse
    ;  texture=fltarr(2,2,2)
    ;  texture[0,*,*]=index[eachlevel]
    ;  texture[1,*,*]=alpha[eachlevel]
      oTexture=obj_new('IDLgrImage',texture,interleave=0)
      if nonequidistant eq 1 then begin
        vert[0,*]=interpol(xx,findgen(n_elements(xx)),vert[0,*])
        vert[1,*]=interpol(yy,findgen(n_elements(yy)),vert[1,*])
        vert[2,*]=interpol(zz,findgen(n_elements(zz)),vert[2,*])
        print,"Attempted to interpolate surface coords to streached grid..."
        oPolygon[eachlevel] = obj_new('IDLgrPolygon', vert, POLYGONS=poly, $
;          COLOR=index[eachlevel], BOTTOM=index[eachlevel], PALETTE=oPalette, $
          COLOR=255, BOTTOM=255, PALETTE=oPalette, $
          /SHADING, TEXTURE_MAP=oTexture, $
          HIDE=1-drawcontent,  $
          XCOORD_CONV=xcc_surfaces, $
          YCOORD_CONV=ycc_surfaces, $
          ZCOORD_CONV=zcc_surfaces)

          if (((size(tmat))[0] eq 2) and n_elements(tmat) eq 16) then begin
            if keyword_set(DEPTH_REORDER) then reorder_by_depth,oPolygon[eachlevel],tmat=tmat
          endif
      endif else begin
        oPolygon[eachlevel] = obj_new('IDLgrPolygon', vert, POLYGONS=poly, $
          COLOR=index[eachlevel], BOTTOM=index[eachlevel], PALETTE=oPalette, $
          /SHADING, TEXTURE_MAP=oTexture, $
          HIDE=1-drawcontent,  $
          XCOORD_CONV=xcc_polygon, $
          YCOORD_CONV=ycc_polygon, $
          ZCOORD_CONV=zcc_polygon)
          if (((size(tmat))[0] eq 2) and n_elements(tmat) eq 16) then begin
            if keyword_set(DEPTH_REORDER) then reorder_by_depth,oPolygon[eachlevel],tmat=tmat
          endif
      endelse
    endfor
  endelse
endelse
  minbboxpal=120
  BBoxPal=intarr(3,256)
;  BBoxPal=spread(bgcolor,1,256)
  for bboxcol=0,255 do begin
    frac=bboxcol/255
    BBoxPal[*,bboxcol]=uint(frac*bboxfgcolor+(1.-frac)*bboxbgcolor)
;    BBoxPal=round(bgcolor+(bboxcolor-bgcolor)*float(bboxcol+1)/(256.-minbboxpal))
  endfor
  pnts=where (BBoxPal gt 255, count)
  if count gt 0 then BBoxPal[pnts]=255
  pnts=where (BBoxPal lt 0, count)
  if count gt 0 then BBoxPal[pnts]=0

  oBBoxPalette = obj_new('IDLgrPalette',BBoxPal[0,*],BBoxPal[1,*],BBoxPal[2,*])
  XBox = box_size[0]*[0,1,0,1,0,1,0,1]+box_xyz0[0]
  YBox = box_size[1]*[0,0,1,1,0,0,1,1]+box_xyz0[1]
  ZBox = box_size[2]*[0,0,0,0,1,1,1,1]+box_xyz0[2]
  BBox = [5,0,1,3,2,0, 5,4,5,7,6,4, 2,0,4, 2,1,5, 2,2,6, 2,3,7, -1]
  oBBox = obj_new('IDLgrPolyLine', XBox, YBox, ZBox, POLYLINES=BBox, $
                   THICK=bboxthick,  $
                   HIDE=1-drawbbox,  $
                   PALETTE=oBBoxPalette ,  $
                   COLOR=bboxfgcolor)
  oDatabox = obj_new('IDLgrPolyLine', XBox, YBox, ZBox, POLYLINES=BBox, $
                   THICK=bboxthick,  $
                   HIDE=1-drawdatabox,  $
;                  PALETTE=oBBoxPalette ,  $
                   COLOR=bboxfgcolor)
  oXTitle= obj_new('IDLgrText',xtitle)
  oYTitle= obj_new('IDLgrText',ytitle)
  oZTitle= obj_new('IDLgrText',ztitle)
  oXAxis = obj_new('IDLgrAxis',0, $
                   XCOORD_CONV=xcc_axis, $
                   YCOORD_CONV=ycc_axis, $
                   ZCOORD_CONV=zcc_axis, $
                   TITLE=oXTitle,HIDE=hideXAxis,   $
                   TICKLEN=(0.01-xcc_axis[0])/xcc_axis[1], $
                   EXTEND=0, $
                   RANGE=xrange, $
                   LOCATION=axes_intersect, $
                   _extra=extra)
  oYAxis = obj_new('IDLgrAxis',1, $
                   XCOORD_CONV=xcc_axis, $
                   YCOORD_CONV=ycc_axis, $
                   ZCOORD_CONV=zcc_axis, $
                   TITLE=oYTitle,HIDE=hideYAxis,   $
                   TICKLEN=(0.01-ycc_axis[0])/ycc_axis[1], $
                   EXTEND=0, $
                   RANGE=yrange, $
                   LOCATION=axes_intersect, $
                   _extra=extra)
  oZAxis = obj_new('IDLgrAxis',2, $
                   XCOORD_CONV=xcc_axis, $
                   YCOORD_CONV=ycc_axis, $
                   ZCOORD_CONV=zcc_axis, $
                   TICKLEN=(0.01-zcc_axis[0])/zcc_axis[1], $
                   TITLE=oZTitle,HIDE=hideZAxis,   $
                   EXTEND=0, $
                   RANGE=zrange, $
                   LOCATION=axes_intersect, $
                   _extra=extra)

  oMainView  = obj_new('IDLgrView', projection=2, color=bgcolor) ;, depth_cue=[0.,0.], zclip=[10.,-10.])

  wBase = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, TITLE="VizIt", $
    /TLB_KILL_REQUEST_EVENTS,/TLB_SIZE_EVENTS)
  wControl = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, TITLE="VizIt Controls", $
    /TLB_KILL_REQUEST_EVENTS)
  wButtonBars = WIDGET_BASE(wControl, /COLUMN, UVALUE='BUTTONS')
  wButtons1 = WIDGET_BASE(wButtonBars, /ROW, UVALUE='BUTTONS')
  wButtons2 = WIDGET_BASE(wButtonBars, /ROW, UVALUE='BUTTONS')
  wLAB  = WIDGET_LABEL (wButtons1, VALUE='snapshots:')
  wGif  = WIDGET_BUTTON(wButtons1, VALUE='idl.gif', UVALUE='GIF')
  wJPEG = WIDGET_BUTTON(wButtons1, VALUE='idl.jpg', UVALUE='JPEG')
  wVRML = WIDGET_BUTTON(wButtons1, VALUE='idl.wrl', UVALUE='VRML')
  wPS   = WIDGET_BUTTON(wButtons1, VALUE='idl.ps', UVALUE='PS')
  wTRANS   = WIDGET_BUTTON(wButtons1, VALUE='TMat', UVALUE='TRANS')
  wREDRAW   = WIDGET_BUTTON(wButtons2, VALUE='Redraw', UVALUE='REDRAW')
  wAUTO   = CW_BGROUP(wButtons2, ['auto'], ROW=1, /NONEXCLUSIVE, /RETURN_NAME)
  wSHOWWHAT   = CW_BGROUP(wButtons2, ['bbox','content'], COLUMN=1, /NONEXCLUSIVE, /RETURN_NAME)
  wAXES   = CW_BGROUP(wButtons2, ['xaxis','yaxis','zaxis'], COLUMN=1, /NONEXCLUSIVE, /RETURN_NAME)
  wDEPTHORDER   = WIDGET_BUTTON(wButtons2, VALUE='Depth Reorder', UVALUE='DEPTHORDER')
  wHIDDEN  =  CW_BGROUP(wButtons2, ['hidden_line'], ROW=1, /NONEXCLUSIVE, /RETURN_NAME)
  wSTYLE  =  CW_BGROUP(wButtons2, ['points','wireframe','solid'], COLUMN=1, /EXCLUSIVE, /RETURN_NAME)
  wSHADING  =  CW_BGROUP(wButtons2, ['flat','Gouraud'], COLUMN=1, /EXCLUSIVE, /RETURN_NAME)
  if nlevels eq 1 then begin
    wLAB  = WIDGET_LABEL (wButtons2, UVALUE='LEVEL', UNAME='LEVEL', $
      VALUE=' level = '+string(level,format='(g10.3)'))
  endif else begin
    strlevels=''
    for eachlevel=0,nlevels-1 do begin
      strlevels=strjoin([strlevels, string(level[eachlevel],format='(g10.3)'),', '])
    endfor
    wLAB  = WIDGET_LABEL (wButtons, UVALUE='LEVEL', UNAME='LEVEL', $
        VALUE=strjoin([' level = ',strlevels]))
  endelse
  wDraw = WIDGET_DRAW(wBase, XSIZE=xsz, YSIZE=ysz, UVALUE='MAIN', $
    RETAIN=0, /EXPOSE_EVENTS, /BUTTON_EVENTS, GRAPHICS_LEVEL=2)
  wZClip0 = CW_FIELD( wControl , /FLOATING , /RETURN_EVENTS, /ROW, STRING=string, TEXT_FRAME=2, TITLE='Near Z Clip', UVALUE='ZCLIP0')
  wZClip1 = CW_FIELD( wControl , /FLOATING , /RETURN_EVENTS, /ROW, STRING=string, TEXT_FRAME=2, TITLE='Far Z Clip', UVALUE='ZCLIP1')

if keyword_set(PNG) then begin
  oMainWindow = obj_new('IDLgrBuffer', $
                   DIMENSIONS=[xsz,ysz])
endif else begin
  WIDGET_CONTROL, wControl, /REALIZE
  WIDGET_CONTROL, wBase, /REALIZE
  WIDGET_CONTROL, wDraw, GET_VALUE=oMainWindow
  WIDGET_CONTROL, wDraw, /DRAW_MOTION
endelse

  oModel = obj_new('IDLgrModel',CLIP_PLANES=clip_planes)
  oDataboxModel = obj_new('IDLgrModel',CLIP_PLANES=clip_planes)
  if nsurfaces gt 0 then begin
    for eachsurface=0,nsurfaces-1 do begin
      oModel->add, oSurfaces[eachsurface]
;      print,"Adding to plot...",eachsurface
    endfor
  endif
if not keyword_set(noiso) then begin
  if nlevels eq 1 then begin
    oModel->add, oPolygon
  endif else begin
    for eachlevel=0,nlevels-1 do oModel->add, oPolygon[eachlevel]
  endelse
endif
  oModel->add, oBBox
  oModel->add, oXAxis
  oModel->add, oYAxis
  oModel->add, oZAxis
  oMainView->add, oModel
  oDataboxModel->add, oDatabox
  oMainView->add, oDataboxModel
  get_bounds, oModel, xr, yr, zr
  oMainView->setproperty, dept=[zr(0),zr(1)+2.5*(zr(1)-zr(0))]
;  oMainView->setproperty, dept=[0.,0.]
;  oMaINvIEw->setproperty, zclip=[20.,-20.]
  set_view, oMainView, oMainWindow, /isotr


  oLight0 = obj_new('IDLgrLight', TYPE=0, INTENSITY=0.50)
  oLight1 = obj_new('IDLgrLight', TYPE=1, INTENSITY=0.75, LOCATION=[ 1, 1, 0])
  oLight2 = obj_new('IDLgrLight', TYPE=1, INTENSITY=0.95, LOCATION=[-1, 1, 1])
  oLights = obj_new('IDLgrModel')
  oLights->add, oLight0
  oLights->add, oLight1
  oLights->add, oLight2
  oMainView->add, oLights

  scl=0.95
  oModel->scale, scl,scl,scl
;  oModel->rotate, [1,1,0], 15
  if keyword_set(tmat) then begin
    if n_elements(tmat) eq 16 then begin
      oModel->SetProperty, TRANSFORM=tmat
    endif else if ((size(tmat))[0]) eq 3 then begin
      if n_elements(tmat[0,*,*]) eq 16 then begin
        oModel->SetProperty, TRANSFORM=reform(tmat[0,*,*])
      endif
    endif
  endif
  yszPlot=64

  if (n_elements(zclip) eq 2) then $
    oMainView->SetProperty, ZCLIP=zclip
  oMainView->GetProperty, ZCLIP=zclipget
  WIDGET_CONTROL, wZClip0, set_value=zclipget[0]
  WIDGET_CONTROL, wZClip1, set_value=zclipget[1]

;  vizit_BBoxFade,oBBox,oModel,bboxfgcolor=bboxfgcolor,bboxbgcolor=bboxbgcolor
  vizit_BBoxFade,oBBox,oModel

if not keyword_set(PNG) then begin
   wPlot = WIDGET_DRAW(wControl, XSIZE=xsz, YSIZE=yszPlot, UVALUE='HIST', $
    RETAIN=0, /EXPOSE_EVENTS, /BUTTON_EVENTS, GRAPHICS_LEVEL=2)
  WIDGET_CONTROL, wPlot, /DRAW_MOTION
  WIDGET_CONTROL, wPlot, GET_VALUE=oPlotWindow

  hist = alog10(histogram(data, min=min, max=max, binsize=(max-min)/255.)+0.1)
  xcc = [0.,1./n_elements(hist)] & ycc = [0.,1./max(hist)]
  oPlot  = obj_new('IDLgrPlot', hist, color=255, XCOORD=xcc, YCOORD=ycc)
  index = byte(255.*plevel)
  oLine  = obj_new('IDLgrPlot', [0,max(hist)], DATAX=[index,index], $
    color=255, XCOORD=xcc, YCOORD=ycc)
  oPlotModel = obj_new('IDLgrModel')
  image = bytscl(rebin(reform(findgen(xsz),xsz,1),xsz,yszPlot))
  oPlotImage = obj_new('IDLgrImage', image, PALETTE=oPalette, $
    XCOORD=[0.,1./xsz], YCOORD=[0.,1./yszPlot])
  oPlotView  = obj_new('IDLgrView', viewplane=[0.,0.,1.,1.], color=[50,50,50])
  oPlotModel->add, oPlotImage
  oPlotModel->add, oPlot
  oPlotModel->add, oLine
  oPlotView->add, oPlotModel
  oPlotWindow->draw, oPlotView
endif else begin
  oPlot=0.
  oPlotView=0.
  oPlotWindow=0.
  oPlotImage=0.
  oPlotModel=0.
  oLine=0.
  wPlot=0.
endelse

  oTrack = obj_new('vizit_trackball', [xsz/2.,ysz/2.], (xsz < ysz)/2.)

  oHolder = obj_new('IDL_Container')
  oHolder->add, oMainView
if not keyword_set(PNG) then oHolder->add, oPlotView
  oHolder->add, oTrack
drawbbox=1
drwacontent=1

  default,xx,0
  default,yy,0
  default,zz,0
  sState = {  oHolder:     oHolder,      $
              oTrack:      oTrack,       $
              oBBox:       oBBox,        $
              oXAxis:      oXAxis,       $
              oYAxis:      oYAxis,       $
              oZAxis:      oZAxis,       $
              oXTitle:      oXTitle,       $
              oYTitle:      oYTitle,       $
              oZTitle:      oZTitle,       $
              oMainView:   oMainView,    $
              oMainWindow: oMainWindow,  $
              wZClip0:     wZClip0,      $
              wZClip1:     wZClip1,      $
              wBase:       wBase,        $
              wControl:    wControl,     $
              wDraw:       wDraw,        $
              oPlot:       oPlot,        $
              oPlotView:   oPlotView,    $
              oPlotImage:  oPlotImage,   $
              oPlotWindow: oPlotWindow,  $
              oPolygon:    oPolygon,     $
              oSurfaces:   oSurfaces,    $
              oLine:       oLine,        $
              oLight0:     oLight0,      $
              oLight1:     oLight1,      $
              oLight2:     oLight2,      $
              oLights:     oLights,      $
              oModel:      oModel,       $
              oPalette:    oPalette,     $
              oPlotModel:  oPlotModel,   $
              oTexture:    oTexture,     $
              wPlot:       wPlot,        $
              xx:          xx,           $
              yy:          yy,           $
              zz:          zz,           $
              nonequidistant: nonequidistant, $
              calcnormals:  calcnormals, $
              xsize:       xsz,          $
              ysize:       yszPlot,      $
              min:         min,          $
              max:         max,          $
              data:        ptr_new(data),$
              level:       level,        $
              index:       index,        $
;              pVector:     pVector,     $
              autodraw:      1,          $
              hideXAxis:     hideXaxis,  $
              hideYAxis:     hideYaxis,  $
              hideZAxis:     hideZaxis,  $
              drawbbox:      drawbbox,    $
              drawcontent:   drawcontent, $
              bboxfgcolor:   bboxfgcolor,  $
              bboxbgcolor:   bboxbgcolor,  $
              noiso:  keyword_set(noiso), $
              nsurfaces:    nsurfaces,   $
              type:        'vizit',      $
              nlevels:      nlevels,     $
              active:      0             $
             }

  if (keyword_set(PNG) or keyword_set(GIF) or keyword_set(JPG)) then begin
    default,filename,'idl'
    if (drawbbox eq 1) then vizit_BBoxFade, oBBox, oModel
    if keyword_set(tmat) then begin
      if n_elements(tmat) eq 16 then begin
        oModel->SetProperty, TRANSFORM=tmat
        oMainWindow->draw, oMainView
        oImage = oMainWindow->read()
        oImage->GetProperty, data=im
        image2D = Color_Quan(im, 1, r, g, b)
        if n_elements(frameno) eq 1 then begin
            frame_filename=filename+string(frameno,format=format_frameno)
            frameno=frameno+1
        endif else begin
            frame_filename=filename
        endelse
        if keyword_set(PNG) then begin
          write_png, frame_filename+'.png', image2d, r, g, b
        endif else if keyword_set(GIF) then begin
          write_gif, frame_filename+'.gif', image2d, r, g, b
        endif else if keyword_set(JPG) then begin
          write_jpeg2000, frame_filename+'.jpg', image2d, r, g, b
        endif
      endif else if ((size(tmat))[0]) eq 3 then begin
        if n_elements(tmat[0,*,*]) eq 16 then begin
          default,frameno,0L
          totalframes=(size(tmat))[1]
          nextinc=10
          for frame=0,totalframes-1 do begin
            oModel->SetProperty, TRANSFORM=reform(tmat[frame,*,*])
            if keyword_set(DEPTH_REORDER) then begin
              if (not keyword_set(noiso)) then begin
                for i=0L,n_elements(oPolygon)-1 do begin
                  reorder_by_depth,oPolygon[i],tmat=reform(tmat[frame,*,*])
                endfor
              endif
              for i=0L,n_elements(oSurfaces)-1 do begin
                reorder_by_depth,oSurfaces[i],tmat=reform(tmat[frame,*,*])
              endfor
            endif
            oMainWindow->draw, oMainView
            oImage = oMainWindow->read()
            oImage->GetProperty, data=im
            image2D = Color_Quan(im, 1, r, g, b)

            frame_filename=filename+string(frameno,format=format_frameno)
            if keyword_set(PNG) then begin
              write_png, frame_filename+'.png', image2d, r, g, b
            endif else if keyword_set(GIF) then begin
              write_gif, frame_filename+'.gif', image2d, r, g, b
            endif else if keyword_set(JPG) then begin
              write_jpeg2000, frame_filename+'.jpg', image2d, r, g, b
            endif
            frameno=frameno+1
            pcent=totalframes/frameno
            if pcent ge nextinc then begin
              print,nextinc,'%'
              nextinc=nextinc+10
            endif
          endfor
        endif
      endif
    endif
    OBJ_DESTROY, oHolder
    WIDGET_CONTROL, wControl, /DESTROY
    WIDGET_CONTROL, wBase, /DESTROY
    if ptr_valid(sState.data) then ptr_free,sState.data
    SAFE_OBJ_DESTROY,oMainWindow
    SAFE_OBJ_DESTROY,oMainView
    SAFE_OBJ_DESTROY,oPolygon
    SAFE_OBJ_DESTROY,oSurfaces
    SAFE_OBJ_DESTROY,oBBox
    SAFE_OBJ_DESTROY,oImage
    SAFE_OBJ_DESTROY,oLight0
    SAFE_OBJ_DESTROY,oLight1
    SAFE_OBJ_DESTROY,oLight2
    SAFE_OBJ_DESTROY,oLights
    SAFE_OBJ_DESTROY,oModel
    SAFE_OBJ_DESTROY,oPalette
    SAFE_OBJ_DESTROY,oPlot
    SAFE_OBJ_DESTROY,oPlotModel
    SAFE_OBJ_DESTROY,oTrack
    SAFE_OBJ_DESTROY,oXAxis
    SAFE_OBJ_DESTROY,oYAxis
    SAFE_OBJ_DESTROY,oZAxis
    SAFE_OBJ_DESTROY,oXTitle
    SAFE_OBJ_DESTROY,oYTitle
    SAFE_OBJ_DESTROY,oZTitle
    SAFE_OBJ_DESTROY,oPlotView
    SAFE_OBJ_DESTROY,oTexture
    SAFE_OBJ_DESTROY,oLine
    SAFE_OBJ_DESTROY,oPlotImage
  endif else begin
    WIDGET_CONTROL, wControl, SET_UVALUE=ptr_new(sState)
    WIDGET_CONTROL, wBase, SET_UVALUE=ptr_new(sState)
    XMANAGER, 'vizit_control', wControl, /NO_BLOCK
    XMANAGER, 'vizit', wBase, /NO_BLOCK
  endelse
END
