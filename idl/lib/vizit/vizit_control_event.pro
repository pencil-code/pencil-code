PRO vizit_set_surface_properties,sState,_EXTRA=_EXTRA
  if (not (*sState).noiso) then begin
    if n_elements((*sState).oPolygon) eq 1 then begin
      (*sState).oPolygon->SetProperty,_EXTRA=_EXTRA
    endif else begin
      for i=0,n_elements((*sState).oPolygon)-1 do begin
        (*sState).oPolygon[i]->SetProperty,_EXTRA=_EXTRA
      endfor
    endelse
  endif

  if (has_tag(*sState, 'oSurfaces')) then begin
    if ((*sState).nsurfaces gt 0) then begin
      for i=0,n_elements((*sState).oSurfaces)-1 do begin
        (*sState).oSurfaces[i]->SetProperty,_EXTRA=_EXTRA
      endfor
    endif
  endif
end

PRO vizit_reorder_by_depth_all_surfaces,sState,_EXTRA=_EXTRA
  if (not (*sState).noiso) then begin
    if n_elements((*sState).oPolygon) eq 1 then begin
      reorder_by_depth,(*sState).oPolygon,_EXTRA=_EXTRA
    endif else begin
      for i=0,n_elements((*sState).oPolygon)-1 do begin
        reorder_by_depth,(*sState).oPolygon[i],_EXTRA=_EXTRA
      endfor
    endelse
  endif

  if (has_tag(*sState, 'oSurfaces')) then begin
    if ((*sState).nsurfaces gt 0) then begin
      for i=0,n_elements((*sState).oSurfaces)-1 do begin
        reorder_by_depth,(*sState).oSurfaces[i],_EXTRA=_EXTRA
      endfor
    endif
  endif
end


PRO vizit_control_event, sEvent
;+
;  Handle events for VIZIT
;-
    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState    ; Get state
    ;print, sEvent

    CASE TAG_NAMES(sEvent, /STRUCTURE_NAME) OF       ; Kill request
    'WIDGET_KILL_REQUEST': BEGIN
        IF not ptr_valid(sState) THEN return
        IF n_elements(*sState) eq 0 THEN return
        OBJ_DESTROY, (*sState).oHolder
;        WIDGET_CONTROL, sEvent.top, /DESTROY
        WIDGET_CONTROL, (*sState).wControl, /DESTROY
        WIDGET_CONTROL, (*sState).wBase, /DESTROY
        if ptr_valid((*sState).data) then ptr_free,(*sState).data
        SAFE_OBJ_DESTROY,(*sState).oMainWindow
        SAFE_OBJ_DESTROY,(*sState).oMainView
        SAFE_OBJ_DESTROY,(*sState).oPolygon
        SAFE_OBJ_DESTROY,(*sState).oSurfaces
        SAFE_OBJ_DESTROY,(*sState).oBBox
        SAFE_OBJ_DESTROY,(*sState).oLight0
        SAFE_OBJ_DESTROY,(*sState).oLight1
        SAFE_OBJ_DESTROY,(*sState).oLight2
        SAFE_OBJ_DESTROY,(*sState).oLights
        SAFE_OBJ_DESTROY,(*sState).oModel
        SAFE_OBJ_DESTROY,(*sState).oPalette
        SAFE_OBJ_DESTROY,(*sState).oPlot
        SAFE_OBJ_DESTROY,(*sState).oPlotModel
        SAFE_OBJ_DESTROY,(*sState).oTrack
        SAFE_OBJ_DESTROY,(*sState).oXAxis
        SAFE_OBJ_DESTROY,(*sState).oYAxis
        SAFE_OBJ_DESTROY,(*sState).oZAxis
        SAFE_OBJ_DESTROY,(*sState).oXTitle
        SAFE_OBJ_DESTROY,(*sState).oYTitle
        SAFE_OBJ_DESTROY,(*sState).oZTitle
        SAFE_OBJ_DESTROY,(*sState).oPlotView
        SAFE_OBJ_DESTROY,(*sState).oTexture
        SAFE_OBJ_DESTROY,(*sState).oLine
        SAFE_OBJ_DESTROY,(*sState).oPlotImage
	PTR_FREE,sState
        RETURN
    END
    'WIDGET_BASE': BEGIN                             ; Resize request
        xsz = sEvent.x
        ysz = sEvent.y

        ;(*sState).oPlotWindow->SetProperty, DIM=[xsz,(*sState).ysize]

        (*sState).oMainWindow->SetProperty, DIM=[xsz,ysz]
        (*sState).oMainView->Remove, (*sState).oLights
        set_view, (*sState).oMainView, (*sState).oMainWindow, /ISOTR
        (*sState).oMainView->Add, (*sState).oLights
        if ((*sState).autodraw eq 1) then  (*sState).oMainWindow->Draw, (*sState).oMainView
        (*sState).oPlotWindow->Draw, (*sState).oPlotView

        (*sState).xsize = xsz
        ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
        RETURN
    END
    ELSE:
    ENDCASE


     if (WIDGET_INFO(sEvent.id,/UNAME) eq 'CW_BGROUP_UNAME') then begin
      IF (sEvent.Value eq 'auto') then begin
        (*sState).autodraw=sEvent.select
      endif else IF (sEvent.Value eq 'bbox') then begin
        (*sState).drawbbox=sEvent.select
      endif else IF (sEvent.Value eq 'content') then begin
        (*sState).drawcontent=sEvent.select
      endif else IF (sEvent.Value eq 'xaxis') then begin
        (*sState).hideXAxis=1-sEvent.select
        (*sState).oXAxis->SetProperty,HIDE=(*sState).hideXAxis
      endif else IF (sEvent.Value eq 'yaxis') then begin
        (*sState).hideYAxis=1-sEvent.select
        (*sState).oYAxis->SetProperty,HIDE=(*sState).hideYAxis
      endif else IF (sEvent.Value eq 'zaxis') then begin
        (*sState).hideZAxis=1-sEvent.select
        (*sState).oZAxis->SetProperty,HIDE=(*sState).hideZAxis
      endif else IF (sEvent.Value eq 'hidden_lines') then begin
        vizit_set_surface_properties,sState,HIDDEN_LINES=1
      endif else IF (sEvent.Value eq 'points') then begin
        if sEvent.select then vizit_set_surface_properties,sState,STYLE=0
      endif else IF (sEvent.Value eq 'wireframe') then begin
        if sEvent.select then vizit_set_surface_properties,sState,STYLE=1
      endif else IF (sEvent.Value eq 'solid') then begin
        if sEvent.select then vizit_set_surface_properties,sState,STYLE=2
      endif else IF (sEvent.Value eq 'flat') then begin
        if sEvent.select then vizit_set_surface_properties,sState,SHADING=0
      endif else IF (sEvent.Value eq 'Gouraud') then begin
        if sEvent.select then vizit_set_surface_properties,sState,SHADING=1
      endif
      ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState

      (*sState).oBBox->SetProperty,HIDE=1-(*sState).drawbbox
      if (has_tag(*sState, 'oVolume')) then begin
        (*sState).oVolume->SetProperty,HIDE=1-(*sState).drawcontent
      endif else if (has_tag(*sState, 'oPolygon')) then begin
        if (not (*sState).noiso) then begin
          if n_elements((*sState).oPolygon) eq 1 then begin
            (*sState).oPolygon->SetProperty,HIDE=1-(*sState).drawcontent
          endif else begin
            for i=0,n_elements((*sState).oPolygon)-1 do begin
              (*sState).oPolygon[i]->SetProperty,HIDE=1-(*sState).drawcontent
            endfor
          endelse
        endif
      endif
      if (has_tag(*sState, 'oSurfaces')) then begin
        if ((*sState).nsurfaces gt 0) then begin
          for i=0,n_elements((*sState).oSurfaces)-1 do begin
            (*sState).oSurfaces[i]->SetProperty,HIDE=1-(*sState).drawcontent
          endfor
        endif
      endif
    endif else begin
    WIDGET_CONTROL, sEvent.id, GET_UVALUE=uval
    CASE uval OF
        'HIST': BEGIN                                    ; Vizit
           CASE sEvent.type OF
           0: (*sState).active = 1
           1: (*sState).active = 0
           ELSE:
           ENDCASE
           IF (*sState).active EQ 1 THEN BEGIN
              if (not (*sState).noiso) then begin
                index = long(255.*((sEvent.x+2) > 0)/(*sState).xsize) < 255
                level = (*sState).min + index/255.*((*sState).max-(*sState).min)
                (*sState).level = level
                (*sState).index = index
;                index = (*sState).index
;                level = (*sState).level
                shade_volume, *((*sState).data), level, vert, poly
                if n_elements(vert) le 1 then begin
                  vert=fltarr(3,3)
                  vert[*,*]=0.
                  poly=[3,0,1,2]
                endif
                if (*sState).calcnormals eq 1 then begin
                  gradvar=grad(*(*sState).data)
                  normals=vert
                  normals[0,*]=interpolate(gradvar[*,*,*,0],vert[0,*],vert[1,*],vert[2,*])
                  normals[1,*]=interpolate(gradvar[*,*,*,1],vert[0,*],vert[1,*],vert[2,*])
                  normals[2,*]=interpolate(gradvar[*,*,*,2],vert[0,*],vert[1,*],vert[2,*])
                endif
                if (*sState).nonequidistant eq 1 then begin
                  vert[0,*]=interpol((*sState).xx,findgen(n_elements((*sState).xx)),reform(vert[0,*]))
                  vert[1,*]=interpol((*sState).yy,findgen(n_elements((*sState).yy)),reform(vert[1,*]))
                  vert[2,*]=interpol((*sState).zz,findgen(n_elements((*sState).zz)),reform(vert[2,*]))
                endif
                (*sState).oPolygon->GetProperty,TEXTURE_MAP=oTexture
                oTexture->GetProperty,DATA=TextImg
                TextImg[0,*,*]=index
                oTexture->SetProperty,DATA=TextImg
                (*sState).oPolygon->SetProperty, DATA=vert, NORMALS=normals, POLYGONS=poly, COLOR=index, BOTTOM=index, TEXTURE_MAP=oTexture
                (*sState).oLine->SetProperty, DATAX=[index,index]
                (*sState).oPlotWindow->draw, (*sState).oPlotView
                if ((*sState).autodraw eq 1) then (*sState).oMainWindow->draw, (*sState).oMainView
                i = WIDGET_INFO(sEvent.top, FIND_BY_UNAME='LEVEL')
                if i gt 0 then WIDGET_CONTROL, i, $
                  set_value='level = '+str(level,format='(g10.3)')
              endif
           ENDIF
           ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
        END

        'OPAC': BEGIN                                    ; VolViz
           CASE sEvent.type OF
             0: BEGIN                                    ; Mouse press
                (*sState).active = 1
                (*sState).xprv = long(255.*(sEvent.x > 0)/(*sState).xsize) < 255
                (*sState).yprv = long(255.*(sEvent.y > 0)/(*sState).ysize) < 255
;                WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
             END
             1: BEGIN                                    ; Mouse release
                (*sState).active=0
                if ((*sState).autodraw eq 1) then (*sState).oMainWindow->draw, (*sState).oMainView
;               WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
             END
             2: BEGIN                                    ; Mouse motion
                if ((*sState).active eq 0) then return
                x = long(255.*(sEvent.x > 0)/(*sState).xsize) < 255
                y = long(255.*(sEvent.y > 0)/(*sState).ysize) < 255
                x0 = (*sState).xprv
                y0 = (*sState).yprv
                i0 = (x0 < x)
                i1 = (x > x0)
                FOR i = i0, i1 DO BEGIN
                    (*sState).opac[i]=y0+((y-y0)*(i-i0))/(i1-i0+1e-6)
                ENDFOR
                (*sState).xprv = x
                (*sState).yprv = y

                (*sState).oPlot->SetProperty, DATAY=(*sState).opac
                (*sState).oVolume->SetProperty, OPACITY_TABLE0=(*sState).opac
                (*sState).oPlotWindow->draw, (*sState).oPlotView

                ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
             END
             ELSE: BEGIN                                 ; Ignore others
             END
           ENDCASE
        END
        'DEPTHORDER': BEGIN
               (*sState).oModel->GetProperty, TRANSFORM=t
               vizit_reorder_by_depth_all_surfaces,sState,tmat=t
        END
        'ZCLIP0': BEGIN
               if (sEvent.update) then begin
                 (*sState).oMainView->GetProperty, ZCLIP=zclip
                 zclip[0]=sEvent.value
                 (*sState).oMainView->SetProperty, ZCLIP=zclip
                 if ((*sState).autodraw eq 1) then  (*sState).oMainWindow->Draw, (*sState).oMainView
               endif
        END
        'ZCLIP1': BEGIN
               if (sEvent.update) then begin
                 (*sState).oMainView->GetProperty, ZCLIP=zclip
                 zclip[1]=sEvent.value
                 (*sState).oMainView->SetProperty, ZCLIP=zclip
                 if ((*sState).autodraw eq 1) then  (*sState).oMainWindow->Draw, (*sState).oMainView
               endif
        END
        'TRANS': BEGIN
               (*sState).oModel->GetProperty, TRANSFORM=t
               print,"Transformation Matrix"
               help,t
               print,'tmat=[[',t[0,0],',',t[1,0],',',t[2,0],',',t[3,0],'], $'
               print,'      [',t[0,1],',',t[1,1],',',t[2,1],',',t[3,1],'], $'
               print,'      [',t[0,2],',',t[1,2],',',t[2,2],',',t[3,2],'], $'
               print,'      [',t[0,3],',',t[1,3],',',t[2,3],',',t[3,3],']]'


        END
        'GIF': BEGIN
               oImage = (*sState).oMainWindow->read()
               oImage->GetProperty, data=im
               write_gif, 'idl.gif', color_quan(im,1,r,g,b), r,g,b
        END
        'JPEG': BEGIN
               oImage = (*sState).oMainWindow->read()
               oImage->GetProperty, data=im
               write_jpeg, 'idl.jpg', im, /true
        END
        'PS':  BEGIN
               oPS = obj_new('IDLgrPrinter')
               diag_result=DIALOG_PRINTERSETUP(oPS)
               diag_result=DIALOG_PRINTJOB(oPS)
               oPS->draw, (*sState).oMainView
               obj_destroy, oPS
        END
        'VRML': BEGIN
               oVRML = obj_new('IDLgrVRML')
               oVRML->draw, (*sState).oMainView
               obj_destroy, oVRML
        END
        'REDRAW': BEGIN
                if ((*sState).drawbbox eq 1) then vizit_BBoxFade, (*sState).oBBox, (*sState).oModel
                (*sState).oPlotWindow->draw, (*sState).oPlotView
                (*sState).oMainWindow->draw, (*sState).oMainView
        END
        ELSE: BEGIN
           print, 'VIZIT_EVENT: Unexpected UVALUE = ',uval
        END
     ENDCASE
     endelse
END
