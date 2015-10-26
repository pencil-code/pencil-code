PRO vizit_event, sEvent
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
        'MAIN': BEGIN

           CASE sEvent.type OF
           0: BEGIN                                    ; Mouse press
                (*sState).active = 1
                if (*sState).type eq 'vizit' then begin
                  loc = [sEvent.x, sEvent.y]
     ;             if (*sState).oMainWindow->PickData((*sState).oMainView,(*sState).oPolygon,loc,xyz) then begin
     ;               print,'PickData: xyz=',xyz
     ;             endif
                  if (not (*sState).noiso) then begin
                    if !d.name eq 'MAC' then $
                      (*sState).oModel->remove, (*sState).oPolygon
                  endif
                endif
                if (*sState).type eq 'volviz' then $
                  (*sState).oModel->remove, (*sState).oVolume
                if ((*sState).autodraw eq 1) then (*sState).oMainWindow->draw, (*sState).oMainView
                ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
           END
           1: BEGIN                                    ; Mouse release
                (*sState).active = 0
                if (not (*sState).noiso) then begin
                  if (*sState).type eq 'vizit' $
                  and !d.name eq 'MAC' then $
                    (*sState).oModel->add, (*sState).oPolygon
                endif
                if (*sState).type eq 'volviz' then $
                  (*sState).oModel->add, (*sState).oVolume
                if ((*sState).autodraw eq 1) then (*sState).oMainWindow->draw, (*sState).oMainView
                ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
           END
           4: BEGIN                                    ; Expose
               if ((*sState).drawbbox eq 1) then vizit_BBoxFade, (*sState).oBBox, (*sState).oModel
               if ((*sState).autodraw eq 1) then (*sState).oMainWindow->Draw, (*sState).oMainView
               (*sState).oPlotWindow->Draw, (*sState).oPlotView
               ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState
               RETURN
           END
           ELSE:
           END

           bHaveTransform = (*sState).oTrack->Update( sEvent, TRANSFORM=qmat )
           IF (bHaveTransform NE 0) THEN BEGIN           ; Trackball updates.
               (*sState).oModel->GetProperty, TRANSFORM=t
               (*sState).oModel->SetProperty, TRANSFORM=t#qmat
               if ((*sState).drawbbox eq 1) then vizit_BBoxFade, (*sState).oBBox, (*sState).oModel
               if ((*sState).autodraw eq 1) then (*sState).oMainWindow->Draw, (*sState).oMainView
           END
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
