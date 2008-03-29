PRO vizit_BBoxFade, oBBox, oModel
;+
;  Fade the intensity of the Bounding Box according to projected z-values
;-
    oBBox->GetProperty, DATA=vert
    oModel->GetProperty, TRANSFORM=t
    nv = (size(vert))[2]
    z = fltarr(nv)
    for i=0,nv-1 do z(i)=([vert[*,i],1.]#t)[2]
    min = 120
    z = min+bytscl(z,top=255-min)
    oBBox->SetProperty, VERT_COLORS=z, /SHADING
END
