function haver,f,xvertical=xvertical,yvertical=yvertical
;
;  Horizontal average of 3-d or 2-d scalar
;  $Id: haver.pro,v 1.2 2003-10-29 09:22:40 theine Exp $
;
        sizef=size(f)
        dim=sizef[0]
        if dim eq 3 then begin
          if keyword_set(xvertical) then begin
            h=fltarr(sizef[1])
            for n=0,sizef[1]-1 do h[n]=aver(f[n,*,*])
          endif else if keyword_set(yvertical) then begin
            h=fltarr(sizef[2])
            for n=0,sizef[2]-1 do h[n]=aver(f[*,n,*])
          endif else begin
            h=fltarr(sizef[3])
            for n=0,sizef[3]-1 do h[n]=aver(f[*,*,n])
          endelse
        endif else if dim eq 2 then begin
          if keyword_set(xvertical) then begin
            h=fltarr(sizef[1])
            for n=0,sizef[1]-1 do h[n]=aver(f[n,*])
          endif else begin
            h=fltarr(sizef[2])
            for n=0,sizef[2]-1 do h[n]=aver(f[*,n])
          endelse
        endif else begin
          print,"error: I don't know how to handle this array"
          return,0
        endelse
        return,h
end
