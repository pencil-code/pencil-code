function haver,f,xvertical=xvertical,yvertical=yvertical
;
;  Horizontal average of 3-d or 2-d scalar; and copy if 1-d scalar
;  $Id$
;
        sizef=size(f)
        dim=sizef[0]
        if dim eq 3 then begin
          if keyword_set(xvertical) then begin
            h=dblarr(sizef[1])
            for n=0,sizef[1]-1 do h[n]=aver(f[n,*,*])
          endif else if keyword_set(yvertical) then begin
            h=dblarr(sizef[2])
            for n=0,sizef[2]-1 do h[n]=aver(f[*,n,*])
          endif else begin
            h=dblarr(sizef[3])
            for n=0,sizef[3]-1 do h[n]=aver(f[*,*,n])
          endelse
        endif else if dim eq 2 then begin
          if keyword_set(xvertical) then begin
            h=dblarr(sizef[1])
            for n=0,sizef[1]-1 do h[n]=aver(f[n,*])
          endif else begin
            h=dblarr(sizef[2])
            for n=0,sizef[2]-1 do h[n]=aver(f[*,n])
          endelse
        endif else if dim eq 1 then begin
          h=f
        endif else begin
          print,"% HAVER: size(f) = ", sizef
          message,"Don't know how to handle array"
        endelse
        return,h
end
