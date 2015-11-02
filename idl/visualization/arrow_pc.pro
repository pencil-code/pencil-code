   pro arrow_pc, x1, y1, x2, y2, _extra=extra, polar=polar

     if keyword_set(polar) then begin
       xa=x1*cos(y1) & ya=x1*sin(y1) & xe=x2*cos(y2) & ye=x2*sin(y2)
     endif else begin
       xa=x1 & ya=y1 & xe=x2 & ye=y2
     endelse

     arrow, xa, ya, xe, ye, _extra=extra

   end
