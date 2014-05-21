function logformat, axis, index, val
 if !x.type then begin
   if val ge 1 and val le 1000 then begin
     if val lt 10 then $
       return, string(val,format='(i1)') $
     else if val lt 100 then $
       return, string(val,format='(i2)') $
     else if val lt 1000 then $
       return, string(val,format='(i3)') $
     else $
       return, string(val,format='(i4)')
   endif else $ 
     return, strtrim(string(val,format='(g10.1)'),2)
 endif else $
   return, strtrim(string(val,format='(g10.1)'),2)
end

