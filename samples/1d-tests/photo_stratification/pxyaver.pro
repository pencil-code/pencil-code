default,iread,0
if iread eq 0 then begin
  pc_read_xyaver,obj=xya
  iread=1
endif
;
plot,xya.z,xya.ttmz[*,0]
s=size(xya.ttmz)
nt=s[2]
;
for it=0,nt-1 do begin
  oplot,xya.z,xya.ttmz[*,it]
endfor
;
END
