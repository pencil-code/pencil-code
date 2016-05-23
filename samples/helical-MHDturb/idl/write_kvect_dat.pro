;
;  procedure to write k and e vectors together with amplitudes and phases
;
ampl=gaunoise(n,seed=seed)
eex=gaunoise(n,seed=seed)
eey=gaunoise(n,seed=seed)
eez=gaunoise(n,seed=seed)
phik=2.*!pi*(random(n,seed=seed)-.5)
close,1
openw,1,'kvect.dat'
fo='(3f7.3,2x,3f7.3,2x,2f11.6)'
fo2='(i6,a,2f7.2)'
;
printf,1,n,' vectors within the range',k1,k2,fo=fo2
printf,1,'    kx     ky     kz       ex     ey     ez        rk        phik'
for i=0,n-1 do begin
  kxe2=(kky(i)*eez(i)-kkz(i)*eey(i))^2 $
      +(kkz(i)*eex(i)-kkx(i)*eez(i))^2 $
      +(kkx(i)*eey(i)-kky(i)*eex(i))^2
  if kxe2 eq 0. then begin
    stop,'bad'
  endif else begin
    printf,1,kkx(i),kky(i),kkz(i),eex(i),eey(i),eez(i),ampl(i),phik(i),fo=fo
  endelse
endfor
close,1
END
