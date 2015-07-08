;
;  procedure to write k and e vectors together with amplitudes and phases
;
ampl=gaunoise(n)
jj=fix(n*random(n))
phik=2.*!pi*(random(n)-.5)
openw,1,'kvect.dat'
fo='(3i5,2x,3i5,2x,2f11.6)'
fo2='(i6,a,2f7.2)'
;
printf,1,n,' vectors within the range',k1,k2,fo=fo2
printf,1,'   kx   ky   kz     ex   ey   ez        rk        phik'
for i=0,n-1 do begin
  j=jj(i)
  kxe2=(kky(i)*kkz(j)-kkz(i)*kky(j))^2 $
      +(kkz(i)*kkx(j)-kkx(i)*kkz(j))^2 $
      +(kkx(i)*kky(j)-kky(i)*kkx(j))^2
  if kxe2 eq 0. then begin
    if kkx(j) ne 0. then begin
      printf,1,kkx(i),kky(i),kkz(i),-kkx(j),kky(j),kkz(j),ampl(i),phik(i),fo=fo
    endif else if kky(j) ne 0. then begin
      printf,1,kkx(i),kky(i),kkz(i),kkx(j),-kky(j),kkz(j),ampl(i),phik(i),fo=fo
    endif else if kkz(j) ne 0. then begin
      printf,1,kkx(i),kky(i),kkz(i),kkx(j),kky(j),-kkz(j),ampl(i),phik(i),fo=fo
    endif
  endif else begin
    printf,1,kkx(i),kky(i),kkz(i),kkx(j),kky(j),kkz(j),ampl(i),phik(i),fo=fo
  endelse
endfor
close,1
END
