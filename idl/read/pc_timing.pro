;$Id$
;
;  read timing file
;
pro pc_timing
;
t=0d0
tprev=0d0
msg=''
counter=0
mul_fac=1
;
close,9
openr,9,'data/timing.dat'
while not eof(9) do begin
  readf,9,t,mul_fac,msg
  if counter eq 0 then begin
    tt=t
    
  endif else begin
    tt=[tt,t] 
  endelse
  fo='(i5,f10.5,2x,a)'
  print,counter,mul_fac*(t-tprev),msg,fo=fo
  counter=counter+1
  tprev=t
endwhile
close,9
;
plot,tt,ps=-6
;
END
