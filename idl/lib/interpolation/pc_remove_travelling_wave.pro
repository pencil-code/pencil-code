;$Id$
pro pc_remove_travelling_wave,f,z,speed,t, $
  average=average,it1=it1
;
;  24-jun-2008/axel: coded
;
;  Assume that we have an array of the form f(z,t,i,j,...)
;  Uses the routine pc_shift_6th of Anders
;
s=size(f)
;
n=s[0]
nz=s[1]
nt=s[2]
;
;  shift by c*t, where c is the speed
;  positive when travelling in the positive z direction
;
print,'n=',n
if n eq 2 then begin
  for it=0,nt-1 do f[*,it]=pc_shift_6th(f[*,it],z,-speed*t(it))
endif else if n eq 3 then begin
  n3=s[3]
  for i=0,n3-1 do begin
    for it=0,nt-1 do f[*,it,i]=pc_shift_6th(f[*,it,i],z,-speed*t(it))
  endfor
endif else if n eq 4 then begin
  n3=s[3]
  n4=s[4]
  print,'n3,n4=',n3,n4
  for j=0,n4-1 do begin
    for i=0,n3-1 do begin
      for it=0,nt-1 do f[*,it,i,j]=pc_shift_6th(f[*,it,i,j],z,-speed*t(it))
    endfor
  endfor
end
;
;  calculate time average (if requested)
;  start with a later timestep (if requested)
;
if arg_present(average) then begin
  if keyword_set(it1) then it1_tmp=it1 else it1_tmp=0
  ntgood=n_elements(t(it1:*))
  if n eq 2 then begin
    average=total(f[*,it1:*],2)/ntgood
  endif else if n eq 3 then begin
    average=total(f[*,it1:*,*],2)/ntgood
  endif else if n eq 4 then begin
    average=total(f[*,it1:*,*,*],2)/ntgood
  end
;
endif
;
end
