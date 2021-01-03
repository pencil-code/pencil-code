;$Id$
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0,nghostx,nghosty,nghostz
;
;  Read xyaver content and write coarsegrain data to save file
;
;  30-mar-08/axel
;
@xder_6th
pc_read_kf,kf
pc_read_xyaver,obj=xyaver
pc_read_param,obj=param,/param2
pc_read_param,obj=param1
z0=param1.xyz0(2)
z1=param1.xyz1(2)
Lz=z1-z0
;
;  put isample=1 to prevent averaging between nevery intervals
;
default,use_grid,1
default,isample,0
default,nevery,50
spawn,'touch parameters.pro'
@parameters
;
if use_grid eq 1 then begin
  pc_read_grid,obj=grid
  pc_read_dim,obj=dim
  nz=dim.nz & n1=dim.n1 & n2=dim.n2
  zzz=grid.z(n1:n2)
endif else begin
  ;
  ;  read in z-array
  ;
  pc_read_dim,obj=dim
  default,nz,dim.nz
  dz=Lz/(nz-1)
  print,'assume non-periodic domain with dz=',dz
  zzz=z0+dz*findgen(nz)
endelse
;
;  check content of xyaver.in and make each variable defined
;
a=''
iread=0
openr,1,'xyaver.in'
openw,2,'.xyaver'
while not eof(1) do begin
  readf,1,a
  if isample eq 0 then begin
    printf,2,'pc_coarsegrain,xyaver.t,transpose(xyaver.'+a+'),nevery,t,'+a+',/aver,/first'
  endif else begin
    printf,2,'pc_coarsegrain,xyaver.t,transpose(xyaver.'+a+'),nevery,t,'+a+',/first'
  endelse
  if iread eq 0 then aa=a else aa=aa+','+a
  iread=iread+1
endwhile
print
print,'print save file line'
printf,2,"save,file='xyaver_coarsegrain_'+str(nevery)+'.sav',zzz,t,nevery,kf,"+aa
print,"save,file='xyaver_coarsegrain_'+str(nevery)+'.sav',zzz,t,nevery,kf,"+aa
close,1
close,2
;print
;print,'wait 3 seconds...'
;print
wait,.3

@.xyaver
print
print,'Note: this procedure just wrote a temporary idl script in .xyaver,'
print,'but after the first execution this new script may not yet be recognized.'
print,"This would be the case if the variable t below doesn't exist."
print,'In that case, just run this script again."
print
print,'t',t
;
;spawn,'mv .xyaver .xyaver_old'
END
