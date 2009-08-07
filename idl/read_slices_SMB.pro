;This is a short a program that reads a slice file. 
;Will be clean up
;A. Svedin Augst 2009

islice=0
file_slice='./data/proc0/slice_uu1.xz'	
global_slice=fltarr(64,64,88)
timer=fltarr(88)
	openr,1,file_slice,/f77	


for i=0,87 do begin & print,islice & plane=fltarr(64,64)*1 & readu,1,plane,t,slice_z2po & 	islice=islice+1 & surface,plane & global_slice(*,*,i)=plane & timer(i)=t

;end 

close,1

