; Check how many particles have collided with the walls

pro channel_deposition,npart_radii,$
                       irmv,trmv,$
                       objpvar=objpvar,startparam=param,$
                       savefile=savefile

dims=size(irmv)
colls=0

xrmv=dblarr(dims[1]) 
yrmv=dblarr(dims[1]) 
zrmv=dblarr(dims[1]) 

for k=0,dims[1]-1 do begin        
    if (dims[0]>0) then begin              
        xrmv(k)=objpvar.xx[irmv[k],0]
        yrmv(k)=objpvar.xx[irmv[k],1]
        zrmv(k)=objpvar.xx[irmv[k],2]
        colls = colls+1
    endif
endfor

print,'Number of collisions with the walls:',colls

;
; Loop over all particle diameters
;

if (savefile) then begin
    for i=0,npart_radii-1 do begin
        filename='./data/rmvdata.sav'
        if (npart_radii gt 1) then begin
            filename='./data/rmvdata'+str(startparam.ap0(i))+'.sav'           
        endif
        save,irmv,trmv,xrmv,yrmv,zrmv,filename=filename
    endfor  

    
endif

END
