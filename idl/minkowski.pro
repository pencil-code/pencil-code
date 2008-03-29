pro shapefinders, arrayin, levelin

;**********************************************************
; a routine to caluclate minkowski funcitonals (see Wilkin 
; et al PRL, 2007) from 3D arrayin  at a particular iso level 
; (levelin) ghost zones are removed and then a binary array
; is created, from this the four minkowski functionals are 
; calculated. 
; A priodic version of this would be useful. 
;**********************************************************

;27/03/08 - abag coded

zero=0.
pi=3.1415926535897931

s=size(arrayin)

nx=s[1]-6
ny=s[2]-6
nz=s[3]-6

field=fltarr(nx,ny,nz)
field=(sqrt(dot2(arrayin)))[2:nx+2,2:ny+2,2:nz+2]
;contour,field[*,*,16],nlevels=64,/fill
print, nx
print, ny
print, nz
print, levelin

n0=0
n1=0
n2=0
n3=0

V0=0.
V1=0.
V2=0.
V3=0.

mask=bytarr(nx,ny,nz)
mask[*,*,*]=0

mask[1:nx-2,1:ny-2,1:nz-2]=field[1:nx-2,1:ny-2,1:nz-2] gt levelin ;this makes an array "mask" with zeros & 1s only

;find n0 - the total number of meshpoints within a structure
n0=TOTAL(mask)
print,"n0=",n0


;find n1 - the total number of edges of each small square inside a structure
for k=0,nz-2 do begin       ;covers full range minus 1
    for j=0,ny-2 do begin   ;because we treat the end point
        for i=0,nx-2 do begin  ;separately
            if TOTAL(mask[i,j,k]+mask[i+1,j,k]) eq 2 then begin
               n1=n1+1 
            endif
            if TOTAL(mask[i,j,k]+mask[i,j+1,k]) eq 2 then begin
               n1=n1+1
            endif
            if TOTAL(mask[i,j,k]+mask[i,j,k+1]) eq 2 then begin
               n1=n1+1
            endif

        end
    end
end

for k=0,nz-2 do begin
    for j=0,ny-2 do begin
        if TOTAL(mask[nx-1,j,k]+mask[nx-1,j+1,k]) eq 2 then begin
           n1=n1+1
        endif
        if TOTAL(mask[nx-1,j,k]+mask[nx-1,j,k+1]) eq 2 then begin
           n1=n1+1
        endif
    end 
end

for k=0,nz-2 do begin
    for i=0,nx-2 do begin
        if TOTAL(mask[i,ny-1,k]+mask[i+1,ny-1,k]) eq 2 then begin
           n1=n1+1
        endif
        if TOTAL(mask[i,ny-1,k]+mask[i,ny-1,k+1]) eq 2 then begin
           n1=n1+1
        endif
    end
end

for j=0,ny-2 do begin
    for i=0,nx-2 do begin
        if TOTAL(mask[i,j,nz-1]+mask[i+1,j,nz-1]) eq 2 then begin
           n1=n1+1
        endif
        if TOTAL(mask[i,j,nz-1]+mask[i,j+1,nz-1]) eq 2 then begin
           n1=n1+1
        endif
    end
end

print,"n1=",n1

;find n2 - the total number of squares inside each structure
;find n2 - this one is for the squares in xy plane
for k=0,nz-1 do begin  ;this was 0,nz-2 but I think we can do nz-1
    for j=0,ny-2 do begin
        for i=0,nx-2 do begin
            if TOTAL(mask[i,j,k]+mask[i+1,j,k]+mask[i,j+1,k]+mask[i+1,j+1,k]) eq 4 then begin
               n2=n2+1
            endif
        end
    end
end

;this one is for squares in the yz plane
for k=0,nz-2 do begin
    for j=0,ny-2 do begin
        for i=0,nx-1 do begin  ;this was nx-2 but I think nx-1 is correct
            if TOTAL(mask[i,j,k]+mask[i,j+1,k]+mask[i,j,k+1]+mask[i,j+1,k+1]) eq 4 then begin
               n2=n2+1
            endif
        end
    end
end

;this one is for squares in the xz plane
for k=0,nz-2 do begin
    for j=0,ny-1 do begin  ;this was ny-2 but I think ny-1 is right
        for i=0,nx-2 do begin
            if TOTAL(mask[i,j,k]+mask[i+1,j,k]+mask[i,j,k+1]+mask[i+1,j,k+1]) eq 4 then begin
               n2=n2+1
            endif
        end
    end
end

print,"n2=",n2


;finally, find n3 - the total number of complete cubes inside a structure
for k=0,nz-2 do begin
    for j=0,ny-2 do begin
        for i=0,nx-2 do begin
            if TOTAL(mask[i,j,k]+mask[i+1,j,k]+mask[i,j+1,k]+mask[i,j,k+1]+mask[i+1,j+1,k]+mask[i,j+1,k+1]+mask[i+1,j,k+1]+mask[i+1,j+1,k+1]) eq 8 then begin
               n3=n3+1
            endif
        end
    end 
end

print,"n3=",n3


a=128.0D

;now use n0,n1,n2 and n3 to estimate the Minkowski functionals
V0=n3
print,"V0=",V0
V1=2.0D/(9.0D*a)*(n2-3.0D*n3)
print,"V1=",V1
V2=2.0D/(9.0D*a*a)*(n1-2.0D*n2+3.0D*n3)
print,"V2=",V2
V3=1.0D/(a*a*a)*(n0-n1+n2-n3)
print,"V3=",V3

;now use the Minkowski functionals to produce the first set of shapefinders
T=V0/(2.0D*V1)  ;thickness
W=2.0D*V1/(V2*pi)  ;width
L=3.0D*V2/(4.0D*V3)  ;length

print,"Thickness=",T
print,"Width=",W
print,"Length=",L

;Finally, calculate the S=(P,F) values
Planarity=(W-T)/(W+T)
Filamentarity=(L-W)/(L+W)

print,""
print,""

print,"Planarity=",Planarity
print,"Filamentarity=",Filamentarity


end
