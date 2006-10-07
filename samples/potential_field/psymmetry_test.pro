pc_potentialfield_exterior_z,aa,aaa=aaa,zzz=zzz
;
arev=aa
for n=n1,n2 do arev(*,*,n,0)=+aa(*,*,n1+n2-n,0)
for n=n1,n2 do arev(*,*,n,1)=+aa(*,*,n1+n2-n,1)
for n=n1,n2 do arev(*,*,n,2)=-aa(*,*,n1+n2-n,2)
pc_potentialfield_exterior_z,arev,aaa=aaarev,zzz=zzz
;
nn1=n1
nn2=n_elements(zzz)-4
;
bbb=curl(aaa)
bbbrev=curl(aaarev)
;
jjj=curl2(aaa)
jjjrev=curl2(aaarev)
;
!p.multi=[0,2,2]
;
j=2
contour,reform(jjj(11,m1:m2,nn1:nn2,j)),nlev=22,x(m1:m2),zzz(nn1:nn2)
contour,reform(jjjrev(11,m1:m2,nn1:nn2,j)),nlev=22,x(m1:m2),zzz(nn1:nn2)
;
j=0
contour,reform(bbb(11,m1:m2,nn1:nn2,j)),nlev=22,x(m1:m2),zzz(nn1:nn2)
contour,reform(bbbrev(11,m1:m2,nn1:nn2,j)),nlev=22,x(m1:m2),zzz(nn1:nn2)
;
END
