;$Id$
;
;  plots result of test field case
;
;  mv idl.ps ~/tex/mhd/alpha_pot/fig/pfield_3d.ps
;
pc_potentialfield_exterior_z,aa,aaa=aaa,zzz=zzz
;
bbb=curl(aaa)
jjj=curl2(aaa)
;
;  ghost zones
;
nn1=n1
nn2=n_elements(zzz)-4
;
;  interfaces
;
n=n_elements(zzz)
nnn1=(n-nz)/2
nnn2=(n+nz)/2-1
;
!p.multi=[0,1,2]
!p.charsize=1.2
;
j=2
!x.title='!6'
!y.title='!8z!6'
!p.title='!8B!dz!n!6'
contour,reform(bbb(11,m1:m2,nn1:nn2,j)),nlev=22,x(m1:m2),zzz(nn1:nn2)
oplot,x,x-x+.5,col=122
oplot,x,x-x-.5,col=122
;
j=0
!x.title='!8y!6'
!y.title='!8z!6'
!p.title='!8J!dx!n!6'
contour,reform(jjj(11,m1:m2,nn1:nn2,j)),nlev=22,x(m1:m2),zzz(nn1:nn2)
oplot,x,x-x+.5,col=122
oplot,x,x-x-.5,col=122
;
END
