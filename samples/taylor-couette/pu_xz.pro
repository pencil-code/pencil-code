;  $Id$
;
;  one should see 2 nice pairs of Taylor vortices after t=50
;
m=ny/2+2
print,y(m)
xxx=x(l1:l2)
zzz=z(n1:n2)
ux=means(reform(uu(l1:l2,m:m+1,n1:n2,0)),2)
uy=means(reform(uu(l1:l2,m:m+1,n1:n2,1)),2)
uz=means(reform(uu(l1:l2,m:m+1,n1:n2,2)),2)
contour,uy,xxx,zzz,/fil,nlev=30
vel_a,ux,uz,xxx,zzz,nvec=1000,/over,len=.02
;
oplot,-par2.r_int*[1,1],[par.xyz0(0),par.xyz1(2)],thick=5
oplot,+par2.r_int*[1,1],[par.xyz0(0),par.xyz1(2)],thick=5
oplot,-par2.r_int*[1,1],[par.xyz0(0),par.xyz1(2)],thick=4,col=0,li=2
oplot,+par2.r_int*[1,1],[par.xyz0(0),par.xyz1(2)],thick=4,col=0,li=2
;
oplot,-par2.r_ext*[1,1],[par.xyz0(0),par.xyz1(2)],thick=5
oplot,+par2.r_ext*[1,1],[par.xyz0(0),par.xyz1(2)],thick=5
oplot,-par2.r_ext*[1,1],[par.xyz0(0),par.xyz1(2)],thick=4,col=0,li=2
oplot,+par2.r_ext*[1,1],[par.xyz0(0),par.xyz1(2)],thick=4,col=0,li=2
;
END
