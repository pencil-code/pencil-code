; $Id$
;
;  expansion shock, to be compared with Fig 1 of S.A.E.G. Falle (2002)
;  ApJL 577, L123-L126 "Rarefaction shocks, shock errors, and low
;  order accuracy in ZEUS"
;
!p.multi=0
!p.charsize=1.5
!p.title='t=100'
!x.title='x+1000'
!y.title='u'
circ_sym,.4,0
;
;  in Falle's paper shock was initially at x=700
;
x0=1000
;
;  exclude ghost zones
;
xxx=x(l1:l2)
uuu=uu(l1:l2,m1:m2,n1:n2,*)
;
plot,xxx+x0,uuu,xr=[250,450],yst=0,ps=8,back=255,col=1
;
;  compare with analytic solution
;
c2_fake=.53
vA2=dot2(bbb)/rho
cA2=vA2+cs2+c2_fake
vA_left=sqrt(vA2(0)) & vA_right=sqrt(vA2(nx-1))
cs_left=sqrt(cs2(0)) & cs_right=sqrt(cs2(nx-1))
cA_left=sqrt(cA2(0)) & cA_right=sqrt(cA2(nx-1))
uu_left=uuu(0) & uu_right=uuu(nx-1)
print,'cA_left,cA_right=',cA_left,cA_right
print,'cs_left,cs_right=',cs_left,cs_right
print,'uu_left,uu_right=',uu_left,uu_right
u1=-cA_left+uu_left
u2=-cA_right+uu_right
x1=u1*t
x2=u2*t
print,x1,x2
print,u1,u2
;
oplot,x0+[xxx(0),x1,x2,xxx(nx-1)],[uu_left,uu_left,uu_right,uu_right],col=122,thick=2
;oplot,xxx+x0,uuu,ps=8,col=1
;
print,'import expans_bfield.gif'
print,'import expans_bfield.jpg'
END
