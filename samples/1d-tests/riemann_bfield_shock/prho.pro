; $Id: prho.pro,v 1.1 2006-08-23 15:13:16 theine Exp $
;
;  expansion shock, to be compared with Fig 1 of S.A.E.G. Falle (2002)
;  ApJL 577, L123-L126 "Rarefaction shocks, shock errors, and low
;  order accuracy in ZEUS"
;
!p.multi=0
!p.charsize=1.5
!p.title='t=30'
!x.title='x'
!y.title='!7q!6'
circ_sym,.4,0
;
;  in Falle's paper shock was initially at x=?
;
x0=0
;
;  exclude ghost zones
;
xxx=x(l1:l2)
rho=exp(lnrho(l1:l2,m1:m2,n1:n2))
;
plot,xxx,rho,yst=3,ps=8,back=255,col=1
;
print,'import riemann_bfield.gif'
END
