; $Id$
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
pc_read_var,obj=data,variables=['rho'],/magic,/trimall
;
plot,data.x,data.rho,yst=3,ps=8,back=255,col=1
;
print,'import riemann_bfield.gif'
END
