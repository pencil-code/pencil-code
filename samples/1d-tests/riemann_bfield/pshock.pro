; $Id$
;
;  expansion shock, to be compared with Fig 1 of S.A.E.G. Falle (2002)
;  ApJL 577, L123-L126 "Rarefaction shocks, shock errors, and low
;  order accuracy in ZEUS"
;
exact=[[0.52,  0.4993], $
[0.5828,  0.4993], $
[0.5828,  0.5414], $
[0.7126,  0.5414], $
[0.7216,  0.5299], $
[0.7843,  0.5299], $
[0.7843,  0.4318], $
[0.8298,  0.4318], $
[0.8298,  0.1961], $
[0.9649,  0.1961], $
[0.9649,  0.1019], $
[1.020,  0.1019]]


!p.multi=0
!p.multi=[0,2,3,0,0]
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
pc_init
pc_read_var,obj=data,variables=['rho','tt','uu','pp','bb'],/magic,/trimall
;
plot,data.x,data.rho,yst=3,ps=8,back=255,col=1,ytitle='rho'
oplot,exact[0,*]*1E3-800,exact[1,*],col=122

plot,data.x,data.tt,yst=3,ps=8,back=255,col=1, ytitle='T'
plot,data.x,data.uu[*,0],yst=3,ps=8,back=255,col=1, ytitle='ux'
plot,data.x,data.pp,yst=3,ps=8,back=255,col=1, ytitle='P'

bb=data.bb
bb[*,0]=bb[*,0]+2.
plot,data.x,dot2(bb),yst=3,ps=8,back=255,col=1, ytitle='B!U2!6'
;
print,'import riemann_bfield.gif'
END
