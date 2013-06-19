;
;  mv idl.ps ../tex/fig/pvar.ps
;
pc_read_var,obj=var,/trimall,variables=['uu','rho','pp','bb'],/magic
velovect,var.bb(*,*,0),var.bb(*,*,1)
bperp=sqrt(var.bb(*,*,0)^2+var.bb(*,*,1)^2)
;plot,var.x,var.bb(*,0,0)
;oplot,var.x,var.uu(*,0,0),col=122,li=2
END
