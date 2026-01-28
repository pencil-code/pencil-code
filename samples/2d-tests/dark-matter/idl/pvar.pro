pc_read_var,obj=var,ivar=ivar,/trimall
pc_read_pvar,obj=pvar,ivar=ivar
;
contour,alog(var.np),/fil,nlev=30,/iso,var.x,var.y
oplot,pvar.xx[*,0],pvar.xx[*,1],ps=3 ;,/iso
;
END
