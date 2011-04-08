pc_read_var,obj=var,/bb,/trimall
;
m=0
bbb=reform(var.bb(*,m,*,1))
contour,bbb,var.x,var.z,/fill,nlev=30,/iso
contour,reform(var.aa(*,m,*,1)),var.x,var.z,/fill,nlev=30,/iso
;
END
