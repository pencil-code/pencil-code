pc_read_var,obj=var,/bb,/trimall
;
m=0
contour,reform(var.aa(*,m,*,1)),var.x,var.z,/fill,nlev=30,/iso
;
END
