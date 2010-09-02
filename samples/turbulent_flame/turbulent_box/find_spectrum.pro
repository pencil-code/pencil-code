;
; Read in data
;
pc_read_var,obj=obj
pc_read_dim,obj=dim
 power3d_v,obj.uu(dim.l1:dim.l2,dim.m1:dim.m2,dim.n1:dim.n2,*)

END
