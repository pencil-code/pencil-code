;$Id$
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=16,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  reads current snapshot and dispays 9 different cross-sections
;  mv idl.ps fig/pn.ps
;
pc_read_var,obj=obj,variable=['bb'],/magic,/trimall
;
l=0
m=0
n=0
!p.multi=[0,3,3]
lev=grange(-1.,1.,21)*.5
for j=0,2 do contour,reform(obj.bb(*,*,n,j)),obj.x,obj.y,/fil,lev=lev,xtit='x',ytit='y'
for j=0,2 do contour,reform(obj.bb(*,m,*,j)),obj.x,obj.z,/fil,lev=lev,xtit='x',ytit='z'
for j=0,2 do contour,reform(obj.bb(l,*,*,j)),obj.y,obj.z,/fil,lev=lev,xtit='y',ytit='z'
!p.multi=0
END
