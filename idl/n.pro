;
;  here we read the rprint files
;  to generate an index catalogue of what is written
;
@tmp/hydro
@tmp/magnetic
print,'nname=',nname
;
file='tmp/n.dat'
a=rtable(file,nname)
if i_t ne 0 then t=reform(a(i_t-1,*))
if i_it ne 0 then it=reform(a(i_it-1,*))
if i_dt ne 0 then dt=reform(a(i_dt-1,*))
if i_u2m ne 0 then u2m=reform(a(i_u2m-1,*))
if i_um2 ne 0 then um2=reform(a(i_um2-1,*))
if i_b2m ne 0 then b2m=reform(a(i_b2m-1,*))
if i_bm2 ne 0 then bm2=reform(a(i_bm2-1,*))
;
plot,t,u2m,yst=0
;save,file='hydro.sav',t,jmax2,j2m,bmax2,b2m
;save,file='magnetic.sav',t,jmax2,j2m,bmax2,b2m
END
