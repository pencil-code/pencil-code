;
;  here we read the rprint files
;  to generate an index catalogue of what is written
;
@tmp/general
@tmp/hydro
@tmp/density
@tmp/entropy
@tmp/magnetic
print,'nname=',nname
;
filen='tmp/n.dat'
a=rtable(filen,nname)
if defined(i_t) ne 0 then t=reform(a(i_t-1,*))
if defined(i_it) ne 0 then it=reform(a(i_it-1,*))
if defined(i_dt) ne 0 then dt=reform(a(i_dt-1,*))
if defined(i_u2m) ne 0 then u2m=reform(a(i_u2m-1,*))
if defined(i_um2) ne 0 then um2=reform(a(i_um2-1,*))
if defined(i_b2m) ne 0 then b2m=reform(a(i_b2m-1,*))
if defined(i_bm2) ne 0 then bm2=reform(a(i_bm2-1,*))
if defined(i_abm) ne 0 then abm=reform(a(i_abm-1,*))
if defined(i_jbm) ne 0 then jbm=reform(a(i_jbm-1,*))
if defined(i_ssm) ne 0 then ssm=reform(a(i_ssm-1,*))
if defined(i_eth) ne 0 then eth=reform(a(i_eth-1,*))
if defined(i_ekin) ne 0 then ekin=reform(a(i_ekin-1,*))
if defined(i_rhom) ne 0 then rhom=reform(a(i_rhom-1,*))
if defined(i_bmz) ne 0 then bmz=reform(a(i_bmz-1,*))
;
;!p.multi=[0,1,2]
if i_um2 ne 0 then plot,t,u2m,yst=0
if i_bm2 ne 0 then oplot,t,b2m,col=122
!p.multi=0
;save,file='hydro.sav',t,jmax2,j2m,bmax2,b2m
;save,file='magnetic.sav',t,jmax2,j2m,bmax2,b2m
END
