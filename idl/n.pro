file='tmp/n.dat'
a=rtable(file,9)
t=reform(a(1,*))
b2m=reform(a(5,*)) & bmax2=reform(a(6,*))
j2m=reform(a(7,*)) & jmax2=reform(a(8,*))
plot,t,1./jmax2,yst=0
save,file='n.sav',t,jmax2,j2m,bmax2,b2m
END
