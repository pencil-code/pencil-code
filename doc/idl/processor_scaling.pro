Nx=1.
Ny=256.
Nz=256.
;
nproc=2.^indgen(11)
!p.multi=[0,2,1]
;
nghost2=6.
nprocy=4. &  nprocz=nproc/nprocy
Work=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)
;
nprocy=8. &  nprocz=nproc/nprocy
Work8=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)
;
nprocy=16. &  nprocz=nproc/nprocy
Work16=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)
;
plot_oo,nproc,Nx*Ny*Nz/nproc,xst=3,yst=3
oplot,nproc,Work,li=1
oplot,nproc,Work8,li=2
oplot,nproc,Work16,li=3
;
nghost2=2.
nprocy=4. &  nprocz=nproc/nprocy
Work=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)
;
nprocy=8. &  nprocz=nproc/nprocy
Work8=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)
;
nprocy=16. &  nprocz=nproc/nprocy
Work16=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)
;
plot_oo,nproc,Nx*Ny*Nz/nproc,xst=3,yst=3
oplot,nproc,Work,li=1
oplot,nproc,Work8,li=2
oplot,nproc,Work16,li=3
;
END
