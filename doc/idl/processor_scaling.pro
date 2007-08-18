if !d.name eq 'PS' then begin
  device,xsize=18,ysize=9,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  mv idl.ps ../figs/processor_scaling.ps
;  convert ../figs/processor_scaling.ps ~/shared/My\ Pictures/Figures/Pencil/processor_scaling.png
;
Nx=1.
Ny=256.
Nz=256.
Nw=Nx*Ny*Nz
;
nproc=2.^indgen(11)
!p.multi=[0,2,1]
!p.charsize=1.6
!x.title='!8N!d!6proc!n'
!y.title='!6Work per processor'
;
nghost2=6.
nprocy=4. &  nprocz=nproc/nprocy
Work4=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
nprocy=8. &  nprocz=nproc/nprocy
Work8=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
nprocy=16. &  nprocz=nproc/nprocy
Work16=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
!p.title='!66th order'
plot_oo,nproc,1./nproc,xst=3,yst=3,li=1
oplot,nproc,Work4,li=2
oplot,nproc,Work8,li=3
oplot,nproc,Work16,li=0
;
xx=10. & dx=30. & siz=.8
legend,xx,dx,10^(-.18),2,'Nprocy=4',siz=siz
legend,xx,dx,10^(-.38),3,'Nprocy=8',siz=siz
legend,xx,dx,10^(-.58),0,'Nprocy=16',siz=siz
;
nghost2=2.
nprocy=4. &  nprocz=nproc/nprocy
Work4=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
nprocy=8. &  nprocz=nproc/nprocy
Work8=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
nprocy=16. &  nprocz=nproc/nprocy
Work16=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
!p.title='!62nd order'
plot_oo,nproc,1./nproc,xst=3,yst=3,li=1
oplot,nproc,Work4,li=2
oplot,nproc,Work8,li=3
oplot,nproc,Work16,li=0
;
END
