if !d.name eq 'PS' then begin
  device,xsize=18,ysize=14,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  plots the work done in the mesh and the ghost zones as a
;  function of total processor number.
;
; mv idl.ps ../figs/processor_scaling_constmem.ps
; convert !$ ~/shared/My\ Pictures/Figures/Pencil/processor_scaling_constmem.png
;
nproc=2.^indgen(22)
;!p.multi=[0,2,1]
!p.multi=0
!p.charsize=1.8
!x.title='!8N!d!6proc!n'
!y.title='!6Work per processor'
MB=1e6
;
nghost2=6.
nprocz=sqrt(nproc) & nprocy=nprocz
Nw=30.*MB*nproc
Nx=Nw^.3333 & Ny=Nx & Nz=Nx
Work=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
!p.title='!66th order'
plot_oo,nproc,1./nproc,xst=3,yst=3,li=1
oplot,nproc,Work,col=122
stop
;
nghost2=2.
nprocz=sqrt(nproc) & nprocy=nprocz
Nw=30.*MB*nproc
Nx=Nw^.3333 & Ny=Nx & Nz=Nx
Work=Nx*(Ny/nprocy+nghost2)*(Nz/nprocz+nghost2)/Nw
;
!p.title='!62nd order'
plot_oo,nproc,1./nproc,xst=3,yst=3,li=1
oplot,nproc,Work,col=122
;
END
