      program todx_magn
! Converts var.dat snapshot to the Dx readable binary data files.
! The data are read in again by .dx routines on regular rectangular
! mesh. The todx_magn.f should be copied to tmp subdirectory of the 
! directory of the run. Compile f77 todx_magn.f and run it in /tmp
! .net files should be linked to tmp directory (from pencil_modular/dx).
! Then, Data Explorer is ready to be started from tmp directory.
      parameter(mx=38,my=70,mz=38,mvar=8)
      parameter(nghost=3)
      dimension f(mx,my,mz,mvar)
      dimension x(mx),y(my),z(mz)
      dimension bbb(mx,my,mz,3)

! Define here the position of the corner of the box
! todx.f is better to be integrated in the main code (with a flag
! of wether one wants to create Dx output).
      x0=-1.
      y0=-4.
      z0=-1.

! Read snapshot file
      open(unit=17, file='./proc0/var.dat', form='unformatted')
      read(17) f
      read(17) t,x,y,z,dx,dy,dz
      close(17)
      print*,'t=',t,' dx=',dx,' dy=',dy,' dz=',dz
      
! Write out velocity field for Dx
      open(unit=18, file='uuu_dx.bin', form='unformatted')
      write(18) ((((real(f(i,j,k,n)),n=1,3),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(18)

! Write out ln density field for Dx
      open(unit=19, file='lnrho_dx.bin', form='unformatted')
      write(19) (((real(f(i,j,k,4)),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(19)

! Write out entropy field for Dx
      open(unit=20, file='ss_dx.bin', form='unformatted')
      write(20) (((real(f(i,j,k,5)),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(20)

! Write out vector potential for Dx
      open(unit=18, file='aaa_dx.bin', form='unformatted')
      write(18) ((((real(f(i,j,k,n)),n=6,8),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(18)

! Calculate magnetic field, using central differences.
      do i=2,mx-1
         do j=2,my-1
            do k=2,mz-1
               bbb(i,j,k,1)=(f(i,j+1,k,8)-f(i,j-1,k,8))/2./dy-
     > (f(i,j,k+1,7)-f(i,j,k-1,7))/2./dz
               bbb(i,j,k,2)=(f(i,j,k+1,6)-f(i,j,k-1,6))/2./dz-
     > (f(i+1,j,k,8)-f(i-1,j,k,8))/2./dx
               bbb(i,j,k,3)=(f(i+1,j,k,7)-f(i-1,j,k,7))/2./dx-
     > (f(i,j+1,k,6)-f(i,j-1,k,6))/2./dy               
            enddo
         enddo
      enddo   

! Write out magnetic field for Dx
      open(unit=18, file='bbb_dx.bin', form='unformatted')
      write(18) ((((real(bbb(i,j,k,n)),n=1,3),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(18)


! Create .dx scripts
      open(unit=21, file='uuu.dx')
      write(21,*) 'object 1 class gridpositions counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(21,*) 'origin  ', x0, y0, z0 
      write(21,*) 'delta  ', dx,'  0   ', '   0    '
      write(21,*) 'delta  ', '  0  ',dy,'  0  '
      write(21,*) 'delta','  0  ','  0  ',dz
      write(21,*) '  '
      write(21,*) 'object 2 class gridconnections counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(21,*) '  '
      write(21,*) 'object 3 class array type float rank 1 shape 3'
      write(21,*) 'items ', (mx-2*nghost)*(my-2*nghost)*(mz-2*nghost),
     > ' lsb binary'
      write(21,*) 'data file "uuu_dx.bin", 4'
      write(21,*) '  '
      write(21,*) 'attribute "dep" string "positions" '
      write(21,*) '  '
      write(21,*) 'object "regular positions regular connections" 
     > class field'
      write(21,*) 'component "positions" value 1'
      write(21,*) 'component "connections" value 2'
      write(21,*) 'component "data" value 3'
      write(21,*) 'end'
      close(21)

      open(unit=24, file='aaa.dx')
      write(24,*) 'object 1 class gridpositions counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(24,*) 'origin  ', x0, y0, z0 
      write(24,*) 'delta  ', dx,'  0   ', '   0    '
      write(24,*) 'delta  ', '  0  ',dy,'  0  '
      write(24,*) 'delta','  0  ','  0  ',dz
      write(24,*) '  '
      write(24,*) 'object 2 class gridconnections counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(24,*) '  '
      write(24,*) 'object 3 class array type float rank 1 shape 3'
      write(24,*) 'items ', (mx-2*nghost)*(my-2*nghost)*(mz-2*nghost),
     > ' lsb binary'
      write(24,*) 'data file "aaa_dx.bin", 4'
      write(24,*) '  '
      write(24,*) 'attribute "dep" string "positions" '
      write(24,*) '  '
      write(24,*) 'object "regular positions regular connections" 
     > class field'
      write(24,*) 'component "positions" value 1'
      write(24,*) 'component "connections" value 2'
      write(24,*) 'component "data" value 3'
      write(24,*) 'end'
      close(24)

      open(unit=25, file='bbb.dx')
      write(25,*) 'object 1 class gridpositions counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(25,*) 'origin  ', x0, y0, z0 
      write(25,*) 'delta  ', dx,'  0   ', '   0    '
      write(25,*) 'delta  ', '  0  ',dy,'  0  '
      write(25,*) 'delta','  0  ','  0  ',dz
      write(25,*) '  '
      write(25,*) 'object 2 class gridconnections counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(25,*) '  '
      write(25,*) 'object 3 class array type float rank 1 shape 3'
      write(25,*) 'items ', (mx-2*nghost)*(my-2*nghost)*(mz-2*nghost),
     > ' lsb binary'
      write(25,*) 'data file "bbb_dx.bin", 4'
      write(25,*) '  '
      write(25,*) 'attribute "dep" string "positions" '
      write(25,*) '  '
      write(25,*) 'object "regular positions regular connections" 
     > class field'
      write(25,*) 'component "positions" value 1'
      write(25,*) 'component "connections" value 2'
      write(25,*) 'component "data" value 3'
      write(25,*) 'end'
      close(25)

      open(unit=22, file='lnrho.dx')
      write(22,*) 'object 1 class gridpositions counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(22,*) 'origin  ', x0, y0, z0 
      write(22,*) 'delta  ', dx,'  0   ', '   0    '
      write(22,*) 'delta  ', '  0  ',dy,'  0  '
      write(22,*) 'delta','  0  ','  0  ',dz
      write(22,*) '  '
      write(22,*) 'object 2 class gridconnections counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(22,*) '  '
      write(22,*) 'object 3 class array type float rank 0'
      write(22,*) 'items ', (mx-2*nghost)*(my-2*nghost)*(mz-2*nghost),
     > ' lsb binary'
      write(22,*) 'data file "lnrho_dx.bin", 4'
      write(22,*) '  '
      write(22,*) 'attribute "dep" string "positions" '
      write(22,*) '  '
      write(22,*) 'object "regular positions regular connections" 
     > class field'
      write(22,*) 'component "positions" value 1'
      write(22,*) 'component "connections" value 2'
      write(22,*) 'component "data" value 3'
      write(22,*) 'end'
      close(22)

      open(unit=23, file='ss.dx')
      write(23,*) 'object 1 class gridpositions counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(23,*) 'origin  ', x0, y0, z0 
      write(23,*) 'delta  ', dx,'  0   ', '   0    '
      write(23,*) 'delta  ', '  0  ',dy,'  0  '
      write(23,*) 'delta','  0  ','  0  ',dz
      write(23,*) '  '
      write(23,*) 'object 2 class gridconnections counts ', mx-2*nghost,
     > my-2*nghost, mz-2*nghost
      write(23,*) '  '
      write(23,*) 'object 3 class array type float rank 0'
      write(23,*) 'items ', (mx-2*nghost)*(my-2*nghost)*(mz-2*nghost),
     > ' lsb binary'
      write(23,*) 'data file "ss_dx.bin", 4'
      write(23,*) '  '
      write(23,*) 'attribute "dep" string "positions" '
      write(23,*) '  '
      write(23,*) 'object "regular positions regular connections" 
     > class field'
      write(23,*) 'component "positions" value 1'
      write(23,*) 'component "connections" value 2'
      write(23,*) 'component "data" value 3'
      write(23,*) 'end'
      close(23)


      stop
      end
