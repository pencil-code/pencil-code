      program todx
! Converts var.dat snapshot to the Dx readable binary data files.
! The data are read in again by .dx routines on regular rectangular
! mesh. The todx.f should be linked in the 
! directory of the run. Compile f77 todx.f -o todx.x and run todx.x
! .net files should be linked to the run directory (from pencil-code/dx).
! Then, Data Explorer is ready to be started from run directory.
      parameter(mx=27,my=27,mz=47,mvar=5)
      parameter(nghost=3)
      dimension f(mx,my,mz,mvar)
      dimension x(mx),y(my),z(mz)

! Define here the position of the corner of the box
! todx.f is better to be integrated in the main code (with a flag
! of wether one wants to create Dx output).
      x0=-0.5
      y0=-0.5
      z0=-0.3

! Read snapshot file
      open(unit=17, file='./tmp/proc0/var.dat', form='unformatted')
      read(17) f
      read(17) t,x,y,z,dx,dy,dz
      close(17)
      print*,'t=',t,' dx=',dx,' dy=',dy,' dz=',dz
      
! Write out velocity field for Dx
      open(unit=18, file='./tmp/uuu_dx.bin', form='unformatted')
      write(18) ((((real(f(i,j,k,n)),n=1,3),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(18)

! Write out ln density field for Dx
      open(unit=19, file='./tmp/lnrho_dx.bin', form='unformatted')
      write(19) (((real(f(i,j,k,4)),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(19)

! Write out entropy field for Dx
      open(unit=20, file='./tmp/ss_dx.bin', form='unformatted')
      write(20) (((real(f(i,j,k,5)),k=nghost+1,mz-nghost),
     > j=nghost+1,my-nghost),i=nghost+1,mx-nghost)
      close(20)

! Create .dx scripts
      open(unit=21, file='./tmp/uuu.dx')
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

      open(unit=22, file='./tmp/lnrho.dx')
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

      open(unit=23, file='./tmp/ss.dx')
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
