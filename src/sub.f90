module Sub 

  use Io

  implicit none

  interface poly              ! Overload the `poly' function
    module procedure poly_1
    module procedure poly_3
  endinterface

  contains

!***********************************************************************
    subroutine max_mn(a,res)
!
!  successively calculate maximum of a, where a is supplied
!  at each call. This routine initializes counter when m=n=1.
!
!   1-apr-01/axel+wolf: coded
!
      use Cdata
!
      real, dimension (nx) :: a
      real :: res
!
!      if (m==m1.and.n==n1) then
      if (lfirstpoint) then
        res=maxval(a)
      else
        res=amax1(res,maxval(a))
      endif
!
    endsubroutine max_mn
!***********************************************************************
    subroutine mean_mn(a,res)
!
!  successively calculate mean of a, where a is supplied
!  at each call. This routine initializes counter when m=n=1.
!
!   17-dec-01/wolf: coded
!
      use Cdata
!
      real, dimension (nx) :: a
      real :: res
!
!      if (m==m1.and.n==n1) then
      if (lfirstpoint) then
        res=sum(a*1.D0)         ! sum at double precision to improve accuracy
      else
        res=res+sum(a*1.D0)
      endif
!
    endsubroutine mean_mn
!***********************************************************************
    subroutine rms_mn(a,res)
!
!  successively calculate rms of a, where a is supplied
!  at each call. This routine initializes counter when m=n=1.
!
!   1-apr-01/axel+wolf: coded
!
      use Cdata
!
      real, dimension (nx) :: a
      real :: res
!
!      if (m==m1.and.n==n1) then
      if (lfirstpoint) then
        res=sum(a**2)
      else
        res=res+sum(a**2)
      endif
!
    endsubroutine rms_mn
!***********************************************************************
    subroutine rms2_mn(a2,res)
!
!  successively calculate rms of a, where a2=a^2 is supplied
!  at each call. This routine initializes counter when m=n=1.
!
!   1-apr-01/axel+wolf: coded
!
      use Cdata
!
      real, dimension (nx) :: a2
      real :: res
!
!      if (m==m1.and.n==n1) then
      if (lfirstpoint) then
        res=sum(a2)
      else
        res=res+sum(a2)
      endif
!
    endsubroutine rms2_mn
!***********************************************************************
    subroutine exps(a,b)
!
!  dot product with itself, to calculate max and rms values of a vector
!  29-sep-97/axel: coded,
!
      use Cdata
!
      real, dimension (mx,my,mz) :: a,b
!
      b=exp(a)
!
    endsubroutine exps
!***********************************************************************
    subroutine dot2(a,b)
!
!  dot product with itself, to calculate max and rms values of a vector
!  29-sep-97/axel: coded,
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a
      real, dimension (nx) :: b
!
      b=a(l1:l2,m,n,1)**2+a(l1:l2,m,n,2)**2+a(l1:l2,m,n,3)**2
!
    endsubroutine dot2
!***********************************************************************
    subroutine dot_mn(a,b,c)
!
!  dot product
!   3-apr-01/axel+gitta: coded
!
      use Cdata
!
      real, dimension (nx,3) :: a,b
      real, dimension (nx) :: c
!
      c=a(:,1)*b(:,1)+a(:,2)*b(:,2)+a(:,3)*b(:,3)
!
    endsubroutine dot_mn
!***********************************************************************
    subroutine dot2_mn(a,b)
!
!  dot product with itself, to calculate max and rms values of a vector
!  29-sep-97/axel: coded
!   1-apr-01/axel: adapted for cache-efficient sub-array formulation
!
      use Cdata
!
      real, dimension (nx,3) :: a
      real, dimension (nx) :: b
!
      b=a(:,1)**2+a(:,2)**2+a(:,3)**2
!
    endsubroutine dot2_mn
!***********************************************************************
    subroutine trace_mn(a,b)
!
!  trace of a matrix
!   3-apr-01/axel+gitta: coded
!
      use Cdata
!
      real, dimension (nx,3,3) :: a
      real, dimension (nx) :: b
!
      b=a(:,1,1)+a(:,2,2)+a(:,3,3)
!
    endsubroutine trace_mn
!***********************************************************************
    subroutine multmv_mn(a,b,c)
!
!  matrix multiplied with vector, gives vector
!   3-apr-01/axel+gitta: coded
!
      use Cdata
!
      real, dimension (nx,3,3) :: a
      real, dimension (nx,3) :: b,c
      real, dimension (nx) :: tmp
      integer :: i,j
!
      do i=1,3
        j=1
        tmp=a(:,i,j)*b(:,j)
        do j=2,3
          tmp=tmp+a(:,i,j)*b(:,j)
        enddo
        c(:,i)=tmp
      enddo
!
    endsubroutine multmv_mn
!***********************************************************************
    subroutine dot2mu(a,b,c)
!
!  dot product with itself, to calculate max and rms values of a vector
!  29-sep-97/axel: coded,
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a
      real, dimension (mx,my,mz) :: b,c
!
      c=b*(a(:,:,:,1)**2+a(:,:,:,2)**2+a(:,:,:,3)**2)
!
    endsubroutine dot2mu
!***********************************************************************
    subroutine dot(a,b,c)
!
!  dot product, c=a.b
!  29-sep-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,b
      real, dimension (mx,my,mz) :: c
!
      c=a(:,:,:,1)*b(:,:,:,1)+a(:,:,:,2)*b(:,:,:,2)+a(:,:,:,3)*b(:,:,:,3)
!
    endsubroutine dot
!***********************************************************************
    subroutine dotneg(a,b,c)
!
!  negative dot product, c=-a.b
!  29-sep-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,b
      real, dimension (mx,my,mz) :: c
!
      c=-a(:,:,:,1)*b(:,:,:,1)-a(:,:,:,2)*b(:,:,:,2)-a(:,:,:,3)*b(:,:,:,3)
!
    endsubroutine dotneg
!***********************************************************************
    subroutine dotadd(a,b,c)
!
!  negative dot product, c=-a.b
!  29-sep-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,b
      real, dimension (mx,my,mz) :: c
!
      c=c+a(:,:,:,1)*b(:,:,:,1)+a(:,:,:,2)*b(:,:,:,2)+a(:,:,:,3)*b(:,:,:,3)
!
    endsubroutine dotadd
!***********************************************************************
    subroutine multsv(a,b,c)
!
!  multiply scalar with a vector
!  29-sep-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: b,c
      real, dimension (mx,my,mz) :: a
      integer :: j
!
      do j=1,3
        c(:,:,:,j)=a*b(:,:,:,j)
      enddo
!
    endsubroutine multsv
!***********************************************************************
    subroutine multsv_add(a,b,c,d)
!
!  multiply scalar with a vector and subtract from another vector
!  29-oct-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,c,d
      real, dimension (mx,my,mz) :: b
      integer :: j
!
      do j=1,3
        d(:,:,:,j)=a(:,:,:,j)+b*c(:,:,:,j)
      enddo
!
    endsubroutine multsv_add
!***********************************************************************
    subroutine multsv_sub(a,b,c,d)
!
!  multiply scalar with a vector and subtract from another vector
!  29-oct-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,c,d
      real, dimension (mx,my,mz) :: b
      integer :: j
!
      do j=1,3
        d(:,:,:,j)=a(:,:,:,j)-b*c(:,:,:,j)
      enddo
!
    endsubroutine multsv_sub
!***********************************************************************
    subroutine cross(a,b,c)
!
!  cross product, c = a x b
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,b,c
!
      c(:,:,:,1)=a(:,:,:,2)*b(:,:,:,3)-a(:,:,:,3)*b(:,:,:,2)
      c(:,:,:,2)=a(:,:,:,3)*b(:,:,:,1)-a(:,:,:,1)*b(:,:,:,3)
      c(:,:,:,3)=a(:,:,:,1)*b(:,:,:,2)-a(:,:,:,2)*b(:,:,:,1)
!
    endsubroutine cross
!***********************************************************************
    subroutine crossp(a,b,c)
!
!  cross product, c = a x b
!
      use Cdata
!
      real, dimension (nx,3) :: a,b,c
!
      c(:,1)=a(:,2)*b(:,3)-a(:,3)*b(:,2)
      c(:,2)=a(:,3)*b(:,1)-a(:,1)*b(:,3)
      c(:,3)=a(:,1)*b(:,2)-a(:,2)*b(:,1)
!
    endsubroutine crossp
!***********************************************************************
    subroutine gij(f,k,g)
!
!  calculate gradient of a vector, return matrix
!   3-apr-01/axel+gitta: coded
!
      use Cdata
      use Deriv
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,3,3) :: g
      real, dimension (nx) :: tmp
      integer :: i,j,k,k1
!
      k1=k-1
      do i=1,3
        do j=1,3
          call der(f,k1+i,tmp,j)
          g(:,i,j)=tmp
        enddo
      enddo
!
    endsubroutine gij
!***********************************************************************
    subroutine grad(f,k,g)
!
!  calculate gradient of a scalar, get vector
!  29-sep-97/axel: coded
!
      use Cdata
      use Deriv
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,3) :: g
      real, dimension (nx) :: tmp
      integer :: k
!
      call der(f,k,tmp,1); g(:,1)=tmp
      call der(f,k,tmp,2); g(:,2)=tmp
      call der(f,k,tmp,3); g(:,3)=tmp
!
    endsubroutine grad
!***********************************************************************
    subroutine curl(f,k,g)
!
!  calculate curl of a vector, get vector
!  12-sep-97/axel: coded
!  10-sep-01/axel: adapted for cache efficiency
!
      use Cdata
      use Deriv
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,3) :: g
      real, dimension (nx) :: tmp1,tmp2
      integer :: k,k1
!
      k1=k-1
!
      call der(f,k1+3,tmp1,2)
      call der(f,k1+2,tmp2,3)
      g(:,1)=tmp1-tmp2
!
      call der(f,k1+1,tmp1,3)
      call der(f,k1+3,tmp2,1)
      g(:,2)=tmp1-tmp2
!
      call der(f,k1+2,tmp1,1)
      call der(f,k1+1,tmp2,2)
      g(:,3)=tmp1-tmp2
!
    endsubroutine curl
!***********************************************************************
    subroutine del2(f,k,del2f)
!
!  calculate del2 of a scalar, get scalar
!  12-sep-97/axel: coded
!
      use Cdata
      use Deriv
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: del2f,dfdx,dfdy,dfdz
      integer :: k
!
      call der2(f,k,dfdx,1)
      call der2(f,k,dfdy,2)
      call der2(f,k,dfdz,3)
      del2f=dfdx+dfdy+dfdz
!
    endsubroutine del2
!***********************************************************************
    subroutine del2v(f,k,del2f)
!
!  calculate del2 of a vector, get vector
!  28-oct-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,3) :: del2f
      real, dimension (nx) :: tmp
      integer :: i,k,k1
!
!  do the del2 diffusion operator
!
      k1=k-1
      do i=1,3
        call del2(f,k1+i,tmp)
        del2f(:,i)=tmp
      enddo
!
    endsubroutine del2v
!***********************************************************************
    subroutine del2v_etc(f,k,del2,graddiv,curlcurl)
!
!  calculates a number of second derivative expressions of a vector
!  outputs a number of different vector fields
!  12-sep-01/axel: coded
!
      use Cdata
      use Deriv
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,3,3) :: fjji,fijj
      real, dimension (nx,3), optional :: del2,graddiv,curlcurl
      real, dimension (nx) :: tmp
      integer :: i,j,k,k1
!
!  do the del2 diffusion operator
!
      k1=k-1
      do i=1,3
      do j=1,3
        call der2 (f,k1+i,tmp,  j); fijj(:,i,j)=tmp  ! f_{i,jj}
        call derij(f,k1+j,tmp,j,i); fjji(:,i,j)=tmp  ! f_{j,ji}
      enddo
      enddo
!
!  the diagonal terms have not been set in derij; do this now
!
      do j=1,3
        fjji(:,j,j)=fijj(:,j,j)
      enddo
!
      if (present(del2)) then
        do i=1,3
          del2(:,i)=fijj(:,i,1)+fijj(:,i,2)+fijj(:,i,3)
        enddo
      endif
!
      if (present(graddiv)) then
        do i=1,3
          graddiv(:,i)=fjji(:,i,1)+fjji(:,i,2)+fjji(:,i,3)
        enddo
      endif
!
      if (present(curlcurl)) then
        curlcurl(:,1)=fjji(:,2,1)-fijj(:,1,2)+fjji(:,3,1)-fijj(:,1,3)
        curlcurl(:,2)=fjji(:,3,2)-fijj(:,2,3)+fjji(:,1,2)-fijj(:,2,1)
        curlcurl(:,3)=fjji(:,1,3)-fijj(:,3,1)+fjji(:,2,3)-fijj(:,3,2)
      endif
!
    endsubroutine del2v_etc
!***********************************************************************
!   subroutine del2v_graddiv(f,del2f,graddiv)
!
!  calculate del2 of a vector, get vector
!  calculate also graddiv of the same vector
!   3-apr-01/axel: coded
!
!     use Cdata
!
!     real, dimension (mx,my,mz,3) :: f
!     real, dimension (mx,my,mz) :: scr
!     real, dimension (mx,3) :: del2f,graddiv
!     real, dimension (mx) :: tmp
!     integer :: j
!
!  do the del2 diffusion operator
!
!     do i=1,3
!       s=0.
!       scr=f(:,:,:,i)
!       do j=1,3
!         call der2(scr,tmp,j)
!tst      if (i==j) graddiv(:,i,j)=tmp
!tst      s=s+tmp
!       enddo
!       del2f(:,j)=s
!     enddo

!     call der2(f,dfdx,1)
!     call der2(f,dfdy,2)
!     call der2(f,dfdz,3)
!     del2f=dfdx+dfdy+dfdz
!
!   endsubroutine del2v_graddiv
!***********************************************************************
    subroutine inpup(file,a,nn)
!
!  read particle snapshot file
!  11-apr-00/axel: adapted from input
!
      use Cdata
!
      integer :: nn
      real, dimension (nn) :: a
      character (LEN=*) :: file
!
      open(1,file=file,form='unformatted')
      read(1) a
      close(1)
    endsubroutine inpup
!***********************************************************************
    subroutine inpui(file,a,nn)
!
!  read particle snapshot file
!  11-apr-00/axel: adapted from input
!
      use Cdata
!
      integer :: nn
      integer, dimension (nn) :: a
      character (LEN=*) :: file
!
      open(1,file=file,form='unformatted')
      read(1) a
      close(1)
    endsubroutine inpui
!***********************************************************************
    subroutine inpuf(file,a,nn)
!
!  read formatted snapshot
!   5-aug-98/axel: coded
!
      use Cdata
!
      integer :: nn
      real, dimension (mx,my,mz,nn) :: a
      character (LEN=*) :: file
!
      open(1,file=file)
      read(1,10) a
      read(1,10) t,x,y,z
      close(1)
!10    format(1p8e10.3)
10    format(8e10.3)
    endsubroutine inpuf
!***********************************************************************
    subroutine outpup(file,a,nn)
!
!  write snapshot file, always write mesh and time, could add other things
!  11-apr-00/axel: adapted from output
!
      integer :: nn
      real, dimension (nn) :: a
      character (LEN=*) :: file
!
      open(1,file=file,form='unformatted')
      write(1) a
      close(1)
    endsubroutine outpup
!***********************************************************************
    subroutine outpui(file,a,nn)
!
!  write snapshot file, always write mesh and time, could add other things
!  11-apr-00/axel: adapted from output
!
      integer :: nn
      integer, dimension (nn) :: a
      character (LEN=*) :: file
!
      open(1,file=file,form='unformatted')
      write(1) a
      close(1)
    endsubroutine outpui
!***********************************************************************
    subroutine outpuf(file,a,nn)
!
!  write formatted snapshot, otherwise like output
!   5-aug-98/axel: coded
!
      use Cdata
!
      integer :: nn
      character (LEN=*) :: file
      real, dimension (mx,my,mz,nn) :: a
!
      open(1,file=file)
      write(1,10) a
      write(1,10) t,x,y,z
      close(1)
!10    format(1p8e10.3)
10    format(8e10.3)
    endsubroutine outpuf
!***********************************************************************
    subroutine wsnap(chsnap,a)
!
!  write snapshot file, labelled consecutively
!  30-sep-97/axel: coded
!
      use Cdata
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar) :: a
      character (LEN=4) :: ch
      character (LEN=9) :: file
      character (LEN=*) :: chsnap
      logical lsnap
      integer, save :: ifirst,nsnap
      real, save :: tsnap
!
!  Output snapshot in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      file='tsnap.dat'
!
!  at first call, need to initialize tsnap
!  tsnap calculated in out1, but only available to root processor
!
      if (ifirst==0) then
        call out1 (file,tsnap,nsnap,dsnap,t)
        ifirst=1
      endif
!
      call out2 (file,tsnap,nsnap,dsnap,t,lsnap,ch,.true.)
      if (lsnap) then
!  update ghost zones for var.dat (cheap, since done infrequently)
        call initiate_isendrcv_bdry(a)
        call finalise_isendrcv_bdry(a)
        call output(chsnap//ch,a,mvar)
      endif
!
    endsubroutine wsnap
!***********************************************************************
    subroutine wdim(file,myout,mzout)
!
!  write dimension to file
!
!   8-sep-01/axel: adapted to take myout,mzout
!
      use Cdata
!
      character (LEN=*) :: file
      integer, optional :: myout,mzout
      integer :: myout1,mzout1
!
!  determine whether myout=my (as on each processor)
!  or whether myout is different (eg when writing out full array)
!
      if(present(myout)) then
        myout1=myout
      else
        myout1=my
      endif
!
!  Do the same for mz.
!
      if(present(mzout)) then
        mzout1=mzout
      else
        mzout1=mz
      endif
!
      open(1,file=file)
      write(1,'(4i7)') mx,myout1,mzout1,mvar
!
!  check for double precision
!
      if (1..eq.1.+1e-10) then
        write(1,'(a)') 'S'
      else
        write(1,'(a)') 'D'
      endif
!
!  write number of ghost cells (could be different in x, y and z)
!
      write(1,*) 0, nghost, nghost
!
      close(1)
    endsubroutine wdim
!***********************************************************************
    subroutine wparam ()
!
!  Write startup parameters
!  21-jan-02/wolf: coded
!
      use Cdata, only: Lx,Ly,Lz,cs0,gamma,gamma1,ss0,grads0
      use Mpicomm
!
      if (lroot) then
        open(1,FILE='tmp/param.dat',FORM='unformatted')
        write(1) Lx,Ly,Lz
        write(1) cs0,gamma,gamma1 ! Write gamma1 here to ensure it is in sync
        write(1) ss0,grads0
      endif
!
    endsubroutine wparam
!***********************************************************************
    subroutine rparam ()
!
!  Read startup parameters
!  21-jan-02/wolf: coded
!
      use Cdata, only: Lx,Ly,Lz,cs0,gamma,gamma1,ss0,grads0
      use Mpicomm
!
!      if (lroot) then
        open(1,FILE='tmp/param.dat',FORM='unformatted')
        read(1) Lx,Ly,Lz
        read(1) cs0,gamma,gamma1
        read(1) ss0,grads0
!      endif
        print*, "Lx,Ly,Lz=", Lx,Ly,Lz
        print*, "cs0,gamma,gamma1=", cs0,gamma,gamma1
!
    endsubroutine rparam
!***********************************************************************
    subroutine out1 (file,tout,nout,dtout,t)
!
      use Mpicomm
!
!  read in output times from control file
!  30-sep-97/axel: coded
!  24-aug-99/axel: allow for logarithmic spacing
!   9-sep-01/axel: adapted for MPI
!
      character (LEN=*) :: file
      integer :: lun,nout
      real :: tout,dtout,ttt,tt,t
      integer, parameter :: nbcast_array=2
      real, dimension(nbcast_array) :: bcast_array
      logical exist
!
!  depending on whether or not file exists, we need to
!  either read or write tout and nout from or to the file
!
      if (lroot) then
        inquire(file=file,exist=exist)
        lun=1
        open(lun,file=file)
        if (exist) then
          read(lun,*) tout,nout
        else
!
!  special treatment when dtout is negative
!
          if (dtout.lt.0.) then
            tout=alog10(t)
          else
            tout=t
          endif
          nout=0
          write(lun,*) tout,nout
        endif
        close(lun)
        bcast_array(1)=tout
        bcast_array(2)=nout
      endif
!
!  broadcast tout and nout, encrypt into floating point array.
!
      call mpibcast_real(bcast_array,nbcast_array)
      tout=bcast_array(1)
      nout=bcast_array(2)
!
!  special treatment when tt is negative
!  this has to do with different integer arithmetic for negative numbers
!
      tt=tout
      if (tt.lt.0.) then
        ttt=tt-1.
      else
        ttt=tt
      endif
!
!  make sure the tout is a good time
!
      if (dtout.ne.0.) tout=ttt+abs(dtout)-amod(tt,abs(dtout))
!
    endsubroutine out1
!***********************************************************************
    subroutine out2 (file,tout,nout,dtout,t,lout,ch,lch)
!
      use General
!
!  check whether we have output, done by all processors
!  30-sep-97/axel: coded
!  24-aug-99/axel: allow for logarithmic spacing
!
      character (LEN=*) :: file
      character (LEN=4) :: ch
      logical lout,lch
      real :: t,tt,tout,dtout
      integer :: lun,nout
!
      if (dtout.lt.0.) then
        tt=alog10(t)
      else
        tt=t
      endif
!
      if (tt.ge.tout) then
        tout=tout+abs(dtout)
        nout=nout+1
        lout=.true.
!
!  write corresponding value of tout to file
!  to make sure we have it, in case the code craches
!  if the disk is full, however, we need to reset the values manually
!
        lun=1
        open(lun,file=file)
        write(lun,*) tout,nout
        close(lun)
      else
        lout=.false.
      endif
!
!  if lch=.false. we don't want to generate a running file number (eg in wvid)
!  if lch=.true. we do want to generate character from nout for file name
!
      if (lch) call chn(nout,ch)
!
    endsubroutine out2
!***********************************************************************
    subroutine debugs (a,label)
!
!  print variable for debug purposes
!  29-oct-97/axel: coded
!
      use Cdata
!
      character (LEN=*) :: label
      real, dimension (mx,my,mz) :: a
!
      if (ip.le.6) then
        print*,'DEBUG: ',label,', min/max=',minval(a),maxval(a)
      endif
!
    endsubroutine debugs
!***********************************************************************
    subroutine debugv (a,label)
!
!  print variable for debug purposes
!  29-oct-97/axel: coded
!
      use Cdata
!
      character (LEN=*) :: label
      real, dimension (mx,my,mz,3) :: a
      integer :: j
!
      if (ip.le.6) then
        do j=1,3
          print*,'DEBUG: ',label,', min/max=',minval(a),maxval(a),j
        enddo
      endif
!
    endsubroutine debugv
!***********************************************************************
    subroutine smooth_3d(ff,nn)
!
!  Smooth scalar vector field FF binomially N times, i.e. with the
!  binomial coefficients (2*N \above k)/2^{2*N}.
!  20-apr-99/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz) :: ff
      integer :: j,nn
!
      do j=1,3
        call smooth_1d(ff,j,nn)
      enddo
!
    endsubroutine smooth_3d
!***********************************************************************
    subroutine smooth_1d(ff,idir,nn)
!
!  Smooth scalar vector field FF binomially N times in direction IDIR.
!  20-apr-99/wolf: coded
!   1-sep-01/axel: adapted for case with ghost layers
!
      use Cdata
!
      real, dimension (mx,my,mz) :: ff,gg
      integer :: idir,i,nn
!
!  don't smooth in directions in which there is no extent
!
      if (idir.eq.1.and.mx.lt.3) return
      if (idir.eq.2.and.my.lt.3) return
      if (idir.eq.3.and.mz.lt.3) return
!
      do i=1,nn
        gg = ff
        select case (idir)
        case (1)                  ! x direction
          ff(2:mx-1,:,:) = (gg(1:mx-2,:,:) + 2*gg(2:mx-1,:,:) + gg(3:mx,:,:))/4.
        case (2)                  ! y direction
          ff(:,2:my-1,:) = (gg(:,1:my-2,:) + 2*gg(:,2:my-1,:) + gg(:,3:my,:))/4.
        case (3)                  ! z direction
          ff(:,:,2:mz-1) = (gg(:,:,1:mz-2) + 2*gg(:,:,2:mz-1) + gg(:,:,3:mz))/4.
        case default
          print*,'Bad call to smooth_1d, idir = ', idir, ' should be 1,2 or 3'
          STOP
        endselect
      enddo
!
    endsubroutine smooth_1d
!***********************************************************************
    subroutine nearmax(f,g)
!
!  extract nearest maxima
!  12-oct-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f,g
!
      g(1     ,:,:)=amax1(f(1     ,:,:),f(2     ,:,:))
      g(2:mx-1,:,:)=amax1(f(1:mx-2,:,:),f(2:mx-1,:,:),f(3:mx,:,:))
      g(  mx  ,:,:)=amax1(              f(  mx-1,:,:),f(  mx,:,:))
!
!  check for degeneracy
!
      if (my.gt.1) then
        f(:,1     ,:)=amax1(g(:,1     ,:),g(:,2     ,:))
        f(:,2:my-1,:)=amax1(g(:,1:my-2,:),g(:,2:my-1,:),g(:,3:my,:))
        f(:,  my  ,:)=amax1(              g(:,  my-1,:),g(:,  my,:))
      else
        f=g
      endif
!
!  check for degeneracy
!
      if (mz.gt.1) then
        g(:,:,1     )=amax1(f(:,:,1     ),f(:,:,2     ))
        g(:,:,2:mz-1)=amax1(f(:,:,1:mz-2),f(:,:,2:mz-1),f(:,:,3:mz))
        g(:,:,  mz  )=amax1(              f(:,:,  mz-1),f(:,:,  mz))
      else
        g=f
      endif
!
    endsubroutine nearmax
!***********************************************************************
    subroutine wmax(lun,f)
!
!  calculate th location of the first few maxima
!   6-jan-00/axel: coded
!
      use Cdata
!
      integer :: lun,l,imax,imax2
      integer, parameter :: nmax=10
      real, dimension (4,nmax) :: fmax
      real, dimension (mx,my,mz) :: f
!
      fmax=0
      do n=1,mz
      do m=1,my
      do l=1,mx
        !
        !  find out whether this f is larger than the smallest max so far
        !
        if (f(l,m,n).gt.fmax(1,1)) then
          !
          !  yes, ok, so now we need to sort it in
          !
          sort_f_in: do imax=nmax,1,-1
            if (f(l,m,n).gt.fmax(1,imax)) then
              !
              !  shift the rest downwards
              !
              do imax2=1,imax-1
                fmax(:,imax2)=fmax(:,imax2+1)
              enddo
              fmax(1,imax)=f(l,m,n)
              fmax(2,imax)=x(l)
              fmax(3,imax)=y(m)
              fmax(4,imax)=z(n)
              exit sort_f_in
!              goto 99
            endif
          enddo sort_f_in
        endif
!99      continue
      enddo
      enddo
      enddo
      write(lun,*) t,fmax
!
    endsubroutine wmax
!***********************************************************************
    subroutine cvs_id(rcsfile, revision, date)
!
!  print CVS Revision info in a compact, yet structured form
!  17-jan-02/wolf: coded
!

      character (len=*) :: rcsfile, revision, date
      integer :: rcsflen, revlen, datelen

      rcsflen=len(rcsfile)
      revlen =len(revision)
      datelen=len(date)
      write(*,'(A,A,T28," version ",A,T50," of ",A)') "CVS: ", &
           rcsfile(10:rcsflen-4), &
           revision(12:revlen-1), &
           date(8:datelen-1)
!
    endsubroutine cvs_id
!***********************************************************************
    function poly_1(coef, x)
!
!  Horner's scheme for polynomial evaluation.
!  Version for 1d array.
!  17-jan-02/wolf: coded 
!
      real, dimension(:) :: coef
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: poly_1
      integer :: Ncoef,Nx,i

      Ncoef = size(coef,1)
      Nx = size(x,1)

      poly_1 = coef(Ncoef)
      do i=Ncoef-1,1,-1
        poly_1 = poly_1*x+coef(i)
      enddo

    endfunction poly_1
!***********************************************************************
    function poly_3(coef, x)
!
!  Horner's scheme for polynomial evaluation.
!  Version for 3d array.
!  17-jan-02/wolf: coded 
!
      real, dimension(:) :: coef
      real, dimension(:,:,:) :: x
      real, dimension(size(x,1),size(x,2),size(x,3)) :: poly_3
      integer :: Ncoef,Nx,i

      Ncoef = size(coef,1)
      Nx = size(x,1)

      poly_3 = coef(Ncoef)
      do i=Ncoef-1,1,-1
        poly_3 = poly_3*x+coef(i)
      enddo

    endfunction poly_3
!***********************************************************************

endmodule Sub
