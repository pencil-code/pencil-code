! $Id: sub.f90,v 1.84 2002-09-21 14:05:53 dobler Exp $ 

module Sub 

!  use Io

  implicit none

  interface poly                ! Overload the `poly' function
    module procedure poly_1
    module procedure poly_3
  endinterface

  interface notanumber          ! Overload the `notanumber' function
    module procedure notanumber_1
    module procedure notanumber_2
    module procedure notanumber_3
  endinterface

  interface cvs_id              ! Overload the cvs_id function
    module procedure cvs_id_1
    module procedure cvs_id_3
  endinterface


  contains

!***********************************************************************
    subroutine save_name(a,iname)
!
!  Lists the value of a (must be treated as real) in gname array
!
!  26-may-02/axel: adapted from max_mn_name
!
      use Cdata
!
      real :: a
      integer :: iname
!
!  Set corresponding entry in itype_name
!  This routine is to be called only once per step
!
      fname(iname)=a
      itype_name(iname)=ilabel_save
!
    endsubroutine save_name
!***********************************************************************
    subroutine max_mn_name(a,iname,lsqrt)
!
!  successively calculate maximum of a, where a is supplied
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!   4-may-02/axel: adapted for fname array
!  23-jun-02/axel: allows for taking square root in the end
!
      use Cdata
!
      real, dimension (nx) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      if (lfirstpoint) then
        fname(iname)=maxval(a)
      else
        fname(iname)=amax1(fname(iname),maxval(a))
      endif
!
!  set corresponding entry in itype_name
!
      if (present(lsqrt)) then
        itype_name(iname)=ilabel_max_sqrt
      else
        itype_name(iname)=ilabel_max
      endif
!
    endsubroutine max_mn_name
!***********************************************************************
    subroutine sum_mn_name(a,iname,lsqrt)
!
!  successively calculate sum of a, where a is supplied
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!   4-may-02/axel: adapted for fname array
!  23-jun-02/axel: allows for taking square root in the end
!
      use Cdata
!
      real, dimension (nx) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      if (lfirstpoint) then
        fname(iname)=sum(a)
      else
        fname(iname)=fname(iname)+sum(a)
      endif
!
!  set corresponding entry in itype_name
!
      if (present(lsqrt)) then
        itype_name(iname)=ilabel_sum_sqrt
      else
        itype_name(iname)=ilabel_sum
      endif
!
    endsubroutine sum_mn_name
!***********************************************************************
    subroutine zsum_mn_name(a,iname)
!
!  successively calculate sum of a, where a is supplied
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   5-jun-02/axel: adapted from sum_mn_name
!
      use Cdata
!
      real, dimension (nx) :: a
      integer :: iname,n_nghost
!
!  always initialize to zero, including other parts of the z-array
!  which later be merged with an mpi reduce command.
!  Thus, it must be zero if is not set otherwise.
!
      if (lfirstpoint) fnamez(:,:,iname)=0.
!
!  n starts with nghost+1=4, so the correct index is n-nghost
!
      n_nghost=n-nghost
      fnamez(n_nghost,ipz+1,iname)=fnamez(n_nghost,ipz+1,iname)+sum(a)
!
    endsubroutine zsum_mn_name
!***********************************************************************
    subroutine xysum_mn_name(a,iname)
!
!  successively calculate sum of a, where a is supplied
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!  19-jun-02/axel: adapted from zsum_mn_name
!
      use Cdata
!
      real, dimension (nx) :: a
      integer :: iname,m_nghost
!
!  always initialize to zero, including other parts of the z-array
!  which later be merged with an mpi reduce command.
!  Thus, it must be zero if is not set otherwise.
!
      if (lfirstpoint) fnamexy(:,:,:,iname)=0.
!
!  m starts with nghost+1=4, so the correct index is m-nghost
!  keep full x-dependence
!
      m_nghost=m-nghost
      fnamexy(:,m_nghost,ipy+1,iname)=fnamexy(:,m_nghost,ipy+1,iname)+a
!
    endsubroutine xysum_mn_name
!***********************************************************************
    subroutine max_mn(a,res)
!
!  successively calculate maximum of a, where a is supplied
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      use Cdata, only: nx,lfirstpoint
!
      real, dimension (nx) :: a
      real :: res
!
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
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   17-dec-01/wolf: coded
!
      use Cdata, only: nx,lfirstpoint
!
      real, dimension (nx) :: a
      real :: res
!
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
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      use Cdata, only: nx,lfirstpoint
!
      real, dimension (nx) :: a
      real :: res
!
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
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      use Cdata, only: nx,lfirstpoint
!
      real, dimension (nx) :: a2
      real :: res
!
      if (lfirstpoint) then
        res=sum(a2)
      else
        res=res+sum(a2)
      endif
!
    endsubroutine rms2_mn
!***********************************************************************
    subroutine sum_mn(a,res)
!
!  successively calculate the sum over all points of a, where a is supplied
!  at each call. This routine initializes counter when lfirstpoint=.true.
!
!   1-apr-01/axel+wolf: coded
!
      use Cdata
!
      real, dimension (nx) :: a
      real :: res
!
      if (lfirstpoint) then
        res=sum(a)
      else
        res=res+sum(a)
      endif
!
    endsubroutine sum_mn
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
      intent(in) :: a
      intent(out) :: b
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
      intent(in) :: a,b
      intent(out) :: c
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
      intent(in) :: a
      intent(out) :: b
!
      b=a(:,1)**2+a(:,2)**2+a(:,3)**2
!
    endsubroutine dot2_mn
!**********************************************************************
    subroutine div(f,k,g)
!
!  calculate divergence of vector, get scalar
!  13-dec-01/nils: coded
!  16-jul-02/nils: adapted from pencil_mpi
!
      use Cdata
      use Deriv
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: g, tmp
      integer :: k,k1
!
      k1=k-1
!
      call der(f,k1+1,tmp,1)
      g=tmp
      call der(f,k1+2,tmp,2)
      g=g+tmp
      call der(f,k1+3,tmp,3)
      g=g+tmp
!
    end subroutine div
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
      intent(in) :: a
      intent(out) :: b
!
      b=a(:,1,1)+a(:,2,2)+a(:,3,3)
!
    endsubroutine trace_mn
!***********************************************************************
    subroutine multvv_mat_mn(a,b,c)
!
!  vector multiplied with vector, gives matrix
!   21-dec-01/nils: coded
!   16-jul-02/nils: adapted from pencil_mpi
!
      use Cdata
!
      real, dimension (nx,3) :: a,b
      real, dimension (nx,3,3) :: c
      integer :: i,j
!
      do i=1,3
         do j=1,3
            c(:,i,j)=a(:,j)*b(:,i)
         end do
      end do
!
    end subroutine multvv_mat_mn
!***********************************************************************
    subroutine multmm_sc_mn(a,b,c)
!
!  matrix multiplied with matrix, gives scalar
!   21-dec-01/nils: coded
!   16-jul-02/nils: adapted from pencil_mpi
!
      use Cdata
!
      real, dimension (nx,3,3) :: a,b
      real, dimension (nx) :: c
      integer :: i,j
!
      c=0
      do i=1,3
         do j=1,3
            c=c+a(:,i,j)*b(:,i,j)
         end do
      end do
!
    end subroutine multmm_sc_mn
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
      intent(in) :: a,b
      intent(out) :: c
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
      intent(in) :: a,b
      intent(out) :: c
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
      intent(in) :: a,b
      intent(out) :: c
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
      intent(in) :: a,b
      intent(out) :: c
!
      c=-a(:,:,:,1)*b(:,:,:,1)-a(:,:,:,2)*b(:,:,:,2)-a(:,:,:,3)*b(:,:,:,3)
!
    endsubroutine dotneg
!***********************************************************************
    subroutine dotadd(a,b,c)
!
!  add dot product, c=c+a.b
!  29-sep-97/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,b
      real, dimension (mx,my,mz) :: c
!
      intent(in) :: a,b
      intent(out) :: c
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
      intent(in) :: a,b
      intent(out) :: c
!
      do j=1,3
        c(:,:,:,j)=a*b(:,:,:,j)
      enddo
!
    endsubroutine multsv
!***********************************************************************
    subroutine multsv_mn(a,b,c)
!
!  vector multiplied with scalar, gives vector
!   22-nov-01/nils erland: coded
!
      use Cdata
!
      intent(in) :: a,b
      intent(out) :: c
!
      real, dimension (nx,3) :: a, c
      real, dimension (nx) :: b
      integer :: i
!
      do i=1,3
        c(:,i)=a(:,i)*b(:)
      enddo
!
    endsubroutine multsv_mn
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
      intent(in) :: a,b,c
      intent(out) :: d
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
      intent(in) :: a,b,c
      intent(out) :: d
!
      do j=1,3
        d(:,:,:,j)=a(:,:,:,j)-b*c(:,:,:,j)
      enddo
!
    endsubroutine multsv_sub
!***********************************************************************
    subroutine multvs_mn(a,b,c)
!
!  vector multiplied with scalar, gives vector
!   22-nov-01/nils erland: coded
!
      use Cdata
!
      real, dimension (nx,3) :: a, c
      real, dimension (nx) :: b
      integer :: i
!
      do i=1,3
        c(:,i)=a(:,i)*b(:)
      enddo
!
    endsubroutine multvs_mn
!***********************************************************************
    subroutine cross(a,b,c)
!
!  cross product, c = a x b
!
      use Cdata
!
      real, dimension (mx,my,mz,3) :: a,b,c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c(:,:,:,1)=a(:,:,:,2)*b(:,:,:,3)-a(:,:,:,3)*b(:,:,:,2)
      c(:,:,:,2)=a(:,:,:,3)*b(:,:,:,1)-a(:,:,:,1)*b(:,:,:,3)
      c(:,:,:,3)=a(:,:,:,1)*b(:,:,:,2)-a(:,:,:,2)*b(:,:,:,1)
!
    endsubroutine cross
!***********************************************************************
    subroutine cross_mn(a,b,c)
!
!  cross product, c = a x b, for stencil variables.
!  Previously called crossp.
!
      use Cdata
!
      real, dimension (nx,3) :: a,b,c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c(:,1)=a(:,2)*b(:,3)-a(:,3)*b(:,2)
      c(:,2)=a(:,3)*b(:,1)-a(:,1)*b(:,3)
      c(:,3)=a(:,1)*b(:,2)-a(:,2)*b(:,1)
!
    endsubroutine cross_mn
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
      intent(in) :: f,k
      intent(out) :: g
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
      intent(in) :: f,k
      intent(out) :: g
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
      intent(in) :: f,k
      intent(out) :: g
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
      intent(in) :: f,k
      intent(out) :: del2f
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: del2f,d2fdx,d2fdy,d2fdz
      integer :: k
!
      call der2(f,k,d2fdx,1)
      call der2(f,k,d2fdy,2)
      call der2(f,k,d2fdz,3)
      del2f=d2fdx+d2fdy+d2fdz
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
      intent(in) :: f,k
      intent(out) :: del2f
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
!  outputs a number of different vector fields.
!  Surprisingly, calling derij only if graddiv or curlcurl are present
!  does not speed up the code on Mephisto @ 32x32x64.
!
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
      intent(in) :: f,k
      intent(out) :: del2,graddiv,curlcurl
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
!  ** They are automatically set above, because derij   **
!  ** doesn't overwrite the value of tmp for i=j!       **
!
!     do j=1,3
!       fjji(:,j,j)=fijj(:,j,j)
!     enddo
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
        curlcurl(:,1)=fjji(:,1,2)-fijj(:,1,2)+fjji(:,1,3)-fijj(:,1,3)
        curlcurl(:,2)=fjji(:,2,3)-fijj(:,2,3)+fjji(:,2,1)-fijj(:,2,1)
        curlcurl(:,3)=fjji(:,3,1)-fijj(:,3,1)+fjji(:,3,2)-fijj(:,3,2)
      endif
!
    endsubroutine del2v_etc
!***********************************************************************
    subroutine g2ij(f,k,g)
!
!  calculates all second derivative of a scalar
!
!  11-jul-02/axel: coded
!
      use Cdata
      use Deriv
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx,3,3) :: g
      real, dimension (nx) :: tmp
      integer :: i,j,k
!
      intent(in) :: f,k
      intent(out) :: g
!
!  run though all 9 possibilities, treat diagonals separately
!
      do j=1,3
        call der2 (f,k,tmp,j); g(:,j,j)=tmp
        do i=j+1,3
          call derij(f,k,tmp,i,j); g(:,i,j)=tmp; g(:,j,i)=tmp
        enddo
      enddo
!
    endsubroutine g2ij
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
    subroutine del6(f,k,del6f)
!
!  calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
!  than del2^3) of a scalar for hyperdiffusion
!  8-jul-02/wolf: coded
!
      use Cdata
      use Deriv
!
      intent(in) :: f,k
      intent(out) :: del6f
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: del6f,d6fdx,d6fdy,d6fdz
      integer :: k
!
      call der6(f,k,d6fdx,1)
      call der6(f,k,d6fdy,2)
      call der6(f,k,d6fdz,3)
      del6f = d6fdx + d6fdy + d6fdz
!
    endsubroutine del6
!***********************************************************************
    subroutine del6_nodx(f,k,del6f)
!
!  calculate something similar to del6, but ignoring the steps dx, dy, dz.
!  Useful for Nyquist filetering, where you just want to remove the
!  Nyquist frequency fully, while retaining the amplitude in small wave
!  numbers.
!  8-jul-02/wolf: coded
!
      use Cdata
      use Deriv
!
      intent(in) :: f,k
      intent(out) :: del6f
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: del6f,d6fdx,d6fdy,d6fdz
      integer :: k
!
      call der6(f,k,d6fdx,1,IGNOREDX=.true.)
      call der6(f,k,d6fdy,2,IGNOREDX=.true.)
      call der6(f,k,d6fdz,3,IGNOREDX=.true.)
      del6f = d6fdx + d6fdy + d6fdz
!
    endsubroutine del6_nodx
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
      character (len=*) :: file
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
      character (len=*) :: file
!
      open(1,file=file,form='formatted')
      read(1,*) a
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
      character (len=*) :: file
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
      character (len=*) :: file
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
      character (len=*) :: file
!
      open(1,file=file,form='formatted')
      write(1,*) a
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
      character (len=*) :: file
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
    subroutine wdim(file,mxout,myout,mzout)
!
!  write dimension to file
!
!   8-sep-01/axel: adapted to take myout,mzout
!
      use Cdata
      use Mpicomm, only: ipx,ipy,ipz
!
      character (len=*) :: file
      integer, optional :: mxout,myout,mzout
      integer :: mxout1,myout1,mzout1
!
!  determine whether mxout=mx (as on each processor)
!  or whether mxout is different (eg when writing out full array)
!
      if(present(mxout)) then; mxout1=mxout; else; mxout1=mx; endif
      if(present(myout)) then; myout1=myout; else; myout1=my; endif
      if(present(mzout)) then; mzout1=mzout; else; mzout1=mz; endif
!
      open(1,file=file)
      write(1,'(4i7)') mxout1,myout1,mzout1,mvar
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
      write(1,'(3i3)') nghost, nghost, nghost
      if (present(mzout)) then
        write(1,'(3i3)') nprocx, nprocy, nprocz
      else
        write(1,'(3i3)') ipx, ipy, ipz
      endif
!
      close(1)
    endsubroutine wdim
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
      character (len=*) :: file
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
      character (len=*) :: file
      character (len=4) :: ch
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
      character (len=*) :: label
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
      character (len=*) :: label
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
!  WARNING: This routine is likely to be broken if you use MPI
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
!  WARNING: This routine is likely to be broken if you use MPI
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
    subroutine cvs_id_1(cvsid)
!
!  print CVS Revision info in a compact, yet structured form
!  Expects the standard CVS Id: line as argument
!  25-jun-02/wolf: coded
!
      character (len=*) :: cvsid
      character (len=20) :: rcsfile, revision, author, date
      character (len=200) :: fmt
      character (len=20) :: tmp1,tmp2,tmp3,tmp4
      integer :: ir0,ir1,iv0,iv1,id0,id1,id2,ia0,ia1
      integer :: rw=18, vw=12, aw=10, dw=19 ! width of individual fields

      !
      !  rcs file name
      !
      ir0 = index(cvsid, ":") + 2
      ir1 = ir0 + index(cvsid(ir0+1:), ",") - 1
      rcsfile = cvsid(ir0:ir1)
      !
      !  version number
      !
      iv0 = ir1 + 4
      iv1 = iv0 + index(cvsid(iv0+1:), " ") - 1
      revision = cvsid(iv0:iv1)
      !
      !  date
      !
      id0 = iv1 + 2             ! first char of date
      id1 = iv1 + 12            ! postition of space
      id2 = iv1 + 20            ! last char of time
      date = cvsid(id0:id2)
      !
      !  author
      !
      ia0 = id2 + 2
      ia1 = ia0 + index(cvsid(ia0+1:), " ") - 1
      author = cvsid(ia0:ia1)
      !
      !  constuct format
      !
      write(tmp1,*) rw
      write(tmp2,*) 6+rw
      write(tmp3,*) 6+rw+4+vw
      write(tmp4,*) 6+rw+4+vw+2+aw
!      fmt = '(A, A' // trim(adjustl(tmp1)) &
      fmt = '(A, A' &
           // ', T' // trim(adjustl(tmp2)) &
           // ', " v. ", A, T' // trim(adjustl(tmp3)) &
           // ', " (", A, T' // trim(adjustl(tmp4)) &
           // ', ") ", A)'
      !
      !  write string
      !
      write(*,fmt) "CVS: ", &
           trim(rcsfile), &
           revision(1:vw), &
           author(1:aw), &
           date(1:dw)
      !write(*,'(A)') '123456789|123456789|123456789|123456789|123456789|12345'
      !write(*,'(A)') '         1         2         3         4         5'
!
    endsubroutine cvs_id_1
!***********************************************************************
    subroutine cvs_id_3(rcsfile, revision, date)
!
!  print CVS revision info in a compact, yet structured form
!  Old version: expects filename, version and date as three separate arguments
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
    endsubroutine cvs_id_3
!***********************************************************************
    subroutine identify_bcs(varname,idx)
!
!  print boundary conditions for scalar field
!
!  19-jul-02/wolf: coded
!
      use Cdata
!
      character (len=*) :: varname
      integer :: idx
!
      write(*,'(A,A6,",  x: <",A6,">, y: <",A6,">,  z: <",A6,">")') &
           'Bcs for ', varname, &
           trim(bcx(idx)), trim(bcy(idx)), trim(bcz(idx))
!
    endsubroutine identify_bcs
!***********************************************************************
    function noform(cname)
!
!  returns the name without format, fills empty space
!  of correct length (depending on format) with dashes
!
!  22-jun-02/axel: coded 
!
      character (len=*) :: cname
      character (len=20) :: noform,cform,cnumber,dash='----------'
      integer :: index_e,index_f,index_g,index_i,index_d,index_r,index1,index2
      integer :: iform0,iform1,iform2,length,number,number1,number2
!
      intent(in)  :: cname
!
!  find position of left bracket to isolate format, cform
!
      iform0=index(cname,' ')
      iform1=index(cname,'(')
      iform2=index(cname,')')
!
!  set format; use default if not given
!  Here we keep the parenthesis in cform
!
      if (iform1>0) then
        cform=cname(iform1:iform2)
        length=iform1-1
      else
        cform='(1p,e10.2,0p)'
        length=iform0-1
      endif
!
!  find length of formatted expression, examples: f10.2, e10.3, g12.1
!  index_1 is the position of the format type (f,e,g), and
!  index_d is the position of the dot
!
      index_e=scan(cform,'eE')
      index_f=scan(cform,'fF')
      index_g=scan(cform,'gG')
      index_i=scan(cform,'iI')
      index_d=index(cform,'.')
      index_r=index(cform,')')
      index1=max(index_e,index_f,index_g,index_i)
      index2=index_d; if(index_d==0) index2=index_r
!
!  calculate the length of the format and assemble expression for legend
!
      cnumber=cform(index1+1:index2-1)
      read(cnumber,'(i4)',err=99) number
10    number1=max(0,(number-length)/2)
      number2=max(0,number-length-number1)
      noform=dash(1:number1)//cname(1:length)//dash(1:number2)
      return
!
! in case of errors:
!
99    print*,'noform: formatting problem'
      print*,'problematic cnumber= <',cnumber,'>'
      number=10
      goto 10     
    endfunction noform
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
    function step(x,x0,width)
!
!  Smooth unit step function centred at x0; implemented as tanh profile
!  23-jan-02/wolf: coded
!
      use Cdata, only: epsi
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: step
      real :: x0,width

        step = 0.5*(1+tanh((x-x0)/(width+epsi)))
!
      endfunction step
!***********************************************************************
    function der_step(x,x0,width)
!
!  Derivative of smooth unit STEP() function given above (i.e. a bump profile).
!  Adapt this if you change the STEP() profile, or you will run into
!  inconsistenies.
!  23-jan-02/wolf: coded
!
      use Cdata, only: epsi
!
      real, dimension(:) :: x
      real, dimension(size(x,1)) :: der_step,arg
      real :: x0,width
!
!  Some argument gymnastics to avoid `floating overflow' for large
!  arguments
!
      arg = abs((x-x0)/(width+epsi))
      arg = min(arg,8.)         ! cosh^2(8) = 3e+27
      der_step = 0.5/(width*cosh(arg)**2)
!
      endfunction der_step
!***********************************************************************
      function notanumber_1(f)
!
!  Check for denormalised floats (in fact NaN or -Inf, Inf).
!  The test used here should work on all architectures even if
!  optimisation is high (something like `if (any(f /= f+1))' would be
!  optimised away).
!  Version for 1d arrays.
!  24-jan-02/wolf: coded
!
        logical :: notanumber_1
        real, dimension(:) :: f
        real, dimension(size(f,1)) :: g

        g = f
        notanumber_1 = (any(f /= g) .or. any(f == g+1))
!
      endfunction notanumber_1
!***********************************************************************
      function notanumber_2(f)
!
!  Check for denormalised floats (in fact NaN or -Inf, Inf).
!  The test used here should work on all architectures even if
!  optimisation is high (something like `if (any(f /= f+1))' would be
!  optimised away).
!  Version for 2d arrays.
!
!  1-may-02/wolf: coded
!d
        logical :: notanumber_2
        real, dimension(:,:) :: f
        real, dimension(size(f,1),size(f,2)) :: g

        g = f
        notanumber_2 = (any(f /= g) .or. any(f == g+1))
!
      endfunction notanumber_2
!***********************************************************************
      function notanumber_3(f)
!
!  Check for denormalised floats (in fact NaN or -Inf, Inf).
!  The test used here should work on all architectures even if
!  optimisation is high (something like `if (any(f /= f+1))' would be
!  optimised away).
!  Version for 3d arrays.
!
!  24-jan-02/wolf: coded
!
        logical :: notanumber_3
        real, dimension(:,:,:) :: f
        real, dimension(size(f,1),size(f,2),size(f,3)) :: g

        g = f
        notanumber_3 = (any(f /= g) .or. any(f == g+1))
!
      endfunction notanumber_3
!***********************************************************************
      subroutine parse_bc(bc,bc1,bc2)
!
!  Parse boundary conditions, which may be in the form `a' (applies to
!  both `lower' and `upper' boundary) or `a:s' (use `a' for lower,
!  `s' for upper boundary.
!
!  24-jan-02/wolf: coded
!
        use Cparam, only: mvar,bclen
        use Mpicomm
!
        character (len=2*bclen+1), dimension(mvar) :: bc
        character (len=bclen), dimension(mvar) :: bc1,bc2
        integer :: j,isep
!
        intent(in) :: bc
        intent(out) :: bc1,bc2
!

        do j=1,mvar
          if (bc(j) == '') then ! will probably never happen due to default='p'
            if (lroot) print*, 'Empty boundary condition No. ', &
                 j, 'in (x, y, or z)'
            call stop_it('PARSE_BC')
          endif
          isep = index(bc(j),':')
          if (isep > 0) then
            bc1(j) = bc(j)(1:isep-1)
            bc2(j) = bc(j)(isep+1:)
          else
            bc1(j) = bc(j)(1:bclen)
            bc2(j) = bc(j)(1:bclen)
          endif
        enddo
!
      endsubroutine parse_bc
!***********************************************************************
      subroutine parse_name(iname,cname,cform,ctest,itest)
!
!  Parse name and format of print variable
!  On output, ITEST is set to INAME if CNAME matches CTEST
!  and CFORM is set to the format given as default.
!  E.g. if CTEST='bmax' *i.e. we are testing input line CNAME for 'bmax',
!  CNAME='bmax' will be parsed to ITEST=INAME, CFORM='(1pe10.2)',
!  CNAME='bmax(G5.1)' to ITEST=INAME, CFORM='G5.1',
!  CNAME='brms' to ITEST=<unchanged, normally 0>, CFORM='(1pe10.2)'
!
      character (len=*) :: cname,cform
      character (len=*) :: ctest
      integer :: iname,itest,iform0,iform1,iform2,length,index_i
!
      intent(in)  :: iname,cname,ctest
      intent(inout) :: itest,cform
!
!  check whether format is given
!
      iform0=index(cname,' ')
      iform1=index(cname,'(')
      iform2=index(cname,')')
!
!  set format; use default if not given
!
      if (iform1>0) then
        cform=cname(iform1+1:iform2-1)
        length=iform1-1
      else
        cform='1p,e10.2,0p'  !!(the nag-f95 compiler requires a comma after 1p)
        length=iform0-1
      endif
!
!  if the name matches, we keep the name and can strip off the format.
!  The remaining name can then be used for the legend.
!
      if (cname(1:length)==ctest .and. itest==0) then
        itest=iname
      endif
!
!  Integer formats are turned into floating point numbers
!
      index_i=index(cform,'i')
      if (index_i/=0) then
        cform(index_i:index_i)='f'
        cform=trim(cform)//'.0'
      endif
!
      endsubroutine parse_name
!***********************************************************************
      subroutine remove_file(fname)
!
!  Remove a file; this variant seems to be portable
!  5-mar-02/wolf: coded
!
        character (len=*) :: fname
!
        open(1,FILE=fname)
        close(1,STATUS="DELETE")
!
      endsubroutine remove_file
!***********************************************************************
      subroutine rmwig0(f)
!
!  There is no diffusion acting on the density, and wiggles in
!  lnrho are not felt in the momentum equation at all (zero gradient).
!  Thus, in order to keep lnrho smooth one needs to smooth lnrho
!  in sporadic time intervals.
!
!  11-Jul-01/axel: adapted from similar version in f77 code
!
!  WARNING: THIS ROUTINE IS LIKELY TO BE BROKEN IF YOU USE MPI
!
      use Cdata
!
      real, dimension (mx,my,mz) :: tmp
      real, dimension (mx,my,mz,mvar) :: f
!
!  copy
!
      print*,'remove wiggles in lnrho, t=',t
      tmp=exp(f(:,:,:,ilnrho))
      call smooth_3d(tmp,1)
      f(:,:,:,ilnrho)=alog(tmp)
!
    endsubroutine rmwig0
!***********************************************************************
    subroutine get_nseed(nseed)
!
!  Get length of state of random number generator. The current seed can
!  be represented by nseed (4-byte) integers.
!  Different compilers have different lengths:
!    SGI: 64, Intel: 47, NAG: 1, Compaq: 2
!
      use Cparam, only: mseed
      use Mpicomm, only: lroot,stop_it      
!
      integer :: nseed
!
      call random_seed(SIZE=nseed)
      !
      ! test whether mseed is large enough for this machine
      !
      if (nseed > mseed) then
        if (lroot) print*, "This machine requires mseed >= ", nseed
        call stop_it("Need to increase mseed")
      endif
!
    endsubroutine get_nseed
!***********************************************************************
    subroutine blob(ampl,f,i,radius,xblob,yblob,zblob)
!
!  single  blob
!
      use Cdata
!
!  27-jul-02/axel: coded
!
      integer :: i
      real, dimension (mx,my,mz,mvar) :: f
      real,optional :: xblob,yblob,zblob
      real :: ampl,radius,x01=0.,y01=0.,z01=0.
!
!  single  blob
!
      if (present(xblob)) x01=xblob
      if (present(yblob)) y01=yblob
      if (present(zblob)) z01=zblob
      if (ampl==0) then
        if (lroot) print*,'ampl=0 in blob'
      else
        if (lroot) print*,'blob: variable i,ampl=',i,ampl
        f(:,:,:,i)=f(:,:,:,i)+ampl*(&
           spread(spread(exp(-((x-x01)/radius)**2),2,my),3,mz)&
          *spread(spread(exp(-((y-y01)/radius)**2),1,mx),3,mz)&
          *spread(spread(exp(-((z-z01)/radius)**2),1,mx),2,my))
      endif
!
    endsubroutine blob
!***********************************************************************

endmodule Sub
