! $Id$
!
!  Calculates 2-point structure functions and/or PDFs
!  and saves them during the run.
!
!  For the time being, the structure functions (or PDFs) are
!  called from power, so the output frequency is set by dspec.
!
!  The save files are under data/proc# under the names
!  sfz1_sum_ or sfz1_sum_transp_ .
!
!-----------------------------------------------------------------------
!   23-dec-02/nils: adapted from postproc/src/struct_func_mpi.f90
!

module struct_func
  !
  implicit none

  include 'struct_func.h'
  !
  contains

!***********************************************************************
    subroutine structure(f,ivec,b_vec,varlabel)
      !
      !  The following parameters may need to be readjusted:
      !  qmax should be set to the largest moment to be calculated
      !  n_pdf gives the number of bins of the PDF
      !
      !   23-dec-02/nils: adapted from postproc/src/struct_func_mpi.f90
      !   28-dec-02/axel: need also n_pdf in normalization
      !
      use Cdata
      use Sub
      use General
      use Messages, only: fatal_error
      use Mpicomm
      !
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivec
      real, dimension (nx,ny,nz) :: b_vec
      character (len=*) :: varlabel
      !
      integer, parameter :: qmax=8+1 ! the extrta 1 is for unsigned 3. moment.
      !integer, parameter :: lb_nxgrid=floor(alog(nxgrid+1.e-6)/lntwo)
      integer, parameter :: imax=lb_nxgrid*2-2
      integer, parameter :: n_pdf=101
      real, dimension (:,:,:,:), allocatable :: flux_vec
      real, dimension (nx,ny,nz) :: vect
      real, dimension (imax,qmax,3) :: sf,sf_sum
      real, dimension (ny,nz,3) :: dvect_flux1,dvect_flux2
      real, dimension (ny,nz) :: dvect1,dvect2,sf3_flux1,sf3_flux2
      real, dimension(n_pdf,imax,3) :: p_du,p_du_sum
      real, dimension(n_pdf) :: x_du
      integer, dimension (ny,nz) :: i_du1,i_du2
      integer :: l,q,direction,ll1,ll2,mtmp,ntmp
      integer :: i,lb_ll,separation,exp1,exp2
      integer(KIND=ikind8) :: ndiv
      real :: pdf_max,pdf_min,normalization,dx_du
      character (len=fnlen):: prefix
      logical :: llsf=.false., llpdf=.false.
      logical, save :: l0=.true.

      if (l0) then
        l0=.false.
        if (2**lb_nxgrid/=nxgrid) then
          print*, 'lb_nxgrid, nxgrid=', lb_nxgrid, nxgrid
          call fatal_error('structure','nxgrid is not a power of 2')
        endif
      endif
      !
      !  Do structure functions
      !
      if (lroot.and.ip<9) print*,'Doing structure functions'
      !
      if (varlabel == 'u') then
        vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
        prefix='/sfu-'
        sf=0.
        llsf=.true.
        llpdf=.false.
      elseif (varlabel == 'b') then
        vect(:,:,:)=b_vec(:,:,:)
        prefix='/sfb-'
        sf=0.
        llsf=.true.
        llpdf=.false.
      elseif (varlabel == 'z1') then
        vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)+b_vec(:,:,:)
        prefix='/sfz1-'
        sf=0.
        llsf=.true.
        llpdf=.false.
      elseif (varlabel == 'z2') then
        vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)-b_vec(:,:,:)
        prefix='/sfz2-'
        sf=0.
        llsf=.true.
        llpdf=.false.
      elseif (varlabel == 'flux') then
        !
        ! Here we calculate the flux like structure functions of
        ! Politano & Pouquet (1998)
        !
        allocate(flux_vec(nx,ny,nz,3))
        mtmp=m
        ntmp=n
        do m=m1,m2
          do n=n1,n2
            call curl(f,iaa,flux_vec(:,m-nghost,n-nghost,:))
          enddo
        enddo
        m=mtmp
        n=ntmp
        vect(:,:,:)  = f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)-b_vec(:,:,:)
        flux_vec=f(l1:l2,m1:m2,n1:n2,iux:iuz)+flux_vec
        prefix='/sfflux-'
        sf=0.
        llsf=.true.
        llpdf=.false.
      endif
      !
      !  Setting some variables depending on wether we want to
      !  calculate pdf or structure functions.
      !
      if (varlabel == 'pdfu') then
        vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)
        prefix='/pdfu-'
        pdf_max= 1.  !(for the time being; assumes |u|<1)
        pdf_min=-pdf_max
        dx_du=(pdf_max-pdf_min)/n_pdf
        do l=1,n_pdf
          x_du(l)=(l-.5)*dx_du+pdf_min
        enddo
        p_du=0.
        llpdf=.true.
        llsf=.false.
      elseif (varlabel == 'pdfb') then
        vect=b_vec
        prefix='/pdfb-'
        pdf_max= 1.  !(for the time being; assumes |u|<1)
        pdf_min=-pdf_max
        dx_du=(pdf_max-pdf_min)/n_pdf
        do l=1,n_pdf
          x_du(l)=(l-.5)*dx_du+pdf_min
        enddo
        p_du=0.
        llpdf=.true.
        llsf=.false.
      elseif (varlabel == 'pdfz1') then
        prefix='/pdfz1-'
        vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)+b_vec(:,:,:)
        pdf_max= 1.  !(for the time being; assumes |u|<1)
        pdf_min=-pdf_max
        dx_du=(pdf_max-pdf_min)/n_pdf
        do l=1,n_pdf
          x_du(l)=(l-.5)*dx_du+pdf_min
        enddo
        p_du=0.
        llpdf=.true.
        llsf=.false.
      elseif (varlabel == 'pdfz2') then
        vect(:,:,:)=f(l1:l2,m1:m2,n1:n2,iuu+ivec-1)-b_vec(:,:,:)
        prefix='/pdfz2-'
        pdf_max= 1.  !(for the time being; assumes |u|<1)
        pdf_min=-pdf_max
        dx_du=(pdf_max-pdf_min)/n_pdf
        do l=1,n_pdf
          x_du(l)=(l-.5)*dx_du+pdf_min
        enddo
        p_du=0.
        llpdf=.true.
        llsf=.false.
      endif
      !
      !  Beginning the loops
      !
      do direction=1,nr_directions
        do l=1,nx
          if (lroot .and. lpostproc) print*,'l=',l
          do lb_ll=1,lb_nxgrid*2-2
            if (lb_ll == 1) then
              exp2=0
            else
              exp2=mod(lb_ll,2)
            endif
            exp1=int(lb_ll/2)-exp2
            separation=(2**exp1)*(3**exp2)
            ll1=mod(l+separation-1,nx)+1
            ll2=mod(l-separation+nx-1,nx)+1
            dvect1=vect(l,:,:)-vect(ll1,:,:)
            dvect2=vect(l,:,:)-vect(ll2,:,:)
            if (varlabel == 'flux') then
              dvect_flux1=flux_vec(l,:,:,:)-flux_vec(ll1,:,:,:)
              dvect_flux2=flux_vec(l,:,:,:)-flux_vec(ll2,:,:,:)
            endif
            if (llpdf) then !if pdf=.true.
              i_du1=1+int((dvect1-pdf_min)*n_pdf/(pdf_max-pdf_min))
              i_du1=min(max(i_du1,1),n_pdf)  !(make sure its inside array bdries)
              i_du2=1+int((dvect2-pdf_min)*n_pdf/(pdf_max-pdf_min))
              i_du2=min(max(i_du2,1),n_pdf)  !(make sure its inside array bdries)
              !
              !  Calculating pdf
              !
              do m=1,ny
                do n=1,nz
                  p_du(i_du1(m,n),lb_ll,direction) &
                       =p_du(i_du1(m,n),lb_ll,direction)+1
                  p_du(i_du2(m,n),lb_ll,direction) &
                       =p_du(i_du2(m,n),lb_ll,direction)+1
                enddo
              enddo
            endif
            !
            if (llsf) then
              !
              !  Calculates sf
              !
              if (varlabel == 'flux') then
                sf3_flux1= &
                     abs(dvect1(:,:))* &
                     (dvect_flux1(:,:,1)**2 &
                     +dvect_flux1(:,:,2)**2 &
                     +dvect_flux1(:,:,3)**2)
                sf3_flux2= &
                     abs(dvect2(:,:))* &
                     (dvect_flux2(:,:,1)**2 &
                     +dvect_flux2(:,:,2)**2 &
                     +dvect_flux2(:,:,3)**2)
              endif
              !
              ! Loop over all q values
              !
              do q=1,qmax-1
                if (varlabel == 'flux') then
                  sf(lb_ll,q,direction) &
                       =sf(lb_ll,q,direction) &
                       +sum(abs(sf3_flux1(:,:))**(q/3.)) &
                       +sum(abs(sf3_flux2(:,:))**(q/3.))
                else
                  sf(lb_ll,q,direction) &
                       =sf(lb_ll,q,direction) &
                       +sum(abs(dvect1(:,:))**q)+sum(abs(dvect2(:,:))**q)
                endif
              enddo
              !
              ! Do unsigned third moment (store in last slot of array)
              !
              if (varlabel == 'flux') then
                sf(lb_ll,qmax,direction) &
                     =sf(lb_ll,qmax,direction) &
                     +sum(sf3_flux1(:,:)) &
                     +sum(sf3_flux2(:,:))
              else
                sf(lb_ll,qmax,direction) &
                     =sf(lb_ll,qmax,direction) &
                     +sum(dvect1(:,:)**3)+sum(dvect2(:,:)**3)
              endif
            endif
          enddo
        enddo
        if (nr_directions > 1) then
          if (direction == 1) then
            !Doing transpose of y direction
            call transp(vect,'y')
          endif
          if (direction == 2) then
            !Doing transpose of z direction
            call transp(vect,'z')
          endif
        endif
      enddo
      !
      !  Collecting all data on root processor and normalizing pdf and sf
      !
      if (llpdf) then
        call mpireduce_sum(p_du,p_du_sum,(/n_pdf,imax,3/))  !Is this safe???
        do i=1,imax
          do direction=1,nr_directions
            normalization=1./(n_pdf*dx_du*sum(p_du_sum(:,i,direction)))
            p_du_sum(:,i,direction)=normalization*p_du_sum(:,i,direction)
          enddo
        enddo
      endif
      !
      if (llsf) then
        call mpireduce_sum(sf,sf_sum,(/imax,qmax,3/))  !Is this safe???
        ndiv=nwgrid*2
        sf_sum=sf_sum/ndiv
      endif
      !
      !  Writing output file
      !
      if (lroot) then
        if (llpdf) then
          if (ip<10) print*,'Writing pdf of variable ',trim(itoa(ivec)), &
               ' to ',trim(datadir)//trim(prefix)//trim(itoa(ivec))//'.dat'
          open(1,file=trim(datadir)//trim(prefix)//trim(itoa(ivec))//'.dat', &
               position='append')
          write(1,*) t,n_pdf
          write(1,'(1p,8e10.2)') p_du_sum(:,:,:)
          write(1,'(1p,8e10.2)') x_du
          close(1)
        endif
        !
        if (llsf) then
          if (ip<10) print*,'Writing structure functions of variable ',&
               trim(itoa(ivec)), &
               ' to ',trim(datadir)//trim(prefix)//trim(itoa(ivec))//'.dat'
          open(1,file=trim(datadir)//trim(prefix)//trim(itoa(ivec))//'.dat', &
               position='append')
          write(1,*) t,qmax
          write(1,'(1p,8e10.2)') sf_sum(:,:,:)
          close(1)
        endif
      endif
      !
    endsubroutine structure
!***********************************************************************

endmodule struct_func
