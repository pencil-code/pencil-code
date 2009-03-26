! $Id: equ.f90 10533 2009-03-26 11:01:45Z ajohan@strw.leidenuniv.nl $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Diagnostics
!
  use Cdata
  use Messages
  use Mpicomm
  use Sub
!
  implicit none
!
  private
!
  public :: diagnostic, initialize_time_integrals
  public :: collect_uumax, xyaverages_z, xzaverages_y, yzaverages_x
  public :: phizaverages_r, yaverages_xz, zaverages_xy
  public :: time_integrals, phiaverages_rz
!
  contains
!***********************************************************************
    subroutine diagnostic
!
!  Calculate diagnostic quantities.
!
!   2-sep-01/axel: coded
!  14-aug-03/axel: began adding surface integrals
!
      real, dimension (mname) :: fmax_tmp, fsum_tmp, fmax, fsum, fweight_tmp
      real :: dv,vol
      integer :: iname,imax_count,isum_count,nmax_count,nsum_count
      logical :: lweight_comm
      logical, save :: first=.true.
      real, save :: dVol_rel1
      real :: intdr_rel,intdtheta_rel,intdphi_rel,intdz_rel

      print*, '111', it, itsub
!
!  Calculate relative volume integral.
!
      if (first) then
        if (lspherical_coords) then
          intdr_rel     =      (xyz1(1)**3-    xyz0(1)**3)/(3.*dx)
          intdtheta_rel = -(cos(xyz1(2))  -cos(xyz0(2)))/dy
          intdphi_rel   =      (xyz1(3)   -    xyz0(3)) /dz
!
!  Prevent zeros from less then 3-dimensional runs
!  (maybe this should be 2pi, but maybe not).
!
          if (nx==1) intdr_rel=1.
          if (ny==1) intdtheta_rel=1.
          if (nz==1) intdphi_rel=1.
          dVol_rel1=1./(intdr_rel*intdtheta_rel*intdphi_rel)
        elseif (lcylindrical_coords) then
          intdr_rel   =      (xyz1(1)**2-    xyz0(1)**2)/(2.*dx)
          intdphi_rel =      (xyz1(2)   -    xyz0(2)) /dy
          intdz_rel   =      (xyz1(3)   -    xyz0(3)) /dz
          dVol_rel1=1./(intdr_rel*intdphi_rel*intdz_rel)
        else
          dVol_rel1=1./(nw*ncpus)
        endif
        first=.false.
        if (lroot.and.ip<=10) print*,'dVol_rel1=',dVol_rel1
        if (lroot) print*,'box volume = ', dx*dy*dz/dVol_rel1
      endif
!
!  Go through all print names, and sort into communicators
!  corresponding to their type.
!
      imax_count=0
      isum_count=0
      lweight_comm=.false.
      do iname=1,nname
        if (itype_name(iname)<0) then
          imax_count=imax_count+1
          fmax_tmp(imax_count)=fname(iname)
        elseif (itype_name(iname)>0) then
          isum_count=isum_count+1
          fsum_tmp(isum_count)=fname(iname)
          if (itype_name(iname)==ilabel_sum_weighted .or. &
              itype_name(iname)==ilabel_sum_weighted_sqrt .or. &
              itype_name(iname)==ilabel_sum_par) then
            fweight_tmp(isum_count)=fweight(iname)
            lweight_comm=.true.
          endif
        endif
      enddo
      nmax_count=imax_count
      nsum_count=isum_count
!
!  Communicate over all processors.
!
      call mpireduce_max(fmax_tmp,fmax,nmax_count)
      call mpireduce_sum(fsum_tmp,fsum,nsum_count)
      if (lweight_comm) call mpireduce_sum(fweight_tmp,fweight,nsum_count)
!
!  The result is present only on the root processor.
!
      if (lroot) then
!
!  Sort back into original array.
!
         imax_count=0
         isum_count=0
         do iname=1,nname
           if (itype_name(iname)<0) then ! max
             imax_count=imax_count+1
!
             if (itype_name(iname)==ilabel_max)            &
                 fname(iname)=fmax(imax_count)
!
             if (itype_name(iname)==ilabel_max_sqrt)       &
                 fname(iname)=sqrt(fmax(imax_count))
!
             if (itype_name(iname)==ilabel_max_dt)         &
                 fname(iname)=fmax(imax_count)
!
             if (itype_name(iname)==ilabel_max_neg)        &
                 fname(iname)=-fmax(imax_count)
!
             if (itype_name(iname)==ilabel_max_reciprocal) &
                 fname(iname)=1./fmax(imax_count)
!
           elseif (itype_name(iname)>0) then ! sum
             isum_count=isum_count+1
!
             if (itype_name(iname)==ilabel_sum)            &
                 fname(iname)=fsum(isum_count)*dVol_rel1
!
             if (itype_name(iname)==ilabel_sum_sqrt)       &
                 fname(iname)=sqrt(fsum(isum_count)*dVol_rel1)
!
             if (itype_name(iname)==ilabel_sum_par)        &
                 fname(iname)=fsum(isum_count)/fweight(isum_count)
!
             if (itype_name(iname)==ilabel_integrate) then
               dv=1.
               if (nxgrid/=1.and.lequidist(1)) dv=dv*dx
               if (nygrid/=1.and.lequidist(2)) dv=dv*dy
               if (nzgrid/=1.and.lequidist(3)) dv=dv*dz
               fname(iname)=fsum(isum_count)*dv
              endif
!
              if (itype_name(iname)==ilabel_surf)          &
                  fname(iname)=fsum(isum_count)
!
              if (itype_name(iname)==ilabel_sum_lim) then
                 vol=1.
                 if (lcylinder_in_a_box)  vol=vol*pi*(r_ext**2-r_int**2)
                 if (nzgrid/=1)           vol=vol*Lz
                 if (lsphere_in_a_box)    vol=1.333333*pi*(r_ext**3-r_int**3)
                 fname(iname)=fsum(isum_count)/vol
              endif
!
             if (itype_name(iname)==ilabel_sum_weighted) then
               if (fweight(isum_count)/=0.0) then
                 fname(iname)=fsum(isum_count)/fweight(isum_count)
               else
                 fname(iname)=0.0
               endif
             endif
!
             if (itype_name(iname)==ilabel_sum_weighted_sqrt) then
               if (fweight(isum_count)/=0.0) then
                 fname(iname)=sqrt(fsum(isum_count)/fweight(isum_count))
               else
                 fname(iname)=0.0
               endif
             endif
!
           endif
!
         enddo
!
      endif
!
    endsubroutine diagnostic
!***********************************************************************
    subroutine collect_UUmax
!
!  Calculate the maximum effective advection velocity in the domain;
!  needed for determining dt at each timestep.
!
!   2-sep-01/axel: coded
!
    real, dimension(1) :: fmax_tmp,fmax
!
!  Communicate over all processors.
!  The result is then present only on the root processor
!  reassemble using old names.
!
    fmax_tmp(1)=UUmax
    call mpireduce_max(fmax_tmp,fmax,1)
    if (lroot) UUmax=fmax(1)
!
    endsubroutine collect_UUmax
!***********************************************************************
    subroutine initialize_time_integrals(f)
!
!  Initialize time_integrals for full chunks.
!
!  28-jun-07/axel+mreinhard: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      if (iuut/=0) f(:,:,:,iuxt:iuzt)=0.
      if (ioot/=0) f(:,:,:,ioxt:iozt)=0.
      if (ibbt/=0) f(:,:,:,ibxt:ibzt)=0.
      if (ijjt/=0) f(:,:,:,ijxt:ijzt)=0.
!
    endsubroutine initialize_time_integrals
!***********************************************************************
    subroutine time_integrals(f,p)
!
!  Calculate time_integrals within each pencil (as long as each
!  pencil case p still contains the current data).
!
!  28-jun-07/axel+mreinhard: coded
!  24-jun-08/axel: moved call to this routine to the individual pde routines
!
      use Hydro, only: time_integrals_hydro
      use Magnetic, only: time_integrals_magnetic
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(in) :: p
!
      if (itsub==itorder) then
        if (lhydro) call time_integrals_hydro(f,p)
        if (lmagnetic) call time_integrals_magnetic(f,p)
      endif
!
    endsubroutine time_integrals
!***********************************************************************
    subroutine xyaverages_z()
!
!  Calculate xy-averages (still depending on z)
!  NOTE: these averages depend on z, so after summation in x and y they
!  are still distributed over nprocz CPUs; hence the dimensions of fsumz
!  (and fnamez).
!  In other words: the whole xy-average is present in one and the same fsumz,
!  but the result is not complete on any of the processors before
!  mpireduce_sum has been called. This is simpler than collecting results
!  first in processors with the same ipz and different ipy, and then
!  assemble result from the subset of ipz processors which have ipy=0
!  back on the root processor.
!
!   6-jun-02/axel: coded
!
      real, dimension (nz,nprocz,mnamez) :: fsumz
!
!  Communicate over all processors.
!  The result is only present on the root processor
!
      if (nnamez>0) then
        call mpireduce_sum(fnamez,fsumz,nz*nprocz*nnamez)
        if (lroot) &
            fnamez(:,:,1:nnamez)=fsumz(:,:,1:nnamez)/(nx*ny*nprocy)
      endif
!
    endsubroutine xyaverages_z
!***********************************************************************
    subroutine xzaverages_y()
!
!  Calculate xz-averages (still depending on y)
!
!  12-oct-05/anders: adapted from xyaverages_z
!
      real, dimension (ny,nprocy,mnamey) :: fsumy
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamey>0) then
        call mpireduce_sum(fnamey,fsumy,ny*nprocy*nnamey)
        if (lroot) &
            fnamey(:,:,1:nnamey)=fsumy(:,:,1:nnamey)/(nx*nz*nprocz)
      endif
!
    endsubroutine xzaverages_y
!***********************************************************************
    subroutine yzaverages_x()
!
!  Calculate yz-averages (still depending on x).
!
!   2-oct-05/anders: adapted from xyaverages_z
!
      real, dimension (nx,mnamex) :: fsumx
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamex>0) then
        call mpireduce_sum(fnamex,fsumx,nx*nnamex)
        if (lroot) &
            fnamex(:,1:nnamex)=fsumx(:,1:nnamex)/(ny*nprocy*nz*nprocz)
      endif
!
    endsubroutine yzaverages_x
!***********************************************************************
    subroutine phizaverages_r()
!
!  Calculate phiz-averages (still depending on r)
!  
!  29-jan-07/wlad: adapted from yzaverages_x and phiaverages_rz
!
      real, dimension (nrcyl,mnamer) :: fsumr
      real, dimension (nrcyl) :: norm
      integer :: in,ir
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamer>0) then
         !the extra slot is where the normalization is stored
         call mpireduce_sum(fnamer,fsumr,nrcyl*(nnamer+1))
         if (lroot) then
            norm=fsumr(:,nnamer+1)
            do in=1,nnamer
               fnamer(:,in)=fsumr(:,in)/norm
            enddo
         endif
      endif
!
    endsubroutine phizaverages_r
!***********************************************************************
    subroutine yaverages_xz()
!
!  Calculate y-averages (still depending on x and z)
!  NOTE: these averages depend on x and z, so after summation in y they
!  are still distributed over nprocy CPUs; hence the dimensions of fsumxz
!  (and fnamexz).
!
!   7-jun-05/axel: adapted from zaverages_xy
!
      real, dimension (nx,nz,nprocz,mnamexz) :: fsumxz
!
!  Communicate over all processors.
!  The result is only present on the root processor.
!
      if (nnamexz>0) then
        call mpireduce_sum(fnamexz,fsumxz,nnamexz*nx*nz*nprocz)
        if (lroot) &
            fnamexz(:,:,:,1:nnamexz)=fsumxz(:,:,:,1:nnamexz)/(ny*nprocy)
      endif
!
    endsubroutine yaverages_xz
!***********************************************************************
    subroutine zaverages_xy()
!
!  Calculate z-averages (still depending on x and y)
!  NOTE: these averages depend on x and y, so after summation in z they
!  are still distributed over nprocy CPUs; hence the dimensions of fsumxy
!  (and fnamexy).
!
!  19-jun-02/axel: coded
!
      real, dimension (nx,ny,nprocy,mnamexy) :: fsumxy
!
!  Communicate over all processors.
!  the result is only present on the root processor.
!
      if (nnamexy>0) then
        call mpireduce_sum(fnamexy,fsumxy,nnamexy*nx*ny*nprocy)
        if (lroot) &
            fnamexy(:,:,:,1:nnamexy)=fsumxy(:,:,:,1:nnamexy)/(nz*nprocz)
      endif
!
    endsubroutine zaverages_xy
!***********************************************************************
    subroutine phiaverages_rz()
!
!  calculate azimuthal averages (as functions of r_cyl,z)
!  NOTE: these averages depend on (r and) z, so after summation they
!  are still distributed over nprocz CPUs; hence the dimensions of fsumrz
!  (and fnamerz).
!
!  9-dec-02/wolf: coded
!
      integer :: i
      real, dimension (nrcyl,0:nz,nprocz,mnamerz) :: fsumrz
!
!  Communicate over all processors.
!  The result is only present on the root processor
!  normalize by sum of unity which is accumulated in fnamerz(:,0,:,1).
!
      if (nnamerz>0) then
        call mpireduce_sum(fnamerz,fsumrz,mnamerz*nrcyl*(nz+1)*nprocz)
        if (lroot) then
          do i=1,nnamerz
            fnamerz(:,1:nz,:,i)=fsumrz(:,1:nz,:,i)/spread(fsumrz(:,0,:,1),2,nz)
          enddo
        endif
      endif
!
    endsubroutine phiaverages_rz
!***********************************************************************
endmodule Diagnostics
