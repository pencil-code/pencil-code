! $Id: wsnaps.f90,v 1.41 2003-11-26 16:11:51 theine Exp $

!!!!!!!!!!!!!!!!!!!!!!!
!!!   wsnaps.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Write snapshot files (variables and powerspectra). Is a separate
!!!  module, so it can be used with whatever IO module we want.

module Wsnaps

  implicit none

contains

!***********************************************************************
    subroutine wsnap(chsnap,a,msnap,enum)
!
!  Write snapshot file, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used for var.dat)
!
!  30-sep-97/axel: coded
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tsnap.dat'
!   5-apr-03/axel: possibility for additional (hard-to-get) output 
!  31-may-03/axel: wsnap can write either w/ or w/o auxiliary variables
!
      use Cdata
      use Mpicomm
      use Boundcond
      use Radiation
!      use Viscosity, only: calc_viscosity
      use Ionization
      use Sub
      use Io
!
!  the dimension msnap can either be mvar+maux (for f-array in run.f90)
!  or just mvar (for f-array in start.f90 or df-array in run.f90
!
      integer :: msnap
      real, dimension (mx,my,mz,msnap) :: a
      character (len=4) :: ch
      character (len=135) :: file
      character (len=*) :: chsnap
      logical lsnap,enum
      integer, save :: ifirst,nsnap
      real, save :: tsnap
!
!  Output snapshot with label in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      if (enum) then
        file=trim(datadir)//'/tsnap.dat'
!
!  at first call, need to initialize tsnap
!  tsnap calculated in read_snaptime, but only available to root processor
!
        if (ifirst==0) then
          call read_snaptime(file,tsnap,nsnap,dsnap,t)
          ifirst=1
        endif
!
!  Check whether we want to output snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently)
!
        call update_snaptime(file,tsnap,nsnap,dsnap,t,lsnap,ch,ENUM=.true.)
        if (lsnap) then
!          call calc_viscosity(a)
          call update_ghosts(a)
          call output(chsnap//ch,a,msnap)
          if(ip<=10.and.lroot) print*,'wsnap: written snapshot ',chsnap//ch
        endif
!
      else
!
!  write snapshot without label (typically, var.dat)
!
!        call calc_viscosity(a)
        call update_ghosts(a)
        call output(chsnap,a,msnap)
      endif
!
    endsubroutine wsnap
!***********************************************************************
   subroutine powersnap(a,lwrite_only)
!
!  Write a snapshot of power spectrum 
!
!  30-sep-97/axel: coded
!  07-oct-02/nils: adapted from wsnap
!  08-oct-02/tony: expanded file to handle 120 character datadir // '/tspec.dat'
!  28-dec-02/axel: call structure from herel; allow optional lwrite_only
!
      use Cdata
      use Mpicomm
      use Boundcond
      use Sub
      use Io
      use Power_spectrum
      use Struct_func
!
      real, dimension (mx,my,mz,mvar+maux) :: a
      real, dimension (nx,ny,nz) :: b_vec
      character (len=135) :: file
      character (len=4) :: ch
      logical :: lspec,llwrite_only=.false.,ldo_all
      logical, optional :: lwrite_only
      integer, save :: ifirst,nspec
      real, save :: tspec
      integer :: ivec,im,in
      real, dimension (nx) :: bb 
!
!  set llwrite_only
!
      if (present(lwrite_only)) llwrite_only=lwrite_only
      ldo_all=.not.llwrite_only
!
!  Output snapshot in 'tpower' time intervals
!  file keeps the information about time of last snapshot
!
      file=trim(datadir)//'/tspec.dat'
!
!  at first call, need to initialize tspec
!  tspec calculated in read_snaptime, but only available to root processor
!
      if(ldo_all.and.ifirst==0) then
         call read_snaptime(file,tspec,nspec,dspec,t)
         ifirst=1
      endif
!
!  Check whether we want to output power snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently)
!
      if(ldo_all) &
           call update_snaptime(file,tspec,nspec,dspec,t,lspec,ch,ENUM=.false.)
      if (lspec.or.llwrite_only) then
         if (ldo_all)  call update_ghosts(a)
         if (vel_spec) call power(a,'u')
         if (mag_spec) call power(a,'b')
         if (vec_spec) call power(a,'a')
         if (ab_spec)  call powerhel(a,'mag')
         if (ou_spec)  call powerhel(a,'kin')
         if (ro_spec)  call powerscl(a,'ro')
         if (ss_spec)  call powerscl(a,'ss')
         if (cc_spec)  call powerscl(a,'cc')
         if (cr_spec)  call powerscl(a,'cr')
         if (oned) then
            if (vel_spec) call power_1d(a,'u',1)
            if (mag_spec) call power_1d(a,'b',1)
            if (vec_spec) call power_1d(a,'a',1)
            if (vel_spec) call power_1d(a,'u',2)
            if (mag_spec) call power_1d(a,'b',2)
            if (vec_spec) call power_1d(a,'a',2)
            if (vel_spec) call power_1d(a,'u',3)
            if (mag_spec) call power_1d(a,'b',3)
            if (vec_spec) call power_1d(a,'a',3)
         endif
         !
         !  Doing structure functions
         !
         do ivec=1,3
            if (lsfb .or. lsfz1 .or. lsfz2 .or. lpdfb .or. lpdfz1 .or. lpdfz2) then
               do n=n1,n2
                  do m=m1,m2
                     call curli(a,iaa,bb,ivec)
                     im=m-nghost
                     in=n-nghost
                     b_vec(:,im,in)=bb
                  enddo
               enddo
               b_vec=b_vec/sqrt(exp(a(l1:l2,m1:m2,n1:n2,ilnrho)))
            endif
            if (lsfu)     call structure(a,ivec,b_vec,'u')
            if (lsfb)     call structure(a,ivec,b_vec,'b')
            if (lsfz1)    call structure(a,ivec,b_vec,'z1')
            if (lsfz2)    call structure(a,ivec,b_vec,'z2')
            if (lpdfu)    call structure(a,ivec,b_vec,'pdfu')
            if (lpdfb)    call structure(a,ivec,b_vec,'pdfb')
            if (lpdfz1)   call structure(a,ivec,b_vec,'pdfz1')
            if (lpdfz2)   call structure(a,ivec,b_vec,'pdfz2')
         enddo
      endif
!
    endsubroutine powersnap
!***********************************************************************
endmodule Wsnaps

