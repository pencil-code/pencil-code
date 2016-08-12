! $Id$
!
!  This module contains useful subroutines for the particle modules.
!  Subroutines that depend on domain decomposition should be put in
!  the Particles_map module.
!
module Particles_sub
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_mpicomm
!
  implicit none
!
  private
!
  public :: init_particle_positions_velocities
  public :: input_particles, output_particles, boundconds_particles
  public :: sum_par_name, max_par_name, sum_par_name_nw, integrate_par_name
  public :: remove_particle, get_particles_interdistance
  public :: count_particles, output_particle_size_dist
  public :: get_rhopswarm, find_grid_volume, find_interpolation_weight
  public :: find_interpolation_indeces, get_gas_density, precalc_cell_volumes
  public :: precalc_weights
!
  interface get_rhopswarm
    module procedure get_rhopswarm_ineargrid
    module procedure get_rhopswarm_point
    module procedure get_rhopswarm_pencil
    module procedure get_rhopswarm_block
  endinterface
!
  contains
!***********************************************************************
    subroutine init_particle_positions_velocities(f,fp,ineargrid,klow,khigh)
!
!  Initial positions and velocities of dust particles.
!
!  29-dec-04/anders: coded
!
      use EquationOfState, only: beta_glnrho_global, cs20
      use General, only: random_number_wrapper, normal_deviate
      use Mpicomm, only: mpireduce_sum, mpibcast_real
      use InitialCondition, only: initial_condition_xxp, initial_condition_vvp
      use Particles_diagnos_dv, only: repeated_init
!
      real, dimension (mx,my,mz,mfarray), intent (out) :: f
      real, dimension (mpar_loc,mparray), intent (out) :: fp
      integer, dimension (mpar_loc,3), intent (out) :: ineargrid
      real, dimension (mpar_loc) :: rr_tmp, az_tmp
!
      real, dimension (3) :: uup, Lxyz_par, xyz0_par, xyz1_par
      real :: vpx_sum, vpy_sum, vpz_sum
      real :: r, p, q, px, py, pz, eps, cs, k2_xxp, rp2
      real :: dim1, npar_loc_x, npar_loc_y, npar_loc_z, dx_par, dy_par, dz_par
      real :: rad,rad_scl,phi,tht,tmp,OO,xx0,yy0,r2
      integer :: l, j, k, ix0, iy0, iz0, n_kill, klow, khigh
      logical :: lequidistant=.false.
      real :: rpar_int,rpar_ext
!
!  Use either a local random position or a global random position for certain
!  initial conditions. The default is a local random position, but the equal
!  number of particles per processors means that this is not completely random.
!
      if (lglobalrandom) then
        Lxyz_par=Lxyz
        xyz0_par=xyz0
        xyz1_par=xyz1
      else
        Lxyz_par=Lxyz_loc
        xyz0_par=xyz0_loc
        xyz1_par=xyz1_loc
      endif
!
!  Initial particle position.
!
      do j=1,ninit
!
        select case (initxxp(j))
!
        case ('nothing')
          if (lroot .and. j==1) print*, 'init_particles: nothing'
!
        case ('origin')
          if (lroot) print*, 'init_particles: All particles at origin'
          fp(klow:khigh,ixp:izp)=0.0
!
        case ('zero-z')
          if (lroot) print*, 'init_particles: Zero z coordinate'
          fp(klow:khigh,izp)=0.0
!
        case ('constant')
          if (lroot) &
              print*, 'init_particles: All particles at x,y,z=', xp0, yp0, zp0
          fp(klow:khigh,ixp)=xp0
          fp(klow:khigh,iyp)=yp0
          fp(klow:khigh,izp)=zp0
!
        case ('constant-1')
          if (lroot) &
              print*, 'init_particles: Particle 1 at x,y,z=', xp1, yp1, zp1
          do k=klow,khigh
            if (ipar(k)==1) then
              fp(k,ixp)=xp1
              fp(k,iyp)=yp1
              fp(k,izp)=zp1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) &
              print*, 'init_particles: Particle 2 at x,y,z=', xp2, yp2, zp2
          do k=klow,khigh
            if (ipar(k)==2) then
              fp(k,ixp)=xp2
              fp(k,iyp)=yp2
              fp(k,izp)=zp2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) &
              print*, 'init_particles: Particle 2 at x,y,z=', xp3, yp3, zp3
          do k=klow,khigh
            if (ipar(k)==3) then
              fp(k,ixp)=xp3
              fp(k,iyp)=yp3
              fp(k,izp)=zp3
            endif
          enddo
!
        case ('random-constz')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=klow,khigh
            if (nxgrid/=1) then
              call random_number_wrapper(r)
              fp(k,ixp)=r
            endif
            if (nygrid/=1) then
              call random_number_wrapper(r)
              fp(k,iyp)=r
            endif
          enddo
          if (nxgrid/=1) &
              fp(klow:khigh,ixp)=xyz0_par(1)+fp(klow:khigh,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(klow:khigh,iyp)=xyz0_par(2)+fp(klow:khigh,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(klow:khigh,izp)=zp0
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=klow,khigh
            if (nxgrid/=1) then
              call random_number_wrapper(r)
              fp(k,ixp)=r
            endif
            if (nygrid/=1) then
              call random_number_wrapper(r)
              fp(k,iyp)=r
            endif
            if (nzgrid/=1) then
              call random_number_wrapper(r)
              fp(k,izp)=r
            endif
          enddo
          if (nxgrid/=1) &
              fp(klow:khigh,ixp)=xyz0_par(1)+fp(klow:khigh,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(klow:khigh,iyp)=xyz0_par(2)+fp(klow:khigh,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(klow:khigh,izp)=xyz0_par(3)+fp(klow:khigh,izp)*Lxyz_par(3)
!
        case ('random-circle')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=klow,khigh
            call random_number_wrapper(r)
            if (zp0>yp0) then
              fp(k,ixp)=xp0*cos((zp0-yp0)*r+yp0)
              fp(k,iyp)=xp0*sin((zp0-yp0)*r+yp0)
            else
              fp(k,ixp)=xp0*cos(2*pi*r)
              fp(k,iyp)=xp0*sin(2*pi*r)
            endif
          enddo
!
        case ('random-sphere')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'in a sphere around (0,0,0) with radius=',rad_sphere
          if (rad_sphere==0) then
            call fatal_error('init_particles','random-sphere '// &
                'radius needs to be larger than zero')
          endif
          if (-rad_sphere+pos_sphere(1)<xyz0(1) .or. &
              +rad_sphere+pos_sphere(1)>xyz1(1) .or. &
              -rad_sphere+pos_sphere(2)<xyz0(2) .or. &
              +rad_sphere+pos_sphere(2)>xyz1(2) .or. &
              -rad_sphere+pos_sphere(3)<xyz0(3) .or. &
              +rad_sphere+pos_sphere(3)>xyz1(3)) then
            call fatal_error('init_particles','random-sphere '// &
                'sphere needs to fit in the box')
          endif
          if (lcartesian_coords) then
            do k=klow,khigh
              rp2=2.*rad_sphere**2
              do while (rp2>rad_sphere**2)
                call random_number_wrapper(r)
                fp(k,ixp)=(r-0.5)*2.*rad_sphere
                call random_number_wrapper(r)
                fp(k,iyp)=(r-0.5)*2.*rad_sphere
                call random_number_wrapper(r)
                fp(k,izp)=(r-0.5)*2.*rad_sphere
                rp2=fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2
              enddo
              fp(k,ixp)=fp(k,ixp)+pos_sphere(1)
              fp(k,iyp)=fp(k,iyp)+pos_sphere(2)
              fp(k,izp)=fp(k,izp)+pos_sphere(3)
            enddo
          else
            call fatal_error('init_particles','random-sphere '// &
                'only implemented for cartesian coordinates')
          endif
!
        case ('random-ellipsoid')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'in an ellipsoid around ', pos_ellipsoid, ' with ' // &
              'semi-principal axes a,b,c =',a_ellipsoid,b_ellipsoid,c_ellipsoid
          if ((a_ellipsoid==0) .or. (b_ellipsoid==0) .or. (c_ellipsoid==0)) then
            call fatal_error('init_particles','random-ellipsoid '// &
                'all semi-principal axes need to be larger than zero')
          endif
          if (-a_ellipsoid+pos_ellipsoid(1)<xyz0(1) .or. &
              +a_ellipsoid+pos_ellipsoid(1)>xyz1(1) .or. &
              -b_ellipsoid+pos_ellipsoid(2)<xyz0(2) .or. &
              +b_ellipsoid+pos_ellipsoid(2)>xyz1(2) .or. &
              -c_ellipsoid+pos_ellipsoid(3)<xyz0(3) .or. &
              +c_ellipsoid+pos_ellipsoid(3)>xyz1(3)) then
            call fatal_error('init_particles','random-ellipsoid '// &
                'ellipsoid needs to fit in the box')
          endif
          if (lcartesian_coords) then
            a_ell2=a_ellipsoid**2
            b_ell2=b_ellipsoid**2
            c_ell2=c_ellipsoid**2
            do k=klow,khigh
              rp2=2.
              do while (rp2>1.)
                call random_number_wrapper(r)
                fp(k,ixp)=(r-0.5)*2.*a_ellipsoid
                call random_number_wrapper(r)
                fp(k,iyp)=(r-0.5)*2.*b_ellipsoid
                call random_number_wrapper(r)
                fp(k,izp)=(r-0.5)*2.*c_ellipsoid
                rp2=fp(k,ixp)**2/a_ell2+fp(k,iyp)**2/b_ell2+fp(k,izp)**2/c_ell2
              enddo
              fp(k,ixp)=fp(k,ixp)+pos_ellipsoid(1)
              fp(k,iyp)=fp(k,iyp)+pos_ellipsoid(2)
              fp(k,izp)=fp(k,izp)+pos_ellipsoid(3)
            enddo
          else
            call fatal_error('init_particles','random-ellipsoid '// &
                'only implemented for cartesian coordinates')
          endif
!
        case ('random-line-x')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=klow,khigh
            if (nxgrid/=1) then
              call random_number_wrapper(r)
              fp(k,ixp)=r
            endif
          enddo
          if (nxgrid/=1) &
              fp(klow:khigh,ixp)=xyz0_par(1)+fp(klow:khigh,ixp)*Lxyz_par(1)
          fp(klow:khigh,iyp)=yp0
          fp(klow:khigh,izp)=zp0
!
        case ('random-line-y')
          if (lroot) print*, 'init_particles: Random particle positions'
          do k=klow,khigh
            if (nygrid/=1) then
              call random_number_wrapper(r)
              fp(k,iyp)=r
            endif
          enddo
          if (nygrid/=1) &
              fp(klow:khigh,iyp)=xyz0_par(2)+fp(klow:khigh,iyp)*Lxyz_par(2)
          fp(klow:khigh,ixp)=xp0
          fp(klow:khigh,izp)=zp0
!
        case ('random-hole')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'with inner hole'
          do k=klow:khigh
            rp2=-1.0
            do while (rp2<rp_int**2)
              rp2=0.0
              if (nxgrid/=1) then
                call random_number_wrapper(r)
                fp(k,ixp)=xyz0(1)+r*Lxyz(1)
                rp2=rp2+fp(k,ixp)**2
              endif
              if (nygrid/=1) then
                call random_number_wrapper(r)
                fp(k,iyp)=xyz0(2)+r*Lxyz(2)
                rp2=rp2+fp(k,iyp)**2
              endif
              if (nzgrid/=1) then
                call random_number_wrapper(r)
                fp(k,izp)=xyz0(3)+r*Lxyz(3)
                rp2=rp2+fp(k,izp)**2
              endif
            enddo
          enddo
!
        case ('random-box')
          if (lroot) print*, 'init_particles: Random particle positions '// &
              'within a box'
          do k=klow,khigh
            if (nxgrid/=1) then
              call random_number_wrapper(r)
              fp(k,ixp)=r
            endif
            if (nygrid/=1) then
              call random_number_wrapper(r)
              fp(k,iyp)=r
            endif
            if (nzgrid/=1) then
              call random_number_wrapper(r)
              fp(k,izp)=r
            endif
            if (lcylindrical_coords) then
              xx0=xp0+fp(k,ixp)*Lx0
              yy0=yp0+fp(k,iyp)*Ly0
              r2=xx0**2+yy0**2
              if (nxgrid/=1) fp(k,ixp)=sqrt(r2)
              if (nygrid/=1) fp(k,iyp)=atan(yy0/xx0)+pi*(xx0/abs(xx0)-1)*0.5
              if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
            else
              if (nxgrid/=1) fp(k,ixp)=xp0+fp(k,ixp)*Lx0
              if (nygrid/=1) fp(k,iyp)=yp0+fp(k,iyp)*Ly0
              if (nzgrid/=1) fp(k,izp)=zp0+fp(k,izp)*Lz0
            endif
          enddo
!
        case ('random-cylindrical','random-cyl')
!
          if (lroot) print*, 'init_particles: Random particle '//&
              'cylindrical positions with power-law = ', dustdensity_powerlaw
!
          do k=klow,khigh
!
! Start the particles obeying a power law
!
            if (lcylindrical_coords.or.lcartesian_coords) then
              tmp=2-dustdensity_powerlaw
            elseif (lspherical_coords) then 
              tmp=3-dustdensity_powerlaw
            else
              call fatal_error("init_particles",&
                   "The world is flat, and we never got here")
            endif
!
            if (lcartesian_coords) then 
              if (nprocx==1) then 
                rpar_int=rp_int
                rpar_ext=rp_ext
              else
                call fatal_error("init_particles",&
                     "random-cyl not yet ready for nprocx/=1 in Cartesian. Parallelize in y or z")
              endif
            else
              rpar_int = xyz0_loc(1)
              rpar_ext = xyz1_loc(1)
            endif
            call random_number_wrapper(rad_scl)
            rad_scl = rpar_int**tmp + rad_scl*(rpar_ext**tmp-rpar_int**tmp)
            rad = rad_scl**(1./tmp)
!
! Random in azimuth
!
            if (lcartesian_coords) then
              call random_number_wrapper(phi)
              phi = 2*pi*phi
              if (nxgrid/=1) fp(k,ixp)=rad*cos(phi)
              if (nygrid/=1) fp(k,iyp)=rad*sin(phi)
              if (nzgrid/=1) then
                call random_number_wrapper(r)
                fp(k,izp)=xyz0_par(3)+r*Lxyz_par(3)
              endif
            elseif (lcylindrical_coords) then
              if (nxgrid/=1) fp(k,ixp)=rad
              if (nygrid/=1) then
                call random_number_wrapper(phi)
                fp(k,iyp) = xyz0_par(2)+phi*Lxyz_par(2)
              endif
              if (nzgrid/=1) then
                call random_number_wrapper(r)
                fp(k,izp)=xyz0_par(3)+r*Lxyz_par(3)
              endif
            elseif (lspherical_coords) then
              if (nxgrid/=1) fp(k,ixp)=rad
              if (nygrid/=1) then
                call random_number_wrapper(tht)
                fp(k,iyp) = xyz0_par(2)+tht*Lxyz_par(2)
              endif
              if (nzgrid/=1) then
                call random_number_wrapper(phi)
                fp(k,izp) = xyz0_par(3)+phi*Lxyz_par(3)
              endif
            endif
!
          enddo
!
        case ('np-constant')
          if (lroot) print*, 'init_particles: Constant number density'
          k=klow
          k_loop: do while (.not. (k>khigh))
          do l=l1,l2
            do m=m1,m2
              do n=n1,n2
                if (nxgrid/=1) then
                  call random_number_wrapper(px)
                  fp(k,ixp)=x(l)+(px-0.5)*dx
                endif
                if (nygrid/=1) then
                  call random_number_wrapper(py)
                  fp(k,iyp)=y(m)+(py-0.5)*dy
                endif
                if (nzgrid/=1) then
                  call random_number_wrapper(pz)
                  fp(k,izp)=z(n)+(pz-0.5)*dz
                endif
                k=k+1
                if (k>khigh) exit k_loop
              enddo
            enddo
          enddo
          enddo k_loop
!
        case ('equidistant')
          if (lroot) print*, 'init_particles: Particles placed equidistantly'
          dim1=1.0/dimensionality
!
!  Number of particles per direction. Found by solving the equation system
!
!    npar_loc_x/npar_loc_y = Lx_loc/Ly_loc
!    npar_loc_x/npar_loc_z = Lx_loc/Lz_loc
!    npar_loc_y/npar_loc_z = Ly_loc/Lz_loc
!    npar_loc_x*npar_loc_y*npar_loc_z = npar_loc
!
!  Found it to be easier to separate in all possible dimensionalities.
!  For a missing direction i, set npar_loc_i=1 in the above equations and
!  ignore any equation that has Li_loc in it.
!
!  Initiate to avoid compiler warnings. Will be overwritten.
          npar_loc_x=1;npar_loc_y=1;npar_loc_z=1
!
          if (dimensionality==3) then
!  3-D
            npar_loc_x=(npar_loc*Lxyz_loc(1)**2/(Lxyz_loc(2)*Lxyz_loc(3)))**dim1
            npar_loc_y=(npar_loc*Lxyz_loc(2)**2/(Lxyz_loc(1)*Lxyz_loc(3)))**dim1
            npar_loc_z=(npar_loc*Lxyz_loc(3)**2/(Lxyz_loc(1)*Lxyz_loc(2)))**dim1
          elseif (dimensionality==2) then
!  2-D
            if (nxgrid==1) then
              npar_loc_x=1
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(3))**dim1
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(2))**dim1
            elseif (nygrid==1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(3))**dim1
              npar_loc_y=1
              npar_loc_z=(npar_loc*Lxyz_loc(3)/Lxyz_loc(1))**dim1
            elseif (nzgrid==1) then
              npar_loc_x=(npar_loc*Lxyz_loc(1)/Lxyz_loc(2))**dim1
              npar_loc_y=(npar_loc*Lxyz_loc(2)/Lxyz_loc(1))**dim1
              npar_loc_z=1
            endif
          elseif (dimensionality==1) then
!  1-D
            if (nxgrid/=1) then
              npar_loc_x=npar_loc
              npar_loc_y=1
              npar_loc_z=1
            elseif (nygrid/=1) then
              npar_loc_x=1
              npar_loc_y=npar_loc
              npar_loc_z=1
            elseif (nzgrid/=1) then
              npar_loc_x=1
              npar_loc_y=1
              npar_loc_z=npar_loc
            endif
          endif
!  Distance between particles.
          dx_par=Lxyz_loc(1)/npar_loc_x
          dy_par=Lxyz_loc(2)/npar_loc_y
          dz_par=Lxyz_loc(3)/npar_loc_z
!  Place first particle.
          fp(1,ixp) = x(l1) ; fp(1,iyp) = y(m1) ; fp(1,izp) = z(n1)
          if (nxgrid/=1) fp(1,ixp) = xyz0_loc(1)+dx_par/2
          if (nygrid/=1) fp(1,iyp) = xyz0_loc(2)+dy_par/2
          if (nzgrid/=1) fp(1,izp) = xyz0_loc(3)+dz_par/2
!  Place all other particles iteratively.
          if (dimensionality==3) then
!  3-D
            do k=2,npar_loc
              fp(k,ixp)=fp(k-1,ixp)+dx_par
              fp(k,iyp)=fp(k-1,iyp)
              fp(k,izp)=fp(k-1,izp)
              if (fp(k,ixp)>xyz1_loc(1)) then
                fp(k,ixp)=fp(1,ixp)
                fp(k,iyp)=fp(k,iyp)+dy_par
              endif
              if (fp(k,iyp)>xyz1_loc(2)) then
                fp(k,iyp)=fp(1,iyp)
                fp(k,izp)=fp(k,izp)+dz_par
              endif
            enddo
          elseif (dimensionality==2) then
!  2-D
            if (nxgrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)+dy_par
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,iyp)>xyz1_loc(2)) then
                  fp(k,iyp)=fp(1,iyp)
                  fp(k,izp)=fp(k,izp)+dz_par
                endif
              enddo
            elseif (nygrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,ixp)>xyz1_loc(1)) then
                  fp(k,ixp)=fp(1,ixp)
                  fp(k,izp)=fp(k,izp)+dz_par
                endif
              enddo
            elseif (nzgrid==1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
                if (fp(k,ixp)>xyz1_loc(1)) then
                  fp(k,ixp)=fp(1,ixp)
                  fp(k,iyp)=fp(k,iyp)+dy_par
                endif
              enddo
            endif
          elseif (dimensionality==1) then
!  1-D
            if (nxgrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)+dx_par
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)
              enddo
            elseif (nygrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)+dy_par
                fp(k,izp)=fp(k-1,izp)
              enddo
            elseif (nzgrid/=1) then
              do k=2,npar_loc
                fp(k,ixp)=fp(k-1,ixp)
                fp(k,iyp)=fp(k-1,iyp)
                fp(k,izp)=fp(k-1,izp)+dz_par
              enddo
            endif
          else
!  0-D
            fp(2:npar_loc,ixp)=fp(1,ixp)
            fp(2:npar_loc,iyp)=fp(1,iyp)
            fp(2:npar_loc,izp)=fp(1,izp)
          endif
          lequidistant=.true.
!
!  Shift particle locations slightly so that a mode appears.
!
        case ('shift')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) then
            call fatal_error('init_particles','must place particles equidistantly before shifting!')
          endif
          k2_xxp=kx_xxp**2+ky_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) then
            call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed!')
          endif
          do k=klow,khigh
            fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,iyp) = fp(k,iyp) - ky_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
            fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
                sin(kx_xxp*fp(k,ixp)+ky_xxp*fp(k,iyp)+kz_xxp*fp(k,izp))
          enddo
!
!  Shift to egg crate mode 2d, cos(x)cos(z)
!
        case ('cosxcosz')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) then
            call fatal_error('init_particles','must place particles equidistantly before shifting!')
          endif
          k2_xxp=kx_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) then
            call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed!')
          endif
          do k=klow,khigh
          fp(k,ixp) = fp(k,ixp) - kx_xxp/k2_xxp*amplxxp* &
              sin(kx_xxp*fp(k,ixp))*cos(kz_xxp*fp(k,izp))
          fp(k,izp) = fp(k,izp) - kz_xxp/k2_xxp*amplxxp* &
              sin(kx_xxp*fp(k,ixp))*cos(kz_xxp*fp(k,izp))
          enddo
!
!  Shift to egg crate mode 2d, sin(x)sin(z)
!
        case ('sinxsinz')
          if (lroot) print*, 'init_particles: shift particle positions'
          if (.not. lequidistant) then
            call fatal_error('init_particles','must place particles equidistantly before shifting!')
          endif
          k2_xxp=kx_xxp**2+kz_xxp**2
          if (k2_xxp==0.0) then
            call fatal_error('init_particles','kx_xxp=ky_xxp=kz_xxp=0.0 is not allowed!')
          endif
          do k=klow,khigh
          fp(k,ixp) = fp(k,ixp) + kx_xxp/k2_xxp*amplxxp* &
              cos(kx_xxp*fp(k,ixp))*sin(kz_xxp*fp(k,izp))
          fp(k,izp) = fp(k,izp) + kz_xxp/k2_xxp*amplxxp* &
              cos(kx_xxp*fp(k,ixp))*sin(kz_xxp*fp(k,izp))
          enddo
!
        case ('gaussian-z')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=klow,khigh
            do while (.true.)
              if (nxgrid/=1) then
                call random_number_wrapper(r)
                fp(k,ixp)=r
              endif
              if (nygrid/=1) then
                call random_number_wrapper(r)
                fp(k,iyp)=r
              endif
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              if (nprocz==2) then
                if (lfirst_proc_z) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
                if (llast_proc_z) fp(k,izp)=abs(zp0*sqrt(-2*alog(r))*cos(2*pi*  p))
              else
                fp(k,izp)= zp0*sqrt(-2*alog(r))*cos(2*pi*p)
              endif
              if ((fp(k,izp)>=xyz0(3)) .and. (fp(k,izp)<=xyz1(3))) exit
            enddo
          enddo
          if (nxgrid/=1) &
              fp(klow:khigh,ixp)=xyz0_par(1)+fp(klow:khigh,ixp)*Lxyz_par(1)
          if (nygrid/=1) &
              fp(klow:khigh,iyp)=xyz0_par(2)+fp(klow:khigh,iyp)*Lxyz_par(2)
!
        case ('gaussian-x')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=klow,khigh
            do while (.true.)
              if (nygrid/=1) then
                call random_number_wrapper(r)
                fp(k,iyp)=r
              endif
              if (nzgrid/=1) then
                call random_number_wrapper(r)
                fp(k,izp)=r
              endif
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              fp(k,ixp)= xp0*sqrt(-2*alog(r))*cos(2*pi*p)
              if ((fp(k,ixp)>=xyz0(1)).and.(fp(k,ixp)<=xyz1(1))) exit
            enddo
          enddo
          if (nygrid/=1) &
              fp(klow:khigh,iyp)=xyz0_par(2)+fp(klow:khigh,iyp)*Lxyz_par(2)
          if (nzgrid/=1) &
              fp(klow:khigh,izp)=xyz0_par(3)+fp(klow:khigh,izp)*Lxyz_par(3)
!
        case ('gaussian-z-pure')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=klow,khigh
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            if (nprocz==2) then
              if (lfirst_proc_z) fp(k,izp)=-abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
              if (llast_proc_z) fp(k,izp)=abs(zp0*sqrt(-2*alog(r))*cos(2*pi*p))
            else
              fp(k,izp)= zp0*sqrt(-2*alog(r))*cos(2*pi*p)
            endif
          enddo
!
        case ('gaussian-r')
          if (lroot) print*, 'init_particles: Gaussian particle positions'
          do k=klow,khigh
            call random_number_wrapper(r)
            call random_number_wrapper(p)
            call random_number_wrapper(q)
            fp(k,ixp)= xp0*sqrt(-2*alog(r))*cos(2*pi*p)*cos(2*pi*q)
            fp(k,iyp)= yp0*sqrt(-2*alog(r))*cos(2*pi*p)*sin(2*pi*q)
          enddo
!
        case ('hole')
!
          call map_nearest_grid(fp,ineargrid)
          call map_xxp_grid(f,fp,ineargrid)
          call sort_particles_imn(fp,ineargrid,ipar)
          do k=k1_imn(imn_array(m_hole+m1-1,n_hole+n1-1)), &
               k2_imn(imn_array(m_hole+m1-1,n_hole+n1-1))
            if (ineargrid(k,1)==l_hole+l1-1) then
              print*, k
              if (nxgrid/=0) fp(k,ixp)=fp(k,ixp)-dx
            endif
          enddo
!
        case ('streaming')
          call streaming(fp,f)
!
        case ('streaming_coldstart')
          call streaming_coldstart(fp,f)
!
        case ('constant-Ri')
          call constant_richardson(fp,f)
!
        case ('birthring')
          if (birthring_width>tini) then
            if (lgaussian_birthring) then
              do k=klow,khigh
                call normal_deviate(rr_tmp(k))
              enddo
            else
              call random_number_wrapper(rr_tmp(klow:khigh))
            endif
            rr_tmp(klow:khigh) = birthring_r+rr_tmp(klow:khigh)*birthring_width
          else
            rr_tmp(klow:khigh) = birthring_r
          endif
          call random_number_wrapper(az_tmp(klow:khigh))
          az_tmp(klow:khigh) = -pi+az_tmp(klow:khigh)*2.0*pi
          if (lcartesian_coords) then
            fp(klow:khigh,ixp) = rr_tmp(klow:khigh)*cos(az_tmp(klow:khigh))
            fp(klow:khigh,iyp) = rr_tmp(klow:khigh)*sin(az_tmp(klow:khigh))
            fp(klow:khigh,izp) = 0.0
          else
            fp(klow:khigh,ixp) = rr_tmp(klow:khigh)
            if (lcylindrical_coords) then
              fp(klow:khigh,iyp) = az_tmp(klow:khigh)
              fp(klow:khigh,izp) = 0.0
            elseif (lspherical_coords) then
              fp(klow:khigh,iyp) = pi/2.0
              fp(klow:khigh,izp) = az_tmp(klow:khigh)
            endif
          endif
          if (lroot .and. nzgrid/=0) print*,"Warning, birthring only implemented for 2D"
!
        case default
          call fatal_error('init_particles','Unknown value initxxp="'//trim(initxxp(j))//'"')
!
        endselect
!
      enddo ! do j=1,ninit
!
!  Interface for user's own initial condition for position
!
      if (linitial_condition) call initial_condition_xxp(f,fp)
!
!  Particles are not allowed to be present in non-existing dimensions.
!  This would give huge problems with interpolation later.
!
      if (nxgrid==1) fp(klow:khigh,ixp)=x(nghost+1)
      if (nygrid==1) fp(klow:khigh,iyp)=y(nghost+1)
      if (nzgrid==1) fp(klow:khigh,izp)=z(nghost+1)
!
      if (init_repeat/=0) call repeated_init(fp,init_repeat)
!
!  Initial particle velocity.
!
      do j=1,ninit
!
        select case (initvvp(j))
!
        case ('nothing')
          if (lroot.and.j==1) print*, 'init_particles: No particle velocity set'
!
        case ('zero')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(klow:khigh,ivpx:ivpz)=0.0
!
        case ('zero-shear')
          if (lroot) print*, 'init_particles: Zero particle velocity'
          fp(klow:khigh,ivpy)=-Sshear*fp(klow:khigh,ixp)
          fp(klow:khigh,ivpx)=0.0
          fp(klow:khigh,ivpz)=0.0
!
        case ('constant')
          if (lroot) print*, 'init_particles: Constant particle velocity'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          if (lcylindrical_coords) then
            fp(klow:khigh,ivpx)=vpx0*cos(fp(k,iyp))+vpy0*sin(fp(k,iyp))
            fp(klow:khigh,ivpy)=vpy0*cos(fp(k,iyp))-vpx0*sin(fp(k,iyp))
            fp(klow:khigh,ivpz)=vpz0
          else
            fp(klow:khigh,ivpx)=vpx0
            fp(klow:khigh,ivpy)=vpy0
            fp(klow:khigh,ivpz)=vpz0
          endif
!
        case ('constant-1')
          if (lroot) &
              print*, 'init_particles: Particle 1 velocity vx,vy,vz=', &
              vpx1, vpy1, vpz1
          do k=klow,khigh
            if (ipar(k)==1) then
              fp(k,ivpx)=vpx1
              fp(k,ivpy)=vpy1
              fp(k,ivpz)=vpz1
            endif
          enddo
!
        case ('constant-2')
          if (lroot) &
              print*, 'init_particles: Particle 2 velocity vx,vy,vz=', &
              vpx2, vpy2, vpz2
          do k=klow,khigh
            if (ipar(k)==2) then
              fp(k,ivpx)=vpx2
              fp(k,ivpy)=vpy2
              fp(k,ivpz)=vpz2
            endif
          enddo
!
        case ('constant-3')
          if (lroot) &
              print*, 'init_particles: Particle 3 velocity vx,vy,vz=', &
              vpx3, vpy3, vpz3
          do k=klow,khigh
            if (ipar(k)==3) then
              fp(k,ivpx)=vpx3
              fp(k,ivpy)=vpy3
              fp(k,ivpz)=vpz3
            endif
          enddo
!
        case ('sinwave-phase')
          if (lroot) print*, 'init_particles: sinwave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=klow,khigh
            fp(k,ivpx)=fp(k,ivpx)+vpx0*sin(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*sin(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*sin(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('coswave-phase')
          if (lroot) print*, 'init_particles: coswave-phase'
          if (lroot) &
              print*, 'init_particles: vpx0, vpy0, vpz0=', vpx0, vpy0, vpz0
          do k=klow,khigh
            fp(k,ivpx)=fp(k,ivpx)+vpx0*cos(kx_vpx*fp(k,ixp)+ky_vpx*fp(k,iyp)+kz_vpx*fp(k,izp)+phase_vpx)
            fp(k,ivpy)=fp(k,ivpy)+vpy0*cos(kx_vpy*fp(k,ixp)+ky_vpy*fp(k,iyp)+kz_vpy*fp(k,izp)+phase_vpy)
            fp(k,ivpz)=fp(k,ivpz)+vpz0*cos(kx_vpz*fp(k,ixp)+ky_vpz*fp(k,iyp)+kz_vpz*fp(k,izp)+phase_vpz)
          enddo
!
        case ('random')
          if (lroot) print*, 'init_particles: Random particle velocities; '// &
              'delta_vp0=', delta_vp0
          do k=klow,khigh
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-x')
          if (lroot) print*, 'init_particles: Random particle x-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=klow,khigh
            call random_number_wrapper(r)
            fp(k,ivpx) = fp(k,ivpx) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-y')
          if (lroot) print*, 'init_particles: Random particle y-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=klow,khigh
            call random_number_wrapper(r)
            fp(k,ivpy) = fp(k,ivpy) + delta_vp0*(2*r-1)
          enddo
!
        case ('random-z')
          if (lroot) print*, 'init_particles: Random particle z-velocity; '// &
              'delta_vp0=', delta_vp0
          do k=klow,khigh
            call random_number_wrapper(r)
            fp(k,ivpz) = fp(k,ivpz) + delta_vp0*(2*r-1)
          enddo
!
        case ('average-to-zero')
          call mpireduce_sum(sum(fp(1:npar_loc,ivpx)),vpx_sum)
          call mpireduce_sum(sum(fp(1:npar_loc,ivpy)),vpy_sum)
          call mpireduce_sum(sum(fp(1:npar_loc,ivpz)),vpz_sum)
          call mpibcast_real(vpx_sum)
          call mpibcast_real(vpy_sum)
          call mpibcast_real(vpz_sum)
          fp(1:npar_loc,ivpx)=fp(1:npar_loc,ivpx)-vpx_sum/npar
          fp(1:npar_loc,ivpy)=fp(1:npar_loc,ivpy)-vpy_sum/npar
          fp(1:npar_loc,ivpz)=fp(1:npar_loc,ivpz)-vpz_sum/npar
!
        case ('follow-gas')
          if (lroot) &
              print*, 'init_particles: Particle velocity equal to gas velocity'
          do k=klow,khigh
            call interpolate_linear(f,iux,iuz,fp(k,ixp:izp),uup, &
                ineargrid(k,:),0,0)
            fp(k,ivpx:ivpz) = uup
          enddo
!
        case ('jeans-wave-dustpar-x')
        ! assumes rhs_poisson_const=1 !
          do k=klow,khigh
            fp(k,ivpx) = fp(k,ivpx) - amplxxp* &
                (sqrt(1+4*1.0*1.0*tausp**2)-1)/ &
                (2*kx_xxp*1.0*tausp)*sin(kx_xxp*(fp(k,ixp)))
          enddo
!
        case ('dragforce_equilibrium','dragforce-equilibrium')
!
!  Equilibrium between drag forces on dust and gas and other forces
!  (from Nakagawa, Sekiya, & Hayashi 1986).
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium'
            print*, 'init_particles: beta_glnrho_global=', beta_glnrho_global
          endif
          cs=sqrt(cs20)
!
          if (ldragforce_equi_global_eps) eps = eps_dtog
!
          if (ldragforce_equi_noback) eps=0.0
!
          if (lroot) print*, 'init_particles: average dust-to-gas ratio=', eps
!  Set gas velocity field.
          if (lhydro) then
            do l=l1,l2; do m=m1,m2; do n=n1,n2
!  Take either global or local dust-to-gas ratio.
              if (.not. ldragforce_equi_global_eps) eps = f(l,m,n,irhop) / get_gas_density(f,l,m,n)
!
              f(l,m,n,iux) = f(l,m,n,iux) - &
                  beta_glnrho_global(1)*eps*Omega*tausp/ &
                  ((1.0+eps)**2+(Omega*tausp)**2)*cs
              f(l,m,n,iuy) = f(l,m,n,iuy) + &
                  beta_glnrho_global(1)*(1+eps+(Omega*tausp)**2)/ &
                  (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
!
            enddo; enddo; enddo
          endif
!  Set particle velocity field.
          do k=klow,khigh
!  Take either global or local dust-to-gas ratio.
            if (ldragforce_equi_noback) then
              eps=0.0
            else
              if (.not. ldragforce_equi_global_eps) then
                ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
                eps = f(ix0,iy0,iz0,irhop) / get_gas_density(f,ix0,iy0,iz0)
              endif
            endif
!
            fp(k,ivpx) = fp(k,ivpx) + &
                beta_glnrho_global(1)*Omega*tausp/ &
                ((1.0+eps)**2+(Omega*tausp)**2)*cs
            fp(k,ivpy) = fp(k,ivpy) + &
                beta_glnrho_global(1)*(1+eps)/ &
                (2*((1.0+eps)**2+(Omega*tausp)**2))*cs
!
          enddo
!
        case ('dragforce_equi_nohydro')
!
          do k=klow,khigh
            fp(k,ivpx) = fp(k,ivpx) - 2*Deltauy_gas_friction* &
                1/(1.0/(Omega*tausp)+Omega*tausp)
            fp(k,ivpy) = fp(k,ivpy) - Deltauy_gas_friction* &
                1/(1.0+(Omega*tausp)**2)
          enddo
!
        case ('dragforce_equi_dust')
!
!  Equilibrium between drag force and Coriolis force on the dust.
!
          if (lroot) then
            print*, 'init_particles: drag equilibrium dust'
            print*, 'init_particles: beta_dPdr_dust=', beta_dPdr_dust
          endif
!  Set particle velocity field.
          cs=sqrt(cs20)
          do k=klow,khigh
            fp(k,ivpx) = fp(k,ivpx) + &
                1/(Omega*tausp+1/(Omega*tausp))*beta_dPdr_dust*cs
            fp(k,ivpy) = fp(k,ivpy) - &
                1/(1.0+1/(Omega*tausp)**2)*beta_dPdr_dust/2*cs
          enddo
!
       case ('Keplerian','keplerian')
!
!  Keplerian velocity based on gravr.
!
          if (lroot) then
            print*, 'init_particles: Keplerian velocity'
            if (lshear) call fatal_error("init_particles",&
                 "Keplerian initial condition is for global disks, not shearing boxes")
          endif
          do k=klow,khigh
            if (lcartesian_coords) then
              rad=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) = -OO*fp(k,iyp)
              fp(k,ivpy) =  OO*fp(k,ixp)
              fp(k,ivpz) =  0.0
            elseif (lcylindrical_coords) then
              rad=fp(k,ixp)
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) =  0.0
              fp(k,ivpy) =  OO*rad
              fp(k,ivpz) =  0.0
            elseif (lspherical_coords) then 
              rad=fp(k,ixp)*sin(fp(k,iyp))
              OO=sqrt(gravr)*rad**(-1.5)
              fp(k,ivpx) =  0.0
              fp(k,ivpy) =  0.0
              fp(k,ivpz) =  OO*rad
            endif
          enddo
!
!  Explosion.
!
       case ('explosion')
         do k=klow,khigh
           rad=sqrt(fp(k,ixp)**2+fp(k,iyp)**2+fp(k,izp)**2)
           fp(k,ivpx) = delta_vp0*fp(k,ixp)/rp_ext
           fp(k,ivpy) = delta_vp0*fp(k,iyp)/rp_ext
           fp(k,ivpz) = delta_vp0*fp(k,izp)/rp_ext
         enddo
!
!
        case default
          call fatal_error('init_particles','Unknown value initvvp="'//trim(initvvp(j))//'"')
!
        endselect
!
      enddo ! do j=1,ninit

    endsubroutine init_particle_positions_velocities
!***********************************************************************
    subroutine input_particles(filename,fp,ipar)
!
!  Read snapshot file with particle data.
!
!  29-dec-04/anders: adapted from input
!
      use Mpicomm, only: mpireduce_max_scl_int
!
      real, dimension (mpar_loc,mparray) :: fp
      character (len=*) :: filename
      integer, dimension (mpar_loc) :: ipar
!
      intent (in) :: filename
      intent (out) :: fp,ipar
!
      open(1,FILE=filename,FORM='unformatted')
!
!  First read the number of particles present at the processor and the index
!  numbers of the particles.
!
        read(1) npar_loc
        if (npar_loc/=0) read(1) ipar(1:npar_loc)
!
!  Then read particle data.
!
        if (npar_loc/=0) read(1) fp(1:npar_loc,:)
!
!  Read snapshot time.
!
!        read(1) t
!
        if (ip<=8) print*, 'input_particles: read ', filename
!
      close(1)
!
!  If we are inserting particles contiuously during the run root must
!  know the total number of particles in the simulation.
!
      if (npar_loc/=0) then
        call mpireduce_max_scl_int(maxval(ipar(1:npar_loc)),npar_total)
      else
        call mpireduce_max_scl_int(npar_loc,npar_total)
      endif
!
    endsubroutine input_particles
!***********************************************************************
    subroutine output_particles(filename,fp,ipar)
!
!  Write snapshot file with particle data.
!
!  29-dec-04/anders: adapted from output
!
      character(len=*) :: filename
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc) :: ipar
      real :: t_sp   ! t in single precision for backwards compatibility
!
      intent (in) :: filename, ipar
!
      t_sp = t
      if (ip<=8.and.lroot) print*,'output_particles: writing snapshot file '// &
          filename
!
      open(lun_output,FILE=filename,FORM='unformatted')
!
!  First write the number of particles present at the processor and the index
!  numbers of the particles.
!
        write(lun_output) npar_loc
        if (npar_loc/=0) write(lun_output) ipar(1:npar_loc)
!
!  Then write particle data.
!
        if (npar_loc/=0) write(lun_output) fp(1:npar_loc,:)
!
!  Write time and grid parameters.
!
        write(lun_output) t_sp, x, y, z, dx, dy, dz
!
      close(lun_output)
!
    endsubroutine output_particles
!***********************************************************************
    subroutine boundconds_particles(fp,ipar,dfp,linsert)
!
!  Global boundary conditions for particles.
!
!  30-dec-04/anders: coded
!
      use Mpicomm
      use General, only: random_number_wrapper
      use Particles_mpicomm
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
      logical, optional :: linsert
!
      real :: xold, yold, rad, r1old, OO, tmp
      integer :: k, ik, k1, k2
      character (len=2*bclen+1) :: boundx, boundy, boundz
      logical :: lnbody
!
      intent (inout) :: fp, ipar, dfp
!
!  Reversing direction of looping is needed for removal.
!
      if ((bcpx=='rmv').or.(bcpy=='rmv').or.(bcpz=='rmv')) then
        k1=npar_loc; k2=1; ik=-1
      else
        k1=1; k2=npar_loc; ik=1
      endif
!
      do k=k1,k2,ik
!
!  Check if we are dealing with a dust or a massive particle
!
        lnbody=(lparticles_nbody.and.any(ipar(k)==ipar_nbody))
        if (.not.lnbody) then
          boundx=bcpx ; boundy=bcpy ; boundz=bcpz
        else
          boundx=bcspx; boundy=bcspy; boundz=bcspz
        endif
!
!  Calculate rad for cylinder-in-a-box calculations
!
        if (lcylinder_in_a_box) then
          rad = sqrt(fp(k,ixp)**2 + fp(k,iyp)**2)
        endif
!
!  Cartesian boundaries: Boundary condition in the x-direction. The physical
!  domain is in the interval
!
!    x \in [x0,x1[
!    y \in [y0,y1[
!    z \in [z0,z1[
!
        if (nxgrid/=1) then
          if (boundx=='p') then
!  xp < x0
            if (fp(k,ixp)< xyz0(1)) then
              fp(k,ixp)=fp(k,ixp)+Lxyz(1)
              if (lshear.and.nygrid/=1) then
                fp(k,iyp)=fp(k,iyp)-deltay
!
!  Keep track of energy gain or lose from shearing boundaries.
!
                if (energy_gain_shear_bcs/=impossible) &
                    energy_gain_shear_bcs=energy_gain_shear_bcs + &
                    (1.0/6.0)*Sshear**2*Omega**2*Lx**2 + &
                    Sshear*(fp(k,ivpy)+Sshear*fp(k,ixp))*Lx - &
                    (4.0/3.0)*Sshear**2*fp(k,ixp)*Lx
              endif
!  Particle position must never need more than one addition of Lx to get back
!  in the box. Often a NaN or Inf in the particle position will show up as a
!  problem here.
              if (fp(k,ixp)< xyz0(1)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
              if (ipyy/=0) fp(k,ipxx)=fp(k,ipxx)+Lxyz(1)
            endif
!  xp > x1
            if (fp(k,ixp)>=xyz1(1)) then
              fp(k,ixp)=fp(k,ixp)-Lxyz(1)
              if (lshear.and.nygrid/=1) then
                fp(k,iyp)=fp(k,iyp)+deltay
!
!  Keep track of energy gain or lose from shearing boundaries.
!
                if (energy_gain_shear_bcs/=impossible) &
                    energy_gain_shear_bcs=energy_gain_shear_bcs + &
                    (1.0/6.0)*Sshear**2*Omega**2*Lx**2 - &
                    Sshear*(fp(k,ivpy)+Sshear*fp(k,ixp))*Lx + &
                    (4.0/3.0)*Sshear**2*fp(k,ixp)*Lx
              endif
              if (fp(k,ixp)>=xyz1(1)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
              if (ipxx/=0) fp(k,ipxx)=fp(k,ipxx)-Lxyz(1)
            endif
          elseif (boundx=='out') then
!
!  Do nothing. A massive particle, can be out of the box. A star, for example,
!  in a cylindrical simulation
!
          elseif (boundx=='flk') then
!
!  Flush-Keplerian - flush the particle to the outer boundary with keplerian
!  speed
!
            if (lcylindrical_coords) then
              if ((fp(k,ixp)< rp_int).or.(fp(k,ixp)>= rp_ext)) then
!   Flush to outer boundary
                fp(k,ixp)  = rp_ext
!   Random new azimuthal y position
                call random_number_wrapper(fp(k,iyp))
                fp(k,iyp)=xyz0_loc(2)+fp(k,iyp)*Lxyz_loc(2)
!   Zero radial, Keplerian azimuthal, velocities
                fp(k,ivpx) = 0.
!   Keplerian azimuthal velocity
                fp(k,ivpy) = fp(k,ixp)**(-1.5)
              endif
!
            elseif (lcartesian_coords) then
!
! The Cartesian case has the option cylinder_in_a_box, sphere_in_a_box
! and nothing (assumed to be the common shearing box. Only the cylinder
! is considered. The code will break otherwise
!
              if (lcylinder_in_a_box) then
!
                if (boundy/='out') then
                  call fatal_error_local("boundconds_particles",&
                       "The radial boundary already does it all"//&
                       "so the y-boundary has to be set to 'out'")
                endif
!
                if ((rad<=rp_int).or.(rad>rp_ext)) then
!  Flush to outer boundary
                  xold=fp(k,ixp); yold=fp(k,iyp); r1old = 1./max(rad,tini)
                  fp(k,ixp) = rp_ext*xold*r1old ! r*cos(theta)
                  fp(k,iyp) = rp_ext*yold*r1old ! r*sin(theta)
!  Set particle velocity to local Keplerian speed.
                  OO=rp_ext**(-1.5)
                  fp(k,ivpx) = -OO*fp(k,iyp)
                  fp(k,ivpy) =  OO*fp(k,ixp)
                endif
              else
                call fatal_error_local('boundconds_particles',&
                     'no clue how to do flush-keplerian in a cartesian box')
              endif
            elseif (lspherical_coords) then
              call fatal_error_local('boundconds_particles',&
                   'flush-keplerian not ready for spherical coords')
            endif
!
          elseif (boundx=='rmv') then
!
!  Remove the particle from the simulation
!
            if (lspherical_coords.or.lcylindrical_coords) then
              if ((fp(k,ixp)< rp_int).or.(fp(k,ixp)>= rp_ext)) then
                if (present(dfp)) then
                  call remove_particle(fp,ipar,k,dfp)
                else
                  call remove_particle(fp,ipar,k)
                endif
              endif
            elseif (lcartesian_coords) then
              if (lcylinder_in_a_box) then
                if ((rad< rp_int).or.(rad>= rp_ext)) then
                  if (present(dfp)) then
                    call remove_particle(fp,ipar,k,dfp)
                  else
                    call remove_particle(fp,ipar,k)
                  endif
                endif
              else
                if (fp(k,ixp)<=xyz0(1) .or. fp(k,ixp)>=xyz1(1)) then
                  if (present(dfp)) then
                    call remove_particle(fp,ipar,k,dfp)
                  else
                    call remove_particle(fp,ipar,k)
                  endif
                endif
              endif
            endif
          elseif (boundx=='hw') then
!
!  Hard wall boundary condition
!
            if (fp(k,ixp)<=xyz0(1)) then
              fp(k,ixp)=2*xyz0(1)-fp(k,ixp)
              fp(k,ivpx)=-fp(k,ivpx)
            elseif (fp(k,ixp)>=xyz1(1)) then
              fp(k,ixp)=2*xyz1(1)-fp(k,ixp)
              fp(k,ivpx)=-fp(k,ivpx)
            endif
          elseif (boundx=='meq') then
            if ((fp(k,ixp).lt.rp_int) .or. (fp(k,ixp).gt.rp_ext)) then
            if (lcylindrical_coords) then
              tmp=2-dustdensity_powerlaw
            elseif (lspherical_coords) then
              tmp=3-dustdensity_powerlaw
            endif
              call random_number_wrapper(fp(k,ixp))
              fp(k,ixp) = (rp_int**tmp + fp(k,ixp)*(rp_ext**tmp-rp_int**tmp))**(1.0/tmp)
              if (lcylindrical_coords) then
                if (nygrid/=1) then
                  call random_number_wrapper(fp(k,iyp))
                  fp(k,iyp) = xyz0(2)+fp(k,iyp)*Lxyz(2)
                endif
                if (nzgrid/=1) then
                  call random_number_wrapper(fp(k,izp))
                  fp(k,izp)=xyz0(3)+fp(k,izp)*Lxyz(3)
                endif
              elseif (lspherical_coords) then
                if (nygrid/=1) then
                  call random_number_wrapper(fp(k,iyp))
                  fp(k,iyp) = xyz0(2)+fp(k,iyp)*Lxyz(2)
                endif
                if (nzgrid/=1) then
                  call random_number_wrapper(fp(k,izp))
                  fp(k,izp) = xyz0(3)+fp(k,izp)*Lxyz(3)
                endif
              endif
            endif
          else
            print*, 'boundconds_particles: No such boundary condition =', boundx
            call stop_it('boundconds_particles')
          endif
        endif
!
!  Boundary condition in the y-direction.
!
        if (nygrid/=1) then
          if (boundy=='p') then
!  yp < y0
            if (fp(k,iyp)< xyz0(2)) then
              fp(k,iyp)=fp(k,iyp)+Lxyz(2)
              if (fp(k,iyp)< xyz0(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
              if (ipyy/=0) fp(k,ipyy)=fp(k,ipyy)+Lxyz(2)
            endif
!  yp > y1
            if (fp(k,iyp)>=xyz1(2)) then
              fp(k,iyp)=fp(k,iyp)-Lxyz(2)
              if (fp(k,iyp)>=xyz1(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
              if (ipyy/=0) fp(k,ipyy)=fp(k,ipyy)-Lxyz(2)
            endif
          elseif (boundy=='out') then
            ! massive particles can be out of the box
            ! the star, for example, in a cylindrical simulation
          elseif (boundy=='rmv') then
            if (lcartesian_coords) then
              if (fp(k,iyp)<=xyz0(2) .or. fp(k,iyp)>=xyz1(2)) then
                if (present(dfp)) then
                  call remove_particle(fp,ipar,k,dfp)
                else
                  call remove_particle(fp,ipar,k)
                endif
              endif
            endif
          elseif (boundy=='hw') then
!
!  Hard wall boundary condition
!
            if (fp(k,iyp)<=xyz0(2)) then
              fp(k,iyp)=2*xyz0(2)-fp(k,iyp)
              fp(k,ivpy)=-fp(k,ivpy)
            elseif (fp(k,iyp)>=xyz1(2)) then
              fp(k,iyp)=2*xyz1(2)-fp(k,iyp)
              fp(k,ivpy)=-fp(k,ivpy)
            endif
          else
            print*, 'boundconds_particles: No such boundary condition =', boundy
            call stop_it('boundconds_particles')
          endif
        endif
!
!  Boundary condition in the z-direction.
!
        if (nzgrid/=1) then
          if (boundz=='p') then
!  zp < z0
            if (fp(k,izp)< xyz0(3)) then
              fp(k,izp)=fp(k,izp)+Lxyz(3)
              if (fp(k,izp)< xyz0(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
              if (ipzz/=0) fp(k,ipzz)=fp(k,ipzz)+Lxyz(3)
            endif
!  zp > z1
            if (fp(k,izp)>=xyz1(3)) then
              fp(k,izp)=fp(k,izp)-Lxyz(3)
              if (fp(k,izp)>=xyz1(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                     ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
              if (ipzz/=0) fp(k,ipzz)=fp(k,ipzz)-Lxyz(3)
            endif
          elseif (boundz=='out') then
            ! massive particles can be out of the box
            ! the star, for example, in a cylindrical simulation
          elseif (boundz=='rmv') then
            if (lcartesian_coords) then
              if (fp(k,izp)<=xyz0(3) .or. fp(k,izp)>=xyz1(3)) then
                if (present(dfp)) then
                  call remove_particle(fp,ipar,k,dfp)
                else
                  call remove_particle(fp,ipar,k)
                endif
              endif
            endif
          elseif (boundz=='ref') then
            if ((fp(k,izp)< xyz0(3)) .or. (fp(k,izp)>=xyz1(3))) then
              fp(k,ivpz)=-fp(k,ivpz)
              if (fp(k,izp)< xyz0(3)) then
                fp(k,izp)=2*xyz0(3)-fp(k,izp)
              endif
              if (fp(k,izp)>=xyz1(3)) then
                fp(k,izp)=2*xyz1(3)-fp(k,izp)
              endif
            endif
          elseif (boundz=='hw') then
!
!  Hard wall boundary condition
!
            if (fp(k,izp)<=xyz0(3)) then
              fp(k,izp)=2*xyz0(3)-fp(k,izp)
              fp(k,ivpz)=-fp(k,ivpz)
            elseif (fp(k,izp)>=xyz1(3)) then
              fp(k,izp)=2*xyz1(3)-fp(k,izp)
              fp(k,ivpz)=-fp(k,ivpz)
            endif
          elseif (boundz=='inc') then
!
!  Move particle from the boundary to the center of the box
!
            if (fp(k,izp)<=xyz0(3) .or. fp(k,izp)>=xyz1(3)) then
              fp(k,izp)=xyz0(3)+(xyz1(3)-xyz0(3))/2.
              fp(k,ivpx:ivpz)=0.
            endif
          else
            print*, 'boundconds_particles: No such boundary condition=', boundz
            call stop_it('boundconds_particles')
          endif
        endif
      enddo
!
!  Redistribute particles among processors (internal boundary conditions).
!
      if (lmpicomm) then
        if (present(dfp)) then
          if (present(linsert).or.(boundx=='meq')) then
            call migrate_particles(fp,ipar,dfp,linsert=.true.)
          else
            call migrate_particles(fp,ipar,dfp)
          endif
        else
          if (present(linsert).or.(boundx=='meq')) then
            call migrate_particles(fp,ipar,linsert=.true.)
          else
            call migrate_particles(fp,ipar)
          endif
        endif
      endif
!
    endsubroutine boundconds_particles
!***********************************************************************
    subroutine sum_par_name(a,iname,lsqrt)
!
!  Successively calculate sum of a, which is supplied at each call.
!  Works for particle diagnostics. The number of particles is stored as
!  a weight used later for normalisation of sum.
!
!  TODO: All processors must enter this subroutine in order to set the
!        diagnostics type right, even processors with no particles.
!        One can do this by calling sum_par_name(fp(1:npar_loc,ixp),...),
!        which sends an array of zero size. This is a bug and needs to be
!        fixed.
!
!  02-jan-05/anders: adapted from sum_mn_name
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      integer, dimension(mname), save :: icount=0
!
      if (iname/=0) then
!
        if (icount(iname)==0) then
          fname(iname)=0.0
          fweight(iname)=0.0
        endif
!
        fname(iname)  =fname(iname)  +sum(a)
        fweight(iname)=fweight(iname)+size(a)
!
!  Set corresponding entry in itype_name
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt_par
        else
          itype_name(iname)=ilabel_sum_par
        endif
!
!  Reset sum when npar_loc particles have been considered.
!
        icount(iname)=icount(iname)+size(a)
        if (icount(iname)==npar_loc) then
          icount(iname)=0
        elseif (icount(iname)>=npar_loc) then
          print*, 'sum_par_name: Too many particles entered this sub.'
          print*, 'sum_par_name: Can only do statistics on npar_loc particles!'
          call fatal_error('sum_par_name','')
        endif
!
      endif
!
    endsubroutine sum_par_name
!***********************************************************************
    subroutine sum_par_name_nw(a,iname,lsqrt)
!
!  successively calculate sum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  22-aug-05/anders: adapted from sum_par_name
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      integer, save :: icount=0
!
      if (iname/=0) then
!
        if (icount==0) fname(iname)=0
!
        fname(iname)=fname(iname)+sum(a)
!
!  Set corresponding entry in itype_name
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt
        else
          itype_name(iname)=ilabel_sum
        endif
!
        icount=icount+size(a)
        if (icount==nw) then
          icount=0
        elseif (icount>=nw) then
          print*, 'sum_par_name_nw: Too many grid points entered this sub.'
          print*, 'sum_par_name_nw: Can only do statistics on nw grid points!'
          call fatal_error('sum_par_name_nw','')
        endif
!
      endif
!
    endsubroutine sum_par_name_nw
!***********************************************************************
    subroutine max_par_name(a,iname,lneg)
!
!  Successively calculate maximum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  28-nov-05/anders: adapted from max_par_name
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lneg
!
      if (iname/=0) then
!
        fname(iname)=0.
        fname(iname)=fname(iname)+maxval(a)
!
!  Set corresponding entry in itype_name.
!
        if (present(lneg)) then
          itype_name(iname)=ilabel_max_neg
        else
          itype_name(iname)=ilabel_max
        endif
!
      endif
!
    endsubroutine max_par_name
!***********************************************************************
    subroutine integrate_par_name(a,iname)
!
!  Calculate integral of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  29-nov-05/anders: adapted from sum_par_name
!
      real, dimension (:) :: a
      integer :: iname
!
      integer, save :: icount=0
!
      if (iname/=0) then
!
        if (icount==0) fname(iname)=0
!
        fname(iname)=fname(iname)+sum(a)
!
!  Set corresponding entry in itype_name.
!
        itype_name(iname)=ilabel_integrate
!
!  Reset sum when npar_loc particles have been considered.
!
        icount=icount+size(a)
        if (icount==npar_loc) then
          icount=0
        elseif (icount>=npar_loc) then
          print*, 'integral_par_name: Too many particles entered this sub.'
          print*, 'integral_par_name: Can only do statistics on npar_loc particles!'
          call fatal_error('integral_par_name','')
        endif
!
      endif
!
    endsubroutine integrate_par_name
!***********************************************************************
    subroutine sharpen_tsc_density(f)
!
!  Sharpen density amplitudes (experimental).
!
!   9-nov-06/anders: coded
!
      use Fourier
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (:,:,:), allocatable :: a1, b1
      integer :: ikx, iky, ikz, stat
      real :: k2, k4
!
!  Allocate memory for large arrays.
!
      allocate(a1(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('sharpen_tsc_density', &
          'Could not allocate memory for a1')
      allocate(b1(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('sharpen_tsc_density', &
          'Could not allocate memory for b1')
!
      a1=f(l1:l2,m1:m2,n1:n2,irhop)
      b1=0.0
      if (lshear) then
        call fourier_transform_shear(a1,b1)
      else
        call fourier_transform(a1,b1)
      endif
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
        k2 = kx_fft(ikx)**2 + ky_fft(iky+ipy*ny)**2 + kz_fft(ikz+ipz*nz)**2
        k4 = kx_fft(ikx)**4 + ky_fft(iky+ipy*ny)**4 + kz_fft(ikz+ipz*nz)**4
        a1(ikx,iky,ikz)=a1(ikx,iky,ikz)/(1-dx**2*k2/8+13*dx**4*k4/1920)
        b1(ikx,iky,ikz)=b1(ikx,iky,ikz)/(1-dx**2*k2/8+13*dx**4*k4/1920)
      enddo; enddo; enddo
      if (lshear) then
        call fourier_transform_shear(a1,b1,linv=.true.)
      else
        call fourier_transform(a1,b1,linv=.true.)
      endif
      f(l1:l2,m1:m2,n1:n2,irhop)=a1
!
!  Deallocate arrays.
!
      if (allocated(a1)) deallocate(a1)
      if (allocated(b1)) deallocate(b1)
!
    endsubroutine sharpen_tsc_density
!***********************************************************************
    subroutine get_particles_interdistance(xx1,xx2,vector,distance,lsquare)
!
!  The name of the subroutine is pretty self-explanatory.
!
!  14-mar-08/wlad: moved here from particles_nbody
!
      real,dimension(3) :: xx1,xx2,evr
      real :: e1,e2,e3,e10,e20,e30
      real,dimension(3),optional :: vector
      real,optional :: distance
      logical,optional :: lsquare
!
      e1=xx1(1);e10=xx2(1)
      e2=xx1(2);e20=xx2(2)
      e3=xx1(3);e30=xx2(3)
!
!  Get the distances in each ortogonal component.
!  These are NOT (x,y,z) for all.
!  For cartesian it is (x,y,z), for cylindrical (s,phi,z)
!  for spherical (r,theta,phi)
!
      if (lcartesian_coords) then
        evr(1) = e1 - e10
        evr(2) = e2 - e20
        evr(3) = e3 - e30
      elseif (lcylindrical_coords) then
        evr(1) = e1 - e10*cos(e2-e20)
        evr(2) = e10*sin(e2-e20)
        evr(3) = e3 - e30
      elseif (lspherical_coords) then
        evr(1) = e1 - e10*sin(e2)*sin(e20)*cos(e3-e30)
        evr(2) = e10*(sin(e2)*cos(e20) - cos(e2)*sin(e20)*cos(e3-e30))
        evr(3) = e10*sin(e20)*sin(e3-e30)
      endif
!
!  Particles relative distance from each other
!
!  r_ij = sqrt(ev1**2 + ev2**2 + ev3**2)
!  invr3_ij = r_ij**(-3)
!
      if (present(distance)) then
        if (present(lsquare)) then
          distance =      sum(evr**2)
        else
          distance = sqrt(sum(evr**2))
        endif
      endif
!
      if (present(vector)) vector=evr
!
    endsubroutine get_particles_interdistance
!***********************************************************************
    subroutine remove_particle(fp,ipar,k,dfp,ineargrid,ks)
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc) :: ipar
      integer :: k
      real, dimension (mpar_loc,mpvar), optional :: dfp
      integer, dimension (mpar_loc,3), optional :: ineargrid
      integer, optional :: ks
!
      logical :: lnbody
      real :: t_sp   ! t in single precision for backwards compatibility
!
      intent (inout) :: fp, dfp, ineargrid
      intent (in)    :: k, ks
!
      t_sp = t
!
!  Check that we are not removing massive particles.
!
      lnbody=(lparticles_nbody.and.any(ipar(k)==ipar_nbody))
!
      if (lnbody) then
        if (present(ks)) then
          print*, 'nbody particle ', ks
          print*, 'is removing the following nbody particle:'
        else
          print*, 'the following nbody particle is being removed'
        endif
        print*, 'ipar(k)=', ipar(k)
        print*, 'xp,yp,zp=' ,fp(k,ixp), fp(k,iyp), fp(k,izp)
        print*, 'vxp,vyp,vzp=', fp(k,ivpx), fp(k,ivpy), fp(k,ivpz)
        call fatal_error("remove_particle","are you sure this should happen?")
      endif
!
!  Write to the respective processor that the particle is removed.
!  We also write the time.
!
!  TODO: It would be better to write this information in binary format to avoid
!  conversion problems when reading t_rmv with pc_read_pvar.
!
      open(20,file=trim(directory)//'/rmv_ipar.dat',position='append')
      if (present(ks)) then
        write(20,*) ipar(k), t_sp, ipar(ks)
      else
        write(20,*) ipar(k), t_sp
      endif
      close(20)
!
      open(20,file=trim(directory)//'/rmv_par.dat', &
          position='append',form='unformatted')
      if (present(ks)) then
        write(20) fp(k,:), fp(ks,:)
      else
        write(20) fp(k,:)
      endif
      close(20)
!
      if (ip<=8) print*, 'removed particle ', ipar(k)
!
!  Switch the removed particle with the last particle present in the processor
!  npar_loc
!
      fp(k,:)=fp(npar_loc,:)
      if (present(dfp)) dfp(k,:)=dfp(npar_loc,:)
      if (present(ineargrid)) ineargrid(k,:)=ineargrid(npar_loc,:)
      ipar(k)=ipar(npar_loc)
      if (lparticles_blocks) inearblock(k)=inearblock(npar_loc)
!
!  Reduce the number of particles by one.
!
      npar_loc=npar_loc-1
!
    endsubroutine remove_particle
!***********************************************************************
    subroutine count_particles(ipar,npar_found)
!
      use Mpicomm, only: mpireduce_sum_int
!
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_found
!
      npar_found=0
!
      call mpireduce_sum_int(npar_loc,npar_found)
      call keep_compiler_quiet(ipar)
!
    endsubroutine count_particles
!***********************************************************************
    subroutine output_particle_size_dist(fp)
!
!  Calculate size distribution of particles and output to file.
!
!   6-feb-11/anders: coded
!
      use Mpicomm, only: mpireduce_sum
!
      real, dimension (mpar_loc,mparray) :: fp
!
      real, dimension (nbin_ap_dist) :: log_ap_loc
      real, dimension (nbin_ap_dist) :: log_ap_loc_low, log_ap_loc_high
      real, dimension (nbin_ap_dist) :: rhop_dist, rhop_dist_sum
      real, save :: delta_log_ap
      integer :: i, k
      logical, save :: lfirstcall=.true.
      logical :: file_exists
!
      intent (in) :: fp
!
!  Define the bins for the size distribution and save to file.
!
      if (lfirstcall) then
        delta_log_ap=(log_ap_max_dist-log_ap_min_dist)/nbin_ap_dist
        do i=1,nbin_ap_dist
          log_ap_loc_low(i)=log_ap_min_dist+(i-1)*delta_log_ap
        enddo
        log_ap_loc_high = log_ap_loc_low + delta_log_ap
        log_ap_loc      = (log_ap_loc_low + log_ap_loc_high)/2
        if (lroot) then
          inquire(file=trim(datadir)//'/particle_size_dist.dat', &
              exist=file_exists)
          if (.not. file_exists) then
            open(20,file=trim(datadir)//'/particle_size_dist.dat', &
                status='new')
            write(20,*) log_ap_min_dist, log_ap_max_dist, nbin_ap_dist
            write(20,*) log_ap_loc
            write(20,*) log_ap_loc_low
            write(20,*) log_ap_loc_high
            close(20)
          endif
        endif
        lfirstcall=.false.
      endif
!
!  Loop over particles and add their mass to the radius bin.
!
      rhop_dist=0.0
      do k=1,npar_loc
        i=(alog10(fp(k,iap))-log_ap_min_dist)/delta_log_ap+1
        if (i>=1 .and. i<=nbin_ap_dist) then
          if (lparticles_number) then
            rhop_dist(i)=rhop_dist(i)+ &
                four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
          else
            rhop_dist(i)=rhop_dist(i)+rhop_const
          endif
        endif
      enddo
!
!  Sum over all processors
!
      call mpireduce_sum(rhop_dist,rhop_dist_sum,nbin_ap_dist)
!
!  Save size distribution function to file.
!
      if (lroot) then
        open(20,file=trim(datadir)//'/particle_size_dist.dat', &
            position='append')
        write(20,*) t, rhop_dist_sum
        close(20)
      endif
!
    endsubroutine output_particle_size_dist
!***********************************************************************
    subroutine get_rhopswarm_ineargrid(mp_swarm_tmp,fp,k,ineark,rhop_swarm_tmp)
!
!  Subroutine to calculate rhop_swarm, the mass density of a single
!  superparticle. The fundamental quantity is actually mp_swarm, the
!  mass of a superparticle. From that one calculates rhop_swarm=mp_swarm/dV,
!  where dV is the volume of a cell. In an equidistant Cartesian grid
!  this is simply a constant, and set in the beginning of the code. In
!  polar and non-equidistant grids dV varies with position such that
!
!    dV = dVol_x(l) * dVol_y(m) * dVol_z(n)
!
!  This subroutine takes care of this variation, also ensuring that the
!  impact on equidistant Cartesian grids is minimal.
!
!  Retrieves a scalar.
!
!  29-apr-11/wlad: coded
!
      real, intent(in)  :: mp_swarm_tmp
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k
      integer, dimension (3), intent(in) :: ineark
      real, intent(out) :: rhop_swarm_tmp
!
      integer  :: il,im,in
!
      if (lcartesian_coords.and.all(lequidist)) then
        if (lparticles_density) then
          rhop_swarm_tmp = fp(k,irhopswarm)
        elseif (lparticles_radius.and.lparticles_number) then
          rhop_swarm_tmp = &
              four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
        else
          rhop_swarm_tmp = rhop_swarm
        endif
      else
        il = ineark(1) ; im = ineark(2) ; in = ineark(3)
        rhop_swarm_tmp = mp_swarm_tmp*dVol1_x(il)*dVol1_y(im)*dVol1_z(in)
      endif
!
    endsubroutine get_rhopswarm_ineargrid
!***********************************************************************
    subroutine get_rhopswarm_point(mp_swarm_tmp,fp,k,il,im,in,rhop_swarm_tmp)
!
!  Same as get_rhopswarm_ineargrid, for general grid points il,im,in.
!  Retrieves a scalar.
!
!  29-apr-11/wlad: coded
!
      real, intent(in)  :: mp_swarm_tmp
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k, il, im, in
      real, intent(out) :: rhop_swarm_tmp
!
      if (lcartesian_coords.and.all(lequidist)) then
        if (lparticles_density) then
          rhop_swarm_tmp = fp(k,irhopswarm)
        elseif (lparticles_radius.and.lparticles_number) then
          rhop_swarm_tmp = &
              four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
        else
          rhop_swarm_tmp = rhop_swarm
        endif
      else
        rhop_swarm_tmp = mp_swarm_tmp*dVol1_x(il)*dVol1_y(im)*dVol1_z(in)
      endif
!
    endsubroutine get_rhopswarm_point
!***********************************************************************
    subroutine get_rhopswarm_block(mp_swarm_tmp,fp,k,il,im,in,ib,rhop_swarm_tmp)
!
!  Same as get_rhopswarm_point, but for the volume elements as block arrays.
!  Retrieves a scalar.
!
!  29-apr-11/wlad: coded
!
      use Particles_mpicomm, only: dVol1xb, dVol1yb, dVol1zb
!
      real, intent(in) :: mp_swarm_tmp
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k, il, im, in, ib
      real, intent(out) :: rhop_swarm_tmp
!
      if (lcartesian_coords.and.all(lequidist)) then
        if (lparticles_density) then
          rhop_swarm_tmp = fp(k,irhopswarm)
        elseif (lparticles_radius.and.lparticles_number) then
          rhop_swarm_tmp = &
              four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
        else
          rhop_swarm_tmp = rhop_swarm
        endif
      else
        rhop_swarm_tmp = mp_swarm_tmp*&
             dVol1xb(il,ib)*dVol1yb(im,ib)*dVol1zb(in,ib)
      endif
!
    endsubroutine get_rhopswarm_block
!***********************************************************************
    subroutine get_rhopswarm_pencil(mp_swarm_tmp,fp,k,im,in,rhop_swarm_tmp)
!
!  Same as get_rhopswarm_ineargrid, but retrieves a pencil.
!
!  29-apr-11/wlad: coded
!
      real, intent(in) :: mp_swarm_tmp
      real, dimension (mpar_loc,mparray), intent(in) :: fp
      integer, intent(in) :: k, im, in
      real, dimension(nx), intent(out) :: rhop_swarm_tmp
!
      if (lcartesian_coords.and.all(lequidist)) then
        if (lparticles_density) then
          rhop_swarm_tmp = fp(k,irhopswarm)
        elseif (lparticles_radius.and.lparticles_number) then
          rhop_swarm_tmp = &
              four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
        else
          rhop_swarm_tmp = rhop_swarm
        endif
      else
        rhop_swarm_tmp = mp_swarm_tmp*dVol1_x(l1:l2)*dVol1_y(im)*dVol1_z(in)
      endif
!
    endsubroutine get_rhopswarm_pencil
!***********************************************************************
    subroutine find_grid_volume(ixx,iyy,izz,volume_cell)
!
!  Find the volume of the grid cell
!
!  24-sep-14/nils: coded
!
      real, intent(out) :: volume_cell
      integer, intent(in) :: ixx,iyy,izz
      real :: dx1, dy1, dz1
!
      if (nxgrid ==1) then
        dx1=1./Lxyz(1)
      else
        dx1=dx_1(ixx)
      endif
      if (nygrid ==1) then
        dy1=1./Lxyz(2)
      else
        dy1=dy_1(iyy)
      endif
      if (nzgrid ==1) then
        dz1=1./Lxyz(3)
      else
        dz1=dz_1(izz)
      endif
!
      volume_cell=1./(dx1*dy1*dz1)
!
      if (lcylindrical_coords .and. nygrid/=1) then
        volume_cell = volume_cell*x(ixx)
      elseif (lspherical_coords) then
        if (nygrid/=1) volume_cell = volume_cell*x(ixx)
        if (nzgrid/=1) volume_cell = volume_cell*sin(y(iyy))*x(ixx)
      endif
!
    end subroutine find_grid_volume
!***********************************************************************
    subroutine find_interpolation_weight(weight,fp,k,ixx,iyy,izz,ix0,iy0,iz0)
!
      real :: weight
      real :: weight_x, weight_y, weight_z
      integer, intent(in) :: k,ixx,iyy,izz,ix0,iy0,iz0
      real, dimension (mpar_loc,mparray), intent(in) :: fp
! 
      if (lparticlemesh_cic) then
!
!  Cloud In Cell (CIC) scheme.
!
!  Particle influences the 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
        weight=1.0
        if (nxgrid/=1) &
            weight=weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
        if (nygrid/=1) &
            weight=weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
        if (nzgrid/=1) &
            weight=weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
      elseif (lparticlemesh_tsc) then
!
!  Triangular Shaped Cloud (TSC) scheme.
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
        if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
          weight_x=1.125-1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
              0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
        else
          if (nxgrid/=1) &
              weight_x=0.75-((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
        endif
        if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
          weight_y=1.125-1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
              0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
        else
          if (nygrid/=1) &
              weight_y=0.75-((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
        endif
        if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
          weight_z=1.125-1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
              0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
        else
          if (nzgrid/=1) &
              weight_z=0.75-((fp(k,izp)-z(izz))*dz_1(izz))**2
        endif
!
        weight=1.0
!
        if (nxgrid/=1) weight=weight*weight_x
        if (nygrid/=1) weight=weight*weight_y
        if (nzgrid/=1) weight=weight*weight_z
      elseif (lparticlemesh_gab) then
!
!  Gaussian box approach. The particles influence is distributed over 
!  the nearest grid point and 3 nodes in each dimension. The weighting of the 
!  source term follows a gaussian distribution
!
        weight = 1.
        if (nxgrid/=1) weight=weight*gab_weights(abs(ixx-ix0)+1)
        if (nygrid/=1) weight=weight*gab_weights(abs(iyy-iy0)+1)
        if (nzgrid/=1) weight=weight*gab_weights(abs(izz-iz0)+1)
      else
!
!  Nearest Grid Point (NGP) scheme.
!
        weight=1.
      endif
!
    end subroutine find_interpolation_weight
!***********************************************************************
    subroutine find_interpolation_indeces(ixx0,ixx1,iyy0,iyy1,izz0,izz1,&
        fp,k,ix0,iy0,iz0)
!
      integer, intent(in)  :: k,ix0,iy0,iz0
      integer, intent(out) :: ixx0,ixx1,iyy0,iyy1,izz0,izz1
      real, dimension (mpar_loc,mparray), intent(in) :: fp
!
!  Cloud In Cell (CIC) scheme.
!
      if (lparticlemesh_cic) then
        ixx0=ix0; iyy0=iy0; izz0=iz0
        ixx1=ix0; iyy1=iy0; izz1=iz0
!
!  Particle influences the 8 surrounding grid points. The reference point is
!  the grid point at the lower left corner.
!
        if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
        if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
        if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
        if (nxgrid/=1) ixx1=ixx0+1
        if (nygrid/=1) iyy1=iyy0+1
        if (nzgrid/=1) izz1=izz0+1
      elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
        if (nxgrid/=1) then
          ixx0=ix0-1; ixx1=ix0+1
        else
          ixx0=ix0  ; ixx1=ix0
        endif
        if (nygrid/=1) then
          iyy0=iy0-1; iyy1=iy0+1
        else
          iyy0=iy0  ; iyy1=iy0
        endif
        if (nzgrid/=1) then
          izz0=iz0-1; izz1=iz0+1
        else
          izz0=iz0  ; izz1=iz0
        endif
      elseif (lparticlemesh_gab) then
!
!  Gaussian box approach. The particles influence is distributed over 
!  the nearest grid point and 3 nodes in each dimension, for a total
!  of 81 influenced points in 3 dimensions
!
        if (nxgrid/=1) then 
          ixx0=ix0-3 
          ixx1=ix0+3
        else
          ixx0=ix0
          ixx1=ix0
        endif
        if (nygrid/=1) then 
          iyy0=iy0-3 
          iyy1=iy0+3
        else
          iyy0=iy0
          iyy1=iy0
        endif
        if (nzgrid/=1) then 
          izz0=iz0-3 
          izz1=iz0+3
        else
          izz0=iz0
          izz1=iz0
        endif
      else
!
!  Nearest Grid Point (NGP) scheme.
!
        ixx0=ix0
        ixx1=ixx0
        iyy0=m
        iyy1=iyy0
        izz0=n
        izz1=izz0
      endif
!
    end subroutine find_interpolation_indeces
!***********************************************************************
    real function get_gas_density(f, ix, iy, iz) result(rho)
!
!  Reads the gas density at location (ix, iy, iz).
!
!  19-jun-14/ccyang: coded.
!
      use EquationOfState, only: get_stratz, rho0
      use DensityMethods , only: getrho_s
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: ix, iy, iz
!
      real, dimension(mz) :: rho0z = 0.0
      logical :: lfirstcall = .true.
!
!  Initialization
!
      init: if (lfirstcall) then
        if (lstratz) call get_stratz(z, rho0z)
        lfirstcall = .false.
      endif init
!
!  Find the gas density at (ix,iy,iz).
!
      if (lstratz) then
        rho = rho0z(iz) * (1.0 + f(ix,iy,iz,irho))
      elseif (ilnrho /= 0) then
        rho = getrho_s(f(ix,iy,iz,ilnrho), ix)
      else
        rho = rho0
      endif
!
    endfunction get_gas_density
!***********************************************************************
    subroutine precalc_cell_volumes(f)
!
      real, dimension (mx,my,mz,mfarray), intent(out) :: f
      integer :: ixx, iyy, izz
      real :: dx1, dy1, dz1, cell_volume1
!
      do ixx=l1,l2
        do iyy=m1,m2
          do izz=n1,n2
            dx1=dx_1(ixx)
            dy1=dy_1(iyy)
            dz1=dz_1(izz)
            cell_volume1 = 1.0
            if (lcartesian_coords) then
              if (nxgrid/=1) cell_volume1 = cell_volume1*dx1
              if (nygrid/=1) cell_volume1 = cell_volume1*dy1
              if (nzgrid/=1) cell_volume1 = cell_volume1*dz1
            elseif (lcylindrical_coords) then
              if (nxgrid/=1) cell_volume1 = cell_volume1*dx1
              if (nygrid/=1) cell_volume1 = cell_volume1*dy1/x(ixx)
              if (nzgrid/=1) cell_volume1 = cell_volume1*dz1
            elseif (lspherical_coords) then
              if (nxgrid/=1) cell_volume1 = cell_volume1*dx1
              if (nygrid/=1) cell_volume1 = cell_volume1*dy1/x(ixx)
              if (nzgrid/=1) cell_volume1 = cell_volume1*dz1/(sin(y(iyy))*x(ixx))
            endif
            f(ixx,iyy,izz,ivol) = 1.0/cell_volume1
          enddo
        enddo
      enddo
    endsubroutine precalc_cell_volumes
!***********************************************************************
    subroutine precalc_weights(weight_array)
!
      real, dimension(:,:,:) :: weight_array
      real, dimension(3) :: tsc_values = (/0.125, 0.75, 0.125/)
      integer :: i,j,k
!
!  This makes sure that the weight array is 1 if the npg approach is chosen
!
      weight_array(:,:,:) = 1.0
!
      if (lparticlemesh_gab) then
        do i=1,7
          do j=1,7
            do k=1,7
              weight_array(i,j,k) = gab_weights(abs(i-4)+1)* &
                  gab_weights(abs(j-4)+1) * gab_weights(abs(k-4)+1)
            enddo
          enddo
        enddo
      endif

      if (lparticlemesh_tsc) then
          do i = 1,3
            do j = 1,3
              do k =1,3
                if (nxgrid/=1) weight_array(i,j,k)=weight_array(i,j,k)*tsc_values(i)
                if (nygrid/=1) weight_array(i,j,k)=weight_array(i,j,k)*tsc_values(j)
                if (nzgrid/=1) weight_array(i,j,k)=weight_array(i,j,k)*tsc_values(k)
            enddo
          enddo
        enddo
      endif

      if (lparticlemesh_cic) then
        weight_array=1.0
        if (nxgrid/=1) &
            weight_array=weight_array*0.5
        if (nygrid/=1) &
            weight_array=weight_array*0.5
        if (nzgrid/=1) &
            weight_array=weight_array*0.5
      endif
      
!
    endsubroutine precalc_weights
!***********************************************************************
endmodule Particles_sub
