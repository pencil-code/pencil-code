! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ywater, lambda; chem_conc(nchemspec)
! PENCILS PROVIDED nucl_rmin, nucl_rate, conc_satm, ff_cond
! PENCILS PROVIDED latent_heat
!
!***************************************************************
module Chemistry
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  real :: Rgas, Rgas_unit_sys
  logical :: lchemistry_diag=.false.
  logical :: lreactions=.false.
  logical, allocatable, dimension(:,:,:) :: lnucleii_generated
!
  include 'chemistry.h'
!
  real, dimension(0,0) :: species_constants

  contains
!***********************************************************************
    subroutine register_chemistry
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine initialize_chemistry(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      Rgas_unit_sys = k_B_cgs/m_u_cgs
      Rgas = Rgas_unit_sys/unit_energy
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chemistry
!*********************************************************************** 
    subroutine init_chemistry(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_chemistry
!***********************************************************************
    subroutine pencil_criteria_chemistry
!
    endsubroutine pencil_criteria_chemistry
!***********************************************************************
    subroutine pencil_interdep_chemistry(lpencil_in)
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine chemistry_before_boundary(f)

      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)

    endsubroutine chemistry_before_boundary
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_chemistry
!***********************************************************************
    subroutine calc_for_chem_mixture(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_for_chem_mixture
!***********************************************************************
    subroutine dchemistry_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dchemistry_dt
!***********************************************************************
    subroutine calc_diagnostics_chemistry(f,p)
!
!  Calculate diagnostic quantities
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_chemistry
!***********************************************************************
    subroutine rprint_chemistry(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_chemistry
!***********************************************************************
    subroutine chemspec_normalization(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine chemspec_normalization
!***********************************************************************
    subroutine get_slices_chemistry(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_chemistry
!***********************************************************************
    subroutine chemistry_clean_up
!
    endsubroutine chemistry_clean_up
!***********************************************************************
    subroutine read_chemistry_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      iostat=0

    endsubroutine read_chemistry_init_pars
!***********************************************************************
    subroutine write_chemistry_init_pars(unit)
!
      integer, intent(in) :: unit
!
    endsubroutine write_chemistry_init_pars
!***********************************************************************
    subroutine read_chemistry_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      iostat=0

    endsubroutine read_chemistry_run_pars
!***********************************************************************
    subroutine write_chemistry_run_pars(unit)
!
      integer, intent(in) :: unit
!
    endsubroutine write_chemistry_run_pars
!***********************************************************************
    subroutine jacobn(f,jacob)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,nchemspec,nchemspec) :: jacob
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(jacob(1,1,1,1,1))
!
    endsubroutine jacobn
!***********************************************************************
    subroutine get_mu1_slice(f,slice,grad_slice,index,sgn,direction)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(ny,nz), intent(out) :: slice, grad_slice
      integer, intent(in) :: index, sgn,direction
!
      call keep_compiler_quiet(slice)
      call keep_compiler_quiet(grad_slice)
      call keep_compiler_quiet(index,sgn,direction)
!
    end subroutine get_mu1_slice
!***********************************************************************
    subroutine get_gamma_slice(f,slice,index,dir)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (:,:), intent(out)  :: slice
      integer, intent(in) :: index,dir
!
      call keep_compiler_quiet(slice)
      call keep_compiler_quiet(index)
      call keep_compiler_quiet(dir)
      !
    endsubroutine get_gamma_slice
!***********************************************************************
    subroutine get_cs2_slice(f,slice,index,dir)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (:,:), intent(out)  :: slice
      integer, intent(in) :: index,dir
!
      call keep_compiler_quiet(slice)
      call keep_compiler_quiet(index)
      call keep_compiler_quiet(dir)
      !
    endsubroutine get_cs2_slice
!***********************************************************************
   subroutine get_cs2_full(cs2_full)
!
      real, dimension (mx,my,mz) :: cs2_full
!
      intent(out) :: cs2_full
!
      call keep_compiler_quiet(cs2_full)
!
    endsubroutine get_cs2_full
!***********************************************************************
    subroutine get_gamma_full(gamma_full)
!
      real, dimension (mx,my,mz) :: gamma_full
!
      intent(out) :: gamma_full
!
      call keep_compiler_quiet(gamma_full)
!
    endsubroutine get_gamma_full
!***********************************************************************
    subroutine get_RHS_Y_full(RHS_Y)
!
      real, dimension (mx,my,mz,nchemspec) :: RHS_Y
!
      intent(out) :: RHS_Y
!
      call keep_compiler_quiet(RHS_Y)
!
    endsubroutine get_RHS_Y_full
!***********************************************************************
    subroutine  write_net_reaction
    endsubroutine  write_net_reaction
!***********************************************************************
    subroutine get_reac_rate(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine get_reac_rate
!!***********************************************************************
    subroutine chemspec_normalization_N2(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine chemspec_normalization_N2
!***********************************************************************
    subroutine chemistry_init_reduc_pointers
!
    endsubroutine chemistry_init_reduc_pointers
!***********************************************************************
    subroutine chemistry_diags_reductions
!
    endsubroutine chemistry_diags_reductions 
!***********************************************************************
   subroutine find_species_index(species_name,ind_glob,ind_chem,found_specie)
!
      integer, intent(out) :: ind_glob
      integer, intent(inout) :: ind_chem
      character (len=*), intent(in) :: species_name
      logical, intent(out) :: found_specie
!
      call keep_compiler_quiet(ind_glob)
      call keep_compiler_quiet(found_specie)

   endsubroutine find_species_index
!***********************************************************************
    subroutine cond_spec_cond(f,df,p,ad,dustbin_width,mfluxcond)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: dustbin_width
      real, dimension (nx) :: mfluxcond
      real, dimension(ndustspec) :: ad
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ad)
      call keep_compiler_quiet(dustbin_width)
      call keep_compiler_quiet(mfluxcond)
!  
    end subroutine cond_spec_cond
!***********************************************************************
    subroutine cond_spec_nucl(f,df,p,kk_vec,ad)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      integer, dimension(nx) :: kk_vec
      real, dimension(ndustspec) :: ad
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ad)
      call keep_compiler_quiet(kk_vec)
!  
    end subroutine cond_spec_nucl
!***********************************************************************
    subroutine condensing_species_rate(p,mfluxcond)
!
      real, dimension (nx) :: mfluxcond
      type (pencil_case) :: p
!
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(mfluxcond)
!      
    end subroutine condensing_species_rate
!***********************************************************************
    subroutine cond_spec_cond_lagr(f,df,p,rp,ix0,ix,np_swarm,dapdt)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: rp,np_swarm,dapdt
      integer :: ix0,ix
 !
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(dapdt)
      call keep_compiler_quiet(ix0)
!
    end subroutine cond_spec_cond_lagr
!***********************************************************************
    subroutine cond_spec_nucl_lagr(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    end subroutine cond_spec_nucl_lagr
!***********************************************************************
    subroutine pushpars2c(p_par)

      use Syscalls, only: copy_addr

      integer, parameter :: n_pars=1
      integer(KIND=ikind8), dimension(n_pars) :: p_par

      call copy_addr(rgas,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
    subroutine make_flame_index(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine make_flame_index
!***********************************************************************
    subroutine make_mixture_fraction(f)
!
! Calculate Bilger mixture fraction and store in f-array. 
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine make_mixture_fraction
!***********************************************************************   
endmodule Chemistry
