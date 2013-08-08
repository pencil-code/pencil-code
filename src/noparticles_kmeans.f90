! $Id: particles_coagulation.f90 19828 2012-11-27 09:58:06Z kalle.jansson.89 $
!
!  Dummy version of k-means vector quantization algorithm using the Forgy
!  initialisation of the codebook.
!
module Particles_kmeans
!
  use Messages
!
  contains
!***********************************************************************
    subroutine ppcvq(measure, dims, dim1, dim2, no_tr_vec, trainv, weights, & 
        cb_size, codebk, codew, do_split, do_norm, do_detail)
!
!  08-aug-13/anders+chao-chin: dummy
!
      integer :: measure, dims, dim1, dim2, no_tr_vec
      real, dimension(dims,no_tr_vec) :: trainv
      real, dimension(no_tr_vec) :: weights
      integer :: cb_size
      real, dimension(dims,cb_size) :: codebk
      real, dimension(cb_size) :: codew
      logical :: do_split, do_norm, do_detail
!
      call fatal_error('ppcvq','you need to use the k-means module')
!
      call keep_compiler_quiet(measure,dims,dim1,dim2,no_tr_vec)
      call keep_compiler_quiet(trainv)
      call keep_compiler_quiet(weights)
      call keep_compiler_quiet(cb_size)
      call keep_compiler_quiet(codebk)
      call keep_compiler_quiet(codew)
      call keep_compiler_quiet(do_split,do_norm,do_detail)
!
    endsubroutine ppcvq
!***********************************************************************
endmodule Particles_kmeans
