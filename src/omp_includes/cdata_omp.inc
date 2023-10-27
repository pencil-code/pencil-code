!$omp threadprivate(dxyz_2)
!$omp threadprivate(dxyz_4)
!$omp threadprivate(dxyz_6)
!$omp threadprivate(dvol)
!$omp threadprivate(dxmax_pencil)
!$omp threadprivate(dxmin_pencil)
!$omp threadprivate(dline_1)
!$omp threadprivate(lcoarse_mn)
!$omp threadprivate(advec2_hypermesh)
!$omp threadprivate(lfirstpoint)
!$omp threadprivate(seed)
!$omp threadprivate(m)
!$omp threadprivate(n)
!$omp threadprivate(itype_name)
!$omp threadprivate(fname)
!$omp threadprivate(fnamer)
!$omp threadprivate(fname_sound)
!$omp threadprivate(fnamex)
!$omp threadprivate(fnamey)
!$omp threadprivate(fnamez)
!$omp threadprivate(fnamexy)
!$omp threadprivate(fnamexz)
!$omp threadprivate(fnamerz)
!$omp threadprivate(ncountsz)
!$omp threadprivate(fname_half)

contains
subroutine copyin_cdata
!$omp parallel copyin( &
!$omp dxyz_2  &
!$omp ,dxyz_4  &
!$omp ,dxyz_6  &
!$omp ,dvol  &
!$omp ,dxmax_pencil  &
!$omp ,dxmin_pencil  &
!$omp ,dline_1  &
!$omp ,lcoarse_mn  &
!$omp ,advec2_hypermesh  &
!$omp ,lfirstpoint  &
!$omp ,seed  &
!$omp ,m  &
!$omp ,n  &
!$omp ,itype_name  &
!$omp ,fname  &
!$omp ,fnamer  &
!$omp ,fname_sound  &
!$omp ,fnamex  &
!$omp ,fnamey  &
!$omp ,fnamez  &
!$omp ,fnamexy  &
!$omp ,fnamexz  &
!$omp ,fnamerz  &
!$omp ,ncountsz  &
!$omp ,fname_half  &
!$omp )
!$omp end parallel
endsubroutine copyin_cdata

