!used to calculate the Morton decomp given the number of processes in cparam.local
!used by simply including this file in cparam.local

integer, parameter :: nprocs_pow_two = int(log(real(ncpus))/log(real(2)))
integer, parameter :: nprocx= 2**(nprocs_pow_two/3)
integer, parameter :: nprocy = 2**((nprocs_pow_two+1)/3)
integer, parameter :: nprocz = 2**((nprocs_pow_two+2)/3)
