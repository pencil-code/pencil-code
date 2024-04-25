!***********************************************************************
   subroutine read_transport_data
!
!  Reading of the chemkin transport data
!
!  01-apr-08/natalia: coded
!  30-jun-17/MR: moved here from eos_chemistry.
!
      use General, only: itoa

      logical :: emptyfile
      logical :: found_specie
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      character (len=10) :: specie_string
      integer :: VarNumber
      integer :: StartInd,StopInd,StartInd_1,StopInd_1
      logical :: tranin=.false.
      logical :: trandat=.false.
!
      emptyFile=.true.
!
      StartInd_1=1; StopInd_1=0

      inquire (file='tran.dat',exist=trandat)
      inquire (file='tran.in',exist=tranin)
      if (tranin .and. trandat) &
        call fatal_error('read_transport_data', &
            'both tran.in and tran.dat found. Please decide which to use.')

      if (tranin) open(file_id,file='tran.in')
      if (trandat) open(file_id,file='tran.dat')
!
      if (lroot) print*, 'the following species are found in tran.in/dat:'
!
      dataloop: do
!
        read(file_id,'(80A)',end=1000) ChemInpLine(1:80)
        emptyFile=.false.
!
        StopInd_1=index(ChemInpLine,' ')
        specie_string=trim(ChemInpLine(1:StopInd_1-1))
!
        call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
!
        if (found_specie) then
          if (lroot) print*,specie_string,' ind_glob=',ind_glob,' ind_chem=',ind_chem
!
          VarNumber=1; StartInd=1; StopInd =0
          do while (VarNumber<7)
!
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
            StartInd=verify(ChemInpLine(StopInd:),' ')+StopInd-1
            StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
!
            if (StopInd==StartInd) then
              StartInd=StartInd+1
            else
              if (VarNumber==1) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E1.0)')  tran_data(ind_chem,VarNumber)
              elseif (VarNumber==2) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') tran_data(ind_chem,VarNumber)
              elseif (VarNumber==3) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') tran_data(ind_chem,VarNumber)
              elseif (VarNumber==4) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') tran_data(ind_chem,VarNumber)
              elseif (VarNumber==5) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') tran_data(ind_chem,VarNumber)
              elseif (VarNumber==6) then
                read (unit=ChemInpLine(StartInd:StopInd),fmt='(E15.8)') tran_data(ind_chem,VarNumber)
              else
                call fatal_error("read_transport_data","no such VarNumber: "//trim(itoa(VarNumber)))
              endif
!
              VarNumber=VarNumber+1
              StartInd=StopInd
            endif
            if (StartInd==80) exit
          enddo
!
        endif
      enddo dataloop
!
! Stop if tran.dat is empty
!
1000  if (emptyFile) call fatal_error('read_transport_data','input file tran.dat is empty')
!
      close(file_id)
!
    endsubroutine read_transport_data
!***********************************************************************
    subroutine find_mass(element_name,MolMass)
!
!  Find mass of element
!
!  05-feb-08/nils: coded
!
      use Mpicomm, only: stop_it
!
      character (len=*), intent(in) :: element_name
      real, intent(out) :: MolMass
!
      select case (element_name)
      case ('H')
        MolMass=1.00794
      case ('C')
        MolMass=12.0107
      case ('N')
        MolMass=14.00674
      case ('O')
        MolMass=15.9994
      case ('Ar','AR')
        MolMass=39.948
      case ('He','HE')
        MolMass=4.0026
      case ('S')
        MolMass=32.0655
      case ('SI')	
        MolMass=28.085
      case ('CLOUD')
        MolMass=0.
      case default
        if (lroot) print*,'element_name=',element_name
        call fatal_error('find_mass','no such element_name: '//trim(element_name))
      end select
!
    endsubroutine find_mass
!***********************************************************************
   subroutine find_species_index(species_name,ind_glob,ind_chem,found_specie)
!
!  Find index in the f array for specie
!
!  05-feb-08/nils: coded
!
      integer, intent(out) :: ind_glob
      integer, intent(inout) :: ind_chem
      character (len=*), intent(in) :: species_name
      integer :: k
      logical, intent(out) :: found_specie
!
      ind_glob=0
    !  ind_chem=0
      do k=1,nchemspec
        if (trim(varname(ichemspec(k)))==species_name) then
          ind_glob=k+ichemspec(1)-1
          ind_chem=k
          exit
        endif
!print*, trim(varname(ichemspec(k))),(species_name)
      enddo
!
!  Check if the species was really found
!
      if ((ind_glob==0)) then
        found_specie=.false.
     !  if (lroot) print*,' no species has been found  ',' species index= ', ind_glob,ind_chem,species_name
        !call fatal_error('find_species_index','index for '//trim(species_name)//'not found')
      else
        found_specie=.true.
    !    if (lroot) print*,species_name,'   species index= ',ind_chem
      endif
!
    endsubroutine find_species_index
!***********************************************************************
    subroutine read_species(input_file)
!
!  This subroutine reads all species information from chem.inp
!  See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  06-mar-08/nils: coded
!
      logical :: IsSpecie=.false., emptyfile
      integer :: k,file_id=123, StartInd, StopInd
      character (len=80) :: ChemInpLine
      character (len=*) :: input_file
!
      emptyFile=.true.
      k=1
      open(file_id,file=input_file)
      do
        read(file_id,'(80A)',end=1000) ChemInpLine(1:80)
        emptyFile=.false.
!
!  Check if we are reading a line within the species section
!
        if (ChemInpLine(1:7)=="SPECIES")            IsSpecie=.true.
        if (ChemInpLine(1:3)=="END" .and. IsSpecie) IsSpecie=.false.
!
!  Read in species
!
        if (IsSpecie) then
          if (ChemInpLine(1:7) /= "SPECIES") then
            StartInd=1; StopInd =0
            do
              StopInd=index(ChemInpLine(StartInd:),' ')+StartInd-1
              if (StopInd==StartInd) then
                StartInd=StartInd+1
              else
                if (k>nchemspec) then
                  print*,'nchemspec=',nchemspec
                  call fatal_error('read_species',"too many species, increase nchemspec")
                endif
                varname(ichemspec(k))=trim(ChemInpLine(StartInd:StopInd-1))
                StartInd=StopInd
                k=k+1
              endif
              if (StartInd==80) exit
            enddo
          endif
        endif
      enddo
!
!  Stop if chem.inp is empty
!
1000  if (emptyFile)  call fatal_error('read_species','input file chem.inp is empty')
!
!  Check if nchemspec where not too large  MR: why fatal?
!
      if (k<nchemspec-1) then
        print*,'nchemspec=',nchemspec
        call fatal_error("read_species","too few species, decrease nchemspec")
      endif
!
      close(file_id)
!
    endsubroutine read_species
!***********************************************************************
   subroutine read_thermodyn(input_file)
!
!  This subroutine reads the thermodynamical data for all species
!  from chem.inp. See the chemkin manual for more information on
!  the syntax of chem.inp.
!
!  06-mar-08/nils: coded
!
      character (len=*), intent(in) :: input_file
      integer :: file_id=123, ind_glob, ind_chem
      character (len=80) :: ChemInpLine
      integer :: In1,In2,In3,In4,In5,iElement,iTemperature,StopInd
      integer :: NumberOfElement_i
      logical :: IsThermo=.false., found_specie, existing_specie
      real, dimension(4) :: MolMass
      real, dimension(3) :: tmp_temp
      character (len=5) :: NumberOfElement_string,element_string
      character (len=10) :: specie_string,TemperatureNr_i
      real :: nne
      integer, dimension(7) :: iaa1,iaa2
!
      integer :: iTemp1=2,iTemp2=3,iTemp3=4
!
      ind_chem=0
!
!  Initialize some index pointers
!
      iaa1(1)=5;iaa1(2)=6;iaa1(3)=7;iaa1(4)=8
      iaa1(5)=9;iaa1(6)=10;iaa1(7)=11
!
      iaa2(1)=12;iaa2(2)=13;iaa2(3)=14;iaa2(4)=15
      iaa2(5)=16;iaa2(6)=17;iaa2(7)=18
!
      open(file_id,file=input_file)
      dataloop2: do
        read(file_id,'(80A)',end=1001) ChemInpLine(1:80)
!
! Check if we are reading a line within the thermo section
!
        if (ChemInpLine(1:6)=="THERMO") IsThermo=.true.
        if (ChemInpLine(1:3)=="END" .and. IsThermo) IsThermo=.false.
!
! Read in thermo data
!
        if (IsThermo) then
          if (ChemInpLine(1:7) /= "THERMO") then
            StopInd=index(ChemInpLine,' ')
            specie_string=trim(ChemInpLine(1:StopInd-1))
!
            call find_species_index(specie_string,ind_glob,ind_chem,found_specie)
!
! Check if we are working with a specie that was found under the SPECIES
! section of chem.inp.
!
            if (ChemInpLine(80:80)=="1") then
              if (found_specie) then
                existing_specie=.true.
              else
                existing_specie=.false.
              endif
            endif
!
! What problems are in the case of  ind_chem=0?
!
            if (ind_chem>0 .and. ind_chem<=nchemspec) then
!
              if (existing_specie.and.found_specie) then
!
! Find molar mass
!
                MolMass=0
                do iElement=1,4
                  In1=25+(iElement-1)*5
                  In2=26+(iElement-1)*5
                  In3=27+(iElement-1)*5
                  In4=29+(iElement-1)*5
                  if (ChemInpLine(In1:In1)==' ') then
                    MolMass(iElement)=0
                  else
                    element_string=trim(ChemInpLine(In1:In2))
                    call find_mass(element_string,MolMass(iElement))
                    In5=verify(ChemInpLine(In3:In4),' ')+In3-1
                    NumberOfElement_string=trim(ChemInpLine(In5:In4))
                    read (unit=NumberOfElement_string,fmt='(I5)') NumberOfElement_i
                    MolMass(iElement)=MolMass(iElement)*NumberOfElement_i
                  endif
                enddo
                species_constants(ind_chem,imass)=sum(MolMass)
!
! Find temperature-ranges for low and high temperature fitting
!
                do iTemperature=1,3
                  In1=46+(iTemperature-1)*10
                  In2=55+(iTemperature-1)*10
                  if (iTemperature==3) In2=73
                  In3=verify(ChemInpLine(In1:In2),' ')+In1-1
                  TemperatureNr_i=trim(ChemInpLine(In3:In2))
                  read (unit=TemperatureNr_i,fmt='(F10.1)') nne
                  tmp_temp(iTemperature)=nne
                enddo
                species_constants(ind_chem,iTemp1)=tmp_temp(1)
                species_constants(ind_chem,iTemp2)=tmp_temp(3)  !MR:indices correct?
                species_constants(ind_chem,iTemp3)=tmp_temp(2)
!
              elseif (ChemInpLine(80:80)=="2") then
                ! Read iaa1(1):iaa1(5)
                read (unit=ChemInpLine(1:75),fmt='(5E15.8)') species_constants(ind_chem,iaa1(1):iaa1(5))
!
              elseif (ChemInpLine(80:80)=="3") then
                ! Read iaa1(6):iaa5(3)
                read (unit=ChemInpLine(1:75),fmt='(5E15.8)') species_constants(ind_chem,iaa1(6):iaa2(3))
              elseif (ChemInpLine(80:80)=="4") then
                ! Read iaa2(4):iaa2(7)
                read (unit=ChemInpLine(1:75),fmt='(4E15.8)') species_constants(ind_chem,iaa2(4):iaa2(7))
              endif
!
            endif  ! if (existing_specie.and.found_specie)
          endif    ! if (ind_chem>0 .and. ind_chem<=nchemspec)
        endif
      enddo dataloop2
1001  continue
      close(file_id)
!
    endsubroutine read_thermodyn
!***********************************************************************
    subroutine write_thermodyn
!
!  This subroutine writes the thermodynamical data for every specie
!  to ./data/chem.out.
!
!  06-mar-08/nils: coded
!
      use General
!
      character (len=fnlen) :: input_file="./data/chem.out"
      character (len=intlen) :: ispec
      integer :: file_id=123,k
      integer, dimension(7) :: iaa1,iaa2
      integer :: iTemp1=2,iTemp3=4
!
!      Initialize some index pointers
!
      iaa1(1)=5;iaa1(2)=6;iaa1(3)=7;iaa1(4)=8
      iaa1(5)=9;iaa1(6)=10;iaa1(7)=11
!
      iaa2(1)=12;iaa2(2)=13;iaa2(3)=14;iaa2(4)=15
      iaa2(5)=16;iaa2(6)=17;iaa2(7)=18
!
      open(file_id,file=input_file)
      write(file_id,*) 'Specie'
      write(file_id,*) 'MolMass Temp1 Temp2 Temp3'
      write(file_id,*) 'a1(1)  a1(2)  a1(3)  a1(4)  a1(5)  a1(6)  a1(7)'
      write(file_id,*) 'a2(1)  a2(2)  a2(3)  a2(4)  a2(5)  a2(6)  a2(7)'
      write(file_id,*) '***********************************************'
      do k=1,nchemspec
        write(file_id,*) varname(ichemspec(k))
        write(file_id,'(F10.2,3F10.2)') species_constants(k,imass),species_constants(k,iTemp1:iTemp3)
        write(file_id,'(7E12.5)') species_constants(k,iaa1)
        write(file_id,'(7E12.5)') species_constants(k,iaa2)
      enddo
!
      close(file_id)
!
      if (lroot) then
        print*,'Write pc_constants.pro in chemistry.f90'
        open (143,FILE=trim(datadir)//'/pc_constants.pro',POSITION="append")
        write (143,*) 'specname=strarr(',nchemspec,')'
        write (143,*) 'specmass=fltarr(',nchemspec,')'
        do k=1,nchemspec
          ispec=itoa(k-1)
          write (143,*) 'specname[',trim(ispec),']=',"'",trim(varname(ichemspec(k))),"'"
          write (143,*) 'specmass[',trim(ispec),']=',species_constants(k,imass)
        enddo
        close (143)
      endif
!
    endsubroutine write_thermodyn
!***********************************************************************
