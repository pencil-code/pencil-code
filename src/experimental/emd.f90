
module EMD
  implicit none
  integer,parameter   :: max_rounds = 50
  real(kind=8),parameter   :: err_coeff = 5.0e-2

  contains

  function analyzer_emd(dataset,xdim1,xdim2,tlen,resultlen,resultlen2) result(analysis)
    implicit none
    integer, intent(in) :: xdim1,xdim2,tlen,resultlen,resultlen2
    integer       :: iIMF,i1,i2,it,iround,crossings,extremas
    logical       :: lIMF, lIMFstop
    real(kind=8), dimension(xdim1,xdim2,tlen), intent(in)       :: dataset
    real(kind=8), dimension(xdim1,xdim2,resultlen,resultlen2)   :: analysis
    real(kind=8), dimension(tlen)          :: data_area,imf
    real(kind=8), dimension(tlen)       :: maxima,minima,yminima,ymaxima,&
                                                   tminima,tmaxima,sminima,smaxima
    real(kind=8)  :: max_mean,env_mean,imf_avg,imf_std
    integer                   :: nminima, nmaxima
    
    ! calculate Hilbert-Huang through IMFs

    analysis(1:xdim1,1:xdim2,1:resultlen,1:resultlen2) = 0
    analysis(1:xdim1,1:xdim2,1:tlen,1) = dataset(1:xdim1,1:xdim2,1:tlen)
    jloop: do i2=1,xdim2
      iloop: do i1=1,xdim1
        data_area(1:tlen) = 0
        data_area(1:tlen) = analysis(i1,i2,1:tlen,1)
        !write(*,*) '1 min:',minval(data_area),' max:',maxval(data_area),i1,i2
        IMFloop: do iIMF=2,resultlen2
          ! Get minima-maxima
          roundloop: do iround=1,max_rounds
            extremas = -1
            crossings = -1
            !do it=1,tlen
            !  write(*,'(100f7.3)') real(it), data_area(it), maxima(it), minima(it)
            !end do
            !write(*,*) '2 min:',minval(data_area),' max:',maxval(data_area),i1,i2,iIMF,iround

            call calc_envelopes()
           ! write(*,*) '3 min:',minval(data_area),' max:',maxval(data_area),i1,i2,iIMF,iround

            if (.not. lIMF) then
              exit roundloop
            end if
            call calc_imf()
            if (lIMF) then
              exit roundloop
            end if
          end do roundloop
          imf(1:tlen) = 0
          imf(1:tlen) = data_area(1:tlen)
          analysis(i1,i2,1:tlen,iIMF) = analysis(i1,i2,1:tlen,iIMF-1)-imf(1:tlen)
          data_area(1:tlen) = analysis(i1,i2,1:tlen,iIMF)
          analysis(i1,i2,1:tlen,iIMF) = imf(1:tlen)
          if (.not. lIMF) then
            write(*,*) 'Warning, IMF ',iIMF-1,' at (',i1,',',i2,') was not good'
            write(*,*) 'nmaxima:',nmaxima,' nminima:',nminima
            write(*,*) 'extrema:',extremas,' crossings:',crossings
            write(*,*) 'imf min:', minval(imf), maxval(imf)
            exit IMFloop
          end if
        end do IMFloop
        exit jloop
      end do iloop
    end do jloop


    !write(*,*) 'Whole results'
    !do iIMF=0,maxIMF-1
    !  write(*,*) 'IMF',iIMF
    !  do i1=1,tlen
    !    do i2=1,xdim2
    !      write(*,'(100f7.3)') analysis(1:xdim1,i2,i1+iIMF*tlen)
    !    end do
    !  end do
    !end do


  contains

  subroutine calc_envelopes
    implicit none
    integer           :: t,t_last
    real(kind=8)       :: y,y0,y1,y2
    
    maxima      = 0
    minima      = 0
    smaxima     = 0
    sminima     = 0 
    ymaxima     = 0
    yminima     = 0
    tminima     = 0
    tmaxima     = 0
    nminima     = 0
    nmaxima     = 0

    t_last = 0

    if (data_area(1)>data_area(2)) then
        nmaxima = nmaxima + 1
        ymaxima(nmaxima) = data_area(1)
        tmaxima(nmaxima) = 1.0
    else
        nminima = nminima + 1
        yminima(nminima) = data_area(1)
        tminima(nminima) = 1.0
    end if

    do t=2,tlen-1
      y0 = data_area(t-1)
      y1 = data_area(t)
      y2 = data_area(t+1)
      if ((y0 - y1 >= 0) .and. (y2 - y1 >= 0)) then
        nminima = nminima + 1
        yminima(nminima) = y1
        tminima(nminima) = real(t,kind=8)
        !nmaxima = nmaxima + 1
        !ymaxima(nmaxima) = 0.5*y1+0.5*ymaxima(nmaxima-1)
        !tmaxima(nmaxima) = real(t)
      else if ((y0 - y1 <= 0) .and. (y2 - y1 <= 0)) then
        nmaxima = nmaxima + 1
        ymaxima(nmaxima) = y1
        tmaxima(nmaxima) = real(t,kind=8)
        !nminima = nminima + 1
        !yminima(nminima) = 0.5*y1+0.5*yminima(nminima-1)
        !tminima(nminima) = real(t)
      end if
    end do

    if (data_area(tlen)>data_area(tlen-1)) then
        nmaxima = nmaxima + 1
        ymaxima(nmaxima) = data_area(tlen)
        tmaxima(nmaxima) = real(tlen)
    else
        nminima = nminima + 1
        yminima(nminima) = data_area(tlen)
        tminima(nminima) = real(tlen)
    end if

    lIMF = ((nmaxima > 3) .and. (nminima > 3))

    if (lIMF) then
      call spline_cubic_set(nmaxima,tmaxima(1:nmaxima),ymaxima(1:nmaxima),2,0.0,2,0.0,smaxima(1:nmaxima))
      do t=1,tlen
        y=real(t)
        call spline_cubic_val(nmaxima,tmaxima(1:nmaxima),ymaxima(1:nmaxima),smaxima(1:nmaxima),y,maxima(t),y1,y2)
      end do
      call spline_cubic_set(nminima,tminima(1:nminima),yminima(1:nminima),2,0.0,2,0.0,sminima(1:nminima))
      do t=1,tlen
        y=real(t)
        call spline_cubic_val(nminima,tminima(1:nminima),yminima(1:nminima),sminima(1:nminima),y,minima(t),y1,y2)
      end do
    end if

    !if (iround == 1) then
     ! write(*,*) '5-data_area min:',minval(data_area),' max:',maxval(data_area),i1,i2
      !write(*,*) '5-yminima min:',minval(yminima),' max:',maxval(yminima),i1,i2
      !write(*,*) '5-ymaxima min:',minval(ymaxima),' max:',maxval(ymaxima),i1,i2
      !write(*,*) '5-minima min:',minval(minima),' max:',maxval(minima),i1,i2
      !write(*,*) '5-maxima min:',minval(maxima),' max:',maxval(maxima),i1,i2
      !write(*,*) '5-extremas nmaxima:',nmaxima,' nminima:',nminima,i1,i2
    !end if

  end subroutine calc_envelopes

  subroutine calc_imf
    implicit none
    integer :: t
    real(kind=8)  :: y0,y1,y2
    
    max_mean  = 0.0
    imf_avg   = 0.0
    imf_std   = 0.0
    crossings = 0
    extremas  = 0

    do t=1,tlen
      y0 = data_area(t)
      y1 = maxima(t)
      y2 = minima(t)
      if ((t==1) .or. (t==tlen)) then
        if (y1 < y0) then
          y1 = y0
        end if
        if (y2 >y0) then
          y2 = y0
        end if
      end if
      env_mean = 0.5*(y1+y2)
      if (abs(env_mean) > max_mean) then
        max_mean = abs(env_mean)
      end if
      imf_avg = imf_avg + y0
      imf_std = imf_std + y0*y0
      data_area(t) = y0 - env_mean
    end do

    imf_avg = imf_avg/real(tlen)
    imf_std = sqrt(imf_std/real(tlen) - imf_avg*imf_avg)

    do t=2,tlen-1
      y0 = data_area(t-1)
      y1 = data_area(t)
      y2 = data_area(t+1)
      if ((y0 - y1 > 0) .and. (y2 - y1 > 0)) then
        extremas = extremas + 1
      else if ((y0 - y1 < 0) .and. (y2 - y1 < 0)) then
        extremas = extremas + 1
      end if
      if ((y0 < 0) .and. (y1 >= 0)) then
        crossings = crossings + 1 
      else if ((y0 > 0) .and. (y1 <= 0)) then
        crossings = crossings + 1
      end if
    end do
    !do t=1,tlen
    !  write(*,'(100f7.3)') real(t), data_area(t)
    !end do

    lIMF = ((max_mean < err_coeff*imf_std) .and. (abs(extremas-crossings) <= 1))
    write(*,*) lIMF, max_mean, err_coeff*imf_std, extremas, crossings

  end subroutine calc_imf

  end function analyzer_emd

end module EMD
