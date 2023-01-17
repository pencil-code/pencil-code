subroutine cfftb ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTB computes the backward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    CFFTB computes a complex periodic sequence from its Fourier coefficients.
!
!    A call of CFFTF followed by a call of CFFTB will multiply the
!    sequence by N.  In other words, the transforms are not normalized.
!
!    The array WSAVE must be initialized by CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N )
!        C_in(K) * exp ( sqrt ( - 1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none
!
  integer n
!
  complex c(n)
  real wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftb1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftb1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTB1 is a lower-level routine used by CFFTB.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  complex c(n)
  complex ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nac
  integer nf
  real wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passb4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passb4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passb2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passb2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passb3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passb3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passb5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passb5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passb ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passb ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cfftb_2d ( ldf, n, f, wsave )
!
!*******************************************************************************
!
!! CFFTB_2D computes a backward two dimensional complex fast Fourier transform.
!
!
!  Discussion:
!
!    The routine computes the backward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to CFFTB_2D followed by a call
!    to CFFTF_2D will return the original data multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) * (-1.)**(I+J),  I,J =0,...,N-1.
!
!    Before calling CFFTF_2D or CFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling CFFTI.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Modified:
!
!    12 March 2001
!
!  Parameters:
!
!    Input, integer LDF, the leading dimension of the matrix.
!
!    Input, integer N, the number of rows and columns in the matrix.
!
!    Input/output, complex F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real WSAVE(4*N+15), a work array whose values depend on N,
!    and which must be initialized by calling CFFTI.
!
  implicit none
!
  integer ldf
  integer n
!
  complex f(ldf,n)
  integer i
  real wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call cfftb ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call cfftb ( n, f(1,i), wsave )
  end do

  return
end
subroutine cfftf ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTF computes the forward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    CFFTF computes the Fourier coefficients of a complex periodic sequence.
!
!    The transform is not normalized.  To obtain a normalized transform,
!    the output must be divided by N.  Otherwise a call of CFFTF
!    followed by a call of CFFTB will multiply the sequence by N.
!
!    The array WSAVE must be initialized by calling CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N )
!        C_in(K) * exp ( - sqrt ( -1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none
!
  integer n
!
  complex c(n)
  real wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftf1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTF1 is a lower level routine used by CFFTF.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  complex c(n)
  complex ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nac
  integer nf
  real wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passf4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passf2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passf2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passf3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passf3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passf5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passf5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passf ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passf ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cfftf_2d ( ldf, n, f, wsave )
!
!*******************************************************************************
!
!! CFFTF_2D computes a two dimensional complex fast Fourier transform.
!
!
!  Discussion:
!
!    The routine computes the forward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to CFFTF_2D,
!    followed by a call to CFFTB_2D will return the original data
!    multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) *(-1.)**(I+J),  I,J =0,...,N-1.
!
!    Before calling CFFTF_2D or CFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling CFFTI.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Modified:
!
!    12 March 2001
!
!  Parameters:
!
!    Input, integer LDF, the leading dimension of the matrix.
!
!    Input, integer N, the number of rows and columns in the matrix.
!
!    Input/output, complex F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real WSAVE(4*N+15), a work array whose values depend on N,
!    and which must be initialized by calling CFFTI.
!
  implicit none
!
  integer ldf
  integer n
!
  complex f(ldf,n)
  integer i
  real wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call cfftf ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call cfftf ( n, f(1,i), wsave )
  end do

  return
end
subroutine cffti ( n, wsave )
!
!*******************************************************************************
!
!! CFFTI initializes WSAVE, used in CFFTF and CFFTB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Output, real WSAVE(4*N+15), contains data, dependent on the value
!    of N, which is necessary for the CFFTF or CFFTB routines.
!
  implicit none
!
  integer n
!
  real wsave(4*n+15)
!
  wsave = 0.0

  if ( n <= 1 ) then
    return
  end if

  call cffti1 ( n, wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! CFFTI1 is a lower level routine used by CFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  real arg
  real argh
  real argld
  real fi
  integer i
  integer i1
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer j
  integer k1
  integer l1
  integer l2
  integer ld
  integer nf
  real r_pi
  real wa(2*n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0E+00 * r_pi ( ) / real ( n )
  i = 2
  l1 = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      i1 = i
      wa(i-1) = 1.0E+00
      wa(i) = 0.0E+00
      ld = ld + l1
      fi = 0.0E+00
      argld = real ( ld ) * argh

      do ii = 4, 2*ido+2, 2
        i = i + 2
        fi = fi + 1.0E+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      if ( ip > 5 ) then
        wa(i1-1) = wa(i-1)
        wa(i1) = wa(i)
      end if

    end do

    l1 = l2

  end do

  return
end
subroutine cosqb ( n, x, wsave )
!
!*******************************************************************************
!
!! COSQB computes the fast cosine transform of quarter wave data.
!
!
!  Discussion:
!
!    COSQB computes a sequence from its representation in terms of a cosine
!    series with odd wave numbers.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!
!        4 * X_in(K) * cos ( ( 2 * K - 1 ) * ( I - 1 ) * PI / ( 2 * N ) )
!
!    COSQB is the unnormalized inverse of COSQF since a call of COSQB
!    followed by a call of COSQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the algorithm.  The WSAVE array must be initialized by
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none
!
  integer n
!
  real, parameter :: tsqrt2 = 2.82842712474619E+00
  real wsave(3*n+15)
  real x(n)
  real x1
!
  if ( n < 2 ) then
    x(1) = 4.0E+00 * x(1)
  else if ( n == 2 ) then
    x1 = 4.0E+00 * ( x(1) + x(2) )
    x(2) = tsqrt2 * ( x(1) - x(2) )
    x(1) = x1
  else
    call cosqb1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqb1 ( n, x, w, xh )
!
!*******************************************************************************
!
!! COSQB1 is a lower level routine used by COSQB.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array.
!
!    Input/output, real X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real W(N).
!
!    Input, real XH(2*N+15).
!
  implicit none
!
  integer n
!
  integer i
  integer k
  integer kc
  integer ns2
  real w(n)
  real x(n)
  real xh(2*n+15)
  real xim1
!
  ns2 = ( n + 1 ) / 2

  do i = 3, n, 2
    xim1 = x(i-1) + x(i)
    x(i) = x(i) - x(i-1)
    x(i-1) = xim1
  end do

  x(1) = x(1) + x(1)

  if ( mod ( n, 2 ) == 0 ) then
    x(n) = 2.0E+00 * x(n)
  end if

  call rfftb ( n, x, xh )

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = w(k-1) * x(kc) + w(kc-1) * x(k)
    xh(kc) = w(k-1) * x(k) - w(kc-1) * x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * ( x(ns2+1) + x(ns2+1) )
  end if

  do k = 2, ns2
    kc = n + 2 - k
    x(k) = xh(k) + xh(kc)
    x(kc) = xh(k) - xh(kc)
  end do

  x(1) = 2.0E+00 * x(1)

  return
end
subroutine cosqf ( n, x, wsave )
!
!*******************************************************************************
!
!! COSQF computes the fast cosine transform of quarter wave data.
!
!
!  Discussion:
!
!    COSQF computes the coefficients in a cosine series representation
!    with only odd wave numbers.
!
!    COSQF is the unnormalized inverse of COSQB since a call of COSQF
!    followed by a call of COSQB will multiply the input sequence X
!    by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + sum ( 2 <= K <= N )
!
!        2 * X_in(K) * cos ( ( 2 * I - 1 ) * ( K - 1 ) * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the algorithm.  The WSAVE array must be initialized by
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N.
!
  implicit none
!
  integer n
!
  real, parameter :: sqrt2 = 1.4142135623731E+00
  real tsqx
  real wsave(3*n+15)
  real x(n)
!
  if ( n < 2 ) then

  else if ( n == 2 ) then
    tsqx = sqrt2 * x(2)
    x(2) = x(1) - tsqx
    x(1) = x(1) + tsqx
  else
    call cosqf1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqf1 ( n, x, w, xh )
!
!*******************************************************************************
!
!! COSQF1 is a lower level routine used by COSQF.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Input/output, real X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real W(N).
!
!    Input, real XH(2*N+15).
!
  implicit none
!
  integer n
!
  integer i
  integer k
  integer kc
  integer ns2
  real w(n)
  real x(n)
  real xh(2*n+15)
  real xim1
!
  ns2 = ( n + 1 ) / 2

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = x(k) + x(kc)
    xh(kc) = x(k) - x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  do k = 2, ns2
    kc = n+2-k
    x(k) = w(k-1) * xh(kc) + w(kc-1) * xh(k)
    x(kc) = w(k-1) * xh(k) - w(kc-1) * xh(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * xh(ns2+1)
  end if

  call rfftf ( n, x, xh )

  do i = 3, n, 2
    xim1 = x(i-1) - x(i)
    x(i) = x(i-1) + x(i)
    x(i-1) = xim1
  end do

  return
end
subroutine cosqi ( n, wsave )
!
!*******************************************************************************
!
!! COSQI initializes WSAVE, used in COSQF and COSQB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The method
!    is more efficient when N is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the COSQB and COSQF algorithms.
!
  implicit none
!
  integer n
!
  real dt
  integer k
  real r_pi
  real wsave(3*n+15)
!
  dt = 0.5E+00 * r_pi ( ) / real ( n )

  do k = 1, n
    wsave(k) = cos ( real ( k ) * dt )
  end do

  call rffti ( n, wsave(n+1) )

  return
end
subroutine cvec_print_some ( n, x, max_print, title )
!
!*******************************************************************************
!
!! CVEC_PRINT_SOME prints some of a complex vector.
!
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, complex X(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none
!
  integer n
!
  integer i
  integer max_print
  character ( len = * ) title
  complex x(n)
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,2g14.6)' ) i, x(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,2g14.6)' ) i, x(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i6,2x,2g14.6)' ) i, x(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,2g14.6)' ) i, x(i)
    end do
    i = max_print
    write ( *, '(i6,2x,2g14.6,2x,a)' ) i, x(i), '...more entries...'

  end if

  return
end
subroutine cvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! CVEC_RANDOM returns a random complex vector in a given range.
!
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, complex A(N), the vector of randomly chosen values.
!
  implicit none
!
  integer n
!
  complex a(n)
  real ahi
  real ai
  real alo
  real ar
  integer i
!
  do i = 1, n

    call r_random ( alo, ahi, ar )
    call r_random ( alo, ahi, ai )

    a(i) = cmplx ( ar, ai )

  end do

  return
end
function d_pi ()
!
!*******************************************************************************
!
!! D_PI returns the value of pi.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision D_PI, the value of PI.
!
  implicit none
!
  double precision d_pi
!
  d_pi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
subroutine d_random ( dlo, dhi, d )
!
!*******************************************************************************
!
!! D_RANDOM returns a random double precision value in a given range.
!
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision DLO, DHI, the minimum and maximum values.
!
!    Output, double precision D, the randomly chosen value.
!
  implicit none
!
  double precision d
  double precision dhi
  double precision dlo
  double precision t
!
!  Pick T, a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set D in ( DLO, DHI ).
!
  d = ( 1.0D+00 - dble ( t ) ) * dlo + dble ( t ) * dhi

  return
end
subroutine dadf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! DADF2 is a lower level routine used by DFFTF1.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  double precision cc(ido,l1,2)
  double precision ch(ido,2,l1)
  integer i
  integer ic
  integer k
  double precision ti2
  double precision tr2
  double precision wa1(ido)
!
  ch(1,1,1:l1)   = cc(1,1:l1,1) + cc(1,1:l1,2)
  ch(ido,2,1:l1) = cc(1,1:l1,1) - cc(1,1:l1,2)

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        tr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ti2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)

        ch(i,1,k) = cc(i,k,1) + ti2
        ch(ic,2,k) = ti2 - cc(i,k,1)
        ch(i-1,1,k) = cc(i-1,k,1) + tr2
        ch(ic-1,2,k) = cc(i-1,k,1) - tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(1,2,1:l1) = -cc(ido,1:l1,2)
  ch(ido,1,1:l1) = cc(ido,1:l1,1)

  return
end
subroutine dadf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! DADF3 is a lower level routine used by DFFTF1.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  double precision cc(ido,l1,3)
  double precision ch(ido,3,l1)
  double precision ci2
  double precision cr2
  double precision di2
  double precision di3
  double precision dr2
  double precision dr3
  integer i
  integer ic
  integer k
  double precision taui
  double precision, parameter :: taur = -0.5D+00
  double precision ti2
  double precision ti3
  double precision tr2
  double precision tr3
  double precision wa1(ido)
  double precision wa2(ido)
!
  taui = sqrt ( 3.0D+00 ) / 2.0D+00

  do k = 1, l1
    cr2 = cc(1,k,2) + cc(1,k,3)
    ch(1,1,k) = cc(1,k,1) + cr2
    ch(1,3,k) = taui * ( cc(1,k,3) - cc(1,k,2) )
    ch(ido,2,k) = cc(1,k,1) + taur * cr2
  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)

      cr2 = dr2 + dr3
      ci2 = di2 + di3

      ch(i-1,1,k) = cc(i-1,k,1) + cr2
      ch(i,1,k)   = cc(i,k,1) + ci2

      tr2 = cc(i-1,k,1) + taur * cr2
      ti2 = cc(i,k,1) + taur * ci2
      tr3 = taui * ( di2 - di3 )
      ti3 = taui * ( dr3 - dr2 )

      ch(i-1,3,k) = tr2 + tr3
      ch(ic-1,2,k) = tr2 - tr3
      ch(i,3,k) = ti2 + ti3
      ch(ic,2,k) = ti3 - ti2

    end do
  end do

  return
end
subroutine dadf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! DADF4 is a lower level routine used by DFFTF1.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  double precision cc(ido,l1,4)
  double precision ch(ido,4,l1)
  double precision ci2
  double precision ci3
  double precision ci4
  double precision cr2
  double precision cr3
  double precision cr4
  double precision hsqt2
  integer i
  integer ic
  integer k
  double precision ti1
  double precision ti2
  double precision ti3
  double precision ti4
  double precision tr1
  double precision tr2
  double precision tr3
  double precision tr4
  double precision wa1(ido)
  double precision wa2(ido)
  double precision wa3(ido)
!
  hsqt2 = sqrt ( 2.0D+00 ) / 2.0D+00

  do k = 1, l1
    tr1 = cc(1,k,2) + cc(1,k,4)
    tr2 = cc(1,k,1) + cc(1,k,3)
    ch(1,1,k) = tr1 + tr2
    ch(ido,4,k) = tr2 - tr1
    ch(ido,2,k) = cc(1,k,1) - cc(1,k,3)
    ch(1,3,k) = cc(1,k,4) - cc(1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        cr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ci2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
        cr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
        ci3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
        cr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
        ci4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)

        tr1 = cr2 + cr4
        tr4 = cr4 - cr2
        ti1 = ci2 + ci4
        ti4 = ci2 - ci4
        ti2 = cc(i,k,1) + ci3
        ti3 = cc(i,k,1) - ci3
        tr2 = cc(i-1,k,1) + cr3
        tr3 = cc(i-1,k,1) - cr3

        ch(i-1,1,k)  = tr1 + tr2
        ch(ic-1,4,k) = tr2 - tr1
        ch(i,1,k)    = ti1 + ti2
        ch(ic,4,k)   = ti1 - ti2
        ch(i-1,3,k)  = ti4 + tr3
        ch(ic-1,2,k) = tr3 - ti4
        ch(i,3,k)    = tr4 + ti3
        ch(ic,2,k)   = tr4 - ti3

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = -hsqt2 * ( cc(ido,k,2) + cc(ido,k,4) )
    tr1 =  hsqt2 * ( cc(ido,k,2) - cc(ido,k,4) )

    ch(ido,1,k) = tr1 + cc(ido,k,1)
    ch(ido,3,k) = cc(ido,k,1) - tr1

    ch(1,2,k) = ti1 - cc(ido,k,3)
    ch(1,4,k) = ti1 + cc(ido,k,3)

  end do

  return
end
subroutine dadf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! DADF5 is a lower level routine used by DFFTF1.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  double precision cc(ido,l1,5)
  double precision ch(ido,5,l1)
  double precision ci2
  double precision ci3
  double precision ci4
  double precision ci5
  double precision cr2
  double precision cr3
  double precision cr4
  double precision cr5
  double precision di2
  double precision di3
  double precision di4
  double precision di5
  double precision dr2
  double precision dr3
  double precision dr4
  double precision dr5
  integer i
  integer ic
  integer k
  double precision, parameter :: ti11 =  0.951056516295154D+00
  double precision, parameter :: ti12 =  0.587785252292473D+00
  double precision ti2
  double precision ti3
  double precision ti4
  double precision ti5
  double precision, parameter :: tr11 =  0.309016994374947D+00
  double precision, parameter :: tr12 = -0.809016994374947D+00
  double precision tr2
  double precision tr3
  double precision tr4
  double precision tr5
  double precision wa1(ido)
  double precision wa2(ido)
  double precision wa3(ido)
  double precision wa4(ido)
!
  do k = 1, l1

    cr2 = cc(1,k,5) + cc(1,k,2)
    ci5 = cc(1,k,5) - cc(1,k,2)
    cr3 = cc(1,k,4) + cc(1,k,3)
    ci4 = cc(1,k,4) - cc(1,k,3)

    ch(1,1,k)   = cc(1,k,1) + cr2 + cr3
    ch(ido,2,k) = cc(1,k,1) + tr11 * cr2 + tr12 * cr3
    ch(1,3,k)   = ti11 * ci5 + ti12 * ci4
    ch(ido,4,k) = cc(1,k,1) + tr12 * cr2 + tr11 * cr3
    ch(1,5,k)   = ti12 * ci5 - ti11 * ci4

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
      dr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
      di4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)
      dr5 = wa4(i-2) * cc(i-1,k,5) + wa4(i-1) * cc(i,k,5)
      di5 = wa4(i-2) * cc(i,k,5)   - wa4(i-1) * cc(i-1,k,5)

      cr2 = dr2 + dr5
      ci5 = dr5 - dr2
      cr5 = di2 - di5
      ci2 = di2 + di5
      cr3 = dr3 + dr4
      ci4 = dr4 - dr3
      cr4 = di3 - di4
      ci3 = di3 + di4

      ch(i-1,1,k) = cc(i-1,k,1) + cr2 + cr3
      ch(i,1,k)   = cc(i,k,1)   + ci2 + ci3

      tr2 = cc(i-1,k,1) + tr11 * cr2 + tr12 * cr3
      ti2 = cc(i,k,1)   + tr11 * ci2 + tr12 * ci3
      tr3 = cc(i-1,k,1) + tr12 * cr2 + tr11 * cr3
      ti3 = cc(i,k,1)   + tr12 * ci2 + tr11 * ci3

      tr5 = ti11 * cr5 + ti12 * cr4
      ti5 = ti11 * ci5 + ti12 * ci4
      tr4 = ti12 * cr5 - ti11 * cr4
      ti4 = ti12 * ci5 - ti11 * ci4

      ch(i-1,3,k)  = tr2 + tr5
      ch(ic-1,2,k) = tr2 - tr5
      ch(i,3,k)    = ti2 + ti5
      ch(ic,2,k)   = ti5 - ti2
      ch(i-1,5,k)  = tr3 + tr4
      ch(ic-1,4,k) = tr3 - tr4
      ch(i,5,k)    = ti3 + ti4
      ch(ic,4,k)   = ti4 - ti3

    end do
  end do

  return
end
subroutine dadfg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! DADFG is a lower level routine used by DFFTF1.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  double precision ai1
  double precision ai2
  double precision ar1
  double precision ar1h
  double precision ar2
  double precision ar2h
  double precision arg
  double precision c1(ido,l1,ip)
  double precision c2(idl1,ip)
  double precision cc(ido,ip,l1)
  double precision ch(ido,l1,ip)
  double precision ch2(idl1,ip)
  double precision d_pi
  double precision dc2
  double precision dcp
  double precision ds2
  double precision dsp
  integer i
  integer ic
  integer idij
  integer ik
  integer ipph
  integer is
  integer j
  integer j2
  integer jc
  integer k
  integer l
  integer lc
  integer nbd
  double precision wa(*)
!
  arg = 2.0D+00 * d_pi ( ) / dble ( ip )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then

    c2(1:idl1,1) = ch2(1:idl1,1)

  else

    ch2(1:idl1,1) = c2(1:idl1,1)
    ch(1,1:l1,2:ip) = c1(1,1:l1,2:ip)

    if ( nbd <= l1 ) then

      is = -ido
      do j = 2, ip
        is = is + ido
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          do k = 1, l1
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    else

      is = -ido

      do j = 2, ip
        is = is + ido
        do k = 1, l1
          idij = is
          do i = 3, ido, 2
            idij = idij + 2
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    end if

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            c1(i-1,k,j)  = ch(i-1,k,j)  + ch(i-1,k,jc)
            c1(i-1,k,jc) = ch(i,k,j)    - ch(i,k,jc)
            c1(i,k,j)    = ch(i,k,j)    + ch(i,k,jc)
            c1(i,k,jc)   = ch(i-1,k,jc) - ch(i-1,k,j)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          c1(i-1,1:l1,j)  = ch(i-1,1:l1,j)  + ch(i-1,1:l1,jc)
          c1(i-1,1:l1,jc) = ch(i,1:l1,j)    - ch(i,1:l1,jc)
          c1(i,1:l1,j)    = ch(i,1:l1,j)    + ch(i,1:l1,jc)
          c1(i,1:l1,jc)   = ch(i-1,1:l1,jc) - ch(i-1,1:l1,j)
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ip + 2 - j
    c1(1,1:l1,j)  = ch(1,1:l1,j)  + ch(1,1:l1,jc)
    c1(1,1:l1,jc) = ch(1,1:l1,jc) - ch(1,1:l1,j)
  end do

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      ch2(ik,l) = c2(ik,1) + ar1 * c2(ik,2)
      ch2(ik,lc) =           ai1 * c2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2 =  dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        ch2(ik,l) =  ch2(ik,l)  + ar2 * c2(ik,j)
        ch2(ik,lc) = ch2(ik,lc) + ai2 * c2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + c2(1:idl1,j)
  end do

  cc(1:ido,1,1:l1) = ch(1:ido,1:l1,1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    cc(ido,j2-2,1:l1) = ch(1,1:l1,j)
    cc(1,j2-1,1:l1)   = ch(1,1:l1,jc)
  end do

  if ( ido == 1 ) then
    return
  end if

  if ( nbd >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do k = 1, l1
        do i = 3, ido, 2
          ic = ido + 2 - i
          cc(i-1,j2-1,k)  = ch(i-1,k,j) + ch(i-1,k,jc)
          cc(ic-1,j2-2,k) = ch(i-1,k,j) - ch(i-1,k,jc)
          cc(i,j2-1,k)    = ch(i,k,j)   + ch(i,k,jc)
          cc(ic,j2-2,k)   = ch(i,k,jc)  - ch(i,k,j)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do i = 3, ido, 2
        ic = ido + 2 - i
        cc(i-1,j2-1,1:l1)  = ch(i-1,1:l1,j) + ch(i-1,1:l1,jc)
        cc(ic-1,j2-2,1:l1) = ch(i-1,1:l1,j) - ch(i-1,1:l1,jc)
        cc(i,j2-1,1:l1)    = ch(i,1:l1,j)   + ch(i,1:l1,jc)
        cc(ic,j2-2,1:l1)   = ch(i,1:l1,jc)  - ch(i,1:l1,j)
      end do
    end do

  end if

  return
end
subroutine dcost ( n, x, wsave )
!
!*******************************************************************************
!
!! DCOST computes the discrete Fourier cosine transform of an even sequence.
!
!
!  Discussion:
!
!    This routine is the unnormalized inverse of itself.  Two successive
!    calls will multiply the input sequence X by 2*(N-1).
!
!    The array WSAVE must be initialized by calling DCOSTI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + (-1) **(I-1) * X_in(N) + sum ( 2 <= K <= N-1 )
!
!        2 * X_in(K) * cos ( ( K - 1 ) * ( I - 1 ) * PI / ( N - 1 ) )
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  The
!    method is more efficient when N-1 is the product of small primes.
!
!    Input/output, double precision X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, double precision WSAVE(3*N+15).
!    The WSAVE array must be initialized by calling DCOSTI.  A different
!    array must be used for each different value of N.
!
  implicit none
!
  integer n
!
  double precision c1
  integer i
  integer k
  integer kc
  integer ns2
  double precision t1
  double precision t2
  double precision tx2
  double precision wsave(3*n+15)
  double precision x(n)
  double precision x1h
  double precision x1p3
  double precision xi
  double precision xim2
!
  ns2 = n / 2

  if ( n <= 1 ) then
    return
  end if

  if ( n == 2 ) then
    x1h = x(1) + x(2)
    x(2) = x(1) - x(2)
    x(1) = x1h
    return
  end if

  if ( n == 3 ) then
    x1p3 = x(1) + x(3)
    tx2 = x(2) + x(2)
    x(2) = x(1) - x(3)
    x(1) = x1p3 + tx2
    x(3) = x1p3 - tx2
    return
  end if

  c1 = x(1) - x(n)
  x(1) = x(1) + x(n)

  do k = 2, ns2
    kc = n + 1 - k
    t1 = x(k) + x(kc)
    t2 = x(k) - x(kc)
    c1 = c1 + wsave(kc) * t2
    t2 = wsave(k) * t2
    x(k) = t1 - t2
    x(kc) = t1 + t2
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  call dfftf ( n-1, x, wsave(n+1) )

  xim2 = x(2)
  x(2) = c1

  do i = 4, n, 2
    xi = x(i)
    x(i) = x(i-2) - x(i-1)
    x(i-1) = xim2
    xim2 = xi
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n) = xim2
  end if

  return
end
subroutine dcosti ( n, wsave )
!
!*******************************************************************************
!
!! DCOSTI initializes WSAVE, used in DCOST.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  The
!    method is more efficient when N-1 is the product of small primes.
!
!    Output, double precision WSAVE(3*N+15), contains data, depending on N, and
!    required by the DCOST algorithm.
!
  implicit none
!
  integer n
!
  double precision d_pi
  double precision dt
  integer k
  double precision wsave(3*n+15)
!
  if ( n <= 3 ) then
    return
  end if

  dt = d_pi ( ) / dble ( n - 1 )

  do k = 2, ( n / 2 )
    wsave(k)     = 2.0D+00 * sin ( dble ( k - 1 ) * dt )
    wsave(n+1-k) = 2.0D+00 * cos ( dble ( k - 1 ) * dt )
  end do

  call dffti ( n-1, wsave(n+1) )

  return
end
subroutine dfftf ( n, r, wsave )
!
!*******************************************************************************
!
!! DFFTF computes the Fourier coefficients of a real periodic sequence.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    The transform is unnormalized.  A call to DFFTF followed by a call
!    to DFFTB will multiply the input sequence by N.
!
!    The transform is defined by:
!
!      R_out(1) = sum ( 1 <= I <= N ) R_in(I)
!
!    Letting L = (N+1)/2, then for K = 2,...,L
!
!      R_out(2*K-2) = sum ( 1 <= I <= N )
!
!        R_in(I) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!      R_out(2*K-1) = sum ( 1 <= I <= N )
!
!        -R_in(I) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    And, if N is even, then:
!
!      R_out(N) = sum ( 1 <= I <= N ) (-1)**(I-1) * R_in(I)
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, double precision R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, double precision WSAVE(2*N+15), a work array.  The WSAVE array
!    must be initialized by calling DFFTI.  A different WSAVE array must be
!    used for each different value of N.
!
  implicit none
!
  integer n
!
  double precision r(n)
  double precision wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call dfftf1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine dfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! DFFTF1 is a lower level routine used by DFFTF and DSINT.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Input/output, double precision C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, double precision CH(N).
!
!    Input, double precision WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  double precision c(n)
  double precision ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer kh
  integer l1
  integer l2
  integer na
  integer nf
  double precision wa(n)
!
  nf = ifac(2)
  na = 1
  l2 = n
  iw = n

  do k1 = 1, nf

    kh = nf - k1
    ip = ifac(kh+3)
    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call dadf4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call dadf4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call dadf2 ( ido, l1, c, ch, wa(iw) )
      else
        call dadf2 ( ido, l1, ch, c, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call dadf3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call dadf3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call dadf5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call dadf5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call dadfg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
        na = 1
      else
        call dadfg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  if ( na /= 1 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine dffti ( n, wsave )
!
!*******************************************************************************
!
!! DFFTI initializes WSAVE, used in DFFTF and DFFTB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Output, double precision WSAVE(2*N+15), contains data, dependent
!    on the value of N, which is necessary for the DFFTF and DFFTB routines.
!
  implicit none
!
  integer n
!
  double precision wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call dffti1 ( n, wsave(n+1), wsave(2*n+1) )

  return
end
subroutine dffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! DFFTI1 is a lower level routine used by DFFTI.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input, double precision WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  double precision arg
  double precision argh
  double precision argld
  double precision d_pi
  double precision fi
  integer i
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer is
  integer j
  integer k1
  integer l1
  integer l2
  integer ld
  integer nf
  double precision wa(n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0D+00 * d_pi ( ) / dble ( n )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      ld = ld + l1
      i = is
      argld = dble ( ld ) * argh
      fi = 0.0D+00

      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0D+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine dsct ( n, x, y )
!
!*******************************************************************************
!
!! DSCT computes a double precision "slow" cosine transform.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!      Y(1) = Sum ( 1 <= J <= N ) X(J)
!
!      For I from 2 to N-1:
!
!        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
!          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
!
!      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
!    and accuracy.
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, double precision X(N), the data sequence.
!
!    Output, double precision Y(N), the transformed data.
!
  implicit none
!
  integer n
!
  double precision d_pi
  integer i
  integer j
  double precision theta
  double precision x(n)
  double precision y(n)
!
  y(1) = sum ( x(1:n) )

  do i = 2, n-1
    y(i) = 0.0D+00
    do j = 1, n
      theta = d_pi ( ) * &
        dble ( mod ( ( j - 1 ) * ( i - 1 ), 2 * ( n - 1 ) ) ) / dble ( n - 1 )
!     theta = dble ( j - 1 ) * dble ( i - 1 ) * d_pi() / dble ( n - 1 )
      y(i) = y(i) + 2.0D+00 * x(j) * cos ( theta )
    end do
  end do

  y(n) = sum ( x(1:n:2) ) - sum ( x(2:n:2) )

  return
end
subroutine dsint ( n, x, wsave )
!
!*******************************************************************************
!
!! DSINT computes the discrete Fourier sine transform of an odd sequence.
!
!
!  Discussion:
!
!    This routine is the unnormalized inverse of itself since two successive
!    calls will multiply the input sequence X by 2*(N+1).
!
!    The array WSAVE must be initialized by calling DSINTI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!        2 * X_in(K) * sin ( K * I * PI / ( N + 1 ) )
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, double precision WSAVE((5*N+30)/2), a work array.  The WSAVE
!    array must be initialized by calling DSINTI.  A different WSAVE array
!    must be used for each different value of N.
!
  implicit none
!
  integer n
!
  integer iw1
  integer iw2
  integer iw3
  double precision wsave((5*n+30)/2)
  double precision x(n)
!
  iw1 = n / 2 + 1
  iw2 = iw1 + n + 1
  iw3 = iw2 + n + 1

  call dsint1 ( n, x, wsave(1), wsave(iw1), wsave(iw2), wsave(iw3) )

  return
end
subroutine dsint1 ( n, war, was, xh, x, ifac )
!
!*******************************************************************************
!
!! DSINT1 is a lower level routine used by DSINT.
!
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input/output, double precision WAR(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, double precision WAS(N/2).
!
!    Input, double precision XH(N).
!
!    Input, double precision X(N+1), ?.
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  integer i
  integer ifac(15)
  integer k
  integer ns2
  double precision sqrt3
  double precision t1
  double precision t2
  double precision war(n)
  double precision was(n/2)
  double precision x(n+1)
  double precision xh(n)
  double precision xhold
!
  sqrt3 = sqrt ( 3.0D+00 )
!
  xh(1:n) = war(1:n)
  war(1:n) = x(1:n)

  if ( n <= 1 ) then
    xh(1) = 2.0D+00 * xh(1)
    return
  end if

  if ( n == 2 ) then
    xhold = sqrt3 * ( xh(1) + xh(2) )
    xh(2) = sqrt3 * ( xh(1) - xh(2) )
    xh(1) = xhold
    return
  end if

  ns2 = n / 2
  x(1) = 0.0D+00

  do k = 1, n/2
    t1 = xh(k) - xh(n+1-k)
    t2 = was(k) * ( xh(k) + xh(n+1-k) )
    x(k+1) = t1 + t2
!
! ??? N+2-K puts us out of the array...DAMN IT, THIS IS AN ERROR.
!
    x(n+2-k) = t2 - t1
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n/2+2) = 4.0D+00 * xh(n/2+1)
  end if
!
!  This call says there are N+1 things in X.
!
  call dfftf1 ( n+1, x, xh, war, ifac )

  xh(1) = 0.5D+00 * x(1)
  do i = 3, n, 2
    xh(i-1) = -x(i)
    xh(i) = xh(i-2) + x(i-1)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(n) = -x(n+1)
  end if

  x(1:n) = war(1:n)
  war(1:n) = xh(1:n)

  return
end
subroutine dsinti ( n, wsave )
!
!*******************************************************************************
!
!! DSINTI initializes WSAVE, used in DSINT.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is a product of small primes.
!
!    Output, double precision WSAVE((5*N+30)/2), contains data, dependent
!    on the value of N, which is necessary for the DSINT routine.
!
  implicit none
!
  integer n
!
  double precision d_pi
  double precision dt
  integer k
  double precision wsave((5*n+30)/2)
!
  if ( n <= 1 ) then
    return
  end if

  dt = d_pi ( ) / dble ( n + 1 )

  do k = 1, n/2
    wsave(k) = 2.0D+00 * sin ( dble ( k ) * dt )
  end do

  call dffti ( n+1, wsave((n/2)+1) )

  return
end
subroutine dsst ( n, x, y )
!
!*******************************************************************************
!
!! DSST computes a double precision "slow" sine transform.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For I from 1 to N,
!
!      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by N / 2.  This is a good check for correctness and accuracy.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, double precision X(N), the data sequence.
!
!    Output, double precision Y(N), the transformed data.
!
  implicit none
!
  integer n
!
  integer i
  double precision d_pi
  double precision theta(n)
  double precision x(n)
  double precision y(n)
!
  call dvec_identity ( n, theta )
  theta(1:n) = theta(1:n) * d_pi() / dble ( n + 1 )

  y(1:n) = 0.0D+00

  do i = 1, n
    y(1:n) = y(1:n) + 2.0D+00 * x(i) * sin ( dble ( i ) * theta(1:n) )
  end do

  return
end
subroutine dvec_identity ( n, a )
!
!*******************************************************************************
!
!! DVEC_IDENTITY sets a double precision vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, double precision A(N), the array to be initialized.
!
  implicit none
!
  integer n
!
  double precision a(n)
  integer i
!
  do i = 1, n
    a(i) = dble ( i )
  end do

  return
end
subroutine dvec_print_some ( n, a, max_print, title )
!
!*******************************************************************************
!
!! DVEC_PRINT_SOME prints "some" of a double precision vector.
!
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, double precision A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none
!
  integer n
!
  double precision a(n)
  integer i
  integer max_print
  character ( len = * ) title
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i6,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(i6,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

  end if

  return
end
subroutine dvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! DVEC_RANDOM returns a random double precision vector in a given range.
!
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, double precision A(N), the vector of randomly chosen values.
!
  implicit none
!
  integer n
!
  double precision a(n)
  double precision ahi
  double precision alo
  integer i
!
  do i = 1, n
    call d_random ( alo, ahi, a(i) )
  end do

  return
end
subroutine ezfftb ( n, r, azero, a, b, wsave )
!
!*******************************************************************************
!
!! EZFFTB computes a real periodic sequence from its Fourier coefficients.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    EZFFTB is a simplified but slower version of RFFTB.
!
!    The transform is defined by:
!
!      R(I) = AZERO + sum ( 1 <= K <= N/2 )
!
!          A(K) * cos ( K * ( I - 1 ) * 2 * PI / N )
!        + B(K) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the output array.  The
!    method is more efficient when N is the product of small primes.
!
!    Output, real R(N), the reconstructed data sequence.
!
!    Input, real AZERO, the constant Fourier coefficient.
!
!    Input, real A(N/2), B(N/2), the Fourier coefficients.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling EZFFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
  implicit none
!
  integer n
!
  real a(n/2)
  real azero
  real b(n/2)
  integer i
  integer ns2
  real r(n)
  real wsave(3*n+15)
!
  if ( n < 2 ) then

    r(1) = azero

  else if ( n == 2 ) then

    r(1) = azero + a(1)
    r(2) = azero - a(1)

  else

    ns2 = ( n - 1 ) / 2

    do i = 1, ns2
      r(2*i) =    0.5E+00 * a(i)
      r(2*i+1) = -0.5E+00 * b(i)
    end do

    r(1) = azero

    if ( mod ( n, 2 ) == 0 ) then
      r(n) = a(ns2+1)
    end if

    call rfftb ( n, r, wsave(n+1) )

  end if

  return
end
subroutine ezfftf ( n, r, azero, a, b, wsave )
!
!*******************************************************************************
!
!! EZFFTF computes the Fourier coefficients of a real periodic sequence.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    EZFFTF is a simplified but slower version of RFFTF.
!
!    The transform is defined by:
!
!      AZERO = sum ( 1 <= I <= N ) R(I) / N,
!
!    and, for K = 1 to (N-1)/2,
!
!      A(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * cos ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      A(N/2) = sum ( 1 <= I <= N ) (-1) **(I-1) * R(I) / N
!
!    For K = 1 to (N-1)/2,
!
!      B(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      B(N/2) = 0.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input, real R(N), the sequence to be transformed.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling EZFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
!    Output, real AZERO, the constant Fourier coefficient.
!
!    Output, real A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none
!
  integer n
!
  real a(n/2)
  real azero
  real b(n/2)
  real cf
  integer i
  integer ns2
  real r(n)
  real wsave(3*n+15)
!
  if ( n < 2 ) then

    azero = r(1)

  else if ( n == 2 ) then

    azero = 0.5E+00 * ( r(1) + r(2) )
    a(1) = 0.5E+00 * ( r(1) - r(2) )

  else

    wsave(1:n) = r(1:n)

    call rfftf ( n, wsave(1), wsave(n+1) )

    cf = 2.0E+00 / real ( n )
    azero = 0.5E+00 * cf * wsave(1)
    ns2 = ( n + 1 ) / 2

    do i = 1, ns2-1
      a(i) = cf * wsave(2*i)
      b(i) = -cf * wsave(2*i+1)
    end do

    if ( mod ( n, 2 ) /= 1 ) then
      a(ns2) = 0.5E+00 * cf * wsave(n)
      b(ns2) = 0.0E+00
    end if

  end if

  return
end
subroutine ezffti ( n, wsave )
!
!*******************************************************************************
!
!! EZFFTI initializes WSAVE, used in EZFFTF and EZFFTB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, dependent on the value
!    of N, which is necessary for the EZFFTF or EZFFTB routines.
!
  implicit none
!
  integer n
!
  real wsave(3*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call ezffti1 ( n, wsave(2*n+1), wsave(3*n+1) )

  return
end
subroutine ezffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! EZFFTI1 is a lower level routine used by EZFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Output, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  real arg1
  real argh
  real ch1
  real ch1h
  real dch1
  real dsh1
  integer i
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer is
  integer j
  integer k1
  integer l1
  integer l2
  integer nf
  real r_pi
  real sh1
  real wa(n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0E+00 * r_pi ( ) / real ( n )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    l2 = l1 * ip
    ido = n / l2
    arg1 = real ( l1 ) * argh
    ch1 = 1.0E+00
    sh1 = 0.0E+00
    dch1 = cos ( arg1 )
    dsh1 = sin ( arg1 )

    do j = 1, ip-1

      ch1h = dch1 * ch1 - dsh1 * sh1
      sh1  = dch1 * sh1 + dsh1 * ch1
      ch1 = ch1h
      i = is + 2
      wa(i-1) = ch1
      wa(i) = sh1

      do ii = 5, ido, 2
        i = i + 2
        wa(i-1) = ch1 * wa(i-3) - sh1 * wa(i-2)
        wa(i)   = ch1 * wa(i-2) + sh1 * wa(i-3)
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine i_factor ( n, ifac )
!
!*******************************************************************************
!
!! I_FACTOR factors an integer.
!
!
!  Modified:
!
!    14 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the number to be factored.
!
!    Output, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer i
  integer ifac(15)
  integer j
  integer n
  integer nf
  integer nl
  integer nq
  integer nr
  integer ntry
!
  ifac(1) = n

  nf = 0
  nl = n

  if ( n == 0 ) then
    nf = 1
    ifac(2) = nf
    ifac(2+nf) = 0
    return
  end if

  if ( n < 1 ) then
    nf = nf + 1
    ifac(2+nf) = -1
    nl = - n
  end if

  if ( nl == 1 ) then
    nf = nf + 1
    ifac(2) = nf
    ifac(2+nf) = 1
    return
  end if

  j = 0

  do while ( nl > 1 )

    j = j + 1
!
!  Choose a trial divisor, NTRY.
!
    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if
!
!  Divide by the divisor as many times as possible.
!
    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nl = nq
      nf = nf + 1
!
!  Make sure factors of 2 appear in the front of the list.
!
      if ( ntry /= 2 ) then

        ifac(2+nf) = ntry

      else

        do i = nf, 2, -1
          ifac(i+2) = ifac(i+1)
        end do
        ifac(3) = 2

      end if

    end do

  end do

  ifac(2) = nf

  return
end
subroutine passb ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSB is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  integer i
  integer idij
  integer idj
  integer idl
  integer idlj
  integer idp
  integer ik
  integer inc
  integer ipph
  integer j
  integer jc
  integer k
  integer l
  integer lc
  integer nac
  integer nt
  real wa(*)
  real wai
  real war
!
  nt = ip * idl1
  ipph = ( ip + 1 ) / 2
  idp = ip * ido

  if ( ido >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      do k = 1, l1
        ch(1:ido,k,j)  = cc(1:ido,j,k) + cc(1:ido,jc,k)
        ch(1:ido,k,jc) = cc(1:ido,j,k) - cc(1:ido,jc,k)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      do i = 1, ido
        ch(i,1:l1,j)  = cc(i,j,1:l1) + cc(i,jc,1:l1)
        ch(i,1:l1,jc) = cc(i,j,1:l1) - cc(i,jc,1:l1)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l) = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =            wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j
      idlj = idlj + inc
      if ( idlj > idp ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  nac = 1

  if ( ido == 2 ) then
    return
  end if

  nac = 0
  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1:2,1:l1,2:ip) = ch(1:2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passb2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSB2 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,2,l1)
  real ch(ido,l1,2)
  integer i
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   + cc(i,2,k)
        ti2         = cc(i,1,k)   - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 + wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 - wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passb3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSB3 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,3,l1)
  real ch(ido,l1,3)
  real ci2
  real ci3
  real cr2
  real cr3
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer k
  real taui
  real, parameter :: taur = -0.5E+00
  real ti2
  real tr2
  real wa1(ido)
  real wa2(ido)
!
  taui = sqrt ( 3.0E+00 ) / 2.0E+00

  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k) - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSB4 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,4,l1)
  real ch(ido,l1,4)
  real ci2
  real ci3
  real ci4
  real cr2
  real cr3
  real cr4
  integer i
  integer k
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,4,k) - cc(2,2,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,2,k) - cc(1,4,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k) - cc(i,3,k)
        ti2 = cc(i,1,k) + cc(i,3,k)
        ti3 = cc(i,2,k) + cc(i,4,k)
        tr4 = cc(i,4,k) - cc(i,2,k)

        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,2,k) - cc(i-1,4,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i,k,1) = ti2 + ti3
        ci3 = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 - wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 + wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 - wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 + wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 - wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 + wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSB5 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,5,l1)
  real ch(ido,l1,5)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer k
  real, parameter :: ti11 = 0.951056516295154E+00
  real, parameter :: ti12 = 0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real, parameter :: tr11 = 0.309016994374947E+00
  real, parameter :: tr12 = -0.809016994374947E+00
  real tr2
  real tr3
  real tr4
  real tr5
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 - wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 + wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 - wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 + wa4(i) * dr5

      end do
    end do

  end if

  return
end
subroutine passf ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSF is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  integer i
  integer idij
  integer idj
  integer idl
  integer idlj
  integer idp
  integer ik
  integer inc
  integer ipph
  integer j
  integer jc
  integer k
  integer l
  integer lc
  integer nac
  integer nt
  real wa(*)
  real wai
  real war
!
  nt = ip * idl1
  ipph = (ip+1) / 2
  idp = ip * ido

  if ( ido >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =           - wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j

      idlj = idlj + inc
      if ( idlj > idp ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) - wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  if ( ido == 2 ) then
    nac = 1
    return
  end if

  nac = 0

  c2(1:idl1,1)    = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)
  c1(2,1:l1,2:ip) = ch(2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) + wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   - wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) + wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   - wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSF2 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,2,l1)
  real ch(ido,l1,2)
  integer i
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)

        ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
        ti2       = cc(i,1,k) - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 - wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 + wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSF3 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,3,l1)
  real ch(ido,l1,3)
  real ci2
  real ci3
  real cr2
  real cr3
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer k
  real taui
  real, parameter :: taur = -0.5E+00
  real ti2
  real tr2
  real wa1(ido)
  real wa2(ido)
!
  taui = - sqrt ( 3.0E+00 ) / 2.0E+00

  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k)   - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSF4 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,4,l1)
  real ch(ido,l1,4)
  real ci2
  real ci3
  real ci4
  real cr2
  real cr3
  real cr4
  integer i
  integer k
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,2,k) - cc(2,4,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,4,k) - cc(1,2,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k)   - cc(i,3,k)
        ti2 = cc(i,1,k)   + cc(i,3,k)
        ti3 = cc(i,2,k)   + cc(i,4,k)
        tr4 = cc(i,2,k)   - cc(i,4,k)
        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,4,k) - cc(i-1,2,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 + wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 - wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 + wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 - wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 + wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 - wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSF5 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,5,l1)
  real ch(ido,l1,5)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer k
  real, parameter :: ti11 = -0.951056516295154E+00
  real, parameter :: ti12 = -0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real tr2
  real tr3
  real tr4
  real tr5
!
!  cos ( 72 ) = +0.3090
!
  real, parameter :: tr11 =  0.309016994374947E+00
!
!  cos ( 36 ) = +0.809016
!
  real, parameter :: tr12 = -0.809016994374947E+00
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 + wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 - wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 + wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 - wa4(i) * dr5

      end do
    end do

  end if

  return
end
function r_cas ( x )
!
!*******************************************************************************
!
!! R_CAS returns the "casine" of a number.
!
!
!  Definition:
!
!    The "casine", used in the discrete Hartley transform, is abbreviated
!    CAS(X), and defined by:
!
!      CAS(X) = cos ( X ) + sin( X )
!             = sqrt ( 2 ) * sin ( X + pi/4 )
!             = sqrt ( 2 ) * cos ( X - pi/4 )
!
!  Modified:
!
!    06 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number whose casine is desired.
!
!    Output, real R_CAS, the casine of X, which will be between
!    plus or minus the square root of 2.
!
  implicit none
!
  real r_cas
  real x
!
  r_cas = cos ( x ) + sin ( x )

  return
end
function r_pi ()
!
!*******************************************************************************
!
!! R_PI returns the value of pi.
!
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R_PI, the value of PI.
!
  implicit none
!
  real r_pi
!
  r_pi = 3.14159265358979323846264338327950288419716939937510E+00

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    06 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  implicit none
!
  real r
  real rhi
  real rlo
  real t
!
!  Pick T, a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R in ( RLO, RHI ).
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP swaps two real values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine radb2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! RADB2 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,2,l1)
  real ch(ido,l1,2)
  integer i
  integer ic
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  ch(1,1:l1,1) = cc(1,1,1:l1) + cc(ido,2,1:l1)
  ch(1,1:l1,2) = cc(1,1,1:l1) - cc(ido,2,1:l1)

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ch(i-1,k,1) = cc(i-1,1,k) + cc(ic-1,2,k)
        tr2         = cc(i-1,1,k) - cc(ic-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   - cc(ic,2,k)
        ti2         = cc(i,1,k)   + cc(ic,2,k)

        ch(i-1,k,2) = wa1(i-2) * tr2 - wa1(i-1) * ti2
        ch(i,k,2)   = wa1(i-2) * ti2 + wa1(i-1) * tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(ido,1:l1,1) =    cc(ido,1,1:l1) + cc(ido,1,1:l1)
  ch(ido,1:l1,2) = -( cc(1,2,1:l1)   + cc(1,2,1:l1) )

  return
end
subroutine radb3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! RADB3 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,3,l1)
  real ch(ido,l1,3)
  real ci2
  real ci3
  real cr2
  real cr3
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer ic
  integer k
  real taui
  real, parameter :: taur = -0.5E+00
  real ti2
  real tr2
  real wa1(ido)
  real wa2(ido)
!
  taui = sqrt ( 3.0E+00 ) / 2.0E+00

  do k = 1, l1

    tr2 = cc(ido,2,k) + cc(ido,2,k)
    cr2 = cc(1,1,k) + taur * tr2
    ch(1,k,1) = cc(1,1,k) + tr2
    ci3 = taui * ( cc(1,3,k) + cc(1,3,k) )

    ch(1,k,2) = cr2 - ci3
    ch(1,k,3) = cr2 + ci3

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      cr2 = cc(i-1,1,k) + taur * tr2
      ch(i-1,k,1) = cc(i-1,1,k) + tr2

      ti2 = cc(i,3,k) - cc(ic,2,k)
      ci2 = cc(i,1,k) + taur * ti2
      ch(i,k,1) = cc(i,1,k) + ti2

      cr3 = taui * ( cc(i-1,3,k) - cc(ic-1,2,k) )
      ci3 = taui * ( cc(i,3,k)   + cc(ic,2,k) )

      dr2 = cr2 - ci3
      dr3 = cr2 + ci3
      di2 = ci2 + cr3
      di3 = ci2 - cr3

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3

    end do
  end do

  return
end
subroutine radb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! RADB4 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,4,l1)
  real ch(ido,l1,4)
  real ci2
  real ci3
  real ci4
  real cr2
  real cr3
  real cr4
  integer i
  integer ic
  integer k
  real, parameter :: sqrt2 = 1.414213562373095E+00
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  do k = 1, l1

    tr1 = cc(1,1,k) - cc(ido,4,k)
    tr2 = cc(1,1,k) + cc(ido,4,k)
    tr3 = cc(ido,2,k) + cc(ido,2,k)
    tr4 = cc(1,3,k) + cc(1,3,k)

    ch(1,k,1) = tr2 + tr3
    ch(1,k,2) = tr1 - tr4
    ch(1,k,3) = tr2 - tr3
    ch(1,k,4) = tr1 + tr4

  end do

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ti1 = cc(i,1,k) + cc(ic,4,k)
        ti2 = cc(i,1,k) - cc(ic,4,k)
        ti3 = cc(i,3,k) - cc(ic,2,k)
        tr4 = cc(i,3,k) + cc(ic,2,k)

        tr1 = cc(i-1,1,k) - cc(ic-1,4,k)
        tr2 = cc(i-1,1,k) + cc(ic-1,4,k)
        ti4 = cc(i-1,3,k) - cc(ic-1,2,k)
        tr3 = cc(i-1,3,k) + cc(ic-1,2,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 - tr4
        cr4 = tr1 + tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-2) * cr2 - wa1(i-1) * ci2
        ch(i,k,2)   = wa1(i-2) * ci2 + wa1(i-1) * cr2
        ch(i-1,k,3) = wa2(i-2) * cr3 - wa2(i-1) * ci3
        ch(i,k,3)   = wa2(i-2) * ci3 + wa2(i-1) * cr3
        ch(i-1,k,4) = wa3(i-2) * cr4 - wa3(i-1) * ci4
        ch(i,k,4)   = wa3(i-2) * ci4 + wa3(i-1) * cr4

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = cc(1,2,k)   + cc(1,4,k)
    ti2 = cc(1,4,k)   - cc(1,2,k)
    tr1 = cc(ido,1,k) - cc(ido,3,k)
    tr2 = cc(ido,1,k) + cc(ido,3,k)

    ch(ido,k,1) = tr2 + tr2
    ch(ido,k,2) = sqrt2 * ( tr1 - ti1 )
    ch(ido,k,3) = ti2 + ti2
    ch(ido,k,4) = -sqrt2 * ( tr1 + ti1 )

  end do

  return
end
subroutine radb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! RADB5 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,5,l1)
  real ch(ido,l1,5)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer ic
  integer k
  real, parameter :: ti11 =  0.951056516295154E+00
  real, parameter :: ti12 =  0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real, parameter :: tr11 =  0.309016994374947E+00
  real, parameter :: tr12 = -0.809016994374947E+00
  real tr2
  real tr3
  real tr4
  real tr5
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  do k = 1, l1

    ti5 = cc(1,3,k) + cc(1,3,k)
    ti4 = cc(1,5,k) + cc(1,5,k)
    tr2 = cc(ido,2,k) + cc(ido,2,k)
    tr3 = cc(ido,4,k) + cc(ido,4,k)

    ch(1,k,1) = cc(1,1,k) + tr2 + tr3
    cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
    cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
    ci5 = ti11 * ti5 + ti12 * ti4
    ci4 = ti12 * ti5 - ti11 * ti4

    ch(1,k,2) = cr2 - ci5
    ch(1,k,3) = cr3 - ci4
    ch(1,k,4) = cr3 + ci4
    ch(1,k,5) = cr2 + ci5

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      ti5 = cc(i,3,k) + cc(ic,2,k)
      ti2 = cc(i,3,k) - cc(ic,2,k)
      ti4 = cc(i,5,k) + cc(ic,4,k)
      ti3 = cc(i,5,k) - cc(ic,4,k)
      tr5 = cc(i-1,3,k) - cc(ic-1,2,k)
      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      tr4 = cc(i-1,5,k) - cc(ic-1,4,k)
      tr3 = cc(i-1,5,k) + cc(ic-1,4,k)

      ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
      ch(i,k,1)   = cc(i,1,k) + ti2 + ti3

      cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
      cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      dr3 = cr3 - ci4
      dr4 = cr3 + ci4
      di3 = ci3 + cr4
      di4 = ci3 - cr4
      dr5 = cr2 + ci5
      dr2 = cr2 - ci5
      di5 = ci2 - cr5
      di2 = ci2 + cr5

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3
      ch(i-1,k,4) = wa3(i-2) * dr4 - wa3(i-1) * di4
      ch(i,k,4)   = wa3(i-2) * di4 + wa3(i-1) * dr4
      ch(i-1,k,5) = wa4(i-2) * dr5 - wa4(i-1) * di5
      ch(i,k,5)   = wa4(i-2) * di5 + wa4(i-1) * dr5

    end do
  end do

  return
end
subroutine radbg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! RADBG is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real ai1
  real ai2
  real ar1
  real ar1h
  real ar2
  real ar2h
  real arg
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  real dc2
  real dcp
  real ds2
  real dsp
  integer i
  integer ic
  integer idij
  integer ik
  integer ipph
  integer is
  integer j
  integer j2
  integer jc
  integer k
  integer l
  integer lc
  integer nbd
  real r_pi
  real wa(*)
!
  arg = 2.0E+00 * r_pi ( ) / real ( ip )
  dcp = cos ( arg )
  dsp = sin ( arg )
  nbd = ( ido - 1 ) / 2
  ipph = ( ip + 1 ) / 2
  ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    ch(1,1:l1,j) =  cc(ido,j2-2,1:l1) + cc(ido,j2-2,1:l1)
    ch(1,1:l1,jc) = cc(1,j2-1,1:l1)   + cc(1,j2-1,1:l1)
  end do

  if ( ido /= 1 ) then

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ic = ido + 2 - i
            ch(i-1,k,j)  = cc(i-1,2*j-1,k) + cc(ic-1,2*j-2,k)
            ch(i-1,k,jc) = cc(i-1,2*j-1,k) - cc(ic-1,2*j-2,k)
            ch(i,k,j)    = cc(i,2*j-1,k)   - cc(ic,2*j-2,k)
            ch(i,k,jc)   = cc(i,2*j-1,k)   + cc(ic,2*j-2,k)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ic = ido + 2 - i
          ch(i-1,1:l1,j)  = cc(i-1,2*j-1,1:l1) + cc(ic-1,2*j-2,1:l1)
          ch(i-1,1:l1,jc) = cc(i-1,2*j-1,1:l1) - cc(ic-1,2*j-2,1:l1)
          ch(i,1:l1,j)    = cc(i,2*j-1,1:l1)   - cc(ic,2*j-2,1:l1)
          ch(i,1:l1,jc)   = cc(i,2*j-1,1:l1)   + cc(ic,2*j-2,1:l1)
        end do
      end do

    end if

  end if

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + ar1 * ch2(ik,2)
      c2(ik,lc) =             ai1 * ch2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2  = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + ar2 * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + ai2 * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    ch(1,1:l1,j)  = c1(1,1:l1,j) - c1(1,1:l1,jc)
    ch(1,1:l1,jc) = c1(1,1:l1,j) + c1(1,1:l1,jc)
  end do

  if ( ido /= 1 ) then

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ch(i-1,k,j)  = c1(i-1,k,j) - c1(i,k,jc)
            ch(i-1,k,jc) = c1(i-1,k,j) + c1(i,k,jc)
            ch(i,k,j)    = c1(i,k,j)   + c1(i-1,k,jc)
            ch(i,k,jc)   = c1(i,k,j)   - c1(i-1,k,jc)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ch(i-1,1:l1,j)  = c1(i-1,1:l1,j) - c1(i,1:l1,jc)
          ch(i-1,1:l1,jc) = c1(i-1,1:l1,j) + c1(i,1:l1,jc)
          ch(i,1:l1,j)    = c1(i,1:l1,j)   + c1(i-1,1:l1,jc)
          ch(i,1:l1,jc)   = c1(i,1:l1,j)   - c1(i-1,1:l1,jc)
        end do
      end do

    end if

  end if

  if ( ido == 1 ) then
    return
  end if

  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)

  if ( nbd <= l1 ) then

    is = -ido

    do j = 2, ip
      is = is + ido
      idij = is
      do i = 3, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    is = -ido
    do j = 2, ip
      is = is + ido
      do k = 1, l1
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine radf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! RADF2 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,l1,2)
  real ch(ido,2,l1)
  integer i
  integer ic
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  ch(1,1,1:l1)   = cc(1,1:l1,1) + cc(1,1:l1,2)
  ch(ido,2,1:l1) = cc(1,1:l1,1) - cc(1,1:l1,2)

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        tr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ti2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)

        ch(i,1,k) = cc(i,k,1) + ti2
        ch(ic,2,k) = ti2 - cc(i,k,1)
        ch(i-1,1,k) = cc(i-1,k,1) + tr2
        ch(ic-1,2,k) = cc(i-1,k,1) - tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(1,2,1:l1) = -cc(ido,1:l1,2)
  ch(ido,1,1:l1) = cc(ido,1:l1,1)

  return
end
subroutine radf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! RADF3 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,l1,3)
  real ch(ido,3,l1)
  real ci2
  real cr2
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer ic
  integer k
  real taui
  real, parameter :: taur = -0.5E+00
  real ti2
  real ti3
  real tr2
  real tr3
  real wa1(ido)
  real wa2(ido)
!
  taui = sqrt ( 3.0E+00 ) / 2.0E+00

  do k = 1, l1
    cr2 = cc(1,k,2) + cc(1,k,3)
    ch(1,1,k) = cc(1,k,1) + cr2
    ch(1,3,k) = taui * ( cc(1,k,3) - cc(1,k,2) )
    ch(ido,2,k) = cc(1,k,1) + taur * cr2
  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)

      cr2 = dr2 + dr3
      ci2 = di2 + di3

      ch(i-1,1,k) = cc(i-1,k,1) + cr2
      ch(i,1,k)   = cc(i,k,1) + ci2

      tr2 = cc(i-1,k,1) + taur * cr2
      ti2 = cc(i,k,1) + taur * ci2
      tr3 = taui * ( di2 - di3 )
      ti3 = taui * ( dr3 - dr2 )

      ch(i-1,3,k) = tr2 + tr3
      ch(ic-1,2,k) = tr2 - tr3
      ch(i,3,k) = ti2 + ti3
      ch(ic,2,k) = ti3 - ti2

    end do
  end do

  return
end
subroutine radf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! RADF4 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,l1,4)
  real ch(ido,4,l1)
  real ci2
  real ci3
  real ci4
  real cr2
  real cr3
  real cr4
  real hsqt2
  integer i
  integer ic
  integer k
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  hsqt2 = sqrt ( 2.0E+00 ) / 2.0E+00

  do k = 1, l1
    tr1 = cc(1,k,2) + cc(1,k,4)
    tr2 = cc(1,k,1) + cc(1,k,3)
    ch(1,1,k) = tr1 + tr2
    ch(ido,4,k) = tr2 - tr1
    ch(ido,2,k) = cc(1,k,1) - cc(1,k,3)
    ch(1,3,k) = cc(1,k,4) - cc(1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        cr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ci2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
        cr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
        ci3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
        cr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
        ci4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)

        tr1 = cr2 + cr4
        tr4 = cr4 - cr2
        ti1 = ci2 + ci4
        ti4 = ci2 - ci4
        ti2 = cc(i,k,1) + ci3
        ti3 = cc(i,k,1) - ci3
        tr2 = cc(i-1,k,1) + cr3
        tr3 = cc(i-1,k,1) - cr3

        ch(i-1,1,k)  = tr1 + tr2
        ch(ic-1,4,k) = tr2 - tr1
        ch(i,1,k)    = ti1 + ti2
        ch(ic,4,k)   = ti1 - ti2
        ch(i-1,3,k)  = ti4 + tr3
        ch(ic-1,2,k) = tr3 - ti4
        ch(i,3,k)    = tr4 + ti3
        ch(ic,2,k)   = tr4 - ti3

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = -hsqt2 * ( cc(ido,k,2) + cc(ido,k,4) )
    tr1 =  hsqt2 * ( cc(ido,k,2) - cc(ido,k,4) )

    ch(ido,1,k) = tr1 + cc(ido,k,1)
    ch(ido,3,k) = cc(ido,k,1) - tr1

    ch(1,2,k) = ti1 - cc(ido,k,3)
    ch(1,4,k) = ti1 + cc(ido,k,3)

  end do

  return
end
subroutine radf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! RADF5 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  implicit none
!
  integer ido
  integer l1
!
  real cc(ido,l1,5)
  real ch(ido,5,l1)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer ic
  integer k
  real, parameter :: ti11 =  0.951056516295154E+00
  real, parameter :: ti12 =  0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real, parameter :: tr11 =  0.309016994374947E+00
  real, parameter :: tr12 = -0.809016994374947E+00
  real tr2
  real tr3
  real tr4
  real tr5
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  do k = 1, l1

    cr2 = cc(1,k,5) + cc(1,k,2)
    ci5 = cc(1,k,5) - cc(1,k,2)
    cr3 = cc(1,k,4) + cc(1,k,3)
    ci4 = cc(1,k,4) - cc(1,k,3)

    ch(1,1,k)   = cc(1,k,1) + cr2 + cr3
    ch(ido,2,k) = cc(1,k,1) + tr11 * cr2 + tr12 * cr3
    ch(1,3,k)   = ti11 * ci5 + ti12 * ci4
    ch(ido,4,k) = cc(1,k,1) + tr12 * cr2 + tr11 * cr3
    ch(1,5,k)   = ti12 * ci5 - ti11 * ci4

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
      dr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
      di4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)
      dr5 = wa4(i-2) * cc(i-1,k,5) + wa4(i-1) * cc(i,k,5)
      di5 = wa4(i-2) * cc(i,k,5)   - wa4(i-1) * cc(i-1,k,5)

      cr2 = dr2 + dr5
      ci5 = dr5 - dr2
      cr5 = di2 - di5
      ci2 = di2 + di5
      cr3 = dr3 + dr4
      ci4 = dr4 - dr3
      cr4 = di3 - di4
      ci3 = di3 + di4

      ch(i-1,1,k) = cc(i-1,k,1) + cr2 + cr3
      ch(i,1,k)   = cc(i,k,1)   + ci2 + ci3

      tr2 = cc(i-1,k,1) + tr11 * cr2 + tr12 * cr3
      ti2 = cc(i,k,1)   + tr11 * ci2 + tr12 * ci3
      tr3 = cc(i-1,k,1) + tr12 * cr2 + tr11 * cr3
      ti3 = cc(i,k,1)   + tr12 * ci2 + tr11 * ci3

      tr5 = ti11 * cr5 + ti12 * cr4
      ti5 = ti11 * ci5 + ti12 * ci4
      tr4 = ti12 * cr5 - ti11 * cr4
      ti4 = ti12 * ci5 - ti11 * ci4

      ch(i-1,3,k)  = tr2 + tr5
      ch(ic-1,2,k) = tr2 - tr5
      ch(i,3,k)    = ti2 + ti5
      ch(ic,2,k)   = ti5 - ti2
      ch(i-1,5,k)  = tr3 + tr4
      ch(ic-1,4,k) = tr3 - tr4
      ch(i,5,k)    = ti3 + ti4
      ch(ic,4,k)   = ti4 - ti3

    end do
  end do

  return
end
subroutine radfg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! RADFG is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer IDO, ?
!
!    Input, integer IP, ?
!
!    Input, integer L1, ?
!
!    Input, integer IDL1, ?
!
!    ?, real CC(IDO,IP,L1), ?
!
!    ?, real C1(IDO,L1,IP), ?
!
!    ?, real C2(IDL1,IP), ?
!
!    ?, real CH(IDO,L1,IP), ?
!
!    ?, real CH2(IDL1,IP), ?
!
!    ?, real WA(*), ?
!
  implicit none
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real ai1
  real ai2
  real ar1
  real ar1h
  real ar2
  real ar2h
  real arg
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  real dc2
  real dcp
  real ds2
  real dsp
  integer i
  integer ic
  integer idij
  integer ik
  integer ipph
  integer is
  integer j
  integer j2
  integer jc
  integer k
  integer l
  integer lc
  integer nbd
  real r_pi
  real wa(*)
!
  arg = 2.0E+00 * r_pi ( ) / real ( ip )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then

    c2(1:idl1,1) = ch2(1:idl1,1)

  else

    ch2(1:idl1,1) = c2(1:idl1,1)
    ch(1,1:l1,2:ip) = c1(1,1:l1,2:ip)

    if ( nbd <= l1 ) then

      is = -ido
      do j = 2, ip
        is = is + ido
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          do k = 1, l1
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    else

      is = -ido

      do j = 2, ip
        is = is + ido
        do k = 1, l1
          idij = is
          do i = 3, ido, 2
            idij = idij + 2
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    end if

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            c1(i-1,k,j)  = ch(i-1,k,j)  + ch(i-1,k,jc)
            c1(i-1,k,jc) = ch(i,k,j)    - ch(i,k,jc)
            c1(i,k,j)    = ch(i,k,j)    + ch(i,k,jc)
            c1(i,k,jc)   = ch(i-1,k,jc) - ch(i-1,k,j)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          c1(i-1,1:l1,j)  = ch(i-1,1:l1,j)  + ch(i-1,1:l1,jc)
          c1(i-1,1:l1,jc) = ch(i,1:l1,j)    - ch(i,1:l1,jc)
          c1(i,1:l1,j)    = ch(i,1:l1,j)    + ch(i,1:l1,jc)
          c1(i,1:l1,jc)   = ch(i-1,1:l1,jc) - ch(i-1,1:l1,j)
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ip + 2 - j
    c1(1,1:l1,j)  = ch(1,1:l1,j)  + ch(1,1:l1,jc)
    c1(1,1:l1,jc) = ch(1,1:l1,jc) - ch(1,1:l1,j)
  end do

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      ch2(ik,l) = c2(ik,1) + ar1 * c2(ik,2)
      ch2(ik,lc) =           ai1 * c2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2 =  dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        ch2(ik,l) =  ch2(ik,l)  + ar2 * c2(ik,j)
        ch2(ik,lc) = ch2(ik,lc) + ai2 * c2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + c2(1:idl1,j)
  end do

  cc(1:ido,1,1:l1) = ch(1:ido,1:l1,1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    cc(ido,j2-2,1:l1) = ch(1,1:l1,j)
    cc(1,j2-1,1:l1)   = ch(1,1:l1,jc)
  end do

  if ( ido == 1 ) then
    return
  end if

  if ( nbd >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do k = 1, l1
        do i = 3, ido, 2
          ic = ido + 2 - i
          cc(i-1,j2-1,k)  = ch(i-1,k,j) + ch(i-1,k,jc)
          cc(ic-1,j2-2,k) = ch(i-1,k,j) - ch(i-1,k,jc)
          cc(i,j2-1,k)    = ch(i,k,j)   + ch(i,k,jc)
          cc(ic,j2-2,k)   = ch(i,k,jc)  - ch(i,k,j)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do i = 3, ido, 2
        ic = ido + 2 - i
        cc(i-1,j2-1,1:l1)  = ch(i-1,1:l1,j) + ch(i-1,1:l1,jc)
        cc(ic-1,j2-2,1:l1) = ch(i-1,1:l1,j) - ch(i-1,1:l1,jc)
        cc(i,j2-1,1:l1)    = ch(i,1:l1,j)   + ch(i,1:l1,jc)
        cc(ic,j2-2,1:l1)   = ch(i,1:l1,jc)  - ch(i,1:l1,j)
      end do
    end do

  end if

  return
end
subroutine random_initialize ( seed )
!
!*******************************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator,
!    and SEED is not changed on output.
!
  implicit none
!
  integer count
  integer count_max
  integer count_rate
  logical, parameter :: debug = .false.
  integer i
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
  real t
!
!  Initialize the random number seed.
!
  call random_seed
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
subroutine rcost ( n, x, wsave )
!
!*******************************************************************************
!
!! RCOST computes the discrete Fourier cosine transform of an even sequence.
!
!
!  Discussion:
!
!    This routine is the unnormalized inverse of itself.  Two successive
!    calls will multiply the input sequence X by 2*(N-1).
!
!    The array WSAVE must be initialized by calling RCOSTI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + (-1) **(I-1) * X_in(N) + sum ( 2 <= K <= N-1 )
!
!        2 * X_in(K) * cos ( ( K - 1 ) * ( I - 1 ) * PI / ( N - 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  The
!    method is more efficient when N-1 is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(3*N+15).
!    The WSAVE array must be initialized by calling RCOSTI.  A different
!    array must be used for each different value of N.
!
  implicit none
!
  integer n
!
  real c1
  integer i
  integer k
  integer kc
  integer ns2
  real t1
  real t2
  real tx2
  real wsave(3*n+15)
  real x(n)
  real x1h
  real x1p3
  real xi
  real xim2
!
  ns2 = n / 2

  if ( n <= 1 ) then
    return
  end if

  if ( n == 2 ) then
    x1h = x(1) + x(2)
    x(2) = x(1) - x(2)
    x(1) = x1h
    return
  end if

  if ( n == 3 ) then
    x1p3 = x(1) + x(3)
    tx2 = x(2) + x(2)
    x(2) = x(1) - x(3)
    x(1) = x1p3 + tx2
    x(3) = x1p3 - tx2
    return
  end if

  c1 = x(1) - x(n)
  x(1) = x(1) + x(n)

  do k = 2, ns2
    kc = n + 1 - k
    t1 = x(k) + x(kc)
    t2 = x(k) - x(kc)
    c1 = c1 + wsave(kc) * t2
    t2 = wsave(k) * t2
    x(k) = t1 - t2
    x(kc) = t1 + t2
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  call rfftf ( n-1, x, wsave(n+1) )

  xim2 = x(2)
  x(2) = c1

  do i = 4, n, 2
    xi = x(i)
    x(i) = x(i-2) - x(i-1)
    x(i-1) = xim2
    xim2 = xi
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n) = xim2
  end if

  return
end
subroutine rcosti ( n, wsave )
!
!*******************************************************************************
!
!! RCOSTI initializes WSAVE, used in RCOST.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  The
!    method is more efficient when N-1 is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the RCOST algorithm.
!
  implicit none
!
  integer n
!
  real dt
  integer k
  real r_pi
  real wsave(3*n+15)
!
  if ( n <= 3 ) then
    return
  end if

  dt = r_pi ( ) / real ( n - 1 )

  do k = 2, ( n / 2 )
    wsave(k)     = 2.0E+00 * sin ( real ( k - 1 ) * dt )
    wsave(n+1-k) = 2.0E+00 * cos ( real ( k - 1 ) * dt )
  end do

  call rffti ( n-1, wsave(n+1) )

  return
end
subroutine rfftb ( n, r, wsave )
!
!*******************************************************************************
!
!! RFFTB computes a real periodic sequence from its Fourier coefficients.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    The transform is unnormalized.  A call to RFFTF followed by a call to
!    RFFTB will multiply the input sequence by N.
!
!    If N is even, the transform is defined by:
!
!      R_out(I) = R_in(1) + (-1)**(I-1) * R_in(N) + sum ( 2 <= K <= N/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    If N is odd, the transform is defined by:
!
!      R_out(I) = R_in(1) + sum ( 2 <= K <= (N+1)/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(2*N+15), a work array.  The WSAVE array must be
!    initialized by calling RFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
  implicit none
!
  integer n
!
  real r(n)
  real wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rfftb1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rfftb1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! RFFTB1 is a lower level routine used by RFFTB.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Input/output, real C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real CH(N).
!
!    Input, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  real c(n)
  real ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nf
  real wa(n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = ido * l1

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radb4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radb4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radb2 ( ido, l1, c, ch, wa(iw) )
      else
        call radb2 ( ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radb3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radb3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radb5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radb5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call radbg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call radbg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( ido == 1 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine rfftf ( n, r, wsave )
!
!*******************************************************************************
!
!! RFFTF computes the Fourier coefficients of a real periodic sequence.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    The transform is unnormalized.  A call to RFFTF followed by a call
!    to RFFTB will multiply the input sequence by N.
!
!    The transform is defined by:
!
!      R_out(1) = sum ( 1 <= I <= N ) R_in(I)
!
!    Letting L = (N+1)/2, then for K = 2,...,L
!
!      R_out(2*K-2) = sum ( 1 <= I <= N )
!
!        R_in(I) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!      R_out(2*K-1) = sum ( 1 <= I <= N )
!
!        -R_in(I) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    And, if N is even, then:
!
!      R_out(N) = sum ( 1 <= I <= N ) (-1)**(I-1) * R_in(I)
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(2*N+15), a work array.  The WSAVE array must be
!    initialized by calling RFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
  implicit none
!
  integer n
!
  real r(n)
  real wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rfftf1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! RFFTF1 is a lower level routine used by RFFTF and RSINT.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Input/output, real C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real CH(N).
!
!    Input, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  real c(n)
  real ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer kh
  integer l1
  integer l2
  integer na
  integer nf
  real wa(n)
!
  nf = ifac(2)
  na = 1
  l2 = n
  iw = n

  do k1 = 1, nf

    kh = nf - k1
    ip = ifac(kh+3)
    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radf4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radf4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radf2 ( ido, l1, c, ch, wa(iw) )
      else
        call radf2 ( ido, l1, ch, c, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radf3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radf3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radf5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radf5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call radfg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
        na = 1
      else
        call radfg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  if ( na /= 1 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine rffti ( n, wsave )
!
!*******************************************************************************
!
!! RFFTI initializes WSAVE, used in RFFTF and RFFTB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Output, real WSAVE(2*N+15), contains data, dependent on the value
!    of N, which is necessary for the RFFTF and RFFTB routines.
!
  implicit none
!
  integer n
!
  real wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rffti1 ( n, wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! RFFTI1 is a lower level routine used by RFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  real arg
  real argh
  real argld
  real fi
  integer i
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer is
  integer j
  integer k1
  integer l1
  integer l2
  integer ld
  integer nf
  real r_pi
  real wa(n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0E+00 * r_pi ( ) / real ( n )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      ld = ld + l1
      i = is
      argld = real ( ld ) * argh
      fi = 0.0E+00

      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0E+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine rsct ( n, x, y )
!
!*******************************************************************************
!
!! RSCT computes a real "slow" cosine transform.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!      Y(1) = Sum ( 1 <= J <= N ) X(J)
!
!      For I from 2 to N-1:
!
!        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
!          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
!
!      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
!    and accuracy.
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real X(N), the data sequence.
!
!    Output, real Y(N), the transformed data.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  real r_pi
  real theta
  real x(n)
  real y(n)
!
  y(1) = sum ( x(1:n) )

  do i = 2, n-1
    y(i) = 0.0E+00
    do j = 1, n
      theta = r_pi ( ) * &
        real ( mod ( ( j - 1 ) * ( i - 1 ), 2 * ( n - 1 ) ) ) / real ( n - 1 )
      y(i) = y(i) + 2.0E+00 * x(j) * cos ( theta )
    end do
  end do

  y(n) = sum ( x(1:n:2) ) - sum ( x(2:n:2) )

  return
end
subroutine rsftb ( n, r, azero, a, b )
!
!*******************************************************************************
!
!! RSFTB computes a "slow" backward Fourier transform of real data.
!
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Output, real R(N), the reconstructed data sequence.
!
!    Input, real AZERO, the constant Fourier coefficient.
!
!    Input, real A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none
!
  integer n
!
  real a(n/2)
  real azero
  real b(n/2)
  integer i
  integer k
  real r(n)
  real r_pi
  real theta
!
  r(1:n) = azero
  do i = 1, n
    do k = 1, n/2
      theta = real ( k * ( i - 1 ) * 2 ) * r_pi() / real ( n )
      r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
    end do
  end do

  return
end
subroutine rsftf ( n, r, azero, a, b )
!
!*******************************************************************************
!
!! RSFTF computes a "slow" forward Fourier transform of real data.
!
!
!  Modified:
!
!    13 March 2001
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real R(N), the data to be transformed.
!
!    Output, real AZERO, = sum ( 1 <= I <= N ) R(I) / N.
!
!    Output, real A(N/2), B(N/2), the Fourier coefficients.
!
  implicit none
!
  integer n
!
  real a(1:n/2)
  real azero
  real b(1:n/2)
  integer i
  integer j
  real r(n)
  real r_pi
  real theta
!
  azero = sum ( r(1:n) ) / real ( n )

  do i = 1, n / 2

    a(i) = 0.0E+00
    b(i) = 0.0E+00

    do j = 1, n
      theta = real ( 2 * i * ( j - 1 ) ) * r_pi() / real ( n )
      a(i) = a(i) + r(j) * cos ( theta )
      b(i) = b(i) + r(j) * sin ( theta )
    end do

    a(i) = a(i) / real ( n )
    b(i) = b(i) / real ( n )

    if ( i /= ( n / 2 ) ) then
      a(i) = 2.0E+00 * a(i)
      b(i) = 2.0E+00 * b(i)
    end if

  end do

  return
end
subroutine rsht ( n, a, b  )
!
!*******************************************************************************
!
!! RSHT computes a "slow" Hartley transform of real data.
!
!
!  Discussion:
!
!    The discrete Hartley transform B of a set of data A is
!
!      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
!
!    Here, the data and coefficients are indexed from 0 to N-1.
!
!    With the above normalization factor of 1/sqrt(N), the Hartley
!    transform is its own inverse.
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines.
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real A(0:N-1), the data to be transformed.
!
!    Output, real B(0:N-1), the transformed data.
!
  implicit none
!
  integer n
!
  real a(0:n-1)
  real b(0:n-1)
  integer i
  integer j
  real r_cas
  real r_pi
  real theta
!
  b(0:n-1) = 0.0E+00

  do i = 0, n-1
    do j = 0, n-1
      theta = 2.0E+00 * r_pi() * real ( mod ( i * j, n ) ) / real ( n )
      b(i) = b(i) + a(j) * r_cas ( theta )
    end do
  end do

  b(0:n-1) = b(0:n-1) / sqrt ( real ( n ) )

  return
end
subroutine rsst ( n, x, y )
!
!*******************************************************************************
!
!! RSST computes a real "slow" sine transform.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For I from 1 to N,
!
!      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
!
!    Applying the routine twice in succession should yield the original data,
!    multiplied by N / 2.  This is a good check for correctness and accuracy.
!
!  Modified:
!
!    18 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real X(N), the data sequence.
!
!    Output, real Y(N), the transformed data.
!
  implicit none
!
  integer n
!
  integer i
  real r_pi
  real theta(n)
  real x(n)
  real y(n)
!
  call rvec_identity ( n, theta )
  theta(1:n) = theta(1:n) * r_pi() / real ( n + 1 )

  y(1:n) = 0.0E+00

  do i = 1, n
    y(1:n) = y(1:n) + 2.0E+00 * x(i) * sin ( real ( i ) * theta(1:n) )
  end do

  return
end
subroutine rsqctb ( n, x, y )
!
!*******************************************************************************
!
!! RSQCTB computes a real "slow" quarter cosine transform backward.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For I from 0 to N-1,
!
!      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
!
!  Reference:
!
!    Briggs and Henson,
!    The Discrete Fourier Transform,
!    SIAM,
!    QA403.5 B75
!
!  Modified:
!
!    21 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real X(0:N-1), the data sequence.
!
!    Output, real Y(0:N-1), the transformed data.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  real r_pi
  real theta
  real x(0:n-1)
  real y(0:n-1)
!
  y(0:n-1) = x(0)

  do i = 0, n-1
    do j = 1, n-1

      theta = 0.5E+00 * r_pi() * real ( j * ( 2 * i + 1 ) ) / real ( n )
      y(i) = y(i) + 2.0E+00 * x(j) * cos ( theta  )

    end do

  end do

  return
end
subroutine rsqctf ( n, x, y )
!
!*******************************************************************************
!
!! RSQCTF computes a real "slow" quarter cosine transform forward.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For I from 0 to N-1,
!
!      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
!
!  Modified:
!
!    07 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Briggs and Henson,
!    The Discrete Fourier Transform,
!    SIAM,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real X(0:N-1), the data sequence.
!
!    Output, real Y(0:N-1), the transformed data.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  real r_pi
  real theta
  real x(0:n-1)
  real y(0:n-1)
!
  y(0:n-1) = 0.0E+00

  do i = 0, n-1
    do j = 0, n-1
      theta = 0.5E+00 * r_pi() * real ( i * ( 2 * j + 1 ) ) / real ( n )
      y(i) = y(i) + x(j) * cos ( theta  )
    end do
  end do

  y(0:n-1) = y(0:n-1) / real ( n )

  return
end
subroutine rsqstb ( n, x, y )
!
!*******************************************************************************
!
!! RSQSTB computes a real "slow" quarter sine transform backward.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For I from 0 to N-1,
!
!      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
!             - X(N) * cos ( pi * I )
!
!  Modified:
!
!    21 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Briggs and Henson,
!    The Discrete Fourier Transform,
!    SIAM,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real X(N), the data sequence.
!
!    Output, real Y(0:N-1), the transformed data.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  real r_pi
  real theta
  real x(1:n)
  real y(0:n-1)
!
  y(0:n-1) = 0.0E+00

  do i = 0, n-1
    do j = 1, n-1

      theta = 0.5E+00 * r_pi() * real ( j * ( 2 * i + 1 ) ) / real ( n )
      y(i) = y(i) - 2.0E+00 * x(j) * sin ( theta  )

    end do

    theta = r_pi() * real ( i )
    y(i) = y(i) - x(n) * cos ( theta )

  end do

  return
end
subroutine rsqstf ( n, x, y )
!
!*******************************************************************************
!
!! RSQSTF computes a real "slow" quarter sine transform forward.
!
!
!  Discussion:
!
!    This routine is provided for illustration and testing.  It is inefficient
!    relative to optimized routines that use fast Fourier techniques.
!
!    For I from 1 to N,
!
!      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
!
!  Modified:
!
!    06 January 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Briggs and Henson,
!    The Discrete Fourier Transform,
!    SIAM,
!    QA403.5 B75
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real X(0:N-1), the data sequence.
!
!    Output, real Y(N), the transformed data.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  real r_pi
  real theta
  real x(0:n-1)
  real y(n)
!
  y(1:n) = 0.0E+00

  do i = 1, n
    do j = 0, n-1
      theta = 0.5E+00 * r_pi() * real ( i * ( 2 * j + 1 ) ) / real ( n )
      y(i) = y(i) + x(j) * sin ( theta  )
    end do
  end do

  y(1:n) = - y(1:n) / real ( n )

  return
end
subroutine rvec_identity ( n, a )
!
!*******************************************************************************
!
!! RVEC_IDENTITY sets a real vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real A(N), the array to be initialized.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
!
  do i = 1, n
    a(i) = real ( i )
  end do

  return
end
subroutine rvec_print_some ( n, a, max_print, title )
!
!*******************************************************************************
!
!! RVEC_PRINT_SOME prints "some" of a real vector.
!
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    10 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer max_print
  character ( len = * ) title
!
  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do

  else if ( max_print >= 3 ) then

    do i = 1, max_print-2
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i6,2x,g14.6)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(i6,2x,g14.6,2x,a)' ) i, a(i), '...more entries...'

  end if

  return
end
subroutine rvec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! RVEC_RANDOM returns a random real vector in a given range.
!
!
!  Modified:
!
!    04 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(N), the vector of randomly chosen values.
!
  implicit none
!
  integer n
!
  real a(n)
  real ahi
  real alo
  integer i
!
  do i = 1, n
    call r_random ( alo, ahi, a(i) )
  end do

  return
end
subroutine rvec_reverse ( n, a )
!
!*******************************************************************************
!
!! RVEC_REVERSE reverses the elements of a real vector.
!
!
!  Example:
!
!    Input:
!
!      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N), the array to be reversed.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
!
  do i = 1, n/2
    call r_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine sinqb ( n, x, wsave )
!
!*******************************************************************************
!
!! SINQB computes the fast sine transform of quarter wave data.
!
!
!  Discussion:
!
!    SINQB computes a sequence from its representation in terms of a sine
!    series with odd wave numbers.
!
!    SINQF is the unnormalized inverse of SINQB since a call of SINQB
!    followed by a call of SINQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!
!        4 * X_in(K) * sin ( ( 2 * K - 1 ) * I * PI / ( 2 * N ) )
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling SINQI.  A different WSAVE array must be used
!    for each different value of N.
!
  implicit none
!
  integer n
!
  real wsave(3*n+15)
  real x(n)
!
  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    x(1) = 4.0E+00 * x(1)
    return
  end if

  x(2:n:2) = -x(2:n:2)

  call cosqb ( n, x, wsave )
!
!  Reverse the X vector.
!
  call rvec_reverse ( n, x )

  return
end
subroutine sinqf ( n, x, wsave )
!
!*******************************************************************************
!
!! SINQF computes the fast sine transform of quarter wave data.
!
!
!  Discussion:
!
!    SINQF computes the coefficients in a sine series representation with
!    only odd wave numbers.
!
!    SINQB is the unnormalized inverse of SINQF since a call of SINQF
!    followed by a call of SINQB will multiply the input sequence X by 4*N.
!
!    The array WSAVE, which is used by SINQF, must be initialized by
!    calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = (-1)**(I-1) * X_in(N) + sum ( 1 <= K <= N-1 )
!        2 * X_in(K) * sin ( ( 2 * I - 1 ) * K * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling SINQI.  A different WSAVE array must be used
!    for each different value of N.
!
  implicit none
!
  integer n
!
  real wsave(3*n+15)
  real x(n)
!
  if ( n <= 1 ) then
    return
  end if
!
!  Reverse the X vector.
!
  call rvec_reverse ( n, x )

  call cosqf ( n, x, wsave )

  x(2:n:2) = -x(2:n:2)

  return
end
subroutine sinqi ( n, wsave )
!
!*******************************************************************************
!
!! SINQI initializes WSAVE, used in SINQF and SINQB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Output, real WSAVE(3*N+15), contains data, dependent on the value
!    of N, which is necessary for the SINQF or SINQB routines.
!
  implicit none
!
  integer n
!
  real wsave(3*n+15)
!
  call cosqi ( n, wsave )

  return
end
subroutine rsint ( n, x, wsave )
!
!*******************************************************************************
!
!! RSINT computes the discrete Fourier sine transform of an odd sequence.
!
!
!  Discussion:
!
!    This routine is the unnormalized inverse of itself since two successive
!    calls will multiply the input sequence X by 2*(N+1).
!
!    The array WSAVE must be initialized by calling RSINTI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!        2 * X_in(K) * sin ( K * I * PI / ( N + 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE((5*N+30)/2), a work array.  The WSAVE array must be
!    initialized by calling RSINTI.  A different WSAVE array must be used
!    for each different value of N.
!
  implicit none
!
  integer n
!
  integer iw1
  integer iw2
  integer iw3
  real wsave((5*n+30)/2)
  real x(n)
!
  iw1 = n / 2 + 1
  iw2 = iw1 + n + 1
  iw3 = iw2 + n + 1

  call rsint1 ( n, x, wsave(1), wsave(iw1), wsave(iw2), wsave(iw3) )

  return
end
subroutine rsint1 ( n, war, was, xh, x, ifac )
!
!*******************************************************************************
!
!! RSINT1 is a lower level routine used by RSINT.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input/output, real WAR(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WAS(N/2).
!
!    Input, real XH(N).
!
!    Input, real X(N+1), ?.
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  implicit none
!
  integer n
!
  integer i
  integer ifac(15)
  integer k
  integer ns2
  real, parameter :: sqrt3 = 1.73205080756888E+00
  real t1
  real t2
  real war(n)
  real was(n/2)
  real x(n+1)
  real xh(n)
  real xhold
!
  xh(1:n) = war(1:n)
  war(1:n) = x(1:n)

  if ( n <= 1 ) then
    xh(1) = 2.0E+00 * xh(1)
    return
  end if

  if ( n == 2 ) then
    xhold = sqrt3 * ( xh(1) + xh(2) )
    xh(2) = sqrt3 * ( xh(1) - xh(2) )
    xh(1) = xhold
    return
  end if

  ns2 = n / 2
  x(1) = 0.0E+00

  do k = 1, n/2
    t1 = xh(k) - xh(n+1-k)
    t2 = was(k) * ( xh(k) + xh(n+1-k) )
    x(k+1) = t1 + t2
!
! ??? N+2-K puts us out of the array...DAMN IT, THIS IS AN ERROR.
!
    x(n+2-k) = t2 - t1
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n/2+2) = 4.0E+00 * xh(n/2+1)
  end if
!
!  This call says there are N+1 things in X.
!
  call rfftf1 ( n+1, x, xh, war, ifac )

  xh(1) = 0.5E+00 * x(1)
  do i = 3, n, 2
    xh(i-1) = -x(i)
    xh(i) = xh(i-2) + x(i-1)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(n) = -x(n+1)
  end if

  x(1:n) = war(1:n)
  war(1:n) = xh(1:n)

  return
end
subroutine rsinti ( n, wsave )
!
!*******************************************************************************
!
!! RSINTI initializes WSAVE, used in RSINT.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations,
!    G. Rodrigue, editor,
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software,
!    W. Cowell, editor,
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is a product of small primes.
!
!    Output, real WSAVE((5*N+30)/2), contains data, dependent on the value
!    of N, which is necessary for the RSINT routine.
!
  implicit none
!
  integer n
!
  real dt
  integer k
  real r_pi
  real wsave((5*n+30)/2)
!
  if ( n <= 1 ) then
    return
  end if

  dt = r_pi ( ) / real ( n + 1 )

  do k = 1, n/2
    wsave(k) = 2.0E+00 * sin ( real ( k ) * dt )
  end do

  call rffti ( n+1, wsave((n/2)+1) )

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function uniform_01_sample ( iseed )
!
!*******************************************************************************
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!
!  Formula:
!
!    ISEED = ISEED * (7**5) mod (2**31 - 1)
!    RANDOM = ISEED * / ( 2**31 - 1 )
!
!  Modified:
!
!    01 March 1999
!
!  Parameters:
!
!    Input/output, integer ISEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real UNIFORM_01_SAMPLE, a random value between 0 and 1.
!
!  Local parameters:
!
!  IA = 7**5
!  IB = 2**15
!  IB16 = 2**16
!  IP = 2**31-1
!
  implicit none
!
  integer, parameter :: ia = 16807
  integer, parameter :: ib15 = 32768
  integer, parameter :: ib16 = 65536
  integer, parameter :: ip = 2147483647
  integer iprhi
  integer iseed
  integer ixhi
  integer k
  integer leftlo
  integer loxa
  real uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real ( iseed ) * 4.656612875E-10

  return
end
