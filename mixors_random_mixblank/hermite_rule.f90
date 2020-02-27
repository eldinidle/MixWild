!  -*-Fortran -*-
!  Time-stamp: "2007-08-02 14:24:52 CU415 kkim32"
! -------------------------------------------------------------------
!
!   subroutines for HERMITE_RULE.
!
!     (c) Michael L. Berbaum, 2007.  All rights Reserved.
!
! -------------------------------------------------------------------
!
!  Discussion:
!
!    Compute a Gauss-Hermite quadrature rule for approximating
!    Integral ( -oo < x < +oo ) f(x) exp ( - b * ( x - a )^2 ) dx
!
!    The user specifies:
!    * the n_quad (number of points) in the rule;
!    * A, the center point;
!    * B, a scale factor;
!
!  Reference:
!
!     Sylvan Elhay, Jaroslav Kautsky,
!     Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!     Interpolatory Quadrature,
!     ACM Transactions on Mathematical Software,
!     Volume 13, Number 4, December 1987, pages 399-415.
!
!     Roger Martin, James Wilkinson,
!     The Implicit QL Algorithm,
!     Numerische Mathematik,
!     Volume 12, Number 5, December 1968, pages 377-383.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2010
!
!  Author:
!
!    John Burkardt
!
! ----------------------------------------------------------------------
!@  subroutine cdgqf

   subroutine cdgqf ( nt, alpha, beta, t, wts )

!  *********************************************************************
!
!! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
!
!  Discussion:
!
!     This routine computes all the knots and weights of a Gauss quadrature
!     formula with a classical weight function with default values for A and B,
!     and only simple knots.
!
!     There are no moments checks and no printing is done.
!     Use routine EIQFS to evaluate a quadrature computed by CGQFS.
!
!  Parameters:
!
!     Input, integer NT, the number of knots.
!     Input, real ( kind=8 ) ALPHA, the value of Alpha, if needed.
!     Input, real ( kind=8 ) BETA, the value of Beta, if needed.
!     Output, real ( kind=8 ) T(NT), the knots.
!     Output, real ( kind=8 ) WTS(NT), the weights.
!
      implicit none
      integer nt
      real(kind=8), dimension(nt) :: aj,bj,t,wts
      real(kind=8) :: alpha,beta
      real(kind=8) :: zemu

!     Get the Jacobi matrix and zero-th moment.
!
      call class_matrix ( nt, alpha, beta, aj, bj, zemu )
!
!     Compute the knots and weights.
!
      call sgqf ( nt, aj, bj, zemu, t, wts )

      return
   end subroutine cdgqf

! ----------------------------------------------------------------------
!@  subroutine cgqf

   subroutine cgqf ( nt, alpha, beta, a, b, t, wts )

!*****************************************************************************
!
!! CGQF computes knots and weights of a Gauss quadrature formula.
!
!  Discussion:
!
!     The user may specify the interval (A,B).
!     Only simple knots are produced.
!     Use routine EIQFS to evaluate this quadrature formula.
!
!  Parameters:
!
!     Input, integer NT, the number of knots.
!     Input, real ( kind=8 ) ALPHA, the value of Alpha, if needed.
!     Input, real ( kind=8 ) BETA, the value of Beta, if needed.
!     Input, real ( kind=8 ) A, B, the interval endpoints, or
!        other parameters.
!     Output, real ( kind=8 ) T(NT), the knots.
!     Output, real ( kind=8 ) WTS(NT), the weights.
!
      implicit none
      integer,intent(in) :: nt
      real(kind=8),intent(in) :: a,b
      real(kind=8),intent(in) :: alpha,beta
      real(kind=8), dimension(nt),intent(out) :: t,wts
      integer :: i
      integer, allocatable, dimension(:) :: mlt,ndx
!
!     Compute the Gauss quadrature formula for default values of A and B.
      call cdgqf ( nt, alpha, beta, t, wts )
!
!     Prepare to scale the quadrature formula to other weight function with 
!     valid A and B.
!
      allocate ( mlt(1:nt) )

      mlt(1:nt) = 1

      allocate ( ndx(1:nt) )

      do i = 1, nt 
         ndx(i) = i
      end do

      call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, alpha, beta, a, b )

      deallocate ( mlt )
      deallocate ( ndx )
      return

   end subroutine cgqf

! ----------------------------------------------------------------------
!@  subroutine class_matrix 

   subroutine class_matrix ( m, alpha, beta, aj, bj, zemu )

!*****************************************************************************
!
!! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
!
!  Discussion:
!
!     This routine computes the diagonal AJ and sub-diagonal BJ
!     elements of the order M tridiagonal symmetric Jacobi matrix
!     associated with the polynomials orthogonal with respect to
!     the weight function.
!
!     For weight function, M elements are defined in BJ even
!     though only M-1 are needed.  For weight function 8, BJ(M) is
!     set to zero.
!
!     The zero-th moment of the weight function is returned in ZEMU.
!
!  Parameters:
!
!     Input, integer M, the order of the Jacobi matrix.
!     Input, real ( kind=8 ) ALPHA, the value of Alpha, if needed.
!     Input, real ( kind=8 ) BETA, the value of Beta, if needed.
!     Output, real ( kind=8 ) AJ(M), BJ(M), the diagonal and subdiagonal
!        of the Jacobi matrix.
!     Output, real ( kind=8 ) ZEMU, the zero-th moment.
!
      implicit none
      integer  :: m,i
      real(kind=8), dimension(m) :: aj,bj
      real(kind=8)  :: alpha,beta
      real(kind=8)  :: r8_gamma
      real(kind=8)  :: temp,temp2
      real(kind=8)  :: zemu
      REAL(kind=8), PARAMETER :: Pi = 3.141592653589793238462643

      temp = sqrt(epsilon ( temp ))

      temp2 = 0.5D+00

      if ( 500.0D+00 * temp < abs ( ( r8_gamma ( temp2 ) )**2 - pi ) ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
         write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
         stop
      end if

      zemu = r8_gamma ( ( alpha + 1.0D+00 ) / 2.0D+00 )

      aj(1:m) = 0.0D+00

      do i = 1, m
         bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
      end do
      bj(1:m) =  sqrt ( bj(1:m) )

      return

   end subroutine class_matrix 

! ----------------------------------------------------------------------
!@  subroutine imtqlx

   subroutine imtqlx ( n, d, e, z )

!*****************************************************************************
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!     This routine is a slightly modified version of the EISPACK routine to 
!     perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!     The authors thank the authors of EISPACK for permission to use this
!     routine. 
!
!     It has been modified to produce the product Q' * Z, where Z is an input 
!     vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!     The changes consist (essentialy) of applying the orthogonal transformations
!     directly to Z as they are generated.
!
!  Parameters:
!
!     Input, integer N, the order of the matrix.
!     Input/output, real ( kind=8 ) D(N), the diagonal entries of the matrix.
!     Input/output, real ( kind=8 ) E(N), the subdiagonal entries of the 
!        matrix, in entries E(1) through E(N-1).
!     Input/output, real ( kind=8 ) Z(N).  On input, a vector.  On output,
!        the value of Q' * Z, where Q is the matrix that diagonalizes the
!        input symmetric tridiagonal matrix.
!
      implicit none
      integer  :: n,i,ii
      integer  :: j,k,l,m,mml
      integer, parameter :: itn = 30
      real(kind=8)  :: b,c,f,g
      real(kind=8), dimension(n) :: d,e,z
      real(kind=8)  :: p,prec,r,s

      prec = epsilon ( prec )

      if ( n == 1 ) then
         return
      end if

      e(n) = 0.0D+00

      do l = 1, n
         j = 0
         do
            do m = l, n
               if ( m == n ) then
                  exit
               end if
               if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
                  exit
               end if
            end do

            p = d(l)

            if ( m == l ) then
               exit
            end if

            if ( itn <= j ) then
               write ( *, '(a)' ) ' '
               write ( *, '(a)' ) 'IMTQLX - Fatal error!'
               write ( *, '(a)' ) '  Iteration limit exceeded.'
               write ( *, '(a,i8)' ) '  J = ', j
               write ( *, '(a,i8)' ) '  L = ', l
               write ( *, '(a,i8)' ) '  M = ', m
               write ( *, '(a,i8)' ) '  N = ', n
               stop
            end if

            j = j + 1
            g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
            r =  sqrt ( g * g + 1.0D+00 )
            g = d(m) - p + e(l) / ( g + sign ( r, g ) )
            s = 1.0D+00
            c = 1.0D+00
            p = 0.0D+00
            mml = m - l

            do ii = 1, mml
               i = m - ii
               f = s * e(i)
               b = c * e(i)
               if ( abs ( g ) <= abs ( f ) ) then
                  c = g / f
                  r =  sqrt ( c * c + 1.0D+00 )
                  e(i+1) = f * r
                  s = 1.0D+00 / r
                  c = c * s
               else
                  s = f / g
                  r =  sqrt ( s * s + 1.0D+00 )
                  e(i+1) = g * r
                  c = 1.0D+00 / r
                  s = s * c
               end if
               g = d(i+1) - p
               r = ( d(i) - g ) * s + 2.0D+00 * c * b
               p = s * r
               d(i+1) = g + p
               g = c * r - b
               f = z(i+1)
               z(i+1) = s * z(i) + c * f
               z(i) = c * z(i) - s * f
            end do
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0D+00
         end do
      end do
!
!     Sorting.
      do ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         do j = ii, n
            if ( d(j) < p ) then
               k = j
               p = d(j)
            end if
         end do
         if ( k /= i ) then
            d(k) = d(i)
            d(i) = p
            p = z(i)
            z(i) = z(k)
            z(k) = p
         end if
      end do

      return

   end subroutine imtqlx

! ----------------------------------------------------------------------
!@  function r8_gamma

   function r8_gamma ( x )

!*****************************************************************************
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!     This routine calculates the gamma function for a real argument X.
!
!     Computation is based on an algorithm outlined in reference 1.
!     The program uses rational functions that approximate the gamma
!     function to at least 20 significant decimal digits.  Coefficients
!     for the approximation over the interval (1,2) are unpublished.
!     Those for the approximation for 12 <= X are from reference 2.
!
!  Reference:
!
!     William Cody,
!     An Overview of Software Development for Special Functions,
!     in Numerical Analysis Dundee, 1975,
!     edited by GA Watson,
!     Lecture Notes in Mathematics 506,
!     Springer, 1976.
!
!     John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!     Charles Mesztenyi, John Rice, Henry Thatcher,
!     Christoph Witzgall,
!     Computer Approximations,
!     Wiley, 1968,
!     LC: QA297.C64.
!
!  Parameters:
!
!     Input, real ( kind=8 ) X, the argument of the function.
!     Output, real ( kind=8 ) R8_GAMMA, the value of the function.
!
      implicit none
      real(kind=8), dimension ( 7 ) :: c = (/ &
         -1.910444077728D-03, &
          8.4171387781295D-04, &
         -5.952379913043012D-04, &
          7.93650793500350248D-04, &
         -2.777777777777681622553D-03, &
          8.333333333333333331554247D-02, &
          5.7083835261D-03 /)
      real(kind=8), parameter :: eps = 2.22D-16
      REAL(kind=8), PARAMETER :: Pi = 3.141592653589793238462643
      REAL(kind=8), PARAMETER :: sqrtPi = sqrt(pi)
      real(kind=8)     :: fact
      integer  :: i,n
      real(kind=8), dimension ( 8 ) :: p = (/ &
         -1.71618513886549492533811D+00, &
          2.47656508055759199108314D+01, &
         -3.79804256470945635097577D+02, &
          6.29331155312818442661052D+02, &
          8.66966202790413211295064D+02, &
         -3.14512729688483675254357D+04, &
         -3.61444134186911729807069D+04, &
          6.64561438202405440627855D+04 /)
      logical  :: parity
      real(kind=8), dimension ( 8 ) :: q = (/ &
         -3.08402300119738975254353D+01, &
          3.15350626979604161529144D+02, &
         -1.01515636749021914166146D+03, &
         -3.10777167157231109440444D+03, &
          2.25381184209801510330112D+04, &
          4.75584627752788110767815D+03, &
         -1.34659959864969306392456D+05, &
         -1.15132259675553483497211D+05 /)
      real(kind=8)  :: r8_gamma
      real(kind=8)  :: res
      real(kind=8)  :: sum
      real(kind=8)  :: x,xden,xnum
      real(kind=8)  :: y,y1,ysq
      real(kind=8)  :: z      
      real(kind=8), parameter :: xbig = 171.624D+00
      real(kind=8), parameter :: xinf = 1.0D+30
      real(kind=8), parameter :: xminin = 2.23D-308

      parity = .false.
      fact = 1.0D+00
      n = 0
      y = x
!
!     Argument is negative.
      if ( y <= 0.0D+00 ) then
         y = - x
         y1 = aint ( y )
         res = y - y1
         if ( res /= 0.0D+00 ) then
            if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
               parity = .true.
            end if
            fact = - pi / sin ( pi * res )
            y = y + 1.0D+00
         else
            res = xinf
            r8_gamma = res
            return
         end if
      end if
!
!     Argument is positive.
      if ( y < eps ) then
!
!     Argument < EPS.
         if ( xminin <= y ) then
            res = 1.0D+00 / y
         else
            res = xinf
            r8_gamma = res
            return
         end if
      else if ( y < 12.0D+00 ) then
         y1 = y
!
!        0.0 < argument < 1.0.
         if ( y < 1.0D+00 ) then
            z = y
            y = y + 1.0D+00
!        1.0 < argument < 12.0.
!        Reduce argument if necessary.
         else
            n = int ( y ) - 1
            y = y - real ( n, kind=8 )
            z = y - 1.0D+00
         end if
!
!        Evaluate approximation for 1.0 < argument < 2.0.
         xnum = 0.0D+00
         xden = 1.0D+00
         do i = 1, 8
            xnum = ( xnum + p(i) ) * z
            xden = xden * z + q(i)
         end do

         res = xnum / xden + 1.0D+00
!
!        Adjust result for case  0.0 < argument < 1.0.
         if ( y1 < y ) then
            res = res / y1
!        Adjust result for case 2.0 < argument < 12.0.
         else if ( y < y1 ) then
            do i = 1, n
               res = res * y
               y = y + 1.0D+00
            end do
         end if
      else
!
!        Evaluate for 12.0 <= argument.
         if ( y <= xbig ) then
            ysq = y * y
            sum = c(7)
            do i = 1, 6
               sum = sum / ysq + c(i)
            end do
            sum = sum / y - y + sqrtpi
            sum = sum + ( y - 0.5D+00 ) * log ( y )
            res = exp ( sum )
         else
            res = xinf
            r8_gamma = res
            return
         end if
      end if
!
!     Final adjustments and return.
      if ( parity ) then
         res = - res
      end if

      if ( fact /= 1.0D+00 ) then
         res = fact / res
      end if

      r8_gamma = res

      return
   
   end function r8_gamma

! ----------------------------------------------------------------------
!@  subroutine scqf

   subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, &
                     st, alpha, beta, a, b )

!*****************************************************************************
!
!! SCQF scales a quadrature formula to a nonstandard interval.
!
!  Discussion:
!
!     The arrays WTS and SWTS may coincide.
!     The arrays T and ST may coincide.
!
!  Parameters:
!
!     Input, integer NT, the number of knots.
!     Input, real ( kind=8 ) T(NT), the original knots.
!     Input, integer MLT(NT), the multiplicity of the knots.
!     Input, real ( kind=8 ) WTS(NWTS), the weights.
!     Input, integer NWTS, the number of weights.
!     Input, integer NDX(NT), used to index the array WTS.  
!        For more details see the comments in CAWIQ.
!     Output, real ( kind=8 ) SWTS(NWTS), the scaled weights.
!     Output, real ( kind=8 ) ST(NT), the scaled knots.
!     Input, real ( kind=8 ) ALPHA, the value of Alpha, if needed.
!     Input, real ( kind=8 ) BETA, the value of Beta, if needed.
!     Input, real ( kind=8 ) A, B, the interval endpoints.
!
      implicit none
      integer  :: nt,nwts
      integer  :: i,k,l
      integer, dimension(nt) :: mlt,ndx
      real(kind=8)  :: a,b
      real(kind=8)  :: al,be
      real(kind=8)  :: alpha,beta
      real(kind=8)  :: p
      real(kind=8)  :: shft,slp
      real(kind=8), dimension(nt)   :: st,t
      real(kind=8), dimension(nwts) :: swts,wts
      real(kind=8)  :: temp,tmp

      temp = epsilon ( temp )

      if ( b <= 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SCQF - Fatal error!'
         write ( *, '(a)' ) '  B <= 0.'
         stop
      end if

      shft = a
      slp = 1.0D+00 / sqrt ( b )
      al = alpha
      be = 0.0D+00
      p = slp**( al + be + 1.0D+00 )

      do k = 1, nt
         st(k) = shft + slp * t(k)
         l = abs ( ndx(k) )
         if ( l /= 0 ) then
            tmp = p
            do i = l, l + mlt(k) - 1
               swts(i) = wts(i) * tmp
               tmp = tmp * slp
            end do
         end if
      end do
      return

   end subroutine scqf

! ----------------------------------------------------------------------
!@  subroutine sgqf

   subroutine sgqf ( nt, aj, bj, zemu, t, wts )

!*****************************************************************************
!
!! SGQF computes knots and weights of a Gauss Quadrature formula.
!
!  Discussion:
!
!     This routine computes all the knots and weights of a Gauss quadrature
!     formula with simple knots from the Jacobi matrix and the zero-th
!     moment of the weight function, using the Golub-Welsch technique.
!
!
!  Reference:
!
!     Sylvan Elhay, Jaroslav Kautsky,
!     Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!     Interpolatory Quadrature,
!     ACM Transactions on Mathematical Software,
!     Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!     Input, integer NT, the number of knots.
!     Input, real ( kind=8 ) AJ(NT), the diagonal of the Jacobi matrix.
!     Input/output, real ( kind=8 ) BJ(NT), the subdiagonal of the Jacobi 
!        matrix, in entries 1 through NT-1.
!     Input, real ( kind=8 ) ZEMU, the zero-th moment of the weight function.
!     Output, real ( kind=8 ) T(NT), the knots.
!     Output, real ( kind=8 ) WTS(NT), the weights.
!
      implicit none
      integer  :: nt
      real(kind=8), dimension(nt)  :: aj,bj,t,wts
      real(kind=8)  :: zemu
!
!     Exit if the zero-th moment is not positive.
      if ( zemu <= 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'SGQF - Fatal error!'
         write ( *, '(a)' ) '  ZEMU <= 0.'
         stop
      end if
!
!     Set up vectors for IMTQLX.
      t(1:nt) = aj(1:nt)
      wts(1) = sqrt ( zemu )
      wts(2:nt) = 0.0D+00
!
!     Diagonalize the Jacobi matrix.
      call imtqlx ( nt, t, bj, wts )
      wts(1:nt) = wts(1:nt)**2

      return

   end subroutine sgqf
   
!call getQuad(dimRE, n_quad, total_quad, pointsF, weightsf)
!dimRE is an integer containing the number of dimensions.
!n_quad is the number of quadrature points in one dimension
!total_quad is n_quad**dimRE.
!weightsF contains the final weights. It is an integer array of length total_quad.
!pointsF contains the final points, in a 2D integer array total_quad by dimRE.
subroutine getQuad(dimRE, n_quad, total_quad, pointsF, weightsF)
  implicit none
  !Input Arguments
  integer, intent(in) :: dimRE, total_quad, n_quad
  !Output Arguments
  real(kind=8), intent(out), dimension (total_quad) :: weightsF
  real(kind=8), intent(out), dimension (total_quad,dimRE) :: pointsF

  !Local variables
  integer  :: i, j, k, counter, skip, j2
  real(kind=8), allocatable, dimension (:,:) :: points,weights
      REAL(kind=8), PARAMETER :: Pi = 3.141592653589793238462643
      REAL(kind=8), PARAMETER :: zero = 0
      REAL(kind=8), PARAMETER :: one = 1

  weightsF(:) = 1
  
  if(dimRE < 1) return
  
  allocate( weights(n_quad, dimRE) )
  allocate( points(n_quad, dimRE) )

  call cgqf(n_quad, zero, zero, zero, one, points(:,1), weights(:,1))
  points(:,1) = points(:,1)*sqrt(2.0)
  weights(:,1) = weights(:,1)/SQRT(PI)
  
  if(dimRE == 1) then
	pointsF(:,1) = points(:,1)
	weightsF(:) = weights(:,1)
	return
  end if
  
  do i=2, dimRE
	points(:,i) = points(:,1)
	weights(:,i) = weights(:,1)
  end do


  do i=1,dimRE
     counter = 1
	 skip = n_quad**(dimRE-i)
     do j=1,total_quad/skip
        do k=1,skip
           call getmod(j2, j, n_quad) 
           weightsF(counter) = weightsF(counter)*weights(j2,i)
           pointsF(counter,i) = points(j2,i)
           counter = counter + 1
        end do
     end do
  end do
  
  IF (ALLOCATED(weights)) deallocate(weights)
  IF (ALLOCATED(points)) deallocate(points)
  
end subroutine getQuad

subroutine getmod(j2, j, n_quad)
  implicit none
  integer  :: j2, j, n_quad

  j2 = mod(j, n_quad)
  if(j2 .eq. 0) j2 = n_quad
end subroutine getmod


! ----------------------------------------------------------------------