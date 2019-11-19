#include "model.h"

double von_mises_cdf ( double x, double a, double b );
double von_mises_cdf_inv ( double cdf, double a, double b );

static double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
/******************************************************************************/

static double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}

static double r8_error_f ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ERROR_F evaluates the error function ERF.

  Discussion:

    Since some compilers already supply a routine named ERF which evaluates
    the error function, this routine has been given a distinct, if
    somewhat unnatural, name.

    The function is defined by:

      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T^2 ) dT.

    Properties of the function include:

      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
                               ERF(0) =           0.0;
                               ERF(0.476936...) = 0.5;
      Limit ( X -> +Infinity ) ERF(X) =          +1.0.

      0.5 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2006

  Author:

    Original FORTRAN77 version by William Cody.
    C version by John Burkardt.

  Reference:

    William Cody,
    "Rational Chebyshev Approximations for the Error Function",
    Mathematics of Computation,
    1969, pages 631-638.

  Parameters:

    Input, double X, the argument of the error function.

    Output, double R8_ERROR_F, the value of the error function.
*/
{
  double a[5] = {
    3.16112374387056560,
    1.13864154151050156E+02,
    3.77485237685302021E+02,
    3.20937758913846947E+03,
    1.85777706184603153E-01 };
  double b[4] = {
    2.36012909523441209E+01,
    2.44024637934444173E+02,
    1.28261652607737228E+03,
    2.84423683343917062E+03 };
  double c[9] = {
    5.64188496988670089E-01,
    8.88314979438837594,
    6.61191906371416295E+01,
    2.98635138197400131E+02,
    8.81952221241769090E+02,
    1.71204761263407058E+03,
    2.05107837782607147E+03,
    1.23033935479799725E+03,
    2.15311535474403846E-08 };
  double d[8] = {
    1.57449261107098347E+01,
    1.17693950891312499E+02,
    5.37181101862009858E+02,
    1.62138957456669019E+03,
    3.29079923573345963E+03,
    4.36261909014324716E+03,
    3.43936767414372164E+03,
    1.23033935480374942E+03 };
  double del;
  double erfxd;
  int i;
  double p[6] = {
    3.05326634961232344E-01,
    3.60344899949804439E-01,
    1.25781726111229246E-01,
    1.60837851487422766E-02,
    6.58749161529837803E-04,
    1.63153871373020978E-02 };
  double q[5] = {
    2.56852019228982242,
    1.87295284992346047,
    5.27905102951428412E-01,
    6.05183413124413191E-02,
    2.33520497626869185E-03 };
  double sqrpi = 0.56418958354775628695;
  double thresh = 0.46875;
  double xabs;
  double xbig = 26.543;
  double xden;
  double xnum;
  double xsmall = 1.11E-16;
  double xsq;

  xabs = fabs ( ( x ) );
/*
  Evaluate ERF(X) for |X| <= 0.46875.
*/
  if ( xabs <= thresh )
  {
    if ( xsmall < xabs )
    {
      xsq = xabs * xabs;
    }
    else
    {
      xsq = 0.0;
    }

    xnum = a[4] * xsq;
    xden = xsq;

    for ( i = 0; i < 3; i++ )
    {
      xnum = ( xnum + a[i] ) * xsq;
      xden = ( xden + b[i] ) * xsq;
    }

    erfxd = x * ( xnum + a[3] ) / ( xden + b[3] );
  }
/*
  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
*/
  else if ( xabs <= 4.0 )
  {
    xnum = c[8] * xabs;
    xden = xabs;
    for ( i = 0; i < 7; i++ )
    {
      xnum = ( xnum + c[i] ) * xabs;
      xden = ( xden + d[i] ) * xabs;
    }

    erfxd = ( xnum + c[7] ) / ( xden + d[7] );
    xsq = ( ( double ) ( ( int ) ( xabs * 16.0 ) ) ) / 16.0;
    del = ( xabs - xsq ) * ( xabs + xsq );
    erfxd = exp ( - xsq * xsq ) * exp ( -del ) * erfxd;

    erfxd = ( 0.5 - erfxd ) + 0.5;

    if ( x < 0.0 )
    {
      erfxd = -erfxd;
    }
  }
/*
  Evaluate ERFC(X) for 4.0 < |X|.
*/
  else
  {
    if ( xbig <= xabs )
    {
      if ( 0.0 < x )
      {
        erfxd = 1.0;
      }
      else
      {
        erfxd = - 1.0;
      }
    }
    else
    {
      xsq = 1.0 / ( xabs * xabs );

      xnum = p[5] * xsq;
      xden = xsq;

      for ( i = 0; i < 4; i++ )
      {
        xnum = ( xnum + p[i] ) * xsq;
        xden = ( xden + q[i] ) * xsq;
      }

      erfxd = xsq * ( xnum + p[4] ) / ( xden + q[4] );
      erfxd = ( sqrpi - erfxd ) / xabs;
      xsq = ( ( double ) ( ( int ) ( xabs * 16.0 ) ) ) / 16.0;
      del = ( xabs - xsq ) * ( xabs + xsq );
      erfxd = exp ( - xsq * xsq ) * exp ( - del ) * erfxd;

      erfxd = ( 0.5 - erfxd ) + 0.5;

      if ( x < 0.0 )
      {
        erfxd = -erfxd;
      }
    }
  }

  return erfxd;
}

static double r8_modp ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MODP returns the nonnegative remainder of R8 division.

  Discussion:

    If
      REM = R8_MODP ( X, Y )
      RMULT = ( X - REM ) / Y
    then
      X = Y * RMULT + REM
    where REM is always nonnegative.

    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360.0) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.

    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.

  Example:

        I         J     MOD R8_MODP  R8_MODP Factorization

      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number to be divided.

    Input, double Y, the number that divides X.

    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
*/
{
  double value;

  if ( y == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_MODP - Fatal error!\n" );
    fprintf ( stderr, "  R8_MODP ( X, Y ) called with Y = %g\n", y );
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( value < 0.0 )
  {
    value = value + fabs ( y );
  }

  return value;
}

double von_mises_cdf ( double x, double a, double b )

/******************************************************************************/
/*
  Purpose:

    VON_MISES_CDF evaluates the von Mises CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2006

  Author:

    Original FORTRAN77 version by Geoffrey Hill
    C version by John Burkardt

  Reference:

    Geoffrey Hill,
    ACM TOMS Algorithm 518,
    Incomplete Bessel Function I0: The von Mises Distribution,
    ACM Transactions on Mathematical Software,
    Volume 3, Number 3, September 1977, pages 279-284.

    Kanti Mardia, Peter Jupp,
    Directional Statistics,
    Wiley, 2000, QA276.M335

  Parameters:

    Input, double X, the argument of the CDF.
    A - PI <= X <= A + PI.

    Input, double A, B, the parameters of the PDF.
    -PI <= A <= PI,
    0.0 < B.

    Output, double VON_MISES_CDF, the value of the CDF.
*/
{
  double a1 = 12.0;
  double a2 = 0.8;
  double a3 = 8.0;
  double a4 = 1.0;
  double arg;
  double c;
  double c1 = 56.0;
  double cdf;
  double ck = 10.5;
  double cn;
  double erfx;
  int ip;
  int n;
  double p;
  const double r8_pi = 3.14159265358979323;
  double r;
  double s;
  double sn;
  double u;
  double v;
  double y;
  double z;
/*
  We expect -PI <= X - A <= PI.
*/
  if ( x - a <= - r8_pi )
  {
    cdf = 0.0;
    return cdf;
  }

  if ( r8_pi <= x - a )
  {
    cdf = 1.0;
    return cdf;
  }
/*
  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
*/
  z = b;

  u = r8_modp ( x - a + r8_pi, 2.0 * r8_pi );

  if ( u < 0.0 )
  {
    u = u + 2.0 * r8_pi;
  }

  y = u - r8_pi;
/*
  For small B, sum IP terms by backwards recursion.
*/
  if ( z <= ck )
  {
    v = 0.0;

    if ( 0.0 < z )
    {
      ip = ( int ) ( z * a2 - a3 / ( z + a4 ) + a1 );
      p = ( double ) ( ip );
      s = sin ( y );
      c = cos ( y );
      y = p * y;
      sn = sin ( y );
      cn = cos ( y );
      r = 0.0;
      z = 2.0 / z;

      for ( n = 2; n <= ip; n++ )
      {
        p = p - 1.0;
        y = sn;
        sn = sn * c - cn * s;
        cn = cn * c + y * s;
        r = 1.0 / ( p * z + r );
        v = ( sn / p + v ) * r;
      }
    }
    cdf = ( u * 0.5 + v ) / r8_pi;
  }
/*
  For large B, compute the normal approximation and left tail.
*/
  else
  {
    c = 24.0 * z;
    v = c - c1;
    r = sqrt ( ( 54.0 / ( 347.0 / v + 26.0 - c ) - 6.0 + c ) / 12.0 );
    z = sin ( 0.5 * y ) * r;
    s = 2.0 * z * z;
    v = v - s + 3.0;
    y = ( c - s - s - 16.0 ) / 3.0;
    y = ( ( s + 1.75 ) * s + 83.5 ) / v - y;
    arg = z * ( 1.0 - s / y / y );
    erfx = r8_error_f ( arg );
    cdf = 0.5 * erfx + 0.5;
  }
  cdf = r8_max ( cdf, 0.0 );
  cdf = r8_min ( cdf, 1.0 );

  return cdf;
}
/******************************************************************************/

double von_mises_cdf_inv ( double cdf, double a, double b )

/******************************************************************************/
/*
  Purpose:

    VON_MISES_CDF_INV inverts the von Mises CDF.

  Discussion:

    A simple bisection method is used on the interval [ A - PI, A + PI ].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double CDF, the value of the CDF.

    Input, double A, B, the parameters of the PDF.
    -PI <= A <= PI,
    0.0 < B.

    Output, double VON_MISES_CDF_INV, the corresponding argument of the CDF.
    A - PI <= X <= A + PI.
*/
{
  double cdf1;
  double cdf3;
  int it;
  int it_max = 100;
  const double r8_pi = 3.14159265358979323;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "VON_MISES_CDF_INV - Fatal error!\n" );
    fprintf ( stderr, "  CDF < 0 or 1 < CDF.\n" );
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = a - r8_pi;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = a + r8_pi;
    return x;
  }
  x1 = a - r8_pi;
  cdf1 = 0.0;

  x2 = a + r8_pi;
/*
  Now use bisection.
*/
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = von_mises_cdf ( x3, a, b );

    if ( fabs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "VON_MISES_CDF_INV - Fatal error!\n" );
      fprintf ( stderr, "  Iteration limit exceeded.\n" );
      exit ( 1 );
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
    }
  }

  return x;
}
/******************************************************************************/
