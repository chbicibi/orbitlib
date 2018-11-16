module orbit_func
  use util
  implicit none

  private
  public :: sgn, euler_angle, cross_product, unit, calc_angle, rotate, mjd
  public :: PI, PI2, RADIANS, DEGREES

  interface sgn
    procedure :: sgn_i, sgn_d
  end interface sgn

  interface calc_angle
    procedure :: calc_angle_acos, calc_angle1, calc_angle2, calc_angle3
  end interface calc_angle

  real(8), parameter :: PI  = acos(-1d0)
  real(8), parameter :: PI2 = 2d0 * PI
  real(8), parameter :: RADIANS = PI / 180d0
  real(8), parameter :: DEGREES = 180d0 / PI

contains

  integer function sgn_i(a)
    implicit none
    integer, intent(in) :: a

    if (a < 0) then
      sgn_i = -1
    else if (a == 0) then
      sgn_i = 0
    else
      sgn_i = 1
    end if
  end function sgn_i

  integer function sgn_d(a)
    implicit none
    real(8), intent(in) :: a

    if (a < 0) then
      sgn_d = -1
    else
      sgn_d = 1
    end if
  end function sgn_d

  ! ----------------------------------------------------------------------------
  ! vector / matrix
  ! ----------------------------------------------------------------------------
  function cross_product(a, b) result(res)
    implicit none
    real(8), intent(in) :: a(:), b(:)
    real(8), allocatable :: res(:)
    integer :: n, i

    n = size(a)
    res = [(a(mod(i, n) + 1) * b(mod(i + 1, n) + 1) &
          - a(mod(i + 1, n) + 1) * b(mod(i, n) + 1), i = 1, n)]
  end function cross_product

  real(8) function triple_product(a, b, c) result(res)
    implicit none
    real(8), intent(in) :: a(:), b(:), c(:)

    res = dot_product(a, cross_product(b, c))
  end function triple_product

  function unit(a)
    implicit none
    real(8), allocatable :: unit(:)
    real(8), intent(in)  :: a(:)

    unit = a / norm(a)
  end function unit

  function rotate(v, axis, theta)
    implicit none
    real(8)             :: rotate(3)
    real(8), intent(in) :: v(3), axis(3), theta
    real(8)             :: i(3, 3), r(3, 3), m(3, 3)

    i      = reshape([1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0], shape(i))
    r      = reshape([0d0, axis(3), -axis(2), -axis(3), 0d0, axis(1), axis(2), -axis(1), 0d0], shape(r))
    m      = i + sin(theta) * r + (1d0 - cos(theta)) * matmul(r, r)
    rotate = matmul(m, v)
  end function rotate

  function euler_angle(theta, axis)
    implicit none
    real(8)              :: euler_angle(3, 3)
    real(8), intent(in)  :: theta(3)
    integer, intent(in)  :: axis(3)
    real(8)              :: s(3), c(3), r(3, 3)
    integer              :: i

    s = sin(theta)
    c = cos(theta)

    do i = 1, 3
      select case (axis(i))
      case (1)
        r = reshape([1d0, 0d0, 0d0, 0d0, c(i), -s(i), 0d0, s(i), c(i)], shape(r))
      case (2)
        r = reshape([c(i), 0d0, s(i), 0d0, 1d0, 0d0, -s(i), 0d0, c(i)], shape(r))
      case (3)
        r = reshape([c(i), -s(i), 0d0, s(i), c(i), 0d0, 0d0, 0d0, 1d0], shape(r))
      end select

      if (i == 1) then
        euler_angle = r
      else
        euler_angle = matmul(euler_angle, r)
      end if
    end do
  end function euler_angle

  ! ----------------------------------------------------------------------------
  real(8) function calc_angle_acos(cos_, sign_) result(rad)
    implicit none
    real(8), intent(in) :: cos_, sign_

    rad = modulo(sgn(sign_) * acos(cos_), PI2)
  end function calc_angle_acos

  real(8) function calc_angle1(a) result(rad)
    implicit none
    real(8), intent(in) :: a

    rad = modulo(a, PI2)
  end function calc_angle1

  real(8) function calc_angle2(a, b) result(rad)
    implicit none
    real(8), intent(in) :: a(3), b(3)

    rad = acos(min(max(dot_product(a, b) / (norm(a) * norm(b)), -1d0), 1d0))
  end function calc_angle2

  real(8) function calc_angle3(a, b, c) result(rad)
    implicit none
    real(8), intent(in) :: a(3), b(3), c(3)

    rad = calc_angle(sgn(triple_product(c, a, b)) * calc_angle(a, b))
  end function calc_angle3


  ! ============================================================================
  ! ============================================================================

  real(8) function mjd(year, day) result(date)
    implicit none
    integer, intent(in) :: year
    real(8), intent(in) :: day
    integer :: y

    y = year - 1
    date = day - y / 100 + (y / 100) / 4 + int(365.25d0 * y) - 678575
  end function mjd
end module orbit_func
