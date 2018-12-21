module orbit_base
  use util
  use orbit_func
  implicit none


  ! ============================================================================
  ! 公開範囲設定
  ! ============================================================================
  private
  ! 構造型
  public :: TOrbit

  ! 定数
  public :: MU, J2, SEC2DAY, AXIS_I, AXIS_J, AXIS_K

  ! 手続き
  public :: get_orbit, calc_period, convert_anomaly, hohmann, lambert
  public :: intersection


  ! ============================================================================
  ! 構造型定義
  ! ============================================================================
  type :: TOrbit
    ! 軌道要素(主)
    real(8) :: epoch        = -1d0 ! 元期 [MJD] : 修正ユリウス日 (Modified Julian Date)
    real(8) :: inclination  = -1d0 ! 軌道傾斜角, i [rad]
    real(8) :: raan         = -1d0 ! 昇交点赤経 (Right Ascension of Ascending Node), Ω [rad]
    real(8) :: eccentricity = -1d0 ! 離心率, e [-]
    real(8) :: arg_perigee  = -1d0 ! 近地点引数, ω (Argument of Perigee) [rad]
    real(8) :: mean_anomaly = -1d0 ! 平均近点角, M [rad]
    real(8) :: mean_motion  = -1d0 ! 平均運動, m [rad/s]

    ! 軌道要素(従)
    real(8) :: semimajor_axis      = -1d0 ! 軌道長半径, a [km]
    real(8) :: semilatus_rectum    = -1d0 ! 半直弦, p [km]
    real(8) :: period              = -1d0 ! 周期, P [s]
    real(8) :: axis_p(3)           = [0d0, 0d0, 0d0] ! 軌道座標ベクトル1 (近点方向)
    real(8) :: axis_q(3)           = [0d0, 0d0, 0d0] ! 軌道座標ベクトル2
    real(8) :: axis_w(3)           = [0d0, 0d0, 0d0] ! 軌道座標ベクトル3 (回転軸方向)
    real(8) :: angular_momentum(3) = [0d0, 0d0, 0d0] ! 角運動量, h

    real(8) :: dt        = -1d0
    real(8) :: eccentric_anomaly = -1d0 ! 離心近点角, E [rad]
    real(8) :: true_anomaly      = -1d0 ! 真近点角, T [rad]

    real(8) :: position(3)       = [0d0, 0d0, 0d0] ! 位置ベクトル, r [km]
    real(8) :: velocity(3)       = [0d0, 0d0, 0d0] ! 速度ベクトル, v [km/s]

    ! 摂動関連
    real(8) :: d_raan        = 0d0 ! 昇交点赤経の時間微分 [rad/s]
    real(8) :: d_arg_perigee = 0d0 ! 近地点引数の時間微分 [rad/s]

    ! その他
    integer :: orbital_shape = 1 ! 軌道形状 (1=>楕円, 0=>放物線, -1=>双曲線)

    contains

    generic :: initialize => initialize_elm, initialize_tle
    generic :: move => move_core, move_default

    procedure :: initialize_elm, initialize_tle
    procedure :: calc_orbit
    procedure :: get_position
    procedure :: get_velocity
    procedure :: move_core, move_default
    procedure :: save
  end type TOrbit

  ! ============================================================================
  ! インターフェイス定義
  ! ============================================================================
  interface get_orbit
    procedure :: create_orbit, position2orbit, peri2orbit
  end interface get_orbit

  interface lambert
    procedure :: lambert_orbit, lambert_core
  end interface lambert


  ! ============================================================================
  ! 定数定義
  ! ============================================================================
  real(8), parameter :: MU        = 3.9860044d5     ! 重力定数µ (== GM) [km^3/s^2]
  real(8), parameter :: J2        = 4.4041964d4     ! 摂動J2項の修正値 (== J2 * re ** 2) [-]
  real(8), parameter :: SEC2DAY   = 1d0 / 86400d0   ! 秒=>日 [day/s]
  real(8), parameter :: AXIS_I(3) = [1d0, 0d0, 0d0] ! 地球座標系ベクトル1 (春分点方向)
  real(8), parameter :: AXIS_J(3) = [0d0, 1d0, 0d0] ! 地球座標系ベクトル2
  real(8), parameter :: AXIS_K(3) = [0d0, 0d0, 1d0] ! 地球座標系ベクトル3 (北極点方向)

contains

  subroutine initialize_elm(this, epc, inc, ran, ecc, ap, ma, mm)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8), intent(in) :: epc, inc, ran, ecc, ap, ma, mm
    real(8) :: matrix(3, 3)

    this%epoch    = epc
    this%inclination    = inc * RADIANS
    this%raan    = ran * RADIANS
    this%eccentricity    = ecc
    this%arg_perigee     = ap * RADIANS
    this%mean_anomaly     = ma * RADIANS
    this%mean_motion     = mm * PI2 * SEC2DAY

    this%semimajor_axis    = (MU / this%mean_motion ** 2) ** (1d0 / 3d0)
    this%semilatus_rectum    = this%semimajor_axis * (1d0 - ecc ** 2)
    this%period    = 86400d0 / mm

    matrix      = euler_angle([this%arg_perigee, this%inclination, this%raan], [3, 1, 3])
    this%axis_p = matrix(1, :)
    this%axis_q = matrix(2, :)
    this%axis_w = matrix(3, :)

    this%d_raan   = -1.5d0 * J2 * this%mean_motion / this%semilatus_rectum ** 2 * cos(this%inclination)
    this%d_arg_perigee    =  1.5d0 * J2 * this%mean_motion / this%semilatus_rectum ** 2 * (2d0 - 2.5d0 * sin(this%inclination) ** 2)
  end subroutine initialize_elm

  subroutine initialize_tle(this, line1, line2)
    implicit none
    class(TOrbit), intent(inout) :: this
    character(69), intent(inout) :: line1, line2
    real(8) :: day, inc, ran, ecc, ap, ma, mm
    integer :: year

    read(line1(19:20), *) year
    read(line1(21:32), *) day

    read(line2( 9:16), *) inc
    read(line2(18:25), *) ran
    line2(25:26) = "0."
    read(line2(25:33), *) ecc
    read(line2(35:42), *) ap
    read(line2(44:51), *) ma
    read(line2(53:63), *) mm

    if (year < 57) then
      year = year + 2000
    else
      year = year + 1900
    end if

    call this%initialize(mjd(year, day), inc, ran, ecc, ap, ma, mm)
  end subroutine initialize_tle


  ! ============================================================================
  ! ============================================================================

  !-----------------------------------------------------------------------------
  ! Constructor
  !-----------------------------------------------------------------------------
  type(TOrbit) function create_orbit(epc, inc, ran, ecc, ap, ma, mm) result (orbit)
    implicit none
    real(8), intent(in) :: epc, inc, ran, ecc, ap, ma, mm
    ! real(8)             :: matrix(3, 3)

    call orbit%initialize(epc, inc, ran, ecc, ap, ma, mm)

    ! orbit%epoch    = epc
    ! orbit%inclination    = inc * RADIANS
    ! orbit%raan    = ran * RADIANS
    ! orbit%eccentricity    = ecc
    ! orbit%arg_perigee     = ap   * RADIANS
    ! orbit%mean_anomaly     = ma   * RADIANS
    ! orbit%mean_motion     = mm   * PI2 * SEC2DAY
    ! orbit%semimajor_axis    = (MU / orbit%mean_motion ** 2) ** (1d0 / 3d0)
    ! orbit%semilatus_rectum    = orbit%semimajor_axis * (1d0 - ecc ** 2)
    ! orbit%period    = 86400d0 / mm

    ! matrix       = euler_angle([orbit%ap, orbit%inc, orbit%ran], [3, 1, 3])
    ! orbit%axis_p = matrix(1, :)
    ! orbit%axis_q = matrix(2, :)
    ! orbit%axis_w = matrix(3, :)

    ! orbit%d_raan   = -1.5d0 * J2 * orbit%mean_motion / orbit%semilatus_rectum ** 2 * cos(orbit%inc)
    ! orbit%d_arg_perigee    =  1.5d0 * J2 * orbit%mean_motion / orbit%semilatus_rectum ** 2 * (2d0 - 2.5d0 * sin(orbit%inc) ** 2)
  end function create_orbit

  type(TOrbit) function position2orbit(pos, vel, date) result (orbit)
    implicit none
    real(8), intent(in) :: pos(3), vel(3), date
    real(8)             :: r, sma, ecc, esin, ecos, ea
    real(8)             :: h(3), axis_p(3), axis_q(3), axis_w(3), axis_n(3)
    integer             :: type

    r      = norm(pos)
    sma    = 1d0 / (2d0 / r - dot_product(vel, vel) / MU)
    h      = cross_product(pos, vel)
    axis_p = unit(-MU * pos / r - cross_product(h, vel))
    axis_w = unit(h)
    axis_q = cross_product(axis_w, axis_p)
    axis_n = unit(cross_product(AXIS_K, axis_w))
    ecc    = sqrt(1d0 - dot_product(h, h) / (MU * sma))
    type   = sgn(1d0 - ecc)
    esin   = dot_product(pos, vel) / sqrt(mu * sma)

    if (ecc > 0) then
      if (type > 0) then
        ecos = 1d0 - r / sma
        ea   = calc_angle(ecos / ecc, esin)
      else
        ea   = asinh(esin / ecc)
      end if
    else
      ea     = 0d0
    end if

    orbit%orbital_shape   = type
    orbit%epoch    = date
    orbit%inclination    = calc_angle(AXIS_K, axis_w)
    orbit%raan    = calc_angle(AXIS_I, axis_n, AXIS_K)
    orbit%eccentricity    = ecc
    orbit%arg_perigee     = calc_angle(axis_n, axis_p, axis_w)
    orbit%mean_anomaly     = type * (ea - esin)
    orbit%mean_motion     = PI2 / calc_period(sma)
    orbit%semimajor_axis    = sma
    orbit%semilatus_rectum    = sma * (1d0 - ecc ** 2)
    orbit%period    = PI2 / orbit%mean_motion
    orbit%axis_p = axis_p
    orbit%axis_q = axis_q
    orbit%axis_w = axis_w
    ! orbit%amom   = h
    ! orbit%dt     = orbit%mean_anomaly / orbit%mean_motion
    ! orbit%eccentric_anomaly     = ea
    orbit%position    = pos
    orbit%velocity   = vel

    orbit%d_raan   = -1.5d0 * J2 * orbit%mean_motion / orbit%semilatus_rectum ** 2 * cos(orbit%inclination)
    orbit%d_arg_perigee    =  1.5d0 * J2 * orbit%mean_motion / orbit%semilatus_rectum ** 2 * (2d0 - 2.5d0 * sin(orbit%inclination) ** 2)
  end function position2orbit

  type(TOrbit) function peri2orbit(pos, sma, axis_w, date) result (orbit)
    implicit none
    real(8), intent(in) :: pos(3), sma, axis_w(3), date
    real(8)  :: r, ecc, axis_p(3), axis_q(3), axis_n(3)
    integer  :: dir

    r      = norm(pos)
    dir    = sgn(sma - r)
    axis_p = dir * unit(pos)
    axis_q = cross_product(axis_w, axis_p)
    axis_n = unit(cross_product(AXIS_K, axis_w))
    ecc    = dir * (1d0 - r / sma)

    orbit%orbital_shape   = sgn(1d0 - ecc)
    orbit%epoch    = date
    orbit%eccentricity    = ecc
    orbit%mean_anomaly     = 0.5d0 * PI * (1 - dir)
    orbit%mean_motion     = PI2 / calc_period(sma)
    orbit%semimajor_axis    = sma
    orbit%semilatus_rectum     = sma * (1d0 - ecc ** 2)
    orbit%period = PI2 / orbit%mean_motion
    orbit%axis_p = axis_p
    orbit%axis_q = axis_q
    orbit%axis_w = axis_w
    orbit%inclination    = calc_angle(AXIS_K, axis_w)
    orbit%raan    = calc_angle(AXIS_I, axis_n, AXIS_K)
    orbit%arg_perigee     = calc_angle(axis_n, axis_p, axis_w)

    orbit%d_raan   = -1.5d0 * J2 * orbit%mean_motion / orbit%semilatus_rectum ** 2 * cos(orbit%inclination)
    orbit%d_arg_perigee    =  1.5d0 * J2 * orbit%mean_motion / orbit%semilatus_rectum ** 2 * (2d0 - 2.5d0 * sin(orbit%inclination) ** 2)
  end function peri2orbit

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  type(TOrbit) function calc_orbit(this, date) result (orbit)
    implicit none
    class(TOrbit), intent(in) :: this
    real(8), intent(in) :: date
    real(8) :: offset, matrix(3, 3)

    orbit%orbital_shape   = this%orbital_shape
    orbit%inclination    = this%inclination
    orbit%eccentricity    = this%eccentricity
    orbit%mean_motion     = this%mean_motion
    orbit%semimajor_axis    = this%semimajor_axis
    orbit%semilatus_rectum    = this%semilatus_rectum
    orbit%period    = this%period
    orbit%d_raan   = this%d_raan
    orbit%d_arg_perigee    = this%d_arg_perigee

    offset       = 86400d0 * (date - this%epoch)
    orbit%epoch    = date
    orbit%raan    = modulo(this%raan + this%d_raan * offset, PI2)
    orbit%arg_perigee     = modulo(this%arg_perigee   + this%d_arg_perigee   * offset, PI2)
    orbit%mean_anomaly     = modulo(this%mean_anomaly   + this%mean_motion    * offset, PI2)
    orbit%eccentric_anomaly     = convert_anomaly(orbit%mean_anomaly, this%eccentricity, from='M', to='E')
    matrix       = euler_angle([orbit%arg_perigee, orbit%inclination, orbit%raan], [3, 1, 3])
    orbit%axis_p = matrix(1, :)
    orbit%axis_q = matrix(2, :)
    orbit%axis_w = matrix(3, :)
    ! orbit%dt    = ma / this%mean_motion
  end function calc_orbit

  function get_position(this) result (pos)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8) :: pos(3), s, c

    if (this%eccentric_anomaly < 0) this%eccentric_anomaly = convert_anomaly(this%mean_anomaly, this%eccentricity, from='M', to='E')

    if (this%orbital_shape > 0) then
      s = sin(this%eccentric_anomaly)
      c = cos(this%eccentric_anomaly)
    else
      s = sinh(this%eccentric_anomaly)
      c = cosh(this%eccentric_anomaly)
    end if
    pos = this%semimajor_axis * (this%orbital_shape * (c - this%eccentricity) * this%axis_p &
                 + sqrt(1d0 - this%eccentricity ** 2) * s * this%axis_q)
  end function get_position

  function get_velocity(this) result (vel)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8) :: vel(3), r, s, c

    if (this%eccentric_anomaly < 0) this%eccentric_anomaly = convert_anomaly(this%mean_anomaly, this%eccentricity, from='M', to='E')

    if (ALL(this%position == 0d0)) then
      r = norm(this%get_position())
    else
      r = norm(this%position)
    end if
    if (this%orbital_shape > 0) then
      s = sin(this%eccentric_anomaly)
      c = cos(this%eccentric_anomaly)
    else
      s = sinh(this%eccentric_anomaly)
      c = cosh(this%eccentric_anomaly)
    end if
    vel = -sqrt(MU) / r * sqrt(this%semimajor_axis) * (s * this%axis_p &
             -  sqrt(1d0 - this%eccentricity ** 2) * c * this%axis_q)
  end function get_velocity

  subroutine move_core(this, arg, mode)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8), intent(in) :: arg
    integer, intent(in) :: mode
    real(8) :: epc, dt, matrix(3, 3)

    select case (mode)
    case (0) ! 時刻 (day)
      epc = arg
      dt = (arg - this%epoch) * 86400d0
    case (1) ! 変動時間 (day)
      epc = this%epoch + arg
      dt = arg * 86400d0
    case (2) !
      epc = this%epoch + this%period * arg / PI2 / 86400d0
      dt = this%period * arg / PI2
    case (-1)
      epc = this%epoch
      dt = 0d0
    case default
      return
    end select

    this%epoch    = epc
    this%raan    = calc_angle(this%raan + this%d_raan * dt)
    this%arg_perigee     = calc_angle(this%arg_perigee  + this%d_arg_perigee  * dt)
    this%mean_anomaly     = calc_angle(this%mean_anomaly  + this%mean_motion   * dt)
    this%eccentric_anomaly     = convert_anomaly(this%mean_anomaly, this%eccentricity, from='M', to='E')
    matrix      = euler_angle([this%arg_perigee, this%inclination, this%raan], [3, 1, 3])
    this%axis_p = matrix(1, :)
    this%axis_q = matrix(2, :)
    this%axis_w = matrix(3, :)
    ! orbit%dt    = ma / this%mean_motion
  end subroutine move_core

  subroutine move_default(this, arg)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8), intent(in) :: arg

    call this%move(arg, 0)
  end subroutine move_default

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  real(8) function convert_anomaly(a, e, from, to) result (z)
    implicit none
    real(8), intent(in) :: a, e
    character, intent(in) :: from, to
    real(8) :: z_next, f, dfdz, c
    integer :: i

    z = a
    select case (from)
    case ('M')
      do i = 1, 100
        if (e < 1) then
          f = z   - e * sin(z) - a
          dfdz = 1d0 - e * cos(z)
        else
          f = e * sinh(z) - z - a
          dfdz = e * cosh(z) - 1d0
        end if
        z_next = z - (f / dfdz)

        if (abs(z_next - z) < 1d-6) exit
        z = z_next
      end do
    case ('T')
      c = cos(z)
      z = calc_angle((c + e) / (1d0 + e * c), PI - z)
    end select

    select case (to)
    case ('M')
      z = z - e * sin(z)
    case ('T')
      c = cos(z)
      z = calc_angle((c - e) / (1d0 - e * c), PI - z)
    end select
  end function convert_anomaly

  real(8) function calc_period(sma) result(period)
    implicit none
    real(8), intent(in) :: sma

    period = PI2 * sqrt(sma ** 3 / MU)
  end function calc_period


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine save(this, date, dt, file, init)
    class(TOrbit), intent(in) :: this
    real(8), intent(in) :: date, dt
    character(*), intent(in) :: file
    logical, intent(in), optional :: init
    type(TOrbit) :: orbit
    ! real(8), allocatable :: data(:, :)
    real(8) :: position(3), ea1, ea2, dea
    integer :: unit, n, i

    orbit = this%calc_orbit(date)

    ea1 = convert_anomaly(orbit%mean_anomaly + orbit%mean_motion * min(dt, 0d0), orbit%eccentricity, from='M', to='E')
    ea2 = convert_anomaly(orbit%mean_anomaly + orbit%mean_motion * max(dt, 0d0), orbit%eccentricity, from='M', to='E')
    dea = ea2 - ea1

    n = int(100 * abs(dt) / orbit%period)
    ! allocate(data(3, 0:n))
    ! print *, n
    ! print *, orbit%period
    ! stop

    if (present(init) .and. init) then
      open(newunit=unit, file=file, status="replace")
    else
      open(newunit=unit, file=file, position="append")
    end if
      do i = 0, n
        orbit%eccentric_anomaly = calc_angle(ea1 + dea * i / dble(n))
        ! orbit%mean_anomaly = convert_anomaly(orbit%eccentric_anomaly, orbit%eccentricity, from='E', to='M')
        ! call orbit%move(0d0, -1)
        position = orbit%get_position()
        write(unit, "(2(f15.5,','),f15.5)") position
      end do
      write(unit, *)
    close(unit)
  end subroutine save


  ! ============================================================================
  ! orbital transfer
  ! ============================================================================

  subroutine hohmann(date, dep, arr, trf, date_d, dt, dv)
    implicit none
    real(8), intent(in) :: date
    type(TOrbit), intent(inout) :: dep, arr
    type(TOrbit), intent(out) :: trf
    real(8), intent(out) :: date_d, dt, dv
    type(TOrbit) :: orbit1, orbit2, orbit_t
    real(8) :: theta, offset, date_a, sma_t, ratio
    real(8) :: axis(3), axis_t(3), pos_d(3), vel_d(3), vel_a(3)
    integer :: i

    i = 1
    ratio = 0.5d0

    orbit1 = dep
    orbit2 = arr

    axis = unit(cross_product(orbit1%axis_w, orbit2%axis_w)) * (-1) ** (i - 1)
    theta = sgn(dot_product(cross_product(orbit1%axis_w, orbit2%axis_w), axis)) &
                  * calc_angle(orbit1%axis_w, orbit2%axis_w)

    offset = calc_angle(convert_anomaly(calc_angle(orbit1%axis_p, axis, orbit1%axis_w), &
                                        orbit1%eccentricity, from='T', to='M') - orbit1%mean_anomaly) / orbit1%mean_motion

    date_d = date + offset * SEC2DAY
    call orbit1%move(date_d, 0)
    pos_d = orbit1%get_position()
    vel_d = orbit1%get_velocity()

    date_a  = date + offset * SEC2DAY
    call orbit2%move(date_a, 0)
    vel_a = orbit2%get_velocity()

    sma_t = norm(pos_d)
    axis_t = Rotate(orbit1%axis_w, axis, ratio * theta)
    orbit_t = get_orbit(pos_d, sma_t, axis_t, date_d)

    trf = orbit_t
    dt = 0.5d0 * trf%period

    dv = norm(vel_d - orbit_t%get_velocity())
    call orbit_t%move(date_a, 0)
    dv = dv + norm(vel_a - orbit_t%get_velocity())
  end subroutine hohmann

  subroutine intersection(dep, arr, axis, theta, offset)
    implicit none
    type(TOrbit), intent(in) :: dep, arr
    real(8), intent(out) :: axis(3), theta, offset
    integer :: i

    i = 1

    axis    = unit(cross_product(dep%axis_w, arr%axis_w)) * (-1) ** (i - 1)
    theta   = sgn(dot_product(cross_product(dep%axis_w, arr%axis_w), axis)) &
                  * calc_angle(dep%axis_w, arr%axis_w)
    offset  = calc_angle(convert_anomaly(calc_angle(dep%axis_p, axis, dep%axis_w), &
                                         dep%eccentricity, from='T', to='M') - dep%mean_anomaly) / dep%mean_motion
  end subroutine intersection

  subroutine lambert_orbit(date, delta_t, dep, arr, delta_v, trf)
    implicit none
    real(8), intent(in) :: date, delta_t
    type(TOrbit), intent(inout) :: dep, arr
    real(8), intent(out) :: delta_v
    type(TOrbit), intent(out), optional :: trf
    real(8) :: pos1(3), pos2(3), vel1(3), vel2(3), axis(3)

    call dep%move(date, 0)
    call arr%move(date + delta_t * SEC2DAY, 0)
    pos1 = dep%get_position()
    pos2 = arr%get_position()
    axis = arr%axis_w
    call lambert(date, delta_t, pos1, pos2, axis, vel1, vel2, trf)

    delta_v = norm(vel1 - dep%get_velocity()) + norm(vel2 - arr%get_velocity())
  end subroutine lambert_orbit

  subroutine lambert_intermediate(date, delta_t, dep, arr, mid, delta_v, trf)
    implicit none
    real(8), intent(in) :: date, delta_t, mid(:, :)
    type(TOrbit), intent(inout) :: dep, arr
    real(8), intent(out) :: delta_v
    type(TOrbit), intent(out), optional :: trf
    real(8) :: pos1(3), pos2(3), vel1(3), vel2(3), axis(3)
    real(8), allocatable :: dv(:)
    integer :: s, i

    ! pos_s = dep%get_position()
    ! pos_d = arr%get_position()
    ! delta_nu = calc_angle(pos1, pos2, axis)

    s = size(mid, dim=2)
    allocate(dv(s + 1))

    do i = 1, s
      call dep%move(date, 0)
      call arr%move(date + delta_t * SEC2DAY, 0)
      pos1 = dep%get_position()
      pos2 = arr%get_position()
      axis = arr%axis_w
      call lambert(date, delta_t, pos1, pos2, axis, vel1, vel2, trf)

      dv(i) = norm(vel1 - dep%get_velocity()) + norm(vel2 - arr%get_velocity())
    end do
  end subroutine lambert_intermediate

  subroutine lambert_core(date, delta_t, pos1, pos2, axis, vel1, vel2, trf)
    implicit none
    real(8), intent(in) :: date, delta_t, pos1(3), pos2(3), axis(3)
    real(8), intent(out) :: vel1(3), vel2(3)
    type(TOrbit), intent(out), optional :: trf
    real(8) :: r1, r2, delta_nu, cos_dnu, sin_dnu, c, s, dtm, dtp
    real(8) :: a, am, t, tm, dt, ddtda, range
    real(8) :: alpha, beta, betam, gamma, delta, f, gi, gd
    real(8) :: cos_nu1, sin_nu1, p, e, ta1, ea1, ma1, axis_n(3)
    integer :: sector, type, round, roundn, i
    logical :: flag

    r1       = norm(pos1)
    r2       = norm(pos2)
    delta_nu = calc_angle(pos1, pos2, axis)
    sin_dnu  = sin(delta_nu)
    cos_dnu  = cos(delta_nu)
    sector   = sgn(PI - delta_nu)

    c        = norm(pos2 - pos1)
    am       = 0.25d0 * (r1 + r2 + c)
    s        = 2d0 * am
    tm       = calc_period(am)
    betam    = 2d0 * asin(sqrt(1d0 - c / s))

    dtm      = tm * (0 + 0.5d0 - sector * (betam - sin(betam)) / PI2)
    dtp      = sqrt(2d0 / MU) * (s ** 1.5d0 - sector * (s - c) ** 1.5d0) / 3d0

    type     = sgn(delta_t - dtp)
    round    = sgn(delta_t - dtm)
    roundn   = (round + 1) / 2 + 0

    a        = 1.1d0 * am
    alpha    = 0d0
    beta     = 0d0
    gamma    = 0d0
    delta    = 0d0
    flag     = .true.
    range    = 0.5d0 * (a - am)

    do i = 1, 1000
      t = calc_period(a)

      if (type > 0) then
        alpha = 2d0 * asin(sqrt(s       / (2d0 * a)))
        beta  = 2d0 * asin(sqrt((s - c) / (2d0 * a)))
        dt    = t * (roundn - (round * (alpha - sin(alpha)) + sector * (beta - sin(beta))) / PI2)
        ! ddtda = 1.5d0 * (roundn * t + dt) / a &
        !       + (round * s ** 2 / sin(alpha) + sector * (s - c) ** 2 / sin(beta)) / sqrt(MU * a ** 3)
        if (round * (delta_t - dt) > 0) then
          if (flag) then
            a     = 2d0 * a
            range = 0.5d0 * a
          else
            a     = a + range
            range = 0.5d0 * range
          end if
        else
          if (flag) flag = .false.
          a     = a - range
          range = 0.5d0 * range
        end if
      else
        gamma = 2d0 * asinh(sqrt(s       / (2d0 * a)))
        delta = 2d0 * asinh(sqrt((s - c) / (2d0 * a)))
        dt    = t / PI2 * (sinh(gamma) - gamma - sector * (sinh(delta) - delta))
        ddtda = 1.5d0 * dt / a &
              - (s ** 2 / sinh(gamma) - sector * (s - c) ** 2 / sinh(delta)) / sqrt(MU * a ** 3)
        a     = max(a + (delta_t - dt) / ddtda, 0.1d0 * a)
      end if

      if (abs(1d0 - dt / delta_t) < 1d-5) exit
    end do

    if (type > 0) then
      p = 4d0 * a * (s - r1) * (s - r2) / c ** 2 * sin(0.5d0 * (alpha - round * sector * beta)) ** 2
      e = sqrt(1d0 - p / a)
    else
      p = 4d0 * a * (s - r1) * (s - r2) / c ** 2 * sinh(0.5d0 * (gamma + sector * delta)) ** 2
      e = sqrt(1d0 + p / a)
    end if

    cos_nu1 = (p - r1) / (e * r1)
    sin_nu1 = (cos_nu1 * cos_dnu - (p - r2) / (e * r2)) / sin_dnu

    ta1     = calc_angle(cos_nu1, sin_nu1)
    ea1     = convert_anomaly(ta1, e, from='T', to='E')
    ma1     = convert_anomaly(ea1, e, from='E', to='M')

    gi      = sqrt(MU * p) / (r1 * r2 * sin_dnu)
    f       = 1d0 - r2 * (1d0 - cos_dnu) / p
    gd      = 1d0 - r1 * (1d0 - cos_dnu) / p

    vel1    = gi * (pos2 - f * pos1)
    vel2    = gi * (gd * pos2 - pos1)

    if (present(trf)) then
      trf%orbital_shape   = type
      trf%epoch    = date
      trf%eccentricity    = e
      trf%semimajor_axis    = a
      trf%semilatus_rectum    = p
      trf%mean_motion     = PI2 / calc_period(a)
      trf%period    = PI2 / trf%mean_motion
      trf%axis_w = sector * unit(cross_product(pos1, pos2))
      trf%axis_p = unit(Rotate(pos1, trf%axis_w, -ta1))
      trf%axis_q = cross_product(trf%axis_w, trf%axis_p)
      axis_n     = unit(cross_product(AXIS_K, trf%axis_w))
      trf%inclination    = calc_angle(AXIS_K, trf%axis_w)
      trf%raan    = calc_angle(AXIS_I, axis_n, AXIS_K)
      trf%arg_perigee     = calc_angle(axis_n, trf%axis_p, trf%axis_w)
      trf%true_anomaly     = ta1
      trf%eccentric_anomaly     = ea1
      trf%mean_anomaly     = ma1
      trf%position    = pos1
      trf%velocity   = vel1
      trf%d_raan   = -1.5d0 * J2 * trf%mean_motion / p ** 2 * cos(trf%inclination)
      trf%d_arg_perigee    =  1.5d0 * J2 * trf%mean_motion / p ** 2 * (2d0 - 2.5d0 * sin(trf%inclination) ** 2)
    end if
  end subroutine lambert_core
end module orbit_base
