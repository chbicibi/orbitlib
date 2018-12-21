module orbit_debri
  use util
  use orbit_func
  use orbit_base
  implicit none


  ! ============================================================================
  ! 公開範囲設定
  ! ============================================================================
  private
  public :: TDebri
  public :: NUM_DEBRIS_ALL, DEBRIS
  public :: load_debri, set_startorbit, transfer


  ! ============================================================================
  ! 構造型定義
  ! ============================================================================
  type :: TDebri
    type(TOrbit) :: orbit
    real(8) :: rcs = 0d0
  end type TDebri


  ! ============================================================================
  ! インターフェイス定義
  ! ============================================================================
  interface load_debri
    procedure :: load_debri1, load_debri2
  end interface load_debri


  ! ============================================================================
  ! 定数定義
  ! ============================================================================
  integer :: NUM_DEBRIS_ALL
  type(TDebri), allocatable :: DEBRIS(:)

contains

  ! ----------------------------------------------------------------------------
  ! Initialize
  ! ----------------------------------------------------------------------------
  subroutine load_debri1(num_debris, file_tle, file_rcs, check_size)
    implicit none
    integer, intent(in) :: num_debris
    character(*), intent(in) :: file_tle, file_rcs
    logical, intent(in) :: check_size
    character(69) :: line1, line2
    integer :: unit, num_lines, i

    num_lines = count_lines(file_tle) / 3

    if (num_lines /= count_lines(file_rcs)) then
      print "(a/a)", "Error: the number of lines mismatch", &
                     "in subroutine 'load_debri'"
      call exit(1)
    end if

    if (num_debris < 0) then
      NUM_DEBRIS_ALL = num_lines
    else
      NUM_DEBRIS_ALL = min(num_debris, num_lines)
    end if

    ALLOCATE(DEBRIS(0:NUM_DEBRIS_ALL))

    open(newunit=unit, file = file_tle)
      do i = 1, NUM_DEBRIS_ALL
        read(unit, *)
        read(unit, "(a)") line1
        read(unit, "(a)") line2
        call DEBRIS(i)%orbit%initialize(line1, line2)
      end do
    close(unit)

    open(newunit=unit, file = file_rcs)
      do i = 1, NUM_DEBRIS_ALL
        read(unit, *) DEBRIS(i)%rcs
      end do
    close(unit)

    call set_startorbit
  end subroutine load_debri1

  subroutine load_debri2(num_debris, file_tle, file_rcs)
    integer, intent(in) :: num_debris
    character(*), intent(in) :: file_tle, file_rcs

    call load_debri(num_debris, file_tle, file_rcs, .true.)
  end subroutine load_debri2


  subroutine set_startorbit!(orbit)
    implicit none
    ! type(TOrbit), intent(in) :: orbit

    ! call DEBRIS(0)%orbit = orbit
    call DEBRIS(0)%orbit%initialize(mjd(2018, 0d0), 99d0, 0d0, 1d-2, 0d0, 0d0, 86400d0 / calc_period(6578d0))
  end subroutine set_startorbit


  subroutine transfer(debri_num, date_s, args, del_v, del_t)
    ! ### 引数 ###
    ! debri_num: [出発デブリ番号, 到達デブリ番号]
    ! date_s: 出発日時
    ! args: 遷移パラメータ []
    ! del_v: 合計速度増分 [km/s]
    ! del_t: 合計遷移時間 [s]

    implicit none
    integer, intent(in)  :: debri_num(2)
    real(8), intent(in)  :: date_s, args(6)
    real(8), intent(out) :: del_v, del_t
    type(TOrbit) :: orbit1, orbit2, orbit_d, orbit_t, orbit_a, orbit1_a, orbit2_a
    real(8) :: rho, dnu, t(3), n(3), date, date_d, date_a, theta, axis(3), axis_t(3), offset
    real(8) :: sma_t, tof, pos_d(3), pos_a(3), vel_d(3), vel_t(3), vel_a(3), v1(3), v2(3)
    integer :: i

    rho    = args(1)
    dnu    = args(2)
    t      = [args(3:4), 1d0]
    n      = [args(5:6), 0d0]

    del_v  = 0d0
    del_t  = 0d0

    if (debri_num(1) == debri_num(2)) return

    date     = date_s
    orbit1   = DEBRIS(debri_num(1))%orbit%calc_orbit(date)
    orbit2   = DEBRIS(debri_num(2))%orbit%calc_orbit(date)

    do i = 1, 3
      axis    = unit(cross_product(orbit1%axis_w, orbit2%axis_w)) * (-1) ** (i - 1)
      theta   = sgn(dot_product(cross_product(orbit1%axis_w, orbit2%axis_w), axis)) &
                    * calc_angle(orbit1%axis_w, orbit2%axis_w)

      offset  = calc_angle(convert_anomaly(calc_angle(orbit1%axis_p, axis, orbit1%axis_w), &
                                      orbit1%eccentricity, from='T', to='M') - orbit1%mean_anomaly) / orbit1%mean_motion

      date_d  = date + offset * SEC2DAY
      orbit_d = orbit1%calc_orbit(date_d)
      pos_d   = orbit_d%get_position()
      vel_d   = orbit_d%get_velocity()

      select case (i)
      case (1)
        sma_t   = 0.5d0 * norm(pos_d) * (1d0 + rho) + 1d-10
      case (2)
        offset  = calc_angle(convert_anomaly(calc_angle(orbit2%axis_p, axis, orbit2%axis_w), &
                                        orbit2%eccentricity, from='T', to='M') - orbit1%mean_anomaly) / orbit1%mean_motion
        date_a  = date + offset * SEC2DAY
        orbit_a = orbit2%calc_orbit(date_a)
        pos_a   = orbit_a%get_position()
        sma_t   = 0.5d0 * (norm(pos_d) + norm(pos_a)) + 1d-10
        ! sma_t   = (norm(pos_d) + orbit_a%semimajor_axis - 30d0) / 2d0
      case (3)
        sma_t   = norm(pos_d) + 1d-10
      end select
      axis_t    = rotate(orbit_d%axis_w, axis, theta * t(i))
      orbit_t   = get_orbit(pos_d, sma_t, axis_t, date_d)

      vel_t     = orbit_t%get_velocity()
      del_v     = del_v + norm(vel_t - vel_d)

      date      = date_d + orbit_t%period * n(i) * SEC2DAY
      orbit1    = orbit_t%calc_orbit(date)
      orbit2    = orbit2%calc_orbit(date)
    end do

    tof      = calc_period(0.5d0 * (orbit1%semimajor_axis + orbit2%semimajor_axis)) * dnu
    orbit1_a = orbit1%calc_orbit(date + tof * SEC2DAY)
    orbit2_a = orbit2%calc_orbit(date + tof * SEC2DAY)
    offset   = calc_angle(sgn(orbit2_a%mean_motion - orbit1_a%mean_motion) &
                 * (convert_anomaly(calc_angle(orbit2_a%axis_p, orbit1_a%get_position(), orbit2_a%axis_w), &
                                    orbit2_a%eccentricity, from='T', to='M') - orbit2_a%mean_anomaly))                  &
               / abs(orbit2_a%mean_motion - orbit1_a%mean_motion)
    date_d   = date + offset * SEC2DAY + 1d-1
    orbit_d  = orbit1%calc_orbit(date_d)
    pos_d    = orbit_d%get_position()
    vel_d    = orbit_d%get_velocity()

    orbit_a  = orbit2%calc_orbit(date_d + tof * SEC2DAY)
    pos_a    = orbit_a%get_position()
    vel_a    = orbit_a%get_velocity()

    call lambert(date_d, tof, pos_d, pos_a, orbit_a%axis_w, v1, v2, orbit1)
    del_v    = del_v + norm(v1 - vel_d)
    del_v    = del_v + norm(vel_a - v2)

    del_t    = date_d - date_s + tof * SEC2DAY
  end subroutine transfer
end module orbit_debri
