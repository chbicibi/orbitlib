module util
  implicit none

  private
  public :: string, serstr
  ! public :: PI
  public :: random, random_array, shuffle, reverse, sort
  public :: filled, integers, mask_index
  public :: norm, distance, normalize, safe_inv
  public :: str, join, nformat, nformat1
  public :: set_randomseed, elapsed_time, passign
  public :: count_lines

  type :: string
    character(:), allocatable :: s

    contains

    generic :: concat => concat_c, concat_s
    procedure :: puts
    procedure :: concat_c
    procedure :: concat_s
  end type string


  ! ******************************************************************
  ! interface
  ! ******************************************************************

  interface random
    procedure :: irandom, drandom
  end interface random

  interface shuffle
    procedure :: ashuffle, ishuffle, lshuffle
  end interface shuffle

  interface reverse
    procedure :: ireverse, dreverse
  end interface reverse

  interface sort
    procedure :: iqsort, dqsort, iqsort_n, dqsort_n
  end interface sort

  interface filled
    procedure :: ifilled, dfilled, lfilled
  end interface filled

  interface str
    procedure :: str_n, str_dig
  end interface str

  interface join
    procedure :: join_c, join_s
  end interface join

  interface passign
    procedure :: passign0, passign1, passign2
  end interface passign

  interface serstr
    procedure :: serstr_n, serstr_a
  end interface serstr

  ! real(8), parameter :: PI = 4d0 * atan(1d0)

  integer, allocatable :: WS_INTEGERS(:)

  contains


  ! ============================================================================
  ! random number
  ! ============================================================================

  subroutine set_randomseed(fixed_num)
    implicit none
    integer, intent(in), optional :: fixed_num
    integer, allocatable :: seed(:)
    integer :: nrand, clock

    call random_seed(size=nrand)
    allocate(seed(nrand))

    if (present(fixed_num) .and. fixed_num > 0) then
      seed = fixed_num
    else
      call system_clock(count=clock)
      seed = clock
    end if

    call random_seed(put=seed)
  end subroutine set_randomseed

  integer function irandom(n) result(res)
    implicit none
    integer, intent(in) :: n
    real(8) :: r

    call random_number(r)
    res = int(r * n) + 1
  end function irandom

  real(8) function drandom(max) result(res)
    implicit none
    real(8), intent(in), optional :: max

    call random_number(res)
    if (present(max)) res = res * max
  end function drandom

  function random_array(size) result(res)
    implicit none
    integer, intent(in) :: size
    real(8), allocatable :: res(:)
    integer :: i

    res = [(random(), i = 1, size)]
  end function random_array

  function ashuffle(array, n) result(res)
    implicit none
    integer, intent(in) :: array(:)
    integer, intent(in), optional :: n
    integer, allocatable :: res(:), a(:)
    integer :: l, s, pos, i

    s = size(array)
    if (present(n) .and. n < s) then
      l = n
    else
      l = s
    end if

    allocate(res(l))
    allocate(a(s), source=array)

    do i = 1, l
      pos = random(s - i + 1)
      res(i) = a(pos)
      a(pos) = a(s - i + 1)
    end do
  end function ashuffle

  function ishuffle(range, n) result(res)
    implicit none
    integer, intent(in) :: range
    integer, intent(in), optional :: n
    integer, allocatable :: res(:)
    ! integer :: size, array(range), pos, i

    ! if (present(n) .and. n < range) then
    !   size = n
    ! else
    !   size = range
    ! end if

    ! res = filled(size, 0)
    ! array = integers(range)

    ! do i = 1, size
    !   pos = random(range - i + 1)
    !   res(i) = array(pos)
    !   array(pos) = array(range - i + 1)
    ! end do
    res = shuffle(integers(range), n)
  end function ishuffle

  function lshuffle(mask, n) result(res)
    implicit none
    logical, intent(in) :: mask(:)
    integer, intent(in) :: n
    logical, allocatable :: res(:)
    integer, allocatable :: offset(:), table(:)
    integer :: s, c, i

    s = size(mask)
    c = count(mask)
    allocate(offset(s))
    offset(1) = 0
    do i = 1, s - 1
      if (mask(i)) then
        offset(i + 1) = offset(i)
      else
        offset(i + 1) = offset(i) + 1
      end if
    end do
    ! table = integers(c) + pack(offset, mask)
    allocate(table(c), source=integers(c) + pack(offset, mask))
    res = filled(s, .false.)
    res(table(shuffle(c, n))) = .true.
  end function lshuffle

  function ireverse(array) result(res)
    implicit none
    integer, intent(in) :: array(:)
    integer, allocatable :: res(:)

    res = array(size(array):1:-1)
  end function ireverse

  function dreverse(array) result(res)
    implicit none
    real(8), intent(in) :: array(:)
    real(8), allocatable :: res(:)

    res = array(size(array):1:-1)
  end function dreverse


  ! ============================================================================
  ! sort
  ! ============================================================================

  recursive function iqsort(array, ord) result(order)
    implicit none
    integer, intent(in) :: array(:)
    integer, intent(in), optional :: ord(:)
    integer, allocatable :: order(:)
    integer :: len, pivot

    if (present(ord)) then
      len = size(ord)
      order = ord
    else
      len = size(array)
      order = integers(len)
    end if
    if (len <= 1) return
    pivot = array(order(len / 2))
    order = [sort(array, pack(order, array(order) < pivot)), &
             pack(order, array(order) == pivot),             &
             sort(array, pack(order, array(order) > pivot))]
  end function iqsort

  recursive function dqsort(array, ord) result(order)
    implicit none
    real(8), intent(in) :: array(:)
    integer, intent(in), optional :: ord(:)
    integer, allocatable :: order(:)
    real(8) :: pivot
    integer :: len

    if (present(ord)) then
      len = size(ord)
      order = ord
    else
      len = size(array)
      order = integers(len)
    end if
    if (len <= 1) return
    pivot = array(order(len / 2))
    order = [sort(array, pack(order, array(order) < pivot)), &
             pack(order, array(order) == pivot),             &
             sort(array, pack(order, array(order) > pivot))]
  end function dqsort

  recursive function iqsort_n(array, n, ord) result(order)
    implicit none
    integer, intent(in) :: array(:)
    integer, intent(in) :: n
    integer, intent(in), optional :: ord(:)
    integer, allocatable :: order(:), left(:), center(:), right(:)
    integer :: len, pivot, lsize, csize

    if (n <= 0) then
      allocate(order(0))
      return
    end if

    if (present(ord)) then
      len = size(ord)
      order = ord
    else
      len = size(array)
      order = integers(len)
    end if
    if (len <= 1) return
    pivot = array(order(len / 2))

    left = sort(array, n, pack(order, array(order) < pivot))
    lsize = size(left)
    if (lsize == n) then
      order = left
      return
    end if

    center = pack(order, array(order) == pivot)
    csize = size(center)
    if (csize >= n - lsize) then
      order = [left, center(1:n-lsize)]
      return
    end if

    right = sort(array, n - lsize - csize, pack(order, array(order) > pivot))
    order = [left, center, right]
  end function iqsort_n

  recursive function dqsort_n(array, n, ord) result(order)
    implicit none
    real(8), intent(in) :: array(:)
    integer, intent(in) :: n
    integer, intent(in), optional :: ord(:)
    integer, allocatable :: order(:), left(:), center(:), right(:)
    real(8) :: pivot
    integer :: len, lsize, csize

    if (n <= 0) then
      allocate(order(0))
      return
    end if

    if (present(ord)) then
      len = size(ord)
      order = ord
    else
      len = size(array)
      order = integers(len)
    end if
    if (len <= 1) return
    pivot = array(order(len / 2))

    left = sort(array, n, pack(order, array(order) < pivot))
    lsize = size(left)
    if (lsize == n) then
      order = left
      return
    end if

    center = pack(order, array(order) == pivot)
    csize = size(center)
    if (csize >= n - lsize) then
      order = [left, center(1:n-lsize)]
      return
    end if

    right = sort(array, n - lsize - csize, pack(order, array(order) > pivot))
    order = [left, center, right]
  end function dqsort_n


  ! ******************************************************************
  ! auxiliary functions
  ! ******************************************************************

  function ifilled(n, a) result(res)
    implicit none
    integer, intent(in) :: n, a
    integer :: res(n)

    res = a
  end function ifilled

  function dfilled(n, a) result(res)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: a
    real(8) :: res(n)

    res = a
  end function dfilled

  function lfilled(n, a) result(res)
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: a
    logical :: res(n)

    res = a
  end function lfilled

  function integers(n) result(res)
    implicit none
    integer, intent(in) :: n
    integer :: res(n), i

    if (.not. allocated(WS_INTEGERS) .or. size(WS_INTEGERS) < n) then
      WS_INTEGERS = [(i, i = 1, n)]
    end if
    res = WS_INTEGERS(1:n)
  end function integers

  function mask_index(mask) result(index)
    implicit none
    logical, intent(in) :: mask(:)
    integer, allocatable :: index(:)
    integer :: i, j

    allocate(index(count(mask)))
    j = 1
    do i = 1, size(mask)
      if (mask(i)) then
        index(j) = i
        j = j + 1
      end if
    end do
    ! index = pack(integers(count(mask)), mask)
  end function mask_index


  ! ============================================================================
  ! vector operations, math
  ! ============================================================================

  real(8) function norm(a) result(res)
    implicit none
    real(8), intent(in) :: a(:)

    res = sqrt(sum((a ** 2)))
  end function norm

  real(8) function distance(a, b) result(res)
    implicit none
    real(8), intent(in) :: a(:), b(:)

    res = sqrt(sum((a - b) ** 2)) ! == norm(a - b)
  end function distance

  function normalize(a) result(res)
    implicit none
    real(8), intent(in) :: a(:)
    real(8), allocatable :: res(:)

    res = a / sqrt(sum((a ** 2))) ! == a / norm(a)
  end function normalize

  real(8) function safe_inv(a) result(res)
    implicit none
    real(8), intent(in) :: a

    if (a == 0d0) then
      res = sign(huge(0d0), a)
    else
      res = 1d0 / a
    end if
  end function safe_inv


  ! ============================================================================
  ! string operations
  ! ============================================================================

  function str_n(n) result(res)
    implicit none
    integer, intent(in) :: n
    character(:), allocatable :: res
    character(32) :: temp

    write(temp, "(i0)") n
    res = trim(temp)
  end function str_n

  function str_dig(n, digits) result(res)
    implicit none
    integer, intent(in) :: n, digits
    character(:), allocatable :: res
    character(32) :: temp

    if (digits == 0) then
      res = str(n)
    else
      write(temp, "(i0." // str(digits) // ")") n
      res = trim(temp)
    end if
  end function str_dig

  function join_c(array, sep) result(res)
    implicit none
    character(*), intent(in) :: array(:)
    character(*), intent(in), optional :: sep
    character(:), allocatable :: res
    integer :: i

    res = array(1)
    if (present(sep)) then
      do i = 2, size(array)
        res = res // sep // array(i)
      end do
    else
      do i = 2, size(array)
        res = res // array(i)
      end do
    end if
  end function join_c

  function join_s(array, sep) result(res)
    implicit none
    type(string), intent(in) :: array(:)
    character(*), intent(in), optional :: sep
    character(:), allocatable :: res
    integer :: i

    res = array(1)%s
    if (present(sep)) then
      do i = 2, size(array)
        res = res // sep // array(i)%s
      end do
    else
      do i = 2, size(array)
        res = res // array(i)%s
      end do
    end if
  end function join_s

  subroutine puts(this)
    implicit none
    class(string), intent(in) :: this

    print "(a)", this%s
  end subroutine puts

  type(string) function concat_c(this, str) result(res)
    implicit none
    class(string), intent(inout) :: this
    character(*), intent(in) :: str

    res = string(this%s // str)
  end function concat_c

  type(string) function concat_s(this, str) result(res)
    implicit none
    class(string), intent(inout) :: this
    class(string), intent(in) :: str

    res = string(this%s // str%s)
  end function concat_s

  type(string) function serstr_n(prefix, n, suffix, digits) result(res)
    implicit none
    character(*), intent(in) :: prefix
    integer, intent(in) :: n
    character(*), intent(in), optional :: suffix
    integer, intent(in), optional :: digits
    integer :: dig

    if (present(digits)) then
      dig = digits
    else
      dig = 0
    end if

    if (present(suffix)) then
      res = string(prefix // str(n, dig) // suffix)
    else
      res = string(prefix // str(n, dig))
    end if
  end function serstr_n

  function serstr_a(prefix, a, suffix, digits) result(res)
    implicit none
    character(*), intent(in) :: prefix
    integer, intent(in) :: a(:)
    character(*), intent(in), optional :: suffix
    integer, intent(in), optional :: digits
    type(string), allocatable :: res(:)
    integer :: i

    res = [(serstr(prefix, a(i), suffix, digits), i = 1, size(a))]
  end function serstr_a

  function nformat(n, s, sep) result(format)
    implicit none
    integer, intent(in) :: n
    character(*), intent(in) :: s
    character(*), intent(in), optional :: sep
    character(:), allocatable :: format

    if (n > 1) then
      if (present(sep)) then
        format = str(n-1) // "(" // s // ",'" // sep // "')" // s
      else
        format = str(n) // "(" // s // ")"
      end if
    else if (n == 1) then
      format = s
    else
      format = ""
    end if
  end function nformat

  function nformat1(n, s, sep) result(format)
    implicit none
    integer, intent(in) :: n
    character(*), intent(in) :: s
    character(*), intent(in), optional :: sep
    character(:), allocatable :: format

    format = "(" // nformat(n, s, sep) // ")"
  end function nformat1


  ! ============================================================================
  ! timer
  ! ============================================================================

  subroutine elapsed_time(start_time)
    implicit none
    integer, intent(in) :: start_time
    integer :: end_time, t_rate, t_max, diff

    call system_clock(end_time, t_rate, t_max)
    if ( end_time < start_time ) then
      diff = (t_max - start_time) + end_time + 1
    else
      diff = end_time - start_time
    endif
    print "(a/, es10.3e2, a)", "Evaluation took:", dble(diff) / t_rate, " seconds of real time"
  end subroutine elapsed_time


  ! ============================================================================
  ! polymorphic assignment
  ! ============================================================================

  subroutine passign0(to, from)
    implicit none
    class(*), intent(out), allocatable :: to
    class(*), intent(in) :: from

    if (allocated(to)) deallocate(to)
    allocate(to, source=from)
  end subroutine passign0

  subroutine passign1(to, from)
    implicit none
    class(*), intent(out), allocatable :: to(:)
    class(*), intent(in) :: from(:)

    if (allocated(to)) deallocate(to)
    allocate(to(size(from)), source=from)
  end subroutine passign1

  subroutine passign2(to, from)
    implicit none
    class(*), intent(out), allocatable :: to(:, :)
    class(*), intent(in) :: from(:, :)
    integer :: s(2)

    s = shape(from)
    if (allocated(to)) deallocate(to)
    allocate(to(s(1), s(2)), source=from)
  end subroutine passign2


  ! ============================================================================
  ! IO
  ! ============================================================================

  integer function count_lines(file) result(n)
    implicit none
    character(*), intent(in) :: file
    character(1000) :: stream
    integer :: unit, stat

    open(newunit=unit, file=file, action="read")
      n = 0
      do
        read(unit, '(a)', iostat=stat) stream
        if (stat < 0) exit
        if (stream(1:1) == "#") then
          cycle
        else
          n = n + 1
        end if
      end do
    close(unit)
  end function count_lines
end module util
