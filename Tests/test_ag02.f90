!
! Verification test for BUG-02, BUG-03, BUG-04 (trace2, diagonal_get_d real return type, sort & quicksort)
! gfortran -fopenmp -I../lib/gnu/ test_ag02.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t
!
program test_bug02_03_04
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use array_mod
  implicit none
  integer, parameter :: dp = 8
  integer :: failures

  failures = 0
  call test_bug02(failures)
  call test_bug03(failures)
  call test_bug04(failures)

  print '(A,I0)', 'test_bug02_03_04 failures = ', failures
  if (failures /= 0) stop 1

contains

  subroutine test_bug02(fail)
    integer, intent(inout) :: fail
    complex(dp) :: Z1(4,4), Z2(4,4), ZH(4,4), Zrect(3,4), ZrectT(4,3), Z00(0,0), zres, zref
    real(dp)    :: D1(4,4), D2(4,4), DS(4,4), Drect(3,4), DrectT(4,3), D00(0,0), dres, dref
    real(dp), parameter :: tol = 1.0e-12_dp
    integer :: i, j

    do j = 1, 4
      do i = 1, 4
        Z1(i,j) = cmplx(real(i - 2*j, dp), real(i + j, dp) * 0.5_dp, kind=dp)
        Z2(i,j) = cmplx(real(3*i + j, dp), real(j - i, dp) * 0.25_dp, kind=dp)
        D1(i,j) = real(i * 2 - j, dp)
        D2(i,j) = real(i + 3 * j, dp)
      end do
    end do
    ZH = Z1
    call hermitian_set(ZH)
    DS = D1
    call symmetric_set(DS)

    ! General complex & real trace2
    zref = sum(Z1 * transpose(Z1))
    if (abs(trace2(Z1) - zref) > tol) fail = fail + 1
    if (abs(trace2(Z1, mtype='') - zref) > tol) fail = fail + 1
    if (abs(trace2(Z1, mtype='G') - zref) > tol) fail = fail + 1   ! 261124: was mtype='H' (BUG-04 of ChangeLog261124)
    ! 261124: one-character 'H' selects the Hermitian shortcut sum(C*conjg(C)), even for a general matrix (BUG-04)
    if (abs(trace2(Z1, mtype='H') - sum(Z1 * conjg(Z1))) > tol) fail = fail + 1

    zref = sum(Z1 * transpose(Z2))
    if (abs(trace2(Z1, Z2) - zref) > tol) fail = fail + 1
    if (abs(trace2(Z1, Z2, mtype='G') - zref) > tol) fail = fail + 1

    dref = sum(D1 * transpose(D1))
    if (abs(trace2(D1) - dref) > tol) fail = fail + 1
    if (abs(trace2(D1, mtype='') - dref) > tol) fail = fail + 1
    if (abs(trace2(D1, mtype='G') - dref) > tol) fail = fail + 1   ! 261124: was mtype='S' (BUG-04 of ChangeLog261124)
    ! 261124: one-character 'S' selects the symmetric shortcut sum(C*C), even for a general matrix (BUG-04)
    if (abs(trace2(D1, mtype='S') - sum(D1 * D1)) > tol) fail = fail + 1

    dref = sum(D1 * transpose(D2))
    if (abs(trace2(D1, D2) - dref) > tol) fail = fail + 1
    if (abs(trace2(D1, D2, mtype='D') - dref) > tol) fail = fail + 1

    ! Hermitian / Symmetric trace2 uppercase and lowercase
    zref = sum(ZH * transpose(ZH))
    if (abs(trace2(ZH, mtype='ZH') - zref) > tol) fail = fail + 1
    if (abs(trace2(ZH, mtype='zh') - zref) > tol) fail = fail + 1
    if (abs(trace2(ZH, mtype='H' ) - zref) > tol) fail = fail + 1   ! 261124: one-character mtype (BUG-04)
    if (abs(trace2(ZH, mtype='h' ) - zref) > tol) fail = fail + 1

    zref = sum(Z1 * transpose(ZH))
    if (abs(trace2(Z1, ZH, mtype='ZH') - zref) > tol) fail = fail + 1
    if (abs(trace2(Z1, ZH, mtype='zh') - zref) > tol) fail = fail + 1
    if (abs(trace2(Z1, ZH, mtype='H' ) - zref) > tol) fail = fail + 1   ! 261124: one-character mtype (BUG-04)

    dref = sum(DS * transpose(DS))
    if (abs(trace2(DS, mtype='DS') - dref) > tol) fail = fail + 1
    if (abs(trace2(DS, mtype='ds') - dref) > tol) fail = fail + 1
    if (abs(trace2(DS, mtype='S' ) - dref) > tol) fail = fail + 1   ! 261124: one-character mtype (BUG-04)
    if (abs(trace2(DS, mtype='s' ) - dref) > tol) fail = fail + 1

    dref = sum(D1 * transpose(DS))
    if (abs(trace2(D1, DS, mtype='DS') - dref) > tol) fail = fail + 1
    if (abs(trace2(D1, DS, mtype='ds') - dref) > tol) fail = fail + 1
    if (abs(trace2(D1, DS, mtype='S' ) - dref) > tol) fail = fail + 1   ! 261124: one-character mtype (BUG-04)

    ! trace2c
    if (abs(trace2c(Z1) - trace2(traceless(Z1))) > tol) fail = fail + 1
    if (abs(trace2c(ZH, mtype='zh') - trace2(traceless(ZH))) > tol) fail = fail + 1
    if (abs(trace2c(D1) - trace2(traceless(D1))) > tol) fail = fail + 1
    if (abs(trace2c(DS, mtype='ds') - trace2(traceless(DS))) > tol) fail = fail + 1

    ! Non-square single matrix -> NaN (no OOB access)
    Zrect  = cmplx(1.0_dp, 2.0_dp, kind=dp)
    ZrectT = cmplx(1.0_dp, 2.0_dp, kind=dp)
    Drect  = 1.0_dp
    DrectT = 1.0_dp

    zres = trace2(Zrect)
    if (.not. isNaN(zres)) fail = fail + 1
    zres = trace2(Zrect, mtype='ZH')
    if (.not. isNaN(zres)) fail = fail + 1
    dres = trace2(Drect)
    if (.not. isNaN(dres)) fail = fail + 1
    dres = trace2(Drect, mtype='DS')
    if (.not. isNaN(dres)) fail = fail + 1

    ! Non-square / mismatched two matrices -> NaN (no OOB access)
    if (.not. isNaN(trace2(Zrect, ZrectT))) fail = fail + 1
    if (.not. isNaN(trace2(Zrect, Zrect))) fail = fail + 1
    if (.not. isNaN(trace2(Zrect, Zrect, mtype='ZH'))) fail = fail + 1
    if (.not. isNaN(trace2(Z1, Z1(1:3,1:3)))) fail = fail + 1
    if (.not. isNaN(trace2(Z1, Z1(1:3,1:3), mtype='ZH'))) fail = fail + 1

    if (.not. isNaN(trace2(Drect, DrectT))) fail = fail + 1
    if (.not. isNaN(trace2(Drect, Drect))) fail = fail + 1
    if (.not. isNaN(trace2(Drect, Drect, mtype='DS'))) fail = fail + 1
    if (.not. isNaN(trace2(D1, D1(1:3,1:3)))) fail = fail + 1
    if (.not. isNaN(trace2(D1, D1(1:3,1:3), mtype='DS'))) fail = fail + 1

    ! 0x0 matrix
    if (trace2(Z00) /= (0.0_dp, 0.0_dp)) fail = fail + 1
    if (trace2(D00) /= 0.0_dp) fail = fail + 1
    if (trace2(Z00, Z00) /= (0.0_dp, 0.0_dp)) fail = fail + 1
    if (trace2(D00, D00) /= 0.0_dp) fail = fail + 1
  end subroutine test_bug02

  subroutine test_bug03(fail)
    integer, intent(inout) :: fail
    real(dp) :: D(4,4), s_expected(4)
    real(dp), allocatable :: diag_v(:), diag_m(:,:), sorted_v(:)
    integer :: i

    D = 0.0_dp
    D(1,1) = -3.0_dp
    D(2,2) =  1.5_dp
    D(3,3) = -0.5_dp
    D(4,4) =  2.0_dp
    D(1,2) =  9.0_dp

    if (storage_size(diagonal(D)) /= 64) then
      print *, 'FAIL BUG-03: storage_size(diagonal(D)) =', storage_size(diagonal(D))
      fail = fail + 1
    end if

    diag_v = diagonal(D)
    if (maxval(abs(diag_v - [-3.0_dp, 1.5_dp, -0.5_dp, 2.0_dp])) > 1.0e-14_dp) fail = fail + 1

    if (storage_size(diagonal(diagonal(D))) /= 64) fail = fail + 1
    diag_m = diagonal(diagonal(D))
    if (maxval(abs(diag_m - diagonalMatrix(D))) > 1.0e-14_dp) fail = fail + 1

    ! Real-only AF and MR keys work on diagonal(D) directly because it is real(8)
    sorted_v = sort(diagonal(D), by='AF')
    s_expected = [-0.5_dp, 1.5_dp, 2.0_dp, -3.0_dp]
    if (maxval(abs(sorted_v - s_expected)) > 1.0e-14_dp) fail = fail + 1

    sorted_v = sort(diagonal(D), by='MR')
    s_expected = [-3.0_dp, 2.0_dp, 1.5_dp, -0.5_dp]
    if (maxval(abs(sorted_v - s_expected)) > 1.0e-14_dp) fail = fail + 1

    if (abs(trace(D) - sum([(D(i,i), i=1,4)])) > 1.0e-14_dp) fail = fail + 1
  end subroutine test_bug03

  subroutine test_bug04(fail)
    integer, intent(inout) :: fail
    real(dp)    :: r_in(10), r0(0), r1(1), r_nan(8), r_allnan(4), r_eq(20)
    complex(dp) :: c_in(10), c0(0), c1(1), c_nan(7), c_allnan(4), c_eq(20)
    real(dp),    allocatable :: r_out(:), r_big(:)
    complex(dp), allocatable :: c_out(:), c_big(:)
    integer,     allocatable :: p(:), p0(:), p1(:), p_nan(:), p_big(:)
    character(2) :: r_keys(15) = ['RF','RR','VF','VR','AF','AR','MF','MR','rf','mr','mR','Mr','mf','mF','Mf']
    character(2) :: c_keys(10) = ['MF','MR','RF','RR','IF','IR','mf','Mr','rF','iR']
    integer :: i, k
    real(dp) :: rnan_val

    rnan_val = NaN(0.0_dp)

    r_in = [3.0_dp, -5.0_dp, 2.0_dp, -2.0_dp, 0.0_dp, 3.0_dp, -1.0_dp, 5.0_dp, -0.5_dp, 2.0_dp]
    c_in = [cmplx(1.0_dp, 2.0_dp, dp), cmplx(-3.0_dp, 1.0_dp, dp), cmplx(1.0_dp, -2.0_dp, dp), &
            cmplx(0.0_dp, 0.0_dp, dp), cmplx(-3.0_dp, -1.0_dp, dp), cmplx(2.0_dp, 0.5_dp, dp), &
            cmplx(1.0_dp, 2.0_dp, dp), cmplx(-1.0_dp, 3.0_dp, dp), cmplx(4.0_dp, -0.5_dp, dp), &
            cmplx(-0.5_dp, -3.0_dp, dp)]

    ! Default real sort ('RF') and all keys with and without P
    allocate(p(10))
    r_out = sort(r_in)
    call check_real_sort(r_in, r_out, 'RF', fail)
    r_out = sort(r_in, P=p)
    call check_real_sort(r_in, r_out, 'RF', fail, p=p)
    do k = 1, size(r_keys)
      r_out = sort(r_in, by=r_keys(k))
      call check_real_sort(r_in, r_out, r_keys(k), fail)
      r_out = sort(r_in, P=p, by=r_keys(k))
      call check_real_sort(r_in, r_out, r_keys(k), fail, p=p)
    end do

    ! Default complex sort ('MF') and all keys with and without P
    c_out = sort(c_in)
    call check_cmplx_sort(c_in, c_out, 'MF', fail)
    c_out = sort(c_in, P=p)
    call check_cmplx_sort(c_in, c_out, 'MF', fail, p=p)
    do k = 1, size(c_keys)
      c_out = sort(c_in, by=c_keys(k))
      call check_cmplx_sort(c_in, c_out, c_keys(k), fail)
      c_out = sort(c_in, P=p, by=c_keys(k))
      call check_cmplx_sort(c_in, c_out, c_keys(k), fail, p=p)
    end do
    deallocate(p)

    ! Empty & 1-element arrays
    allocate(p0(0), p1(1))
    r_out = sort(r0, P=p0, by='MR')
    if (size(r_out) /= 0 .or. size(p0) /= 0) fail = fail + 1
    c_out = sort(c0, P=p0, by='IF')
    if (size(c_out) /= 0 .or. size(p0) /= 0) fail = fail + 1

    r1 = [-4.25_dp]
    r_out = sort(r1, P=p1, by='MR')
    if (size(r_out) /= 1 .or. r_out(1) /= -4.25_dp .or. p1(1) /= 1) fail = fail + 1
    c1 = [cmplx(2.0_dp, -3.0_dp, dp)]
    c_out = sort(c1, P=p1, by='IR')
    if (size(c_out) /= 1 .or. c_out(1) /= c1(1) .or. p1(1) /= 1) fail = fail + 1
    deallocate(p0, p1)

    ! Arrays with NaN (start, middle, end) & all-NaN
    r_nan = [rnan_val, 3.0_dp, -2.0_dp, rnan_val, 0.0_dp, 3.0_dp, -5.0_dp, rnan_val]
    allocate(p_nan(size(r_nan)))
    do k = 1, 8
      r_out = sort(r_nan, P=p_nan, by=r_keys(k))
      call check_real_sort(r_nan, r_out, r_keys(k), fail, p=p_nan, num_nan=3)
    end do
    deallocate(p_nan)

    r_allnan = rnan_val
    allocate(p_nan(size(r_allnan)))
    r_out = sort(r_allnan, P=p_nan, by='RF')
    call check_real_sort(r_allnan, r_out, 'RF', fail, p=p_nan, num_nan=4)
    deallocate(p_nan)

    c_nan = [NaN(), cmplx(2.0_dp, -1.0_dp, dp), cmplx(-3.0_dp, 4.0_dp, dp), NaN(), &
             cmplx(0.0_dp, 0.0_dp, dp), cmplx(-3.0_dp, -2.0_dp, dp), NaN()]
    allocate(p_nan(size(c_nan)))
    do k = 1, 6
      c_out = sort(c_nan, P=p_nan, by=c_keys(k))
      call check_cmplx_sort(c_nan, c_out, c_keys(k), fail, p=p_nan, num_nan=3)
    end do
    deallocate(p_nan)

    c_allnan = NaN()
    allocate(p_nan(size(c_allnan)))
    c_out = sort(c_allnan, P=p_nan, by='MF')
    call check_cmplx_sort(c_allnan, c_out, 'MF', fail, p=p_nan, num_nan=4)
    deallocate(p_nan)

    ! All-equal arrays
    r_eq = 7.5_dp
    allocate(p(size(r_eq)))
    r_out = sort(r_eq, P=p, by='RF')
    call check_real_sort(r_eq, r_out, 'RF', fail, p=p)
    c_eq = cmplx(1.5_dp, -2.5_dp, dp)
    c_out = sort(c_eq, P=p, by='MF')
    call check_cmplx_sort(c_eq, c_out, 'MF', fail, p=p)
    deallocate(p)

    ! Large N=50000 random, already-sorted, and reverse-sorted arrays
    allocate(r_big(50000), c_big(50000), p_big(50000))
    call random_number(r_big, sigma=10.0_dp)
    call random_number(c_big, sigma=10.0_dp)
    r_out = sort(r_big, P=p_big, by='MR')
    call check_real_sort(r_big, r_out, 'MR', fail, p=p_big)
    r_out = sort(r_out, P=p_big, by='MR')
    call check_real_sort(r_out, r_out, 'MR', fail)
    r_out = sort(r_out, P=p_big, by='MF')
    call check_real_sort(r_out, r_out, 'MF', fail)

    c_out = sort(c_big, P=p_big, by='RF')
    call check_cmplx_sort(c_big, c_out, 'RF', fail, p=p_big)
    c_out = sort(c_out, P=p_big, by='RR')
    call check_cmplx_sort(c_out, c_out, 'RR', fail)
    deallocate(r_big, c_big, p_big)
  end subroutine test_bug04

  function to_upper2(s) result(u)
    character(2), intent(in) :: s
    character(2) :: u
    integer :: i
    u = s
    do i = 1, 2
      if (iachar(u(i:i)) >= iachar('a') .and. iachar(u(i:i)) <= iachar('z')) &
        u(i:i) = achar(iachar(u(i:i)) - 32)
    end do
  end function to_upper2

  subroutine check_perm(n, p, fail)
    integer, intent(in)    :: n, p(:)
    integer, intent(inout) :: fail
    logical :: seen(n)
    integer :: i
    if (size(p) /= n) then
      fail = fail + 1
      return
    end if
    seen = .false.
    do i = 1, n
      if (p(i) < 1 .or. p(i) > n) then
        fail = fail + 1
        return
      end if
      if (seen(p(i))) then
        fail = fail + 1
        return
      end if
      seen(p(i)) = .true.
    end do
  end subroutine check_perm

  subroutine check_real_sort(cin, cout, by, fail, p, num_nan)
    real(dp), intent(in) :: cin(:), cout(:)
    character(2), intent(in) :: by
    integer, intent(inout) :: fail
    integer, optional, intent(in) :: p(:), num_nan
    character(2) :: u
    integer :: i, m, n, expected_nan
    real(dp) :: k1, k2

    n = size(cin)
    expected_nan = 0
    if (present(num_nan)) expected_nan = num_nan
    m = n - expected_nan

    if (size(cout) /= n) then
      fail = fail + 1
      return
    end if
    u = to_upper2(by)

    do i = 1, m
      if (ieee_is_nan(cout(i))) fail = fail + 1
    end do
    do i = m + 1, n
      if (.not. ieee_is_nan(cout(i))) fail = fail + 1
    end do

    do i = 1, m - 1
      select case (u(1:1))
      case('V','R')
        k1 = cout(i); k2 = cout(i+1)
      case('A','M')
        k1 = abs(cout(i)); k2 = abs(cout(i+1))
      end select
      if (u(2:2) == 'F') then
        if (k1 > k2) fail = fail + 1
      else
        if (k1 < k2) fail = fail + 1
      end if
    end do

    if (present(p)) then
      call check_perm(n, p, fail)
      do i = 1, m
        if (cout(i) /= cin(p(i))) fail = fail + 1
      end do
      do i = m + 1, n
        if (.not. ieee_is_nan(cin(p(i)))) fail = fail + 1
      end do
    end if
  end subroutine check_real_sort

  subroutine check_cmplx_sort(cin, cout, by, fail, p, num_nan)
    complex(dp), intent(in) :: cin(:), cout(:)
    character(2), intent(in) :: by
    integer, intent(inout) :: fail
    integer, optional, intent(in) :: p(:), num_nan
    character(2) :: u
    integer :: i, m, n, expected_nan
    real(dp) :: k1, k2

    n = size(cin)
    expected_nan = 0
    if (present(num_nan)) expected_nan = num_nan
    m = n - expected_nan

    if (size(cout) /= n) then
      fail = fail + 1
      return
    end if
    u = to_upper2(by)

    do i = 1, m
      if (isNaN(cout(i))) fail = fail + 1
    end do
    do i = m + 1, n
      if (.not. isNaN(cout(i))) fail = fail + 1
    end do

    do i = 1, m - 1
      select case (u(1:1))
      case('M')
        k1 = abs(cout(i)); k2 = abs(cout(i+1))
      case('R')
        k1 = real(cout(i), dp); k2 = real(cout(i+1), dp)
      case('I')
        k1 = aimag(cout(i)); k2 = aimag(cout(i+1))
      end select
      if (u(2:2) == 'F') then
        if (k1 > k2) fail = fail + 1
      else
        if (k1 < k2) fail = fail + 1
      end if
    end do

    if (present(p)) then
      call check_perm(n, p, fail)
      do i = 1, m
        if (cout(i) /= cin(p(i))) fail = fail + 1
      end do
      do i = m + 1, n
        if (.not. isNaN(cin(p(i)))) fail = fail + 1
      end do
    end if
  end subroutine check_cmplx_sort

end program test_bug02_03_04
