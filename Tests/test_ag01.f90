!
! Verification test for BUG-01 (real mvmult/lmatmul dispatch to dsymv, 0-size & empty type guards)
! gfortran -fopenmp -I../lib/gnu/ test_ag01.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t
!
program test_bug01
  use array_mod
  implicit none
  integer, parameter :: dp = 8
  real(dp) :: A_sym(4,4), A_upper(4,4), A_rect(3,4)
  real(dp) :: v4(4), v3(3), w4(4), w3(3), ref4(4), ref3(3), empty_v(0), empty_m00(0,0), empty_m04(0,4), empty_m40(4,0)
  real(dp), allocatable :: res(:)
  real(dp) :: err, tol
  integer :: i, j, failures
  character(len=1) :: types(10) = ['N','n','T','t','C','c','S','s','H','h']
  character(len=1) :: sym_types(4) = ['S','s','H','h']

  tol = 1.0e-14_dp
  failures = 0

  do j = 1, 4
    do i = 1, 4
      A_sym(i,j) = real(i + 2*j, dp)
    end do
  end do
  call symmetric_set(A_sym)
  v4 = [1.0_dp, -2.0_dp, 3.5_dp, -0.5_dp]
  ref4 = matmul(A_sym, v4)

  ! Matrix with matching upper triangle and poisoned strictly lower triangle
  A_upper = A_sym
  do j = 1, 4
    do i = j + 1, 4
      A_upper(i,j) = 999.0_dp
    end do
  end do

  ! Rectangular matrix for N, T, C and error cases
  do j = 1, 4
    do i = 1, 3
      A_rect(i,j) = real(3*i - j, dp)
    end do
  end do
  v3 = [2.0_dp, -1.0_dp, 4.0_dp]

  ! 1. Symmetric A_sym with all valid types (lmatmul & mvmult)
  do i = 1, size(types)
    res = lmatmul(A_sym, v4, type=types(i))
    err = maxval(abs(res - ref4))
    if (err > tol .or. isNaN(res)) then
      print *, 'FAIL 1a lmatmul type=', types(i), ' err=', err
      failures = failures + 1
    end if

    w4 = 99.0_dp
    call mvmult(A_sym, v4, w4, type=types(i))
    err = maxval(abs(w4 - ref4))
    if (err > tol .or. isNaN(w4)) then
      print *, 'FAIL 1b mvmult type=', types(i), ' err=', err
      failures = failures + 1
    end if
  end do

  ! 2. Proof of dsymv('U', ...) dispatch: poisoned lower triangle ignored for S, s, H, h
  do i = 1, size(sym_types)
    res = lmatmul(A_upper, v4, type=sym_types(i))
    err = maxval(abs(res - ref4))
    if (err > tol .or. isNaN(res)) then
      print *, 'FAIL 2a dsymv proof lmatmul type=', sym_types(i), ' err=', err
      failures = failures + 1
    end if

    w4 = 99.0_dp
    call mvmult(A_upper, v4, w4, type=sym_types(i))
    err = maxval(abs(w4 - ref4))
    if (err > tol .or. isNaN(w4)) then
      print *, 'FAIL 2b dsymv proof mvmult type=', sym_types(i), ' err=', err
      failures = failures + 1
    end if
  end do

  ! 3. Rectangular N, T, C
  ref3 = matmul(A_rect, v4)
  res = lmatmul(A_rect, v4, type='N')
  if (maxval(abs(res - ref3)) > tol) failures = failures + 1
  call mvmult(A_rect, v4, w3, type='n')
  if (maxval(abs(w3 - ref3)) > tol) failures = failures + 1

  ref4 = matmul(transpose(A_rect), v3)
  do i = 3, 6
    res = lmatmul(A_rect, v3, type=types(i))
    if (maxval(abs(res - ref4)) > tol) then
      print *, 'FAIL 3a rect lmatmul type=', types(i)
      failures = failures + 1
    end if
    w4 = 99.0_dp
    call mvmult(A_rect, v3, w4, type=types(i))
    if (maxval(abs(w4 - ref4)) > tol) then
      print *, 'FAIL 3b rect mvmult type=', types(i)
      failures = failures + 1
    end if
  end do

  ! 4. Non-square matrix with S, s, H, h -> NaN
  do i = 1, size(sym_types)
    res = lmatmul(A_rect, v4, type=sym_types(i))
    if (.not. isNaN(res)) then
      print *, 'FAIL 4a nonsquare lmatmul type=', sym_types(i)
      failures = failures + 1
    end if
    w3 = 0.0_dp
    call mvmult(A_rect, v4, w3, type=sym_types(i))
    if (.not. isNaN(w3)) then
      print *, 'FAIL 4b nonsquare mvmult type=', sym_types(i)
      failures = failures + 1
    end if
  end do

  ! 5. Invalid type='X' and empty type='' -> NaN
  res = lmatmul(A_sym, v4, type='X')
  if (.not. isNaN(res)) failures = failures + 1
  res = lmatmul(A_sym, v4, type='')
  if (.not. isNaN(res)) failures = failures + 1
  w4 = 0.0_dp
  call mvmult(A_sym, v4, w4, type='X')
  if (.not. isNaN(w4)) failures = failures + 1
  w4 = 0.0_dp
  call mvmult(A_sym, v4, w4, type='')
  if (.not. isNaN(w4)) failures = failures + 1

  ! 6. Dimension mismatches -> NaN
  res = lmatmul(A_sym, v3, type='N')
  if (.not. isNaN(res)) failures = failures + 1
  res = lmatmul(A_sym, v3, type='T')
  if (.not. isNaN(res)) failures = failures + 1
  res = lmatmul(A_sym, v3, type='S')
  if (.not. isNaN(res)) failures = failures + 1
  res = lmatmul(A_sym, v3, type='H')
  if (.not. isNaN(res)) failures = failures + 1
  call mvmult(A_sym, v3, w4, type='N')
  if (.not. isNaN(w4)) failures = failures + 1
  call mvmult(A_sym, v4, w3, type='N')
  if (.not. isNaN(w3)) failures = failures + 1
  call mvmult(A_sym, v3, w4, type='T')
  if (.not. isNaN(w4)) failures = failures + 1
  call mvmult(A_sym, v4, w3, type='T')
  if (.not. isNaN(w3)) failures = failures + 1
  call mvmult(A_sym, v3, w4, type='S')
  if (.not. isNaN(w4)) failures = failures + 1
  call mvmult(A_sym, v4, w3, type='S')
  if (.not. isNaN(w3)) failures = failures + 1

  ! 7. 0-sized matrices (0x0, 0x4, 4x0)
  res = lmatmul(empty_m00, empty_v, type='N')
  if (size(res) /= 0) failures = failures + 1
  res = lmatmul(empty_m00, empty_v, type='S')
  if (size(res) /= 0) failures = failures + 1
  call mvmult(empty_m00, empty_v, empty_v, type='S')
  res = lmatmul(empty_m04, v4, type='N')
  if (size(res) /= 0) failures = failures + 1
  res = lmatmul(empty_m40, empty_v, type='N')
  if (size(res) /= 4 .or. any(res /= 0.0_dp)) failures = failures + 1
  w4 = 99.0_dp
  call mvmult(empty_m40, empty_v, w4, type='N')
  if (any(w4 /= 0.0_dp)) failures = failures + 1
  res = lmatmul(empty_m40, v4, type='T')
  if (size(res) /= 0) failures = failures + 1
  res = lmatmul(empty_m04, empty_v, type='T')
  if (size(res) /= 4 .or. any(res /= 0.0_dp)) failures = failures + 1
  w4 = 99.0_dp
  call mvmult(empty_m04, empty_v, w4, type='T')
  if (any(w4 /= 0.0_dp)) failures = failures + 1

  print '(A,I0)', 'test_bug01 failures = ', failures
  if (failures /= 0) stop 1
end program test_bug01
