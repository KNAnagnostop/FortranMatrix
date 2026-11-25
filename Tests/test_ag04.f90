!
! Verification test for BUG-08, BUG-09, BUG-10, BUG-11, MINOR-01, MINOR-03, MINOR-04
! gfortran -fopenmp -I../lib/gnu/ test_ag04.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t
!
program test_bug08_11_minor
  use array_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  integer, parameter :: dp = 8
  integer :: nfail, i, j, fd1, fd2, un, ios
  character(256) :: line
  logical :: found_imag

  ! BUG-08 & BUG-10 arrays
  complex(dp) :: z35(3,5), z53(5,3), z33_a(3,3), z33_b(3,3), z55_a(5,5), z55_b(5,5)
  complex(dp) :: z03(0,3), z30(3,0), z00(0,0)
  complex(dp), allocatable :: zres(:,:), zbad(:,:)
  real(dp)    :: d35(3,5), d53(5,3), d33_a(3,3), d33_b(3,3), d55_a(5,5), d55_b(5,5)
  real(dp)    :: d03(0,3), d30(3,0), d00(0,0), dbad(2,2)
  real(dp), allocatable    :: dres(:,:)
  real(dp)    :: rnd1, rnd2

  ! BUG-11 arrays
  complex(dp) :: zh(4,4), zg(4,4), z1(1,1)
  real(dp)    :: ds(4,4), dg(4,4), d1(1,1)
  complex(dp), allocatable :: ev1d(:), ev2d(:,:), zin(:,:), z0v(:), zw(:)
  real(dp), allocatable    :: din(:,:), d0v(:), dw(:)
  real(dp)    :: z33r(3,3), v3(3), err

  ! MINOR-04 arrays
  complex(dp) :: zsq(3,3), zrect53(5,3), zrect35(3,5)
  real(dp)    :: dsq(3,3), drect53(5,3), drect35(3,5)
  complex(dp) :: z1_0(0), z3_0(0,0,0)
  real(dp)    :: d1_0(0), d3_0(0,0,0)

  nfail = 0

  ! ------------------------------------------------------------------
  ! 1. BUG-09: matrix_random_init FD leak & signed overflow safety
  ! ------------------------------------------------------------------
  call execute_command_line('ls -1 /proc/self/fd | wc -l > fd_before.txt')
  do i = 1, 50
    call matrix_random_init()
  end do
  call execute_command_line('ls -1 /proc/self/fd | wc -l > fd_after.txt')
  open(newunit=un, file='fd_before.txt', status='old', action='read')
  read(un, *) fd1
  close(un)
  open(newunit=un, file='fd_after.txt', status='old', action='read')
  read(un, *) fd2
  close(un)
  call check(fd1 == fd2, 'BUG-09: matrix_random_init does not leak FDs')
  call random_number(rnd1)
  call matrix_random_init()
  call random_number(rnd2)
  call check(rnd1 >= 0.0_dp .and. rnd1 < 1.0_dp .and. rnd2 >= 0.0_dp .and. rnd2 < 1.0_dp, &
             'BUG-09: random_number works after matrix_random_init')

  ! ------------------------------------------------------------------
  ! 2. BUG-08: zherk interface (real alpha, beta) & BLAS lmatmul/mmmult
  ! ------------------------------------------------------------------
  do j = 1, 5
    do i = 1, 3
      z35(i,j) = cmplx(real(i + 2*j, dp), real(3*i - j, dp), kind=dp)
      d35(i,j) = real(i + 2*j, dp) - 0.5_dp * real(j, dp)
    end do
  end do
  z53 = transpose(z35) + cmplx(0.0_dp, 1.0_dp, kind=dp)
  d53 = transpose(d35) + 1.25_dp

  z33_a = lmatmul(z35)
  z33_b = matmul(z35, conjg(transpose(z35)))
  call check(maxval(abs(z33_a - z33_b)) < 1.0e-12_dp, 'BUG-08: lmatmul(z35) == Z * Z^H')

  z55_a = lmatmul(z35, side='H')
  z55_b = matmul(conjg(transpose(z35)), z35)
  call check(maxval(abs(z55_a - z55_b)) < 1.0e-12_dp, 'BUG-08: lmatmul(z35, side=H) == Z^H * Z')

  call mmmult(z53, z55_a)
  z55_b = matmul(z53, conjg(transpose(z53)))
  call check(maxval(abs(z55_a - z55_b)) < 1.0e-12_dp, 'BUG-08: mmmult(z53, C) == Z * Z^H')

  call mmmult(z53, z33_a, side='h')
  z33_b = matmul(conjg(transpose(z53)), z53)
  call check(maxval(abs(z33_a - z33_b)) < 1.0e-12_dp, 'BUG-08: mmmult(z53, C, side=h) == Z^H * Z')

  ! Empty side string safety & mismatch check
  z33_a = lmatmul(z35, side='')
  call check(maxval(abs(z33_a - matmul(z35, conjg(transpose(z35))))) < 1.0e-12_dp, &
             'BUG-08: lmatmul(z35, side="") defaults safely')
  allocate(zbad(2,2))
  zbad = (0.0_dp, 0.0_dp)
  call mmmult(z35, zbad)
  call check(isNaN(zbad) .and. all(ieee_is_nan(real(zbad, dp))), 'BUG-08: mmmult complex shape mismatch gives NaN')

  ! 0-sized lmatmul/mmmult (BUG-08 & BUG-11)
  zres = lmatmul(z03)
  call check(all(shape(zres) == [0, 0]), 'BUG-08/11: lmatmul(0x3) shape [0,0]')
  zres = lmatmul(z03, side='H')
  call check(all(shape(zres) == [3, 3]) .and. maxval(abs(zres)) == 0.0_dp, &
             'BUG-08/11: lmatmul(0x3, side=H) gives 3x3 zero')
  zres = lmatmul(z30)
  call check(all(shape(zres) == [3, 3]) .and. maxval(abs(zres)) == 0.0_dp, &
             'BUG-08/11: lmatmul(3x0) gives 3x3 zero')
  zres = lmatmul(z00)
  call check(all(shape(zres) == [0, 0]), 'BUG-08/11: lmatmul(0x0) shape [0,0]')

  ! ------------------------------------------------------------------
  ! 3. BUG-10: real mmmult(D, R, side='S'/'T'/'C') (opA='T' vs 'C')
  ! ------------------------------------------------------------------
  d33_a = lmatmul(d35)
  call mmmult(d35, d33_b)
  call check(maxval(abs(d33_a - matmul(d35, transpose(d35)))) < 1.0e-12_dp .and. &
             maxval(abs(d33_b - d33_a)) < 1.0e-12_dp, 'BUG-10: mmmult(d35, R) == D * D^T')

  d55_a = lmatmul(d35, side='S')
  call mmmult(d35, d55_b, side='C')
  call check(maxval(abs(d55_a - matmul(transpose(d35), d35))) < 1.0e-12_dp .and. &
             maxval(abs(d55_b - d55_a)) < 1.0e-12_dp, 'BUG-10: mmmult(d35, R, side=C) == D^T * D')

  call mmmult(d35, d55_b, side='T')
  call check(maxval(abs(d55_b - d55_a)) < 1.0e-12_dp, 'BUG-10: mmmult(d35, R, side=T) == D^T * D')

  dres = lmatmul(d03, side='T')
  call check(all(shape(dres) == [3, 3]) .and. maxval(abs(dres)) == 0.0_dp, &
             'BUG-10/11: lmatmul(d 0x3, side=T) gives 3x3 zero')
  call mmmult(d30, dres)
  call check(all(shape(dres) == [3, 3]) .and. maxval(abs(dres)) == 0.0_dp, &
             'BUG-10/11: mmmult(d 3x0, R) gives 3x3 zero')
  dbad = 0.0_dp
  call mmmult(d35, dbad)
  call check(isNaN(dbad) .and. all(ieee_is_nan(dbad)), 'BUG-10/MINOR-03: mmmult real shape mismatch gives NaN')

  ! ------------------------------------------------------------------
  ! 4. BUG-11: eigenvalues, eigenvectors, inverse, 0-sized BLAS wrappers
  ! ------------------------------------------------------------------
  ! 4a. Real general with conjugate pairs AND real eigenvalues:
  !     Block 1: [[1,-2],[2,1]] -> 1 +/- 2i
  !     Block 2: [[3,-4],[4,3]] -> 3 +/- 4i
  dg = 0.0_dp
  dg(1,1) = 1.0_dp; dg(1,2) = -2.0_dp; dg(2,1) = 2.0_dp; dg(2,2) = 1.0_dp
  dg(3,3) = 3.0_dp; dg(3,4) = -4.0_dp; dg(4,3) = 4.0_dp; dg(4,4) = 3.0_dp
  dg(1,3) = 0.25_dp; dg(2,4) = -0.1_dp

  ev1d = sort(eigenvalues(dg, mtype='DG'), by='RF')
  ev2d = eigenvectors(dg, mtype='DG', sortby='RF')
  call check(maxval(abs(ev1d - ev2d(:,1))) < 1.0e-12_dp, 'BUG-11: DG eigenvalues match sorted eigenvectors(:,1)')
  err = maxval(abs(matmul(cmplx(dg, 0.0_dp, kind=dp), ev2d(:,2:)) - matmul(ev2d(:,2:), diagonal(ev2d(:,1)))))
  call check(err < 1.0e-12_dp, 'BUG-11: DG A*V == V*diag(E) with complex pairs')

  ! Also test DG ending in a conjugate pair after real eigenvalues: [5, -2, 0+/-3i]
  dg = 0.0_dp
  dg(1,1) = 5.0_dp; dg(2,2) = -2.0_dp; dg(3,4) = -3.0_dp; dg(4,3) = 3.0_dp
  dg(1,4) = 0.3_dp
  ev2d = eigenvectors(dg, mtype='DG', sortby='IF')
  err = maxval(abs(matmul(cmplx(dg, 0.0_dp, kind=dp), ev2d(:,2:)) - matmul(ev2d(:,2:), diagonal(ev2d(:,1)))))
  call check(err < 1.0e-12_dp, 'BUG-11: DG mixed real + trailing conjugate pair A*V == V*diag(E)')

  ! 4b. Real symmetric (DS)
  ds = reshape([4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                1.0_dp, 3.0_dp, 1.0_dp, 0.0_dp, &
                0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, 1.0_dp, 5.0_dp], [4,4])
  ev1d = eigenvalues(ds, mtype='DS', sortby='RF')
  ev2d = eigenvectors(ds, mtype='DS', sortby='RF')
  err = maxval(abs(matmul(cmplx(ds, 0.0_dp, kind=dp), ev2d(:,2:)) - matmul(ev2d(:,2:), diagonal(ev2d(:,1)))))
  call check(maxval(abs(ev1d - ev2d(:,1))) < 1.0e-12_dp .and. err < 1.0e-12_dp, 'BUG-11: DS eigensystem')

  ! 4c. Complex general (ZG) & Hermitian (ZH)
  zg = cmplx(dg, 0.5_dp * transpose(dg), kind=dp)
  ev1d = sort(eigenvalues(zg, mtype='ZG'), by='MF')
  ev2d = eigenvectors(zg, mtype='ZG', sortby='MF')
  err = maxval(abs(matmul(zg, ev2d(:,2:)) - matmul(ev2d(:,2:), diagonal(ev2d(:,1)))))
  call check(maxval(abs(ev1d - ev2d(:,1))) < 1.0e-12_dp .and. err < 1.0e-12_dp, 'BUG-11: ZG eigensystem')

  zh = zg
  call hermitian_set(zh)
  ev1d = eigenvalues(zh, mtype='ZH', sortby='RF')
  ev2d = eigenvectors(zh, mtype='ZH', sortby='RF')
  err = maxval(abs(matmul(zh, ev2d(:,2:)) - matmul(ev2d(:,2:), diagonal(ev2d(:,1)))))
  call check(maxval(abs(ev1d - ev2d(:,1))) < 1.0e-12_dp .and. err < 1.0e-12_dp, 'BUG-11: ZH eigensystem')

  ! 4d. 1x1 edge case (zheev RWORK=max(1, 3*n-2)=1)
  z1(1,1) = (3.5_dp, 0.0_dp)
  d1(1,1) = -2.25_dp
  ev2d = eigenvectors(z1, mtype='ZH')
  call check(abs(ev2d(1,1) - (3.5_dp, 0.0_dp)) < 1.0e-13_dp .and. abs(abs(ev2d(1,2)) - 1.0_dp) < 1.0e-13_dp, &
             'BUG-11: 1x1 ZH eigenvectors')
  ev1d = eigenvalues(d1, mtype='DG')
  call check(abs(ev1d(1) - (-2.25_dp, 0.0_dp)) < 1.0e-13_dp, 'BUG-11: 1x1 DG eigenvalues')

  ! 4e. 0x0 eigensolvers & inverse
  ev1d = eigenvalues(z00, mtype='ZG', sortby='RF')
  call check(size(ev1d) == 0, 'BUG-11: eigenvalues(z00, ZG) size 0')
  ev1d = eigenvalues(z00, mtype='ZH')
  call check(size(ev1d) == 0, 'BUG-11: eigenvalues(z00, ZH) size 0')
  ev1d = eigenvalues(d00, mtype='DG')
  call check(size(ev1d) == 0, 'BUG-11: eigenvalues(d00, DG) size 0')
  ev1d = eigenvalues(d00, mtype='DS', sortby='RF')
  call check(size(ev1d) == 0, 'BUG-11: eigenvalues(d00, DS) size 0')

  ev2d = eigenvectors(z00, mtype='ZG', sortby='RF')
  call check(all(shape(ev2d) == [0, 1]), 'BUG-11: eigenvectors(z00, ZG) shape [0,1]')
  ev2d = eigenvectors(z00, mtype='ZH')
  call check(all(shape(ev2d) == [0, 1]), 'BUG-11: eigenvectors(z00, ZH) shape [0,1]')
  ev2d = eigenvectors(d00, mtype='DG', sortby='RF')
  call check(all(shape(ev2d) == [0, 1]), 'BUG-11: eigenvectors(d00, DG) shape [0,1]')
  ev2d = eigenvectors(d00, mtype='DS')
  call check(all(shape(ev2d) == [0, 1]), 'BUG-11: eigenvectors(d00, DS) shape [0,1]')

  zin = inverse(z00)
  call check(all(shape(zin) == [0, 0]), 'BUG-11: inverse(z00) shape [0,0]')
  din = inverse(d00)
  call check(all(shape(din) == [0, 0]), 'BUG-11: inverse(d00) shape [0,0]')
  zin = inverse(zh)
  call check(maxval(abs(matmul(zh, zin) - cidentitymatrix(4))) < 1.0e-12_dp, 'BUG-11: inverse(4x4 complex) accurate')
  din = inverse(ds)
  call check(maxval(abs(matmul(ds, din) - didentitymatrix(4))) < 1.0e-12_dp, 'BUG-11: inverse(4x4 real) accurate')

  ! 4f. 0-sized matrix/vector BLAS wrappers
  zres = lmatmul(z30, z03)
  call check(all(shape(zres) == [3, 3]) .and. maxval(abs(zres)) == 0.0_dp, 'BUG-11: lmatmul(z 3x0, 0x3) -> 3x3 zero')
  dres = lmatmul(d30, d03)
  call check(all(shape(dres) == [3, 3]) .and. maxval(abs(dres)) == 0.0_dp, 'BUG-11: lmatmul(d 3x0, 0x3) -> 3x3 zero')
  allocate(z0v(0), zw(3), d0v(0), dw(3))
  zw = lmatmul(z30, z0v, type='N')
  call check(size(zw) == 3 .and. maxval(abs(zw)) == 0.0_dp, 'BUG-11: lmatmul(z 3x0, v0) -> zero')
  call mvmult(z03, z0v, zw, type='T')
  call check(maxval(abs(zw)) == 0.0_dp, 'BUG-11: mvmult(z 0x3, v0, type=T) -> zero')
  call vmmult(z0v, z03, zw)
  call check(maxval(abs(zw)) == 0.0_dp, 'BUG-11: vmmult(v0, z 0x3) -> zero')
  dw = lmatmul(d30, d0v, type='')
  call check(isNaN(dw) .and. all(ieee_is_nan(dw)), 'BUG-11: lmatmul(d, v, type="") -> NaN')
  zw = lmatmul(z30, z0v, type='')
  call check(isNaN(zw) .and. all(ieee_is_nan(real(zw, dp))), 'BUG-11: lmatmul(z, v, type="") -> NaN')
  call vmmult(d0v, d03, dw)
  call check(maxval(abs(dw)) == 0.0_dp, 'BUG-11: vmmult(v0, d 0x3) -> zero')
  deallocate(z0v, zw, d0v, dw)

  ! ------------------------------------------------------------------
  ! 5. MINOR-01: array2_print uses standard AIMAG cleanly
  ! ------------------------------------------------------------------
  zbad = reshape([(1.5_dp, -2.75_dp), (0.0_dp, 3.25_dp), (-4.5_dp, 1.0_dp), (2.0_dp, 0.0_dp)], [2,2])
  open(newunit=un, file='print_complex.txt', status='replace', action='readwrite')
  call print(zbad, unit=un, form='Mathematica', name='zm')
  rewind(un)
  found_imag = .false.
  do
    read(un, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (index(line, '-2.75') > 0) found_imag = .true.
  end do
  close(un)
  call check(found_imag, 'MINOR-01: array2_print Mathematica output contains AIMAG values')

  ! ------------------------------------------------------------------
  ! 6. MINOR-04: 0-sized & non-square guards in pure basic array procs
  ! ------------------------------------------------------------------
  ! 6a. 0-sized
  call check(norm(z1_0) == 0.0_dp .and. norm(z00) == 0.0_dp .and. norm(z3_0) == 0.0_dp, &
             'MINOR-04: complex norm on 0-sized arrays is 0.0')
  call check(norm(d1_0) == 0.0_dp .and. norm(d00) == 0.0_dp .and. norm(d3_0) == 0.0_dp, &
             'MINOR-04: real norm on 0-sized arrays is 0.0')
  call check(abs(trace(z00)) == 0.0_dp .and. trace(d00) == 0.0_dp, 'MINOR-04: trace on 0x0 is 0')
  call check(abs(trace2(z00)) == 0.0_dp .and. trace2(d00) == 0.0_dp, 'MINOR-04: trace2 on 0x0 is 0')
  call check(abs(trace2c(z00)) == 0.0_dp .and. trace2c(d00) == 0.0_dp, 'MINOR-04: trace2c on 0x0 is 0')
  call check(all(shape(traceless(z00)) == [0,0]) .and. all(shape(traceless(d00)) == [0,0]), &
             'MINOR-04: traceless on 0x0 gives [0,0]')
  call traceless_set(z00)
  call traceless_set(d00)
  call check(size(diagonal(z00)) == 0 .and. size(diagonal(d00)) == 0, 'MINOR-04: diagonal on 0x0 size 0')
  call check(all(shape(diagonalMatrix(z00)) == [0,0]) .and. all(shape(diagonalMatrix(d00)) == [0,0]), &
             'MINOR-04: diagonalMatrix(0x0 matrix) shape [0,0]')
  call check(all(shape(diagonal(1.0_dp, 0)) == [0,0]) .and. all(shape(diagonal(1.0_dp, -3)) == [0,0]), &
             'MINOR-04: diagonal(r, n<=0) shape [0,0]')
  call check(all(shape(diagonal((1.0_dp, 2.0_dp), -2)) == [0,0]), &
             'MINOR-04: diagonal(z, n<=0) shape [0,0]')
  call check(all(shape(cidentitymatrix(-1)) == [0,0]) .and. all(shape(didentitymatrix(0)) == [0,0]), &
             'MINOR-04: cidentitymatrix/didentitymatrix(n<=0) shape [0,0]')
  call hermitian_set(z00); call symmetric_set(z00); call antisymmetric_set(z00)
  call symmetric_set(d00); call antisymmetric_set(d00)
  call check(isHermitian(z00) == 0.0_dp .and. isSymmetric(z00) == 0.0_dp .and. isAntiSymmetric(z00) == 0.0_dp, &
             'MINOR-04: isHermitian/isSymmetric/isAntiSymmetric on z00 is 0')
  call check(isSymmetric(d00) == 0.0_dp .and. isAntiSymmetric(d00) == 0.0_dp, &
             'MINOR-04: isSymmetric/isAntiSymmetric on d00 is 0')

  ! 6b. Non-square 5x3 & 3x5
  do j = 1, 3
    do i = 1, 5
      zrect53(i,j) = cmplx(real(10*i + j, dp), real(i - j, dp), kind=dp)
      drect53(i,j) = real(10*i + j, dp)
    end do
  end do
  zrect35 = z35
  drect35 = d35

  call check(size(diagonal(zrect53)) == 3 .and. &
             maxval(abs(diagonal(zrect53) - [zrect53(1,1), zrect53(2,2), zrect53(3,3)])) == 0.0_dp, &
             'MINOR-04: diagonal(z 5x3) size 3 & exact elements')
  call check(size(diagonal(drect53)) == 3 .and. &
             maxval(abs(diagonal(drect53) - [drect53(1,1), drect53(2,2), drect53(3,3)])) == 0.0_dp, &
             'MINOR-04: diagonal(d 5x3) size 3 & exact elements')
  zres = diagonalMatrix(zrect53)
  dres = diagonalMatrix(drect53)
  call check(all(shape(zres) == [5,3]) .and. all(shape(dres) == [5,3]) .and. &
             abs(zres(2,2) - zrect53(2,2)) == 0.0_dp .and. abs(zres(5,1)) == 0.0_dp .and. &
             abs(dres(3,3) - drect53(3,3)) == 0.0_dp .and. abs(dres(4,2)) == 0.0_dp, &
             'MINOR-04: diagonalMatrix(5x3) preserves shape & zeros off-diagonal')

  call check(isNaN(trace(zrect53)) .and. ieee_is_nan(trace(drect53)), &
             'MINOR-04: trace(non-square) returns NaN')
  call check(isNaN(trace2c(zrect53)) .and. ieee_is_nan(trace2c(drect53)), &
             'MINOR-04: trace2c(non-square) returns NaN')
  call check(ieee_is_nan(isHermitian(zrect53)) .and. ieee_is_nan(isSymmetric(zrect53)) .and. &
             ieee_is_nan(isAntiSymmetric(zrect53)) .and. ieee_is_nan(isSymmetric(drect35)) .and. &
             ieee_is_nan(isAntiSymmetric(drect35)), &
             'MINOR-04: isHermitian/isSymmetric/isAntiSymmetric(non-square) returns NaN')

  zres = traceless(zrect53)
  dres = traceless(drect35)
  call check(all(shape(zres) == [5,3]) .and. all(ieee_is_nan(real(zres, dp))) .and. &
             all(shape(dres) == [3,5]) .and. all(ieee_is_nan(dres)), &
             'MINOR-04: traceless(non-square) returns NaN array of same shape')

  call traceless_set(zrect53); call check(all(ieee_is_nan(real(zrect53, dp))), 'MINOR-04: traceless_set(z non-square) -> NaN')
  call traceless_set(drect35); call check(all(ieee_is_nan(drect35)), 'MINOR-04: traceless_set(d non-square) -> NaN')

  call hermitian_set(zrect53);     call check(all(ieee_is_nan(real(zrect53, dp))), 'MINOR-04: hermitian_set(z non-square) -> NaN')
  call symmetric_set(zrect35);     call check(all(ieee_is_nan(real(zrect35, dp))), 'MINOR-04: symmetric_set(z non-square) -> NaN')
  call antisymmetric_set(zrect35); call check(all(ieee_is_nan(real(zrect35, dp))), 'MINOR-04: antisymmetric_set(z non-square) -> NaN')
  call symmetric_set(drect53);     call check(all(ieee_is_nan(drect53)), 'MINOR-04: symmetric_set(d non-square) -> NaN')
  call antisymmetric_set(drect53); call check(all(ieee_is_nan(drect53)), 'MINOR-04: antisymmetric_set(d non-square) -> NaN')

  ! 6c. Empty uplo='' on square matrix
  zsq = reshape([(1.0_dp, 2.0_dp), (3.0_dp, 4.0_dp), (5.0_dp, -1.0_dp), &
                 (-2.0_dp, 1.0_dp), (0.0_dp, 2.0_dp), (4.0_dp, 0.5_dp), &
                 (1.0_dp, -3.0_dp), (2.0_dp, 1.0_dp), (6.0_dp, -2.0_dp)], [3,3])
  call hermitian_set(zsq, uplo='')
  call check(isHermitian(zsq) < 1.0e-13_dp, 'MINOR-04: hermitian_set(zsq, uplo="") succeeds')
  call symmetric_set(zsq, uplo='')
  call check(isSymmetric(zsq) < 1.0e-13_dp, 'MINOR-04: symmetric_set(zsq, uplo="") succeeds')
  call antisymmetric_set(zsq, uplo='')
  call check(isAntiSymmetric(zsq) < 1.0e-13_dp, 'MINOR-04: antisymmetric_set(zsq, uplo="") succeeds')

  dsq = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp], [3,3])
  call symmetric_set(dsq, uplo='')
  call check(isSymmetric(dsq) < 1.0e-13_dp, 'MINOR-04: symmetric_set(dsq, uplo="") succeeds')
  call antisymmetric_set(dsq, uplo='')
  call check(isAntiSymmetric(dsq) < 1.0e-13_dp, 'MINOR-04: antisymmetric_set(dsq, uplo="") succeeds')

  ! Keep dummy v3/z33r unused check clean
  z33r = 0.0_dp; v3 = 0.0_dp
  if (sum(z33r) + sum(v3) > 1.0_dp) nfail = nfail + 1

  write(*,'(A,I0)') 'TOTAL_FAILURES = ', nfail
  if (nfail > 0) stop 1

contains
  subroutine check(cond, label)
    logical, intent(in) :: cond
    character(*), intent(in) :: label
    if (.not. cond) then
      write(*,'(A,A)') 'FAIL: ', label
      nfail = nfail + 1
    end if
  end subroutine check
end program test_bug08_11_minor
