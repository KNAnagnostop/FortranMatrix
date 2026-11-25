!
! Verification test for BUG-05, BUG-06, MINOR-02 (determinant/lndet and pfaffian/lnPfaffian singular & odd-N handling)
! gfortran -fopenmp -I../lib/gnu/ test_ag03.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t
!
program test_bug05_06
 use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
 use, intrinsic :: ieee_arithmetic
 use :: matrix_mod_common, only: ONE, ZERO, IMU
 use :: array_mod
 use :: matrix_mod_array, only: array2_pfaffian2, array2_log_pfaffian2
 implicit none

 integer :: failures

 failures = 0

 call matrix_random_init()

 call test_det_and_lndet()
 call test_pfaffian_and_lnpf()

 if (failures == 0) then
  print '(A)', 'PASS: all BUG-05 and BUG-06 checks passed (0 failures).'
 else
  print '(A,I0,A)', 'FAIL: ', failures, ' checks failed.'
  stop 1
 end if

contains

 subroutine fail_msg(msg)
  character(*), intent(in) :: msg
  failures = failures + 1
  write(error_unit, '(A,A)') 'FAIL: ', trim(msg)
 end subroutine fail_msg

 subroutine check_singular_complex_pair(label, pair)
  character(*), intent(in) :: label
  complex(dp), intent(in) :: pair(2)
  complex(dp) :: recon

  recon = exp(pair(1)) * pair(2)
  if (ieee_is_nan(real(pair(1), dp)) .or. ieee_is_nan(aimag(pair(1))) .or. &
      ieee_is_nan(real(pair(2), dp)) .or. ieee_is_nan(aimag(pair(2))) .or. &
      ieee_is_nan(real(recon, dp))   .or. ieee_is_nan(aimag(recon))) then
   call fail_msg(label // ' produced NaN')
   return
  end if
#ifndef FAST_MATH_TEST
  if (ieee_class(real(pair(1), dp)) /= ieee_negative_inf) then
   call fail_msg(label // ' real(res(1)) is not -Inf')
  end if
#endif
  if (aimag(pair(1)) /= 0.0_dp) call fail_msg(label // ' aimag(res(1)) /= 0')
  if (pair(2) /= ZERO)          call fail_msg(label // ' res(2) /= (0,0)')
  if (recon /= ZERO)            call fail_msg(label // ' exp(res(1))*res(2) /= (0,0)')
 end subroutine check_singular_complex_pair

 subroutine check_singular_real_pair(label, pair)
  character(*), intent(in) :: label
  real(dp), intent(in) :: pair(2)
  real(dp) :: recon

  recon = exp(pair(1)) * pair(2)
  if (ieee_is_nan(pair(1)) .or. ieee_is_nan(pair(2)) .or. ieee_is_nan(recon)) then
   call fail_msg(label // ' produced NaN')
   return
  end if
#ifndef FAST_MATH_TEST
  if (ieee_class(pair(1)) /= ieee_negative_inf) then
   call fail_msg(label // ' res(1) is not -Inf')
  end if
#endif
  if (pair(2) /= 0.0_dp) call fail_msg(label // ' res(2) /= 0')
  if (recon /= 0.0_dp)   call fail_msg(label // ' exp(res(1))*res(2) /= 0')
 end subroutine check_singular_real_pair

 subroutine check_singular_complex_scalar(label, val)
  character(*), intent(in) :: label
  complex(dp), intent(in) :: val

  if (ieee_is_nan(real(val, dp)) .or. ieee_is_nan(aimag(val))) then
   call fail_msg(label // ' produced NaN')
   return
  end if
  if (val /= ZERO) call fail_msg(label // ' /= (0,0)')
 end subroutine check_singular_complex_scalar

 subroutine test_det_and_lndet()
  integer, parameter :: sizes(4) = [1, 2, 5, 16]
  integer :: k, n
  complex(dp), allocatable :: C(:,:), zblk(:,:)
  real(dp), allocatable :: D(:,:), dblk(:,:)
  complex(dp) :: zdet, zldet(2), zrecon, z0(0,0)
  real(dp) :: rdet, rldet(2), rrecon, r0(0,0), tol
  real(dp) :: dneg(2,2)

  ! 1. 0x0 matrices
  zdet = determinant(z0)
  if (zdet /= ONE) call fail_msg('determinant(z 0x0) /= (1,0)')
  zldet = lndet(z0)
  if (zldet(1) /= ZERO .or. zldet(2) /= ONE) call fail_msg('lndet(z 0x0) /= [(0,0),(1,0)]')

  rdet = determinant(r0)
  if (rdet /= 1.0_dp) call fail_msg('determinant(d 0x0) /= 1.0')
  rldet = lndet(r0)
  if (rldet(1) /= 0.0_dp .or. rldet(2) /= 1.0_dp) call fail_msg('lndet(d 0x0) /= [0.0, 1.0]')

  ! 2. Non-singular complex & real matrices
  do k = 1, size(sizes)
   n = sizes(k)
   allocate(C(n,n), D(n,n))
   call random_number(C, sigma=1.0_dp)
   call random_number(D, sigma=1.0_dp)

   zdet = determinant(C)
   zldet = lndet(C)
   zrecon = exp(zldet(1)) * zldet(2)
   tol = 1.0e-11_dp * max(1.0_dp, abs(zdet))
   if (abs(zdet - zrecon) > tol) call fail_msg('complex nonsing det vs lndet mismatch')
   if (abs(abs(zldet(2)) - 1.0_dp) > 1.0e-12_dp) call fail_msg('complex nonsing |lndet(2)| /= 1')
   if (aimag(zldet(1)) /= 0.0_dp) call fail_msg('complex nonsing aimag(lndet(1)) /= 0')

   rdet = determinant(D)
   rldet = lndet(D)
   rrecon = exp(rldet(1)) * rldet(2)
   tol = 1.0e-11_dp * max(1.0_dp, abs(rdet))
   if (abs(rdet - rrecon) > tol) call fail_msg('real nonsing det vs lndet mismatch')
   if (abs(abs(rldet(2)) - 1.0_dp) > 1.0e-12_dp) call fail_msg('real nonsing |lndet(2)| /= 1')

   deallocate(C, D)
  end do

  ! Negative real determinant check: [[1, 2], [3, 4]] -> det = -2
  dneg = reshape([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], [2,2])
  rdet = determinant(dneg)
  rldet = lndet(dneg)
  if (abs(rdet + 2.0_dp) > 1.0e-13_dp) call fail_msg('real 2x2 neg det /= -2')
  if (abs(rldet(2) + 1.0_dp) > 1.0e-13_dp) call fail_msg('real 2x2 neg lndet(2) /= -1')
  if (abs(exp(rldet(1))*rldet(2) + 2.0_dp) > 1.0e-13_dp) call fail_msg('real 2x2 neg lndet recon /= -2')

  ! Negative complex entry check via random scaling: i * I_2 -> det = -1
  allocate(C(2,2))
  C = ZERO
  C(1,1) = IMU; C(2,2) = IMU
  zdet = determinant(C)
  zldet = lndet(C)
  if (abs(zdet + ONE) > 1.0e-13_dp) call fail_msg('complex i*I2 det /= -1')
  if (abs(exp(zldet(1))*zldet(2) + ONE) > 1.0e-13_dp) call fail_msg('complex i*I2 lndet recon /= -1')
  deallocate(C)

  ! 3. Singular zero matrices (1x1, 2x2, 4x4)
  do k = 1, 3
   if (k == 1) n = 1
   if (k == 2) n = 2
   if (k == 3) n = 4
   allocate(C(n,n), D(n,n))
   C = ZERO
   D = 0.0_dp

   zdet = determinant(C)
   call check_singular_complex_scalar('determinant(zero complex)', zdet)
   call check_singular_complex_pair('lndet(zero complex)', lndet(C))

   rdet = determinant(D)
   if (ieee_is_nan(rdet) .or. rdet /= 0.0_dp) call fail_msg('determinant(zero real) /= 0')
   call check_singular_real_pair('lndet(zero real)', lndet(D))

   deallocate(C, D)
  end do

  ! 4. Exact singular 5x5 matrices:
  ! 4a. Zero row / zero column (exact info > 0 in LAPACK LU)
  n = 5
  allocate(C(n,n), D(n,n))
  call random_number(C, sigma=1.0_dp)
  call random_number(D, sigma=1.0_dp)
  C(4,:) = ZERO
  D(:,4) = 0.0_dp

  call check_singular_complex_scalar('determinant(zero row complex)', determinant(C))
  call check_singular_complex_pair('lndet(zero row complex)', lndet(C))
  rdet = determinant(D)
  if (ieee_is_nan(rdet) .or. rdet /= 0.0_dp) call fail_msg('determinant(zero col real) /= 0')
  call check_singular_real_pair('lndet(zero col real)', lndet(D))

  ! 4b. Exact dyadic rank-deficient 5x5 (row 4 = 2 * row 2 with exact powers of 2, so exact LU zero pivot)
  D = reshape([ 1.0_dp,  2.0_dp,  0.0_dp,  4.0_dp,  0.0_dp, &
                2.0_dp,  4.0_dp,  1.0_dp,  8.0_dp,  0.0_dp, &
                4.0_dp,  0.0_dp,  2.0_dp,  0.0_dp,  1.0_dp, &
                0.0_dp,  2.0_dp,  4.0_dp,  4.0_dp,  2.0_dp, &
                1.0_dp,  1.0_dp,  0.0_dp,  2.0_dp,  4.0_dp ], [5,5])
  C = cmplx(D, 0.0_dp, kind=dp)
  C(1,5) = IMU; C(2,3) = cmplx(0.0_dp, 1.0_dp, kind=dp); C(4,3) = cmplx(0.0_dp, 2.0_dp, kind=dp)
  call check_singular_complex_scalar('determinant(dyadic dup row complex)', determinant(C))
  call check_singular_complex_pair('lndet(dyadic dup row complex)', lndet(C))
  rdet = determinant(D)
  if (ieee_is_nan(rdet) .or. rdet /= 0.0_dp) call fail_msg('determinant(dyadic dup row real) /= 0')
  call check_singular_real_pair('lndet(dyadic dup row real)', lndet(D))
  deallocate(C, D)

  ! 5. Block-diagonal 5x5: nonsing 3x3 + zero 2x2, and reversed (zero 2x2 + nonsing 3x3)
  allocate(C(5,5), zblk(3,3), D(5,5), dblk(3,3))
  call random_number(zblk, sigma=1.0_dp)
  call random_number(dblk, sigma=1.0_dp)

  C = ZERO; C(1:3,1:3) = zblk
  D = 0.0_dp; D(1:3,1:3) = dblk
  call check_singular_complex_scalar('determinant(blk+0 complex)', determinant(C))
  call check_singular_complex_pair('lndet(blk+0 complex)', lndet(C))
  if (ieee_is_nan(determinant(D)) .or. determinant(D) /= 0.0_dp) call fail_msg('determinant(blk+0 real) /= 0')
  call check_singular_real_pair('lndet(blk+0 real)', lndet(D))

  C = ZERO; C(3:5,3:5) = zblk
  D = 0.0_dp; D(3:5,3:5) = dblk
  call check_singular_complex_scalar('determinant(0+blk complex)', determinant(C))
  call check_singular_complex_pair('lndet(0+blk complex)', lndet(C))
  if (ieee_is_nan(determinant(D)) .or. determinant(D) /= 0.0_dp) call fail_msg('determinant(0+blk real) /= 0')
  call check_singular_real_pair('lndet(0+blk real)', lndet(D))

  deallocate(C, zblk, D, dblk)
 end subroutine test_det_and_lndet

 subroutine make_random_skew(A)
  complex(dp), intent(out) :: A(:,:)
  call random_number(A, sigma=0.5_dp)
  call antisymmetric_set(A)
 end subroutine make_random_skew

 subroutine check_pf_all_variants_zero(label, A)
  character(*), intent(in) :: label
  complex(dp), contiguous, intent(in) :: A(:,:)

  call check_singular_complex_scalar(label // ' pfaffian(default)', pfaffian(A))
  call check_singular_complex_scalar(label // ' pfaffian(P)',       pfaffian(A, alg='P'))
  call check_singular_complex_scalar(label // ' pfaffian(C)',       pfaffian(A, alg='C'))
  call check_singular_complex_scalar(label // ' pfaffian("")',      pfaffian(A, alg=''))
  call check_singular_complex_scalar(label // ' array2_pfaffian2',  array2_pfaffian2(A))

  call check_singular_complex_pair(label // ' lnPfaffian(default)', lnPfaffian(A))
  call check_singular_complex_pair(label // ' lnPfaffian(P)',       lnPfaffian(A, alg='P'))
  call check_singular_complex_pair(label // ' lnPfaffian(C)',       lnPfaffian(A, alg='C'))
  call check_singular_complex_pair(label // ' lnPfaffian("")',      lnPfaffian(A, alg=''))
  call check_singular_complex_pair(label // ' array2_log_pfaffian2', array2_log_pfaffian2(A))
 end subroutine check_pf_all_variants_zero

 subroutine test_pfaffian_and_lnpf()
  integer, parameter :: even_sizes(4) = [2, 4, 8, 16]
  integer, parameter :: odd_sizes(4)  = [1, 3, 5, 7]
  integer :: k, n
  complex(dp), allocatable :: A(:,:), blk2(:,:), blk4(:,:)
  complex(dp) :: z0(0,0)
  complex(dp) :: pf_def, pf_p, pf_c, pf_empty, pf_c2, detA, zrecon
  complex(dp) :: lpf_def(2), lpf_p(2), lpf_c(2), lpf_empty(2), lpf_c2(2), ldetA(2)
  real(dp) :: tol, tol_det

  ! 1. 0x0 matrix
  if (pfaffian(z0)               /= ONE) call fail_msg('pfaffian(0x0 default) /= 1')
  if (pfaffian(z0, alg='P')      /= ONE) call fail_msg('pfaffian(0x0 P) /= 1')
  if (pfaffian(z0, alg='C')      /= ONE) call fail_msg('pfaffian(0x0 C) /= 1')
  if (pfaffian(z0, alg='')       /= ONE) call fail_msg('pfaffian(0x0 empty) /= 1')
  if (array2_pfaffian2(z0)       /= ONE) call fail_msg('array2_pfaffian2(0x0) /= 1')

  lpf_p  = lnPfaffian(z0, alg='P')
  lpf_c  = lnPfaffian(z0, alg='C')
  lpf_c2 = array2_log_pfaffian2(z0)
  if (any(lpf_p  /= [ZERO, ONE])) call fail_msg('lnPfaffian(0x0 P) /= [0,1]')
  if (any(lpf_c  /= [ZERO, ONE])) call fail_msg('lnPfaffian(0x0 C) /= [0,1]')
  if (any(lpf_c2 /= [ZERO, ONE])) call fail_msg('array2_log_pfaffian2(0x0) /= [0,1]')

  ! 2. Non-singular even skew-symmetric matrices
  do k = 1, size(even_sizes)
   n = even_sizes(k)
   allocate(A(n,n))
   call make_random_skew(A)

   detA     = determinant(A)
   pf_def   = pfaffian(A)
   pf_p     = pfaffian(A, alg='P')
   pf_c     = pfaffian(A, alg='C')
   pf_empty = pfaffian(A, alg='')
   pf_c2    = array2_pfaffian2(A)

   tol     = 1.0e-10_dp * max(1.0_dp, abs(pf_p))
   tol_det = 1.0e-10_dp * max(1.0_dp, abs(detA))

   if (abs(pf_def - pf_p) > tol)       call fail_msg('pfaffian default vs P mismatch')
   if (abs(pf_empty - pf_p) > tol)     call fail_msg('pfaffian empty vs P mismatch')
   if (abs(pf_c - pf_p) > tol)         call fail_msg('pfaffian C vs P mismatch')
   if (abs(pf_c2 - pf_p) > tol)        call fail_msg('array2_pfaffian2 vs P mismatch')
   if (abs(pf_p*pf_p - detA) > tol_det) call fail_msg('pfaffian**2 vs det mismatch')

   lpf_def   = lnPfaffian(A)
   lpf_p     = lnPfaffian(A, alg='P')
   lpf_c     = lnPfaffian(A, alg='C')
   lpf_empty = lnPfaffian(A, alg='')
   lpf_c2    = array2_log_pfaffian2(A)

   if (abs(lpf_def(1) - lpf_p(1)) > 1.0e-10_dp .or. abs(lpf_def(2) - lpf_p(2)) > 1.0e-10_dp) then
    call fail_msg('lnPfaffian default vs P mismatch')
   end if
   if (abs(lpf_empty(1) - lpf_p(1)) > 1.0e-10_dp .or. abs(lpf_empty(2) - lpf_p(2)) > 1.0e-10_dp) then
    call fail_msg('lnPfaffian empty vs P mismatch')
   end if
   if (abs(lpf_c(1) - lpf_p(1)) > 1.0e-10_dp .or. abs(lpf_c(2) - lpf_p(2)) > 1.0e-10_dp) then
    call fail_msg('lnPfaffian C vs P mismatch')
   end if
   if (abs(lpf_c2(1) - lpf_p(1)) > 1.0e-10_dp .or. abs(lpf_c2(2) - lpf_p(2)) > 1.0e-10_dp) then
    call fail_msg('array2_log_pfaffian2 vs P mismatch')
   end if

   zrecon = exp(lpf_p(1)) * lpf_p(2)
   if (abs(zrecon - pf_p) > tol) call fail_msg('lnPfaffian(P) recon vs pfaffian mismatch')
   zrecon = exp(lpf_c(1)) * lpf_c(2)
   if (abs(zrecon - pf_c) > tol) call fail_msg('lnPfaffian(C) recon vs pfaffian mismatch')

   deallocate(A)
  end do

  ! 3. Large even N = 200 scaled so |Pf| > 1e308 (tests overflow avoidance in lnPfaffian P & C)
  n = 200
  allocate(A(n,n))
  call make_random_skew(A)
  A = A * 1000.0_dp
  lpf_p = lnPfaffian(A, alg='P')
  lpf_c = lnPfaffian(A, alg='C')
  ldetA = lndet(A)
  if (ieee_is_nan(real(lpf_p(1), dp)) .or. ieee_is_nan(real(lpf_c(1), dp))) then
   call fail_msg('N=200 overflow-avoidance lnPfaffian produced NaN')
  end if
  if (abs(lpf_p(1) - lpf_c(1)) > 1.0e-10_dp * abs(lpf_p(1)) .or. &
      abs(lpf_p(2) - lpf_c(2)) > 1.0e-9_dp) then
   call fail_msg('N=200 lnPfaffian P vs C mismatch')
  end if
  if (abs(2.0_dp * real(lpf_p(1), dp) - real(ldetA(1), dp)) > 1.0e-10_dp * abs(real(ldetA(1), dp))) then
   call fail_msg('N=200 2*Re(lnPfaffian) vs Re(lndet) mismatch')
  end if
  deallocate(A)

  ! 4. Odd-dimensional skew-symmetric matrices (N = 1, 3, 5, 7)
  do k = 1, size(odd_sizes)
   n = odd_sizes(k)
   allocate(A(n,n))
   call make_random_skew(A)
   call check_pf_all_variants_zero('odd skew-symmetric', A)
   deallocate(A)
  end do

  ! 5. Singular even skew-symmetric matrices
  ! 5a. Zero matrices (N = 2, 4, 6)
  do k = 1, 3
   n = 2 * k
   allocate(A(n,n))
   A = ZERO
   call check_pf_all_variants_zero('zero even skew', A)
   deallocate(A)
  end do

  ! 5b. 6x6 skew with row 1 & col 1 zero (zero pivot at i = 1)
  allocate(A(6,6))
  call make_random_skew(A)
  A(1,:) = ZERO
  A(:,1) = ZERO
  call check_pf_all_variants_zero('6x6 zero row1/col1', A)
  deallocate(A)

  ! 5c. 6x6 block-diagonal: nonsing 4x4 skew + zero 2x2 (zero pivot at i = 5)
  allocate(A(6,6), blk4(4,4))
  call make_random_skew(blk4)
  A = ZERO
  A(1:4,1:4) = blk4
  call check_pf_all_variants_zero('6x6 block4 + 0x2', A)

  ! 5d. Exact integer skew block4 moved to rows/cols (1,3,5,6) and zeros on (2,4)
  !     (tests pivot swap + zero pivot at later elimination step with exact integer arithmetic)
  A = ZERO
  A(1,3) = cmplx( 2.0_dp,  1.0_dp, kind=dp); A(3,1) = -A(1,3)
  A(5,6) = cmplx(-1.0_dp,  3.0_dp, kind=dp); A(6,5) = -A(5,6)
  A(1,5) = cmplx( 1.0_dp, -1.0_dp, kind=dp); A(5,1) = -A(1,5)
  A(3,6) = cmplx( 0.0_dp,  2.0_dp, kind=dp); A(6,3) = -A(3,6)
  call check_pf_all_variants_zero('6x6 permuted exact zero block', A)
  deallocate(A, blk4)

  ! 5e. 6x6 block-diagonal: nonsing 2x2 skew + zero 4x4 (zero pivot at i = 3)
  allocate(A(6,6), blk2(2,2))
  call make_random_skew(blk2)
  A = ZERO
  A(1:2,1:2) = blk2
  call check_pf_all_variants_zero('6x6 block2 + 0x4', A)
  deallocate(A, blk2)
 end subroutine test_pfaffian_and_lnpf

end program test_bug05_06
