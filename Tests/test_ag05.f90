!
! Verification test for BUG-07, BUG-12, BUG-13 (print/printna, timer_mod, save/read)
! gfortran -fopenmp -I../lib/gnu/ test_ag05.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t
!
program test_bug07_12_13
  use array_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  complex(dp), parameter :: ZERO = cmplx(0.0_dp, 0.0_dp, kind=dp)
  complex(dp), parameter :: ONE  = cmplx(1.0_dp, 0.0_dp, kind=dp)
  complex(dp), parameter :: IMU  = cmplx(0.0_dp, 1.0_dp, kind=dp)

  integer :: nfail, un, ios, nlines, row_idx
  character(4096) :: bigline
  character(:), allocatable :: long_fmt, comp_fmt
  logical :: ok1, ok2, ok3, ok4, in_gamma

  ! Test arrays for BUG-12
  real(dp) :: ta(40,40), tb(40,40)

  ! Test arrays for BUG-07
  complex(dp), allocatable :: c2(:,:), c1(:), c00(:,:), c30(:,:), c03(:,:), c0(:), cgam(:,:)
  real(dp),    allocatable :: d2(:,:), d1(:), d00(:,:), d30(:,:), d03(:,:), d0(:)
  complex(dp) :: cstatic(2,3), cvstatic(3)
  real(dp)    :: dstatic(2,3), dvstatic(3)

  ! Test arrays for BUG-13
  complex(dp), allocatable :: sz1(:), sz1r(:), sz2(:,:), sz2r(:,:), sz3(:,:,:), sz3r(:,:,:)
  real(dp),    allocatable :: sd1(:), sd1r(:), sd2(:,:), sd2r(:,:), sd3(:,:,:), sd3r(:,:,:)
  complex(dp) :: sz0(0), sz20(0,2), sz30(2,0,3)
  real(dp)    :: sd0(0), sd20(2,0), sd30(0,2,2)

  integer :: i, j, k

  nfail = 0
  long_fmt = '2X,ES24.16E3,1X,1X,1X,1X'
  comp_fmt = '1X,G28.17'

  ! ==========================================================================
  ! 1. BUG-12: timer_mod BEFORE timer_mod_init (must run first!)
  ! ==========================================================================
  call timer_mod_now()
  call check(timer_count == 1, 'BUG-12: pre-init timer_mod_now increments count to 1')
  call check(ieee_is_finite(timer_cputot) .and. ieee_is_finite(timer_cpudel) .and. &
             ieee_is_finite(timer_wtimetot) .and. ieee_is_finite(timer_wtimedel) .and. &
             timer_cputot >= 0.0_dp .and. timer_cputot < 10.0_dp .and. &
             timer_wtimetot >= 0.0_dp .and. timer_wtimetot < 10.0_dp, &
             'BUG-12: pre-init timer_mod_now produces finite bounded times')
  call check(timer_date(1:6) == 'Date: ', 'BUG-12: pre-init timer_mod_now formats timer_date')

  open(newunit=un, file='timer_out.txt', status='replace', action='readwrite')
  call timer_mod_print(un=un, messg='pre-init')
  call check(timer_count == 2, 'BUG-12: pre-init timer_mod_print increments count to 2')

  call timer_mod_init(un=un)
  call check(timer_count == 0 .and. timer_cputot == 0.0_dp .and. timer_cpudel == 0.0_dp .and. &
             timer_wtimetot == 0.0_dp .and. timer_wtimedel == 0.0_dp, &
             'BUG-12: timer_mod_init resets count and all timing vars to 0')

  ta = 1.0_dp
  do i = 1, 100
    tb = matmul(ta, ta)
  end do
  if (tb(1,1) < 0.0_dp) nfail = nfail + 1

  call timer_mod_now()
  call timer_mod_print(un=un, messg='post-init')
  call check(timer_count == 2 .and. timer_cputot >= 0.0_dp .and. timer_wtimetot >= 0.0_dp, &
             'BUG-12: post-init timer_mod_now and timer_mod_print track count and times')

  rewind(un)
  ok1 = .false.; ok2 = .false.; ok3 = .false.
  do
    read(un, '(A)', iostat=ios) bigline
    if (ios /= 0) exit
    if (index(bigline, '#TIMERMOD') > 0 .and. index(bigline, 'pre-init') > 0) ok1 = .true.
    if (index(bigline, '#TIMERINIT::') > 0) ok2 = .true.
    if (index(bigline, '#TIMERMOD') > 0 .and. index(bigline, 'post-init') > 0) ok3 = .true.
  end do
  close(un)
  call check(ok1 .and. ok2 .and. ok3, 'BUG-12: timer_mod output file contains expected lines')

  ! ==========================================================================
  ! 2. BUG-07: print & printna
  ! ==========================================================================
  allocate(c2(-2:0, 3:5), d2(-1:0, 2:4), c1(-2:0), d1(3:5))
  do j = 3, 5
    do i = -2, 0
      c2(i,j) = cmplx(real(i + j, dp), real(i - j, dp), kind=dp)
    end do
  end do
  do j = 2, 4
    do i = -1, 0
      d2(i,j) = real(10*i + j, dp)
    end do
  end do
  c1 = [(cmplx(real(i, dp), -real(i, dp), kind=dp), i = -2, 0)]
  d1 = [(real(i, dp) * 1.5_dp, i = 3, 5)]

  open(newunit=un, file='print_modes.txt', status='replace', action='readwrite')
  ! Absent optionals & custom bounds across modes
  call print(c2, unit=un)
  call print(d2, unit=un)
  call print(c1, unit=un)
  call print(d1, unit=un)
  call print(c2, unit=un, fmt=long_fmt, form='Mathematica', name='C2M')
  call print(d2, unit=un, fmt=long_fmt, form='Mathematica', name='D2M')
  call print(c1, unit=un, fmt=long_fmt, form='Mathematica', name='C1M')
  call print(d1, unit=un, fmt=long_fmt, form='Mathematica', name='D1M')
  call print(c2, unit=un, fmt='', form='', name='EMPTY_FMT')
  call print(d2, unit=un, fmt='', form='', name='EMPTY_FMT_D')
  call print(c1, unit=un, fmt='', form='', name='EMPTY_FMT_V')
  call print(d1, unit=un, fmt='', form='', name='EMPTY_FMT_DV')
  call print(c2, unit=un, form='Info', name='INFO_C2')
  call print(d2, unit=un, form='Info', name='INFO_D2')

  ! Non-allocatable printna checks with long format
  cstatic = reshape([(cmplx(real(i, dp), 0.25_dp, kind=dp), i = 1, 6)], [2,3])
  dstatic = reshape([(real(i, dp), i = 1, 6)], [2,3])
  cvstatic = [(1.0_dp, -1.0_dp), (2.0_dp, 0.5_dp), (3.0_dp, 0.0_dp)]
  dvstatic = [1.0_dp, 2.0_dp, 3.0_dp]
  call printna(cstatic, unit=un, fmt=long_fmt, form='Default', name='CS')
  call printna(dstatic, unit=un, fmt=long_fmt, form='Mathematica', name='DS')
  call printna(cvstatic, unit=un, fmt=long_fmt, form='Default', name='CVS')
  call printna(dvstatic, unit=un, fmt=long_fmt, form='Mathematica', name='DVS')

  ! Reversed / out-of-range slice bounds (previously read OOB in Mathematica j=js and Default)
  call print(c2, unit=un, form='Default',     ips=5)
  call print(c2, unit=un, form='Mathematica', ips=5)
  call print(c2, unit=un, form='Mathematica', jps=10)
  call print(c2, unit=un, form='Mathematica', jpe=0)
  call print(c2, unit=un, form='Mathematica', ips=-1, ipe=-2, jps=5, jpe=3)
  call print(c2, unit=un, form='Gamma',       jps=5, jpe=3)
  call print(d2, unit=un, form='Mathematica', jps=10)
  call print(d2, unit=un, form='Mathematica', jpe=0)
  call print(d2, unit=un, form='Default',     jps=5, jpe=3)
  call print(c1, unit=un, form='Default',     ips=10)
  call print(c1, unit=un, form='Mathematica', ips=10)
  call print(d1, unit=un, form='Default',     ipe=-5)
  call print(d1, unit=un, form='Mathematica', ipe=-5)

  ! 0-sized allocated arrays
  allocate(c00(0,0), c30(3,0), c03(0,3), c0(0))
  allocate(d00(0,0), d30(3,0), d03(0,3), d0(0))
  call print(c00, unit=un, form='Default', name='C00')
  call print(c00, unit=un, form='Mathematica', name='C00M')
  call print(c30, unit=un, form='Mathematica', name='C30M')
  call print(c03, unit=un, form='Default', name='C03D')
  call print(c0,  unit=un, form='Default', name='C0D')
  call print(c0,  unit=un, form='Mathematica', name='C0M')
  call print(d00, unit=un, form='Default', name='D00')
  call print(d30, unit=un, form='Mathematica', name='D30M')
  call print(d03, unit=un, form='Mathematica', name='D03M')
  call print(d0,  unit=un, form='Mathematica', name='D0M')
  close(un)
  call check(.true., 'BUG-07: print/printna handle slices, long/empty fmt, 0-size without crashing')

  ! Specific test for Gamma mode: unit respected & no 2000-char row truncation
  allocate(cgam(2, 700))
  cgam = ZERO
  cgam(:,1) = ZERO; cgam(:,2) = ONE; cgam(:,3) = -ONE; cgam(:,4) = IMU; cgam(:,5) = -IMU
  cgam(1,700) = cmplx(2.0_dp, 3.0_dp, kind=dp)
  cgam(2,700) = -ONE
  open(newunit=un, file='print_gamma.txt', status='replace', action='readwrite')
  call print(cgam, unit=un, form='Gamma', name='BIGGAM')
  rewind(un)
  ok1 = .false.; ok2 = .false.; ok3 = .false.; ok4 = .false.
  in_gamma = .false.; row_idx = 0
  do
    read(un, '(A)', iostat=ios) bigline
    if (ios /= 0) exit
    if (index(bigline, '#pmstart:') > 0) then
      in_gamma = .true.
      ok1 = .true.
    else if (index(bigline, '#pmend:') > 0) then
      in_gamma = .false.
      ok2 = .true.
    else if (in_gamma .and. index(bigline, '#---') == 0) then
      row_idx = row_idx + 1
      if (row_idx == 1) then
        ok3 = (len_trim(bigline) == 2100 .and. bigline(1:15) == '  0  1 -1  i -i' .and. &
               bigline(2098:2100) == ' NN')
      else if (row_idx == 2) then
        ok4 = (len_trim(bigline) == 2100 .and. bigline(1:15) == '  0  1 -1  i -i' .and. &
               bigline(2098:2100) == ' -1')
      end if
    end if
  end do
  close(un)
  call check(ok1 .and. ok2 .and. row_idx == 2 .and. ok3 .and. ok4, &
             'BUG-07: Gamma mode writes untruncated 2100-char rows to unit')

  ! ==========================================================================
  ! 3. BUG-13: save & read round-trips
  ! ==========================================================================
  allocate(sz1(4), sz1r(4), sz2(3,4), sz2r(3,4), sz3(2,3,2), sz3r(2,3,2))
  allocate(sd1(4), sd1r(4), sd2(3,4), sd2r(3,4), sd3(2,3,2), sd3r(2,3,2))
  do i = 1, 4
    sz1(i) = cmplx(real(i, dp) * 0.375_dp, -real(i, dp) * 1.125_dp, kind=dp)
    sd1(i) = real(i, dp) * 2.75_dp - 5.0_dp
  end do
  do j = 1, 4
    do i = 1, 3
      sz2(i,j) = cmplx(real(10*i + j, dp), real(i - 2*j, dp), kind=dp) * 1.0e-3_dp
      sd2(i,j) = real(10*i - j, dp) * 0.125_dp
    end do
  end do
  do k = 1, 2
    do j = 1, 3
      do i = 1, 2
        sz3(i,j,k) = cmplx(real(i + j*k, dp), real(k - i, dp), kind=dp)
        sd3(i,j,k) = real(i*3 + j + k*5, dp) * 0.25_dp
      end do
    end do
  end do

  ! 3a. Composite fmt='1X,G28.17' (verifies line count == size and exact round-trip)
  open(newunit=un, file='save_comp_c.txt', status='replace', action='readwrite')
  call save(sz1, unit=un, fmt=comp_fmt)
  call save(sz2, unit=un, fmt=comp_fmt)
  call save(sz3, unit=un, fmt=comp_fmt)
  rewind(un)
  nlines = 0
  do
    read(un, '(A)', iostat=ios) bigline
    if (ios /= 0) exit
    nlines = nlines + 1
  end do
  call check(nlines == size(sz1) + size(sz2) + size(sz3), &
             'BUG-13: save(complex, fmt="1X,G28.17") writes 1 line per element')
  rewind(un)
  call read(sz1r, unit=un); call read(sz2r, unit=un); call read(sz3r, unit=un)
  close(un)
  call check(maxval(abs(sz1 - sz1r)) < 1.0e-13_dp .and. &
             maxval(abs(sz2 - sz2r)) < 1.0e-13_dp .and. &
             maxval(abs(sz3 - sz3r)) < 1.0e-13_dp, &
             'BUG-13: save/read complex composite fmt round-trip')

  open(newunit=un, file='save_comp_d.txt', status='replace', action='readwrite')
  call save(sd1, unit=un, fmt=comp_fmt)
  call save(sd2, unit=un, fmt=comp_fmt)
  call save(sd3, unit=un, fmt=comp_fmt)
  rewind(un)
  call read(sd1r, unit=un); call read(sd2r, unit=un); call read(sd3r, unit=un)
  close(un)
  call check(maxval(abs(sd1 - sd1r)) < 1.0e-13_dp .and. &
             maxval(abs(sd2 - sd2r)) < 1.0e-13_dp .and. &
             maxval(abs(sd3 - sd3r)) < 1.0e-13_dp, &
             'BUG-13: save/read real composite fmt round-trip')

  ! 3b. Long fmt > 20 chars
  open(newunit=un, file='save_long.txt', status='replace', action='readwrite')
  call save(sz1, unit=un, fmt=long_fmt)
  call save(sz2, unit=un, fmt=long_fmt)
  call save(sz3, unit=un, fmt=long_fmt)
  call save(sd1, unit=un, fmt=long_fmt)
  call save(sd2, unit=un, fmt=long_fmt)
  call save(sd3, unit=un, fmt=long_fmt)
  rewind(un)
  sz1r = ZERO; sz2r = ZERO; sz3r = ZERO
  sd1r = 0.0_dp; sd2r = 0.0_dp; sd3r = 0.0_dp
  call read(sz1r, unit=un); call read(sz2r, unit=un); call read(sz3r, unit=un)
  call read(sd1r, unit=un); call read(sd2r, unit=un); call read(sd3r, unit=un)
  close(un)
  call check(maxval(abs(sz1 - sz1r)) < 1.0e-13_dp .and. &
             maxval(abs(sz2 - sz2r)) < 1.0e-13_dp .and. &
             maxval(abs(sz3 - sz3r)) < 1.0e-13_dp .and. &
             maxval(abs(sd1 - sd1r)) < 1.0e-13_dp .and. &
             maxval(abs(sd2 - sd2r)) < 1.0e-13_dp .and. &
             maxval(abs(sd3 - sd3r)) < 1.0e-13_dp, &
             'BUG-13: save/read long (>20 chars) fmt round-trip')

  ! 3c. Blank fmt='' & 0-sized save
  open(newunit=un, file='save_blank.txt', status='replace', action='readwrite')
  call save(sz1, unit=un, fmt='')
  call save(sd2, unit=un, fmt='')
  call save(sz0, unit=un); call save(sz20, unit=un); call save(sz30, unit=un)
  call save(sd0, unit=un); call save(sd20, unit=un); call save(sd30, unit=un)
  rewind(un)
  sz1r = ZERO; sd2r = 0.0_dp
  call read(sz1r, unit=un); call read(sd2r, unit=un)
  close(un)
  call check(maxval(abs(sz1 - sz1r)) < 1.0e-13_dp .and. &
             maxval(abs(sd2 - sd2r)) < 1.0e-13_dp, &
             'BUG-13: save/read blank fmt="" fallback and 0-sized save')

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
end program test_bug07_12_13
