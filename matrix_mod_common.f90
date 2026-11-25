!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_common.f90 ............................................................................!:.:.:
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
module              matrix_mod_common                                                                                               !:.:.:
 use, intrinsic ::  iso_fortran_env, only: real64, compiler_version, compiler_options, &
                                           stdin => input_unit, stdout => output_unit, stderr => error_unit 
 use, intrinsic ::  ieee_arithmetic   
 implicit none
 save
!-----------------------------------------------------------------------------------------------------------------------------------
 integer    , parameter                           :: dp  = real64
!-----------------------------------------------------------------------------------------------------------------------------------
 real   (dp), parameter                           :: PI  = atan2(0.0_dp,-1.0_dp        ), TWOPI = 2.0_dp * PI
 complex(dp), parameter                           :: ONE = CMPLX(1.0_dp, 0.0_dp,KIND=dp), ZERO  = CMPLX(0.0_dp,0.0_dp,KIND=dp)
 complex(dp), parameter                           :: IMU = CMPLX(0.0_dp, 1.0_dp,KIND=dp)
!-----------------------------------------------------------------------------------------------------------------------------------
 integer    , parameter                           :: mtype_len = 4, name_len = 20                                                   ! size of character variables MAT % mtype, MAT % name
!-----------------------------------------------------------------------------------------------------------------------------------
 integer                                          :: f_mout    = stdout, f_minput = stdin                                           ! Set this to your unit of preference for module output/input
!-----------------------------------------------------------------------------------------------------------------------------------
! Used for printing fortran numbers to Mathematica numbers: print *,Mathematica_nleft, r, '"]' must be a valid Mathematica number
! Try also: nleft='Internal`StringToDouble['; nright=']'
 integer                , parameter               :: math_len=25
 character(len=math_len), parameter               :: Mathematica_nleft ='(Interpreter["Number"]["'
 character(len=math_len), parameter               :: Mathematica_nright='"])'

!-----------------------------------------------------------------------------------------------------------------------------------
 interface          NaN                                                                                                             !:.:.:
  module procedure                                :: matrix_return_nan_complex8
  module procedure                                :: matrix_return_nan_complex8_z   , matrix_return_nan_real8_d
 end interface      NaN
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isNaN                                                                                                           !:.:.:
  module procedure                                :: matrix_return_is_nan_complex8_z, matrix_return_is_nan_real8_d
 end interface      isNaN
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_random_init                                                                                              !:.:.:
  integer                                         :: un,nseeds,ios
  integer, allocatable                            :: seeds(:)
  
  write(f_mout,'(A)')'# matrix_random_init: Initializing random_number from /dev/urandom'
  
  call random_seed(size=nseeds)
  ALLOCATE(  seeds(nseeds))
  ! 261124 (BUG-06): do not crash when /dev/urandom is missing or unreadable (Windows, sandboxes, containers)
  open(newunit=un,file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=ios)
  if(ios == 0)then
   read(un, iostat=ios) seeds
   close(un)
  end if
  if(ios /= 0)then
   write(f_mout,'(A)')'# matrix_random_init: WARNING: cannot read /dev/urandom, using call random_seed() instead'
   call random_seed()                  ! processor-dependent seeding (non-repeatable with gfortran >= 7)
   DEALLOCATE(seeds)
   return
  end if
  seeds    = iand(seeds, huge(seeds))  ! clear the sign bit (no overflow, unlike ABS(-huge-1))
  seeds(1) = ior (seeds(1), 1)         ! guarantee a seed that is not all zeros
  call random_seed(put = seeds)
  DEALLOCATE(seeds) 
! The Fortran 2018 call random_init(.false.,.false.) is an alternative (gfortran >= 9, ifx); /dev/urandom is kept for portability
! to older compilers and for full-entropy seeding.

 end subroutine     matrix_random_init
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_error                                           (message)                                                !:.:.:
  character(*), intent(in)                        :: message
  write(stderr,'(A)')                                                                                 NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
       '# ERROR: matrix module error: '// TRIM(message)                                             //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//NEW_LINE('A')
  flush(stderr)
  error stop 1                                                                                     ! 261124: non-zero exit status (was: stop)
 end subroutine     matrix_error
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_return_nan_complex8                             ()                       result(r)                       !:.:.:
  complex(dp)                                     :: r
  real   (dp)                                     :: x

  x = IEEE_VALUE (0.0_dp,IEEE_QUIET_NAN)
  r = CMPLX(x,x,kind=dp)

 end function       matrix_return_nan_complex8
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_return_nan_complex8_z                           (z)                      result(r)                       !:.:.:
  complex(dp), intent(in)                         :: z
  complex(dp)                                     :: r
  real   (dp)                                     :: x

  x = IEEE_VALUE (REAL(z,kind=dp),IEEE_QUIET_NAN)
  r = CMPLX(x,x,kind=dp)

 end function       matrix_return_nan_complex8_z
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_return_nan_real8_d                              (z)                      result(r)                       !:.:.:
  real   (dp), intent(in)                         :: z
  real   (dp)                                     :: r

  r = IEEE_VALUE (z,IEEE_QUIET_NAN)

 end function       matrix_return_nan_real8_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_return_is_nan_complex8_z                        (z)                      result(itis)                    !:.:.:
  complex(dp), intent(in)                         :: z
  logical                                         :: itis
  logical                                         :: isx, isy

  isx = IEEE_IS_NAN (REAL (z,kind=dp))
  isy = IEEE_IS_NAN (AIMAG(z        ))

  if( isx .or. isy )then
   itis = .true.
  else
   itis = .false.
  end if
  
 end function       matrix_return_is_nan_complex8_z
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_return_is_nan_real8_d                           (z)                      result(itis)                    !:.:.:
  real(dp),    intent(in)                         :: z
  logical                                         :: itis

  itis = IEEE_IS_NAN (z)
  
 end function       matrix_return_is_nan_real8_d


!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module          matrix_mod_common
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================



!-----------------------------------------------------------------------------------------------------------------------------------
!  Copyright by Konstantinos N. Anagnostopoulos, Physics Department, National Technical University of Athens, 2020
!  konstant@mail.ntua.gr, www.physics.ntua.gr/konstant
!  
!  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
!  published by the Free Software Foundation, version 3 of the License.
!  
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public Liense along with this program. If not, see http://www.gnu.org/licenses
!-----------------------------------------------------------------------------------------------------------------------------------
