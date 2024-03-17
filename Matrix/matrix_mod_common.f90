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
  integer                                         :: un,nseeds
  integer, allocatable                            :: seeds(:)
  
  write(f_mout,'(A)')'# matrix_random_init: Initializing random_number from /dev/urandom'
  
  call random_seed(size=nseeds)
  ALLOCATE(  seeds(nseeds))
  open(newunit=un,file="/dev/urandom", access="stream", form="unformatted")
  read(un)   seeds
  seeds=ABS (seeds)
  call random_seed(put = seeds)
  DEALLOCATE(seeds) 
! call random_init(.false.,.false.)    ! fortran2018 function, not available with gfortran

 end subroutine     matrix_random_init
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_error                                           (message)                                                !:.:.:
  character(*) message
  write(stderr,'(A)')                                                                                 NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
       '# ERROR: matrix module error: '// TRIM(message)                                             //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//NEW_LINE('A')
  stop
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
