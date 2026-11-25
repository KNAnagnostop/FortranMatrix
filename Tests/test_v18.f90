!
! gfortran -I../lib/gnu/ test_v17.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t
!
program             testme
 use, intrinsic         :: iso_fortran_env, stdin => input_unit, stdout => output_unit, stderr => error_unit 
 use, intrinsic         :: ieee_arithmetic  
 use                    :: array_mod 
 implicit none
!-----------------------------------------------------------------------------------------------------------------------------------
 integer   ,parameter   :: dp  = real64
!-----------------------------------------------------------------------------------------------------------------------------------
 real   (dp), parameter :: PI  = atan2(0.0_dp,-1.0_dp        ), TWOPI = 2.0_dp * PI
 complex(dp), parameter :: ONE = CMPLX(1.0_dp, 0.0_dp,KIND=dp), ZERO  = CMPLX(0.0_dp,0.0_dp,KIND=dp)
 complex(dp), parameter :: IMU = CMPLX(0.0_dp, 1.0_dp,KIND=dp)
!-----------------------------------------------------------------------------------------------------------------------------------
 integer   ,parameter   :: f_out = 13, f_math = 14
 integer                :: i,j,k,m,n,l,itest
 real   (dp),allocatable:: d1(:,:),d2(:,:),d3(:,:),d4(:,:),v1(:),v2(:)
 complex(dp),allocatable:: m1(:,:),m2(:,:),m3(:,:),m4(:,:),w1(:),w2(:)
 character(1000)        :: string
 character(1)           :: form
!-----------------------------------------------------------------------------------------------------------------------------------
 complex(dp)            :: trz,trz1
 real   (dp)            :: trd,trd1
 real   (dp)            :: tol,err
!-----------------------------------------------------------------------------------------------------------------------------------
 call matrix_random_init
 f_mout = stdout        ! module prints messages to stdout
!-----------------------------------------------------------------------------------------------------------------------------------
!m=6;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!allocate(m1(is:is+m-1,js:js+n-1)); call random_number(m1,sigma=PI); call print(m1,stdout,fmt='F12.3',form=form,name='m1')
!-----------------------------------------------------------------------------------------------------------------------------------
tol=-1.0e-15_dp
m=10;l=5;n=12
allocate(m1(m,l),m2(l,n),m3(m,n),m4(m,n),w1(n))
allocate(d1(m,l),d2(l,n),d3(m,n),d4(m,n),v1(n))

!-----------------------------------------------------------------------------------------------------------------------------------
!Matrix-matrix multiplication:
call random_number(m1);call random_number(m2);call random_number(m3);call random_number(w1);
call random_number(d1);call random_number(d2);call random_number(d3);call random_number(v1);

itest=itest+1;err=norm(lmatmul(m1,m2)-matmul(m1,m2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltz: ',err,norm(m1)

!-----------------------------------------------------------------------------------------------------------------------------------
10 continue
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
