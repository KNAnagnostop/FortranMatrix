!
! gfortran -I../lib/gnu/ test_v16.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t,
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
 integer                :: i,j,k,m,n,itest
 real   (dp)            :: s
 real   (dp),allocatable:: d1(:),d2(:,:),d3(:,:,:)
 complex(dp),allocatable:: m1(:),m2(:,:),m3(:,:,:)
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
m=10;n=m;itest=0;tol=-1.0e-15_dp
allocate(m1(m),m2(m,n),m3(m,m,n))
allocate(d1(m),d2(m,n),d3(m,m,n))

call random_number(m1);call random_number(m2);call random_number(m3);
call random_number(d1);call random_number(d2);call random_number(d3);

itest=itest+1;trz=trace2c(m2);call traceless_set(m2);trz1=trace2(m2);err=abs(trz-trz1);call random_number(m2);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' trc2z: ',err,abs(trz)
itest=itest+1;trd=trace2c(d2);call traceless_set(d2);trd1=trace2(d2);err=abs(trd-trd1);call random_number(d2);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' trc2d: ',err,abs(trd)

itest=itest+1;trz=trace2c(m2);err=abs(trace2c(m2)-trace2(traceless(m2)));call random_number(m2);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' trc2z: ',err,abs(trz)
itest=itest+1;trd=trace2c(d2);err=abs(trace2c(d2)-trace2(traceless(d2)));call random_number(d2);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' trc2d: ',err,abs(trd)

itest=itest+1;do i=1,n;m1(i)=m2(i,i);end do;err=norm(diagonal(m2)-m1);call random_number(m1);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagz: ',err,norm(m1)
itest=itest+1;do i=1,n;d1(i)=d2(i,i);end do;err=norm(diagonal(d2)-d1);call random_number(d1);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagd: ',err,norm(d1)

itest=itest+1;m2=diagonal(m1);err=norm(m1-diagonal(m2));call random_number(m2);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagz: ',err,norm(m1)
itest=itest+1;d2=diagonal(d1);err=norm(d1-diagonal(d2));call random_number(d2);
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagd: ',err,norm(d1)

itest=itest+1;m1=3;m2=diagonal(3.0_dp,n);err=norm(m2-diagonal(m1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagz: ',err,norm(m1);call random_number(m1);call random_number(m2);
itest=itest+1;d1=3;d2=diagonal(3.0_dp,n);err=norm(d2-diagonal(d1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagd: ',err,norm(d1);call random_number(d1);call random_number(d2);

itest=itest+1;call random_number(trz);m1=trz;m2=diagonal(trz,n);err=norm(m2-diagonal(m1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagz: ',err,norm(m1);call random_number(m1);call random_number(m2);
itest=itest+1;call random_number(trz);d1=trz;d2=diagonal(trz,n);err=norm(d2-diagonal(d1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' diagz: ',err,norm(d1);call random_number(d1);call random_number(d2);

itest=itest+1;err=norm(hermitian(m2)-conjg(transpose(m2)))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' hermz: ',err,norm(m2);
itest=itest+1;err=norm(symmetric(m2)-transpose(m2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' symmz: ',err,norm(m2);
itest=itest+1;err=norm(symmetric(d2)-transpose(d2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' symmd: ',err,norm(d2);

itest=itest+1;err=norm(antisymmetric(m2)+transpose(m2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' asymz: ',err,norm(m2);
itest=itest+1;err=norm(antisymmetric(d2)+transpose(d2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' asymd: ',err,norm(d2);

!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
