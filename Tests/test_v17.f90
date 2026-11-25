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
itest=itest+1;err=norm(lmatmul(d1,d2)-matmul(d1,d2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltd: ',err,norm(d1)

itest=itest+1;call mmmult(m1,m2,m3);err=norm(m3-matmul(m1,m2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltz: ',err,norm(m1)
itest=itest+1;call mmmult(d1,d2,d3);err=norm(d3-matmul(d1,d2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltd: ',err,norm(d1)

!-----------------------------------------------------------------------------------------------------------------------------------
!Square matrices:
deallocate(m1,m2,m3,m4,d1,d2,d3,d4,v1,w1);
allocate(m1(n,n),m2(n,n),m3(n,n),m4(n,n),w1(n),w2(n))
allocate(d1(n,n),d2(n,n),d3(n,n),d4(n,n),v1(n),v2(n))
call random_number(m1);call random_number(m2);call random_number(m3);call random_number(w1);
call random_number(d1);call random_number(d2);call random_number(d3);call random_number(v1);
!-----------------------------------------------------------------------------------------------------------------------------------

itest=itest+1;err=norm(lmatmul(m1)         -matmul(m1,hermitian(m1)))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmmdz: ',err,norm(m1)
itest=itest+1;err=norm(lmatmul(m1,side='H')-matmul(hermitian(m1),m1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmmdz: ',err,norm(m1)

itest=itest+1;err=norm(lmatmul(d1)         -matmul(d1,symmetric(d1)))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmmdd: ',err,norm(d1)
itest=itest+1;err=norm(lmatmul(d1,side='S')-matmul(symmetric(d1),d1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmmdd: ',err,norm(d1)

itest=itest+1;call mmmult(m1,m4);err=norm(m4-matmul(m1,hermitian(m1)))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltz: ',err,norm(m1)
itest=itest+1;call mmmult(d1,d4);err=norm(d4-matmul(d1,symmetric(d1)))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltd: ',err,norm(d1)


itest=itest+1;call mmmult(m1,m4,side='H');err=norm(m4-matmul(hermitian(m1),m1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltz: ',err,norm(m1)
itest=itest+1;call mmmult(d1,d4,side='S');err=norm(d4-matmul(symmetric(d1),d1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltd: ',err,norm(d1)

!-----------------------------------------------------------------------------------------------------------------------------------
!Symmetric matrix multiplication:
call hermitian_set(m1);call symmetric_set(d1)

itest=itest+1;err=norm(lmatmul(m1,m2,mtype='ZH')-matmul(m1,m2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmzhz: ',err,norm(m1)
itest=itest+1;err=norm(lmatmul(d1,d2,mtype='DS')-matmul(d1,d2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmdsd: ',err,norm(d1)

itest=itest+1;call mmmult(m1,m2,m3,mtype='ZH');err=norm(m3-matmul(m1,m2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltz: ',err,norm(m1)
itest=itest+1;call mmmult(d1,d2,d3,mtype='DS');err=norm(d3-matmul(d1,d2))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mmltd: ',err,norm(d1)

itest=itest+1;err=norm(lmatmul(m1,w1,type='H')-matmul(m1,w1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mvHEz: ',err,norm(m1)

itest=itest+1;call mvmult(m1,w1,w2,type='H');err=norm(w2-matmul(m1,w1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mvhez: ',err,norm(m1)

!-----------------------------------------------------------------------------------------------------------------------------------
!Matrix-vector multiplication: non-square matrices
deallocate(m1,m2,m3,m4,d1,d2,d3,d4,v1,w1,v2,w2); 
allocate(m1(m,l),m2(l,n),m3(m,n),m4(n,m),w1(n),w2(m))
allocate(d1(m,l),d2(l,n),d3(m,n),d4(n,m),v1(n),v2(m))
call random_number(m1);call random_number(m2);call random_number(m3);call random_number(w1);
call random_number(d1);call random_number(d2);call random_number(d3);call random_number(v1);
!-----------------------------------------------------------------------------------------------------------------------------------

itest=itest+1;err=norm(lmatmul(  m3,  w1         )-matmul(m3,w1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mv31z: ',err,norm(m3)
itest=itest+1;err=norm(lmatmul(v=w1,A=m4,type='T')-matmul(w1,m4))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mvT3z: ',err,norm(m4)
itest=itest+1;err=norm(lmatmul(v=w1,A=m4,type='C')-matmul(w1,conjg(m4)))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mvH3z: ',err,norm(m4)

itest=itest+1;err=norm(lmatmul(  d3,  v1         )-matmul(d3,v1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mv31d: ',err,norm(d3)
itest=itest+1;err=norm(lmatmul(v=v1,A=d4,type='T')-matmul(v1,d4))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mvT3d: ',err,norm(d4)

itest=itest+1;err=norm(lmatmul(v=w1,A=m4,type='T')-matmul(w1,m4))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' vm14z: ',err,norm(m4)
itest=itest+1;err=norm(lmatmul(v=v1,A=d4,type='T')-matmul(v1,d4))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' vm14d: ',err,norm(d4)

itest=itest+1;call mvmult(m3,w1,w2);err=norm(w2-matmul(m3,w1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mv31Z: ',err,norm(m3)
itest=itest+1;call mvmult(d3,v1,v2);err=norm(v2-matmul(d3,v1))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mv31D: ',err,norm(d3)

itest=itest+1;call mvmult(m3,w2,w1,type='T');err=norm(w1-matmul(w2,m3))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mv23Z: ',err,norm(m3)
itest=itest+1;call mvmult(m3,w2,w1,type='C');err=norm(w1-matmul(w2,conjg(m3)))
if(err>tol) print '(I0.3,A,1000G16.3)',itest,' mv23C: ',err,norm(m3)


!-----------------------------------------------------------------------------------------------------------------------------------
10 continue
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
