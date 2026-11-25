!
! gfortran -I../lib/gnu/ test_v15.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t > hi
!
! gnuplot:
! s    = 1  # sigma: make sure this is correct
! g(x) = sqrt(1.0/(2.0*pi*s**2)) * exp(- x**2 / (2.0 * s**2))
!
! s=1;plot "<awk '/^R1/{print $2}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R1/{print $3}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R1/{print $4}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R2/{print $2}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R2/{print $3}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R2/{print $4}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R3/{print $2}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R3/{print $3}' hi| hist -b 0.1" u 1:4 w p,g(x)
! s=1;plot "<awk '/^R3/{print $4}' hi| hist -b 0.1" u 1:4 w p,g(x)
!
! set log y  # and repeat the above to check tails
! s=1;plot for [n in "2 3 4"] "<awk '/^R1/{print $".n."}' hi| hist -b 0.1 -d 10000000" u 1:4 w p t n,g(x)
! s=1;plot for [n in "2 3 4"] "<awk '/^R2/{print $".n."}' hi| hist -b 0.1 -d 10000000" u 1:4 w p t n,g(x)
! s=1;plot for [n in "2 3 4"] "<awk '/^R3/{print $".n."}' hi| hist -b 0.1 -d 10000000" u 1:4 w p t n,g(x)
!
! Check accuracy: set log y
! s=1;plot for [n in "2 3 4"] "<awk '/^R3/{print $".n."}' hi| hist -b 0.1 -d 10000000" u 1:(abs($4-g($1))) w p t n

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
 integer                :: i,j,k,imeas,nmeas,m,n
 real   (dp)            :: s
 real   (dp),allocatable:: d1(:),d2(:,:),d3(:,:,:)
 complex(dp),allocatable:: m1(:),m2(:,:),m3(:,:,:)
 character(1000)        :: string
 character(1)           :: form
!-----------------------------------------------------------------------------------------------------------------------------------
 call matrix_random_init
 f_mout = stdout        ! module prints messages to stdout
!-----------------------------------------------------------------------------------------------------------------------------------
!m=6;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!allocate(m1(is:is+m-1,js:js+n-1)); call random_number(m1,sigma=PI); call print(m1,stdout,fmt='F12.3',form=form,name='m1')
!-----------------------------------------------------------------------------------------------------------------------------------
m=10;n=m+8
allocate(m1(m),m2(m,n),m3(m,m,n))
allocate(d1(m),d2(m,n),d3(m,m,n))

s = 1.0_dp

nmeas    = 100000
do imeas = 1, nmeas
 call random_number(m1,sigma=s); call random_number(d1,sigma=s); 
 do i = 1,size(m1,1)
  print '(A,20G28.16)','R1 ',d1(i),m1(i)
 end do
end do

nmeas    = 10000
do imeas = 1, nmeas
 call random_number(m2,sigma=s); call random_number(d2,sigma=s); 
 do j = 1,size(m2,2); do i = 1,size(m2,1)
  print '(A,20G28.16)','R2 ',d2(i,j),m2(i,j)
 end do             ; end do
end do

nmeas    = 1000
do imeas = 1, nmeas
 call random_number(m3,sigma=s); call random_number(d3,sigma=s); 
 do k = 1,size(m3,3); do j = 1,size(m3,2); do i = 1,size(m3,1)
  print '(A,20G28.16)','R3 ',d3(i,j,k),m3(i,j,k)
 end do             ; end do             ; end do
end do


!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
