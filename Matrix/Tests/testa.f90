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
 integer                :: i,j
 real   (dp)            :: r,r1,r2,r3,r4,r5,r6,r7,r8,r9,sigma, ldetd(2),lPfd(2)
 complex(dp)            :: z,z1,z2,z3,z4,z5,z6,z7,z8,z9,       ldet (2),lPf (2)
 integer                :: m,n,is,ie,js,je
 real   (8),allocatable :: d1(:,:  ),d2(:,:  ),d3(:,:  ),d4(:,:  ),d5(:,:),d6(:,:),d7(:,:),d8(:,:),d9(:,:),d10(:,:)
 complex(8),allocatable :: m1(:,:  ),m2(:,:  ),m3(:,:  ),m4(:,:  ),m5(:,:),m6(:,:),m7(:,:),m8(:,:),m9(:,:),m10(:,:)
 real   (8),allocatable :: u1(:    ),u2(:    ),u3(:    ),u4(:    ),u5(:  ),u6(:  ),u7(:  ),u8(:  ),u9(:  ),u10(:  )
 complex(8),allocatable :: v1(:    ),v2(:    ),v3(:    ),v4(:    ),v5(:  ),v6(:  ),v7(:  ),v8(:  ),v9(:  ),v10(:  ),w(:,:)
 character(1000)        :: string
 character(1)           :: form
!-----------------------------------------------------------------------------------------------------------------------------------
 call matrix_random_init
 f_mout = stdout        ! module prints messages to stdout
!-----------------------------------------------------------------------------------------------------------------------------------
!m=6;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!allocate(m1(is:is+m-1,js:js+n-1)); call random_number(m1,sigma=PI); call print(m1,stdout,fmt='F12.3',form=form,name='m1')
!-----------------------------------------------------------------------------------------------------------------------------------
m=48;n=m+8;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
allocate(m1(m,n),m2(n,m),m3(m,m),m4(n,n),v1(m),v2(n),v3(m),v4(n))
allocate(d1(m,n),d2(n,m),d3(m,m),d4(n,n),u1(m),u2(n),u3(m),u4(n))
allocate(m5(n,n),m6(n,n),m7(n,n),d5(n,n),d6(n,n),d7(n,n))
call random_number(m1,sigma=PI); call random_number(m2,sigma=PI); 
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
