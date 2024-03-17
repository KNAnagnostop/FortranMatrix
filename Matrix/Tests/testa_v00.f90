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
call random_number(m5,sigma=PI); call random_number(m6,sigma=PI); call hermitian_set(m5); call hermitian_set(m6)
call random_number(v1,sigma=PI); call random_number(v2,sigma=PI);
call random_number(d1,sigma=PI); call random_number(d2,sigma=PI); 
call random_number(d5,sigma=PI); call random_number(d6,sigma=PI); call symmetric_set(d5); call symmetric_set(d6)
call random_number(u1,sigma=PI); call random_number(u2,sigma=PI);
!-----------------------------------------------------------------------------------------------------------------------------------
!call print(m1,stdout,fmt='F12.3',form=form,name='m1')
call mmmult(m1,m2,m3)
call mmmult(m2,m1,m4)
print *,'m1*m2: ', norm( m3 - MATMUL(m1, m2)) , norm( m3 - (m1 .mm. m2)) 
print *,'m2*m1: ', norm( m4 - MATMUL(m2, m1)) 
!-----------------------------------------------------------------------------------------------------------------------------------
call mvmult(m1,v2,v3)
call mvmult(m2,v1,v4)
print *,'m1*v2: ', norm(v3-MATMUL(m1,v2)),norm(v3 - (m1 .mm. v2) )
print *,'m2*v1: ', norm(v4-MATMUL(m2,v1)),norm(v4 - (m2 .mm. v1) )
!-----------------------------------------------------------------------------------------------------------------------------------
call vmmult(v1,m1,v4)
call vmmult(v2,m2,v3)
print *,'v1*m1: ', norm(v4-MATMUL(v1,m1)),norm(v4 - (v1 .mm. m1) )
print *,'v2*m2: ', norm(v3-MATMUL(v2,m2)),norm(v3 - (v2 .mm. m2) )
!-----------------------------------------------------------------------------------------------------------------------------------
call mmmult(d1,d2,d3)
call mmmult(d2,d1,d4)
print *,'d1*d2: ', norm( d3 - MATMUL(d1, d2)) , norm( d3 - (d1 .mm. d2)) 
print *,'d2*d1: ', norm( d4 - MATMUL(d2, d1)) 
!-----------------------------------------------------------------------------------------------------------------------------------
call mvmult(d1,u2,u3)
call mvmult(d2,u1,u4)
print *,'d1*u2: ', norm(u3-MATMUL(d1,u2)),norm(u3 - (d1 .mm. u2) )
print *,'d2*u1: ', norm(u4-MATMUL(d2,u1)),norm(u4 - (d2 .mm. u1) )
!-----------------------------------------------------------------------------------------------------------------------------------
call vmmult(u1,d1,u4)
call vmmult(u2,d2,u3)
print *,'u1*d1: ', norm(u4-MATMUL(u1,d1)),norm(u4 - (u1 .mm. d1) )
print *,'u2*d2: ', norm(u3-MATMUL(u2,d2)),norm(u3 - (u2 .mm. d2) )
!-----------------------------------------------------------------------------------------------------------------------------------
call mmmult(m1,m3);
print *,'m1*m1^+:', norm(m3 - MATMUL(m1, hermitian(m1)))
call mmmult(m1,m4,side='H');
print *,'m1^+*m1:', norm(m4 - MATMUL(hermitian(m1), m1))
!-----------------------------------------------------------------------------------------------------------------------------------
call mmmult(d1,d3);
print *,'d1*d1^+:', norm(d3 - MATMUL(d1, symmetric(d1)))
call mmmult(d1,d4,side='S');
print *,'d1^+*d1:', norm(d4 - MATMUL(symmetric(d1), d1))
!-----------------------------------------------------------------------------------------------------------------------------------
call mmmult(m5,m6,m7,mtype='ZH')
print *,'m5*m6: (ZH)  ', norm( m7 - MATMUL(m5, m6)) , norm( m7 - (m5 .mm. m6))
!-----------------------------------------------------------------------------------------------------------------------------------
call mmmult(d5,d6,d7,mtype='DS')
print *,'d5*d6: (DS)  ', norm( d7 - MATMUL(d5, d6)) , norm( d7 - (d5 .mm. d6))
!-----------------------------------------------------------------------------------------------------------------------------------
m7 = inverse(m5)
print *,'m5^-1:  ', norm(diagonal(ONE,n) - MATMUL(m5,m7)),norm(diagonal(ONE,n) - MATMUL(m7,m5))
!-----------------------------------------------------------------------------------------------------------------------------------
d7 = inverse(d5)
print *,'d5^-1:  ', norm(diagonal(1.0_dp,n) - MATMUL(d5,d7)),norm(diagonal(1.0_dp,n) - MATMUL(d7,d5))
!-----------------------------------------------------------------------------------------------------------------------------------
call random_number(m5,sigma=1.0_dp);call random_number(d5,sigma=1.0_dp)
print *,'det(m5):', abs(determinant(m5)-PRODUCT(eigenvalues(m5)))/abs(determinant(m5))
print *,'det(d5):', abs(determinant(d5)-PRODUCT(eigenvalues(d5)))/abs(determinant(d5))
call antisymmetric_set(m5)
print *,'Pf(m5) :', abs(determinant(m5)-Pfaffian(m5)**2)/abs(determinant(m5))
!-----------------------------------------------------------------------------------------------------------------------------------
call random_number(m5,sigma=1.0_dp);
w = eigenvectors(m5);
m7=w(:,2:);v7=w(:,1)
print *,'evec:   ',norm(MATMUL(m5,m7)-MATMUL(m7,diagonal(v7)))

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
