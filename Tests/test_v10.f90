program             testme
 use                matrix_mod_common
 use                matrix_mod_matrix
 use                matrix_mod_array
 implicit none
 integer   ,parameter   :: f_out = 13, f_math = 14
 integer                :: i,j
 real   (dp)            :: r,r1,r2,r3,r4,r5,r6,r7,r8,r9,sigma, ldetd(2),lPfd(2)
 complex(dp)            :: z,z1,z2,z3,z4,z5,z6,z7,z8,z9,       ldet (2),lPf (2)
 type( Matrix)          :: m1, m2, m3, m4, m5, m6, m7, m8, m9,II
 type(DMatrix)          :: d1, d2, d3, d4, d5, d6, d7, d8, d9,ID
 type( Vector)          :: v1, v2, v3, v4, v5, v6, v7, v8, v9
 type(DVector)          :: u1, u2, u3, u4, u5, u6, u7, u8, u9
 integer                :: m,n,is,ie,js,je
 real   (8),allocatable :: ad1(:,:  ),ad2(:,:  ),ad3(:,:  ),ad4(:,:  ),ad5(:,:)
 complex(8),allocatable :: am1(:,:  ),am2(:,:  ),am3(:,:  ),am4(:,:  ),am5(:,:),am6(:,:),am7(:,:),am8(:,:),am9(:,:),am10(:,:)
 real   (8),allocatable :: vd1(:    ),vd2(:    ),vd3(:    ),vd4(:    )
 complex(8),allocatable :: vm1(:    ),vm2(:    ),vm3(:    ),vm4(:    ),vm5(:  ),vm6(:  ),vm7(:  ),vm8(:  ),vm9(:  ),vm10(:  )
 character(1000)        :: string
 character(1)           :: form
!-----------------------------------------------------------------------------------------------------------------------------------
 open(unit=f_out ,file='test.dat');  open(unit=f_math,file='test.m'); f_mout   =f_out ! output unit for matrix modules set to f_out
 call matrix_random_init
 f_mout = stdout
!-----------------------------------------------------------------------------------------------------------------------------------
!m=4;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call d1 % print(stdout,fmt='F12.3',form=form);
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!-----------------------------------------------------------------------------------------------------------------------------------
! m=612;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
! allocate(am1(is:is+m-1,js:js+n-1)); call random_number(am1,sigma=PI); call print(am1,stdout,fmt='F12.3',form=form,name='am1')
! allocate(am2(js:js+n-1,is:is+m-1)); call random_number(am2,sigma=PI); call print(am2,stdout,fmt='F12.3',form=form,name='am2')
!!if we allocate am3, then is:is+m-1 remains after the multiplication:
!allocate(am3(is:is+m-1,is:is+m-1)); call random_number(am3,sigma=PI); call print(am3,stdout,fmt='F12.3',form=form,name='am3')
!am3 = am1 .mm. am2                                                  ; call print(am3,stdout,fmt='F12.3',form=form,name='am1*am2') 
!am4 = am2 .mm. am1                                                  ; call print(am4,stdout,fmt='F12.3',form=form,name='am2*am1')
!allocate(vm1(is:is+m-1          )); call random_number(vm1,sigma=PI); call print(vm1,stdout,fmt='F12.3',form=form,name='vm1')
!allocate(vm2(js:js+n-1          )); call random_number(vm2,sigma=PI); call print(vm2,stdout,fmt='F12.3',form=form,name='vm2')
!allocate(ad1(is:is+m-1,js:js+n-1)); call random_number(ad1,sigma=PI); call print(ad1,stdout,fmt='F12.3',form=form,name='ad1')
!allocate(ad2(js:js+n-1,is:is+m-1)); call random_number(ad2,sigma=PI); call print(ad2,stdout,fmt='F12.3',form=form,name='ad2')
!allocate(vd1(is:is+m-1          )); call random_number(vd1,sigma=PI); call print(vd1,stdout,fmt='F12.3',form=form,name='vd1')
!allocate(vd2(js:js+n-1          )); call random_number(vd2,sigma=PI); call print(vd2,stdout,fmt='F12.3',form=form,name='vd2')
!-----------------------------------------------------------------------------------------------------------------------------------
!call random_number(am1,sigma=PI); call print(am1,stdout,fmt='F12.3',form=form,name='am1')
!call hermitian_set(am1,'l');      call print(am1,stdout,fmt='F12.3',form=form,name='am1')
!print *,'am1: ',isHermitian(am1)
!call random_number(am1,sigma=PI); call print(am1,stdout,fmt='F12.3',form=form,name='am1')
!call symmetric_set(am1,'l');      call print(am1,stdout,fmt='F12.3',form=form,name='am1')
!print *,'am1: ',isSymmetric(am1)
!call random_number(am1,sigma=PI); call print(am1,stdout,fmt='F12.3',form=form,name='am1')
!call antisymmetric_set(am1,'l');  call print(am1,stdout,fmt='F12.3',form=form,name='am1')
!print *,'am1: ',isAntiSymmetric(am1)
!call random_number(ad1,sigma=PI); call print(ad1,stdout,fmt='F12.3',form=form,name='ad1')
!call symmetric_set(ad1,'l');      call print(ad1,stdout,fmt='F12.3',form=form,name='ad1')
!print *,'ad1: ',isSymmetric(ad1)
!call random_number(ad1,sigma=PI); call print(ad1,stdout,fmt='F12.3',form=form,name='ad1')
!call antisymmetric_set(ad1,'l');  call print(ad1,stdout,fmt='F12.3',form=form,name='ad1')
!print *,'ad1: ',isAntiSymmetric(ad1)
!-----------------------------------------------------------------------------------------------------------------------------------
!am2 = traceless(am1);  call print(am2,stdout,fmt='F12.3',form=form,name='am2')
!am3 = traceless(am1);  call print(am3,stdout,fmt='F12.3',form=form,name='am3')
!print *,'am2: ', trace(am2)/m, abs(trace2(am2)-trace2c(am2))/m
!print *,'am3: ', trace(am3)/m, abs(trace2(am3)-trace2c(am3))/m
!m1 = Matrix('g',m,n,is,js,name='d1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!m2 = traceless(m1)
!print *,' m2: ', trace(m2)/m, abs(trace2(m2)-trace2c(m2))/m
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call d1 % print(stdout,fmt='F12.3',form=form);
!d2 = traceless(d1)
!print *,' d2: ', trace(d2)/m, abs(trace2(d2)-trace2c(d2))/m
!-----------------------------------------------------------------------------------------------------------------------------------
! m=612;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!call hermitian_set(m1,'l')           ;call m1 % print(stdout,fmt='F12.3',form=form);
!print *,' m1: ', isHermitian(m1)
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!call symmetric_set(m1,'l')           ;call m1 % print(stdout,fmt='F12.3',form=form);
!print *,' m1: ', isSymmetric(m1)
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!call antisymmetric_set(m1,'l')       ;call m1 % print(stdout,fmt='F12.3',form=form);
!print *,' m1: ', isAntiSymmetric(m1)
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call d1 % print(stdout,fmt='F12.3',form=form);
!call symmetric_set(d1,'u')           ;call d1 % print(stdout,fmt='F12.3',form=form);
!print *,' d1: ', isSymmetric(d1)
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call d1 % print(stdout,fmt='F12.3',form=form);
!call antisymmetric_set(d1,'u')       ;call d1 % print(stdout,fmt='F12.3',form=form);
!print *,' d1: ', isAntiSymmetric(d1)
!-----------------------------------------------------------------------------------------------------------------------------------
!m=612;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js,name='m1') ;call hermitian_set(m1); call m1 % print(stdout,fmt='F12.3',form=form);
!m2 = sin(m1)                         ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'sin (m2)? ', isHermitian(m2), norm(m2 % v - sin(m1 % v))
!m2 = cos(m1)                         ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'cos (m2)? ', isHermitian(m2), norm(m2 % v - cos (m1 % v))
!m2 = exp(m1)                         ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'exp (m2)? ', isHermitian(m2), norm(m2 % v - exp (m1 % v))
!m2 = log(m1)                         ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'log (m2)? ', isHermitian(m2), norm(m2 % v - log (m1 % v))
!m2 = sqrt(m1)                        ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'sqrt(m2)? ', isHermitian(m2), norm(m2 % v - sqrt(m1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!m=4;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js,name='m1') ;call hermitian_set(m1); call m1 % print(stdout,fmt='F12.3',form=form);
!d2 = abs(m1)                         ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'abs (m1)? ',  issymmetric(d2), norm(d2 % v - abs(m1 % v))
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call symmetric_set(d1); call d1 % print(stdout,fmt='F12.3',form=form);
!d2 = abs(d1)                         ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'abs (d2)? ', issymmetric(d2), norm(d2 % v - abs(d1 % v))
!d2 = sin(d1)                         ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'sin (d2)? ', issymmetric(d2), norm(d2 % v - sin(d1 % v))
!d2 = cos(d1)                         ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'cos (d2)? ', issymmetric(d2), norm(d2 % v - cos (d1 % v))
!d2 = exp(d1)                         ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'exp (d2)? ', issymmetric(d2), norm(d2 % v - exp (d1 % v))
!d1 = abs(d1)
!d2 = log(d1)                         ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'log (d2)? ', issymmetric(d2), norm(d2 % v - log (d1 % v))
!d2 = sqrt(d1)                        ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'sqrt(d2)? ', issymmetric(d2), norm(d2 % v - sqrt(d1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!m=6;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!v1 = Vector('g',m,is,name='v1') ; call v1 % print(stdout,fmt='F12.3',form=form);
!v2 = sin(v1)                         ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'sin (v2)? ', norm(v2 % v - sin(v1 % v))
!v2 = cos(v1)                         ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'cos (v2)? ', norm(v2 % v - cos (v1 % v))
!v2 = exp(v1)                         ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'exp (v2)? ', norm(v2 % v - exp (v1 % v))
!v2 = log(v1)                         ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'log (v2)? ', norm(v2 % v - log (v1 % v))
!v2 = sqrt(v1)                        ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'sqrt(v2)? ', norm(v2 % v - sqrt(v1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!m=48;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!v1 = Vector('g',m,is,name='v1') ;call v1 % print(stdout,fmt='F12.3',form=form);
!u2 = abs(v1)                         ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'abs (v1)? ',   norm(u2 % v - abs(v1 % v))
!u1 =DVector('g',m,is,name='u1') ;call u1 % print(stdout,fmt='F12.3',form=form);
!u2 = abs(u1)                         ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'abs (u2)? ',  norm(u2 % v - abs(u1 % v))
!u2 = sin(u1)                         ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'sin (u2)? ',  norm(u2 % v - sin(u1 % v))
!u2 = cos(u1)                         ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'cos (u2)? ',  norm(u2 % v - cos (u1 % v))
!u2 = exp(u1)                         ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'exp (u2)? ',  norm(u2 % v - exp (u1 % v))
!u1 = abs(u1)
!u2 = log(u1)                         ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'log (u2)? ',  norm(u2 % v - log (u1 % v))
!u2 = sqrt(u1)                        ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'sqrt(u2)? ',  norm(u2 % v - sqrt(u1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!m  =49;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js,name='m1') ;call hermitian_set(m1); call m1 % print(stdout,fmt='F12.3',form=form);
!m2 =     m1**12                      ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'m1**12  ? ', isHermitian(m2), norm(m2 % v - (m1 % v)**12 )
!m2 =     m1**PI                      ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'m1**PI  ? ', isHermitian(m2), norm(m2 % v - (m1 % v)**PI )
!call random_number(z1)
!m2 =     m1**z1                      ;call m2 % print(stdout,fmt='F12.3',form=form);
!print *,'m1**z1  ? ', isHermitian(m2), norm(m2 % v - (m1 % v)**z1 )
!-----------------------------------------------------------------------------------------------------------------------------------
!m  =49;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call symmetric_set(d1); call d1 % print(stdout,fmt='F12.3',form=form);
!d2 =     d1**12                      ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'d1**12  ? ', isSymmetric(d2), norm(d2 % v - (d1 % v)**12 )
!d1 = abs(d1)
!d2 =     d1**PI                      ;call d2 % print(stdout,fmt='F12.3',form=form);
!print *,'d1**PI  ? ', isSymmetric(d2), norm(d2 % v - (d1 % v)**PI )
!-----------------------------------------------------------------------------------------------------------------------------------
!m  =4;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!v1 = Vector('g',m,is,name='v1')      ;call v1 % print(stdout,fmt='F12.3',form=form);
!v2 =     v1**12                      ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'v1**12  ? ', norm(v2 % v - (v1 % v)**12 )
!v2 =     v1**PI                      ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'v1**PI  ? ', norm(v2 % v - (v1 % v)**PI )
!call random_number(z1)
!v2 =     v1**z1                      ;call v2 % print(stdout,fmt='F12.3',form=form);
!print *,'v1**z1  ? ', norm(v2 % v - (v1 % v)**z1 )
!-----------------------------------------------------------------------------------------------------------------------------------
!m  =4;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!u1 =DVector('g',m,is,name='u1')      ;call u1 % print(stdout,fmt='F12.3',form=form);
!u2 =     u1**12                      ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'u1**12  ? ', norm(u2 % v - (u1 % v)**12 )
!u1 = abs(u1)
!u2 =     u1**PI                      ;call u2 % print(stdout,fmt='F12.3',form=form);
!print *,'u1**PI  ? ', norm(u2 % v - (u1 % v)**PI )
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
