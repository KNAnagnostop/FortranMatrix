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
!m1 = Matrix('g',m,n,is,js,name='d1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!-----------------------------------------------------------------------------------------------------------------------------------
!m=71;n=m+22;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!allocate(am1(is:is+m-1,js:js+n-1)); call random_number(am1,sigma=PI); call print(am1,stdout,fmt='F12.3',form=form,name='am1')
!allocate(am2(js:js+n-1,is:is+m-1)); call random_number(am2,sigma=PI); call print(am2,stdout,fmt='F12.3',form=form,name='am2')
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
!print *,'----------------------------------------------------------'
!print *,'am1 .mm. am2', norm((am1 .mm. am2)  - MATMUL(am1,am2) )
!print *,'am2 .mm. am1', norm((am2 .mm. am1)  - MATMUL(am2,am1) )
!print *,'----------------------------------------------------------'
!print *,'am1 .mm. ad2', norm((am1 .mm. ad2)  - MATMUL(am1,ad2) )
!print *,'ad2 .mm. am1', norm((ad2 .mm. am1)  - MATMUL(ad2,am1) )
!print *,'----------------------------------------------------------'
!print *,'ad1 .mm. am2', norm((ad1 .mm. am2)  - MATMUL(ad1,am2) )
!print *,'am2 .mm. ad1', norm((am2 .mm. ad1)  - MATMUL(am2,ad1) )
!print *,'----------------------------------------------------------'
!print *,'ad1 .mm. ad2', norm((ad1 .mm. ad2)  - MATMUL(ad1,ad2) )
!print *,'ad2 .mm. ad1', norm((ad2 .mm. ad1)  - MATMUL(ad2,ad1) )
!print *,'----------------------------------------------------------'
!print *,'am1 .mm. vm2', norm((am1 .mm. vm2)  - MATMUL(am1,vm2) )
!print *,'vm1 .mm. am1', norm((vm1 .mm. am1)  - MATMUL(vm1,am1) )
!print *,'----------------------------------------------------------'
!print *,'am1 .mm. vd2', norm((am1 .mm. vd2)  - MATMUL(am1,vd2) )
!print *,'vd1 .mm. am1', norm((vd1 .mm. am1)  - MATMUL(vd1,am1) )
!print *,'----------------------------------------------------------'
!print *,'ad1 .mm. vm2', norm((ad1 .mm. vm2)  - MATMUL(ad1,vm2) )
!print *,'vm1 .mm. ad1', norm((vm1 .mm. ad1)  - MATMUL(vm1,ad1) )
!print *,'----------------------------------------------------------'
!print *,'ad1 .mm. vd2', norm((ad1 .mm. vd2)  - MATMUL(ad1,vd2) )
!print *,'vd1 .mm. ad1', norm((vd1 .mm. ad1)  - MATMUL(vd1,ad1) )
!print *,'----------------------------------------------------------'
!-----------------------------------------------------------------------------------------------------------------------------------
! recheck after MATMUL -> LMATMUL
!m=256;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call d1 % print(stdout,fmt='F12.3',form=form);
!d2 =DMatrix('g',n,m,js,is,name='d2') ;call d2 % print(stdout,fmt='F12.3',form=form);
!d3 = d1 * d2                         ;call d3 % print(stdout,fmt='F12.3',form=form);
!d4 = d2 * d1                         ;call d4 % print(stdout,fmt='F12.3',form=form);
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!m2 = Matrix('g',n,m,js,is,name='m2') ;call m2 % print(stdout,fmt='F12.3',form=form);
!m3 = m1 * m2                         ;call m3 % print(stdout,fmt='F12.3',form=form);
!m4 = m2 * m1                         ;call m4 % print(stdout,fmt='F12.3',form=form);
!m5 = d1 * m2                         ;call m5 % print(stdout,fmt='F12.3',form=form);
!m6 = m1 * d2                         ;call m6 % print(stdout,fmt='F12.3',form=form);
!print * ,'d1*d2= ', norm(  d3 % v - MATMUL(d1 % v , d2 % v)  )
!print * ,'d2*d1= ', norm(  d4 % v - MATMUL(d2 % v , d1 % v)  )
!print * ,'m1*m2= ', norm(  m3 % v - MATMUL(m1 % v , m2 % v)  )
!print * ,'m2*m1= ', norm(  m4 % v - MATMUL(m2 % v , m1 % v)  )
!print * ,'d1*m2= ', norm(  m5 % v - MATMUL(d1 % v , m2 % v)  )
!print * ,'m1*d2= ', norm(  m6 % v - MATMUL(m1 % v , d2 % v)  )
!-----------------------------------------------------------------------------------------------------------------------------------
! recheck after MATMUL -> LMATMUL
!m=214;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!d1 =DMatrix('g',m,n,is,js,name='d1') ;call d1 % symmetric_set ; call d1 % print(stdout,fmt='F12.3',form=form);
!d2 =DMatrix('g',n,m,js,is,name='d2') ;call d2 % symmetric_set ; call d2 % print(stdout,fmt='F12.3',form=form);
!d3 = d1 * d2                         ;call d3 % print(stdout,fmt='F12.3',form=form);
!d4 = d2 * d1                         ;call d4 % print(stdout,fmt='F12.3',form=form);
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % hermitian_set ; call m1 % print(stdout,fmt='F12.3',form=form);
!m2 = Matrix('g',n,m,js,is,name='m2') ;call m2 % hermitian_set ; call m2 % print(stdout,fmt='F12.3',form=form);
!m3 = m1 * m2                         ;call m3 % print(stdout,fmt='F12.3',form=form);
!m4 = m2 * m1                         ;call m4 % print(stdout,fmt='F12.3',form=form);
!m5 = d1 * m2                         ;call m5 % print(stdout,fmt='F12.3',form=form);
!m6 = m1 * d2                         ;call m6 % print(stdout,fmt='F12.3',form=form);
!print * ,'d1*d2= ', norm(  d3 % v - MATMUL(d1 % v , d2 % v)  )
!print * ,'d2*d1= ', norm(  d4 % v - MATMUL(d2 % v , d1 % v)  )
!print * ,'m1*m2= ', norm(  m3 % v - MATMUL(m1 % v , m2 % v)  )
!print * ,'m2*m1= ', norm(  m4 % v - MATMUL(m2 % v , m1 % v)  )
!print * ,'d1*m2= ', norm(  m5 % v - MATMUL(d1 % v , m2 % v)  )
!print * ,'m1*d2= ', norm(  m6 % v - MATMUL(m1 % v , d2 % v)  )
!-----------------------------------------------------------------------------------------------------------------------------------
m=412;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if

!m1 =  Matrix('g',m,n,is,js,name='m1') ;call m1 % hermitian_set ; call m1 % print(stdout,fmt='F12.3',form=form);
!print *, isHermitian(m1) 
!m1 =  Matrix('g',m,n,is,js,name='m1') ;call symmetric_set(m1) ; call m1 % print(stdout,fmt='F12.3',form=form);
!print *, isSymmetric(m1) 
!m1 =  Matrix('g',m,n,is,js,name='m1') ;call antisymmetric_set(m1) ; call m1 % print(stdout,fmt='F12.3',form=form);
!print *, isAntiSymmetric(m1) 
!m1 =  Matrix('g',m,n,is,js,name='m1') ;call traceless_set(m1) ; call m1 % print(stdout,fmt='F12.3',form=form);
!print *, trace(m1) , abs(trace2(m1)-trace2c(m1))

!d1 = DMatrix('g',m,n,is,js,name='d1') ;call symmetric_set(d1) ; call d1 % print(stdout,fmt='F12.3',form=form);
!print *, isSymmetric(d1) 

!d1 = DMatrix('g',m,n,is,js,name='d1') ;call antisymmetric_set(d1) ; call d1 % print(stdout,fmt='F12.3',form=form);
!print *, isAntiSymmetric(d1) 

!d1 = DMatrix('g',m,n,is,js,name='d1') ;call traceless_set(d1) ; call d1 % print(stdout,fmt='F12.3',form=form);
!print *, trace(d1), abs(trace2(d1)-trace2c(d1)) 

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
