program             testme
 use                matrix_mod_common
 use                matrix_mod_matrix
 use                matrix_mod_array
 implicit none
 integer   ,parameter   :: f_out = 13, f_math = 14
 integer                :: i,j
 real   (dp)            :: r,r1,r2,r3,r4,r5,r6,r7,r8,r9,sigma
 complex(dp)            :: z,z1,z2,z3,z4,z5,z6,z7,z8,z9
 type( Matrix)          :: m1, m2, m3
 type(DMatrix)          :: d1, d2, d3
 type( Vector)          :: v1, v2, v3, v4, v5, v6, v7
 type(DVector)          :: u1, u2, u3, u4, u5, u6, u7
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
 m=10;  n=m;  is=13;js=21            ; ie = is + m - 1; je = js + n -1
 if(m<=10) then; form='d';else;form='i';end if ! call v1 % print(stdout,fmt='F12.3',form=form)
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!v1 = Vector(n       ,is=is,name='v1'); if(n<20) call v1 % print(stdout,fmt='F12.3',form='d')
!v1 = Vector(n,PI    ,is=is,name='v2'); if(n<20) call v1 % print(stdout,fmt='F12.3',form='d')
!v1 = Vector(n,PI+IMU,is=is,name='v3'); if(n<20) call v1 % print(stdout,fmt='F12.3',form='d'); 
!allocate(vm1(n)) ; call random_number(vm1)
!v1 = Vector(vm1     ,is=is,name='v4'); if(n<20) call v1 % print(stdout,fmt='F12.3',form='d')
!v1 = Vector('u',n   ,is=is,name='v5'); if(n<20) call v1 % print(stdout,fmt='F12.3',form='d')
!v1 = Vector('g',n   ,is=is,name='v5'); if(n<20) call v1 % print(stdout,fmt='F12.3',form='d')
!n  = 1000000;v1 = Vector('g',n   ,is=is,name='v5',sigma=PI);  call v1 % print(20)
!gnuplot> s = pi; f(x) = 1.0/sqrt(2*pi*s**2)*exp(-0.5*x**2/s**2);plot [:][0:] "< awk '!/^#/{print $1}' fort.20|hist -d 20000000 -b 0.05" u 1:4 w l,f(x)
!-----------------------------------------------------------------------------------------------------------------------------------
!allocate(vm1(n));call random_number(vm1); v1 = Vector(vm1,is=is,name='v5');call v1 % print(stdout,fmt='F12.3',form='i')
!v1 = Vector('u',n   ,is=is,name='v5'); if(n<20) call v1 % print(stdout,fmt='F12.3',form='i')
!open(unit=11,file='test.con');call v1 % save(11);close(11);open(unit=11,file='test.con');v2 = Vector(n,is=225,name='v2');call v2 % read(11);print *,'|v1-v2|= ',norm(v1 % v - v2 % v)
!-----------------------------------------------------------------------------------------------------------------------------------
!u1 = DVector(n       ,is=is,name='u1'); if(n<20) call u1 % print(stdout,fmt='F12.3',form='d')
!u1 = DVector(n,PI    ,is=is,name='u2'); if(n<20) call u1 % print(stdout,fmt='F12.3',form='d')
!u1 = DVector(n,PI+IMU,is=is,name='u3'); if(n<20) call u1 % print(stdout,fmt='F12.3',form='d')
!allocate(vd1(n)) ; call random_number(vd1)
!u1 = DVector(vd1     ,is=is,name='u4'); if(n<20) call u1 % print(stdout,fmt='F12.3',form='d')
!u1 = DVector('u',n   ,is=is,name='u5'); if(n<20) call u1 % print(stdout,fmt='F12.3',form='d')
!u1 = DVector('g',n   ,is=is,name='u6'); if(n<20) call u1 % print(stdout,fmt='F12.3',form='d')
!n=1000000; u1=DVector('g',n   ,is=is,name='u5',sigma=PI); call u1 % save(20); call u1 % print(stdout,form='I')
!-----------------------------------------------------------------------------------------------------------------------------------
!open(unit=11,file='test.con');call u1  % save(11);close(11);open(unit=11,file='test.con');u2 = DVector(n,is=225,name='u2'); call u2 % read(11);print *,'|u1-u2|= ',norm(u1 % v - u2 % v)
!u1 = DVector(n, is=is, name='u1'); call u1 % random();call u1 % gaussian(sigma=PI);call u1 % save(20)
!v1 =  Vector(n, is=is, name='v1'); call v1 % random();call v1 % gaussian(sigma=PI);call v1 % save(21)
!-----------------------------------------------------------------------------------------------------------------------------------
!allocate(vm1(12:12+n-1));call random_number(vm1);v1 = Vector(vm1,is=lbound(vm1,1));print *,'norm:',norm(v1)-sum(abs(v1 % v))/v1 % n
!u1 = v1 % re();u2 = v1 % im();print *, norm(v1 % v - CMPLX(u1 % v,u2 % v,kind=dp))
!v2 = v1 % conjg();print * , norm( conjg(v1 % v) - v2 % v)
!allocate(vd1(12:12+n-1));call random_number(vd1);u1 =DVector(vd1,is=lbound(vd1,1));print *,'norm:',norm(u1)-sum(abs(u1 % v))/u1 % n
!-----------------------------------------------------------------------------------------------------------------------------------
!n=622;is=258;if(n<10) then; form='d';else;form='i';end if
!v1 = Vector('u',n,is=is,name='v1'); call v1  % print(stdout,fmt='F12.3',form=form)
!v2 = v1 ;                           call v2  % print(stdout,fmt='F12.3',form=form); 
!print *,'|v2-v1|=',norm(v2 % v - v1 % v)
!allocate(vm3(123:123+n-1)); call random_number(vm3)
!v3 = v2
!v3 = vm3;                           call v3  % print(stdout,fmt='F12.3',form=form); 
!print *,'|v3-vm3|=',norm(v3 % v - vm3)
!v4 = v2
!v4 = PI;                            call v4  % print(stdout,fmt='F12.3',form=form); 
!print *,'|v4- r|=',norm(v4 % v - PI)
!v4 = PI+IMU;                        call v4  % print(stdout,fmt='F12.3',form=form); 
!print *,'|v4- z|=',norm(v4 % v - (PI+IMU))
!u1 = DVector('u',n,is=is+4,name='u1');
!v5 = u1;                            call v5  % print(stdout,fmt='F12.3',form=form); 
!print *,'|v5-u1|=',norm(v5 % v - u1 % v)
!u2 = u1;                            call u2  % print(stdout,fmt='F12.3',form=form); 
!print *,'|u2-u1|=',norm(u2 % v - u1 % v)
!allocate(vd3(223:223+n-1)); call random_number(vd3)
!u3 = u2
!u3 = vd3;                           call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3-vm3|=',norm(u3 % v - vd3)
!u4 = u2
!u4 = PI;                            call u4  % print(stdout,fmt='F12.3',form=form); 
!print *,'|u4- r|=',norm(u4 % v - PI)
!u4 = PI+IMU;                        call u4  % print(stdout,fmt='F12.3',form=form); 
!print *,'|u4- z|=',norm(u4 % v - PI)
!-----------------------------------------------------------------------------------------------------------------------------------
!n=622;is=258;if(n<10) then; form='d';else;form='i';end if
!v1 =  Vector('u',n,is=is+4,name='u1');call v1  % print(stdout,fmt='F12.3',form=form);
!v2 =  Vector('u',n,is=is  ,name='u2');call v2  % print(stdout,fmt='F12.3',form=form);  
!allocate(vm1(223:223+n-1)); call random_number(vm1)
!v3 = v1 + v2                         ;call v3  % print(stdout,fmt='F12.3',form=form);
!print *,'|v1+v2|',norm( v3 % v - (v1 % v + v2 % v))
!v4 = v3 + vm1                        ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3+vm|',norm( v4 % v - (v3 % v + vm1))
!v4 = vm1+v3                          ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|vm+v3|',norm( v4 % v - (v3 % v + vm1))
!v5 = v3 +PI                          ;call v5  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3+r |',norm( v5 % v - (v3 % v + PI))
!v5 = PI +v3                          ;call v5  % print(stdout,fmt='F12.3',form=form);
!print *,'|r +v3|',norm( v5 % v - (v3 % v + PI))
!v5 = v3 +(PI+IMU)                    ;call v5  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3+z |',norm( v5 % v - (v3 % v + (PI+IMU)))
!v5 = (PI+IMU) +v3                    ;call v5  % print(stdout,fmt='F12.3',form=form);
!print *,'|z +v3|',norm( v5 % v - (v3 % v + PI))
!u1 = DVector('u',n,is=12,name='u1')
!v6 = v3 + u1                         ;call v6  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3+u1|',norm( v6 % v - (v3 % v + u1 % v))
!v7 = u1 + v3                         ;call v7  % print(stdout,fmt='F12.3',form=form);
!print *,'|u1+v7|',norm( v7 % v - (v3 % v + u1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!u1 = DVector('u',n,is=is+4,name='u1');call u1  % print(stdout,fmt='F12.3',form=form);
!u2 = DVector('u',n,is=is  ,name='u2');call u2  % print(stdout,fmt='F12.3',form=form);  
!allocate(vd1(223:223+n-1)); call random_number(vd1)
!u3 = u1 + u2                         ;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|u1+u2|',norm( u3 % v - (u1 % v + u2 % v))
!u4 = u3 + vd1                        ;call u4  % print(stdout,fmt='F12.3',form=form);
!print *,'|u3+vd|',norm( u4 % v - (u3 % v + vd1))
!u4 = vd1+u3                          ;call u4  % print(stdout,fmt='F12.3',form=form);
!print *,'|vd+u3|',norm( u4 % v - (u3 % v + vd1))
!u5 = u3 +PI                          ;call u5  % print(stdout,fmt='F12.3',form=form);
!print *,'|u3+r |',norm( u5 % v - (u3 % v + PI))
!u5 = PI +u3                          ;call u5  % print(stdout,fmt='F12.3',form=form);
!print *,'|r +u3|',norm( u5 % v - (u3 % v + PI))
!!test dvector_assignFrom_vector
!u6 = v7                              ;call u6  % print(stdout,fmt='F12.3',form=form);
!print *,'|u6-Re(v7)|',norm(u6 % v - Real(v7 % v, kind=dp) )
!-----------------------------------------------------------------------------------------------------------------------------------
!n=6;is=258;if(n<10) then; form='d';else;form='i';end if
!v1 =  Vector('u',n,is=is+4   ,name='v1');call v1  % print(stdout,fmt='F12.3',form=form);
!v2 =  Vector('u',n,is=is     ,name='v2');call v2  % print(stdout,fmt='F12.3',form=form);
!u1 = DVector('u',n,is=is+9321,name='u2');
!allocate(vm1(223:223+n-1)); call random_number(vm1)
!v3 = v1 -v2                         ;call v3  % print(stdout,fmt='F12.3',form=form);
!print *,'|v1-v2|', norm(v3 % v - (v1 % v - v2 % v) )
!v3 = v2 -v1                         ;call v3  % print(stdout,fmt='F12.3',form=form);
!print *,'|v2-v1|', norm(v3 % v + (v1 % v - v2 % v) )
!v4 = vm1-v3                         ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|vm-v3|',norm( v4 % v - (-v3 % v + vm1))
!v4 = v3 - vm1                       ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3-vm|',norm( v4 % v + (-v3 % v + vm1))
!v4 = v3 - u1                        ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3-u1|',norm( v4 % v - ( v3 % v - u1 % v))
!v4 = u1 - v3                        ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|u1-v3|',norm( v4 % v + ( v3 % v - u1 % v))
!v4 = v3 - PI                        ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3- r|',norm( v4 % v - ( v3 % v - PI))
!v4 = PI - v3                       ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|r -v3|',norm( v4 % v + ( v3 % v - PI))
!v4 = v3 - (PI+IMU)                 ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|v3- z|',norm( v4 % v - ( v3 % v - (PI+IMU)))
!v4 = (PI+IMU) - v3                 ;call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|z -v3|',norm( v4 % v + ( v3 % v - (PI+IMU)))
!v4 = -v1                           ;call v1  % print(stdout,fmt='F12.3',form=form);call v4  % print(stdout,fmt='F12.3',form=form);
!print *,'|-v1  |',norm( v4 % v + v1 % v)
!-----------------------------------------------------------------------------------------------------------------------------------
!n=6;is=258;if(n<10) then; form='d';else;form='i';end if
!u1 = DVector('u',n,is=is+4   ,name='u1');call u1  % print(stdout,fmt='F12.3',form=form);
!u2 = DVector('u',n,is=is+4   ,name='u2');call u2  % print(stdout,fmt='F12.3',form=form);
!call random_number(r1); allocate(vd1(223:223+n-1)); call random_number(vd1)
!u3 = u1 - u2;                           ;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|u1-u2|', norm(u3 % v - (u1 % v - u2 % v) )
!u3 = u1 - vd1;                          ;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|u1-v |', norm(u3 % v - (u1 % v - vd1) )
!u3 = vd1 - u1;                          ;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|v -u1|', norm(u3 % v + (u1 % v - vd1) )
!u3 = u1  - r1;                          ;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|u1- r|', norm(u3 % v - (u1 % v - r1) )
!u3 = r1  - u1;                          ;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|r -u1|', norm(u3 % v + (u1 % v - r1) )
!u3 = -u1                                ;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'|-u1|'  , norm(u3 % v + (u1 % v)      )
!-----------------------------------------------------------------------------------------------------------------------------------
!n=6;is=258;if(n<10) then; form='d';else;form='i';end if
!v1 =  Vector('u',n,is=is+4   ,name='v1');call v1  % print(stdout,fmt='F12.3',form=form);
!call random_number(r1); call random_number(z1)
!v2 = v1 / r1                            ;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'v2/r', norm( v2 % v - v1 % v / r1)
!v2 = v1 / z1                            ;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'v2/z', norm( v2 % v - v1 % v / z1)
!v2 = v1 * r1                            ;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'v2*r', norm( v2 % v - v1 % v * r1)
!v2 = v1 * z1                            ;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'v2*z', norm( v2 % v - v1 % v * z1)
!v2 = r1 * v1                            ;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'r*v2', norm( v2 % v - v1 % v * r1)
!v2 = z1 * v1                            ;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'z*v2', norm( v2 % v - v1 % v * z1)
!-----------------------------------------------------------------------------------------------------------------------------------
!n=6;is=258;if(n<10) then; form='d';else;form='i';end if
!u1 = DVector('u',n,is=is+4   ,name='u1');call u1  % print(stdout,fmt='F12.3',form=form);
!call random_number(r1);
!u2 = u1/r1                              ;call u2  % print(stdout,fmt='F12.3',form=form);
!print *,'u2/r', norm( u2 % v - u1 % v / r1)
!u2 = u1*r1                              ;call u2  % print(stdout,fmt='F12.3',form=form);
!print *,'u2*r', norm( u2 % v - u1 % v * r1)
!u2 = r1*u1                              ;call u2  % print(stdout,fmt='F12.3',form=form);
!print *,'r*u2', norm( u2 % v - u1 % v * r1)
!-----------------------------------------------------------------------------------------------------------------------------------
!m=623;n=m;is=258;js=123;if(n<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1  % print(stdout,fmt='F12.3',form=form);
!v1 = Vector('g',n,is=is+83,name='v1');call v1  % print(stdout,fmt='F12.3',form=form);
!v2 = m1 * v1;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'m1*v1:(G)', norm( v2 % v - MATMUL(m1 % v, v1 % v))
!v3 = v1 * m1;call v3  % print(stdout,fmt='F12.3',form=form);
!print *,'v1*m1:(G)', norm( v3 % v - MATMUL(v1 % v, m1 % v))
!call m1 % hermitian_set ;call m1  % print(stdout,fmt='F12.3',form=form);
!v2 = m1 * v1;call v2  % print(stdout,fmt='F12.3',form=form);
!print *,'m1*v1:(H)', norm( v2 % v - MATMUL(m1 % v, v1 % v))
!v3 = v1 * m1;call v3  % print(stdout,fmt='F12.3',form=form);
!print *,'v1*m1:(H)', norm( v3 % v - MATMUL(v1 % v, m1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!d1 = DMatrix('g',m,n,is,js ,name='d1') ;call d1  % print(stdout,fmt='F12.3',form=form);
!u1 = DVector('g',n,is=is+83,name='u1') ;call u1  % print(stdout,fmt='F12.3',form=form);
!u2 = d1 * u1;call u2  % print(stdout,fmt='F12.3',form=form);
!print *,'d1*u1:   ', norm( u2 % v - MATMUL(d1 % v, u1 % v))
!u2 = u1 * d1;call u2  % print(stdout,fmt='F12.3',form=form);
!print *,'u1*d1:   ', norm( u2 % v - MATMUL(u1 % v, d1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!m=512;n=m+37;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js ,name='m1');call m1  % print(stdout,fmt='F12.3',form=form);
!v1 = Vector('g',n,is=is+83,name='v1');call v1  % print(stdout,fmt='F12.3',form=form);
!v2 = Vector('g',m,is=js+33,name='v2');call v2  % print(stdout,fmt='F12.3',form=form);
!v3 = m1 * v1;call v3  % print(stdout,fmt='F12.3',form=form);
!print *,'m1*v1:', norm( v3 % v - MATMUL(m1 % v, v1 % v))
!v3 = v2 * m1;call v3  % print(stdout,fmt='F12.3',form=form);
!print *,'v2*m1:', norm( v3 % v - MATMUL(v2 % v, m1 % v))
!d1 =DMatrix('g',m,n,is,js ,name='d1');call d1  % print(stdout,fmt='F12.3',form=form);
!u1 =DVector('g',n,is=is+83,name='u1');call u1  % print(stdout,fmt='F12.3',form=form);
!u2 =DVector('g',m,is=js+33,name='u2');call u2  % print(stdout,fmt='F12.3',form=form);
!u3 = d1 * u1;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'d1*u1:', norm( u3 % v - MATMUL(d1 % v, u1 % v))
!u3 = u2 * d1;call u3  % print(stdout,fmt='F12.3',form=form);
!print *,'u2*d1:', norm( u3 % v - MATMUL(u2 % v, d1 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!m=615;n=m+37;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m);call m1  % print(stdout,fmt='F12.3',form=form);
!d1 =DMatrix('g',m);call d1  % print(stdout,fmt='F12.3',form=form);
!v1 = Vector('g',m);call v1  % print(stdout,fmt='F12.3',form=form);
!u1 = Vector('g',m);call u1  % print(stdout,fmt='F12.3',form=form);
!print *,'m1,d1,v1,u1:', isNaN(m1 % v), isNaN(d1 % v), isNaN(v1 % v), isNaN(u1 % v)
!print *,'m1,d1,v1,u1:', isNaN(m1)    , isNaN(d1)    , isNaN(v1)    , isNaN(u1)
!m1 % v(1,1) = nan(z);call m1  % print(stdout,fmt='F12.3',form=form);
!d1 % v(1,1) = nan(z);call d1  % print(stdout,fmt='F12.3',form=form);
!v1 % v(1)   = nan( );call v1  % print(stdout,fmt='F12.3',form=form);
!u1 % v(1)   = nan( );call u1  % print(stdout,fmt='F12.3',form=form);
!print *,'m1,d1,v1,u1:', isNaN(m1)    , isNaN(d1)    , isNaN(v1)    , isNaN(u1)

!call u1 % print(stdout,fmt='F12.3',form=form)
!if(m<7) call print(am5,stdout,fmt='F12.3',name='diag')
!-----------------------------------------------------------------------------------------------------------------------------------
!m=6;n=m-2;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!v1 = Vector(m,is=is,name='v1');call random_number(v1,sigma=PI);call v1  % print(stdout,fmt='F12.3',form=form);
!u1 =DVector(m,is=is,name='u1');call random_number(u1,sigma=PI);call u1  % print(stdout,fmt='F12.3',form=form);
!m1 = Matrix(m,n,name='m1')    ;call random_number(m1,sigma=PI);call m1  % print(stdout,fmt='F12.3',form=form);
!d1 =DMatrix(m,n,name='d1')    ;call random_number(d1,sigma=PI);call d1  % print(stdout,fmt='F12.3',form=form);
!-----------------------------------------------------------------------------------------------------------------------------------
!m=6;n=m-2;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!v1 = Vector(m,is=is,name='v1');call random_number(v1,sigma=PI);call v1  % print(stdout,fmt='F12.3',form=form);
!v2 = sort(v1,by='IF') ;call v2  % print(stdout,fmt='F12.3',form=form);
!u1 = DVector(m,is=is,name='u1');call random_number(u1,sigma=PI);call u1  % print(stdout,fmt='F12.3',form=form);
!u2 = sort(u1,by='AF') ; call u2  % print(stdout,fmt='F12.3',form=form);
!-----------------------------------------------------------------------------------------------------------------------------------
!m=6;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!v1 = Vector('g',m);call v1  % print(stdout,fmt='F12.3',form=form);
!v2 = Vector('g',m);call v2  % print(stdout,fmt='F12.3',form=form);
!print *,dot_product(v1,v2), abs(dot_product(v1,v2) - dot_product(v1%v,v2%v))
!u1 =DVector('g',m);call u1  % print(stdout,fmt='F12.3',form=form);
!u2 =DVector('g',m);call u2  % print(stdout,fmt='F12.3',form=form);
!print *,dot_product(u1,u2), SUM(u1%v*u2%v),dot_product(u1%v,u2%v),abs(dot_product(u1,u2) - dot_product(u1%v,u2%v))
!print *, maxval(u1), maxval(u2, u2 % v < 0.0_dp)
!print *, minval(u1), minval(u2, u2 % v > 0.0_dp)
!-----------------------------------------------------------------------------------------------------------------------------------
!m=4;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1= Matrix('g',m,n,is,js,name='m1'); call m1  % print(stdout,fmt='F12.3',form=form);
!v1= Vector('g',m,is,name='v1')     ; call v1  % print(stdout,fmt='F12.3',form=form);
!u1=DVector('g',m,is,name='d1')     ; call u1  % print(stdout,fmt='F12.3',form=form);
!print *, abs(trace(m1)-trace(m1 % v)), abs(trace2(m1)-trace2(m1 % v)), abs(trace2c(m1)-trace2c(m1 % v))
!call m1 % hermitian_set;  call m1  % print(stdout,fmt='F12.3',form=form); !check hermitian flag
!print *, abs(trace(m1)-trace(m1 % v)), abs(trace2(m1)-trace2(m1 % v)), abs(trace2c(m1)-trace2c(m1 % v))
!call m1 % traceless_set;  call m1  % print(stdout,fmt='F12.3',form=form); !check hermitian flag
!print *, abs(trace(m1)),abs(trace(m1)-trace(m1 % v)), abs(trace2(m1)-trace2(m1 % v)), abs(trace2c(m1)-trace2c(m1 % v))
!v1=diagonal(m1)                    ; call v1  % print(stdout,fmt='F12.3',form=form);
!print *, norm(v1 % v - diagonal(m1 % v))
!m1=diagonal(v1)                    ; call m1  % print(stdout,fmt='F12.3',form=form);
!m2= Matrix('g',m,n,is,js,name='m2'); call m2  % print(stdout,fmt='F12.3',form=form);
!m1=diagonalMatrix(m2)              ; call m1  % print(stdout,fmt='F12.3',form=form);
!call traceless(m1)
!print *, trace(m1)
!-----------------------------------------------------------------------------------------------------------------------------------
!m=4;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!d1=DMatrix('u',m,n,is,js,name='d1'); call d1  % print(stdout,fmt='F12.3',form=form);
!print *, trace(d1), trace2(d1),trace2c(d1)
!print *,abs(trace(d1)-trace(d1 % v)), abs(trace2(d1)-trace2(d1 % v)), abs(trace2c(d1)-trace2c(d1 % v))
!call d1 % traceless_set !traceless(d1)
!print *,trace(d1), trace2(d1),trace2c(d1)
!print *,abs(trace(d1)-trace(d1 % v)), abs(trace2(d1)-trace2c(d1))
!u1 = diagonal(d1)                  ; call u1  % print(stdout,fmt='F12.3',form=form);
!print *, norm(u1 % v - diagonal(d1 % v))
!d2 = diagonal(u1)                  ; call d2  % print(stdout,fmt='F12.3',form=form);
!d2 = diagonalMatrix(d1)            ; call d2  % print(stdout,fmt='F12.3',form=form);
!call d1 % symmetric_set
!d2 = diagonalMatrix(d1)            ; call d2  % print(stdout,fmt='F12.3',form=form);
!d1=DMatrix('u',m,n,is,js,name='d1'); call d1  % symmetric_set; call d1  % print(stdout,fmt='F12.3',form=form);
!print *,trace(d1), trace2(d1),trace2c(d1)
!print *,abs(trace(d1)-trace(d1 % v)), abs(trace2(d1)-trace2(d1 % v)), abs(trace2c(d1)-trace2c(d1 % v))
!call d1 % traceless_set !traceless(d1)
!print *,trace(d1), trace2(d1),trace2c(d1)
!print *,abs(trace(d1)-trace(d1 % v)), abs(trace2(d1)-trace2c(d1))
!-----------------------------------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
