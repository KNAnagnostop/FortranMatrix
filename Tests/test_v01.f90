program             testme
 use                matrix_mod_common
 use                matrix_mod_matrix
 use                matrix_mod_array
 implicit none
 integer   ,parameter   :: f_out = 13, f_math = 14
 integer                :: i,j
 real   (dp)            :: r,r1,r2,r3,sigma
 complex(dp)            :: z,z1,z2,z3
 type( Matrix)          :: mm, m1, m2, m3
 type(DMatrix)          :: dd, d1, d2, d3
 integer                :: m,n,is,ie,js,je
 real   (8),allocatable :: ad(:,:),ad1(:,:),ad2(:,:),ad3(:,:)
 complex(8),allocatable :: am(:,:),am1(:,:),am2(:,:),am3(:,:)
 complex(8)             :: sm(3:5,6:8)
 real   (8)             :: sd(3:5,6:8)
 character(1000)        :: string

 open(unit=f_out ,file='test.dat');  open(unit=f_math,file='test.m')
 f_mout   =f_out ! output unit for matrix modules set to f_out
 call matrix_random_init
!-----------------------------------------------------------------------------------------------------------------------------------
! mm =  Matrix('g',4)
! mm % v = 0.23e-12 * mm % v;     call mm % print(f_math,fmt='G12.4',form='Mathematica')
! m1 =  Matrix('g',4,name='m1');  call m1 % print(f_math,fmt='G12.4',form='Mathematica')
! d1 = DMatrix('g',4,name='d1');  call d1 % print(f_math,fmt='G12.4',form='Mathematica')
!-----------------------------------------------------------------------------------------------------------------------------------
!test constructors:
!
! m2 = Matrix('g',4000,sigma=PI); call m2 % save(20) !also uniform has been checked
!gnuplot> s = pi; f(x) = 1.0/sqrt(2*pi*s**2)*exp(-0.5*x**2/s**2);plot [:][0:] "< head -n 10000000 fort.20|awk '{print $2}'|hist -d 20000000 -b 0.05" u 1:4 w l,f(x) 
! m2 = Matrix(1000);call m2 % gaussian(sigma=PI);call m2 % save(unit=20) ! checked both $1,$2
! m2 = Matrix('g',1000,mtype='ZH');print *,'ZH= ',SUM(ABS(m2 % v - TRANSPOSE( CONJG(m2 % v))))/(m2 % n * m2 % m) ! check hermitian up to 1000x1000
! m2 = Matrix(4,mtype='ZH'); call m2 % random;   call m2 % print(stdout,fmt='F8.3') !returns hermitian
! m2 = Matrix(4,mtype='ZH'); call m2 % gaussian; call m2 % print(stdout,fmt='F8.3') !returns hermitian
! m1 = Matrix('u',1400,232,12,24,mtype='ZG');allocate(am,source= m1 % v);m2 = Matrix(am,is=21,js=34); !check constructor from array
! print *,SUM(ABS(m1 % v-m2 % v))/(m1 % n * m1 % m) ! check constructor by array
! call m1  % print(stdout,fmt='F8.3',form='Default',ips=14,ipe=17,jps=55,jpe=58)
! print *,m2 % is,m2 % ie,m2 % js,m2 % je; print *,m1 % is,m1 % ie,m1 % js,m1 % je
! m2 = Matrix(1254); call m2 % hermitian; print *,'ZH= ',SUM(ABS(m2 % v - TRANSPOSE( CONJG(m2 % v))))/(m2 % n * m2 % m) ! check hermitian up to 1000x1000
!-----------------------------------------------------------------------------------------------------------------------------------
! d1 = DMatrix(REAL(IMU*IMU,kind=dp),4,js=6,mtype='DS',name='DD')
! d1 = DMatrix('u',1000,mtype='DS'); ! print *,'DS= ',SUM(ABS(d1 % v - TRANSPOSE( d1 % v)))/(d1 % n * d1 % m) ! check symmetric up to 1000x1000
! call d1  % print(stdout,fmt='F8.3',form='Default')
! d1 = DMatrix('u',1400,232,12,24,mtype='DG');allocate(ad,source= d1 % v);d2 = DMatrix(ad,is=21,js=34);
! print *,SUM(ABS(d1 % v-d2 % v))/(d1 % n * d1 % m) ! check constructor by array
! call d1  % print(stdout,fmt='F8.3',form='Default',ips=14,ipe=17,jps=55,jpe=58)
! print *,d2 % is,d2 % ie,d2 % js,d2 % je;  print *,d1 % is,d1 % ie,d1 % js,d1 % je
!-----------------------------------------------------------------------------------------------------------------------------------
! d1 = DMatrix('g',1000,sigma=PI); call d1 % save(20); d1 = DMatrix(1000); call d1 % gaussian(sigma=PI);  call d1 % save(20) ! check random start - checked both uniform + gaussian
!gnuplot> s = pi; f(x) = 1.0/sqrt(2*pi*s**2)*exp(-0.5*x**2/s**2);plot [:][0:] "< head -n 10000000 fort.20|hist -d 20000000 -b 0.05" u 1:4 w l,f(x) 
! d1 = DMatrix(1021); call d1 % symmetric; print *,'DS= ',SUM(ABS(d1 % v - TRANSPOSE( d1 % v)))/(d1 % n * d1 % m) ! check symmetric up to 1000x1000
!-----------------------------------------------------------------------------------------------------------------------------------
!!  assignment and (+) operator check:
!m1 = Matrix('u',3,  name='m1',mtype='ZH',is=12,js=85); m2 = Matrix('g',4,2,name='m2')
!call m1 % print(stdout,fmt='F8.3');                                      call m2 % print(stdout,fmt='F8.3') 
!m1 % name = 'n1 new';                         m2 = m1;                   call m2 % print(stdout,fmt='F8.3')   ! Matrix           assignment - overrides previous definition
!r = PI      ;                                 m2 = r ; m2 % name = 'r' ; call m2 % print(stdout,fmt='F8.3')   ! real    constant assignment
!z = PI+2*IMU;                                 m2 = z ; m2 % name = 'z' ; call m2 % print(stdout,fmt='F8.3')   ! complex constant assignment
!allocate(am(4,5)); call array2_gauss_set(am); m2 = am; m2 % name = 'am'; call m2 % print(stdout,fmt='F8.3')   ! array2           assignment
!m1 = 1000.0_dp + m2 + 2000.0_dp                      ; m1 % name = 'p1'; call m1 % print(stdout,fmt='F14.3')  ! add a  real
!m1 = 1000.0_dp + 2000._dp*IMU + m2 +10000._dp*IMU    ; m1 % name = 'p1'; call m1 % print(stdout,fmt='F14.3')  ! add a  complex
!m1 = (1000.0_dp + am) + m2 + IMU*(2000._dp+am)       ; m1 % name = 'p1'; call m1 % print(stdout,fmt='F14.3')  ! add an array2
!-----------------------------------------------------------------------------------------------------------------------------------
!m=4;  n=m;  is=1;ie=m;js=1;je=n
!m1 = Matrix('u',m,n, name='m1',mtype='ZH',is=12,js=85);d1 = DMatrix('u',m,n, name='m1',mtype='DS',is=12,js=85);
!m2 = Matrix('g',n,m, name='m2',mtype='ZH',is=-5,js=-8);d2 = DMatrix('g',n,m, name='m2',mtype='DS',is=-5,js=-8);
!allocate(am(5:m+4,2:n+1));allocate(am1(m,n));allocate(am2(n,m));allocate(am3(m,n));
!allocate(ad(5:m+4,2:n+1));allocate(ad1(m,n));allocate(ad2(n,m));allocate(ad3(m,n));
!call random_number(am);call random_number(am1,2*PI);call random_number(am2,PI);call random_number(am3,PI);
!call random_number(r );call random_number(r1 ,2*PI);call random_number(r2 ,PI);call random_number(r3 ,PI);
!call random_number(z );call random_number(z1 ,2*PI);call random_number(z2 ,PI);call random_number(z3 ,PI);
!-----------------------------------------------------------------------------------------------------------------------------------
!!multiplication operator check:
!-----------------------------------------------------------------------------------------------------------------------------------
!am  = IMU + PI; m2 = am; m3 =  IMU * 4.0_dp * m2 + 2.0_dp * m2 + m1;m3 =   m1 / (2*PI)
!m3  = Matrix(m,n,3,4);
!m1  = Matrix('u',m,n, name='m1',mtype='ZH',is=5,js=8);
!m2  = Matrix('g',m,n, name='m2',mtype='ZH',is=5,js=8);
!m3  =  z1 *  m1       -        m1    * m2       +  m1      / z3 +        m1    *am2  / z2 +        am1* m2      * r2
!am3 =  z1 * (m1 % v)  - MATMUL(m1 % v, m2 % v)  + (m1 % v) / z3 + MATMUL(m1 % v,am2) / z2 + MATMUL(am1, m2 % v) * r2 
!call m3 % print(stdout,fmt='F8.3')
! print '(A,10I6)  ', 'size: ',size(m3 % v), shape(m3 % v),lbound(m3 % v), ubound(m3 % v)
! print '(A,G28.17)', 'dif='  ,SUM(ABS(am3 - m3 % v))/(m3 % m * m3 % n)
! print '(A,G28.17)', 'dif='  ,SUM(ABS(MATMUL(m1 % v, m2 % v) - m3 % v))/(m3 % m * m3 % n)
!-----------------------------------------------------------------------------------------------------------------------------------
!!transpose/conjg/hermitian check:
!m3 = m2 % hermitian();!m3 = m2 % transpose();!m3 = m2 % conjg(); m3 = hermitian(m2); m3 = transpose(m2); m3=conjg(m2)
!call m2 % print(stdout,fmt='F8.3')
!call m3 % print(stdout,fmt='F8.3')
!print '(A,10I6)  ', 'size m2: ',size(m2 % v), shape(m2 % v),lbound(m2 % v), ubound(m2 % v)
!print '(A,10I6)  ', 'size m3: ',size(m3 % v), shape(m3 % v),lbound(m3 % v), ubound(m3 % v)
!print *, 'dif=',SUM(ABS((conjg(transpose(m2%v))) - m3 % v))/(m3 % m * m3 % n)
!-----------------------------------------------------------------------------------------------------------------------------------
!Dmatrix tests:
!m=613;n=m;is=-20;js=1022
!d1 = DMatrix('g',m,n,is=is,js=js,mtype='DS',sigma=PI)
!d1 = DMatrix(    m,n,is=is,js=js,mtype='DS');call random_number(d1);!call d1 % random
!d1 = DMatrix(    m,n,is=is,js=js,mtype='DS');!call random_number(d1,sigma=PI);!call d1 % gaussian(sigma=PI);!call random_number(d1,sigma=PI);
!m1 =  Matrix('u',m,n,is=is,js=js,mtype='ZG',sigma=PI)
!m1 =  Matrix(    m,n,is=is,js=js,mtype='ZH');call m1 % random!call random_number(m1);!call m1 % random! the imaginary part has a peak at 0 due to real diagonal for ZH
!m1 =  Matrix(    m,n,is=is,js=js,mtype='ZH');call random_number(m1,sigma=PI)!;call m1 % gaussian(sigma=PI)!call m1 % gaussian(sigma=PI);! the imaginary part has a peak at 0 due to real diagonal for ZH
!call d1 % print(stdout,'F8.3',form='d')
!call d1 % save(20);
!print '(A,10G28.16)','diff= ',SUM(ABS( m1 % v - conjg(transpose(m1 % v))))/(m1 % m * m1 % m)
!do  i = is, is+m-1;  do j = js, js+n-1; write(21,*) d1 % v(i,j); end do; end do               ! lower diagonal
! call m1 % print(stdout,'F8.3',form='d')
! call m1 % save(20);
!do  i = is, is+m-1;  do j = js, js+m-1; write(21,'(2G28.16)') m1 % v(i,j); end do; end do     ! lower diagonal
!do i = 1, 100
! d1 = DMatrix('g',m,mtype='DS',sigma=PI)
! do j = 1, m; write(22,'(2G28.16)') d1 % v(j,j); end do                                       ! test the diagonal (special for hermitian/symmetric)
! m1 =  Matrix('g',m,mtype='ZG',sigma=PI)
! do j = 1, m; write(22,'(2G28.16)') m1 % v(j,j); end do
!end do
!gnuplot> s = pi; f(x) = 1.0/sqrt(2*pi*s**2)*exp(-0.5*x**2/s**2);plot [:][0:] "< head -n 10000000 fort.20|awk '{print $1}'|hist -d 20000000 -b 0.05" u 1:4 w l,f(x) 
!-----------------------------------------------------------------------------------------------------------------------------------
!m=4;  n=6;  is=2;js=4;
!m1 = Matrix('u',m,n, name='m1',mtype='ZH',is=is,js=js);d1 = DMatrix('u',m,n, name='m1',mtype='DS',is=is,js=js);
!m2 = Matrix('g',n,m, name='m2',mtype='ZH',is=is,js=js);d2 = DMatrix('g',n,m, name='m2',mtype='DS',is=is,js=js);
!allocate(am(5:m+4,2:n+1));allocate(am1(m,n));allocate(am2(n,m));allocate(am3(m,n));
!allocate(ad(5:m+4,2:n+1));allocate(ad1(m,n));allocate(ad2(n,m));allocate(ad3(m,n));
!call random_number(am);call random_number(am1,2*PI);call random_number(am2,PI);call random_number(am3,PI);
!call random_number(r );call random_number(r1 ,2*PI);call random_number(r2 ,PI);call random_number(r3 ,PI);
!call random_number(z );call random_number(z1 ,2*PI);call random_number(z2 ,PI);call random_number(z3 ,PI);
!-----------------------------------------------------------------------------------------------------------------------------------
!m=4;  n=4;  is=2;js=4;
!d1 = DMatrix(2*n)
!allocate(ad(is+3:is+m-1+3,js-5:js+n-1-5)); call random_number(ad)
!d1 = DMatrix(ad,is=is,js=js,name='d1') ! get arbitrary bounds
!d1 = DMatrix(ad, is=lbound(ad,1),js=lbound(ad,2),mtype='DS') ! get bounds from ad
!call array2_print_d(ad,stdout,'F8.3',name='ad'); call d1 % print(stdout,'F8.3')
!-----------------------------------------------------------------------------------
!m1 = Matrix(2*n)
!allocate(am(is+3:is+m-1+3,js-5:js+n-1-5)); call random_number(am)
!m1 = Matrix(am,name='m1') ! get is=js=1
!m1 = Matrix(am,is=is,js=js,name='m1') ! get arbitrary bounds
!m1 = Matrix(am, is=lbound(am,1),js=lbound(am,2),mtype='ZG') ! get bounds from ad
!m1 = Matrix(am, is=lbound(am,1),mtype='ZG') ! get bounds from ad
!call array2_print(am,stdout,'F8.3',name='am'); call m1 % print(stdout,'F8.3')
!-----------------------------------------------------------------------------------
!m=6;  n=m;  is=2;js=4;
!d1 = DMatrix('u',m,n,is,js)
!call d1 % print(stdout,'F8.3')
!call d1 % save(20);close(20)
!d1 = DMatrix('u',m,n,is,js)
!call d1 % print(stdout,'F8.3')
!open(unit=20,file='fort.20'); call d1 % read(20)
!call d1 % print(stdout,'F8.3')
!m1 = Matrix('u',m,n,is,js)
!call m1 % print(stdout,'F8.3')
!call m1 % save(20);close(20)
!m1 = Matrix('u',m,n,is,js)
!call m1 % print(stdout,'F8.3')
!open(unit=20,file='fort.20'); call m1 % read(20)
!call m1 % print(stdout,'F8.3')
!m=6;  n=m-1;  is=2;js=4;
!d1 = DMatrix('u',m,n,is,js); call d1 % print(stdout,'F8.3')!; d2 = d1 % symmetric()
!d2 = d1 % transpose()      ; call d2 % print(stdout,'F8.3')
! call d1 % symmetric_set; print  '(A,10G28.16)','diff= ',SUM(ABS( d1 % v - (transpose(d1 % v))))/(d1 % m * d1 % m)
!-----------------------------------------------------------------------------------
! Check DMatrix matmul, random_number
!m=33;  n=m+2;  is=2;js=4; 
!d1 = DMatrix('g',m,n,is=is  ,js=js  ,mtype='DG',name='d1'); allocate(ad1(is:is+m-1,js:js+n-1)); ad1 = d1 % v
!d2 = DMatrix('g',n,m,js=js-6,is=is+5,mtype='DG',name='d2'); allocate(ad2(js:js+n-1,is:is+m-1)); ad2 = d2 % v
!r  = PI; d1 = r; ad = 2*r; d2 = -d1 - ad;                       
!d1 = DMatrix('g',m,is=is,js=js,mtype='DS',name='d1');     ad1 = d1 % v ; d2 = d1 / PI
!d3 =         d2 * d1  ;d3 =         ad2 * d1  ; ad3 = MATMUL(ad2 ,ad1) ; print *,norm(d3 - ad3)
! m=544;  n=m;  is=2;js=4;
!d1 = DMatrix(m,n,is=is  ,js=js  ,mtype='DS',name='d1')
!call random_number(d1,sigma=PI)
!call d1 % print(stdout,'F8.3');
!call d1 % save(20)
!gnuplot> s = pi; f(x) = 1.0/sqrt(2*pi*s**2)*exp(-0.5*x**2/s**2);plot [:][0:] "< head -n 10000000 fort.20|awk '{print $1}'|hist -d 20000000 -b 0.05" u 1:4 w l,f(x)
!m=343;  n=m+2;  is=2;js=4; 
!d1 = DMatrix('g',m,n,is=is  ,js=js  ,mtype='DG',name='d1'); d2 = d1 % transpose();print *,norm(d2 % v - transpose(d1 % v)) 
!call d1 % print(stdout,'F8.3');d2 % name = 'd2';call d2 % print(stdout,'F8.3');
!-----------------------------------------------------------------------------------
!Combined operations:
!m=4;  n=m;  is=2;js=4;
!call random_number(r );call random_number(r1 ,2*PI);call random_number(r2 ,PI);call random_number(r3 ,PI);
!call random_number(z );call random_number(z1 ,2*PI);call random_number(z2 ,PI);call random_number(z3 ,PI);
!d1  = DMatrix('g',m,n,is=is  ,js=js  ,mtype='DS',name='d1'); allocate(ad1(is:is+m-1,js:js+n-1)); ad1 = d1 % v
!m1  =  Matrix('g',n,m,is=is  ,js=js  ,mtype='ZH',name='m1'); allocate(am1(is:is+m-1,js:js+n-1)); am1 = m1 % v
!d2  = DMatrix('g',n,m,js=js-6,is=is+5,mtype='DG',name='d2'); allocate(ad2(js:js+n-1,is:is+m-1)); ad2 = d2 % v
!d3  = d1 / r3 * d2 * r2          - d2 * r3 * r1 *z d1 + d2/r1 ;
!ad3 = (r2/r3) *  MATMUL(ad1,ad2)   - (r3*r1) * MATMUL(ad2,ad1) + ad2/r1 + IMU * ad2
!-----------------------------------------------------------------------------------
!m3  =         d1 *  m1;am3 = MATMUL(ad1 , am1)
!m3  =  d1 -  m1; am3 = ad1 - am1
!if(m < 5) call   m3 % print(stdout,'F8.3')   ;   if(m > 4) call   m3 % print(stdout,'F8.3',form='i')
!print *,norm(m3 - am3) , norm(m3 - m3 % hermitian() ), m3 % mtype
!-----------------------------------------------------------------------------------
m=4;  n=m;  is=2;js=4;
m1  =  Matrix('g',n,m,is=is  ,js=js  ,mtype='ZG',name='m1'); allocate(am1(is:is+m-1,js:js+n-1))     ; am1 = m1 % v
d1  = DMatrix('g',n,m,is=is  ,js=js  ,mtype='DG',name='d1'); allocate(ad1(is:is+m-1,js:js+n-1))     ; ad1 = d1 % v
d2  = DMatrix('g',n,m,is=is  ,js=js  ,mtype='DG',name='d2'); allocate(ad2(is:is+m-1,js:js+n-1))     ; ad2 = d2 % v

call m1 % hermitian_set

!if(m < 5) call   d1 % print(stdout,'F8.3')   ;   if(m > 4) call   d1 % print(stdout,'F8.3',form='i')
!if(m < 5) call   d2 % print(stdout,'F8.3')   ;   if(m > 4) call   d2 % print(stdout,'F8.3',form='i')
 if(m < 5) call   m1 % print(stdout,'F8.3')   ;   if(m > 4) call   m1 % print(stdout,'F8.3',form='i')
!print *, norm(m1- z * ad1)


!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
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
