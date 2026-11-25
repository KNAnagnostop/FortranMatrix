program             testme
 use                matrix_mod_common
 use                matrix_mod_matrix
 use                matrix_mod_array
 implicit none
 integer   ,parameter   :: f_out = 13, f_math = 14
 integer                :: i,j
 real   (dp)            :: r,r1,r2,r3,sigma
 complex(dp)            :: z,z1,z2,z3
 type( Matrix)          :: m1, m2, m3
 type(DMatrix)          :: d1, d2, d3
 integer                :: m,n,is,ie,js,je
 real   (8),allocatable :: ad1(:,:),ad2(:,:),ad3(:,:),ad4(:,:)
 complex(8),allocatable :: am1(:,:),am2(:,:),am3(:,:),am4(:,:)
 real   (8),allocatable :: vd1(:  ),vd2(:  ),vd3(:  ),vd4(:  )
 complex(8),allocatable :: vm1(:  ),vm2(:  ),vm3(:  ),vm4(:  )
 complex(8)             :: sm(3:5,6:8)
 real   (8)             :: sd(3:5,6:8)
 real   (dp)            :: dtr
 complex(dp)            :: ctr
 character(1000)        :: string
!-----------------------------------------------------------------------------------------------------------------------------------
 open(unit=f_out ,file='test.dat');  open(unit=f_math,file='test.m'); f_mout   =f_out ! output unit for matrix modules set to f_out
 call matrix_random_init
!-----------------------------------------------------------------------------------------------------------------------------------
 m=96;  n=m;  is=4;js=8            ; ie = is + m - 1; je = js + n -1
 allocate(ad1(is:ie,js:je),ad2(is:ie,js:je),ad3(is:ie,js:je),ad4(is:ie,js:je))
 allocate(am1(is:ie,js:je),am2(is:ie,js:je),am3(is:ie,js:je),am4(is:ie,js:je))
 allocate(vd1(      js:je),vd2(      js:je),vd3(      js:je),vd4(      js:je))
 allocate(vm1(      js:je),vm2(      js:je),vm3(      js:je),vm4(      js:je))
!-----------------------------------------------------------------------------------------------------------------------------------
!Test trace2(am1,am2) and trace2(ad1,ad2)
 call random_number(am1,sigma=PI); call random_number(am2,sigma=PI); 
 call random_number(ad1,sigma=PI); call random_number(ad2,sigma=PI); 
 am3 = matmul(am1,am2); ad3 = matmul(ad1,ad2);
 print *,'G complextrace: ',trace2(am1,am2)-trace(am3)
 print *,'G    realtrace: ',trace2(ad1,ad2)-trace(ad3)
! Test Hermitian calculation:
 call hermitian_set(am2);call symmetric_set(ad2); 
 am3 = matmul(am1,am2); ad3 = matmul(ad1,ad2);
 print *,'S complextrace: ',trace2(am1,am2,mtype='ZH')-trace(am3)
 print *,'S    realtrace: ',trace2(ad1,ad2,mtype='DS')-trace(ad3)
!-----------------------------------------------------------------------------------------------------------------------------------
!call random_number(am1,sigma=PI)
!am2  = hermitian(am1);! call hermitian_set(am1);print *,norm(am1-am2),norm(am1 - CONJG(TRANSPOSE(am1)))
!call hermitian_set(am1);!call traceless(am1)
!if(m < 7) call print(am1,stdout,'F8.3',name='m1'); !if(m < 7) call print(am2,stdout,'F8.3',name='m2')
!z=ZERO;do i = 0, m-1; z = z + am1(is+i,js+i);end do;print *,trace(am1)-z
!z=ZERO;do i = 0, m-1;do j = 0, n-1; z = z + am1(is+i,js+j) * am1(is+j,js+i);end do ;end do; print *,trace2(am1,mtype='ZH') -z
!print *, diagonal(am1)
!z1 = trace(am1); vm1 = diagonal(am1)-z1/m; call traceless(am1); z2 = trace(am1); vm2 = diagonal(am1);print *,z2,SUM(ABS(vm1-vm2))/m
!call hermitian_set(am1); am2 = am1; call traceless(am2);print *,trace2c(am1)-trace2(am2)
!z = ZERO; do i = 0, m-1;do j = 0, n-1;z = z + am1(is+i,js+j) * am1(is+j,js+i);end do;end do;print *,trace2(am1)-z
!z=trace(am1);print *,(trace2c(am1)-(trace2(am1)-z*z/m))/m
!z=trace(am2);print *,(trace2c(am2,mtype='ZH')-(trace2(am2,mtype='ZH')-z*z/m))/m
!if(m < 7) call print(am1,stdout,'F8.3',name='m1'); !if(m < 7) call print(am2,stdout,'F8.3',name='m2')
!-----------------------------------------------------------------------------------------------------------------------------------
!call random_number(ad1,sigma=PI); ad2 = (ad1+TRANSPOSE(ad1))/2 !ad2 is symmetric
!print *,diagonal(ad1)
!r = 0.0_dp;do i = 0, m-1; r = r + ad1(is+i,js+i);end do;print *,trace(ad1)-r,SUM(diagonal(ad1))-r
!print *, trace2(ad2,mtype='DS')-trace2(ad2)
!r=0.0_dp;do i = 0, m-1;               r = r + ad2(is+i,js+i)                 ;end do ;print *,trace(ad2)-r
!r=0.0_dp;do i = 0, m-1;do j = 0, n-1; r = r + ad1(is+i,js+j) * ad1(is+j,js+i);end do ;end do; print *,trace2(ad1) -r
!r=0.0_dp;do i = 0, m-1;do j = 0, n-1; r = r + ad2(is+i,js+j) * ad2(is+j,js+i);end do ;end do; print *,trace2 (ad2,mtype='DS') -r
!r=trace(ad1); print *,(trace2c(ad1) - trace2(ad1) + r*r/m)/m
!call traceless(ad2);call traceless(ad1)
!r=0.0_dp;do i = 0, m-1;do j = 0, n-1; r = r + ad2(is+i,js+j) * ad2(is+j,js+i);end do ;end do; print *,trace2c(ad2,mtype='DS') -r
!r=0.0_dp;do i = 0, m-1;do j = 0, n-1; r = r + ad1(is+i,js+j) * ad1(is+j,js+i);end do ;end do; print *,trace2c(ad1           ) -r
!-----------------------------------------------------------------------------------------------------------------------------------
!call random_number(ad1,sigma=PI); ad2 = symmetric(ad1); vd1=diagonal(ad1); vd2=diagonal(ad2);ad3=diagonal(vd1);ad4=diagonal(vd2)
!if(m < 7) call print(ad1,stdout,'F8.3',name='d1');if(m < 7) call print(ad2,stdout,'F8.3',name='d2');
!if(m < 7) call print(ad3,stdout,'F8.3',name='v1');if(m < 7) call print(ad4,stdout,'F8.3',name='v2')
!print *,norm(ad3-ad4)
!-----------------------------------------------------------------------------------------------------------------------------------
!deallocate(ad1)
!ad1 = diagonal(1.0_dp,n); am1 = diagonal(IMU,n)
!call random_number(ad2);call random_number(am2);
!ad1 = diagonal_matrix(ad2)-diagonal(diagonal(ad2));am1 = diagonal_matrix(am2)-diagonal(diagonal(am2));print *,norm(ad1),norm(am1)
!-----------------------------------------------------------------------------------------------------------------------------------
!call random_number(ad1);call random_number(am1); call random_number(sm); call random_number(sd)
!call print  (ad1,stdout,'F8.3',name='d1')
!call print  (ad1,stdout,'F8.3',name='d1',ips=5,ipe=6,jps=12,jpe=13)
!call printna(sd ,stdout,'F8.3',name='d1')
!call printna(sd ,stdout,'F8.3',name='d1',ips=2,jps=2)
!call print  (am1,stdout,'F8.3',name='m1');call print(ad1,stdout,'F8.3',name='d1')
!call save   (am1,20);call save(ad1,21)
!close(20);close(21)
!open(unit=20,file='fort.20');open(unit=21,file='fort.21'); 
!deallocate(am2);allocate(am2(222:222+m-1,627:627+n-1))
!call read(am2,20);call read(ad2,21)
!call print(am2,stdout,'F8.3',name='m2');call print(ad2,stdout,'F8.3',name='d2')
!print *,norm(am2-am1),norm(ad2-ad1)
!-----------------------------------------------------------------------------------------------------------------------------------
!call random_number(ad1,sigma=PI); ad2=symmetric(ad1); print *,isSymmetric(ad1),isSymmetric(ad2)
!call random_number(am1,sigma=PI); am2=hermitian(am1); am3 = 0.5_dp * (am1 + TRANSPOSE(am1))
!print *, isSymmetric(am1),isSymmetric(am3),isHermitian(am2)

!if(m < 7) call print(ad1,stdout,'F8.3',name='d1')!;if(m < 7) call print(am1,stdout,'F8.3',name='m1');

!if(m < 5) call   d1 % print(stdout,'F8.3')   ;   if(m > 4) call   d1 % print(stdout,'F8.3',form='i')
!if(m < 5) call   d2 % print(stdout,'F8.3')   ;   if(m > 4) call   d2 % print(stdout,'F8.3',form='i')
!if(m < 5) call   m1 % print(stdout,'F8.3')   ;   if(m > 4) call   m1 % print(stdout,'F8.3',form='i')
!print *, norm(m1- z * ad1)
!-----------------------------------------------------------------------------------------------------------------------------------
!r = nan();r1 = nan(r1);z = nan();z1 = nan(z);  z2=CMPLX(nan(r1),0.0_dp); z3=CMPLX(0.0_dp,nan(r1));
!print *,r,r1,isNaN(r),isNaN(r1); print *,z,z1,isNaN(z),isNaN(z1);print *,z2,z3,isNaN(z2),isNaN(z3)
!r=0;   print *,r,isNaN(r); z=IMU; print *,z,isNaN(z); z=PI;  print *,z,isNaN(z); z=PI+IMU;print *,z,isNaN(z)
!ad1 = nan(); ad1(5,9) = nan(r);print *, isNaN(ad1)
!am1 = nan(); am1 = cmplx(PI,nan(r),kind=dp);am1 = nan(r);
!am1(5,9) = nan();am1(5,9) = cmplx(PI,nan(r),kind=dp);am1(5,9) = nan(r); print *, isNaN(am1)
!-----------------------------------------------------------------------------------------------------------------------------------
!deallocate(tm2);allocate(tm2,mold=tm1);if(m<7) call print(tm1,stdout,fmt='G15.3',name='v1')
!open(unit=20,file="test.con");call save(tm1,unit=20); close(20); open(unit=20,file="test.con"); call read(tm2,unit=20);print *,norm(tm2-tm1);close(20)
!vm1(js) = nan() !nan() !cmplx(nan(r),0.0_dp) !cmplx(0.0_dp,nan(r)),print *,isNaN(vm1),vm1(js),vm1(je) !also tested tm1,td1
!-----------------------------------------------------------------------------------------------------------------------------------
! allocate(tm1(is:ie,js:je,js:je));allocate(td1(is:ie,js:je,js:je));call random_number(tm1);call random_number(td1)
! deallocate(vm1);allocate(vm1(1000000));deallocate(vd1);allocate(vd1(1000000));call random_number(vm1);call random_number(vd1);call save(vd1,unit=20)
!gnuplot> s = pi; f(x) = 1.0/sqrt(2*pi*s**2)*exp(-0.5*x**2/s**2);plot [:][0:] "< head -n 10000000 fort.20|awk '{print $1}'|hist -d 20000000 -b 0.05" u 1:4 w l,f(x)
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
