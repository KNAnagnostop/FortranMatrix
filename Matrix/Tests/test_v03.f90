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
 real   (8),allocatable :: td1(:,:,:),td2(:,:,:),td3(:,:,:),td4(:,:,:)
 complex(8),allocatable :: tm1(:,:,:),tm2(:,:,:),tm3(:,:,:),tm4(:,:,:)
 real   (8),allocatable :: ad1(:,:  ),ad2(:,:  ),ad3(:,:  ),ad4(:,:  )
 complex(8),allocatable :: am1(:,:  ),am2(:,:  ),am3(:,:  ),am4(:,:  )
 real   (8),allocatable :: vd1(:    ),vd2(:    ),vd3(:    ),vd4(:    )
 complex(8),allocatable :: vm1(:    ),vm2(:    ),vm3(:    ),vm4(:    )
 complex(8)             :: sm(3:5,6:8)
 real   (8)             :: sd(3:5,6:8)
 character(1000)        :: string
!-----------------------------------------------------------------------------------------------------------------------------------
 open(unit=f_out ,file='test.dat');  open(unit=f_math,file='test.m'); f_mout   =f_out ! output unit for matrix modules set to f_out
 call matrix_random_init
 f_mout = stdout
!-----------------------------------------------------------------------------------------------------------------------------------
 m=100;  n=m+23;  is=4;js=8            ; ie = is + m - 1; je = js + n -1
 allocate(ad1(is:ie,js:je),ad2(js:je,is:ie),ad3(is:ie,is:ie),ad4(js:je,js:je)) ! ad3 = MATMUL(ad1,ad2) ad4 = MATMUL(ad2,ad1)
 allocate(am1(is:ie,js:je),am2(js:je,is:ie),am3(is:ie,is:ie),am4(js:je,js:je)) ! am3 = MATMUL(am1,am2) am4 = MATMUL(am2,am1)
 allocate(vd1(      js:je),vd2(      is:ie),vd3(      is:ie),vd4(      js:je)) 
 allocate(vm1(      js:je),vm2(      is:ie),vm3(      is:ie),vm4(      js:je)) ! vm3 = MATMUL(am1,vm1) vm4 = MATMUL(vm2,am1)
!-----------------------------------------------------------------------------------------------------------------------------------
! MATRIX MULTIPLICATION:   lmatmul
!-----------------------------------------------------------------------------------------------------------------------------------
! am3 = LMATMUL(am1,am2);print *, norm(am3 - MATMUL(am1,am2))
! call hermitian_set(am1); am3 = LMATMUL(am1,am2,mtype='ZH'); print *, norm(am3 - MATMUL(am1,am2))
! ad3 = LMATMUL(ad1,ad2) ; print *, norm(ad3 - MATMUL(ad1,ad2))
! call symmetric_set(ad1); ad3 = LMATMUL(ad1,ad2,mtype='DS'); print *, norm(ad3 - MATMUL(ad1,ad2))
!-----------------------------------------------------------------------------------------------------------------------------------
!am3 = lmatmul(am1)         ;print *, norm(am3 - MATMUL(am1,CONJG(TRANSPOSE(am1)) ))
!am4 = lmatmul(am1,side='H');print *, norm(am4 - MATMUL(CONJG(TRANSPOSE(am1)),am1 ))
!ad3 = lmatmul(ad1)         ;print *, norm(ad3 - MATMUL(ad1,(TRANSPOSE(ad1)) ))
!ad4 = lmatmul(ad1,side='S');print *, norm(ad4 - MATMUL((TRANSPOSE(ad1)),ad1 ))
!-----------------------------------------------------------------------------------------------------------------------------------
! matrix-vector multiplication:
!call random_number(am1,sigma=PI); call random_number(vm1,sigma=PI); call random_number(vm2,sigma=PI)
!vm3 = lmatmul(A=am1,v=vm1);! vm3 = lmatmul(v=vm1,A=am1);lmatmul(am1,vm1)
!vm4 = lmatmul(v=vm2,A=am1,type='T')!lmatmul(am1,vm2,type='T');
!print *, norm(vm3 - matmul(am1,vm1)),size(vm3),lbound(vm3),ubound(vm3)
!print *, norm(vm4 - matmul(vm2,am1)),size(vm4),lbound(vm4),ubound(vm4)
!print *, norm(vm4 - matmul(TRANSPOSE(am1),vm2)),size(vm4),lbound(vm4),ubound(vm4)
!vm4 = lmatmul(v=vm2,A=am1,type='H')
!print *, norm(vm4 - matmul(vm2,CONJG(am1))),size(vm4),lbound(vm4),ubound(vm4)
!print *, norm(vm4 - matmul(hermitian(am1),vm2)),size(vm4),lbound(vm4),ubound(vm4)
!call random_number(ad1,sigma=PI); call random_number(vd1,sigma=PI); call random_number(vd2,sigma=PI)
!vd3 = lmatmul(A=ad1,v=vd1);! vd3 = lmatmul(v=vd1,A=ad1);lmatmul(ad1,vd1)
!vd4 = lmatmul(v=vd2,A=ad1,type='T')!lmatmul(ad1,vd2,type='T');
!print *, norm(vd3 - matmul(ad1,vd1)),size(vd3),lbound(vd3),ubound(vd3)
!print *, norm(vd4 - matmul(vd2,ad1)),size(vd4),lbound(vd4),ubound(vd4)
!-----------------------------------------------------------------------------------------------------------------------------------
!if(m < 7) call print(ad1,stdout,'F8.3',name='d1')!;if(m < 7) call print(am1,stdout,'F8.3',name='m1');
!if(m < 5) call   d1 % print(stdout,'F8.3')   ;   if(m > 4) call   d1 % print(stdout,'F8.3',form='i')
!if(m < 5) call   m1 % print(stdout,'F8.3')   ;   if(m > 4) call   m1 % print(stdout,'F8.3',form='i')
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
