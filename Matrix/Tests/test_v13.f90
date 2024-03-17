!gfortran -ffree-line-length-0 -I ../lib/gnu/ test_v13.f90 -L ../lib/gnu/ -lMatrix -llapack -lblas -o t
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
 integer   ,allocatable :: pm1(:  ),pm2(:  ),pm3(:  ),pm4(:  )
 complex(8)             :: sm(3:5,6:8)
 real   (8)             :: sd(3:5,6:8)
 real   (dp)            :: dtr
 complex(dp)            :: ctr
 character(1000)        :: string
!-----------------------------------------------------------------------------------------------------------------------------------
 complex(8),allocatable :: zw(:,:), zev(:),zevec(:,:)
 real   (8),allocatable :: dw(:,:), dev(:),devec(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
 open(unit=f_out ,file='test.dat');  open(unit=f_math,file='test.m'); f_mout   =f_out ! output unit for matrix modules set to f_out
 call matrix_random_init
!-----------------------------------------------------------------------------------------------------------------------------------
 m=512;  n=m;  is=1;js=1            ; ie = is + m - 1; je = js + n -1
 allocate(ad1(is:ie,js:je),ad2(is:ie,js:je),ad3(is:ie,js:je),ad4(is:ie,js:je))
 allocate(am1(is:ie,js:je),am2(is:ie,js:je),am3(is:ie,js:je),am4(is:ie,js:je))
 allocate(vd1(      js:je),vd2(      js:je),vd3(      js:je),vd4(      js:je))
 allocate(vm1(      js:je),vm2(      js:je),vm3(      js:je),vm4(      js:je))
 allocate(pm1(      js:je),pm2(      js:je),pm3(      js:je),pm4(      js:je))
!-----------------------------------------------------------------------------------------------------------------------------------
!Test trace2(am1,am2) and trace2(ad1,ad2)
 call random_number(am1,sigma=PI); call random_number(am2,sigma=PI); 
 call random_number(ad1,sigma=PI); call random_number(ad2,sigma=PI); 

 allocate(zev(m),zevec(m,m)) ;  allocate(dev(m),devec(m,m))


 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(am1)             ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','ZEV_US: eval: ',zev
 print '(A,G28.16  )','ZEV_US: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(am1,sortby='RF') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','ZEV_RF: eval: ',zev
 print '(A,G28.16  )','ZEV_RF: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(am1,sortby='RR') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','ZEV_RR: eval: ',zev
 print '(A,G28.16  )','ZEV_RR: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(am1,sortby='IF') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','ZEV_IF: eval: ',zev
 print '(A,G28.16  )','ZEV_IF: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(am1,sortby='IR') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','ZEV_IR: eval: ',zev
 print '(A,G28.16  )','ZEV_IR: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(am1,sortby='MF') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','ZEV_MF: eval: ',zev
 print '(A,G28.16  )','ZEV_MF: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(am1,sortby='MR') ; zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','ZEV_MR: eval: ',zev
 print '(A,G28.16  )','ZEV_MR: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

!-----------------------------------------------------------------------------------------------------------------------------------
! Double arrays
!-----------------------------------------------------------------------------------------------------------------------------------
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 !a general real matrix has complex eigenvalues:
 zw = eigenvectors(ad1)             ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','DEV_US: eval: ',zev
 print '(A,G28.16  )','DEV_US: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(ad1,sortby='RF') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','DEV_RF: eval: ',zev
 print '(A,G28.16  )','DEV_RF: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(ad1,sortby='RR') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','DEV_RR: eval: ',zev
 print '(A,G28.16  )','DEV_RR: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(ad1,sortby='IF') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','DEV_IF: eval: ',zev
 print '(A,G28.16  )','DEV_IF: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(ad1,sortby='IR') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','DEV_IR: eval: ',zev
 print '(A,G28.16  )','DEV_IR: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(ad1,sortby='MF') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','DEV_MF: eval: ',zev
 print '(A,G28.16  )','DEV_MF: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw = eigenvectors(ad1,sortby='MR') ;  zev = zw(:,1); zevec = zw(:,2:)
 print '(A,1000F8.3)','DEV_MR: eval: ',zev
 print '(A,G28.16  )','DEV_MR: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

!-----------------------------------------------------------------------------------------------------------------------------------
! Hermitian arrays
!-----------------------------------------------------------------------------------------------------------------------------------
 call hermitian_set(am1)
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  = eigenvectors(am1,mtype='ZH')             ;  zev = zw(:,1); zevec = zw(:,2:)
 dev = eigenvalues (am1,mtype='ZH')             
 print '(A,1000F8.3)','HEV_US: eval: ',zev
 print '(A,1000F8.3)','HEV_US: eval: ',dev
 print '(A,G28.16  )','HEV_US: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','HEV_US: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - dev(i) * zevec(:,i)),i=1,m)])
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  =      eigenvectors(am1,mtype='ZH',sortby='RF') ;  zev = zw(:,1); zevec = zw(:,2:)
 dev = sort(eigenvalues (am1,mtype='ZH'),   by='RF')
 print '(A,1000F8.3)','HEV_RF: eval: ',zev
 print '(A,1000F8.3)','HEV_RF: eval: ',dev
 print '(A,G28.16  )','HEV_RF: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','HEV_RF: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - dev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  =      eigenvectors(am1,mtype='ZH',sortby='RR') ;  zev = zw(:,1); zevec = zw(:,2:)
 dev = sort(eigenvalues (am1,mtype='ZH'),   by='RR')
 print '(A,1000F8.3)','HEV_RR: eval: ',zev
 print '(A,1000F8.3)','HEV_RR: eval: ',dev
 print '(A,G28.16  )','HEV_RR: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','HEV_RR: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - dev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  =      eigenvectors(am1,mtype='ZH',sortby='MF') ;  zev = zw(:,1); zevec = zw(:,2:)
 dev = sort(eigenvalues (am1,mtype='ZH'),   by='MF')
 print '(A,1000F8.3)','HEV_MF: eval: ',zev
 print '(A,1000F8.3)','HEV_MF: eval: ',dev
 print '(A,G28.16  )','HEV_MF: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','HEV_MF: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - dev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  =      eigenvectors(am1,mtype='ZH',sortby='MR') ;  zev = zw(:,1); zevec = zw(:,2:)
 dev = sort(eigenvalues (am1,mtype='ZH'),   by='MR')
 print '(A,1000F8.3)','HEV_MR: eval: ',zev
 print '(A,1000F8.3)','HEV_MR: eval: ',dev
 print '(A,G28.16  )','HEV_MR: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','HEV_MR: norm: ',norm([(norm( matmul(am1,zevec(:,i))  - dev(i) * zevec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')

!-----------------------------------------------------------------------------------------------------------------------------------
! Symmetric arrays
!-----------------------------------------------------------------------------------------------------------------------------------
 call symmetric_set(ad1)
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  = eigenvectors(ad1,mtype='DS')             ;  zev = zw(:,1); zevec = zw(:,2:)
 dw  = eigenvectors(ad1,mtype='DS')             ;  dev = dw(:,1); devec = dw(:,2:)
 print '(A,1000F8.3)','SEV_US: eval: ',zev
 print '(A,1000F8.3)','SEV_US: eval: ',dev
 print '(A,G28.16  )','SEV_US: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','SEV_US: norm: ',norm([(norm( matmul(ad1,devec(:,i))  - dev(i) * devec(:,i)),i=1,m)])
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  = eigenvectors(ad1,mtype='DS',sortby='RF') ;  zev = zw(:,1); zevec = zw(:,2:)
 dw  = eigenvectors(ad1,mtype='DS',sortby='RF') ;  dev = dw(:,1); devec = dw(:,2:)
 print '(A,1000F8.3)','SEV_RF: eval: ',zev
 print '(A,1000F8.3)','SEV_RF: eval: ',dev
 print '(A,G28.16  )','SEV_RF: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','SEV_RF: norm: ',norm([(norm( matmul(ad1,devec(:,i))  - dev(i) * devec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
!call print(devec,unit=6,fmt='F8.3',name='devec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  = eigenvectors(ad1,mtype='DS',sortby='RR') ;  zev = zw(:,1); zevec = zw(:,2:)
 dw  = eigenvectors(ad1,mtype='DS',sortby='RR') ;  dev = dw(:,1); devec = dw(:,2:)
 print '(A,1000F8.3)','SEV_RR: eval: ',zev
 print '(A,1000F8.3)','SEV_RR: eval: ',dev
 print '(A,G28.16  )','SEV_RR: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','SEV_RR: norm: ',norm([(norm( matmul(ad1,devec(:,i))  - dev(i) * devec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
!call print(devec,unit=6,fmt='F8.3',name='devec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  = eigenvectors(ad1,mtype='DS',sortby='MF') ;  zev = zw(:,1); zevec = zw(:,2:)
 dw  = eigenvectors(ad1,mtype='DS',sortby='MF') ;  dev = dw(:,1); devec = dw(:,2:)
 print '(A,1000F8.3)','SEV_MF: eval: ',zev
 print '(A,1000F8.3)','SEV_MF: eval: ',dev
 print '(A,G28.16  )','SEV_MF: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','SEV_MF: norm: ',norm([(norm( matmul(ad1,devec(:,i))  - dev(i) * devec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
!call print(devec,unit=6,fmt='F8.3',name='devec')
 print '(A)','---------------------------------------------------------------------------------------------------------------------'
 zw  = eigenvectors(ad1,mtype='DS',sortby='MR') ;  zev = zw(:,1); zevec = zw(:,2:)
 dw  = eigenvectors(ad1,mtype='DS',sortby='MR') ;  dev = dw(:,1); devec = dw(:,2:)
 print '(A,1000F8.3)','SEV_MR: eval: ',zev
 print '(A,1000F8.3)','SEV_MR: eval: ',dev
 print '(A,G28.16  )','SEV_MR: norm: ',norm([(norm( matmul(ad1,zevec(:,i))  - zev(i) * zevec(:,i)),i=1,m)])
 print '(A,G28.16  )','SEV_MR: norm: ',norm([(norm( matmul(ad1,devec(:,i))  - dev(i) * devec(:,i)),i=1,m)])
!call print(zevec,unit=6,fmt='F8.3',name='zevec')
!call print(devec,unit=6,fmt='F8.3',name='devec')
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
