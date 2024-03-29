!gfortran -ffree-line-length-0 -I ../lib/gnu/ test_v14.f90 -L ../lib/gnu/ -lMatrix -llapack -lblas -o t
program             testme
 use                array_mod
 implicit none
 integer   ,parameter   :: dp = 8
 integer   ,parameter   :: f_out = 13, f_math = 14
 integer                :: i,j
 real   (dp)            :: r,r1,r2,r3,sigma
 complex(dp)            :: z,z1,z2,z3
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
 m=5;  n=m;  is=1;js=1            ; ie = is + m - 1; je = js + n -1
 allocate(ad1(is:ie,js:je),ad2(is:ie,js:je),ad3(is:ie,js:je),ad4(is:ie,js:je))
 allocate(am1(is:ie,js:je),am2(is:ie,js:je),am3(is:ie,js:je),am4(is:ie,js:je))
 allocate(vd1(      js:je),vd2(      js:je),vd3(      js:je),vd4(      js:je))
 allocate(vm1(      js:je),vm2(      js:je),vm3(      js:je),vm4(      js:je))
 allocate(pm1(      js:je),pm2(      js:je),pm3(      js:je),pm4(      js:je))
!-----------------------------------------------------------------------------------------------------------------------------------
!Test trace2(am1,am2) and trace2(ad1,ad2)

 pm1 = [(i,i=is,ie)]
 vd1 = [(i,i=is,ie)]
 vm1 = [(i,i=is,ie)]

 !print '(1000I6  )',random_sort(m)
 !print '(1000I6  )',random_sort(is,ie)
 !print '(1000I6  )',random_sort(pm1)
 !print '(1000F6.2)',random_sort(vd1)
 !print '(1000F6.2)',random_sort(vm1)
 !print '(A       )',random_sort('Superstrings! E8 x E8')

 call random_sort_array(pm1); print '(1000I6  )',pm1
 call random_sort_array(vd1); print '(1000F6.2)',vd1
 call random_sort_array(vm1); print '(1000F6.2)',vm1

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
