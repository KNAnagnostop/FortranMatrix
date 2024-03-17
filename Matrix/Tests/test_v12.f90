!gfortran -ffree-line-length-0 -I ../lib/gnu/ test_v12.f90 -L ../lib/gnu/ -lMatrix -llapack -lblas -o t
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
 open(unit=f_out ,file='test.dat');  open(unit=f_math,file='test.m'); f_mout   =f_out ! output unit for matrix modules set to f_out
 call matrix_random_init
!-----------------------------------------------------------------------------------------------------------------------------------
 m=256;  n=m;  is=1;js=1            ; ie = is + m - 1; je = js + n -1
 allocate(ad1(is:ie,js:je),ad2(is:ie,js:je),ad3(is:ie,js:je),ad4(is:ie,js:je))
 allocate(am1(is:ie,js:je),am2(is:ie,js:je),am3(is:ie,js:je),am4(is:ie,js:je))
 allocate(vd1(      js:je),vd2(      js:je),vd3(      js:je),vd4(      js:je))
 allocate(vm1(      js:je),vm2(      js:je),vm3(      js:je),vm4(      js:je))
 allocate(pm1(      js:je),pm2(      js:je),pm3(      js:je),pm4(      js:je))
!-----------------------------------------------------------------------------------------------------------------------------------
!Test trace2(am1,am2) and trace2(ad1,ad2)
 call random_number(vm1,sigma=PI); call random_number(vm2,sigma=PI); 
 call random_number(vd1,sigma=PI); call random_number(vd2,sigma=PI); 

if( .true.)then

 if( m>= 4) vm1(4) = vm1(2 ) ! create degeneracy
 if( m>= 8) vm1(1) = vm1(2 ) ! create degeneracy
 if( m>=12) vm1(9) = vm1(11) ! create degeneracy
 if( m>=16) vm1(3) = vm1(15) 
 if( m>=64) then ;vm1(31) = -vm1(28); vm1(48) = -vm1(21); vm1(41) = -vm1(21);vm1(54) = vm1(2); vm1(61) = vm1(2); vm1(1) = vm1(2);end if 

 if( m>= 4) vd1(4) = vd1(2 ) ! create degeneracy
 if( m>= 8) vd1(1) = vd1(2 ) ! create degeneracy
 if( m>=12) vd1(9) = vd1(11) ! create degeneracy
 if( m>=16) vd1(3) = vd1(15) 
 if( m>=64) then ;vd1(31) = -vd1(28); vd1(48) = -vd1(21); vd1(41) = -vd1(21);vd1(54) = vd1(2); vd1(61) = vd1(2); vd1(1) = vd1(2);end if 

end if

!-------------------------------------------------------------------
!Complex sorting
 print '(A          )','--------------------------------------------'
 print '(A          )','............................................'
 print '(A,1000F10.4)','US: ',            vm1
 print '(A          )','............................................'
 print '(A,1000F10.4)','US: ', abs  (     vm1 )
 print '(A          )','............................................'
                                          vm2 = sort(vm1,P=pm1,by='MF')
 print '(A,1000F10.4)','MF: ', abs  (     vm2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vm3(i) = vm1(pm1(i)) ; end do
 print '(A,G28.16   )','ZMF norm= ',norm(vm2-vm3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
                                          vm2 = sort(vm1,P=pm1,by='MR')
 print '(A,1000F10.4)','MR: ', abs  (     vm2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vm3(i) = vm1(pm1(i)) ; end do
 print '(A,G28.16   )','ZMR norm= ',norm(vm2-vm3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
 print '(A,1000F10.4)','US: ',       dble(vm1)
 print '(A          )','............................................'
                                          vm2 = sort(vm1,P=pm1,by='RF')
 print '(A,1000F10.4)','RF: ', dble (     vm2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vm3(i) = vm1(pm1(i)) ; end do
 print '(A,G28.16   )','ZRF norm= ',norm(vm2-vm3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
                                          vm2 = sort(vm1,P=pm1,by='RR')
 print '(A,1000F10.4)','RR: ', dble (     vm2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vm3(i) = vm1(pm1(i)) ; end do
 print '(A,G28.16   )','ZRR norm= ',norm(vm2-vm3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
 print '(A,1000F10.4)','US: ',      aimag(vm1)
 print '(A          )','............................................'
                                          vm2 = sort(vm1,P=pm1,by='IF')
 print '(A,1000F10.4)','IF: ', aimag(     vm2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vm3(i) = vm1(pm1(i)) ; end do
 print '(A,G28.16   )','ZIF norm= ',norm(vm2-vm3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
                                          vm2 = sort(vm1,P=pm1,by='IR')
 print '(A,1000F10.4)','IR: ', aimag(     vm2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vm3(i) = vm1(pm1(i)) ; end do
 print '(A,G28.16   )','ZIR norm= ',norm(vm2-vm3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
!-------------------------------------------------------------------
!Double sorting
 print '(A          )','--------------------------------------------'
 print '(A          )','--------------------------------------------'
 print '(A,1000F10.4)','US: ',            vd1
 print '(A          )','............................................'
                                          vd2 = sort(vd1,P=pm1,by='RF')
 print '(A,1000F10.4)','RF: ', dble (     vd2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vd3(i) = vd1(pm1(i)) ; end do
 print '(A,G28.16   )','DRF norm= ',norm(vd2-vd3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
                                          vd2 = sort(vd1,P=pm1,by='RR')
 print '(A,1000F10.4)','RR: ', dble (     vd2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vd3(i) = vd1(pm1(i)) ; end do
 print '(A,G28.16   )','DRR norm= ',norm(vd2-vd3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
 print '(A,1000F10.4)','US: ',        abs(vd1)
 print '(A          )','............................................'
                                          vd2 = sort(vd1,P=pm1,by='AF')
 print '(A,1000F10.4)','AF: ', abs  (     vd2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vd3(i) = vd1(pm1(i)) ; end do
 print '(A,G28.16   )','DAF norm= ',norm(vd2-vd3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
                                          vd2 = sort(vd1,P=pm1,by='AR')
 print '(A,1000F10.4)','AR: ', abs  (     vd2 )
 print '(A,1000I10  )','PO: ',                             pm1
 do i = 1, m; vd3(i) = vd1(pm1(i)) ; end do
 print '(A,G28.16   )','DAF norm= ',norm(vd2-vd3)
 do i = 1, m; do j = 1, m; if(i /= j .and. pm1(i) == pm1(j) ) print '(A,4I6)','DBG: oops same indices',i,j,pm1(i),pm1(j)  ; end do ;end do;
 print '(A          )','............................................'
!-------------------------------------------------------------------
 
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
