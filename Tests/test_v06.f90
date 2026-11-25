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
 integer                :: m,n,is,ie,js,je
 real   (8),allocatable :: td1(:,:,:),td2(:,:,:),td3(:,:,:),td4(:,:,:)
 complex(8),allocatable :: tm1(:,:,:),tm2(:,:,:),tm3(:,:,:),tm4(:,:,:)
 real   (8),allocatable :: ad1(:,:  ),ad2(:,:  ),ad3(:,:  ),ad4(:,:  ),ad5(:,:)
 complex(8),allocatable :: am1(:,:  ),am2(:,:  ),am3(:,:  ),am4(:,:  ),am5(:,:),am6(:,:),am7(:,:),am8(:,:),am9(:,:),am10(:,:)
 real   (8),allocatable :: vd1(:    ),vd2(:    ),vd3(:    ),vd4(:    )
 complex(8),allocatable :: vm1(:    ),vm2(:    ),vm3(:    ),vm4(:    ),vm5(:  ),vm6(:  ),vm7(:  ),vm8(:  ),vm9(:  ),vm10(:  )
 complex(8)             :: sm(3:5,6:8)
 real   (8)             :: sd(3:5,6:8)
 character(1000)        :: string
!-----------------------------------------------------------------------------------------------------------------------------------
 open(unit=f_out ,file='test.dat');  open(unit=f_math,file='test.m'); f_mout   =f_out ! output unit for matrix modules set to f_out
 call matrix_random_init
 f_mout = stdout
!-----------------------------------------------------------------------------------------------------------------------------------
 m=531;  n=m;  is=13;js=21            ; ie = is + m - 1; je = js + n -1
!-----------------------------------------------------------------------------------------------------------------------------------
 allocate(vm1(is:ie)); call random_number(vm1,sigma=PI)
 allocate(vd1(is:ie)); call random_number(vd1,sigma=PI)
 allocate(vm2(is:ie))
 allocate(vd2(is:ie))

!-----------------------------------------------------------------------------------------------------------------------------------
if(m <= 20)then
 call print(vm1,stdout,fmt='F12.3',name='m1')
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     MF'
 vm2 = sort(vm1,by='MF') ! MR/MF   modulus   reverse/modulus   forward
 do i = is,ie
  print '(20F15.5)', Real(vm2(i),kind=dp), AIMAG(vm2(i)), ABS(vm2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     MR'
 vm2 = sort(vm1,by='MR') ! MR/MF   modulus reverse/modulus forward
 do i = is,ie
  print '(20F15.5)', Real(vm2(i),kind=dp), AIMAG(vm2(i)), ABS(vm2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     RF'
 vm2 = sort(vm1,by='RF') ! MR/MF   modulus reverse/modulus forward
 do i = is,ie
  print '(20F15.5)', Real(vm2(i),kind=dp), AIMAG(vm2(i)), ABS(vm2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     RR'
 vm2 = sort(vm1,by='RR') ! MR/MF   modulus reverse/modulus forward
 do i = is,ie
  print '(20F15.5)', Real(vm2(i),kind=dp), AIMAG(vm2(i)), ABS(vm2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     IF'
 vm2 = sort(vm1,by='IF') ! MR/MF   modulus reverse/modulus forward
 do i = is,ie
  print '(20F15.5)', Real(vm2(i),kind=dp), AIMAG(vm2(i)), ABS(vm2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     IR'
 vm2 = sort(vm1,by='IR') ! MR/MF   modulus reverse/modulus forward
 do i = is,ie
  print '(20F15.5)', Real(vm2(i),kind=dp), AIMAG(vm2(i)), ABS(vm2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
!-----------------------------------------------------------------------------------------------------------------------------------

 call print(vd1,stdout,fmt='F12.3',name='d1')

 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     RF'
 vd2 = sort(vd1,by='RF') 
 do i = is,ie
  print '(20F15.5)', vd2(i), abs(vd2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     RR'
 vd2 = sort(vd1,by='RR') 
 do i = is,ie
  print '(20F15.5)', vd2(i), abs(vd2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     AF'
 vd2 = sort(vd1,by='AF') 
 do i = is,ie
  print '(20F15.5)', vd2(i), abs(vd2(i))
 end do
 print '(A)','#------------------------------------------------------------------------------------------'
 print '(A)','#     AR'
 vd2 = sort(vd1,by='AR') 
 do i = is,ie
  print '(20F15.5)', vd2(i), abs(vd2(i))
 end do
end if  ! if( m < 20)
!-----------------------------------------------------------------------------------------------------------------------------------

print '(A)','#------------------------------------------------------------------------------------------'
vm2 = sort(vm1,by='MF'); j = 0; do i= is,ie-1;  if(ABS  (vm2(i))>ABS  (vm2(i+1))) j = j +1; end do; print '(A,I6)','#     MF: ',j
vm2 = sort(vm1,by='MR'); j = 0; do i= is,ie-1;  if(ABS  (vm2(i))<ABS  (vm2(i+1))) j = j +1; end do; print '(A,I6)','#     MR: ',j
vm2 = sort(vm1,by='RF'); j = 0; do i= is,ie-1;  if(REAL (vm2(i))>REAL (vm2(i+1))) j = j +1; end do; print '(A,I6)','#     RF: ',j
vm2 = sort(vm1,by='RR'); j = 0; do i= is,ie-1;  if(REAL (vm2(i))<REAL (vm2(i+1))) j = j +1; end do; print '(A,I6)','#     RR: ',j
vm2 = sort(vm1,by='IF'); j = 0; do i= is,ie-1;  if(AIMAG(vm2(i))>AIMAG(vm2(i+1))) j = j +1; end do; print '(A,I6)','#     IF: ',j
vm2 = sort(vm1,by='IR'); j = 0; do i= is,ie-1;  if(AIMAG(vm2(i))<AIMAG(vm2(i+1))) j = j +1; end do; print '(A,I6)','#     IR: ',j
print '(A)','#------------------------------------------------------------------------------------------'
vd2 = sort(vd1,by='AF'); j = 0; do i= is,ie-1;  if(ABS  (vd2(i))>ABS  (vd2(i+1))) j = j +1; end do; print '(A,I6)','#     AF: ',j
vd2 = sort(vd1,by='AR'); j = 0; do i= is,ie-1;  if(ABS  (vd2(i))<ABS  (vd2(i+1))) j = j +1; end do; print '(A,I6)','#     AR: ',j
vd2 = sort(vd1,by='RF'); j = 0; do i= is,ie-1;  if(     (vd2(i))>     (vd2(i+1))) j = j +1; end do; print '(A,I6)','#     RF: ',j
vd2 = sort(vd1,by='RR'); j = 0; do i= is,ie-1;  if(     (vd2(i))<     (vd2(i+1))) j = j +1; end do; print '(A,I6)','#     RR: ',j

!if(m<7) call print(am5,stdout,fmt='F12.3',name='diag')
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
