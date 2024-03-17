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
 m=14;  n=m;  is=1;js=1            ; ie = is + m - 1; je = js + n -1

 allocate(am1(is:ie,js:je),vm1(is:ie)); call random_number(am1,sigma=PI)     ;  ! call hermitian_set(am1)
 !---------------------------------
 vm1 = eigenvalues(am1)
 z   = SUM(LOG(vm1)); r1 = Real(z,kind=dp)     ; z1 = exp(IMU * AIMAG(z)); r2 = cos(AIMAG(z))       ; r3 = sin(AIMAG(z))
 !---------------------------------
 z2  = determinant(am1)
 !---------------------------------
 vm2 = lndet(am1)   ; r4 = real(vm2(1),kind=dp);                           r5 = real(vm2(2),kind=dp); r6 = aimag(vm2(2))

 if(m<7) call print(am1,stdout,fmt='F12.3',name='m1')
 if(m<7) call print(vm1,stdout,fmt='F12.3',name='v1')
!if(m<193)print *,'det  =',product(vm1), determinant(am1), exp(vm2(1))*vm2(2)
 if(m<193)print *,'det  =',abs(product(vm1)-z2)/abs(product(vm1)),abs(product(vm1)-exp(vm2(1))*vm2(2))/abs(product(vm1))

!print  *        ,'lndet=', r1, r2, r3!, r4, r5,r6
 print  *        ,'lndet=', ABS(r1-r4)/ABS(r1), ABS(r2-r5), ABS(r3-r6)
 print  *        ,'----------------------------------------------------------------------------------------------------------------'

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 m=512;  n=m;  is=1;js=1            ; ie = is + m - 1; je = js + n -1
 allocate(ad1(is:ie,js:je),vd1(is:ie)); call random_number(ad1,sigma=PI); if(m<7  ) call print(ad1,stdout,fmt='F12.3',name='d1')
 !---------------------------------
 vm1 = eigenvalues(ad1) ; if(m<7  ) call print(vm1,stdout,fmt='F12.3',name='v1')
 !---------------------------------
 r7  = product(vm1)
 z   = SUM(LOG(vm1));  r1 = Real(z,kind=dp) ; r2 = cos(AIMAG(z))       ; r3 = sin(AIMAG(z))
 !---------------------------------
 r4  = determinant(ad1)
 !---------------------------------
 vd2 = lndet(ad1)  ;   r5 = vd2(1)          ; r6 = vd2(2)
 

!if(m<193)print *,'detval   =',r7 , r4, exp(r5)*r6
 if(m<193)print *,'det      =',ABS(r7-r4)/r7,  abs( r7 - exp(r5)*r6 )/r7

!print  *        ,'lndetval =',r1,r2,r3
!print  *        ,'lndetval =',r5,r6
 print  *        ,'lndet    =',ABS(r5 - r1)/r5,ABS(r6 - r2),ABS(r3)
 
!if(m<7) call print(am5,stdout,fmt='F12.3',name='diag')
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
