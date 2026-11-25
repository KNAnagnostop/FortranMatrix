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
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
