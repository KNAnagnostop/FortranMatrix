! gfortran -I../lib/gnu/ test.f90 -L../lib/gnu -lMatrix -llapack -o t
program testme
 use array_mod
 integer    , parameter   :: dp = 8
 complex(dp), allocatable :: s0(:,:), s1(:,:), s2(:,:), s3(:,:)
 complex(dp), allocatable :: g1(:,:), g3(:,:)
 complex(dp), parameter   :: ZERO=(0._dp,0._dp),ONE=(1._dp,0._dp),IMU=(0._dp,1._dp)

 s0 = PauliMatrix(0); s1 = PauliMatrix(1); s2 = PauliMatrix(2); s3 = PauliMatrix(3);

 g1 = IMU * tensorprod(s2,s2,s2,s2); call print(g1,form='Gamma',name='g1')
 g3 = IMU * tensorprod(s2,s2,s0,s3); call print(g3,form='Gamma',name='g3')

end program testme
