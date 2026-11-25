!
! gfortran -I../lib/gnu/ test_basic.f90 -L../lib/gnu -lMatrix -llapack -o t
program     testme
 use        array_mod
 implicit none
!------------------------------------
 integer    , parameter   :: dp = 8
 complex(dp), allocatable :: s0(:,:), s1(:,:), s2(:,:), s3(:,:)
 complex(dp), allocatable :: g1(:,:), g3(:,:)
 complex(dp), parameter   :: ZERO=(0.0_dp,0.0_dp),ONE=(1.0_dp,0.0_dp),IMU=(0.0_dp,1.0_dp)
!------------------------------------
 call init
!------------------------------------

 s0 = PauliMatrix(0); s1 = PauliMatrix(1); s2 = PauliMatrix(2); s3 = PauliMatrix(3); 
 print *,'-----------------------------------------------'
 call printgamma(s2)
 call print(s2,form='Gamma',name='s2')
 print *,'-----------------------------------------------'
 g1 = IMU * tensorprod(s2,s2,s2,s2)
 g3 = IMU * tensorprod(s2,s2,s0,s3)
 call print(g1,form='Gamma',name='g1')!; call printgamma(g1)
 print *,'-----------------------------------------------'
 call print(g3,form='Gamma',name='g3')!; call printgamma(g3)
 print *,'-----------------------------------------------'
!------------------------------------
contains
!------------------------------------
 subroutine     init

!------------------------------------
 call matrix_random_init
 
 end subroutine init
!------------------------------------
 subroutine     printgamma(g)
  complex  (dp)   :: g(:,:)
  integer         :: N, i,j
  character(2000) :: s
  character(3)    :: e
  real     (dp)   :: cut=1.0e-7_dp

  N = size(g,1)

  do i  = 1, N
   s    = ''
   do j = 1, N
    if     ( abs(g(i,j)-ZERO) < cut)then
     e = '  0'
    else if( abs(g(i,j)- ONE) < cut)then
     e = '  1'
    else if( abs(g(i,j)+ ONE) < cut)then
     e = ' -1'
    else if( abs(g(i,j)- IMU) < cut)then
     e = '  i'
    else if( abs(g(i,j)+ IMU) < cut)then
     e = ' -i'
    else
     e = ' NN'
    end if
    s = trim(s) // e
   end do
   print '(A)', trim(s)
  end do
 end subroutine printgamma
!------------------------------------
end program testme
