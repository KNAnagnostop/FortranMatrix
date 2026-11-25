!
! gfortran -I../lib/gnu/ test_basic.f90 -L../lib/gnu -lMatrix -llapack -o t
program     testme
 use      array_mod
 implicit none
!------------------------------------
 integer    , parameter   :: dp = 8
 complex(dp), allocatable :: X (:,:), Y (:,:), Z (:,:), ZZ(:,:)
 complex(dp), allocatable :: Xv(:  )
 real   (dp), allocatable :: A (:,:), B (:,:), C (:,:), D (:,:), CC (:,:)
 real   (dp), allocatable :: Av(:  ), Bv(:  ), Cv(:  ), Dv(:  ), CCv(:  )
 real   (dp), allocatable ::          Bt(:,:)
 integer                  :: Ni , Nj , Nm , Nn , Np , Nq
 integer                  :: Ni1, Nj1, Nm1, Nn1, Np1, Nq1
 integer                  :: Ni2, Nj2, Nm2, Nn2, Np2, Nq2
 integer                  ::  i,  j,  m,  n,  p,  q
!------------------------------------
!Test 
 Nm = 3; Nn = 9
 Np = 2; Nq = 2
 Nm1= 12; Nm2= 7; Nn2=6
!------------------------------------
 call init
 print *,''
 print *,'---------------------------------------------------------------------'

 !-------------------------------------------
 ! Check identity:
 ! A . D . B = C <=> (Bt x A).vec(D) = vec(C)   x:outer product, .: matrix product vec(.): vectorization of matrix
 ! A(Nm1,Nn) D(Nn,Nm2) B(Nm2,Nn2)  C(Nm1,Nn2)
 !-------------------------------------------
 ! LHS:
 CC=matmul(D,B); C = matmul(A,CC); Cv=vectorize(C)
 !RHS:
 Dv=vectorize(D)
 CCv = matmul(tensorprod(transpose(B),A),Dv)

 print *,'ADBC: ',norm(CCv-Cv)
 




!------------------------------------
 stop
!------------------------------------
 call print(X ,fmt='F8.3',name='X' )
 Xv = vectorize(X)
 call print(Xv,fmt='F8.3',name='Xv')

 print *,''
 print *,'---------------------------------------------------------------------'
!------------------------------------
 call print(A ,fmt='F8.3',name='A' )
 Av = vectorize(A)
 call print(Av,fmt='F8.3',name='Av')

 print *,''
!------------------------------------
contains
!------------------------------------
 subroutine     init

!------------------------------------
 allocate(X(Nm ,Nn),Y(Np,Nq )           )
 allocate(A(Nm1,Nn),D(Nn,Nm2),B(Nm2,Nn2))

 
!------------------------------------
 call matrix_random_init
 call random_number(X); call random_number(Y); 
 call random_number(A); call random_number(B);  call random_number(D); 

 end subroutine init
!------------------------------------
end program testme
