!
! gfortran -I../lib/gnu/ test_basic.f90 -L../lib/gnu -lMatrix -llapack -o t
program     testme
 use      array_mod
 implicit none
!------------------------------------
 integer    , parameter   :: dp = 8
 complex(dp), allocatable :: H1(:,:), H2(:,:), H3(:,:) ! Hermitian Matrices
 complex(dp), allocatable :: X (:,:), Y (:,:), Z (:,:), ZZ(:,:)
 complex(dp), allocatable :: X1(:,:), Y1(:,:), Z1(:,:)
 complex(dp), allocatable :: X2(:,:), Y2(:,:), Z2(:,:)
 real   (dp), allocatable :: A (:,:), B (:,:), C (:,:), CC(:,:)
 real   (dp), allocatable :: A1(:,:), B1(:,:), C1(:,:)
 real   (dp), allocatable :: A2(:,:), B2(:,:), C2(:,:)
 integer                  :: Ni , Nj , Nm , Nn , Np , Nq
 integer                  :: Ni1, Nj1, Nm1, Nn1, Np1, Nq1
 integer                  :: Ni2, Nj2, Nm2, Nn2, Np2, Nq2
 integer                  ::  i,  j,  m,  n,  p,  q
!------------------------------------
!Test 
!Test matrix multiplication: C1 = A1 x B1   C2 = A1 x B2 => C1.C2 = A1.A2 x B1.B2
 Nm1 = 3; Nn1 = 5; Nm2 = Nn1 ; Nn2 = 11 ; Nm = Nm1 ; Nn = Nn2 ! X1(Nm1,Nn1) X2(Nm2,Nn2)   => X(Nm,Nn) = X1.X2(Nm1,Nn2)
 Np1 = 7; Nq1 = 8; Np2 = Nq1 ; Nq2 = 9  ; Np = Np1 ; Nq = Nq2 ! Y1(Np1,Nq1) Y2(Np2,Nq2)   => Y(Np,Nq) = Y1.Y2(Np1,Nq2)
!------------------------------------
 call init
 print *,''
 print *,'---------------------------------------------------------------------'
!------------------------------------
!Test the property:
!Z1 = X1 x Y1 Z2 = X2 x Y2 => Z1.Z2 = X1.X2  x  Y1.Y2  (.: matrix multiplication  x: tensor product)
 X = matmul    (X1,X2) ; Y  = matmul    (Y1,Y2); ZZ = tensorprod(X ,Y )
 Z1= tensorprod(X1,Y1) ; Z2 = tensorprod(X2,Y2); Z  = matmul    (Z1,Z2)
 print *,'Z  : ',norm(Z-ZZ)


 print *,''
 print *,'---------------------------------------------------------------------'

!------------------------------------
!Test the property:
!C1 = A1 x B1 C2 = A2 x B2 => C1.C2 = A1.A2  x  B1.B2  (.: matrix multiplication  x: tensor product)
 A = matmul    (A1,A2) ; B  = matmul    (B1,B2); CC = tensorprod(A ,B )
 C1= tensorprod(A1,B1) ; C2 = tensorprod(A2,B2); C  = matmul    (C1,C2)
 print *,'C  : ',norm(C-CC)

 print *,''
 print *,'---------------------------------------------------------------------'

!Test associativity: X1 x (Y2 x Z) = (X1 x Y2) x Z
 X = tensorprod(Y2,Z) ; X = tensorprod(X1,X)
 Y = tensorprod(X1,Y2); Y = tensorprod(Y ,Z)
 print *,'XYZ: ',norm(X-Y)

 print *,''
!------------------------------------
contains
!------------------------------------
 subroutine     init

!------------------------------------
 allocate(X1(Nm1,Nn1),Y1(Np1,Nq1))
 allocate(X2(Nm2,Nn2),Y2(Np2,Nq2))

 allocate(A1(Nm1,Nn1),B1(Np1,Nq1))
 allocate(A2(Nm2,Nn2),B2(Np2,Nq2))
 
!------------------------------------
 call matrix_random_init
 call random_number(X1); call random_number(X2); 
 call random_number(Y1); call random_number(Y2); 
 call random_number(A1); call random_number(A2); 
 call random_number(B1); call random_number(B2); 

 end subroutine init
!------------------------------------
end program testme
