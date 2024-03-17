!----------------------------------------------------------------------------------------------------------------------
!
! 1. Compile library:     make gnu
! 2. Compile file:        gfortran -I./lib/gnu/ hello.f90 -L./lib/gnu/ -lMatrix -llapack -lblas -o hello
!
!----------------------------------------------------------------------------------------------------------------------
program hello_world
 use :: array_mod
 complex(8), allocatable :: M(:,:),v(:), v2(:), w(:,:)
 complex(8), allocatable :: evec(:,:), eval(:)
 complex(8), allocatable :: P(:,:), Pinv(:,:)
 complex(8), allocatable :: ldet(:),lpf(:)
 complex(8)              :: z
 complex(8), allocatable :: s0(:,:), s1(:,:), s2(:,:), s3(:,:), gamma(:,:)

 n = 8

 allocate(m(n,n),v(n))

 !create random entries in the arrays:  overloaded random_number() for complex arrays, and for Gaussian random numbers
 call matrix_random_init                          ! seed random number from /dev/urandom
 call random_number(M,sigma=1._8); call random_number(v,sigma=1._8) 
 
 z    = determinant(M)
 ldet = lndet      (M)                            !logarithm of determinant for large matrices

 v2   = lmatmul(M,v)                              !use BLAS for for matrix multiplication
 
 w    = eigenvectors(M,sortby='RF')               !Eigenvalues+right eigenvectors of general complex matrix, ordered by ascending real part
 eval = w(:,1 )                                   !first column of w stores eigenvalues
 evec = w(:,2:)                                   !the rest, the eigenvectors, ordered as the eigenvalues
 
 i=3;
 v2 = lmatmul(M,evec(:,i)) - eval(i) * evec(:,i)  !M.v - λ v = 0
 call print( v2 )                                 !The library's printing subroutine for allocatable arrays

 !Diagonalize m:
 P    = evec ; Pinv = inverse(P)
 M    = lmatmul(Pinv,lmatmul(M,P))                ! P^{-1} . M . P
 
 call print(M,FMT='F8.3',name='Mdiag')

 !Compute Pfaffian:
 call random_number    (M)
 call antisymmetric_set(M)
 z    =   Pfaffian(M)                             !Pfaffian of antisymmetric matrix
 lpf  = lnPfaffian(M)                             !log(Pfaffian) for large matrices
 call print(lpf,form='Mathematica',name='lpf')    !copy/paste to define in Mathematica

 !Tensor product. Use Pauli matrices:
 s0    = PauliMatrix(0); s1 = PauliMatrix(1); s2 = PauliMatrix(2); s3 = PauliMatrix(3);
 gamma = tensorprod(s2,s2,s0,s3)
 call print(gamma,form='Gamma',name='gamma')

end program hello_world
