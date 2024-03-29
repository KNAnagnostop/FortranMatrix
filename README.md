# Fortran Matrix Modules

The modules of the library Matrix provide convenient interfaces to many frequently used linear algebra calculations for real(8) and complex(8) arrays. 

These include matrix multiplication, trace, hermitianization, special matices, Pauli matrices, sorting, as well as more complicated calculations like computing eigenvalues, eigenvectors, determinants, Pfaffians, inverse of a matrix, tensor products. It overloads the random_number() subroutine for complex(8) arrays, as well as for generating Gaussian random numbers. 

The programming is aiming at performance, and the fortran 77 version of Lapack and BLAS is used under the hood. Some performance critical routines, like matrix-matrix and matrix-vector multiplication, are given in the form of subroutine for optimal performance, as well as as functions and operator overloading for convenience. 

## Documentation

The up to date documentation can be found in the [array  module documentation](https://docs.google.com/document/d/19KoFFvpxTcm9FN1zGqdGaBopyJhl2vYA0D1mEHbTUZ8/edit?usp=sharing).

The [matrix module documentation](https://docs.google.com/document/d/1PLhbGWSkTO2lGfq7dNz5SUUcmbekR8NTH-dGsxdLsfs/edit?usp=sharing). The matrix module is mostly written for educational reasons, and is not to be used for HPC applications.


## Downloading

You may download the [source code also from here.](https://physics.ntua.gr/konstant/PUB/Matrix.tgz)


## Hello World Program

See the included file hello.f90 for details on how to compile the code below, and how to run it. 

``` fortran
program hello_world
 use :: array_mod
 complex(8), allocatable :: M(:,:),v(:), v2(:), w(:,:)
 complex(8), allocatable :: evec(:,:), eval(:)
 complex(8), allocatable :: P(:,:), Pinv(:,:)
 complex(8), allocatable :: ldet(:),lpf(:)
 complex(8)              :: z
 complex(8), allocatable :: s0(:,:), s1(:,:), s2(:,:), s3(:,:), gamma(:,:)

 n = 8  ; allocate(m(n,n),v(n))

 !overloaded random_number() for complex arrays, and for Gaussian random numbers
 !create random entries in the arrays in a Gaussian distribution with σ=1:
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
 call print( v2 )                                 !The library printing subroutine for allocatable arrays

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
```
 
 
