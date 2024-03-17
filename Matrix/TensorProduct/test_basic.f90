!
! gfortran -I../lib/gnu/ test_basic.f90 -L../lib/gnu -lMatrix -llapack -o t
program     testme
 use      array_mod
 implicit none
!------------------------------------
 integer    , parameter   :: dp = 8
 complex(dp), allocatable :: H1(:,:), H2(:,:), H3(:,:) ! Hermitian Matrices
 complex(dp), allocatable :: X (:,:), Y (:,:), Z (:,:)
 real   (dp), allocatable :: A (:,:), B (:,:), C (:,:)
 complex(dp), allocatable :: Xh(:,:), Yh(:,:), Zh(:,:) ! Hermitian Conjgugates
 complex(dp), allocatable :: Xt(:,:), Yt(:,:), Zt(:,:) ! transposes
 complex(dp), allocatable :: Xi(:,:), Yi(:,:), Zi(:,:) ! inverses
 real   (dp), allocatable :: At(:,:), Bt(:,:), Ct(:,:) ! transposes
 real   (dp), allocatable :: Ai(:,:), Bi(:,:), Ci(:,:) ! inverses
 complex(dp), allocatable :: e1(:  ), e2(:  ), e3(:  ) ! eigenvalues
 real   (dp), allocatable :: ev(:  )
 integer                  :: Ni, Nj, Nm, Nn, Np, Nq
 integer                  ::  i,  j,  m,  n,  p,  q
!------------------------------------
!Test only square matrices here:
 Nm= 22; Nn = Nm  ! Matrix sizes: X(Nm,Nn),Y(Np,Nq)   Z(Ni,Nj) = X \otimes Y
 Np= 12; Nq = Np  !               A(Nm,Nn),B(Np,Nq)   C(Ni,Nj) = A \otimes B
!------------------------------------
 call init
!------------------------------------

 Z  = tensorprod(X ,Y ) ; C  = tensorprod(A ,B ) ; Ni = size(Z,1) ; Nj = size(Z,2)
 Zt = tensorprod(Xt,Yt) ; Ct = tensorprod(At,Bt) ; 
 Zi = tensorprod(Xi,Yi) ; Ci = tensorprod(Ai,Bi)
 Zh = tensorprod(Xh,Yh) ; 
                 
!call print(Z,fmt='F8.3',name='Z')
 print *,''
 print *,'---------------------------------------------------------------------------'
 print *,'trZ   : ',  trace         (Z) - trace      (X)     *  trace      (Y)
 print *,'trZ2  : ',  trace2        (Z) - trace2     (X)     *  trace2     (Y)
 print *,'Zt    : ',  norm(transpose(Z) - Zt)
 print *,'Zh    : ',  norm(hermitian(Z) - Zh)
 if(Ni < 500) &
 print *,'detZ  : ', (determinant   (Z) - determinant(X)**Np *  determinant(Y)**Nm)/abs(determinant   (Z))
 print *,'lndetZ: ',  lndet         (Z) - lndet      (X)* Np -  lndet      (Y)* Nm 
 print *,''
 print *,'---------------------------------------------------------------------------'
 print *,'trC   : ',  trace         (C) - trace      (A)     *  trace      (B)
 print *,'trC2  : ',  trace2        (C) - trace2     (A)     *  trace2     (B)
 print *,'Ct    : ',  norm(transpose(C) - Ct)
 print *,'Ci    : ',  norm(inverse  (C) - Ci)
 if(Ni < 500) &
 print *,'deC   : ', (determinant   (C) - determinant(A)**Np *  determinant(B)**Nm)/abs(determinant   (C))
 print *,'lndetC: ',  lndet         (C) - lndet      (A)* Np -  lndet      (B)* Nm 
 print *,''

!Test if eigenvalues of the tensor product, are the product of eigenvalues: use Hermitian matrices
 H1   = X; call hermitian_set(H1);
 H2   = Y; call hermitian_set(H2);
 H3   = tensorprod (H1,H2)
 e1   = eigenvalues(H1,mtype='ZH'); e2 = eigenvalues(H2,mtype='ZH') ! ZH: Hermitian alogirhtm => real eigenvalues, sorted by value 
 e3   = eigenvalues(H3,mtype='ZH'); allocate(ev(Ni))
 i    = 0
 do m = 1, Nm; do p = 1, Np
  i   = i + 1; ev(i) = DBLE(e1(m)) * DBLE(e2(p))
 end do      ; end do
 ev = sort(ev)        ! eigenvalues need to be sorted to be compared
 print *,'---------------------------------------------------------------------------'
 print '(A,1000F10.6)','EV1    : ',norm(ev-DBLE(e3))
 print *,''


!Z = tensorprod(PauliMatrix(1),PauliMatrix(1),PauliMatrix(1))
!call print(Z,fmt='F8.3',name='Z')

!Z = tensorprod(PauliMatrix(1),PauliMatrix(1),PauliMatrix(1),PauliMatrix(1))
!call print(Z,fmt='F8.3',name='Z')

!Z = tensorprod(PauliMatrix(1),PauliMatrix(1),PauliMatrix(1),PauliMatrix(1),PauliMatrix(1))
!call print(Z,fmt='F8.3',name='Z')

!Z = tensorprod(PauliMatrix(1),PauliMatrix(1),PauliMatrix(1),PauliMatrix(1),PauliMatrix(1),PauliMatrix(1))
!call print(Z,fmt='F8.3',name='Z')

!A = PauliMatrix(1) ; B = PauliMatrix(3)
!call print(A,fmt='F8.3',name='A')

!C = tensorprod(A,B)
!call print(C,fmt='F8.3',name='AB')

!C = tensorprod(A,B)
!call print(C,fmt='F8.3',name='BA')

!C = tensorprod(B,B)
!call print(C,fmt='F8.3',name='BB')


!C = tensorprod(A,A)
!call print(C,fmt='F8.3',name='C2')

!C = tensorprod(A,A,A)
!call print(C,fmt='F8.3',name='C3')

!C = tensorprod(A,A,A,A)
!call print(C,fmt='F8.3',name='C4')

!C = tensorprod(A,A,A,A,A)
!call print(C,fmt='F8.3',name='C5')

!C = tensorprod(A,A,A,A,A,A)
!call print(C,fmt='F8.3',name='C6')

!------------------------------------
contains
!------------------------------------
 subroutine     init

!------------------------------------
 allocate(X(Nm,Nn),Y(Np,Nq))
 allocate(A(Nm,Nn),B(Np,Nq))
!------------------------------------
 X = 0.0_dp ; Y = 1.0_dp
 call matrix_random_init
 call random_number(X)
 call random_number(Y)
 call random_number(A)
 call random_number(B)
!Y = identityMatrix(Np)
 Xh = hermitian(X) ; Yh = hermitian(Y)
 Xt = transpose(X) ; Yt = transpose(Y)
 Xi = inverse  (X) ; Yi = inverse  (Y)
 At = transpose(A) ; Bt = transpose(B)
 Ai = inverse  (A) ; Bi = inverse  (B)
 end subroutine init
!------------------------------------
end program testme
