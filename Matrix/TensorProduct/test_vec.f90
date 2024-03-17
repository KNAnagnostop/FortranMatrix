!
! gfortran -I../lib/gnu/ test_vec.f90 -L../lib/gnu -lMatrix -llapack -o t
program     testme
 use      array_mod
 implicit none
!------------------------------------
 integer    , parameter   :: dp = 8
 complex(dp), allocatable :: u (:  ), v (:  )
 complex(dp), allocatable :: uc(:,:), vc(:,:) ! vector column
 complex(dp), allocatable :: ur(:,:), vr(:,:) ! vector row
 complex(dp), allocatable :: Z (:,:), Z1(:,:), Z2(:,:)
 real   (dp), allocatable :: e (:  ), g (:  )
 real   (dp), allocatable :: ec(:,:), gc(:,:)
 real   (dp), allocatable :: er(:,:), gr(:,:)
 real   (dp), allocatable :: C (:,:), C1(:,:), C2(:,:)
 integer                  :: Ni, Nj, Nm, Nn, Np, Nq
 integer                  ::  i,  j,  m,  n,  p,  q
!------------------------------------
 Nm= 24; Nn = Nm  ! Sizes: u(Nm),v(Np)
 Np= 32; Nq = Np  ! 
!------------------------------------
 call init
!------------------------------------

 print *,''
 print *,'------------------------------------------------------------------------------'

!call print(v,fmt='F8.3',name='v')
 uc=vec2col(u); ur=vec2row(u)
 vc=vec2col(v); vr=vec2row(v)
 print *,'u  : ',norm(u-uc(:,1)),norm(u-ur(1,:)) 
 print *,'v  : ',norm(v-vc(:,1)),norm(v-vr(1,:)) 

 Z = tensorprod(u,v); Z1 = tensorprod(ur,vc) ; Z2 = tensorprod(vec2row(u),vec2col(v))! must give the same results
!print *,'nu : ', size(u)               ,' nur: ',size(ur,1), size(ur,2), ' nuc: ',size(uc,1), size(uc,2)
!print *,'nv : ', size(v)               ,' nvr: ',size(vr,1), size(vr,2), ' nvc: ',size(vc,1), size(vc,2)
!print *,'nuv: ', size(Z, 1), size(Z, 2),' nZ1: ',size(Z1,1), size(Z1,2), ' nZ2: ',size(Z2,1), size(Z2,2)
 print *,' uv: ', norm(Z-Z1), norm(Z-Z2), norm(Z1-Z2)

 Z = tensorprod(v,u); Z1 = tensorprod(vr,uc) ; Z2 = tensorprod(vec2row(v),vec2col(u))! must give the same results
!print *,'nvu: ', size(Z, 1), size(Z, 2),' nZ1: ',size(Z1,1), size(Z1,2), ' nZ2: ',size(Z2,1), size(Z2,2)
 print *,' vu: ', norm(Z-Z1), norm(Z-Z2), norm(Z1-Z2)
 print *,''
 print *,'------------------------------------------------------------------------------'

 ec=vec2col(e); er=vec2row(e)
 gc=vec2col(g); gr=vec2row(g)
 print *,'e  : ',norm(e-ec(:,1)),norm(e-er(1,:)) 
 print *,'g  : ',norm(g-gc(:,1)),norm(g-gr(1,:)) 
 C = tensorprod(e,g); C1 = tensorprod(er,gc) ; C2 = tensorprod(vec2row(e),vec2col(g))! must give the same results
!print *,'neg: ', size(C, 1), size(C, 2),' nC1: ',size(C1,1), size(C1,2), ' nC2: ',size(C2,1), size(C2,2)
 print *,' eg: ', norm(C-C1), norm(C-C2), norm(C1-C2)

 C = tensorprod(g,e); C1 = tensorprod(gr,ec) ; C2 = tensorprod(vec2row(g),vec2col(e))! must give the same results
!print *,'nge: ', size(C, 1), size(C, 2),' nC1: ',size(C1,1), size(C1,2), ' nC2: ',size(C2,1), size(C2,2)
 print *,' ge: ', norm(C-C1), norm(C-C2), norm(C1-C2)





 !call print(u ,fmt='F8.3',name='u' )
 !call print(uc,fmt='F8.3',name='uc')
 !call print(ur,fmt='F8.3',name='ur')
 !call print(v ,fmt='F8.3',name='v' )
 !call print(vc,fmt='F8.3',name='vc')
 !call print(vr,fmt='F8.3',name='vr')


 !Z = tensorprod(u,v)
 !call print(Z,fmt='F8.3',name='Z')

 !Z = tensorprod(v,u)
 !call print(Z,fmt='F8.3',name='Z')

 !call print(e,fmt='F8.3',name='e')
 !call print(g,fmt='F8.3',name='g')

 !C = tensorprod(e,g)
 !call print(C,fmt='F8.3',name='C')

 !C = tensorprod(g,e)
 !call print(C,fmt='F8.3',name='C')
!------------------------------------
contains
!------------------------------------
 subroutine     init

!------------------------------------
 allocate(u(Nm),v(Np))
 allocate(e(Nm),g(Np))
!------------------------------------
 call matrix_random_init
 call random_number(u)
 call random_number(v)
 call random_number(e)
 call random_number(g)

 end subroutine init
!------------------------------------
end program testme
