program             testme
 use                matrix_mod_common
 use                matrix_mod_matrix
 use                matrix_mod_array
 implicit none
 integer   ,parameter   :: f_out = 13, f_math = 14
 integer                :: i,j
 real   (dp)            :: r,r1,r2,r3,sigma
 complex(dp)            :: z,z1,z2,z3
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
 m=4;  n=m;  is=1;js=1            ; ie = is + m - 1; je = js + n -1
! allocate(ad1(is:ie,js:je),ad2(js:je,is:ie),ad3(is:ie,is:ie),ad4(js:je,js:je),ad5(js:je,js:je)) ! ad3 = MATMUL(ad1,ad2) ad4 = MATMUL(ad2,ad1)
! allocate(am1(is:ie,js:je),am2(js:je,is:ie),am3(is:ie,is:ie),am4(js:je,js:je),am5(js:je,js:je)) ! am3 = MATMUL(am1,am2) am4 = MATMUL(am2,am1)
! allocate(vd1(      js:je),vd2(      is:ie),vd3(      is:ie),vd4(      js:je)) 
! allocate(vm1(      js:je),vm2(      is:ie),vm3(      is:ie),vm4(      js:je)) ! vm3 = MATMUL(am1,vm1) vm4 = MATMUL(vm2,am1)
!-----------------------------------------------------------------------------------------------------------------------------------

 m = 417; n = m

 allocate(ad1(n,n)) ; allocate(am1(n,n)); allocate(am7(n,n))
 call random_number(ad1);call random_number(am7,sigma=PI);
!call print(ad1,20,form='m',name='m')
!call symmetric_set(ad1)
 if(m<7) call print(ad1,stdout,fmt='F12.3',name='A')
 am1=ad1

 vm8 =eigenvalues (am7)
 am8 =eigenvectors(am7)
 vm9 =am8(:,1)
 am10=am8(:,2:)
!print *,'Deval89:', norm(vm9-vm8),norm(vm8) !wrong: the eigenvalues come out at a different order for large n. In order to compare:
 do i = 1, m
  write(20,'(1000G28.16)') vm9(i)
  write(21,'(1000G28.16)') vm8(i)
 end do     
!Then: sort fort.20 > hi1;sort fort.21 > hi2; paste hi[12] |awk -v OFMT="%.17g" '{print $1,$2,$3,$4,sqrt(($1-$3)^2+($2-$4)^2)}'|ct|lesss

 vm2=eigenvalues (ad1,mtype='DG')
 am5=eigenvectors(ad1,mtype='DG')
 vm3=am5(:,1 )
 am3=am5(:,2:)

 vm1=eigenvalues (ad1)
 am6=eigenvectors(ad1)
 vm4=am6(:,1)
 am4=am6(:,2:)

 if(m<7) call print(vm2,stdout,fmt='F12.3',name='evalS2')
!if(m<7) call print(vm3,stdout,fmt='F12.3',name='evalS3')
!if(m<7) call print(vm1,stdout,fmt='F12.3',name='evalG1')
!if(m<7) call print(vm4,stdout,fmt='F12.3',name='evalG4')
!if(m<7) call print(am3,stdout,fmt='F12.3',name='evecS')
 if(m<7) call print(am4,stdout,fmt='F12.3',name='evecG')

!print *,'Deval32:', norm(vm3-vm2),norm(vm3)!wrong: the eigenvalues come out at a different order for large n. In order to compare:
 do i = 1, m
  write(22,'(1000G28.16)') vm3(i)
  write(23,'(1000G28.16)') vm2(i)
 end do     
!Then: sort fort.22 > hi1;sort fort.23 > hi2; paste hi[12] |awk -v OFMT="%.17g" '{print $1,$2,$3,$4,sqrt(($1-$3)^2+($2-$4)^2)}'|ct|lesss

!print *,'Deval41:', norm(vm4-vm1),norm(vm4)
!print *,'Deval43:', norm(vm4-vm3)
!print *,'Devec43:', norm(am4-am3)

 r    = 0.0_dp
 do i = 1, m
  r1  = norm( MATMUL(am1,am3(:,i)) - vm3(i) * am3(:,i) )
  r   = r + r1
! print *,'r1= ',i,r1
 end do
 print *, 'r = ',r/n

 ! if V = (v1 ... vn) Λ = diag(λ1,...,λn) then  A.V = V.Λ
 am7 = diagonal(vm3)
!if(m<7) call print(am7,stdout,fmt='F12.3',name='diagEV')
 print *, 'rS = ', norm( MATMUL(am1,am3) - MATMUL(am3,am7) )
 print *, 'rG = ', norm( MATMUL(am1,am4) - MATMUL(am4,diagonal(vm4)) )





!-----------------------------------------------------------------------------------------------------------------------------------
!am2 = inverse(am1);am3 = MATMUL(am1,am2)-diagonal(ONE   ,m);am4 = MATMUL(am2,am1)-diagonal(ONE   ,m)
!ad2 = inverse(ad1);ad3 = MATMUL(ad1,ad2)-diagonal(1.0_dp,m);ad4 = MATMUL(ad2,ad1)-diagonal(1.0_dp,m)

! call random_number(am1);
! call hermitian_set(am1)
! if(m<7) call print(am1,stdout,fmt='F12.3',name='m1') 
! call array2_zgeev(am1,vm2,am2,job='V') !arrays vm2, am2 get lbound = 1
! vm4=eigenvalues(am1,mtype='ZH')
! deallocate(am5); deallocate(vm3,am3);allocate(vm3(n),am3(n,n))
! am5=eigenvectors(am1,mtype='ZH')
! vm3=am5(:,1 )
! am3=am5(:,2:)
! if(m<7) call print(vm2,stdout,fmt='F12.3',name='evalz')
! if(m<7) call print(vm4,stdout,fmt='F12.3',name='evale')
! if(m<7) call print(am2,stdout,fmt='F12.3',name='evecz')
! if(m<7) call print(am3,stdout,fmt='F12.3',name='evech')
!!print *,'Deval: zheev-zgeev (N):', norm(vm4-vm2) ! we cannot compare, the ordering of eigenvalues is different in zgeev and zheev
!!print *,'Deval: zheev-zgeev (V):', norm(vm3-vm2)
! print *,'Deval: zheev-zheev    :', norm(vm3-vm4)
! print *,'Devec3 :', norm(am3-am2)
! ! if V = (v1 ... vn) Λ = diag(λ1,...,λn) then  A.V = V.Λ
! print *, 'rG = ', norm( MATMUL(am1,am2) - MATMUL(am2,diagonal(vm2)) )
! print *, 'rH = ', norm( MATMUL(am1,am3) - MATMUL(am3,diagonal(vm3)) )
!-----------------------------------------------------------------------------------------------------------------------------------
!r    = 0.0_dp
!do i = 1, m
! r1  = norm( MATMUL(am1,am2(:,i)) - vm2(i) * am2(:,i) )
! r   = r + r1
! !print *,'r1= ',i,r1
!end do
!print *, 'r = ',r/n,lbound(am2),lbound(vm2)

!if(m<7) call print(am5,stdout,fmt='F12.3',name='diag')
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
