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
!m1 = Matrix('g',m,n,is,js,name='m1') ;  call m1 % print(stdout,fmt='F12.3',form=form);
!m2 = inverse(m1);    m2 % name='mi'  ;  call m2 % print(stdout,fmt='F12.3',form=form);
!m3 = Matrix('1',m,n,is,js,name='II') ;  call m3 % print(stdout,fmt='F12.3',form=form);
!m4 = m1 * m2                         ;  call m4 % print(stdout,fmt='F12.3',form=form);
!print *,'|m1 * m1^-1|= ',norm(m4 - m3)
!m4 = m2 * m1                         ;  call m4 % print(stdout,fmt='F12.3',form=form);
!print *,'|m1^-1 * m1|= ',norm(m4 - m3)
!-----------------------------------------------------------------------------------------------------------------------------------
!m=80;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js,name='m1') ;call m1 % print(stdout,fmt='F12.3',form=form);
!v1 = eigenvalues (m1)                ;call v1 % print(stdout,fmt='F12.3',form=form);
!m2 = eigenvectors(m1,v2)             ;call m2 % print(stdout,fmt='F12.3',form=form);call v2 % print(stdout,fmt='F12.3',form=form);
!! if V = (v1 ... vn) Λ = diag(λ1,...,λn) then  A.V = V.Λ    -    see test_v04.f90 for array2 tests
!! !!! CAREFUL: v1 and v2 have the eigenvalues possibly sorted differently !!!
!m3 = diagonal(v2)                    ;call m3 % print(stdout,fmt='F12.3',form=form);
!print *,'|A.v - l.v|=', norm(m1 * m2 - m2 * m3)
!! determinant:
!z1 = determinant(m1)
!print *,'|det - Π λ|=', abs(z1-PRODUCT(v2%v))/abs(z1), abs(z1-PRODUCT(v1%v)) /abs(z1)
!ldet = lndet(m1)
!print *,'|ln det   |=', abs(ldet(1)), ldet(2), log(ldet(2)), abs(log(z1)-ldet(1) - log(ldet(2)))
!! The 1st term is not always zero. The 2nd term is aways zero: |ln det|= Σ ln ρ_j  , λ_j = ρ_j exp(i φ_j)
!print *,'|ldet -Σ λ|=', abs(ldet(1) + log(ldet(2)) - SUM(log(v2%v)) )/abs(ldet(1)), &
!                        abs(ldet(1) - SUM(log(ABS(v2%v))))/abs(ldet(1))
!-------------------------------------------------------------------------
! Compute Pfaffian:
!m1 % v = 0.5_dp * ( m1 % v - TRANSPOSE( m1 % v) )
!z1     = determinant(m1)
!z2     = Pfaffian(m1)
!print *,'|det- Pf^2|=', abs(z1 - z2 *z2)/abs(z1)
!lPf    = lnPfaffian(m1)
!ldet   = lnDet(m1)
!print *,'|ldet-lPf^2|=',abs((ldet(1)+log(ldet(2))) - 2.0_dp*(lPf(1)+log(lPf(2))))/abs(ldet(1)), &
!                        abs( ldet(1) - 2.0_dp * lPf(1) )/abs(ldet(1))
!-----------------------------------------------------------------------------------------------------------------------------------
!Use the Hermitian version for the computation:
!m=74;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!m1 = Matrix('g',m,n,is,js,name='m1')!;call m1 % print(stdout,fmt='F12.3',form=form);
!call m1 % hermitian_set              ;call m1 % print(stdout,fmt='F12.3',form=form);
!v1 = eigenvalues (m1)                ;call v1 % print(stdout,fmt='F12.3',form=form);
!m2 = eigenvectors(m1,v2)             ;call m2 % print(stdout,fmt='F12.3',form=form);call v2 % print(stdout,fmt='F12.3',form=form);
!m3 = diagonal(v2)                    ;call m3 % print(stdout,fmt='F12.3',form=form);
!print *,'|A.v - l.v|=', norm(m1 * m2 - m2 * m3), norm(aimag(v1 % v)), norm(aimag(v2 % v))
!-----------------------------------------------------------------------------------------------------------------------------------
!m=412;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
!d1 =DMatrix('g',m,n,is,js,name='m1') ;  call d1 % print(stdout,fmt='F12.3',form=form);
!d2 = inverse(d1);    d2 % name='mi'  ;  call d2 % print(stdout,fmt='F12.3',form=form);
!d3 =DMatrix('1',m,n,is,js,name='II') ;  call d3 % print(stdout,fmt='F12.3',form=form);
!d4 = d1 * d2                         ;  call d4 % print(stdout,fmt='F12.3',form=form);
!print *,'|d1 * d1^-1|= ',norm(d4 - d3)
!d4 = d2 * d1                         ;  call d4 % print(stdout,fmt='F12.3',form=form);
!print *,'|d1^-1 * d1|= ',norm(d4 - d3)
!-----------------------------------------------------------------------------------------------------------------------------------
m=220;n=m;is=258;js=123;if(m<10) then; form='d';else;form='i';end if
d1 =DMatrix('g',m,n,is,js,name='d1') ;call d1 % print(stdout,fmt='F12.3',form=form);
v1 = eigenvalues(d1)                 ;call v1 % print(stdout,fmt='F12.3',form=form);
m2 = eigenvectors(d1,v2)             ;call m2 % print(stdout,fmt='F12.3',form=form);call v2 % print(stdout,fmt='F12.3',form=form);
m3 = diagonal(v2)                    ;call m3 % print(stdout,fmt='F12.3',form=form);
print *,'|A.v - l.v|=', norm(d1 * m2 - m2 * m3)
r1 = determinant(d1)
print *,'|det - Π λ|=', abs(r1-PRODUCT(v2%v))/abs(r1), abs(r1-PRODUCT(v1%v)) /abs(r1)
ldetd = lndet(d1)
print *,'lndet= ',abs(r1 - exp(ldetd(1)) * ldetd(2))/abs(r1), ldetd(1), ldetd(2)

v3 = sort(v1) ; v4 = sort(v2); v3 % name = 'v3';v4 % name = 'v4'
call v3 % print(stdout,fmt='F12.3',form=form);call v3 % print(stdout,fmt='F12.3',form=form);

!call v3 % save(20); call v4 % save(21)
print *,'|v1-v2|= ', norm(v3-v4)
!does not work: for large matrices, complex conjugate pairs might be sorted differently: 
!1.5255e-13   -13.506843823629497          1.7380992486099580                 -13.506843823629648          1.7380992486099798    
!1.5255e-13   -13.506843823629497         -1.7380992486099580                 -13.506843823629648         -1.7380992486099798   
!2.73559e-13  -13.763856087733195          0.0000000000000000                 -13.763856087733469          0.0000000000000000    
!7.43231       5.8184618482390764          3.7161536472709207                  5.8184618482390835         -3.7161536472709140    
!7.43231       5.8184618482390764         -3.7161536472709207                  5.8184618482390835          3.7161536472709140    
!11.5462      -2.1360745132469088          5.7731031869837954                 -2.1360745132469097         -5.7731031869837794    
!11.5462      -2.1360745132469088         -5.7731031869837954                 -2.1360745132469097          5.7731031869837794    

!-----------------------------------------------------------------------------------------------------------------------------------
call d1 % symmetric_set              ;call d1 % print(stdout,fmt='F12.3',form=form);
v1 = eigenvalues(d1)                 ;call v1 % print(stdout,fmt='F12.3',form=form);
m2 = eigenvectors(d1,v2)             ;call m2 % print(stdout,fmt='F12.3',form=form);call v2 % print(stdout,fmt='F12.3',form=form);
m3 = diagonal(v2)                    ;call m3 % print(stdout,fmt='F12.3',form=form);
print *,'|A.v - l.v|=', norm(d1 * m2 - m2 * m3)
r1 = determinant(d1)
print *,'|det - Π λ|=', abs(r1-PRODUCT(v2%v))/abs(r1), abs(r1-PRODUCT(v1%v)) /abs(r1)
ldetd = lndet(d1)
print *,'lndet= ',abs(r1 - exp(ldetd(1)) * ldetd(2))/abs(r1), ldetd(1), ldetd(2)

v3 = sort(v1) ; v4 = sort(v2); 
print *,'|v1-v2|= ', norm(v3-v4) ! here it is ~ 0 , only real eigenvalues
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end program         testme
!===================================================================================================================================
