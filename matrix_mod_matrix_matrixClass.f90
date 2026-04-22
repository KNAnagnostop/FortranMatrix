!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_matrix_matrixClass.f90 ................................................................!:.:.:
!...................................................................................................................................
!...................................................................................................................................
!.......................... type/class MatrixClass procedures ......................................................................!:.:.:
!...................................................................................................................................
!...................................................................................................................................

!...................................................................................................................................
!.......................... Used by Constructors: ..................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
!Usage:  call matrix_gaussian_set(MAT,sigma)        Default: sigma=1          (from array2_gauss_set/array2_gauss_set_d)
!Result: MAT % v is filled with gaussianly distributed random numbers by mod_array: array2_gauss_set/array2_gauss_set_d
!        If MAT % mtype is 'ZH'/'DS', MAT % v is symmetrized by a call to MAT % hermitian/MAT % symmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_gaussian_set                                    (MAT,sigma)                                              !:.:.:
  class(MatrixClass)                    :: MAT
  real (dp)      , optional, intent(in) :: sigma
  integer                               :: i,is,js
  real (dp)                             :: s
  real           , parameter            :: sqrt2 = sqrt(2.0_dp), isqrt2 = 1.0_dp/sqrt2

  is = MAT % is;       js = MAT % js;               s = 1.0_dp
  if(present(sigma))                                s = sigma

  if( MAT % mtype == 'ZH' .or. MAT % mtype == 'DS') s = s * sqrt2                                                                   ! (X1+X2)/2 have sigma2 = sigma * sqrt(2), if independent (the diagonal isn't, X1=X2)

  !--------------------------------------------------------------------------------------------------------------------------------
  select type ( MAT   )
  !--------------------------------------------------------------------------------------------------------------------------------
  class is    ( Matrix)
  !--------------------------------------------------------------------------------------------------------------------------------

   if( .not.    allocated  (MAT % v)) call matrix_error('matrix_gaussian_set: MAT % v not allocated. (Z)')

   call array2_gauss_set   (MAT % v, sigma=s)

   if( MAT % mtype(2:2) == 'H') then                                                                                                ! if(hermitian): X -> (X1+X2)/2, where X1,X2 independent random vars, so sigma->sigma * sqrt(2)
                                                                                                                                    !                The diagonal is exception, since X1=X2
    call MAT % hermitian_set                                                                                                        ! X = (X1+X2)/2

    do i = 0 , MAT % n - 1
     MAT % v(is+i,js+i) = MAT % v(is+i,js+i) * isqrt2                                                                               ! the diagonal X = (X1+X2)/2, where X1=X2, so we have to scale back to the original sigma
    end do

   end if
  !--------------------------------------------------------------------------------------------------------------------------------
  type  is    (DMatrix)
  !--------------------------------------------------------------------------------------------------------------------------------

   if( .not.    allocated  (MAT % v)) call matrix_error('matrix_gaussian_set: MAT % v not allocated. (D)')

   call array2_gauss_set_d (MAT % v, sigma=s)

   if( MAT % mtype(2:2) == 'S') then                                                                                                ! if(symmetric): X -> (X1+X2)/2, where X1,X2 independent random vars, so sigma->sigma * sqrt(2)
                                                                                                                                    !                The diagonal is exception, since X1=X2

    call MAT % symmetric_set

    do i = 0 , MAT % n - 1
     MAT % v(is+i,js+i) = MAT % v(is+i,js+i) * isqrt2                                                                               ! the diagonal X = (X1+X2)/2, where X1=X2, so we have to scale back to the original sigma
    end do

   end if
  !--------------------------------------------------------------------------------------------------------------------------------
  class default
    call                                   matrix_error('matrix_gaussian_set: MAT is of unknown class.')
  !--------------------------------------------------------------------------------------------------------------------------------
  end select
  !--------------------------------------------------------------------------------------------------------------------------------

 end subroutine     matrix_gaussian_set
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage:  call matrix_random_set(MAT)
!Result: MAT % v is filled with uniformly distributed random numbers by mod_array: array2_random_set/array2_random_set_d
!        If MAT % mtype is 'ZH'/'DS', MAT % v is symmetrized by a call to MAT % hermitian_set/MAT % symmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_random_set                                      (MAT)                                                    !:.:.:
  class(MatrixClass)                    :: MAT
  integer                               :: i, j, is, js, m, n

  m  = MAT % m  ; n  = MAT % n
  is = MAT % is ; js = MAT % js 

  !--------------------------------------------------------------------------------------------------------------------------------
  select type ( MAT   )
  !--------------------------------------------------------------------------------------------------------------------------------
  class is    ( Matrix)
  !--------------------------------------------------------------------------------------------------------------------------------
   if( .not.    allocated  (MAT % v)) call matrix_error('matrix_random_set: MAT % v not allocated. (Z)')

   call array2_random_set  (MAT % v)                                                                                               ! X1,X2 are independent random variables uniformely distributed, 
                                                                                                                                   ! X=(X1+X2)/2 has f(x)= {4x 0<x<1/2, 4(1-x) 1/2<x<1}
   if( MAT % mtype(2:2) == 'H') then                                                                                               ! We cannot use MAT % hermitian_set because X=(X1+X2)/2 is not uniformely distributed

    do  j = 0, n - 1
     MAT  % v(is+j,js+j) = REAL (MAT %  v(is+j,js+j) , kind=dp )                                                                   ! diagonal
     do i = 0, j - 1                                    
      MAT % v(is+i,js+j) = CONJG(MAT %  v(is+j,js+i)           )                                                                   ! one before the diagonal and conjugate
     end do
    end do

   end if

  !--------------------------------------------------------------------------------------------------------------------------------
  type  is    (DMatrix)
  !--------------------------------------------------------------------------------------------------------------------------------
   if( .not.    allocated  (MAT % v)) call matrix_error('matrix_random_set: MAT % v not allocated. (D)')
 
   call array2_random_set_d(MAT % v)                                                                                               ! X1,X2 are independent random variables uniformely distributed, 
                                                                                                                                   ! X=(X1+X2)/2 has f(x)= {4x 0<x<1/2, 4(1-x) 1/2<x<1}
   if( MAT % mtype(2:2) == 'S') then                                                                                               ! We cannot use MAT % symmetric_set because X=(X1+X2)/2 is not uniformely distributed

    do  j = 0, n - 1
     do i = 0, j - 1                                    
      MAT % v(is+i,js+j) = MAT %  v(is+j,js+i)                                                                                     ! one before the diagonal and symmetric
     end do
    end do

   end if

  !--------------------------------------------------------------------------------------------------------------------------------
  class default
    call                                   matrix_error('matrix_random_set: MAT is of unknown class.')
  !--------------------------------------------------------------------------------------------------------------------------------
  end select
  !--------------------------------------------------------------------------------------------------------------------------------

 end subroutine     matrix_random_set
!-----------------------------------------------------------------------------------------------------------------------------------
! Used by constructors: sets the MatrixClass metadata and makes the necessary tests:
! input: MAT, m
! optional: n = m is = 1          js = 1             If js not provided, js=is If is not provided, is=1
! sets:           ie = is + m - 1 je = js + n - 1        
!           mtype='DG'/'ZG'      name=''
! mtype: allowed: 'ZG'/'ZH' and 'DG'/'DH' for Matrix/DMatrix. For 'ZH'/'DS', only a square matrix is allowed
! MAT % v is not allocated
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_metadata_put                                    (MAT,m,n,is,js,mtype,name)                               !:.:.:
  class(MatrixClass)                    :: MAT
  integer               , intent(in)    :: m
  integer     , optional, intent(in)    :: n,is,js
  character(*), optional, intent(in)    :: mtype,name
  character(1000)                       :: err, info
  integer                               :: mm,nm,ism,iem,jsm,jem
  character(mtype_len)                  :: mtypem
  character( name_len)                  ::  namem

  err =''
  mm  = m;           nm = m;                                                                                                        ! default: m x m square matrix
  if(present(n))     nm = n

  ism = 1; iem = mm; 
  jsm = 1; jem = nm;                                                                                                                ! default values

  if(mm <= 0 .or. nm <= 0) call matrix_error('matrix_metadata_put: m,n<= 0')

  if(present(is)) then
   ism = is; iem = is + mm - 1
   jsm = is; jem = is + nm - 1
  end if
  if(present(js))then
   jsm = js
   jem = js + nm - 1
  end if

  if(iem - ism /= mm - 1) write(err,*) trim(err) //NEW_LINE('A')// 'MAT bounds incorrect m,is,ie: ',mm,ism,iem
  if(jem - jsm /= nm - 1) write(err,*) trim(err) //NEW_LINE('A')// 'MAT bounds incorrect n,js,je: ',nm,jsm,jem
  if(err       /= ''    ) call matrix_error('matrix_metadata_put: ERROR'//NEW_LINE('A')//trim(err))

  namem = ''
  if(present( name))  namem = name

  write(info,*)     'matrix_metadata_put: info: '//trim(namem)//NEW_LINE('A')//                  &
                    'n,m,is,js,ie,je: ',nm,mm,ism,jsm,iem,jem,  NEW_LINE('A')
  !---------------------------------------------------------------------------------------------------------------------------------
  select type(MAT)
   !--------------------------------------------------------------------------------------------------------------------------------
   class is( Matrix)                                                                                                                ! Matrix, CMatrix
    mtypem= 'ZG';
    if(present(mtype)) mtypem = mtype
    if(mtypem(1:1)/='Z'                       )call matrix_error(trim(info)//'wrong mtype (Z): '//trim(mtypem))
    if(mtypem(2:2)/='G' .and. mtypem(2:2)/='H')call matrix_error(trim(info)//'wrong mtype (Z): '//trim(mtypem))                     ! ZG, ZH available only
    if(mtypem(2:2)=='H' .and. mm         /=nm )call matrix_error(trim(info)//'type (H) must be a square matrix')                    ! ZH: only square matrix allowed
   !--------------------------------------------------------------------------------------------------------------------------------
   type  is(DMatrix)                                                                                                                 
    mtypem= 'DG';
    if(present(mtype)) mtypem = mtype
    if(mtypem(1:1)/='D'                       )call matrix_error(trim(info)//'wrong mtype (D): '//trim(mtypem))
    if(mtypem(2:2)/='G' .and. mtypem(2:2)/='S')call matrix_error(trim(info)//'wrong mtype (D): '//trim(mtypem))                     ! DG, DS available only
    if(mtypem(2:2)=='S' .and. mm         /=nm )call matrix_error(trim(info)//'type (S) must be a square matrix')                    ! DS: only square matrix allowed
   !--------------------------------------------------------------------------------------------------------------------------------
   class default
    call matrix_error('matrix_metadata_put: MAT is of unknown class.')
   end select
  !---------------------------------------------------------------------------------------------------------------------------------

   MAT % m     =  mm; MAT % n  =  nm
   MAT % is    = ism; MAT % ie = iem; 
   MAT % js    = jsm; MAT % je = jem;
   MAT % mtype = mtypem
   MAT % name  =  namem

 end subroutine     matrix_metadata_put
!-----------------------------------------------------------------------------------------------------------------------------------
! Careful: if the metadata of MATA makes sense, it is the responsibility of the user
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    matrix_metadata_copy(MATA,MATB)                                                                                 !:.:.:
  class(MatrixClass), intent(in)        :: MATA
  class(MatrixClass), intent(in out)    :: MATB

  MATB % m     = MATA % m
  MATB % n     = MATA % n
  MATB % is    = MATA % is
  MATB % ie    = MATA % ie
  MATB % js    = MATA % js
  MATB % je    = MATA % je
  MATB % mtype = MATA % mtype
  MATB % name  = MATA % name
  !---------------------------------------------------------------------------------------------------------------------------------
  select type(MATB)
   !--------------------------------------------------------------------------------------------------------------------------------
   class is( Matrix)  
    if( MATB % mtype(1:1) == 'D') MATB % mtype = 'ZG'
   !--------------------------------------------------------------------------------------------------------------------------------
   type  is(DMatrix)   
    if( MATB % mtype(1:1) == 'Z') MATB % mtype = 'DG'       
   end select
 end subroutine     matrix_metadata_copy
!...................................................................................................................................
!.......................... Components   ...........................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call MAT % read(unit=f_mout)        ! read raw values of matrix from unit
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_read                                            (MAT,unit)                                               !:.:.:
  class(MatrixClass)    , intent(in out):: MAT
  integer     , optional, intent(in)    :: unit
  integer                               :: n,m,i,j,is,ie,js,je,un
  real(dp)                              :: x,y

  un    = f_minput;
  if(present(unit)) un     = unit

  m  = MAT % m ; n  = MAT % n
  is = MAT % is; ie = MAT % ie
  js = MAT % js; je = MAT % je
  
  select type (MAT)
   class is (Matrix) ! Matrix, CMatrix
    if( .not. allocated(MAT % v) ) return
    do j = js, je;do i = is, ie
     read (un,*) x, y
     MAT % v(i,j) = CMPLX(x,y,kind=dp)
    end do;        end do
   type is (DMatrix)
    if( .not. allocated(MAT % v) ) return
    do j = js, je;do i = is, ie
     read (un,*) MAT % v(i,j)
    end do;        end do
   class default
    call matrix_error('matrix_read: MAT is of unknown class')
   end select

  end subroutine matrix_read
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call MAT % save(unit=f_mout,fmt='G28.16')                                                                                   ! print raw values of matrix to unit
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_save                                            (MAT,unit,fmt)                                           !:.:.: 
  class(MatrixClass)    , intent(in)    :: MAT
  integer     , optional, intent(in)    :: unit
  character(*), optional, intent(in)    :: fmt
  character(20)                         :: fmtout
  integer                               :: n,m,i,j,is,ie,js,je,un

  un    = f_mout;
  fmtout='G28.17'
  if(present(unit)) un     = unit
  if(present(fmt )) fmtout = fmt
  m  = MAT % m ; n  = MAT % n
  is = MAT % is; ie = MAT % ie
  js = MAT % js; je = MAT % je
  
  select type (MAT)
   class is (Matrix) ! Matrix, CMatrix
    if( .not. allocated(MAT % v) ) return
    do j = js, je;do i = is, ie
     write(un,'(2'//fmtout//')') MAT % v(i,j)
    end do;        end do
   type is (DMatrix)
    if( .not. allocated(MAT % v) ) return
    do j = js, je;do i = is, ie
     write(un,'( '//fmtout//')') MAT % v(i,j)
    end do;        end do
   class default
    call matrix_error('matrix_save: MAT is of unknown class')
   end select

 end subroutine matrix_save
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call MAT % print(unit=f_mout,fmt='G28.17',form='Mathematica'/'Default'/'Info',ips=12,ipe=24,jps=55,jpe=65)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_print                                           (MAT,unit,fmt,form,ips,ipe,jps,jpe)                      !:.:.: 
  class(MatrixClass)    , intent(in)    :: MAT
  integer     , optional, intent(in)    :: unit
  character(*), optional, intent(in)    :: fmt
  character(*), optional, intent(in)    :: form
  integer     , optional, intent(in)    :: ips,ipe,jps,jpe
  character(25)                         :: fmts,fmtn,myform,mname,eol
  character(math_len)                   :: nleft,nright
  character(mtype_len)                  :: mtype
  character( name_len)                  :: name
  character(1000)                       :: message,ibounds
  integer                               :: n,m,i,j,is,ie,js,je,un

  if( MAT % m == 0              ) call matrix_error('matrix_print: MAT not initialized: MAT % m = 0')

  un    = f_mout;   fmtn    = 'G28.17'
  if(present(unit)) un      = unit
  if(present(fmt )) fmtn    = fmt

  m     = MAT % m    ;  n   = MAT % n
  is    = MAT % is   ; ie   = MAT % ie; 
  js    = MAT % js   ; je   = MAT % je
  mtype = MAT % mtype; name = MAT % name

  if(present(ips) .and. ips >= is) is = ips  ! if asked, print only the ips:ipe, jps:jpe section of the arrays
  if(present(ipe) .and. ipe <= ie) ie = ipe
  if(present(jps) .and. jps >= js) js = jps
  if(present(jpe) .and. jpe <= je) je = jpe

  myform = 'default'
  if(present(form)) myform = form

  select case (myform(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('d','D')
   write(un,'(A,I6,A3,I6,A,I6,A3,I6,A2,I6,A3,I6,A2)') &
        '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
        '#pmstart: ' // TRIM(name) // ' =[                    (is a ',m,' x ',n,'  '//mtype//' matrix) (printing bounds: ',    &
        is,' : ',ie,'  ',js,' : ',je,' )'

   select type (MAT)
   !-------------------------------------------------
   class is (Matrix) ! Matrix, CMatrix
    if( .not. allocated(MAT % v) ) call matrix_error('matrix_print: vec % v not allocated')
    write(fmts,'(A,I12,2A)') '(',2*n, TRIM(fmtn),')'
    do i = is, ie
     write(un,fmts) (MAT % v(i,j), j = js, je)
    end do;
   !-------------------------------------------------
   type is (DMatrix)
    if( .not. allocated(MAT % v) ) call matrix_error('matrix_print: vec % v not allocated')
    write(fmts,'(A,I12,2A)') '(',n, TRIM(fmtn),')'
    do i = is, ie;
     write(un,fmts) (MAT % v(i,j), j = js, je)
    end do; 
   !-------------------------------------------------
   class default
    call matrix_error('matrix_print: MAT is of unknown class')
   !-------------------------------------------------
   end select

  write(un,'(A           )') &
       '#pmend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('i','I')
   write(message,'(A,I12,A3,I12,A,I12,A3,I12,A2,I12,A3,I12,A)') &
        '#info: matrix ' // TRIM(name) // ' is a ',m,' x ',n,'  '//mtype//' matrix, bounds: ',    &
        is,' : ',ie,'  ',js,' : ',je,'    :: v(is:ie,js:je) ::'

   select type (MAT)
   !-------------------------------------------------
   class is (Matrix) ! Matrix, CMatrix
    if( allocated(MAT % v) ) &
    write(ibounds,'(I12,A,I12,A,I12,A,I12,A)') &
    lbound(MAT % v,1),' : ',ubound(MAT % v,1),' , ',lbound(MAT % v,2),' : ',ubound(MAT % v,2)
   !-------------------------------------------------
   type is (DMatrix)
    if( allocated(MAT % v) ) &
    write(ibounds,'(I12,A,I12,A,I12,A,I12,A)') &
    lbound(MAT % v,1),' : ',ubound(MAT % v,1),' , ',lbound(MAT % v,2),' : ',ubound(MAT % v,2)
   !-------------------------------------------------
   class default
    call matrix_error('matrix_print: MAT is of unknown class')
   !-------------------------------------------------
   end select

   write(un,'(A           )') trim(message) // trim(ibounds)
  !---------------------------------------------------------------------------------------------------------------------------------
 case ('m','M') ! Mathematica Format
  nleft=Mathematica_nleft; nright=Mathematica_nright ! from matrix_mod_common, a mathematica number in exponential notation must be between nleft and nright
  mname='m'
  if(len_trim(MAT % name) > 0) mname= MAT % name
  write(un,'(A)') trim(mname) // '={'
  select type(MAT)
   !-------------------------------------------------
  class is (Matrix)
   do  i = is, ie
    write(un,'(A)') '{'
    do j = js, je - 1
     write(un,'(A,'//trim(fmtn)//',A,'//trim(fmtn)//',A)')  &
          trim(nleft),REAL(MAT % v(i,j), kind=dp),trim(nright)//' + I ('//trim(nleft),IMAG(MAT % v (i,j)),trim(nright)//'),'
    end do
    eol             = ') },'
    if(i == ie) eol = ') }'
    write (un,'(A,'//trim(fmtn)//',A,'//trim(fmtn)//',A)')  &
         trim(nleft), REAL(MAT % v(i,j), kind=dp),trim(nright)//' + I ('//trim(nleft),IMAG(MAT % v (i,j)),trim(nright)//trim(eol)
   end do !i = is, ie
   !-------------------------------------------------
  type  is (DMatrix)
   do  i = is, ie
    write(un,'(A)') '{'
    do j = js, je - 1
     write(un,'(A,'//trim(fmtn)//                  ',A)')  &
          trim(nleft),     MAT % v(i,j)          ,trim(nright)//','
    end do
    eol             = '},'
    if(i == ie) eol = '}'
    write (un,'(A,'//trim(fmtn)//                  ',A)')  & 
         trim(nleft),      MAT % v(i,j)          ,trim(nright)//trim(eol)
   end do !i = is, ie
   !-------------------------------------------------
  class default
   !-------------------------------------------------
  end select
  write(un,'(A)') '};'
  !---------------------------------------------------------------------------------------------------------------------------------
 case default
  call matrix_error('matrix_print: unknown form for printing MAT')
  !---------------------------------------------------------------------------------------------------------------------------------
 end select

 end subroutine matrix_print
!-----------------------------------------------------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------------------------------------------------
!  Copyright by Konstantinos N. Anagnostopoulos, Physics Department, National Technical University of Athens, 2020
!  konstant@mail.ntua.gr, www.physics.ntua.gr/konstant
!  
!  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
!  published by the Free Software Foundation, version 3 of the License.
!  
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public Liense along with this program. If not, see http://www.gnu.org/licenses
!-----------------------------------------------------------------------------------------------------------------------------------
