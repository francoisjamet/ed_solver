MODULE impurity_class
  use common_def

  !$$$$$$$$$$$$$$$$$$$$
  !$$ IMPURITY CLASS $$
  !$$$$$$$$$$$$$$$$$$$$

  use eigen_sector_class 
  use mask_class, only : new_diag

  !use latticerout    , only : web,T1_T2_connect_unitcells
  !use quantum_hamilt , only : Hamiltonian
  !use matrix         , only : diag

  IMPLICIT NONE

INTERFACE swap
 module procedure swapi,swaps,swapr,swapa,swapc,swapqc,swapqr,zswap_,zswap__,swapcs
END INTERFACE

INTERFACE diag
 MODULE PROCEDURE diagr,diagc,diagi,diagrr,diagr_,diagc_,diagc__,diagr__
END INTERFACE

 !---------------------------------------------------!
TYPE arrayr
 real(8),dimension(:),allocatable :: subarray
END TYPE
 !---------------------------------------------------!
TYPE arrayi
 integer,dimension(:),allocatable :: subarray
END TYPE
 !---------------------------------------------------!
TYPE arrayc
 complex(8),dimension(:),allocatable :: subarray
END TYPE
 !---------------------------------------------------!


TYPE symmetry

 ! ======= groupe ponctuel ======== !

 integer,dimension(:,:),allocatable      ::  rot,trans,alltrans,flip,invtrans
 real(8),dimension(:,:),allocatable      ::  ttrans
 integer                                 ::  nmax,nnmax !max numb of sym.
 real(8)                                 ::  BZVEC(3),BZVEC2(3)
 complex(8),dimension(:),allocatable     ::  rotchar,flipchar,transchar
 complex(8)                              ::  chirot,chiflip
 integer                                 ::  ntrans,nrot,nflip
 real(8),dimension(:,:,:),allocatable    ::  rotmat
 real(8)                                 ::  centerrot(3),diagonale(3)
 integer                                 ::  flag(20),flageff(20)
 integer,dimension(:,:,:),allocatable    ::  symtab
 integer,dimension(:,:),allocatable      ::  sympoint,whichk
 type(arrayi),dimension(:),allocatable   ::  tabid
 integer                                 ::  tabidli
 integer,dimension(:),allocatable        ::  tabidcol

 ! ======= groupe ponctuel ======== !

END TYPE

TYPE unitcell
   integer                                   :: q1,q1small,q1mag,ncor,npk
   real(8)                                   :: band_plot_shift
   real(8),dimension(:,:),allocatable        :: insidepos,pk
   real(8),dimension(:,:,:),allocatable      :: clock
   integer,dimension(:,:),allocatable        :: half
   integer,dimension(:),allocatable          :: nneigh
   integer,dimension(:,:,:),allocatable      :: invlink
   real(8)                                   :: a(3),b(3),c(3),as(3),bs(3),cs(3),VBZ
   integer,dimension(:),allocatable          :: smallsite,smallsite2,orbitals
   integer,dimension(:,:),allocatable        :: phase2,phase3
   integer                                   :: struct_nm,struct_nlink
   integer,dimension(:,:),allocatable        :: struct
   logical                                   :: FLAG_MATRIX_ELEMENTS
   complex(8),dimension(:,:,:),allocatable   :: matkxky,matkxkyshift
   real(8), dimension(:,:),allocatable       :: dispersion
   complex(8),dimension(:,:,:),allocatable   :: eigenstates
   integer                                   :: k_zero
   real(8)                                   :: reciproc_base(3,3),P(3,3),T1(3),T2(3),T3(3)
   real(8)                                   :: reciproc_boundary(3,3)
   real(8), dimension(:),allocatable         :: xr,yr,zr,xrpp,yrpp,zrpp
   integer                                   :: nL,nLz
   integer                                   :: nfermiNnumb=2000
   real(8),dimension(:,:),allocatable        :: nfermiN,Akwfree

   real(8),dimension(:),allocatable          :: polygonx,polygony
   integer                                   :: Npoly
   integer,dimension(:,:),allocatable        :: site_connected_to_bath
   integer,dimension(:,:,:),allocatable      :: site_connected_to_cor
   logical,dimension(:),allocatable          :: sitecor
   integer,dimension(:),allocatable          :: reduced_space
   integer                                   :: mbath=0,mdata=0,nmanybody=0
   logical                                   :: FLAG_ONLY_NN
   integer                                   :: FLAG_MAG_FIELD
   integer                                   :: FLAG_HORSDIAG_SIGMA
END TYPE


TYPE web

 type(unitcell)                            :: cell
 type(symmetry)                            :: sym

 !web:
 logical                                   :: symmetries
 real(8), dimension(:,:),allocatable       :: nlink_
 integer                                   :: mini,maxi,Ncu,maxconn
 integer,dimension(:,:,:,:),allocatable    :: struct_link
 real(8),dimension(:,:,:),allocatable      :: struct_pos
 integer                                   :: plaquettecentre(100),plaquetteN
 integer,dimension(:),allocatable          :: latticemapi
 real(8)                                   :: P(3,3),teta(2)
 real(8)                                   :: T1(3),T2(3),T3(3)
 logical                                   :: FLAG_NO_SYM_AT_ALL
 logical                                   :: FLAG_OPEN_BC,FLAG_OPEN_BUT_PERIODIC_LINKED,FLAG_UNIT_CELL_INSIDE,FLAG_PLOT_3D
 real(8),dimension(:,:),allocatable        :: x,xma
 real(8)                                   :: randvec(3)
 integer,dimension(:),allocatable          :: site,site2,site3
 integer                                   :: open1,open2,open3,nL,nLz
 complex(8),dimension(:,:),allocatable     :: phase,phase4
 integer                                   :: centresite
 integer,dimension(:,:),allocatable        :: vecin
 real(8),dimension(:),allocatable          :: latticemap
 real(8),dimension(:,:),allocatable        :: distbord,periodist,periodist_open_bc
 real(8),dimension(:),allocatable          :: distances,angles
 real(8),dimension(:,:),allocatable        :: distancesp
 integer,dimension(:),allocatable          :: bordure,whichbord
 integer,dimension(:,:),allocatable        :: nombrevoisin
 integer                                   :: maxd,maxd_real,maxvoisins,bordcount
 integer,dimension(:,:),allocatable        :: maxdistsite
 integer,dimension(:,:,:),allocatable      :: longvoisin
 real(8)                                   :: centre(3)
 integer,dimension(:,:),allocatable        :: ineigh,cadran,direc,centerCu,centreangle
 integer,dimension(:,:),allocatable        :: links
 integer                                   :: nlinks,N
 integer,dimension(:),allocatable          :: conv_unitcell_site_,centredist,k_inv
 logical,dimension(:),allocatable          :: major_site
 integer                                   :: normalisation,nlink_cell
 complex(8),dimension(:,:),allocatable     :: unitary_transform

 !--------------------------------!
 ! copy of the unitcell variables !
 !--------------------------------!

 real(8), dimension(:,:),allocatable       :: dispersion
 real(8)                                   :: reciproc_base(3,3)
 real(8)                                   :: reciproc_boundary(3,3)
 real(8), dimension(:),allocatable         :: xr,yr,zr,xrpp,yrpp,zrpp
 complex(8),dimension(:,:,:),allocatable   :: matkxky,matkxkyshift
 integer                                   :: q1,q1small,q1mag
 real(8),dimension(:,:),allocatable        :: insidepos
 real(8),dimension(:,:,:),allocatable      :: clock
 integer,dimension(:,:),allocatable        :: half
 integer,dimension(:),allocatable          :: nneigh
 real(8)                                   :: a(3),b(3),c(3)
 integer,dimension(:),allocatable          :: smallsite,smallsite2
 integer,dimension(:,:),allocatable        :: phase2,phase3
 integer                                   :: struct_nm,struct_nlink
 integer,dimension(:,:),allocatable        :: struct

END TYPE


TYPE Hamiltonian
  character(22)                               :: hamiltonian_title
  integer                                     :: RANGE_TETA=1
  logical                                     :: spin_sector_not_coupled,NO_DC=.false.,DC_NO_SAME_ORBITALS=.false.
  integer                                     :: ordertype,ncor,mbath,mdata  ! 1=para, 2=spinupdn, 3=bcs, 4=currents
  character(22)                               :: orderlabel(100)
  complex(8),dimension(:,:),allocatable       :: teta,delta,tetadn,deltap,sigma_hf
  complex(8),dimension(:,:,:,:,:),allocatable :: teta_,tetadn_
  complex(8),dimension(:,:),allocatable       :: delta_mat,teta_mat,teta_mat_dn,teta_mat0
  complex(8),dimension(:,:),allocatable       :: delta_mat_p
  real(8),dimension(:,:),allocatable          :: rrrrsign,rrrrrsign
  real(8),dimension(:),allocatable            :: rrsign,rsign
  real(8),dimension(:),allocatable            :: eps_mat
  logical,dimension(:),allocatable            :: full_proj
  real(8),dimension(:),allocatable            :: eps,dU
  complex(8),dimension(:,:),allocatable       :: epsk
  integer,dimension(:,:),allocatable          :: type
  real(8),dimension(:,:),allocatable          :: Vrep,Jterm,field,field_cell
  integer                                     :: subE,ntype,subDiag,q1
  logical,dimension(:),allocatable            :: diagornot,forbidden
  character(20),dimension(:),allocatable      :: label
  integer                                     :: nH
  character(22)                               :: labH(100)
  real(8)                                     :: min(100),max(100),hund,disorder,a_Angstrom,b_Angstrom,c_Angstrom,dist_plane
  real(8)                                     :: cor1,cor2,cor3,cor4,quanta,temperature
  real(8)                                     :: rrrsign
end type



  TYPE(masked_matrix_type), PUBLIC, SAVE :: Eccc

  REAL(DBL),  PARAMETER, PRIVATE         :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,    PARAMETER, PRIVATE         :: F=.FALSE.,T=.TRUE.
  

  TYPE impurity_type
    !$$$$$$$$$$$$$$$$$$$
    !$$ IMPURITY TYPE $$
    !$$$$$$$$$$$$$$$$$$$
    ! NUMBER OF SITES
    INTEGER :: Nc      = 0
    ! NUMBER OF 1-PARTICLE ORBITALS 
    INTEGER :: norbs   = 0
    ! SIZE OF REDUCED HILBERT SPACE OF THE IMPURITY
    INTEGER :: nstates = 0
    ! RANK OF ORBITALS
    INTEGER, POINTER :: iorb(:,:) => NULL()
    ! IMPURITY QUADRATIC ENERGY 
    TYPE(masked_matrix_type), POINTER :: Ec(:) => NULL()   ! in (site,site) basis for a given spin
    ! IMPURITY QUARTIC ENERGY 
    TYPE(masked_real_matrix_type) :: U  ! in (site,site) basis 
  END TYPE

  PRIVATE :: invmat_jordan

CONTAINS

elemental subroutine swapi(a,b)
integer,intent(inout) :: a,b
integer :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swaps(a,b)
real(4),intent(inout) :: a,b
real(4) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapr(a,b)
real(8),intent(inout) :: a,b
real(8) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapqr(a,b)
real(16),intent(inout) :: a,b
real(16) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapqc(a,b)
complex(16),intent(inout) :: a,b
complex(16) :: c
c=a
a=b
b=c
end subroutine
!*****************
subroutine swapa(a,b)
character*(*),intent(inout) :: a,b
character*20000:: c
integer :: maxab
maxab=LEN_TRIM(a)
if(LEN_TRIM(b)>maxab) maxab=LEN_TRIM(b)
c(1:maxab)=a(1:maxab)
a(1:maxab)=b(1:maxab)
b(1:maxab)=c(1:maxab)
end subroutine
!*****************
elemental subroutine swapc(a,b)
complex(8),intent(inout) :: a,b
complex(8) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapcs(a,b)
complex(4),intent(inout) :: a,b
complex(4) :: c
c=a
a=b
b=c
end subroutine
!*****************
    SUBROUTINE ZSWAP__(N,ZX,INCX,ZY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX(16) ZX(*),ZY(*)
      COMPLEX(16) ZTEMP
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZTEMP = ZX(IX)
          ZX(IX) = ZY(IY)
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
          ZTEMP = ZX(I)
          ZX(I) = ZY(I)
          ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END subroutine

!*****************

      SUBROUTINE ZSWAP_(N,ZX,INCX,ZY,INCY)
      INTEGER     :: INCX,INCY,N
      COMPLEX(8)  :: ZX(*),ZY(*)
      COMPLEX(8)  :: ZTEMP
      INTEGER     :: I,IX,IY

      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZTEMP = ZX(IX)
          ZX(IX) = ZY(IY)
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
          ZTEMP = ZX(I)
          ZX(I) = ZY(I)
          ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END subroutine



SUBROUTINE invmat_jordan(nn,a)
  IMPLICIT NONE
  REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER, DIMENSION(SIZE(a,1))          :: ipiv,indxr,indxc
  INTEGER                                :: nn
  LOGICAL, DIMENSION(SIZE(a,1))          :: lpiv
  REAL(8)                                :: pivinv
  REAL(8), DIMENSION(SIZE(a,1))          :: dumc
  INTEGER, TARGET                        :: irc(2)
  INTEGER                                :: i,l,n
  INTEGER, POINTER                       :: irow,icol

  n=SIZE(a,1)

  irow => irc(1)
  icol => irc(2)

  ipiv=0

  DO i=1,n
     !Main loop over columns to be reduced. 
     lpiv = (ipiv == 0)
     !Begin search for a pivot element. 
     irc=MAXLOC(ABS(a),outerand(lpiv,lpiv))
     ipiv(icol)=ipiv(icol)+1
     IF (ipiv(icol) > 1) STOP 'gaussj:singular matrix (1)'

     !We now have the pivot element, so we interchange
     !rows, if needed, to put the pivot element on the diagonal. The columns
     !are not physically interchanged, only relabeled:
     !indxc(i),the column of the ith pivot element, is the ith column that is
     !reduced, while indxr(i) is the row in which that pivot element was
     !originally located. If indxr(i) = indxc(i) there is an implied column
     !interchange. With this form of bookkeeping, the inverse matrix will be
     !scrambled by
     !columns. 

     IF (irow /= icol) CALL swap(a(irow,:),a(icol,:))

     indxr(i)=irow !We are now ready to divide the pivot row by the pivot element, 
                   !located at irow and icol.
     indxc(i)=icol

     IF (a(icol,icol) == zero) STOP 'gaussj:singular matrix (2)'
     pivinv=one/a(icol,icol)
     a(icol,icol)=CMPLX(one,zero)
     a(icol,:)=a(icol,:)*pivinv
     dumc=a(:,icol)

     !Next, we reduce the rows, except for the pivot one, of course. 
     a(:,icol)     = CMPLX(zero,zero)
     a(icol,icol)  = pivinv
     a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
     a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
  END DO

  !It only remains to unscramble the solution in view of the column
  !interchanges. 
  !We do this by interchanging pairs of columns in the reverse order that the
  !permutation 
  !was built up. 
  DO l=n,1,-1
     CALL swap(a(:,indxr(l)),a(:,indxc(l)))
  END DO

CONTAINS

  FUNCTION outerprod(a,b)
    REAL(8), DIMENSION(:), INTENT(in) :: a,b
    REAL(8), DIMENSION(SIZE(a),SIZE(b)) :: outerprod
    outerprod=SPREAD(a,dim=2,ncopies=SIZE(b))* SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION

 FUNCTION outerand(a,b)
   IMPLICIT NONE
   LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b
   LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand
   outerand = SPREAD(a,dim=2,ncopies=SIZE(b)).AND.SPREAD(b,dim=1,ncopies=SIZE(a))
 END FUNCTION

END SUBROUTINE


     subroutine decomposevec(coef,vec,base)
     implicit none
     real(8),intent(inout)   :: vec(:),coef(:)
     real(8)                 :: matrice(size(vec),size(vec)),base(:,:)
     integer                 :: sizemat
      sizemat=size(vec(:))
      if(sizemat/=size(matrice(1,:))) stop 'error decomposevec : dimensions'
      matrice=transpose(base)
      call invmat_jordan(sizemat,matrice)
      coef=matmul(matrice,vec)
     end subroutine

 function diagr__(mat)
 implicit none
 real(8)::mat(:,:,:,:)
 real(8)::diagr__(size(mat(:,1,1,1)),size(mat(1,:,1,1)),size(mat(1,1,1,:)))
 integer::k,i,j
  do j=1,size(mat(:,1,1,1))
   do i=1,size(mat(1,1,1,:))
     diagr__(j,:,i)=(/(mat(j,k,k,i),k=1,size(mat(1,:,1,1)))/)
   enddo
  enddo
 end function

      !------------------------------------!

 function diagc__(mat)
 implicit none
 complex(8)::mat(:,:,:,:)
 complex(8)::diagc__(size(mat(:,1,1,1)),size(mat(1,:,1,1)),size(mat(1,1,1,:)))
 integer::k,i,j
  do j=1,size(mat(:,1,1,1))
   do i=1,size(mat(1,1,1,:))
     diagc__(j,:,i)=(/(mat(j,k,k,i),k=1,size(mat(1,:,1,1)))/)
   enddo
  enddo
 end function

      !------------------------------------!

 function diagrr(mat)
 implicit none
 real(4)::mat(:,:)
 real(4)::diagrr(size(mat(1,:)))
 integer::k
   diagrr=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function

      !------------------------------------!

 function diagc(mat)
 implicit none
 complex(8)::mat(:,:)
 complex(8)::diagc(size(mat(1,:)))
 integer::k
   diagc=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function

      !------------------------------------!

 function diagi(mat)
 implicit none
 integer::mat(:,:)
 integer::diagi(size(mat(1,:)))
 integer::k
   diagi=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function

      !------------------------------------!

 function diagr(mat)
 implicit none
 real(8)::mat(:,:)
 real(8)::diagr(size(mat(1,:)))
 integer::k
   diagr=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function

      !------------------------------------!

 function diagr_(mat)
 implicit none
 real(8) :: mat(:,:,:)
 real(8) :: diagr_(size(mat(:,1,1)),size(mat(1,1,:)))
 integer :: k,i
  do i=1,size(mat(1,1,:))
   diagr_(:,i)=(/( mat(k,k,i), k=1,size(mat(:,1,1))  )/)
  enddo
 end function

      !------------------------------------!

 function diagc_(mat,i1,i2)
 implicit none
 complex(8) :: mat(:,:,:)
 complex(8) :: diagc_(size(mat(:,1,1)),size(mat(1,1,:)))
 integer    :: k,i,i1,i2
  do i=1,size(mat(1,1,:))
   diagc_(:,i)=(/(mat(k,k,i),k=1,size(mat(:,1,1)))/)
  enddo
 end function


 subroutine T1_T2_connect_unitcells(xx,i,j,n1,n2,n3,cadran)
 implicit none
 type(web)        :: xx
 integer          :: i,j,k,l,m,n1,n2,n3
 real(8)          :: base(3,3),Tdecomp(3),vv(3),aa,bb,cc
 logical,optional :: cadran

  if(.not.present(cadran))then
    base(1,:) = xx%T1
    base(2,:) = xx%T2
    base(3,:) = xx%T3
    vv        = xx%xma(j,:) - xx%xma(i,:)
    call decomposevec(Tdecomp,vv,base)
    n1=NINT(Tdecomp(1))
    n2=NINT(Tdecomp(2))
    n3=NINT(Tdecomp(3))
  else
    call cells(xx%cadran(i,j),aa,bb,cc)
    n1=NINT(aa)
    n2=NINT(bb)
    n3=NINT(cc)
  endif

 return
 end subroutine

subroutine cells(cadran,aa,bb,cc)
 integer   :: cadran
 real(8)   :: aa,bb,cc

 aa=0.d0;bb=0.d0;cc=0.d0

 SELECT CASE (cadran)
 !----------------------------!
 CASE(2)
  aa=1.d0 ;  bb=1.d0
 CASE(3)
  aa=0.d0 ;  bb=1.d0
 CASE(4)
  aa=-1.d0;  bb= 0.d0
 CASE(5)
  aa=0.d0 ;  bb=0.d0
 CASE(6)
  aa=1.d0 ;  bb=0.d0
 CASE(7)
  aa= 0.d0;  bb=-1.d0
 CASE(8)
  aa=-1.d0;  bb=-1.d0
 CASE(9)
  aa= 1.d0;  bb=-1.d0
 CASE(1)
  aa=-1.d0;  bb= 1.d0
 !----------------------------!
 CASE(20)
  aa=1.d0 ;  bb=1.d0  ; cc=-1.d0
 CASE(21)
  aa=0.d0 ;  bb=1.d0  ; cc=-1.d0
 CASE(22)
  aa=-1.d0;  bb=0.d0  ; cc=-1.d0
 CASE(23)
  aa=0.d0 ;  bb=0.d0  ; cc=-1.d0
 CASE(24)
  aa=1.d0 ;  bb=0.d0  ; cc=-1.d0
 CASE(25)
  aa= 0.d0;  bb=-1.d0 ; cc=-1.d0
 CASE(26)
  aa=-1.d0;  bb=-1.d0 ; cc=-1.d0
 CASE(27)
  aa= 1.d0;  bb=-1.d0 ; cc=-1.d0
 CASE(19)
  aa=-1.d0;  bb= 1.d0 ; cc=-1.d0
 !----------------------------!
 CASE(11)
  aa=1.d0 ;  bb=1.d0  ; cc=1.d0
 CASE(12)
  aa=0.d0 ;  bb=1.d0  ; cc=1.d0
 CASE(13)
  aa=-1.d0;  bb= 0.d0 ; cc=1.d0
 CASE(14)
  aa=0.d0 ;  bb=0.d0  ; cc=1.d0
 CASE(15)
  aa=1.d0 ;  bb=0.d0  ; cc=1.d0
 CASE(16)
  aa= 0.d0;  bb=-1.d0 ; cc=1.d0
 CASE(17)
  aa=-1.d0;  bb=-1.d0 ; cc=1.d0
 CASE(18)
  aa= 1.d0;  bb=-1.d0 ; cc=1.d0
 CASE(10)
  aa=-1.d0;  bb= 1.d0 ; cc=1.d0
 !----------------------------!
 END SELECT

 return
 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_impurity(impurity,Nc,IMASKE,IMASKU) 
    TYPE(impurity_type), INTENT(INOUT) :: impurity
    INTEGER,             INTENT(IN)    :: Nc
    INTEGER, OPTIONAL,   INTENT(IN)    :: IMASKU(Nc,Nc),IMASKE(Nc,Nc,2)
    INTEGER                            :: spin

    CALL delete_impurity(impurity)

    ! NUMBER OF SITES
    impurity%Nc    = Nc

    ! NUMBER OF 1-PARTICLE ORBITALS
    impurity%norbs = Nc * 2
   
    ! NUMBER OF IMPURITY STATES
    impurity%nstates = 2**impurity%norbs
   
    ! ORDERING OF ORBITALS WITH INCREASING RANK  = |(site,up)> |(site,down)>

    IF(ASSOCIATED(impurity%iorb)) DEALLOCATE(impurity%iorb,STAT=istati) ; ALLOCATE(impurity%iorb(Nc,2)) 

    CALL ramp(impurity%iorb(:,1))

    impurity%iorb(:,2) = impurity%iorb(:,1) + Nc

    IF(PRESENT(IMASKE))THEN
      ! QUADRATIC ENERGY
      if(ASSOCIATED(impurity%Ec)) DEALLOCATE(impurity%Ec,STAT=istati) ; ALLOCATE(impurity%Ec(SIZE(IMASKE,3)))
      DO spin=1,SIZE(IMASKE,3)
       CALL new_masked_matrix(impurity%Ec(spin),"Ec(sz="//TRIM(cspin(spin))//")",Nc,Nc,IMASK=IMASKE(:,:,spin),IS_HERM=T)
      ENDDO
    ENDIF

    IF(PRESENT(IMASKU))THEN
      ! QUARTIC ENERGY
      CALL new_masked_real_matrix(impurity%U,"U",Nc,Nc,IMASK=IMASKU,IS_HERM=T)
    ENDIF

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_impurity(IMP)
    TYPE(impurity_type), INTENT(INOUT) :: IMP
    INTEGER                            :: spin

    IF(ASSOCIATED(IMP%iorb)) DEALLOCATE(IMP%iorb,STAT=istati) 
    IF(ASSOCIATED(IMP%Ec))THEN
      DO spin=1,SIZE(IMP%Ec)
       CALL delete_masked_matrix(IMP%Ec(spin))
      ENDDO
     DEALLOCATE(IMP%Ec,STAT=istati)
    ENDIF

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE copy_impurity(IMPOUT,IMPIN)

    TYPE(impurity_type), INTENT(INOUT) :: IMPOUT
    TYPE(impurity_type), INTENT(IN)    :: IMPIN
    INTEGER                            :: spin

    IF(.NOT.ASSOCIATED(IMPIN%Ec))  STOP "ERROR IN copy_impurity: INPUT  ISNT ALLOCATED!"
    IF(.NOT.ASSOCIATED(IMPOUT%Ec)) STOP "ERROR IN copy_impurity: OUTPUT ISNT ALLOCATED!"
    IMPOUT%Nc      = IMPIN%Nc
    IMPOUT%norbs   = IMPIN%norbs
    IMPOUT%nstates = IMPIN%nstates
    DO spin=1,SIZE(IMPIN%Ec)
      CALL copy_masked_matrix(IMPOUT%Ec(spin),IMPIN%Ec(spin))
    ENDDO
    CALL copy_masked_real_matrix(IMPOUT%U,IMPIN%U)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE define_impurity(impurity,mmu,impurity_,Himp,Eimp)
  
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ READ IMPURITY PARAMETERS $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    TYPE(impurity_type), INTENT(INOUT) :: impurity
    type(web)                          :: impurity_
    type(hamiltonian)                  :: Himp
    integer                            :: jj,ff,k,i,j,kkk,ijk,ki,kj,kki,kkj,ii
    REAL(DBL)                          :: rval,mmu,val(impurity_%N**2)
    INTEGER                            :: mu,Nc,spin,iind_,iind
    INTEGER, ALLOCATABLE               :: IMASKU(:,:)

#ifdef _complex
    COMPLEX(DBL),OPTIONAL              :: Eimp(:,:)
#else
    REAL(DBL),OPTIONAL                 :: Eimp(:,:)
#endif

    CALL delete_impurity(impurity)
    Nc=impurity_%N

    write(log_unit,*) '......  starting impurity problem with number of site  :  ', Nc
    CALL new_impurity(impurity,Nc)

    !======================!
    ! NON-QUADRATIC ENERGY 
    !======================!

    if(allocated(IMASKU)) deallocate(IMASKU,STAT=istati); ALLOCATE(IMASKU(Nc,Nc))

    IMASKU = 0; kkk = 0 ; val = 0.d0
    
    do jj=1,Nc
      kkk           = kkk + 1
      IMASKU(jj,jj) = kkk
      if(.not.allocated(UUmatrix))then
        val(kkk)      = Himp%dU(impurity_%site(jj))
      else
        val(kkk)      = UUmatrix(jj,jj)
      endif
    enddo

    kkk=jj-1

    do jj=1,Nc
    if(.not.allocated(UUmatrix))then
     do i=1,impurity_%nneigh(impurity_%site(jj))
      if(impurity_%cadran(jj,i)==5)then
        kkk                               = kkk+1
        if(kkk>size(val)) then
         do j=1,Nc
         do ii=1,impurity_%nneigh(impurity_%site(j))
           if(impurity_%cadran(j,ii)==5)then
            write(log_unit,*) ' site, neighbor : ', j, impurity_%ineigh(j,ii)
           endif
         enddo
         enddo
         write(*,*) 'size val   : ', size(val)
         write(*,*) 'impurity%N : ', impurity_%N
         write(*,*) 'Nc         : ', Nc
         stop 'error build impurity Vrep matrix, too much elements, a geometry problem?'
        endif
        IMASKU(jj,impurity_%ineigh(jj,i)) = kkk
        val(kkk)                          = Himp%Vrep(impurity_%site(jj),i)
      endif
     enddo
    else
     do i=jj+1,Nc
       kkk=kkk+1
       if(kkk>size(val)) then
          write(*,*) 'stop error in building impurity'
          stop
       endif
       IMASKU(jj,i) = kkk
       val(kkk)     = UUmatrix(jj,i)
     enddo
    endif
    enddo
 
    CALL new_masked_real_matrix(impurity%U,"U",Nc,Nc,IMASK=IMASKU,IS_HERM=T)

    write(145+rank,*) '===== DIAGONALIZING IMP WITH ONSITE REPULSION ====='
    DO iind=1,impurity%U%MASK%nind
      rval=val(iind)
      CALL fill_masked_real_matrix(impurity%U,iind,rval)
    ENDDO
    call write_array(impurity%U%mat, ' Coulomb repulsion ', unit=145+rank, short=.true.)
    write(145+rank,*) '==================================================='


    ! TEST HERMITICITY
    write(log_unit,*) 'test hermiticity impurity U'
    CALL test_masked_real_matrix_symmetric(impurity%U)

    !======================!
    ! QUADRATIC ENERGY
    !======================!

    if(associated(impurity%Ec)) deallocate(impurity%Ec,STAT=istati); ALLOCATE(impurity%Ec(2))

    call update_impurity(impurity_%N,impurity,mmu,impurity_,Himp,Eimp=Eimp)

    if(allocated(IMASKU)) deallocate(IMASKU,STAT=istati)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE update_impurity(Nc,impurity,mmu,impurity_,Himp,Eimp)
  
    TYPE(impurity_type), INTENT(INOUT) :: impurity
    type(web)                          :: impurity_
    type(hamiltonian)                  :: Himp
    integer                            :: jj,ff,k,i,j

#ifdef _complex
    COMPLEX(DBL)                       :: val(Nc**2*4)
    COMPLEX(DBL),OPTIONAL              :: Eimp(:,:)
#else
    REAL(DBL)                          :: val(Nc**2*4)
    REAL(DBL),OPTIONAL                 :: Eimp(:,:)
#endif

    REAL(DBL)                          :: rval,mmu
    INTEGER                            :: mu,Nc,spin,iind_,iind
    INTEGER                            :: IMASKE(Nc,Nc,2)

   IMASKE=0; k=0; val=0.d0 

   if(.not.present(Eimp))then

    do i=1,impurity_%N
     k=k+1
     if(allocated(Himp%eps)) val(k)=-mmu+Himp%eps(impurity_%site(i))
     IMASKE(i,i,1:2)=k
     if(allocated(impurity_%nneigh))then
     do j=1,impurity_%nneigh(impurity_%site(i))
      if(allocated(impurity_%ineigh).and.impurity_%cadran(i,j)==5)then
       k=k+1
       val(k)=Himp%teta(impurity_%site(i),j)
       IMASKE(i,impurity_%ineigh(i,j),1:2)=k
      endif
     enddo
     endif
    enddo

   else

    if(Nc/=impurity_%N)then
       write(*,*) 'something wrong in impurity class, size Nc'
       stop
    endif
    if(size(Eimp,1)/=2*impurity_%N)then
       write(*,*) 'something wrong in impurity class, size Eimp and Nc do not match'
       write(*,*) 'shape Eimp : ', shape(Eimp)
       write(*,*) 'impurity%N : ', impurity_%N 
      stop
    endif

    do i=1,impurity_%N
     do j=1,impurity_%N
      if(abs(Eimp(i,j))>1.d-7.or.i==j)then
       k=k+1
                val(k) =   Eimp(i,j)
       if(i==j) val(k) =   val(k) - mmu
       IMASKE(i,j,1)=k
      endif
      if(abs(Eimp(j+Nc,i+Nc))>1.d-7.or.i==j)then
       k=k+1
                val(k) = - Eimp(j+Nc,i+Nc)   !   Eimp in the supra form, here we want the TB form
       if(i==j) val(k) =   val(k) - mmu
       IMASKE(i,j,2)=k
      endif
     enddo
    enddo

   endif

    DO spin=1,2
      CALL new_masked_matrix(impurity%Ec(spin),"Ec(sz="//TRIM(cspin(spin))//")",Nc,Nc,IMASK=IMASKE(:,:,spin),IS_HERM=T)
    ENDDO
    CALL clean_redundant_imask(impurity%Ec)

    DO spin=1,SIZE(impurity%Ec)
      do i=1,k
       CALL fill_masked_matrix(impurity%Ec(spin),i,val(i)) 
      enddo
      call write_array( impurity%Ec(spin)%rc%mat , ' Ec(spin) ' , unit=log_unit, short=.true. )
      call write_array( IMASKE(:,:,spin), ' mask spin ', unit=log_unit )
    enddo

    if(present(Eimp))then
      call write_array( Eimp, ' real Eimp ' , unit=log_unit, short=.true.)
    endif

    ! TEST HERMITICITY
    DO spin=1,SIZE(impurity%Ec)
      write(log_unit,*) 'test hermiticity Ec spin ', spin
      CALL test_masked_matrix_hermitic(impurity%Ec(spin))
    ENDDO

    CALL delete_masked_matrix(Eccc)
    CALL Nambu_Ec(Eccc,impurity%Ec)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_impurity(impurity,UNIT)

    TYPE(impurity_type), INTENT(IN) :: impurity
    INTEGER, OPTIONAL,   INTENT(IN) :: UNIT
    INTEGER                         :: unit_,spin

    IF(.NOT.ASSOCIATED(impurity%U%mat)) STOP "ERROR IN write_impurity: INPUT  ISNT ALLOCATED!"

    CALL dump_message(UNIT=UNIT,TEXT="################")
    CALL dump_message(UNIT=UNIT,TEXT="### IMPURITY ###")
    CALL dump_message(UNIT=UNIT,TEXT="################")

                      unit_ = log_unit 
    IF(PRESENT(UNIT)) unit_ = UNIT

    WRITE(unit_,'(a,I0)') "# Nb of sites in the impurity : Nc = ",impurity%Nc

    DO spin=1,SIZE(impurity%Ec)
      write(unit_,*) ' =================================== '
      CALL write_masked_matrix(impurity%Ec(spin),UNIT=UNIT,SHORT=T)
    ENDDO

    write(unit_,*) ' =================================== '
    CALL write_masked_real_matrix(impurity%U,UNIT=UNIT,SHORT=T)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE Nambu_Ec(EcNambu,Ec)

    ! EMBED IMPURITY QUADRATIC ENERGY IN NAMBU MATRICES 

    TYPE(masked_matrix_type), INTENT(INOUT) :: EcNambu
    TYPE(masked_matrix_type), INTENT(IN)    :: Ec(:)
    INTEGER                                 :: Nc ! for clarity only

    IF(SIZE(Ec)==0) STOP "ERROR IN Nambu_Ec: INPUT Ec ISNT ALLOCATED!"

    Nc = Ec(1)%rc%n1
    CALL new_masked_matrix(EcNambu,"EcNambu",Nc*2,Nc*2,IS_HERM=T)

    ! UPPER LEFT BLOCK (SPIN UP)
    EcNambu%rc%mat(   1:Nc,     1:Nc)   =             Ec(1)%rc%mat

    ! LOWER RIGHT BLOCK (SPIN DOWN)
    EcNambu%rc%mat(Nc+1:Nc*2,Nc+1:Nc*2) = - TRANSPOSE(Ec(2)%rc%mat)

    energy_global_shift =  sum(diag( Ec(2)%rc%mat ))

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION average_chem_pot(impurity) 
  TYPE(impurity_type), INTENT(IN) :: impurity
  REAL(DBL)                       :: average_chem_pot
  INTEGER                         :: spin,nspin
  LOGICAL                         :: is_diag(impurity%Nc,impurity%Nc)

    CALL new_diag(is_diag,impurity%Nc)

    ! COMPUTE AVERAGE CHEMICAL POTENTIAL
    nspin = SIZE(impurity%Ec)

    DO spin=1,nspin
     average_chem_pot = average_chem_pot + SUM(impurity%Ec(spin)%rc%mat,is_diag)
    ENDDO
    average_chem_pot = average_chem_pot / ( nspin * impurity%Nc ) 

  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE shift_average_chem_pot(new_chem_pot,impurity) 
    TYPE(impurity_type)                :: impurity
    REAL(DBL),           INTENT(IN)    :: new_chem_pot
    REAL(DBL)                          :: mean_chem_pot
    INTEGER                            :: spin,nspin
    LOGICAL                            :: is_diag(impurity%Nc,impurity%Nc)

    CALL new_diag(is_diag,impurity%Nc)

    ! SHIFT AVERAGE CHEMICAL POTENTIAL TO new_chem_pot
    nspin = SIZE(impurity%Ec)

    ! COMPUTE OLD AVERAGE CHEMICAL POTENTIAL
    mean_chem_pot = average_chem_pot(impurity)

    DO spin=1,nspin
      WHERE(is_diag)
       impurity%Ec(spin)%rc%mat = impurity%Ec(spin)%rc%mat + new_chem_pot - mean_chem_pot
      END WHERE
    ENDDO

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE 
