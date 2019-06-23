MODULE density_matrix

  use HAIMupdo_class


  IMPLICIT NONE

INTERFACE diagonalize
    MODULE PROCEDURE diagonalize_comp,diagonalize_real
END INTERFACE

  REAL(DBL), PARAMETER, PRIVATE :: zero=0.0_DBL,one=1.0_DBL
  LOGICAL,   PARAMETER, PRIVATE :: F=.FALSE.,T=.TRUE.

  !----------------------------------------------------!
  ! COMPUTE THE REDUCED DENSITY MATRIX OF THE IMPURITY !
  !----------------------------------------------------!

CONTAINS

 SUBROUTINE diagonalize_comp(matrix,VALP,VECP,EIGENVAL_ONLY)
    implicit none
  !--------------------------------------------------------------------------------------! 
  ! CONVENIENT WRAPPER OF DSYEVD/ZHEEVD (DIVIDE & CONQUER OF SYMMETRIC/HERMITIAN MATRIX) !
  !--------------------------------------------------------------------------------------!

    COMPLEX(8)                         :: matrix(:,:)
    COMPLEX(8),       INTENT(INOUT)    :: VECP(:,:)
    REAL(8),          INTENT(INOUT)    :: VALP(:)
    LOGICAL,  OPTIONAL, INTENT(IN)     :: EIGENVAL_ONLY
    INTEGER                            :: LRWORK
    CHARACTER(LEN=1)                   :: FLAG
    INTEGER                            :: LWORK,LIWORK,N,info
    complex(8), allocatable            :: work(:)
    real(8), allocatable               :: rwork(:)
    INTEGER, allocatable               :: iwork(:)

    N = size(VALP)

    write(*,*) '......... DIAGONALIZE complex matrix ..........'
    if(any(shape(matrix)-shape(VECP)>0)) then
      write(*,*) 'error diagonalize_comp shape do not match'
      write(*,*) 'shape matrix : ', shape(matrix)
      write(*,*) 'shape VECP   : ', shape(VECP)
      stop 'critical'
    endif
    IF(N/=SIZE(matrix,1).OR.N/=SIZE(matrix,2).OR.N/=SIZE(VECP,1).OR.N/=SIZE(VECP,2)) &
   & STOP "ERROR IN diagonalize_dsyevd: INCONSISTENT DIMENSIONS!"

                        FLAG = 'V'
    IF(PRESENT(EIGENVAL_ONLY))THEN
      IF(EIGENVAL_ONLY) FLAG = 'N'
    ENDIF

    IF(N<=1)THEN
      LWORK  = 1
      LIWORK = 1
      LRWORK = 1
    ELSE
      SELECT CASE(FLAG)
        CASE('N')
          LWORK  = N + 1
          LIWORK = 1
          LRWORK = N
        CASE('V')
          LWORK  = 2*N + N*N + 5
          LRWORK = 5*N + 2*N*N + 5
          LIWORK = 3 + 5*N + 5
      END SELECT
    ENDIF

    ALLOCATE(work(LWORK), rwork(LRWORK), iwork(LIWORK))
    VECP=matrix;
    write(*,*) '.... call ZHEEVD ....',FLAG
    CALL ZHEEVD(FLAG,'U',N,VECP,N,VALP,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,info)
    DEALLOCATE(work, rwork, iwork)

    IF(info>0)THEN
      SELECT CASE(FLAG)
        CASE('N')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO CONVERGE")
          CALL dump_message(TEXT=c2s(i2c(info))//" ELEMENTS OF AN INTERMEDIATE TRIDIAGONAL FORM DIDN'T CONVERGE TO 0!")
        CASE('V')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO COMPUTE AN EIGENVALUE")
          CALL dump_message(TEXT="WHILE WORKING ON SUBMATRIX "//c2s(i2c(info/(N+1)))//" x "//c2s(i2c(MOD(info,N+1))))
      END SELECT
    ENDIF
    IF(info<0) CALL dump_message(TEXT=c2s(i2c(info))//"-TH ARGUMENT HAD AN ILLEGAL VALUE")


  END SUBROUTINE

 SUBROUTINE diagonalize_real(matrix,VALP,VECP,EIGENVAL_ONLY)
    implicit none

    !--------------------------------------------------------------------------------------!
    ! CONVENIENT WRAPPER OF DSYEVD/ZHEEVD (DIVIDE & CONQUER OF SYMMETRIC/HERMITIAN MATRIX) !
    !--------------------------------------------------------------------------------------!

    REAL(8)                         :: matrix(:,:)
    REAL(8),          INTENT(INOUT) :: VECP(:,:)
    REAL(8),          INTENT(INOUT) :: VALP(:)
    LOGICAL,  OPTIONAL, INTENT(IN)  :: EIGENVAL_ONLY
    CHARACTER(LEN=1)                :: FLAG
    INTEGER                         :: LWORK,LIWORK,N,info
    real(8), allocatable            :: work(:)
    INTEGER, allocatable            :: iwork(:)

    N = SIZE(VALP)
    IF(N/=SIZE(matrix,1).OR.N/=SIZE(matrix,2).OR.N/=SIZE(VECP,1).OR.N/=SIZE(VECP,2)) &
        &  STOP "ERROR IN diagonalize_dsyevd: INCONSISTENT DIMENSIONS!"

    if(any(shape(matrix)-shape(VECP)>0)) then
        write(*,*) 'error diagonalize_comp shape do not match'
        write(*,*) 'shape matrix : ', shape(matrix)
        write(*,*) 'shape VECP   : ', shape(VECP)
        stop 'critical'
    endif

                        FLAG = 'V'
    IF(PRESENT(EIGENVAL_ONLY))THEN
      IF(EIGENVAL_ONLY) FLAG = 'N'
    ENDIF

    IF(N<=1)THEN
      LWORK  = 1
      LIWORK = 1
    ELSE
      SELECT CASE(FLAG)
        CASE('N')
          LWORK  = 2 * N + 1
          LIWORK = 1
        CASE('V')
          LWORK  = 2 * N**2 + 6 * N + 1
          LIWORK =            5 * N + 3
      END SELECT
    ENDIF

    ALLOCATE(work(LWORK),iwork(LIWORK))
    VECP=matrix
    write(*,*) '.... call DSYEVD ....', FLAG
    CALL DSYEVD(FLAG,'U',N,VECP,N,VALP,WORK,LWORK,IWORK,LIWORK,info)
    DEALLOCATE(work,iwork)

    IF(info>0)THEN
      SELECT CASE(FLAG)
        CASE('N')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO CONVERGE")
          CALL dump_message(TEXT=c2s(i2c(info))//" ELEMENTS OF AN INTERMEDIATE TRIDIAGONAL FORM DIDN'T CONVERGE TO 0!")
        CASE('V')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO COMPUTE AN EIGENVALUE")
          CALL dump_message(TEXT="WHILE WORKING ON SUBMATRIX "//c2s(i2c(info/(N+1)))//" x "//c2s(i2c(MOD(info,N+1))))
      END SELECT
    ENDIF
    IF(info<0) CALL dump_message(TEXT=c2s(i2c(info))//"-TH ARGUMENT HAD AN ILLEGAL VALUE")


  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE compute_density_matrix(dmat,AIM,beta,GS)

    !-----------------------------------------------------!
    ! COMPUTE THE REDUCED DENSITY MATRIX OF THE IMPURITY  !
    !-----------------------------------------------------!

    TYPE(rcmatrix_type),        INTENT(INOUT) :: dmat
    TYPE(AIM_type),             INTENT(IN)    :: AIM
    REAL(DBL),                  INTENT(IN)    :: beta
    TYPE(eigensectorlist_type), INTENT(IN)    :: GS      ! list of lowest eigenstates
    REAL(DBL)                                 :: boltz,Zpart,E0
    INTEGER                                   :: AIMrank1,AIMrank2,IMPrank1,IMPrank2,BATHrank
    INTEGER                                   :: AIMstate1,AIMstate2
    INTEGER                                   :: nIMPstates,nBATHstates,Nb,Nc
    INTEGER                                   :: start_compute
    INTEGER                                   :: IMPchunk(nproc),IMPstatemin(nproc),IMPstatemax(nproc)
#ifdef _complex 
    COMPLEX(DBL)                              :: coeff1,coeff2 
    COMPLEX(DBL), ALLOCATABLE                 :: dmat_vec(:),dmat_vec_tot(:) 
#else
    REAL(DBL)                                 :: coeff1,coeff2 
    REAL(DBL),    ALLOCATABLE                 :: dmat_vec(:),dmat_vec_tot(:) 
#endif
    INTEGER,      ALLOCATABLE                 :: rankmin(:),rankmax(:),rankchunk(:) 
    INTEGER                                   :: thisrank,isector,ieigen 
    TYPE(eigensector_type), POINTER           :: es    => NULL()
    TYPE(eigen_type),       POINTER           :: eigen => NULL()
    
    CALL dump_message(TEXT="# START COMPUTING THE REDUCED DENSITY MATRIX... ")
    CALL reset_timer(start_compute)
    
    nIMPstates  = AIM%impurity%nstates
    Nc          = AIM%impurity%Nc
    nBATHstates = AIM%bath    %nstates
    Nb          = AIM%bath    %Nb

    CALL split(nIMPstates,IMPstatemin,IMPstatemax,IMPchunk)

    Zpart       = partition(beta,GS)
    E0          = GSenergy(GS)
    dmat%rc     = zero

    !=========================================================================================!

    DO isector=1,GS%nsector

      es => GS%es(isector)

      !----------------------------------------------!
      ! PARSE THE LIST OF EIGENSTATES IN THIS SECTOR !
      !----------------------------------------------!

      DO ieigen=1,es%lowest%neigen

        eigen => es%lowest%eigen(ieigen)

        boltz   = DEXPc(-beta*(eigen%val-E0)) ! Boltzman factor

        !----------------------------------------------!
        ! LOOP OVER MATRIX ELEMENTS WE WANT TO COMPUTE !
        !----------------------------------------------!

!$OMP PARALLEL PRIVATE(IMPrank1,IMPrank2,BATHrank,AIMstate1,AIMrank1,AIMstate2,AIMrank2,coeff1,coeff2)
!$OMP DO
        DO IMPrank1=IMPstatemin(iproc),IMPstatemax(iproc)
          DO IMPrank2=1,IMPrank1

            !---------------------------------------!
            ! TRACE OUT THE BATH DEGREES OF FREEDOM !
            !---------------------------------------!

        !----------------------------------------------------------------------------------------------!
            DO BATHrank=1,nBATHstates
              CALL IMPBATH2AIMstate(AIMstate1,IMPrank1-1,BATHrank-1,Nc,Nb) 
              IF(is_in_sector(AIMstate1,es%sector))THEN
               AIMrank1 = rank_func(AIMstate1,es%sector)
               coeff1   = eigen%vec%rc(AIMrank1)
                IF(coeff1/=zero)THEN
                  CALL IMPBATH2AIMstate(AIMstate2,IMPrank2-1,BATHrank-1,Nc,Nb)
                  IF(is_in_sector(AIMstate2,es%sector))THEN
                   AIMrank2 = rank_func(AIMstate2,es%sector)
                   coeff2   = eigen%vec%rc(AIMrank2)
                    IF(coeff2/=zero)THEN
                      dmat%rc(IMPrank1,IMPrank2) = dmat%rc(IMPrank1,IMPrank2) +      coeff1  * conj(coeff2) * boltz
                      IF(IMPrank1/=IMPrank2) &
                      dmat%rc(IMPrank2,IMPrank1) = dmat%rc(IMPrank2,IMPrank1) + conj(coeff1) *       coeff2  * boltz
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO ! end loop over bath states
        !----------------------------------------------------------------------------------------------!

          ENDDO 
        ENDDO ! loop over impurity states
!$OMP END DO
!$OMP END PARALLEL 
      ENDDO ! loop over eigenstates

    ENDDO ! loop over eigensectors
    !=========================================================================================!


    write(*,*) '............ end main loops density matrix ............'

    if(size2>1.and..not.no_mpi) call mpi_collect_

 
    !-------------------------------------------------------!
    ! RENORMALIZE USING PARTITION FUNCTION                  !
    ! 1 if T=0, 1 also at T>0 if we had all the eigenstates !
    !-------------------------------------------------------!

    dmat%rc = dmat%rc / Zpart 

    CALL timer_fortran(start_compute,"# ... TOOK ")



  contains

  !--------------!
  !--------------!
  !--------------!
  !--------------!

   subroutine mpi_collect_

    write(*,*) 'start collecting MPI chunks'

    ALLOCATE(rankmin(nproc),rankmax(nproc),rankchunk(nproc))
    rankmin   = IMPstatemin * ( IMPstatemin - 1 ) / 2 + 1
    rankmax   = IMPstatemax * ( IMPstatemax + 1 ) / 2
    rankchunk = rankmax - rankmin + 1

    ALLOCATE(dmat_vec(rankchunk(iproc)),dmat_vec_tot(rankmax(nproc)))

    DO IMPrank1=IMPstatemin(iproc),IMPstatemax(iproc)
      thisrank = IMPrank1*(IMPrank1-1)/2
      DO IMPrank2=1,IMPrank1
        dmat_vec_tot(thisrank+IMPrank2)                  = dmat%rc(IMPrank1,IMPrank2)
        dmat_vec    (thisrank+IMPrank2-rankmin(iproc)+1) = dmat%rc(IMPrank1,IMPrank2)
      ENDDO
    ENDDO

write(*,*) 'COLLECTING DATA IN DENSITY MATRIX'

if(size2>1.and..not.no_mpi)then
#ifdef _complex
      CALL MPI_ALLGATHERV( &
      & dmat_vec,rankchunk(iproc),MPI_DOUBLE_COMPLEX,dmat_vec_tot,rankchunk,rankmin-1,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
#else
      CALL MPI_ALLGATHERV( &
      & dmat_vec,rankchunk(iproc),MPI_DOUBLE_PRECISION,dmat_vec_tot,rankchunk,rankmin-1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#endif
endif

    IF(ierr/=0.and.size2>1.and..not.no_mpi) CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)

    DO IMPrank1=1,nIMPstates; IF(IMPrank1<IMPstatemin(iproc).OR.IMPrank1>IMPstatemax(iproc))THEN
      thisrank = IMPrank1*(IMPrank1-1)/2
      dmat%rc(IMPrank1,IMPrank1)   = dmat_vec_tot(thisrank+IMPrank1)
      DO IMPrank2=1,IMPrank1-1
        dmat%rc(IMPrank1,IMPrank2) = dmat_vec_tot(thisrank+IMPrank2)
        dmat%rc(IMPrank2,IMPrank1) = conj(dmat%rc(IMPrank1,IMPrank2))
      ENDDO
    ENDIF; ENDDO

    if(allocated(dmat_vec)    ) deallocate(dmat_vec)
    if(allocated(dmat_vec_tot)) deallocate(dmat_vec_tot)
    if(allocated(rankmin))      deallocate(rankmin)
    if(allocated(rankmax))      deallocate(rankmax)
    if(allocated(rankchunk))    deallocate(rankchunk)

  end subroutine

  !--------------!
  !--------------!
  !--------------!
  !--------------!

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

      SUBROUTINE DSORT (DX, DY, N, KFLAG)

      INTEGER KFLAG, N
      DOUBLE PRECISION DX(*), DY(*)
      DOUBLE PRECISION R, T, TT, TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
      INTEGER IL(21), IU(21)
      EXTERNAL XERMSG
      INTRINSIC ABS, INT
      NN = N
      IF (NN .LT. 1) THEN
       write(*,*) 'The number of values to be sorted is not positive.'
         RETURN
      ENDIF
    KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
     write(*,*)      'The sort control parameter, K, is not 2, 1, -1, or -2.'
      return
      endif

      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            DX(I) = -DX(I)
   10    CONTINUE
      ENDIF

      IF (KK .EQ. 2) GO TO 100

      M = 1
      I = 1
      J = NN
      R = 0.375D0

   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF

   30 K = I

      IJ = I + INT((J-I)*R)
      T = DX(IJ)

      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
      ENDIF
      L = J

      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)

         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
         ENDIF
      ENDIF

   40 L = L-1
      IF (DX(L) .GT. T) GO TO 40

   50 K = K+1
      IF (DX(K) .LT. T) GO TO 50

      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         GO TO 40
      ENDIF
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70

   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1

   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = DX(I+1)
      IF (DX(I) .LE. T) GO TO 80
      K = I

   90 DX(K+1) = DX(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 90
      DX(K+1) = T
      GO TO 80

  100 M = 1
      I = 1
      J = NN
      R = 0.375D0

  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF

  120 K = I

      IJ = I + INT((J-I)*R)
      T = DX(IJ)
      TY = DY(IJ)

      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J

      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)

         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF

  130 L = L-1
      IF (DX(L) .GT. T) GO TO 130


  140 K = K+1
      IF (DX(K) .LT. T) GO TO 140

      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
         GO TO 130
      ENDIF

      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160

  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1

  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = DX(I+1)
      TY = DY(I+1)
      IF (DX(I) .LE. T) GO TO 170
      K = I

  180 DX(K+1) = DX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 180
      DX(K+1) = T
      DY(K+1) = TY
      GO TO 170

  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            DX(I) = -DX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE analyze_density_matrix(dmat,IMPiorb,NAMBU) 

    !--------------------------------------------------------------!
    ! COMPUTE & DISPLAY THE FIRST EIGENPAIRS OF THE DENSITY MATRIX !
    ! COMPUTE THE ENTANGLEMENT ENTROPY:                            !
    ! S = - SUM_i ( lambda_i * log(lambda_i) )                     !
    ! WHERE lambda_i ARE THE EIGENVALUES OF THE DENSITY MATRIX     !
    !--------------------------------------------------------------!

    INTEGER,      INTENT(IN) :: IMPiorb(:,:)
    LOGICAL,      INTENT(IN) :: NAMBU
#ifdef _complex
    COMPLEX(DBL), INTENT(IN) :: dmat(:,:)
    COMPLEX(DBL)             :: VECP(SIZE(dmat,1),SIZE(dmat,1))
#else
    REAL(DBL),    INTENT(IN) :: dmat(:,:)
    REAL(DBL)                :: VECP(SIZE(dmat,1),SIZE(dmat,1))
#endif
    REAL(DBL)                :: VALP(SIZE(dmat,1)),absvec(SIZE(dmat,1)),states(SIZE(dmat,1)) 
    REAL(DBL)                :: ENTROPY,TRACE
    INTEGER                  :: ivp,nvp,icomp
    CHARACTER(LEN=1000)      :: cvec,prefix
    INTEGER                  :: npairs,ncomp,ncompi,i,j
    character(100)           :: intwri

    !-------------------------------------------------------------------------------!
    ! NUMBER OF EIGENVALUES TO DISPLAY, NUMBER OF EIGENVECTOR COMPONENTS TO DISPLAY !
    !-------------------------------------------------------------------------------!

    ncomp=min(size(dmat,1),16) ; npairs=ncomp 

    if(size(dmat,1)<=16) call write_array(dmat," DENSITY MATRIX ", unit=log_unit, short=.true.)

    open(unit=1111,file='density_matrix_before_diag')
    do i=1,size(dmat,1)
     do j=1,size(dmat,2)
      write(1111,*) i,j,real(dmat(i,j))
     enddo
    enddo
    close(1111)
    open(unit=1111,file='density_matrix_ket_basis')
    do i=1,size(dmat,1)
    write(1111,*) cket_from_state(i-1,IMPiorb,NAMBU)
    enddo
    close(1111)


    write(*,*) 'analyse density matrix'
 
    nvp = SIZE(dmat,1)
    CALL diagonalize(dmat,VALP,VECP,EIGENVAL_ONLY=F)
    CALL dump_message(TEXT="# "//c2s(i2c(npairs))//" LARGEST EIGENVALUES OF THE DENSITY MATRIX:") 

    !===============================================================================================!
    DO ivp=nvp,nvp-(npairs-1),-1 ! loop over the npairs largest eigenvalues
      absvec = ABS(VECP(:,ivp)) 
      states = DBLE((/(icomp,icomp=1,nvp)/))
      CALL dsort(absvec,states,nvp,-2) ! sort components of eigenvector in decreasing order and sort state accordingly

      ncompi=ncomp
      do icomp=ncomp,1,-1
       if(abs(VECP(INT(states(icomp)),ivp))>1.d-10)then
        ncompi=icomp
        exit
       endif
      enddo

      write(log_unit,*) '======================================================='
#ifdef _complex
      intwri = '('//c2s(i2c(ncompi))//'(2f9.6,a))'
      write(log_unit,*) 'non zerp elements in state vector  : ', ncompi
      write(log_unit,*)  (VECP(INT(states(icomp)),ivp)," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
      WRITE(cvec,ADJUSTL(TRIM(intwri))) &
      (real(VECP(INT(states(icomp)),ivp)),aimag(VECP(INT(states(icomp)),ivp))," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
#else
      intwri = '('//c2s(i2c(ncompi))//'(f9.6,a))'
      write(log_unit,*) 'non zero elements in state vector  : ', ncompi
      write(log_unit,*)  (VECP(INT(states(icomp)),ivp)," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
      WRITE(cvec,ADJUSTL(TRIM(intwri))) &
      (VECP(INT(states(icomp)),ivp)," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
#endif
      write(log_unit,*) '-------------------------------------------------------'
      prefix = "# "//c2s(i2c(nvp-ivp+1))//"   "
      WRITE(log_unit,'(a'//c2s(i2c(LEN_TRIM(prefix)))//',f9.6,a)') TRIM(prefix),VALP(ivp)," -> "//TRIM(cvec)
      write(log_unit,*) '======================================================='

    ENDDO
    !===============================================================================================!

    write(*,*) '.... start entropy ....'

    write(log_unit,*) '-----------------------------------------'
    write(log_unit,*) 'EIGENVALUES ARE : '
    write(log_unit,'(1000f6.3)') VALP
    write(log_unit,*) '-----------------------------------------'

    ENTROPY = zero; TRACE = zero
    DO ivp = nvp,1,-1
      IF(VALP(ivp)>1.d-16) ENTROPY = ENTROPY - VALP(ivp) * LOG(VALP(ivp))
      TRACE = TRACE + VALP(ivp)
    ENDDO
    write(*,*) '.... done ....'
    WRITE(log_unit,'(a,f9.6)') "# CHECK TRACE          = ",TRACE 
    WRITE(log_unit,'(a,f9.6)') "# ENTANGLEMENT ENTROPY = ",ENTROPY 
    CALL flush(log_unit)

  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
 
END MODULE 
