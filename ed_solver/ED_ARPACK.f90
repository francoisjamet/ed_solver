MODULE ED_ARPACK

  use Lanczos_fast
#ifdef _GPU
  use fortran_cuda
#endif
  IMPLICIT NONE

 !-------------------------------!
 ! DIAGONALIZATION FULL SPETRUM  !
 !-------------------------------!

  REAL(DBL), PARAMETER, PRIVATE  ::  zero=0.0_DBL,one=1.0_DBL

CONTAINS

#ifdef _arpack
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

      subroutine arpack_eigenvector_sym_matrix(verbose,amode,tol,maxn,rvec,values,v,maxncv,nconv_,av,sigma_,mode_)
      implicit none

      !---------------------------------------------!
      !     Solve A*x = lambda*x in regular mode    !
      !---------------------------------------------!

      integer              :: maxn, maxnev, maxncv,nconv_
      Double precision     :: values(maxncv),v(maxn,maxncv), workl(maxncv*(maxncv+8))
      Double precision     :: workd(3*maxn), d(maxncv,2), resid(maxn), ax(maxn)
      logical              :: select(maxncv)
      integer              :: iparam(11), ipntr(11)
      character            :: bmat*1, which*2
      integer              :: ido, n, nev, ncv, lworkl, info, ierr, j, nconv, maxitr, mode, ishfts
      logical              :: rvec
      Double precision     :: tol, sigma,dnrm2
      external             :: dnrm2
      Double precision     :: zero
      parameter        (zero = 0.0D+0)
      character(2)         :: amode
      logical              :: verbose
      real(8),optional     :: sigma_
      integer,optional     :: mode_
     !----------------------------------!
     ! amode : SM around 0              !
     ! amode : BE extremas eigenvalues  !
     !----------------------------------!
      interface
       subroutine av(n,w,v)
        integer          :: n
        Double precision,intent(in)    :: v(n)
        Double precision,intent(inout) :: w(n)
       end subroutine
      end interface


      if(present(sigma_))then
        sigma=sigma_
      else
        sigma=0.d0
      endif

      maxnev=maxncv-1
      n=maxn
      nev=maxnev
      ncv=maxncv

!     %----------------------------------------------------%
!     | A standard eigenvalue                              |
!     | problem is solved (BMAT = 'I'). NEV is the number  |
!     | of eigenvalues to be approximated.  The user can   |
!     | modify NEV, NCV, WHICH to solve problems of        |
!     | different sizes, and to get different parts of the |
!     | spectrum.  However, The following conditions must  |
!     | be satisfied:                                      |
!     |                   N <= MAXN,                       | 
!     |                 NEV <= MAXNEV,                     |
!     |             NEV + 1 <= NCV <= MAXNCV               | 
!     %----------------------------------------------------% 
      if ( n .gt. maxn ) then
         write(*,*) 'n,maxn : ',n,maxn
         print *, ' ERROR with _SDRV1: N is greater than MAXN '
         stop 'arpack error'
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
         stop 'arpack error'
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
         stop 'arpack error'
      end if
      bmat = 'I'
      which = amode
!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%
      lworkl = ncv*(ncv+8)
      info = 0
      ido = 0
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%
      ishfts = 1
      maxitr = 300
      mode   = 1

      if(present(mode_))then
       mode=mode_
       if(mode/=1)then
          write(*,*) 'ARPACK routine needs to be modified for mode/=1';stop
       endif
      else
       mode=1
      endif

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

 10   continue
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, info )
         if (ido .eq. -1 .or. ido .eq. 1) then
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           %--------------------------------------%
            call av (n, workd(ipntr(2)), workd(ipntr(1)))
            go to 10
         end if
      if ( info .lt. 0 ) then
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
      else
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        | Computed eigenvalues may be extracted.    |  
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
         call dseupd ( rvec, 'All', select, d, v, maxn, sigma, bmat, n, which, &
             & nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, ierr )
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
         if ( ierr .ne. 0) then
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
         else
             nconv =  iparam(5)
             do 20 j=1, nconv
!               %---------------------------%
!               | Compute the residual norm |
!               |   ||  A*x - lambda*x ||   |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
                call av(n, ax, v(1,j))
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
 20          continue
            if(verbose)  call dmout(6, nconv, 2, d, maxncv, -6, 'Ritz values and relative residuals')
         end if
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
         if ( info .eq. 1) then
           if(verbose)then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
           endif
         else if ( info .eq. 3) then
           if(verbose)then
            print *, ' '
            print *, ' No shifts could be applied during implicit', ' Arnoldi update, try increasing NCV.'
            print *, ' '
           endif
         end if
        if(verbose)then
         print *, ' '
         print *, ' _SDRV1 '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',   nconv
         print *, ' The number of Implicit Arnoldi update',  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
        endif
      end if
      nconv_=nconv
      values=d(:,1)

 9000 continue

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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

      subroutine arpack_eigenvector_sym_matrix_(verbose,amode,tol,maxn,rvec,values,v,maxncv,nconv_,av,sigma_,mode_)
      implicit none
      !---------------------------------------------!
      !     Solve A*x = lambda*x in regular mode    !
      !---------------------------------------------!
      integer           maxn, maxnev, maxncv, nconv_
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Complex(8)        ax(maxn),d(maxncv),v(maxn,maxncv),workd(3*maxn), workev(3*maxncv), resid(maxn), workl(3*maxncv*maxncv+5*maxncv)
      Double precision  rwork(maxncv), rd(maxncv,3),values(:)
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j, ierr, nconv, maxitr, ishfts, mode
      Complex(8)        sigma
      Double precision  tol
      logical           rvec
      Double precision  dznrm2 , dlapy2
      external          dznrm2 , dlapy2
      character(2)      amode
      logical           verbose
     !----------------------------------!
     ! amode : SM around 0              !
     ! amode : BE extremas eigenvalues  !
     !----------------------------------!

      interface
       subroutine av(n,w,v)
        integer                  :: n
        complex(8),intent(in)    :: v(n)
        complex(8),intent(inout) :: w(n)
       end subroutine
      end interface

      complex(8),optional :: sigma_
      integer,optional    :: mode_


      if(present(sigma_))then
        sigma=sigma_
      else
        sigma=0.d0
      endif

      maxnev=maxncv-1
      n=maxn
      nev=maxnev
      ncv=maxncv

      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if

      bmat  = 'I'
      which = amode
      lworkl  = 3*ncv**2+5*ncv
      ido    = 0
      info   = 0
      ishfts = 1
      maxitr = 300
      mode   = 1

      if(present(mode_))then
       mode=mode_
       if(mode/=1)then
          write(*,*) 'ARPACK routine needs to be modified for mode/=1';stop
       endif
      else
       mode=1
      endif

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

 10   continue

      call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, rwork,info )
      if (ido .eq. -1 .or. ido .eq. 1) then
         call av (n, workd(ipntr(2)), workd(ipntr(1)))
         go to 10
      end if
      if ( info .lt. 0 ) then
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
      else
         call zneupd  (rvec, 'A', select, d, v, maxn, sigma, workev, bmat, n, which, &
            & nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, rwork, ierr)
         if ( ierr .ne. 0) then
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
         else
             nconv = iparam(5)
             do 20 j=1, nconv
                call av(n, ax,v(1,j))
                call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble (d(j))
                rd(j,2) = dimag (d(j))
                rd(j,3) = dznrm2 (n, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
 20          continue
             if(verbose) call dmout (6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and relative residuals')
          end if
        if(verbose)then
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit', ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
         print *, ' '
         print *, '_NDRV1'
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', nconv
         print *, ' The number of Implicit Arnoldi update', ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
        endif
      end if
      nconv_=nconv
      values=real(d(:))

 9000 continue

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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



#endif

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

!subroutine eigenvector_matrix_r(lsize,mat,vaps,eigenvec)
!implicit none
!integer                   :: lsize,i
!real(8)                   :: WORK(3*lsize),RWORK(3*lsize)
!real(8)                   :: mat(lsize,lsize)
!real(8)                   :: eigenvec(lsize,lsize)
!real(8)                   :: rrr
!integer                   :: INFO,order(lsize)
!real(8)                   :: vaps(lsize),temp(1,lsize)
!
!   if(lsize<1) stop 'error eigenvector_matrix, 0 dim'
!
!   eigenvec=mat
!
!   call DSYEV('V','U',lsize,eigenvec,lsize,vaps,WORK,3*lsize,INFO)
!
!   if(INFO/=0)then
!    write(*,*) 'BAD eigenvalue calculation , info = :', INFO
!    write(*,*) 'stop calculations...'
!    stop
!   endif
!
!end subroutine


  SUBROUTINE ARPACK_diago(lowest)
    implicit none
    TYPE(eigenlist_type), INTENT(INOUT) :: lowest
#ifdef _complex
    COMPLEX(DBL),ALLOCATABLE :: VECP(:,:)
#else
    REAL(DBL),ALLOCATABLE    :: VECP(:,:)
#endif
    REAL(DBL), ALLOCATABLE              :: VALP(:)
    INTEGER                             :: start_diagH
    REAL(DBL)                           :: a
    INTEGER                             :: jj,ii,i,j
    REAL(DBL)                           :: coeff
    INTEGER                             :: dimenvec,Neigen_,nconv

    CALL reset_timer(start_diagH)
    if(dimen_H()==0) stop 'error Hilbert space has 0 dimension'
    dimenvec=dimen_H()

#ifdef _complex
    Neigen_=min(dimen_H()-1,  Neigen_arpack)
#else
    Neigen_=min(dimen_H()-1,2*Neigen_arpack)
#endif

    if(Neigen_arpack<2*Neigen)then
      write(*,*) 'WARNING WARNING WARNING : Neigen arpack should be at least twice as large as Neigen, as Arnoldi converges the extremum of the spectra'
    endif

    ALLOCATE(VECP(dimen_H(),Neigen_),VALP(Neigen_))

    !----------------------------------------------------------------------------------------------!
     if(Neigen_>0)then
     write(*,*) '======= RUNNING ARPACK WITH NEIGEN = ', Neigen_ , '========'
#ifdef _arpack
#ifdef _complex
     call arpack_eigenvector_sym_matrix_(.true.,'SR',tolerance,dimenvec,.true.,VALP(1:Neigen_),VECP(1:dimenvec,1:Neigen_),Neigen_,nconv,av=Hmultc)
#else
     call arpack_eigenvector_sym_matrix (.true.,'BE',tolerance,dimenvec,.true.,VALP(1:Neigen_),VECP(1:dimenvec,1:Neigen_),Neigen_,nconv,av=Hmultr)
#endif
#else
     write(*,*) 'error ARPACK library not present, recompile with -D_arpack'
     stop
#endif

#ifndef _complex
     !if(mod(nconv,2)==0)then
     ! nconv=nconv/2
     !else
     ! if(nconv/=1)then
     !  nconv=(nconv+1)/2
     ! endif
     !endif
#endif


     if(nconv>Neigen_) nconv=Neigen_

    if(nconv>0)then
      j=0
      do i=2,nconv
       if(.not.is_eigen_in_window(VALP(i),[VALP(1),VALP(1)+dEmax]))then
        j=i
        exit
       endif
      enddo
      if(j/=0) nconv=j   
     else
      nconv=0
     endif


     write(*,*) 'ARPACK converged ---> : ', nconv
     CALL timer_fortran(start_diagH,"# BUILDING OF "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
     if(nconv>0)then
      CALL reset_timer(start_diagH)
      CALL add_eigen(nconv,VECP,VALP,lowest)
      CALL timer_fortran(start_diagH,"# STORING "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
      if(rank==0)then
       write(*,*) '====================================================='
       write(*,*) '  GETTING THE GROUD STATE                            '
       write(*,*) '  N eigenvalues    : ', lowest%neigen
       write(*,*) '  eigenvalues      : ', lowest%eigen(:)%val
       write(*,*) '  SECTOR dimension : ', lowest%eigen(1)%dim_space
       write(*,*) '====================================================='
      endif
     endif
    else
      nconv=0
     endif
    !----------------------------------------------------------------------------------------------!

    IF(ALLOCATED(VALP)) DEALLOCATE(VALP); IF(ALLOCATED(VECP)) DEALLOCATE(VECP)
 
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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE ED_diago(lowest)
    implicit none
    TYPE(eigenlist_type), INTENT(INOUT) :: lowest
#ifdef _complex
    COMPLEX(DBL), ALLOCATABLE           :: VECP(:,:)
#else
    REAL(DBL),    ALLOCATABLE           :: VECP(:,:)
#endif
    REAL(DBL),    ALLOCATABLE           :: VALP(:)
    TYPE(rcvector_type)                 :: invec,outvec
    INTEGER                             :: start_diagH
    REAL(DBL)                           :: a
    INTEGER                             :: jj,ii,i,j
    REAL(DBL)                           :: coeff
    INTEGER                             :: dimenvec
    
    CALL reset_timer(start_diagH)
    ALLOCATE(VECP(dimen_H(),dimen_H()),VALP(dimen_H()))
    if(dimen_H()==0) stop 'error Hilbert space has 0 dimension'
    dimenvec=dimen_H()
    CALL new_rcvector(invec,dimenvec); CALL new_rcvector(outvec,dimenvec);
    invec%rc=0.d0


    !----------------------------------------------------------------------------------------------!
    DO ii=1,dimenvec
     ! H |v_i> to build H_ij !
      invec%rc(ii)=1.d0
      CALL Hmult__(dimenvec,outvec%rc,invec%rc)
      VECP(:,ii)=outvec%rc
      invec%rc(ii)=0.d0
    ENDDO
 
    CALL timer_fortran(start_diagH,"# BUILDING OF "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
    CALL reset_timer(start_diagH)

    if(dimenvec/=1)then
#ifdef _complex
    call eigenvector_matrix(lsize=dimenvec,mat=VECP,vaps=VALP,eigenvec=VECP)
#else
#ifdef _GPU
    call choose_gpu(rank)
#ifdef _GPU_SINGLEPREC
    call eigenvector_cuda_r(dimenvec,VECP,VALP,VECP,.true.) 
#else
    call eigenvector_cuda_r(dimenvec,VECP,VALP,VECP,.false.)
#endif
#else
   call eigenvector_matrix_r(lsize=dimenvec,mat=VECP,vaps=VALP,eigenvec=VECP)
#endif
#endif
    else
      VALP(1)=VECP(1,1) !BUGGG
      VECP=1.  
    endif

    CALL timer_fortran(start_diagH,"# DIAGONALIZATION OF "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
    CALL reset_timer(start_diagH)

    if(.not.FLAG_FULL_ED_GREEN)then
     j=1
     do i=2,min(dimenvec,Neigen) !size(VALP) !BUG CW
     !write(*,*) 'diff val, dEmax : ', VALP(i)-VALP(1),dEmax,is_eigen_in_window(VALP(i),[VALP(1),VALP(1)+dEmax])
      if(.not.is_eigen_in_window(VALP(i),[VALP(1),VALP(1)+dEmax]))then
       j=i
       goto 333
      endif
     enddo
     j=min(dimenvec,Neigen) !size(VALP)
    else
     j=dimenvec
    endif
333 continue

    !if(rank==0) write(*,*) 'ADDING X EIGENVALUES : ', j

    CALL add_eigen(j,VECP,VALP,lowest)
    !----------------------------------------------------------------------------------------------!

    IF(ALLOCATED(VALP)) DEALLOCATE(VALP); IF(ALLOCATED(VECP)) DEALLOCATE(VECP) 

    CALL delete_rcvector(invec) 
    CALL delete_rcvector(outvec)
    CALL timer_fortran(start_diagH,"# STORING "//TRIM(ADJUSTL(title_H_()))//" TOOK ")

    if(rank==0)then
      write(*,*) '====================================================='
      write(*,*) '  GETTING THE GROUD STATE                            '
      write(*,*) '  N eigenvalues    : ', lowest%neigen
      write(*,*) '  dEmax            : ', dEmax
      write(*,*) '  j                : ', j
     !write(*,*) '  eigenvalues      : ', lowest%eigen(:)%val
      write(*,*) '  SECTOR dimension : ', lowest%eigen(1)%dim_space
      write(*,*) '====================================================='
    endif

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

END MODULE 
