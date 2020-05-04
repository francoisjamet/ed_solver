module vertex

  use correlations
  use namelistmod
  use Lanczos_Cullum_wrapper

PRIVATE

PUBLIC :: four_leg_vertex_matrices_routine

  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE cdagger_copy
  END INTERFACE

  TYPE cdagger_mat
   integer    :: k_p,l_p
   integer    :: k_m,l_m
#ifdef _complex
   complex(8),allocatable :: c_p(:,:) , c_m(:,:)
#else
   real(8),   allocatable :: c_p(:,:) , c_m(:,:)
#endif
  END TYPE

  logical,parameter :: verbose=.false.

contains

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

         !-------------------------!

  subroutine cdagger_copy(cout,cin)
  implicit none
  Type(cdagger_mat),intent(inout) :: cout
  Type(cdagger_mat),intent(in)    :: cin
   call kill_cdagger(cout)
   cout%k_p=cin%k_p
   cout%l_p=cin%l_p
   cout%k_m=cin%k_m
   cout%l_m=cin%l_m
   call allocate_dagger(cout)
   cout%c_p=cin%c_p
   cout%c_m=cin%c_m
  end subroutine

         !-------------------------!

  subroutine kill_cdagger(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
      if(allocated(cdagger%c_p)) deallocate(cdagger%c_p)
      if(allocated(cdagger%c_m)) deallocate(cdagger%c_m)
  end subroutine

         !-------------------------!

  subroutine allocate_dagger(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
   call allocate_dagger_p(cdagger)
   call allocate_dagger_m(cdagger)
  end subroutine

         !-------------------------!

  subroutine allocate_dagger_p(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
    if(allocated(cdagger%c_p)) deallocate(cdagger%c_p)
    allocate(cdagger%c_p(cdagger%k_p,cdagger%l_p))
  end subroutine

         !-------------------------!

  subroutine allocate_dagger_m(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
    if(allocated(cdagger%c_m)) deallocate(cdagger%c_m)
    allocate(cdagger%c_m(cdagger%k_m,cdagger%l_m))
  end subroutine

         !-------------------------!

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

  subroutine four_leg_vertex_matrices_routine(AIM,GS)
  implicit none
  TYPE(eigensectorlist_type)                :: GS
  TYPE(cdagger_mat)                         :: cup(GS%nsector),cdn(GS%nsector)
  TYPE(AIM_type)                            :: AIM
  integer                                   :: kkk,itot,i,j,k,nup,ndn,ii
  integer                                   :: min_up,max_up,min_dn,max_dn
  real(8)                                   :: beta,ZZ

     CALL init_apply_C (AIM)
     itot = AIM%bath%Nb+AIM%impurity%Nc
     ZZ   = partition(beta_ED,GS)
     min_up=10000; max_up=-10000; min_dn=10000; max_dn=-10000;
      do i=1,GS%nsector
        nup=GS%es(i)%sector%updo%up%npart
        ndn=GS%es(i)%sector%updo%down%npart
        if(nup>max_up) max_up=nup; if(ndn>max_dn) max_dn=ndn; if(nup<min_up) min_up=nup; if(ndn<min_dn) min_dn=ndn
      enddo
     write(*,*) 'nup ranging : ',min_up,max_up
     write(*,*) 'nup ranging : ',min_dn,max_dn
     if(verbose) write(*,*) 'writing chiloc'
     if(rank==0) open(unit=1414,file='c_transition',form='unformatted')
     if(Verbose) write(*,*) AIM%impurity%Nc,ZZ,beta_ED,itot,min_up,max_up,min_dn,max_dn,GS%nsector
     if(rank==0) write(1414) AIM%impurity%Nc,ZZ,beta_ED,itot,min_up,max_up,min_dn,max_dn,GS%nsector

     call mpibarrier

     do kkk=1,AIM%impurity%Nc

     call four_leg_vertex_matrices(AIM,GS,Cup_sector,apply_Cup,cup,kkk)
     call mpibarrier
     call four_leg_vertex_matrices(AIM,GS,Cdo_sector,apply_Cdo,cdn,kkk)
     call mpibarrier

     do i=1,GS%nsector
       nup=GS%es(i)%sector%updo%up%npart
       ndn=GS%es(i)%sector%updo%down%npart
       if(FLAG_MPI_GREENS>0) then
         ii=mod(i-1,size2)
         call mpibcast(cup(i)%k_p,iii=ii)
         call mpibcast(cup(i)%l_p,iii=ii)
         call mpibcast(cup(i)%k_m,iii=ii)
         call mpibcast(cup(i)%l_m,iii=ii)
         call mpibcast(cdn(i)%k_p,iii=ii)
         call mpibcast(cdn(i)%l_p,iii=ii)
         call mpibcast(cdn(i)%k_m,iii=ii)
         call mpibcast(cdn(i)%l_m,iii=ii)
         if(rank/=ii)then
          call allocate_dagger_p(cup(i))
          call allocate_dagger_m(cup(i))
          call allocate_dagger_p(cdn(i))
          call allocate_dagger_m(cdn(i))
         endif
         call mpibarrier
         call mpibcast(cup(i)%c_m,iii=ii)
         call mpibcast(cup(i)%c_p,iii=ii)
         call mpibcast(cdn(i)%c_m,iii=ii)
         call mpibcast(cdn(i)%c_p,iii=ii)
       endif
       call mpibarrier
       if(verbose) write(*,*) kkk,nup,ndn
       if(rank==0) write(1414) kkk,nup,ndn
       if(Verbose) write(*,*) GS%es(i)%lowest%neigen
       if(rank==0) write(1414) GS%es(i)%lowest%neigen
       if(rank==0) write(1414) GS%es(i)%lowest%eigen(:)%val
       if(Verbose) write(*,*) 'cup p'
       if(rank==0) write(1414) cup(i)%k_p,cup(i)%l_p
       if(rank==0) write(1414) cup(i)%c_p
       if(Verbose) write(*,*) 'cup '
       if(rank==0) write(1414) cup(i)%k_m,cup(i)%l_m
       if(rank==0) write(1414) cup(i)%c_m
       if(Verbose) write(*,*) 'cdn p'
       if(rank==0) write(1414) cdn(i)%k_p,cdn(i)%l_p
       if(rank==0) write(1414) cdn(i)%c_p
       if(Verbose) write(*,*) 'cdn '
       if(rank==0) write(1414) cdn(i)%k_m,cdn(i)%l_m
       if(rank==0) write(1414) cdn(i)%c_m
     enddo

     do i=1,GS%nsector
       call kill_cdagger(cup(i))
       call kill_cdagger(cdn(i))
     enddo

    enddo

    if(rank==0) close(1414)

  return
  end subroutine

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

  SUBROUTINE four_leg_vertex_matrices(AIM,GS,Asector,applyA,cdagger,kkk)
  implicit none
    TYPE(eigensectorlist_type)                :: GS
    TYPE(cdagger_mat)                         :: cdagger(GS%nsector)
    TYPE(AIM_type)                            :: AIM
    TYPE(sector_type)                         :: Asec
    TYPE(eigensector_type)                    :: Apm_es(2)
    TYPE(eigensector_type), POINTER           :: es    => NULL()
    TYPE(eigen_type),       POINTER           :: eigen => NULL(), OPeigen => NULL()
    INTEGER                                   :: i1,j1,isec_back,ii,i,i_size,j_size,kkk
    INTEGER                                   :: ipm,jpm,kpm,iph,iorb,jorb,jjj1,jjj2,j
    LOGICAL,allocatable                       :: orb(:)
    INTEGER                                   :: isector,ieigen,uup,ddn,itot,i_,v(2),issz
    INTEGER                                   :: ktot,iiorb,jjorb,kk,jj,iorb_f,jorb_f,jsector
    LOGICAL, PARAMETER                        :: linecolumn = .true. !default should be true
    INTEGER                                   :: rank_,size2_
   !----------------------------------------------------!
    INTERFACE
    ! RETURNS SECTOR OF A,A^+|0>
      SUBROUTINE Asector(Asec,pm,sector)
        use sector_class, only : sector_type
       TYPE(sector_type),      INTENT(INOUT) :: Asec
       TYPE(sector_type),      INTENT(IN)    :: sector
       CHARACTER(LEN=1),       INTENT(IN)    :: pm
      END SUBROUTINE
   !----------------------------------------------------!
    ! COMPUTES A,A^+|0>
      SUBROUTINE applyA(Aes,pm,MASK,es,rankr)
        use eigen_sector_class , only : eigensector_type
       TYPE(eigensector_type), INTENT(INOUT) :: Aes
       CHARACTER(LEN=1),       INTENT(IN)    :: pm
       TYPE(eigensector_type), INTENT(IN)    :: es
       INTEGER,                INTENT(IN)    :: rankr
       LOGICAL,                INTENT(IN)    :: MASK(:)
      END SUBROUTINE
    END INTERFACE
   !----------------------------------------------------!

    allocate(orb(AIM%impurity%Nc))
    itot = AIM%bath%Nb+AIM%impurity%Nc

!=====================================================================!
!=====================================================================!
!=====================================================================!

    orb=.false.;orb(kkk)=.true.
    write(*,*) 'ORBITAL MASK : ', orb

    if(FLAG_MPI_GREENS>0)then
      rank_=rank; size2_=size2
     else
      rank_=0; size2_=1
     endif

    DO isector=rank_+1,GS%nsector,size2_
         write(*,*) 'RANK SECTOR ',rank,isector
         es   => GS%es(isector)
         uup  =  GS%es(isector)%sector%updo%up%npart
         ddn  =  GS%es(isector)%sector%updo%down%npart

         !write(*,*) 'es,uup,ddn : ', isector,uup,ddn

         if(verbose) write(*,*) 'BUILDING NEW A SECTORS'
         DO i=1,2
             CALL Asector(Asec,pm(i),es%sector)
             CALL new_eigensector(Apm_es(i),Asec)
         ENDDO

     !----------------------------------------------------!
         DO kpm=1,2

            if(Apm_es(kpm)%sector%updo%not_physical) then
             if(kpm==1)then
              cdagger(isector)%k_p=1
              cdagger(isector)%l_p=1
              call allocate_dagger_p(cdagger(isector))
              cdagger(isector)%c_p=0.d0
             else
              cdagger(isector)%k_m=1
              cdagger(isector)%l_m=1
              call allocate_dagger_m(cdagger(isector))
              cdagger(isector)%c_m=0.d0
             endif
             cycle
            endif

            if(kpm==1.and.verbose) write(*,*) '-------------C^+----------------'
            if(kpm==2.and.verbose) write(*,*) '--------------C-----------------'

            CALL delete_eigenlist(Apm_es(kpm)%lowest)
            if(verbose) write(*,*) 'starting loop over states in sector, there are [x] states : ', GS%es(isector)%lowest%neigen

            DO ieigen=1,GS%es(isector)%lowest%neigen
               if(verbose) write(*,*) 'number of states,ieigen: ', ieigen,GS%es(isector)%lowest%neigen
               if(verbose) write(*,*) 'orb,ieigen,kpm,pm(kpm)   ', orb,ieigen,kpm,pm(kpm)
               if(verbose) write(*,*) 'orbital number           ', kkk
               if(.not.orb(kkk))Then
                 write(*,*) 'orb is false : ', orb
                 write(*,*) 'orb number   : ', kkk
                 stop
               endif
               CALL applyA(Apm_es(kpm),pm(kpm),orb,es,ieigen)
               if(verbose) write(*,*) 'applying A, neigen = ',Apm_es(kpm)%lowest%neigen
            ENDDO

           if(Apm_es(kpm)%lowest%neigen==0)then
            write(*,*) 'error zero eigenstates : ', Apm_es(kpm)%lowest%neigen
            stop
           endif

            jsector=0
            do ii=1,GS%nsector
              if(GS%es(ii)%sector%updo%up%npart   == Apm_es(kpm)%sector%updo%up%npart.and. &
               & GS%es(ii)%sector%updo%down%npart == Apm_es(kpm)%sector%updo%down%npart)then
               jsector=ii
              exit
             endif
            enddo

            if(jsector==0)Then
              write(*,*) 'missing sector in four leg Chi : ', jsector
              stop
            endif

            i_size=GS%es(jsector)%lowest%neigen
            j_size=GS%es(isector)%lowest%neigen

            if(verbose) then
               write(*,*) 'applyA is in sector and max is         : ',ii,GS%nsector
               write(*,*) 'in target sector there are [x] states  : ',GS%es(ii)%lowest%neigen
               write(*,*) 'up,dn in original sector               : ',GS%es(isector)%sector%updo%up%npart,GS%es(isector)%sector%updo%down%npart
               write(*,*) 'up,dn in target sector                 : ',GS%es(ii)%sector%updo%up%npart,GS%es(ii)%sector%updo%down%npart
               write(*,*) 'there are [x] collected states from A  : ',Apm_es(kpm)%lowest%neigen
               write(*,*) 'up,dn in Apm target sector             : ',Apm_es(kpm)%sector%updo%up%npart,GS%es(ii)%sector%updo%down%npart
               write(*,*) 'destroy or creates                     : ',kpm,pm(kpm)
               write(*,*) 'target/initial sector sizes            : ',i_size,j_size
               write(*,*) 'size sector GS                         : ',size(GS%es(ii)%lowest%eigen(1)%vec%rc)
               write(*,*) 'size sector Apm                        : ',size(Apm_es(kpm)%lowest%eigen(1)%vec%rc)
               write(*,*) 'dim space Gs                           : ',GS%es(ii)%lowest%eigen(1)%dim_space
               write(*,*) 'dim space Apm                          : ',Apm_es(kpm)%lowest%eigen(1)%dim_space
               write(*,*) 'Gs lanczos iter                        : ',GS%es(ii)%lowest%eigen(1)%lanczos_iter
               write(*,*) 'Apm lanczos iter                       : ',Apm_es(kpm)%lowest%eigen(1)%lanczos_iter

            endif

            if(Apm_es(kpm)%lowest%neigen/=GS%es(isector)%lowest%neigen)then
             write(*,*) 'obtained c|i> and initial sector size : ',Apm_es(kpm)%lowest%neigen,GS%es(isector)%lowest%neigen
             write(*,*) 'they should be the same, error'
             stop
            endif

            if(kpm==1)then

             if(linecolumn)then
              cdagger(isector)%k_p=i_size
              cdagger(isector)%l_p=j_size
             else
              cdagger(isector)%k_p=j_size
              cdagger(isector)%l_p=i_size
             endif

              call allocate_dagger_p(cdagger(isector))
              if(verbose) write(*,*) 'cdagger matrix, kpm:', kpm
              do i=1,i_size
               do j=1,j_size
                if(verbose) write(*,*) i,j,size(GS%es(ii)%lowest%eigen(i)%vec%rc)
                if(verbose) write(*,*) i,j,size(Apm_es(kpm)%lowest%eigen(j)%vec%rc)
                if(linecolumn)then
                 cdagger(isector)%c_p(i,j)= scalprod( GS%es(ii)%lowest%eigen(i)%vec%rc , Apm_es(kpm)%lowest%eigen(j)%vec%rc)
                else
                 cdagger(isector)%c_p(j,i)= scalprod( GS%es(ii)%lowest%eigen(i)%vec%rc , Apm_es(kpm)%lowest%eigen(j)%vec%rc)
                endif
               enddo
              enddo
              if(verbose) write(*,*) 'maxval(cdagger) : ', maxval(abs(cdagger(isector)%c_p))
            else

             if(linecolumn)then
              cdagger(isector)%k_m=i_size
              cdagger(isector)%l_m=j_size
             else
              cdagger(isector)%k_m=j_size
              cdagger(isector)%l_m=i_size
             endif

              call allocate_dagger_m(cdagger(isector))
              if(verbose) write(*,*) 'cdagger matrix, kpm:', kpm
              do i=1,i_size
               do j=1,j_size
                if(verbose) write(*,*) i,j,size(GS%es(ii)%lowest%eigen(i)%vec%rc)
                if(verbose) write(*,*) i,j,size(Apm_es(kpm)%lowest%eigen(j)%vec%rc)
                if(linecolumn)Then
                 cdagger(isector)%c_m(i,j)= scalprod( GS%es(ii)%lowest%eigen(i)%vec%rc , Apm_es(kpm)%lowest%eigen(j)%vec%rc)
                else
                 cdagger(isector)%c_m(j,i)= scalprod( GS%es(ii)%lowest%eigen(i)%vec%rc , Apm_es(kpm)%lowest%eigen(j)%vec%rc)
                endif
               enddo
              enddo

              if(verbose) write(*,*) 'maxval(cdagger) : ', maxval(abs(cdagger(isector)%c_m))
            endif

         ENDDO
     !----------------------------------------------------!

    CALL delete_sector(Asec)
    DO iph=1,2
      CALL delete_eigensector(Apm_es(iph))
    ENDDO

    ENDDO
!=====================================================================!
!=====================================================================!
!=====================================================================!

    deallocate(orb)

    write(*,*) 'RANK DONE FOUR LEG : ', rank
    call mpibarrier

return

contains

real(8) function scalprod(v1,v2)
implicit none
real(8) :: v1(:),v2(:)
integer::i
 scalprod = 0.d0
 do i=1,size(v1)
  scalprod=scalprod+v1(i)*v2(i)
 enddo
end function

end subroutine

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

end module
