 module liste3
 use genvar
 use common_def
 
  implicit none

  type node
    type(node), pointer :: next, prev
    real(8)             :: a
    integer             :: info
    integer             :: istart=1
    integer             :: itot
  end type


 type(node), pointer, private    :: ll, cur
 REAL(DBL),ALLOCATABLE, PRIVATE  :: rstates(:)
 INTEGER  ,ALLOCATABLE, PRIVATE  :: states(:)
 logical, private                :: verbose=.true.

  contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine kill_list(tot)
  integer :: tot,ii

  cur => ll
  do ii=1,tot
   cur     => cur%next
   write(*,*) 'cur%itot : ', cur%itot
   deallocate(cur%prev)
  enddo
  deallocate(cur)

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

  subroutine initialize_list
   allocate(ll)
   ll%next => ll
   ll%prev => ll
   ll%istart=-1
   ll%itot=0
   cur     => ll
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

  subroutine add_to_list(aa,info)
   real(8) :: aa
   integer :: info
     allocate(cur%next)
     cur%next%a    =  aa
     cur%next%info =  info
     cur%next%itot =  cur%itot+1
     cur%next%prev => cur
     cur%next%next => ll
     cur           => cur%next
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

  subroutine sweep_list
   if(verbose) write(*,*) 'there are ... elements in list : ', cur%itot
   do while(cur%istart==1)
        if(verbose) write(*,*) cur%a
        cur => cur%prev
   enddo
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

 subroutine mpi_reduce_list(vec_tot_out)
  REAL(DBL)              :: vec_tot_out(:)
  REAL(DBL),ALLOCATABLE  :: rtmp(:)
  INTEGER  ,ALLOCATABLE  :: tmp(:)
  INTEGER                :: counti,i,j,kk

     vec_tot_out = 0.d0

     kk     = cur%itot

     if(verbose) write(*,*) 'list contains ... elements : ', kk

     if(allocated(states)) then
      deallocate(states)  
     endif
     if(allocated(rstates)) then
      deallocate(rstates)
     endif

     allocate(rstates(kk))
     allocate(states(kk))

     do i=kk,1,-1
       rstates(i) = cur%a
        states(i) = cur%info
        write(*,*) ' a, info : ', i,kk,cur%a,cur%info
        cur => cur%prev
     enddo   

     call kill_list(kk)

     do i=1,size2
  
     counti=0
     if(i==rank+1) counti =  kk
     if(size2>1.and..not.no_mpi) CALL MPI_BCAST(counti,1,MPI_INTEGER,i-1,MPI_COMM_WORLD,ierr)

     if(counti>0)then
      if(allocated(tmp)) deallocate(tmp)
      if(allocated(rtmp)) deallocate(rtmp)
      allocate(tmp(counti),rtmp(counti))
      if(i==rank+1)then
        tmp(1:counti)   =  states(1:counti)
        rtmp(1:counti)  = rstates(1:counti)
      endif
      if(size2>1.and..not.no_mpi) CALL MPI_BCAST(tmp(1:counti),counti,MPI_INTEGER,i-1,MPI_COMM_WORLD,ierr)
      if(size2>1.and..not.no_mpi) CALL MPI_BCAST(rtmp(1:counti),counti,MPI_DOUBLE_PRECISION,i-1,MPI_COMM_WORLD,ierr)
       do j=1,counti
        vec_tot_out(tmp(j)) = vec_tot_out(tmp(j)) + rtmp(j)
       enddo
      endif
     enddo

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

 subroutine test_list3
 integer :: i
 real(8) :: final_vec(20)

 call initialize_list

 do i=1,10
  write(*,*) 'add an element'
  call add_to_list(dble(i),i+rank*10)
 enddo

 call mpi_reduce_list(final_vec)
 write(*,*) 'final vec : ', final_vec 
 write(*,*) ' kill list ... '
 write(*,*) 'done'

 call initialize_list

 do i=1,10
  write(*,*) 'add an element'
  call add_to_list(dble(i),i+rank*10)
 enddo

 call mpi_reduce_list(final_vec)
 write(*,*) 'final vec : ', final_vec
 
call MPI_FINALIZE(ierr)
 stop



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

 end module
