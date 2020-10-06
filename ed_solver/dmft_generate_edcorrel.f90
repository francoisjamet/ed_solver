program generateedcorrel
implicit none
integer                :: ii,i,j,k,l,ien,spins
integer                :: num,len__,status__,paramagnetic
character(20)         :: value
integer,allocatable    :: dd(:,:)
logical                :: check_mask

call get_command_argument(1,value)
read(value,*) k
call get_command_argument(2,value)
read(value,*) spins

 write(*,*) '================ SIZE OF CORREL : ',k
 write(*,*) '================ SPIN           : ',spins

open(unit=140,file='ed_correl1')
allocate(dd(2*k,2*k))
dd=0.0
do i=1,2*k
 dd(i,i)=i
enddo
inquire(file='mask_sym_green_ed',exist=check_mask)
if(check_mask)then
open(unit=2222,file='mask_sym_green_ed')
 do i=1,2*k
  read(2222,*) (dd(i,j),j=1,2*k)
 enddo
close(2222)
endif

write(140,*) '###############################################'
write(140,*) '### ELECTRONIC CORRELATIONS ON THE IMPURITY ###'
write(140,*) '###############################################'
write(140,*) '########################################'
write(140,*) '###     MASK FOR GREENS FUNCTION    ###'
write(140,*) '### [ <c[up]*C[up]>  <c[up]*c[do]> ] ###'
write(140,*) '### [ <C[do]*C[up]>  <c[do]*C[do]> ] ###'
write(140,*) '########################################'
do i=1,2*k
 write(140,'(200i4)') (dd(i,j),j=1,2*k)
enddo
write(140,*) '###########################################'
write(140,*) '### MASKS FOR SPIN/DENSITY CORRELATIONS ###'
write(140,*) '###########################################'
write(140,*) 'F'
write(140,*) 'F'
write(140,*) 'F'


inquire(file='mask_sym_green_ed',exist=check_mask)
if(check_mask)then
   open(unit=2222,file='mask_sym_green_ed')
    do i=1,2*k
     read(2222,*) (dd(i,j),j=1,2*k)
    enddo
   close(2222)
else
   dd=0
   if(k==1.or.spins==0)then
    do i=1,k
     dd(i,i)=i
    enddo
   else
    ii=0
    do i=1,k
     ii=ii+1
     dd(ii,ii)=ii
    enddo
    do i=1,k-1
     do j=i+1,k
       ii=ii+1
       dd(i,j)=ii
       dd(j,i)=ii
     enddo
    enddo
   endif
endif

do i=1,k
 write(140,'(200i4)') (dd(i,j),j=1,k)
enddo

write(140,*) 2.  ,'# wmax = real freq. range [-wmax,wmax] for Sz'
write(140,*) 2.  ,'# wmax = real freq. range [-wmax,wmax] for S+-'
write(140,*) 30. ,'# wmax = real freq. range [-wmax,wmax] for N'
write(140,*) '##################################'
write(140,*) '### MASK FOR Pijk CORRELATIONS ###'
write(140,*) '##################################'
write(140,*) 'F'
write(140,*) 0
write(140,*) '# LIST OF TRIPLETS'
write(140,*) '# MASK OF CORRELATIONS'
write(140,*) 4.
write(140,*) '###################################'
write(140,*) '### MASK FOR Pijkl CORRELATIONS ###'
write(140,*) '###################################'
write(140,*) 'F'
write(140,*) 0
write(140,*)
write(140,*)
write(140,*) 4.
write(140,*) '##################################################'
write(140,*) '### LIST OF PROJECTION VECTORS IN BASIS |.ud2> ###'
write(140,*) '##################################################'
write(140,*) 0

close(140)
end program
