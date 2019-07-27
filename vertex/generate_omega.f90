program gen
real(8) :: w1,w2,pi
integer :: i,j,k,op,nw
logical :: symmetry

!symmetry=.true.
symmetry=.false.

op=0
pi=dacos(-1.d0)


open(file='omega_list_path',unit=30)

nw=120
do i=-nw,nw
  if(i==0) cycle
  if(i<0) then
  w1=(2*i+1)
  else
  w1=(2*i-1)
  endif
  w2=w1
  write(30,*) w1,w2
end do

do i=-nw,nw
  if(i==0) cycle
  if(i<0) then
  w1=(2*i+1)
  else
  w1=(2*i-1)
  endif
  w2=-w1
  write(30,*) w1,w2
end do

close(30)
open(file='omega_list_area',unit=30)

nw=80
do i=-nw,nw
 do j=-nw,nw
  if(i==0) cycle
  if(j==0) cycle
  if(i<0) then
  w1=(2*i+1)
  else
  w1=(2*i-1)
  endif
  if(j<0) then
  w2=(2*j+1)
  else
  w2=(2*j-1)
  endif
  if(symmetry.and.w2<-w1) cycle
  write(30,*) w1,w2
 end do
end do

close(30)

end program
