module splines

  use genvar
  use random, only: drand1

 !===============================================================!
   type spline
     real,dimension(:),allocatable   :: om,Frc,w,cc,tt,smooth
     real,dimension(:,:),allocatable :: der
     integer                         :: nn,kk,N,initialized=0
     integer                         :: order,nest,ierr
   end type
 !===============================================================!

 !---------------------------------------------------!
  INTERFACE resampleit
      MODULE PROCEDURE resampleit___b,resampleit__b
  END INTERFACE
 !---------------------------------------------------!

INTERFACE qsort_adj_array
 MODULE PROCEDURE qsort_adj_array_c,qsort_adj_array_r,qsort_adj_array_c_s,qsort_adj_array_rs &
               & ,qsort_adj_array_veccs,qsort_adj_array_vecrs,qsort_adj_array_vecc,qsort_adj_array_vecr,&
               &  qsort_adj_array_i
END INTERFACE

INTERFACE qsort_array
 MODULE PROCEDURE qsort_array_r,qsort_array_rs,qsort_array_i
END INTERFACE

contains

function compare(f,g)
implicit none
real(8) :: f,g
integer :: compare
 if(f<g) then
  compare=-1
 else
  compare=1
 endif
end function

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
integer function qsort_rand( lower, upper )
implicit none
    integer                       :: lower, upper
    real(4)                       :: r
    r=drand1()
    qsort_rand =  lower + nint(r * (upper-lower))
end function qsort_rand
!***********************************************
!***********************************************
!***********************************************
!***********************************************

function compare_r(f,g)
real :: f,g
integer :: compare_r
 if(f<g) then
  compare_r=-1
 else
  compare_r=1
 endif
end function

recursive subroutine qsort_sort( array, order, left, right )
implicit none
    real(8), dimension(:)         :: array
    integer, dimension(:)         :: order
    integer                       :: left
    integer                       :: right
    integer                       :: i
    integer                       :: last
    if ( left .ge. right ) return
    call qsort_swap( order, left, qsort_rand(left,right) )
    last = left
    do i = left+1, right
        if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
        endif
    enddo
    call qsort_swap( order, left, last )
    call qsort_sort( array, order, left, last-1 )
    call qsort_sort( array, order, last+1, right )

end subroutine qsort_sort

subroutine qsort_swap( order, first, second )
implicit none
    integer, dimension(:)         :: order
    integer                       :: first, second
    integer                       :: tmp
    tmp           = order(first)
    order(first)  = order(second)
    order(second) = tmp
end subroutine

recursive subroutine qsort_sort_r( array, order, left, right )
implicit none
    real,dimension(:)             :: array
    integer, dimension(:)         :: order
    integer                       :: left
    integer                       :: right
    integer                       :: i
    integer                       :: last
    if ( left .ge. right ) return
    call qsort_swap( order, left, qsort_rand(left,right) )
    last = left
    do i = left+1, right
        if ( compare_r(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
        endif
    enddo
    call qsort_swap( order, left, last )
    call qsort_sort_r( array, order, left, last-1 )
    call qsort_sort_r( array, order, last+1, right )
end subroutine

subroutine qsort_adj_array_i(array,order)
implicit none
   integer, dimension(:)            :: array
   integer, dimension(size(array))  :: backup
   integer, dimension(size(array))  :: order
   integer                          :: i,j
   backup=array
   do j=1,size(array)
     array(j)=backup(order(j))
   enddo
end subroutine

           !-------------!

subroutine qsort_adj_array_veccs(array,order)
implicit none
   complex(4), dimension(:,:)                          :: array
   complex(4), dimension(size(array,1),size(array,2))  :: backup
   integer, dimension(size(array,2))                   :: order
   integer                                             :: i,j
   backup=array
   do j=1,size(array,2)
     array(:,j)=backup(:,order(j))
   enddo
end subroutine

           !-------------!

subroutine qsort_adj_array_vecrs(array,order)
implicit none
   real(4), dimension(:,:)                          :: array
   real(4), dimension(size(array,1),size(array,2))  :: backup
   integer, dimension(size(array,2))                :: order
   integer                                          :: i,j
   backup=array
   do j=1,size(array,2)
     array(:,j)=backup(:,order(j))
   enddo
end subroutine

           !-------------!

subroutine qsort_adj_array_vecc(array,order)
implicit none
   complex(8), dimension(:,:)                          :: array
   complex(8), dimension(size(array,1),size(array,2))  :: backup
   integer, dimension(size(array,2))                :: order
   integer                                          :: i,j
   backup=array
   do j=1,size(array,2)
     array(:,j)=backup(:,order(j))
   enddo
end subroutine

           !-------------!

subroutine qsort_adj_array_vecr(array,order)
implicit none
   real(8), dimension(:,:)                          :: array
   real(8), dimension(size(array,1),size(array,2))  :: backup
   integer, dimension(size(array,2))                :: order
   integer                                          :: i,j
   backup=array
   do j=1,size(array,2)
     array(:,j)=backup(:,order(j))
   enddo
end subroutine

subroutine qsort_adj_array_r(array,order)
implicit none
   real(8), dimension(:)            :: array
   real(8), dimension(size(array))  :: backup
   integer, dimension(size(array)) :: order
   integer                         :: i,j
   backup=array
   do j=1,size(array)
     array(j)=backup(order(j))
   enddo
end subroutine

           !-------------!

subroutine qsort_adj_array_rs(array,order)
implicit none
   real, dimension(:)            :: array
   real, dimension(size(array))  :: backup
   integer, dimension(size(array)) :: order
   integer                         :: i,j
   backup=array
   do j=1,size(array)
     array(j)=backup(order(j))
   enddo
end subroutine

           !-------------!

subroutine qsort_adj_array_c(array,order)
implicit none
   complex(8), dimension(:)            :: array
   complex(8), dimension(size(array))  :: backup
   integer, dimension(size(array))     :: order
   integer                             :: i,j
   backup=array
   do j=1,size(array)
     array(j)=backup(order(j))
   enddo
end subroutine

           !-------------!

subroutine qsort_adj_array_c_s(array,order)
implicit none
   complex, dimension(:)            :: array
   complex, dimension(size(array))  :: backup
   integer, dimension(size(array))  :: order
   integer                          :: i,j
   backup=array
   do j=1,size(array)
     array(j)=backup(order(j))
   enddo
end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

subroutine qsort_array_r(array,order2)
implicit none
    real(8),dimension(:)                    :: array
    real(8),dimension(size(array))          :: backup
    integer,dimension(size(array))          :: order
    integer,dimension(size(array)),optional :: order2
    integer                                 :: i
    do i=1,size(order)
      order(i)=i
    enddo
    call qsort_sort( array, order, 1, size(array) )
    do i=1,size(order)
       backup(i)=array(order(i))
    enddo
    array=backup
    if(present(order2)) order2=order
end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

subroutine qsort_array_rs(array,order2)
implicit none
    real,dimension(:)                       :: array
    real,dimension(size(array))             :: backup
    integer,dimension(size(array))          :: order
    integer,dimension(size(array)),optional :: order2
    integer                                 :: i
    do i=1,size(order)
      order(i)=i
    enddo
    call qsort_sort_r( array, order, 1, size(array) )
    do i=1,size(order)
       backup(i)=array(order(i))
    enddo
    array=backup
    if(present(order2)) order2=order
end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

subroutine qsort_array_i(array_,order2)
implicit none
    integer,dimension(:)                    :: array_
    real,dimension(size(array_))            :: array
    real,dimension(size(array))             :: backup
    integer,dimension(size(array))          :: order
    integer,dimension(size(array)),optional :: order2
    integer                                 :: i
    do i=1,size(order)
      order(i)=i
    enddo
    array=float(array_)
    call qsort_sort_r( array, order, 1, size(array) )
    do i=1,size(order)
       backup(i)=array(order(i))
    enddo
    array=backup
    if(present(order2)) order2=order
    array_=NINT(array)
end subroutine




subroutine group_data_rrr(n1,array,array2,array2b,lbin)
implicit none
integer     :: n1
real(8)     :: array(n1),dist,temp,good
real(8)     :: array2(n1),array2b(n1)
integer     :: lbin,i,j,k,l,m,n,siz,count
real(8)     :: mean(n1)
integer     :: howmany(n1)

   siz=n1; temp=0.;count=0; good=0.d0; howmany=0; mean=0.

   do i=1,siz
    temp=array(i); dist=abs(temp-good)
    if(dist>3.d-3.or.count==0)then
     count=count+1
     howmany(count)=1
     mean(count)=array2(i)
     array2b(count)=array(i)
     good=temp
    else
     howmany(count)=howmany(count)+1
     mean(count)=mean(count)+array2(i)
    endif
   enddo

   lbin=count
   do i=1,lbin
   temp=dble(howmany(i))
    if(abs(temp)>1.d-4) then
     array2(i)=mean(i)/temp
    else
     array2(i)=0.
    endif
   enddo

return
end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

      SUBROUTINE CUBSPL ( TAU, C, N, IBCBEG, IBCEND)
!     N = NUMBER OF DATA POINTS. ASSUMED TO BE .GE. 2.
!     (TAU(I), C(1,I), I=1,...,N) = ABSCISSAE AND ORDINATES OF THE
!     DATA POINTS. TAU IS ASSUMED TO BE STRICTLY INCREASING.
!     IBCBEG, IBCEND = BOUNDARY CONDITION INDICATORS, AND
!     C(2,1), C(2,N) = BOUNDARY CONDITION INFORMATION. SPECIFICALLY,
!        IBCBEG = 0  MEANS NO BOUNDARY CONDITION AT TAU(1) IS GIVEN.
!           IN THIS CASE, THE NOT-A-KNOT CONDITION IS USED, I.E. THE
!           JUMP IN THE THIRD DERIVATIVE ACROSS TAU(2) IS FORCED TO
!           ZERO, THUS THE FIRST AND THE SECOND CUBI! POLYNOMIAL PIECES
!           ARE MADE TO COINCIDE.)
!        IBCBEG = 1  MEANS THAT THE SLOPE AT TAU(1) IS MADE TO EQUAL
!           C(2,1), SUPPLIED BY INPUT.
!        IBCBEG = 2  MEANS THAT YHE SECOND DERIVATIVE TAU(1) IS 
!           MADE TO EQUAL C(2,1), SUPPLIED BY INPUT.
!        IBCEND = 0, 1, OR 2 HAS ANALOGOUS MEANING CONCERNING THE
!           BOUNDARY CONDITION AT TAU(N), WITH THE ADDITIONAL INFOR-
!           MATION TAKEN FROM C(2,N).
!     C(J,I), J=1,...,4; I=1,...,L (= N-1) = THE POLYNOMIAL COEFFICIENTS 
!        OF THE CUBI! INTERPOLATING SPLINE WITH INTERIOR KNOTS (OR
!        JOINTS) TAU(2),...,TAU(N-1). PRECISELY, IN THE 
!        INTERVAL (TAU(I), TAU(I+1)), THE SPLINE F IS GIVEN BY
!           F(X) = C(1,I)+H*(C(2,I)+H*(C(3,I)+H*C(4,I)/3.)/2.)
!        WHERE H = X - TAU(I). THE FUNCTION PROGRAM *PPVALU* MAY BE 
!        USED TO EVALUATE F OR ITS DERIVATIVES FROM TAU, C, L = N-1,
!        AND K=4.
      IMPLICIT NONE
      INTEGER IBCBEG, IBCEND, N, I, J, L, M
      DOUBLE PRECISION C(4,N), TAU(N), DIVDF1, DIVDF3, DTAU, G
!     A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SPLOPES S(I) OF
!     F AT TAU(I), I=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS 
!     ELIMINATION, WITH S(I) ENDING UP IN C(2,I), ALL I.
!     C(3,.) AND C(4,.) ARE USED INITIALLY FOR TEMPORARY STORAGE.
      L = N - 1
!     COMPUTE FIRST DIFFERENCES OF TAU SEQUENCE AND STORE IN C(3,.). ALSO,
!     COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN C(4,.).
      DO 10 M = 2, N
         C(3,M) = TAU(M) - TAU(M-1)
   10    C(4,M) = (C(1,M) - C(1,M-1))/C(3,M)
!     CONSTRUCT FIRST EQUATION FROM THE BOUNDARY CONDITION, OF THE FORM
!             C(4,1)*S(1) + C(3,1)*S(2) = C(2,1)
      IF (IBCBEG-1)                     11,15,16
   11 IF (N .GT. 2)                     GO TO 12
!     NO CONDITION AT LEFT END AND N = 2.
      C(4,1) = 1.D0
      C(3,1) = 1.D0
      C(2,1) = 2.D0*C(4,2)
                                        GO TO 25
!     NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
   12 C(4,1) = C(3,3)
      C(3,1) = C(3,2) + C(3,3)
      C(2,1) = ((C(3,2)+2.D0*C(3,1))*C(4,2)*C(3,3)+C(3,2)**2*C(4,3))/C(3,1)
                                        GO TO 19
!     SLOPE PRESCRIBED AT LEFT END.
   15 C(4,1) = 1.D0
      C(3,1) = 0.D0
                                        GO TO 18
!     SECOND DERIVATIVE PRESCRIBED AT LEFT END.
   16 C(4,1) = 2.D0
      C(3,1) = 1.D0
      C(2,1) = 3.D0*C(4,2) - C(3,2)/2.D0*C(2,1)
   18 IF (N .EQ. 2)                     GO TO 25
!  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESP. EQUATIONS AND CAR-
!  RY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE M-TH
!  EQUATION READS    C(4,M)*S(M) + C(3,M)*S(M+1) = C(2,M).
   19 DO 20 M=2,L
         G = -C(3,M+1)/C(4,M-1)
         C(2,M) = G*C(2,M-1) + 3.D0*(C(3,M)*C(4,M+1)+C(3,M+1)*C(4,M))
   20    C(4,M) = G*C(3,M-1) + 2.D0*(C(3,M) + C(3,M+1))
!     CONSTRUCT LAST EQUATION FROM THE SECOND BOUNDARY CONDITION, OF THE FORM
!           (-G*C(4,N-1))*S(N-1) + C(4,N)*S(N) = C(2,N)
!     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!     SUBSTITUTION, SINCE C ARRAY HAPPENS TO BE SET UP JUST RIGHT FOR IT
!     AT THIS POINT.
      IF (IBCEND-1)                     21,30,24
   21 IF (N .EQ. 3 .AND. IBCEND .EQ. 0) GO TO 22
!     NOT-A-KNOT AND N .GE. 3, AND EITHER N .GT. 3 OR ALSO NOT-A-KNOT AT
!     LEFT END POINT.
      G = C(3,N-1) + C(3,N)
      C(2,N) = ((C(3,N)+2.D0*G)*C(4,N)*C(3,N-1) + C(3,N)**2*(C(1,N-1)-C(1,N-2))/C(3,N-1))/G
      G = -G/C(4,N-1)
      C(4,N) = C(3,N-1)
                                        GO TO 29
!     EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND NOT-A-
!     KNOT AT LEFT END POINT).
   22 C(2,N) = 2.D0*C(4,N)
      C(4,N) = 1.D0
                                        GO TO 28
!     SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
   24 C(2,N) = 3.D0*C(4,N) + C(3,N)/2.D0*C(2,N)
      C(4,N) = 2.D0
                                        GO TO 28
   25 IF (IBCEND-1)                     26,30,24
   26 IF (IBCBEG .GT. 0)                GO TO 22
!     NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
      C(2,N) = C(4,N)
                                        GO TO 30
   28 G = -1.D0/C(4,N-1)
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
   29 C(4,N) = G*C(3,N-1) + C(4,N)
      C(2,N) = (G*C(2,N-1) + C(2,N))/C(4,N)
!  CARRY OUT BACK SUBSTITUTION
   30 DO 40 J=L,1,-1
   40    C(2,J) = (C(2,J) - C(3,J)*C(2,J+1))/C(4,J)
!  GENERATE CUBIC COEFFICIENTS IN EACH INTERVAL, I.E., THE DERIV.S
!  AT ITS LEFT ENDPOINT, FROM VALUE AND SLOPE AT ITS ENDPOINTS.
      DO 50 I=2,N
         DTAU = C(3,I)
         DIVDF1 = (C(1,I) - C(1,I-1))/DTAU
         DIVDF3 = C(2,I-1) + C(2,I) - 2.D0*DIVDF1
         C(3,I-1) = 2.D0*(DIVDF1 - C(2,I-1) - DIVDF3)/DTAU
   50    C(4,I-1) = (DIVDF3/DTAU)*(6.D0/DTAU)
         RETURN
   END subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!

!--------------------------------------------------------------------!
!  Extracted from "A practical guide to splines," 1st edition,       !
!  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.   !
!--------------------------------------------------------------------!

   SUBROUTINE INTERV ( XT, LXT, X, LEFT, MFLAG )
!  COMPUTES LEFT = MAX( I , 1 .LE. I .LE. LXT .AND. XT(I) .LE. X ).
!
!  XT.....A REAL SEQUENCE, OF LENGHT LXT, ASSUMED TO BE NONDECREASING
!  LXT.....NUMBER OF TERMS IN THE SEQUENCE XT.
!  X.....THE POINT WHOSE LOCATION WITH RESPECT TO THE SEQUENCE XT IS 
!        TO BE DETERMINED.
!  LEFT, MFLAG.....BOTH INTEGERS, WHOSE VALUE IS
!
!   1     -1      IF               X .LT. XT(1)
!   I      0      IF   XT(I)  .LE. X .LT. XT(I+1)
!  LXT     1      IF  XT(LXT) .LE. X
!         IN PARTICULAR, MFLAG = 0 IS THE 'USUAL' CASE. MFLAG .NE. 0
!         INDICATES THAT X LIES OUTSIDE THE HALFOPEN INTERVAL
!         XT(1) .LE. Y .LT. XT(LXT). THE ASYMETRIC TREATMENT OF THE 
!         INTERVAL IS DUE TO THE DECISION TO MAKE ALL PP FUNCTIONS CONT-
!         INUOUS FROM THE RIGHT.
!  THE PROGRAM IS DESIGNED TO BE EFFICIENT IN THE COMMON SITUATION THAT
!  IT IS CALLED REPEATEDLY, WITH X TAKEN FROM AN INCREASING OR DECREA-
!  SING SEQUENCE. THIS WILL HAPPEN, E.G., WHEN A PP FUNCTION IS TO BE
!  GRAPHED. THE FIRST GUESS FOR LEFT IS THEREFORE TAKEN TO BE THE VAL-
!  UE RETURNED AT THE PREVIOUS CALL AND STORED IN THE L O C A L VARIA-
!  BLE  ILO . A FIRST CHECK ASCERTAINS THAT  ILO .LT. LXT (THIS IS NEC-
!  ESSARY SINCE THE PRESENT CALL MAY HAVE NOTHING TO DO WITH THE PREVI- 
!  OUS CALL). THEN, IF  XT(ILO) .LE. X .LT. XT(ILO+1), WE SET  LEFT = 
!  ILO  AND ARE DONE AFTER JUST THREE COMPARISONS.
!     OTHERWISE, WE REPEATEDLY DOUBLE THE DIFFERENCE  ISTEP = IHI - ILO
!  WHILE ALSO MOVING  ILO  AND  IHI  IN THE DIRECTION OF  X, UNTIL
!                      XT(ILO) .LE. X .LT. XT(IHI) ,
!  AFTER WHICH WE USE BISECTION TO GET, IN ADDITION, ILO+1 = IHI .
!  LEFT = ILO  IS THE RETURNED.

      IMPLICIT NONE
      INTEGER  LEFT, LXT,MFLAG, IHI, ILO, ISTEP, MIDDLE
      REAL(8) X, XT(LXT)
      DATA ILO /1/
!     SAVE ILO  (A VALID FORTRAN STATEMENT IN THE NEW 1977 STANDARD)
      IHI = ILO + 1
      IF (IHI .LT. LXT)                 GO TO 20
         IF (X .GE. XT(LXT))            GO TO 110
         IF (LXT .LE. 1)                GO TO 90
         ILO = LXT - 1
         IHI = LXT
  20  IF (X .GE. XT(IHI))               GO TO 40
      IF (X .GE. XT(ILO))               GO TO 100
!     NOW X .LT. XT(ILO). DECREASE  ILO  TO CAPTURE X .
      ISTEP = 1
  31     IHI = ILO
         ILO = IHI - ISTEP
         IF (ILO .LE. 1)                GO TO 35
         IF (X .GE. XT(ILO))            GO TO 50
         ISTEP = ISTEP*2
                                        GO TO 31
  35  ILO = 1
      IF (X .LT. XT(1))                 GO TO 90
                                        GO TO 50
!              **** NOW X .GE. XT(IHI). INCREASE  IHI  TO CAPTURE  X . 
  40  ISTEP = 1
  41     ILO = IHI
         IHI = ILO + ISTEP
         IF (IHI .GE. LXT)              GO TO 45
         IF (X .LT. XT(IHI))            GO TO 50
         ISTEP = ISTEP*2
                                        GO TO 41
  45  IF (X .GE. XT(LXT))               GO TO 110
      IHI = LXT
!           **** NOW  XT(ILO) .LE. X .LT. XT(IHI). NARROW THE INTERVAL.
  50  MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO)              GO TO 100
!     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1 .
      IF (X .LT. XT(MIDDLE))            GO TO 53
         ILO = MIDDLE
                                        GO TO 50
  53     IHI = MIDDLE
                                        GO TO 50
  90  MFLAG = -1
      LEFT = 1
                                        RETURN
 100  MFLAG = 0
      LEFT = ILO
                                        RETURN
 110  MFLAG = 1
      LEFT = LXT
                                        RETURN
    END SUBROUTINE

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

!---------------------------------------------------------------------!
!  Extracted from "A practical guide to splines," 1st edition,        !
!  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.    !
!---------------------------------------------------------------------!

   REAL(8) FUNCTION PPVALU (BREAK, COEF, L, K, X, JDERIV)
   IMPLICIT NONE
!  CALCULATES VALUE AT X OF JDERIV-TH DERIVATIVE OF PP FCT FROM PP-REPR
!  TO BE EVALUATED. SPECIFICALLY, THE J-TH DERIVATIVE OF F IS      
!  GIVEN BY
!  (D**J)F(X) = COEF(J+1,I) + H*(COEF(J+2,I) + H*( ... (COEF(K-1,I) +
!                             + H*COEF(K,I)/K-J-I))/(K-J-2) ... )/2)/1
!  WITH  H = X - BREAK(I),  AND
!  I = MAX( 1 , MAX( J , BREAK(J) .LE. X , 1 .LE. J .LE. L ) ).
!  X.....THE POINT AT WHICH TO EVALUATE.
!  JDERIV.....INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUAT-
!  ED. ASSUMED TO BE ZERO OR POSITIVE.
!  PPVALU.....THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF F AT X.
!  THE INTERVAL INDEX I, APPROPRIATE FOR X, IS FOUND THROUGHT A  
!  CALL TO INTERV. THE FORMULA ABOVE FOR THE JDERIV-YH DERIVATIVE
!  OF F IS THEN EVALUATED (BY NESTED MULTIPLICATION).
      INTEGER JDERIV, K, L, I, M, NDUMMY
      REAL(8) BREAK(L), COEF(K,L), X, FMMJDR, H
      PPVALU = 0.D0
      FMMJDR = K - JDERIV
!     DERIVATIVES OF ORDER K OR HIGHER ARE IDENTICALLY ZERO.
      IF (FMMJDR .LE. 0.D0)           GO TO 99
!
!              FIND INDEX I OF LARGEST BREAKPOINT TO THE LEFT OF X.
!     CALL INTERV ( BREAK, L, X, I, NDUMMY )
!     EVALUATE JDERIV-TH DERIVATIVE OF I-TH POLYNOMIAL PIECE AT X.
      H = X - BREAK(I)
      DO 10 M=K,JDERIV+1,-1
         PPVALU = (PPVALU/FMMJDR)*H + COEF(M,I)
  10     FMMJDR = FMMJDR - 1.D0
  99  RETURN
      END function

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 !-------------------------------------------------------------------------------------!
 ! interpolate a G^qmc(-Lfak:Lfak) to G^smooth(-L:L) using the cubic-spline algorithm. !
 !-------------------------------------------------------------------------------------!

      subroutine interp(gtmp,g,Lfak,Lfak1,L)
      implicit real*8(a-h,o-z)
      integer L,Lfak,Lfak1
      double precision gtmp(-Lfak:Lfak),g(-L:L)
      double precision xa(Lfak1),ya(4,Lfak1)
      do 10 i=1,Lfak1
         xa(i)=float(i-1)/float(Lfak)
         ya(1,i) = gtmp(i-1)
 10   continue
      call CUBSPL(xa,ya,Lfak1,0,0)
      do 20 i=1,L
         x=float(i)/float(L)
         g(i) = PPVALU(xa,ya,Lfak,4,x,0)
 20   continue
      g(0)=gtmp(0)
      do 40 i=1,L
         g(-i)=-g(L-i)
 40   continue
      return
      end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine resampleit__b(x,y,xx,yy,ss)
 implicit none
 real(8)      :: x(:),y(:)
 real(8)      :: xx(:),yy(:),ss
 integer      :: i,j,k
 type(spline) :: rspline

  if(size(x)==size(xx))then
  if(maxval(abs(x-xx))<1.d-12) then
   yy=y
   return
  endif
  endif

 if(maxval(abs(y))<1.d-15) then
   yy=0.d0
   return
 endif

  call init_spline(rspline,x,y,4,ss)
  call evaluate_splineb(rspline,xx,yy)
  call kill_spline(rspline)

 return
 end subroutine

  !--------------------------------!
  !--------------------------------!

 subroutine resampleit___b(x,y,xx,yy,ss)
 implicit none
 complex(8)   :: y(:),yy(:)
 real(8)      :: x(:),xx(:)
 real(8)      :: ss,yyy(size(yy))
 integer      :: i,j,k
 type(spline) :: rspline

  if(size(x)==size(xx))then
  if(maxval(abs(x-xx))<1.d-12) then
   yy=y
   return
  endif
  endif
 if(maxval(abs(y))<1.d-15) then
   yy=0.d0
   return
 endif

  call init_spline(rspline,x,real(y),4,ss)
  call evaluate_splineb(rspline,xx,yyy)
  yy=yyy
  call kill_spline(rspline)
  call init_spline(rspline,x,aimag(y),4,ss)
  call evaluate_splineb(rspline,xx,yyy)
  yy=yy+imi*yyy
  call kill_spline(rspline)

 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine init_spline(rspline,omi,Frc,kkk2,s2,errors)
  implicit none
    !--------------------------------------------------------------!
    ! omi : mesh x, Frc : array to fit , kk : order of the splines
    ! nn  : number of knots
    ! tt  : knots
    ! cc  : coefficient of spline
    ! w   : weights
    !--------------------------------------------------------------!
  real(8)               :: omi(:),Frc(size(omi))
  type(spline)          :: rspline 
  integer               :: kk,kkk,i,j,kkk2
  real(8),optional      :: s2
  real(8),optional      :: errors(size(omi))
  real(4)               :: s
  integer               :: ier,nn,lwrk,k,tryit
  real(4)               :: fp,xl,xr

  rspline%ierr=0
  kk=size(omi)

  if(rspline%initialized>0) then
     if(testing) write(*,*) 'danger :  spline already initialized'
     goto 21
  endif
  rspline%initialized=1

  if(allocated(rspline%om))     deallocate(rspline%om)
  if(allocated(rspline%Frc))    deallocate(rspline%Frc)
  if(allocated(rspline%w ))     deallocate(rspline%w)
  if(allocated(rspline%cc ))    deallocate(rspline%cc)
  if(allocated(rspline%tt ))    deallocate(rspline%tt)
  if(allocated(rspline%der))    deallocate(rspline%der)
  if(allocated(rspline%smooth)) deallocate(rspline%smooth)

  kkk           = kkk2
  rspline%N     = kk
  rspline%order = kkk
  rspline%nest  = max(kk+kkk+1,2*kkk+2)
  if(kkk>5) then
    if(testing) write(*,*) 'error spline of too high order'
    rspline%ierr=1
    return
  endif
  if(kk<=kkk)    kkk=kk-1

  allocate(rspline%om(kk))
  allocate(rspline%Frc(kk))
  allocate(rspline%cc(rspline%nest))
  allocate(rspline%tt(rspline%nest))
  allocate(rspline%der(kkk,kk))
  allocate(rspline%smooth(kk))
  allocate(rspline%w(kk))


  21 continue

  rspline%om=0
  rspline%Frc=0
  rspline%cc=0
  rspline%tt=0
  rspline%der=0
  rspline%smooth=0
  rspline%w=0


  rspline%om=omi
  rspline%Frc=Frc
  rspline%w=(/( 1.d0,i=1,kk )/)

  if(present(errors))then
   rspline%w=abs(errors)
   rspline%w=(10.d0*rspline%w/maxval(rspline%w))**2.d0
   where(rspline%w<1.d-2) rspline%w=1.d-2
   rspline%w=1.d0/rspline%w
  endif

  tryit=0
  20 continue

  if(present(s2)) then
    call define_spline(s2=s2)
  else
    call define_spline
  endif

  if(ier==10) then
   if(tryit==0)then
    if(testing) write(*,*) 'error init spline, try to sort it out....'
    call sortitout
    tryit=1
    goto 20
   else
    rspline%ierr=1
    write(*,*) 'error init spline, nope, nothing possible....'
    if(strongstop) stop
   endif
  endif


 return
  
 contains
 
   !-----------------------!
   !-----------------------!
   !-----------------------!

 subroutine define_spline(s2) 
  implicit none
  real(8),optional  :: s2
  real(4)           :: wrk(rspline%N*(rspline%order+3)+rspline%nest*(7+4*rspline%order)+10),tt(rspline%nest),cc(rspline%nest)
  integer           :: iwrk(rspline%nest)

   rspline%ierr=0
   lwrk=size(wrk)
   if(present(s2))then
    s=real(s2)
   else
    s=0.
   endif
   if(rspline%order>=rspline%N) then
    write(*,*) 'spline not enough points'
    write(*,*) 'number of points : ' , rspline%N
    write(*,*) ' order of splines: ', rspline%order
    rspline%ierr=1
    return
   endif
   if(rspline%nest<=2*rspline%order+2) then
     write(*,*) 'spline nest too small'
     rspline%ierr=1
     return
   endif
   xl=rspline%om(1)-1.e-7
   xr=rspline%om(rspline%N)+1.e-7
   call curfit(0,rspline%N,rspline%om,rspline%Frc,rspline%w,xl,xr, &
            & rspline%order,s,rspline%nest,rspline%nn,rspline%tt,rspline%cc,fp,wrk,lwrk,iwrk,ier)

 end subroutine
 
   !-----------------------!
   !-----------------------!
   !-----------------------!

 subroutine sortitout
 implicit none
 integer :: order(rspline%N)
 real(8) :: iiFrc(rspline%N),iiom(rspline%N)
  iiFrc=rspline%Frc
  iiom=rspline%om
  call qsort_array(iiom,order)
  call qsort_adj_array(iiFrc,order)
  k=rspline%N
  call group_data_rrr(k,iiom,iiFrc,iiom,k)
  call kill_spline(rspline)
  rspline%ierr=0
  kk=k
  if(allocated(rspline%om))     deallocate(rspline%om)
  if(allocated(rspline%Frc))    deallocate(rspline%Frc)
  if(allocated(rspline%w ))     deallocate(rspline%w)
  if(allocated(rspline%cc ))    deallocate(rspline%cc)
  if(allocated(rspline%tt ))    deallocate(rspline%tt)
  if(allocated(rspline%der))    deallocate(rspline%der)
  if(allocated(rspline%smooth)) deallocate(rspline%smooth)
  kkk           = kkk2
  rspline%N     = kk
  rspline%order = kkk
  rspline%nest  = max(kk+kkk+1,2*kkk+2)
  if(kk<=kkk)    kkk=kk-1
  allocate(rspline%om(kk))
  allocate(rspline%Frc(kk))
  allocate(rspline%cc(rspline%nest))
  allocate(rspline%tt(rspline%nest))
  allocate(rspline%der(kkk,kk))
  allocate(rspline%smooth(kk))
  allocate(rspline%w(kk))
  rspline%om=0
  rspline%Frc=0
  rspline%cc=0
  rspline%tt=0
  rspline%der=0
  rspline%smooth=0
  rspline%w=0
  rspline%om=iiom
  rspline%Frc=iiFrc
  rspline%w=(/( 1.d0,i=1,kk )/)
 end subroutine

   !-----------------------!
   !-----------------------!
   !-----------------------!

 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine kill_spline(rspline)
  implicit none
  type(spline)    :: rspline
  if(allocated(rspline%om))     deallocate(rspline%om)
  if(allocated(rspline%Frc))    deallocate(rspline%Frc)
  if(allocated(rspline%w ))     deallocate(rspline%w)
  if(allocated(rspline%cc ))    deallocate(rspline%cc)
  if(allocated(rspline%tt ))    deallocate(rspline%tt)
  if(allocated(rspline%der))    deallocate(rspline%der)
  if(allocated(rspline%smooth)) deallocate(rspline%smooth)
  rspline%initialized=0
end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine evaluate_splineb(rspline,xx,yy)
 implicit none
 type(spline)      :: rspline
 integer           :: kk,nn,iii,mm,siz
 real(8)           :: xx(:),yy(:)
 real(4)           :: yy2(size(yy)),xxx(size(xx))
 integer           :: i,ier

   if(rspline%ierr>0) return
   mm=size(xx(:))
   kk=rspline%order
   xxx(1:size(xx))=xx
   call splev(rspline%tt,rspline%nn,rspline%cc,rspline%order,xxx,yy2,mm,ier)
   yy=yy2(1:size(yy))

 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine evaluate_splinebr(rspline,xx,yy)
 implicit none
 type(spline)      :: rspline
 integer           :: kk,nn,iii,mm,siz
 real(4)           :: xx(:),yy(:)
 integer           :: i,ier

   if(rspline%ierr>0) return
   mm=size(xx(:))
   kk=rspline%order
   call splev(rspline%tt,rspline%nn,rspline%cc,rspline%order,xx,yy,mm,ier)

 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

end module
