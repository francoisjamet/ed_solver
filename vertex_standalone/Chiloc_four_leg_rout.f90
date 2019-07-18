module chitilde
implicit none

logical,parameter:: bypass=.false.

contains

      complex(8) function  PhiM_ii(beta,Ei,Ej,Ek,El,w1,w2,w3,PHI_EPS)
      implicit none
      real(8)    :: beta,Ei,Ej,Ek,El,PHI_EPS
      complex(8) :: w1,w2,w3,I_w3_Ekl,w23_Ejl,I_w1_Eij,T24o6,w12_Eik,T11o12
        I_w3_Ekl    = 1.d0/(w3+(Ek-El))
        w23_Ejl     =   w2+w3+(Ej-El)
        I_w1_Eij    = 1.d0/(w1+(Ei-Ej))
        if(abs(w23_Ejl) > PHI_EPS)then
            T24o6 = 1.d0/w23_Ejl *  ( I_w1_Eij - 1.d0/(w1+w2+w3+Ei-El) )
        else
            T24o6 = I_w1_Eij * I_w1_Eij
        endif
        w12_Eik   = w1+w2+(Ei-Ek)
        if(abs(w12_Eik) > PHI_EPS)then
            T11o12 = -1.d0/w12_Eik
        else
            T11o12 = beta
        endif
        PhiM_ii= I_w3_Ekl *  (  T24o6 - 1.d0/(w2+Ej-Ek) * ( I_w1_Eij  +  T11o12 )  )
      end function

      complex(8) function PhiM_ji(beta,Ei,Ej,Ek,El,w1,w2,w3,PHI_EPS)
      implicit none
      real(8) :: beta,Ei,Ej,Ek,El,PHI_EPS
      complex(8) :: w1,w2,w3,I_w3_Ekl,w23_Eil,I_w1_Eji,T1o57
        I_w3_Ekl    = 1.d0/(w3+(Ek-El));
        w23_Eil     =   (w2+w3+Ei-El);
        I_w1_Eji    = 1.d0/(w1+Ej-Ei);
        if(abs(w23_Eil)>PHI_EPS)then
            T1o57 = 1.d0/w23_Eil * I_w1_Eji
        else
            T1o57 = I_w1_Eji * I_w1_Eji - beta*I_w1_Eji
        endif
        PhiM_ji= I_w3_Ekl * ( T1o57 - 1.d0/(w2+Ei-Ek) * I_w1_Eji )
      end function

      complex(8) function PhiM_ki(beta,Ei,Ej,Ek,El,w1,w2,w3,PHI_EPS)
      implicit none
      real(8) :: beta,Ei,Ej,Ek,El,PHI_EPS
      complex(8) :: w1,w2,w3,w12_Eki
        w12_Eki   = (w1+w2+Ek-Ei)
        if(abs(w12_Eki)>PHI_EPS)then
            PhiM_ki= -1.d0/(w3+Ei-El) * 1.d0/w12_Eki * 1.d0/(w2+Ej-Ei)
        else
            PhiM_ki=0.0
        endif
      end function

      complex(8) function PhiM_li(beta,Ei,Ej,Ek,El,w1,w2,w3,PHI_EPS)
      implicit none
      real(8) :: beta,Ei,Ej,Ek,El,PHI_EPS
      complex(8) :: w1,w2,w3,w23_Eji
        w23_Eji   = w2 + w3 + Ej - Ei
        if(abs(w23_Eji)>PHI_EPS) then
            PhiM_li = - 1.d0 / ( (w3+Ek-Ei) * w23_Eji * (w1+w2+w3+El-Ei) )
        else
            PhiM_li=    0.0
        endif
      end function

subroutine  chi_tilde_loc( &
 &    k_,k__,l_,l__,cutoff,op, w1,  w2,  w3, PHI_EPS, beta, Z, gsE, sites, nup, ndn, &
 &    cp_i_E,            dim_E_i,          cp_pup_E,          dim_E_pup,     &
 &    cp_pdn_E,          dim_E_pdn,        cp_mup_E,          dim_E_mup,     &
 &    cp_mdn_E,          dim_E_mdn,        cp_p2dn_E,         dim_E_p2dn,    &
 &    cp_m2dn_E,         dim_E_m2dn,       cp_puppdn_E,       dim_E_puppdn,  &
 &    cp_muppdn_E,       dim_E_muppdn,     cp_pupmdn_E,       dim_E_pupmdn,  &
 &    cp_mupmdn_E,       dim_E_mupmdn,     cp_i_cdup,         cp_i_cddn,     &
 &    cp_pup_cddn,       cp_pdn_cdup,      cp_pdn_cddn,       cp_mup_cdup,   &
 &    cp_mup_cddn,       cp_mdn_cddn,      cp_mdn_cdup,       cp_muppdn_cdup,&
 &    cp_pupmdn_cddn,    cp_mupmdn_cdup,   cp_mupmdn_cddn,    cp_m2dn_cddn,  pDNDN,pUPDN )
implicit none
integer    :: op,sites,nup,ndn,k_,l_,k__,l__,k1,l1,k2,l2
complex(8) :: w1,w2,w3
real(8)    :: PHI_EPS,beta,Z,gsE
complex(8) :: pDNDN,pUPDN
integer    :: dim_E_i
integer    :: dim_E_pup
integer    :: dim_E_pdn
integer    :: dim_E_mup
integer    :: dim_E_mdn
integer    :: dim_E_p2dn
integer    :: dim_E_m2dn
integer    :: dim_E_puppdn
integer    :: dim_E_muppdn
integer    :: dim_E_pupmdn
integer    :: dim_E_mupmdn
real(8)    ::           cp_i_E(:)
real(8)    ::         cp_pup_E(:)
real(8)    ::         cp_pdn_E(:)
real(8)    ::         cp_mup_E(:)
real(8)    ::         cp_mdn_E(:)
real(8)    ::        cp_p2dn_E(:)
real(8)    ::        cp_m2dn_E(:)
real(8)    ::      cp_puppdn_E(:)
real(8)    ::      cp_muppdn_E(:)
real(8)    ::      cp_pupmdn_E(:)
real(8)    ::      cp_mupmdn_E(:)
real(8)    ::      cp_i_cdup(:,:,:)
real(8)    ::      cp_i_cddn(:,:,:)
real(8)    ::    cp_pup_cddn(:,:,:)
real(8)    ::    cp_pdn_cdup(:,:,:)
real(8)    ::    cp_pdn_cddn(:,:,:)
real(8)    ::    cp_mup_cdup(:,:,:)
real(8)    ::    cp_mup_cddn(:,:,:)
real(8)    ::    cp_mdn_cddn(:,:,:)
real(8)    ::    cp_mdn_cdup(:,:,:)
real(8)    :: cp_muppdn_cdup(:,:,:)
real(8)    :: cp_pupmdn_cddn(:,:,:)
real(8)    :: cp_mupmdn_cdup(:,:,:)
real(8)    :: cp_mupmdn_cddn(:,:,:)
real(8)    ::   cp_m2dn_cddn(:,:,:)
real(8)    :: boltzZ
complex(8) :: cDNDN,cUPDN,xi1,yi1,xi2,yi2
integer    :: dimi,diml,dimk,dimj,i,j,k,l,stati,stati_,u_,d_
integer    :: mu,pu,md,pd
real(8)    :: cutoff
real(8)    :: d1,d2
real(8)    :: cc_cutoff,cc

!various notation for the same thing
u_=k_;
d_=l_;
k1=k_
l1=l_
k2=k__
l2=l__
mu=k1
pu=k2
md=l1
pd=l2
cc_cutoff = 1d-4
pDNDN = CMPLX(0.d0,0.d0);pUPDN = CMPLX(0.d0,0.d0)

! call run_tests(k_) ; if(k_/=k_) call run_tests(l_) ;  return

   do stati=1,dim_E_i

    boltzZ = exp(-(cp_i_E(stati)-gsE)*beta)/Z;

    if(boltzZ<cutoff) cycle

    dimi = dim_E_i

    if(op == 1 .and. (k2==k1.and.l2==l1.and.k1==l2))then

       cDNDN = CMPLX(0.d0,0.d0)

        if(bypass.or.ndn/=sites) then
            diml = dim_E_pdn ;
            dimk = dim_E_i;
            dimj = dim_E_pdn;
            do j =1, dimj;
               do l =1, diml;
               cc =  cp_i_cddn(k_,l,+stati) *  cp_i_cddn(k_,j,+stati)
               if( abs(cc) < cc_cutoff ) cycle
               do k =1, dimk;
                xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w1,w2,w3,PHI_EPS);
                xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w3,w2,w1,PHI_EPS);
                cDNDN =cDNDN+  (xi1+xi2)* cp_i_cddn(k_,j,+stati)* cp_i_cddn(k_,j,+k) * cp_i_cddn(k_,l,+k) * cp_i_cddn(k_,l,+stati);
            enddo;enddo;enddo
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_i_cddn(k_,j,+stati) *  cp_i_cddn(k_,l,+stati)
                  if( abs(cc) < cc_cutoff ) cycle
                  do k =1, dimk;
                yi1 = PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w1,w2,w3,PHI_EPS);
                yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w3,w2,w1,PHI_EPS);

                cDNDN =cDNDN+  (yi1+yi2)* cp_i_cddn(k_,l,+stati)* cp_i_cddn(k_,l,+k) * cp_i_cddn(k_,j,+k) * cp_i_cddn(k_,j,+stati);
            enddo;enddo;enddo
        endif


        if(bypass.or.(ndn/=sites .and. ndn/=0))then
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_i_cddn(k_,l,+stati) *  cp_mdn_cddn(k_,stati,+j)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
                     xi1 =-PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w1,w3,PHI_EPS);
                     xi2 = PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w3,w1,PHI_EPS);
                     cDNDN = cDNDN +  (xi1+xi2)* cp_mdn_cddn(k_,stati,+j)* cp_mdn_cddn(k_,k,+j) * cp_i_cddn(k_,l,+k) * cp_i_cddn(k_,l,+stati);
            enddo;enddo;enddo

            diml = dim_E_mdn;
            dimk = dim_E_pdn;
            dimj = dim_E_i;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_i_cddn(k_,k,+stati) *  cp_mdn_cddn(k_,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
                     yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),w3,w1,w2,PHI_EPS);
                     yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),w1,w3,w2,PHI_EPS);
                     cDNDN = cDNDN+  (yi1+yi2)* cp_mdn_cddn(k_,stati,+l)* cp_mdn_cddn(k_,j,+l) * cp_i_cddn(k_,k,+j) * cp_i_cddn(k_,k,+stati);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;
            do j =1, dimj;
               do l =1, diml;

                  cc =  cp_mdn_cddn(k_,stati,+j) *  cp_i_cddn(k_,l,+stati)
                  if (abs(cc) < cc_cutoff) cycle
               do k =1, dimk;

                xi1 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w1,w3,PHI_EPS);
                xi2 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w3,w1,PHI_EPS);
                cDNDN =cDNDN+  (xi1+xi2)* cp_i_cddn(k_,l,+stati)* cp_i_cddn(k_,l,+k) * cp_mdn_cddn(k_,k,+j) * cp_mdn_cddn(k_,stati,+j);
            enddo;enddo;enddo

            diml = dim_E_i;
            dimk = dim_E_pdn;
            dimj = dim_E_mdn;
            do j =1, dimj;
               do k =1, dimk;
                  cc =  cp_mdn_cddn(k_,stati,+j) *  cp_i_cddn(k_,k,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do l =1, diml;
                yi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),w3,w1,w2,PHI_EPS);
                yi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),w1,w3,w2,PHI_EPS);
                cDNDN =cDNDN+  (yi1+yi2)* cp_i_cddn(k_,k,+stati)* cp_i_cddn(k_,k,+l) * cp_mdn_cddn(k_,l,+j) * cp_mdn_cddn(k_,stati,+j);
            enddo;enddo;enddo
        endif

        if(bypass.or.ndn<sites-1)then

            diml = dim_E_p2dn;
            dimk = dim_E_pdn;
            dimj = dim_E_pdn;
            do j =1, dimj;
               do k =1, dimk;
                  cc =  cp_i_cddn(k_,j,+stati) *  cp_i_cddn(k_,k,+stati)
                  if (abs(cc) < cc_cutoff) cycle
            do l =1, diml;
                yi1 =-PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),w2,w1,w3,PHI_EPS);
                yi2 = PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),w2,w3,w1,PHI_EPS);
                cDNDN =cDNDN+  (yi1+yi2)* cp_i_cddn(k_,k,+stati)* cp_pdn_cddn(k_,l,+k) * cp_pdn_cddn(k_,l,+j) * cp_i_cddn(k_,j,+stati);
            enddo;enddo;enddo
            diml = dim_E_pdn;
            dimk = dim_E_p2dn;
            dimj = dim_E_pdn;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_i_cddn(k_,l,+stati) *  cp_i_cddn(k_,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;

                xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),w3,w1,w2,PHI_EPS);
                xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),w1,w3,w2,PHI_EPS);
                cDNDN =cDNDN+  (xi1+xi2)* cp_i_cddn(k_,j,+stati)* cp_pdn_cddn(k_,k,+j) * cp_pdn_cddn(k_,k,+l) * cp_i_cddn(k_,l,+stati);
            enddo;enddo;enddo
        endif

        if(bypass.or.ndn>1)then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_m2dn;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mdn_cddn(k_,stati,+k) *  cp_mdn_cddn(k_,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
                     yi1 =-PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),w2,w1,w3,PHI_EPS);
                     yi2 = PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),w2,w3,w1,PHI_EPS);
                     cDNDN =cDNDN+  (yi1+yi2)* cp_mdn_cddn(k_,stati,+l)* cp_m2dn_cddn(k_,l,+j) * cp_m2dn_cddn(k_,k,+j) * cp_mdn_cddn(k_,stati,+k);
                  enddo;enddo;enddo

            diml = dim_E_mdn;
            dimk = dim_E_m2dn;
            dimj = dim_E_mdn;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_mdn_cddn(k_,stati,+j) *  cp_mdn_cddn(k_,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
            do k =1, dimk;
                yi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),w3,w1,w2,PHI_EPS);
                yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),w1,w3,w2,PHI_EPS);
                cDNDN =cDNDN+  (yi1+yi2)* cp_mdn_cddn(k_,stati,+l)* cp_m2dn_cddn(k_,l,+k) * cp_m2dn_cddn(k_,j,+k) * cp_mdn_cddn(k_,stati,+j);
            enddo;enddo;enddo
        endif

        if(bypass.or.ndn/=0)then

            diml = dim_E_i;
            dimk = dim_E_mdn;
            dimj = dim_E_mdn;
            do j =1, dimj;
            do k =1, dimk;

                  cc =  cp_mdn_cddn(k_,stati,+j) *  cp_mdn_cddn(k_,stati,+k)
                  if (abs(cc) < cc_cutoff) cycle
                  do l =1, diml;
                xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),w1,w2,w3,PHI_EPS);
                xi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),w3,w2,w1,PHI_EPS);
                cDNDN =cDNDN+  (xi1+xi2)* cp_mdn_cddn(k_,stati,+k) * cp_mdn_cddn(k_,l,+k) * cp_mdn_cddn(k_,l,+j) * cp_mdn_cddn(k_,stati,+j);
            enddo;enddo;enddo

            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_i;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mdn_cddn(k_,stati,+k) *  cp_mdn_cddn(k_,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
                yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),w1,w2,w3,PHI_EPS);
                yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),w3,w2,w1,PHI_EPS);
                cDNDN =cDNDN+  (yi1+yi2)* cp_mdn_cddn(k_,stati,+l) * cp_mdn_cddn(k_,j,+l) * cp_mdn_cddn(k_,j,+k) * cp_mdn_cddn(k_,stati,+k);
           enddo;enddo;enddo
        endif

        pDNDN = pDNDN+cDNDN*boltzZ;

    endif

   if(op == 1 .and. (k2/=k1.or.l2/=l1.or.k1/=l2))then

       cDNDN = CMPLX(0.d0,0.d0)

! For UPDN *all* terms are important:
! tot  perm    w   w   w         matrix elements      fermi       //  Indices:
!                                k1  k2  l1  l2
!   0   0      1   2   3         -u  +u  -d  +d   ii    1         //    ii : i j k l i
!   1   0      1   2   3         +u  -d  +d  -u   ji    1         //    ji : i k l j i
!   2   0      1   2   3         -d  +d  -u  +u   ki    1         //    ki : i l k j i
!   3   0      1   2   3         +d  -u  +u  -d   li    1         //    li : i l j k i
!
!   4   1      2   1   3         +u  -u  -d  +d   ii   -1
!   5   1      2   1   3         -u  -d  +d  +u   ji   -1
!   6   1      2   1   3         -d  +d  +u  -u   ki   -1
!   7   1      2   1   3         +d  +u  -u  -d   li   -1
!
!   8   2      3   1   2         -d  -u  +u  +d   ii    1
!   9   2      3   1   2         -u  +u  +d  -d   ji    1
!  10   2      3   1   2         +u  +d  -d  -u   ki    1
!  11   2      3   1   2         +d  -d  -u  +u   li    1
!
!  12   3      1   3   2         -u  -d  +u  +d   ii   -1
!  13   3      1   3   2         -d  +u  +d  -u   ji   -1
!  14   3      1   3   2         +u  +d  -u  -d   ki   -1
!  15   3      1   3   2         +d  -u  -d  +u   li   -1
!
!  16   4      2   3   1         +u  -d  -u  +d   ii    1
!  17   4      2   3   1         -d  -u  +d  +u   ji    1
!  18   4      2   3   1         -u  +d  +u  -d   ki    1
!  19   4      2   3   1         +d  +u  -d  -u   li    1
!
!  20   5      3   2   1         -d  +u  -u  +d   ii   -1
!  21   5      3   2   1         +u  -u  +d  -d   ji   -1
!  22   5      3   2   1         -u  +d  -d  +u   ki   -1
!  23   5      3   2   1         +d  -d  +u  -u   li   -1

        if(bypass.or.ndn/=sites) then
            diml = dim_E_pdn ;
            dimk = dim_E_i;
            dimj = dim_E_pdn;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cddn(l2,l,+stati) *  cp_i_cddn(k1,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
            do k =1, dimk;
               xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w1,w2,w3,PHI_EPS);
                cDNDN =cDNDN+  (xi1)* cp_i_cddn(k1,j,+stati)* cp_i_cddn(k2,j,+k) * cp_i_cddn(l1,l,+k) * cp_i_cddn(l2,l,+stati);

                xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w3,w2,w1,PHI_EPS);
                cDNDN =cDNDN+  (xi2)* cp_i_cddn(l1,j,+stati)* cp_i_cddn(k2,j,+k) * cp_i_cddn(k1,l,+k) * cp_i_cddn(l2,l,+stati);

            enddo;enddo;enddo
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cddn(k1,j,+stati) *  cp_i_cddn(k1,l,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  cc =  cp_i_cddn(k2,j,+stati) *  cp_i_cddn(l1,l,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               yi1 = PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w1,w2,w3,PHI_EPS);
                cDNDN =cDNDN+  ( yi1)* cp_i_cddn(l1,l,+stati)* cp_i_cddn(l2,l,+k) * cp_i_cddn(k1,j,+k) * cp_i_cddn(k2,j,+stati);

                yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),w3,w2,w1,PHI_EPS);
                cDNDN =cDNDN+  (+yi2)* cp_i_cddn(k1,l,+stati)* cp_i_cddn(l2,l,+k) * cp_i_cddn(l1,j,+k) * cp_i_cddn(k1,j,+stati);

            enddo;enddo;enddo
        endif


        if(bypass.or.(ndn/=sites .and. ndn/=0))then
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_mdn_cddn(pu,stati,+j)
                  if (abs(cc) < cc_cutoff) cycle
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_mdn_cddn(pu,stati,+j)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
                     xi1 =-PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w1,w3,PHI_EPS);
                     cDNDN = cDNDN + (xi1)* cp_mdn_cddn(pu,stati,+j)* cp_mdn_cddn(mu,k,+j) * cp_i_cddn(md,l,+k) * cp_i_cddn(pd,l,+stati);

                     xi2 = PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w3,w1,PHI_EPS);
                     cDNDN = cDNDN + (xi2)* cp_mdn_cddn(pu,stati,+j)* cp_mdn_cddn(md,k,+j) * cp_i_cddn(mu,l,+k) * cp_i_cddn(pd,l,+stati);

                  enddo;enddo;enddo

            diml = dim_E_mdn;
            dimk = dim_E_pdn;
            dimj = dim_E_i;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_i_cddn(pu,k,+stati) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  cc =  cp_i_cddn(pu,k,+stati) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
                     yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),w3,w1,w2,PHI_EPS);
                     cDNDN = cDNDN+  (yi1)* cp_mdn_cddn(pd,stati,+l)* cp_mdn_cddn(md,j,+l) * cp_i_cddn(mu,k,+j) * cp_i_cddn(pu,k,+stati);
                     yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),w1,w3,w2,PHI_EPS);
                     cDNDN = cDNDN+  (+yi2)* cp_mdn_cddn(pd,stati,+l)* cp_mdn_cddn(mu,j,+l) * cp_i_cddn(md,k,+j) * cp_i_cddn(pu,k,+stati);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;
            do j =1, dimj;
            do k =1, dimk;
            do l =1, diml;
                xi1 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w1,w3,PHI_EPS);
                cDNDN =cDNDN+  (xi1)* cp_i_cddn(md,l,+stati)* cp_i_cddn(pd,l,+k) * cp_mdn_cddn(pu,k,+j) * cp_mdn_cddn(mu,stati,+j);

                xi2 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),w2,w3,w1,PHI_EPS);
                cDNDN =cDNDN+  (+xi2)* cp_i_cddn(mu,l,+stati)* cp_i_cddn(pd,l,+k) * cp_mdn_cddn(pu,k,+j) * cp_mdn_cddn(md,stati,+j);

            enddo;enddo;enddo

            diml = dim_E_i;
            dimk = dim_E_pdn;
            dimj = dim_E_mdn;
            do j =1, dimj;
            do k =1, dimk;
               cc =  cp_mdn_cddn(mu,stati,+j) *  cp_i_cddn(md,k,+stati)
               if (abs(cc) < cc_cutoff) cycle
               cc =  cp_mdn_cddn(md,stati,+j) *  cp_i_cddn(mu,k,+stati)
               if (abs(cc) < cc_cutoff) cycle
               do l =1, diml;
                yi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),w3,w1,w2,PHI_EPS);
                cDNDN =cDNDN+  (yi1)* cp_i_cddn(mu,k,+stati)* cp_i_cddn(pu,k,+l) * cp_mdn_cddn(pd,l,+j) * cp_mdn_cddn(md,stati,+j);
                yi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),w1,w3,w2,PHI_EPS);
                cDNDN =cDNDN+  (yi2)* cp_i_cddn(md,k,+stati)* cp_i_cddn(pu,k,+l) * cp_mdn_cddn(pd,l,+j) * cp_mdn_cddn(mu,stati,+j);
            enddo;enddo;enddo
        endif

        if(bypass.or.ndn<sites-1)then

            diml = dim_E_p2dn;
            dimk = dim_E_pdn;
            dimj = dim_E_pdn;
            do j =1, dimj;
            do k =1, dimk;

               cc =  cp_i_cddn(pu,j,+stati) *  cp_i_cddn(mu,k,+stati)
               if (abs(cc) < cc_cutoff) cycle
               cc =  cp_i_cddn(pu,j,+stati) *  cp_i_cddn(md,k,+stati)
               if (abs(cc) < cc_cutoff) cycle
               do l =1, diml;
               yi1 =-PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),w2,w1,w3,PHI_EPS);
               cDNDN =cDNDN+  ( yi1) * cp_i_cddn(mu,k,+stati)* cp_pdn_cddn(md,l,+k) * cp_pdn_cddn(pd,l,+j) * cp_i_cddn(pu,j,+stati);
               yi2 = PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),w2,w3,w1,PHI_EPS);
               cDNDN =cDNDN+  (+yi2) * cp_i_cddn(md,k,+stati)* cp_pdn_cddn(mu,l,+k) * cp_pdn_cddn(pd,l,+j) * cp_i_cddn(pu,j,+stati);
            enddo;enddo;enddo
            diml = dim_E_pdn;
            dimk = dim_E_p2dn;
            dimj = dim_E_pdn;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_i_cddn(md,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_i_cddn(mu,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
                     xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),w3,w1,w2,PHI_EPS);
                     cDNDN =cDNDN+  (xi1)* cp_i_cddn(md,j,+stati)* cp_pdn_cddn(mu,k,+j) * cp_pdn_cddn(pu,k,+l) * cp_i_cddn(pd,l,+stati);
                     xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),w1,w3,w2,PHI_EPS);
                     cDNDN =cDNDN+  (xi2)* cp_i_cddn(mu,j,+stati)* cp_pdn_cddn(md,k,+j) * cp_pdn_cddn(pu,k,+l) * cp_i_cddn(pd,l,+stati);
                  enddo;enddo;enddo
        endif

        if(bypass.or.ndn>1)then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_m2dn;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mdn_cddn(md,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  cc =  cp_mdn_cddn(mu,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
            do j =1, dimj;
               yi1 =-PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),w2,w1,w3,PHI_EPS);
                cDNDN =cDNDN+  (yi1)* cp_mdn_cddn(pd,stati,+l)* cp_m2dn_cddn(pu,l,+j) * cp_m2dn_cddn(mu,k,+j) * cp_mdn_cddn(md,stati,+k);
                yi2 = PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),w2,w3,w1,PHI_EPS);
                cDNDN =cDNDN+  (yi2)* cp_mdn_cddn(pd,stati,+l)* cp_m2dn_cddn(pu,l,+j) * cp_m2dn_cddn(md,k,+j) * cp_mdn_cddn(mu,stati,+k);
            enddo;enddo;enddo

            diml = dim_E_mdn;
            dimk = dim_E_m2dn;
            dimj = dim_E_mdn;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_mdn_cddn(mu,stati,+j) *  cp_mdn_cddn(pu,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  cc =  cp_mdn_cddn(md,stati,+j) *  cp_mdn_cddn(pu,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
                     yi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),w3,w1,w2,PHI_EPS);
                     cDNDN = cDNDN+  ( yi1)* cp_mdn_cddn(pu,stati,+l)* cp_m2dn_cddn(pd,l,+k) * cp_m2dn_cddn(md,j,+k) * cp_mdn_cddn(mu,stati,+j);
                     yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),w1,w3,w2,PHI_EPS);
                     cDNDN = cDNDN+  (+yi2)* cp_mdn_cddn(pu,stati,+l)* cp_m2dn_cddn(pd,l,+k) * cp_m2dn_cddn(mu,j,+k) * cp_mdn_cddn(md,stati,+j);
                  enddo;enddo;enddo
        endif

        if(bypass.or.ndn/=0)then

            diml = dim_E_i;
            dimk = dim_E_mdn;
            dimj = dim_E_mdn;
            do j =1, dimj;
               do k =1, dimk;

                  cc =  cp_mdn_cddn(mu,stati,+j) *  cp_mdn_cddn(pu,stati,+k)
                  if (abs(cc) < cc_cutoff) cycle
                  cc =  cp_mdn_cddn(md,stati,+j) *  cp_mdn_cddn(pu,stati,+k)
                  if (abs(cc) < cc_cutoff) cycle
                  do l =1, diml;
                     xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),w1,w2,w3,PHI_EPS);
                     cDNDN =cDNDN+  (xi1)* cp_mdn_cddn(pu,stati,+k) * cp_mdn_cddn(md,l,+k) * cp_mdn_cddn(pd,l,+j) * cp_mdn_cddn(mu,stati,+j);
                     xi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),w3,w2,w1,PHI_EPS);
                     cDNDN =cDNDN+  (xi2)* cp_mdn_cddn(pu,stati,+k) * cp_mdn_cddn(mu,l,+k) * cp_mdn_cddn(pd,l,+j) * cp_mdn_cddn(md,stati,+j);
                  enddo;enddo;enddo

            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_i;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mdn_cddn(mu,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                if (abs(cc) < cc_cutoff) cycle
                cc =  cp_mdn_cddn(md,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                if (abs(cc) < cc_cutoff) cycle
                do j =1, dimj;
                yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),w1,w2,w3,PHI_EPS);
                cDNDN =cDNDN+  ( yi1)* cp_mdn_cddn(pd,stati,+l) * cp_mdn_cddn(mu,j,+l) * cp_mdn_cddn(pu,j,+k) * cp_mdn_cddn(md,stati,+k);
                yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),w3,w2,w1,PHI_EPS);
                cDNDN =cDNDN+  (+yi2)* cp_mdn_cddn(pd,stati,+l) * cp_mdn_cddn(md,j,+l) * cp_mdn_cddn(pu,j,+k) * cp_mdn_cddn(mu,stati,+k);
           enddo;enddo;enddo
        endif
        pDNDN = pDNDN+cDNDN*boltzZ;

    endif


! **** UPDN *****************************************************************************

    if(op == 2)then
        cUPDN = CMPLX(0.d0,0.d0)
        if(bypass.or.(ndn/=sites .and. nup /= sites))then
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_pup;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_i_cdup(mu,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pup_E(j),cp_i_E(k),cp_pdn_E(l),w1,w2,w3,PHI_EPS);
              cUPDN =cUPDN + xi1 * cp_i_cdup(mu,j,+stati)* cp_i_cdup(pu,j,+k) * cp_i_cddn(md,l,+k) * cp_i_cddn(pd,l,+stati);
            enddo;enddo;enddo

            dimj = dim_E_pup;
            dimk = dim_E_i;
            diml = dim_E_pdn;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cdup(pu,j,+stati) *  cp_i_cddn(md,l,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1 = PhiM_ki(beta,cp_i_E(stati), cp_pup_E(j),cp_i_E(k),cp_pdn_E(l),w1,w2,w3,PHI_EPS);
             cUPDN =cUPDN + xi1 * cp_i_cddn(md,l,+stati)* cp_i_cddn(pd,l,+k) * cp_i_cdup(mu,j,+k) * cp_i_cdup(pu,j,+stati);
            enddo;enddo;enddo

            diml = dim_E_puppdn;
            dimk = dim_E_pup;
            dimj = dim_E_pup;
            do j =1, dimj;
            do k =1, dimk;
               cc =  cp_i_cdup(pu,j,+stati) *  cp_i_cdup(mu,k,+stati)
               if (abs(cc) < cc_cutoff) cycle
               do l =1, diml;
               xi1 = PhiM_ji(beta,cp_i_E(stati), cp_pup_E(j),cp_pup_E(k),cp_puppdn_E(l),w2,w1,w3,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_i_cdup(mu,k,+stati)* cp_pup_cddn(md,l,+k) * cp_pup_cddn(pd,l,+j) * cp_i_cdup(pu,j,+stati);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_puppdn;
            dimj = dim_E_pdn;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_i_cddn(md,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1   = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_puppdn_E(k),cp_pdn_E(l),w3,w1,w2,PHI_EPS);
                cUPDN = cUPDN + xi1 * cp_i_cddn(md,j,+stati)* cp_pdn_cdup(mu,k,+j) * cp_pdn_cdup(pu,k,+l) * cp_i_cddn(pd,l,+stati);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_puppdn;
            dimj = dim_E_pup;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_i_cdup(mu,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pup_E(j),cp_puppdn_E(k),cp_pdn_E(l),w1,w3,w2,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_i_cdup(mu,j,+stati)* cp_pup_cddn(md,k,+j) * cp_pdn_cdup(pu,k,+l) * cp_i_cddn(pd,l,+stati);
            enddo;enddo;enddo

            diml = dim_E_puppdn;
            dimk = dim_E_pdn;
            dimj = dim_E_pup;
            do j =1, dimj;
            do k =1, dimk;

                  cc =  cp_i_cdup(pu,j,+stati) *  cp_i_cddn(md,k,+stati)
                  if (abs(cc) < cc_cutoff) cycle
               do l =1, diml;
               xi1 = PhiM_ji(beta,cp_i_E(stati), cp_pup_E(j),cp_pdn_E(k),cp_puppdn_E(l),w2,w3,w1,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_i_cddn(md,k,+stati)* cp_pdn_cdup(mu,l,+k) * cp_pup_cddn(pd,l,+j) * cp_i_cdup(pu,j,+stati);
            enddo;enddo;enddo
        endif

        if(bypass.or.(ndn/=sites .and. nup /= 0)) then
            diml = dim_E_muppdn;
            dimk = dim_E_mup;
            dimj = dim_E_mup;
            do j =1, dimj;
            do k =1, dimk;

                  cc =  cp_mup_cdup(mu,stati,+j) *  cp_mup_cdup(pu,stati,+k)
                  if (abs(cc) < cc_cutoff) cycle
               do l =1, diml;
               xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mup_E(j),cp_mup_E(k),cp_muppdn_E(l),w1,w2,w3,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_mup_cdup(pu,stati,+k)* cp_mup_cddn(md,l,+k) * cp_mup_cddn(pd,l,+j) * cp_mup_cdup(mu,stati,+j);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mup;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_mup_cdup(pu,stati,+j)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1 = PhiM_ii(beta,cp_i_E(stati), cp_mup_E(j),cp_i_E(k),cp_pdn_E(l),w2,w1,w3,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_mup_cdup(pu,stati,+j)* cp_mup_cdup(mu,k,+j) * cp_i_cddn(md,l,+k) * cp_i_cddn(pd,l,+stati);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mup;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_mup_cdup(mu,stati,+j) *  cp_i_cddn(md,l,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mup_E(j),cp_i_E(k),cp_pdn_E(l),w2,w1,w3,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_i_cddn(md,l,+stati)* cp_i_cddn(pd,l,+k) * cp_mup_cdup(pu,k,+j) * cp_mup_cdup(mu,stati,+j);
            enddo;enddo;enddo
            diml = dim_E_muppdn;
            dimk = dim_E_pdn;
            dimj = dim_E_mup;
            do j =1, dimj;
            do k =1, dimk;

                  cc =  cp_mup_cdup(mu,stati,+j) *  cp_i_cddn(md,k,+stati)
                  if (abs(cc) < cc_cutoff) cycle
               do l =1, diml;
               xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mup_E(j),cp_pdn_E(k),cp_muppdn_E(l),w1,w3,w2,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_i_cddn(md,k,+stati)* cp_muppdn_cdup(pu,k,+l) * cp_mup_cddn(pd,l,+j) * cp_mup_cdup(mu,stati,+j);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_muppdn;
            dimj = dim_E_mup;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_mup_cdup(pu,stati,+j)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1 = PhiM_ii(beta,cp_i_E(stati), cp_mup_E(j),cp_muppdn_E(k),cp_pdn_E(l),w2,w3,w1,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_mup_cdup(pu,stati,+j)* cp_mup_cddn(md,k,+j) * cp_muppdn_cdup(mu,l,+k) * cp_i_cddn(pd,l,+stati);
            enddo;enddo;enddo

            diml = dim_E_pdn;
            dimk = dim_E_muppdn;
            dimj = dim_E_pdn;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_i_cddn(pd,l,+stati) *  cp_i_cddn(md,j,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
                     xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_muppdn_E(k),cp_pdn_E(l),w3,w2,w1,PHI_EPS);
                     cUPDN =cUPDN - xi1 * cp_i_cddn(md,j,+stati)* cp_muppdn_cdup(pu,j,+k) * cp_muppdn_cdup(mu,l,+k) * cp_i_cddn(pd,l,+stati);
            enddo;enddo;enddo
        endif

        if(bypass.or.(ndn/=0 .and. nup /= sites))then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_pupmdn;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mdn_cddn(md,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
               xi1 = PhiM_li(beta,cp_i_E(stati), cp_pupmdn_E(j),cp_mdn_E(k),cp_mdn_E(l),w1,w2,w3,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_mdn_cddn(pd,stati,+l)* cp_mdn_cdup(mu,j,+l) * cp_mdn_cdup(pu,j,+k) * cp_mdn_cddn(md,stati,+k);
            enddo;enddo;enddo
            diml = dim_E_i;
            dimk = dim_E_pup;
            dimj = dim_E_mdn;
            do j =1, dimj;
            do k =1, dimk;

                  cc =  cp_mdn_cddn(md,stati,+j) *  cp_i_cdup(mu,k,+stati)
                  if (abs(cc) < cc_cutoff) cycle
               do l =1, diml;
               xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pup_E(k),cp_i_E(l),w3,w1,w2,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_i_cdup(mu,k,+stati)* cp_i_cdup(pu,k,+l) * cp_mdn_cddn(pd,l,+j) * cp_mdn_cddn(md,stati,+j);
            enddo;enddo;enddo
            diml = dim_E_mdn;
            dimk = dim_E_pup;
            dimj = dim_E_i;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_i_cdup(pu,k,+stati) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
               xi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pup_E(k),cp_mdn_E(l),w3,w1,w2,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_mdn_cddn(pd,stati,+l)* cp_mdn_cddn(md,j,+l) * cp_i_cdup(mu,k,+j) * cp_i_cdup(pu,k,+stati);
            enddo;enddo;enddo
            diml = dim_E_mdn;
            dimk = dim_E_pup;
            dimj = dim_E_pupmdn;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_i_cdup(pu,k,+stati) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
               xi1 = PhiM_li(beta,cp_i_E(stati), cp_pupmdn_E(j),cp_pup_E(k),cp_mdn_E(l),w1,w3,w2,PHI_EPS);
              cUPDN =cUPDN - xi1 * cp_mdn_cddn(pd,stati,+l)* cp_mdn_cdup(mu,j,+l) * cp_pupmdn_cddn(md,k,+j) * cp_i_cdup(pu,k,+stati);
            enddo;enddo;enddo
            diml = dim_E_pup;
            dimk = dim_E_pupmdn;
            dimj = dim_E_mdn;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_mdn_cddn(md,stati,+j) *  cp_i_cdup(mu,l,+stati)
                  if (abs(cc) < cc_cutoff) cycle
                  do k =1, dimk;
               xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_pupmdn_E(k),cp_pup_E(l),w2,w3,w1,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_i_cdup(mu,l,+stati)* cp_pupmdn_cddn(pd,l,+k) * cp_mdn_cdup(pu,k,+j) * cp_mdn_cddn(md,stati,+j);
            enddo;enddo;enddo
            diml = dim_E_pup;
            dimk = dim_E_pupmdn;
            dimj = dim_E_pup;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_i_cdup(pu,j,+stati) *  cp_i_cdup(mu,l,+stati)
                  if (abs(cc) < cc_cutoff) cycle
            do k =1, dimk;
               xi1 = PhiM_ki(beta,cp_i_E(stati), cp_pup_E(j),cp_pupmdn_E(k),cp_pup_E(l),w3,w2,w1,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_i_cdup(mu,l,+stati)* cp_pupmdn_cddn(pd,l,+k) * cp_pupmdn_cddn(md,j,+k) * cp_i_cdup(pu,j,+stati);
            enddo;enddo;enddo
        endif

        if(bypass.or.(ndn/=0 .and. nup /= 0))then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_mupmdn;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mdn_cddn(md,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
               xi1 = PhiM_li(beta,cp_i_E(stati), cp_mupmdn_E(j),cp_mdn_E(k),cp_mdn_E(l),w2,w1,w3,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_mdn_cddn(pd,stati,+l)* cp_mupmdn_cdup(pu,l,+j) * cp_mupmdn_cdup(mu,k,+j) * cp_mdn_cddn(md,stati,+k);
            enddo;enddo;enddo
            diml = dim_E_mup;
            dimk = dim_E_mupmdn;
            dimj = dim_E_mup;
            do j =1, dimj;

               do l =1, diml;
                  cc =  cp_mup_cdup(mu,stati,+j) *  cp_mup_cdup(pu,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
            do k =1, dimk;
               xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mup_E(j),cp_mupmdn_E(k),cp_mup_E(l),w3,w1,w2,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_mup_cdup(pu,stati,+l)* cp_mupmdn_cddn(pd,l,+k) * cp_mupmdn_cddn(md,j,+k) * cp_mup_cdup(mu,stati,+j);
           enddo;enddo;enddo
            diml = dim_E_mup;
            dimk = dim_E_mupmdn;
            dimj = dim_E_mdn;
            do j =1, dimj;
               do l =1, diml;
                  cc =  cp_mdn_cddn(md,stati,+j) *  cp_mup_cdup(pu,stati,+l)
                 if (abs(cc) < cc_cutoff) cycle
            do k =1, dimk;
               xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_mupmdn_E(k),cp_mup_E(l),w1,w3,w2,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_mup_cdup(pu,stati,+l)* cp_mupmdn_cddn(pd,l,+k) * cp_mupmdn_cdup(mu,j,+k) * cp_mdn_cddn(md,stati,+j);
            enddo;enddo;enddo
            diml = dim_E_mdn;
            dimk = dim_E_mup;
            dimj = dim_E_mupmdn;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mup_cdup(mu,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
               xi1 = PhiM_li(beta,cp_i_E(stati), cp_mupmdn_E(j),cp_mup_E(k),cp_mdn_E(l),w2,w3,w1,PHI_EPS);
                cUPDN =cUPDN + xi1 * cp_mdn_cddn(pd,stati,+l)* cp_mupmdn_cdup(pu,l,+j) * cp_mupmdn_cddn(md,k,+j) * cp_mup_cdup(mu,stati,+k);
           enddo;enddo;enddo
            diml = dim_E_i;
            dimk = dim_E_mup;
            dimj = dim_E_mdn;
            do j =1, dimj;
            do k =1, dimk;

                  cc =  cp_mdn_cddn(md,stati,+j) *  cp_mup_cdup(pu,stati,+k)
                  if (abs(cc) < cc_cutoff) cycle
                  do l =1, diml;
               xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mup_E(k),cp_i_E(l),w3,w2,w1,PHI_EPS);
                cUPDN =cUPDN - xi1 * cp_mup_cdup(pu,stati,+k)* cp_mup_cdup(mu,l,+k) * cp_mdn_cddn(pd,l,+j) * cp_mdn_cddn(md,stati,+j);
            enddo;enddo;enddo
            diml = dim_E_mdn;
            dimk = dim_E_mup;
            dimj = dim_E_i;

            do k =1, dimk;
               do l =1, diml;
                  cc =  cp_mup_cdup(mu,stati,+k) *  cp_mdn_cddn(pd,stati,+l)
                  if (abs(cc) < cc_cutoff) cycle
                  do j =1, dimj;
               xi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mup_E(k),cp_mdn_E(l),w3,w2,w1,PHI_EPS);
                cUPDN = cUPDN - xi1 * cp_mdn_cddn(pd,stati,+l)* cp_mdn_cddn(md,j,+l) * cp_mup_cdup(pu,j,+k) * cp_mup_cdup(mu,stati,+k);
            enddo;enddo;enddo
        endif

        pUPDN = pUPDN + cUPDN*boltzZ;
    endif
    enddo

contains

subroutine run_tests(k_)
implicit none
integer :: k_

   write(*,*) 'test operator for obtaining density'
   write(*,*) 'orbital : ', k_
   do stati=1,dim_E_i
    boltzZ = exp(-(cp_i_E(stati)-gsE)*beta)/Z
    if(nup>0)then
    if(size(cp_mup_cdup,3)/=dim_E_mup.or.size(cp_mup_cdup,2)/=dim_E_i) then
      write(*,*) 'error shape'
      write(*,*) 'shape cp_i_cdup : ', shape(cp_mup_cdup)
      write(*,*) 'dimEi ' , dim_E_i
      write(*,*) 'dimE mup ', dim_E_mup
      stop
    endif
    endif
    if(op==2.and.nup>0.and.dim_E_mup>0) pUPDN=pUPDN + &
   & boltzZ*  sum(  (/( cp_mup_cdup(k_,stati,j) * cp_mup_cdup(k_,stati,j) , j=1,dim_E_mup  )/) )
    if(op==1.and.ndn>0.and.dim_E_mdn>0) pDNDN=pDNDN + &
   & boltzZ*  sum(  (/( cp_mdn_cddn(k_,stati,j) * cp_mdn_cddn(k_,stati,j) , j=1,dim_E_mdn  )/) )
   enddo
   if(op==1) write(*,*) 'up',nup,ndn,minval(cp_i_E(:))-gsE,real(pDNDN)
   if(op==2) write(*,*) 'dn',nup,ndn,minval(cp_i_E(:))-gsE,real(pUPDN)

 write(*,*) 'end test density'
 write(*,*) 'test commutators up/do'

 do stati=1,dim_E_i
  do stati_=1,dim_E_pupmdn
    if(op==2.and.ndn>0.and.nup<sites) then
       d1=sum(  (/( cp_mdn_cdup(k_,stati_,j)    * cp_mdn_cddn(k_,stati,j) , j=1,dim_E_mdn  )/) )
       d2=sum(  (/( cp_pupmdn_cddn(k_,j,stati_) * cp_i_cdup(k_,j,stati)   , j=1,dim_E_pup  )/) )
       write(*,*) 'commutator', nup , ndn , d1+d2
       pUPDN=pUPDN + d1 + d2
    endif
    if(op==1.and.ndn>0.and.dim_E_mdn>0) pDNDN=0.
  enddo
 enddo

 write(*,*) 'end test commutators'

end subroutine

end subroutine
end module
