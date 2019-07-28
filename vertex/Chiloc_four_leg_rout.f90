module chitilde
Use openmpmod
implicit none



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
 &    cutoff,cccc_cutoff,op,  PHI_EPS, beta, Z, gsE, sites, nup, ndn, &
 &    cp_i_E,            dim_E_i,          cp_pup_E,          dim_E_pup,     &
 &    cp_pdn_E,          dim_E_pdn,        cp_mup_E,          dim_E_mup,     &
 &    cp_mdn_E,          dim_E_mdn,        cp_p2dn_E,         dim_E_p2dn,    &
 &    cp_m2dn_E,         dim_E_m2dn,       cp_puppdn_E,       dim_E_puppdn,  &
 &    cp_muppdn_E,       dim_E_muppdn,     cp_pupmdn_E,       dim_E_pupmdn,  &
 &    cp_mupmdn_E,       dim_E_mupmdn,     cp_i_cdup,         cp_i_cddn,     &
 &    cp_pup_cddn,       cp_pdn_cdup,      cp_pdn_cddn,       cp_mup_cdup,   &
 &    cp_mup_cddn,       cp_mdn_cddn,      cp_mdn_cdup,       cp_muppdn_cdup,&
 &    cp_pupmdn_cddn,    cp_mupmdn_cdup,   cp_mupmdn_cddn,    cp_m2dn_cddn,  &
 &    norb,nomg, frequ_,rank,size2,chi_loc )

implicit none
integer    :: op,sites,nup,ndn,k_,l_,k__,l__
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
real(8)    :: cutoff
real(8)    :: d1,d2
real(8)    :: cccc, cccc_cutoff, cccct(norb,norb,norb,norb)
integer :: norb, nomg, iomg
complex(8) :: frequ_(nomg,3)
complex(8) :: chi_loc(norb,norb,norb,norb,nomg,2)
real(8) :: st,ft
logical :: bypass=.false., connect
integer :: rank,size2
integer :: iomg_mpi(0:size2)
!various notation for the same thing
u_=k_;
d_=l_;
iomg_mpi(0) = 1
k_ = nomg / size2
k__ = mod(nomg,size2)
do i =1,size2
   if(i<=k__) then
      iomg_mpi(i)  = iomg_mpi(i-1) +k_ +1
   else
      iomg_mpi(i)  = iomg_mpi(i-1)  +k_
   endif
enddo
          !$OMP PARALLEL DEFAULT(NONE)   SHARED(cp_i_E,cp_mdn_cddn,cp_mdn_E,nomg,cccc_cutoff,beta,&
          !$OMP dim_E_pdn,dim_E_i,cp_pdn_E,frequ_,PHI_EPS,boltzZ,cp_i_cddn,op,ndn,sites, cp_p2dn_E,dim_E_m2dn,&
          !$OMP dim_E_mdn,dim_E_p2dn,cp_pdn_cddn,cp_m2dn_cddn,cp_m2dn_E,nup,dim_e_pup,cp_i_cdup,cp_pup_e,dim_e_puppdn,&
          !$OMP cp_pup_cddn,cp_puppdn_e,cp_pdn_cdup,dim_e_muppdn,dim_e_mup,cp_mup_cdup,cp_mup_cddn,cp_muppdn_e,cp_mup_e,&
          !$OMP cp_muppdn_cdup,dim_e_pupmdn,cp_mdn_cdup,cp_pupmdn_e,cp_pupmdn_cddn,dim_e_mupmdn,cp_mupmdn_cdup,cp_mupmdn_e,&
          !$OMP  cp_mupmdn_cddn,gse,Z,cutoff,dimi,norb,rank,size2, iomg_mpi)&
          !$OMP reduction(+:chi_loc) &
          !$OMP PRIVATE(xi1,xi2,yi1,yi2,cccc,diml,dimk,dimj,iomg,k_,k__,l_,l__,cccct,connect,st,ft,j,k,l,stati)


   do stati=1,dim_E_i

    boltzZ = exp(-(cp_i_E(stati)-gsE)*beta)/Z;

    if(boltzZ<cutoff) cycle
    dimi = dim_E_i

    if(op == 1 )then

            if(ndn/=sites) then


            diml = dim_E_pdn
            dimk = dim_E_i;
            dimj = dim_E_pdn;
            !$OMP DO  SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
                  connect = .false.
                   do l__ = 1,norb
                      do l_ = 1, norb
                         do k__ = 1,norb
                            do k_ =1,norb
                               if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                               cccct(k_,k__,l_,l__)  = boltzZ * cp_i_cddn(k_,j,+stati)* cp_i_cddn(k_,j,+k) * cp_i_cddn(k_,l,+k) * cp_i_cddn(k_,l,+stati);
                               if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                            enddo
                         enddo
                      enddo
                   enddo

                   if (.not. connect) cycle
                   do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                      xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                      xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);

                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  * (xi1+xi2)
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
          enddo;enddo;enddo
          !$OMP END DO NOWAIT

          !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                         cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,l,+stati)* cp_i_cddn(k_,l,+k) * cp_i_cddn(k_,j,+k) * cp_i_cddn(k_,j,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle
                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                  yi1 = PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                  yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1+yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT

        endif


        if((ndn/=sites .and. ndn/=0))then
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;

            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
                  connect = .false.
                  do l__ = 1,norb
                     do l_ = 1, norb
                        do k__ = 1,norb
                           do k_ =1,norb
                              if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                              cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k_,stati,+j)* cp_mdn_cddn(k_,k,+j) * cp_i_cddn(k_,l,+k) * cp_i_cddn(k_,l,+stati)
                              if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                           enddo
                        enddo
                     enddo
                  enddo
                  if(.not. connect) cycle

                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                  xi1 =-PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                  xi2 = PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1) + cccct(k_,k__,l_,l__)  *  (xi1+xi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT

            diml = dim_E_mdn;
            dimk = dim_E_pdn;
            dimj = dim_E_i;


            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k_,stati,+l)* cp_mdn_cddn(k_,j,+l) * cp_i_cddn(k_,k,+j) * cp_i_cddn(k_,k,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                        yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);

                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                    chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1+yi2)
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo
             enddo;enddo;enddo
             !$OMP END DO NOWAIT

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,l,+stati)* cp_i_cddn(k_,l,+k) * cp_mdn_cddn(k_,k,+j) * cp_mdn_cddn(k_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                        xi2 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);

                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                    chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi1+xi2)
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo
             enddo;enddo;enddo
             !$OMP END DO NOWAIT

            diml = dim_E_i;
            dimk = dim_E_pdn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do l =1, diml;
                  do k =1, dimk;

                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,k,+stati)* cp_i_cddn(k_,k,+l) * cp_mdn_cddn(k_,l,+j) * cp_mdn_cddn(k_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        yi1 =PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                        yi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);

                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                    chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1+yi2)
                                 enddo
                              enddo
                           enddo
                        enddo

                enddo
             enddo;enddo;enddo
             !$OMP END DO NOWAIT

        endif

        if(ndn<sites-1)then
            diml = dim_E_p2dn;
            dimk = dim_E_pdn;
            dimj = dim_E_pdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,k,+stati)* cp_pdn_cddn(k_,l,+k) * cp_pdn_cddn(k_,l,+j) * cp_i_cddn(k_,j,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        yi1 =-PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                        yi2 = PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);

                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                    chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1+yi2)
                                 enddo
                              enddo
                           enddo
                        enddo

                enddo
             enddo;enddo;enddo
             !$OMP END DO NOWAIT

            diml = dim_E_pdn;
            dimk = dim_E_p2dn;
            dimj = dim_E_pdn;


            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
                  connect = .false.
                  do l__ = 1,norb
                     do l_ = 1, norb
                        do k__ = 1,norb
                           do k_ =1,norb
                              if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                              cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,j,+stati)* cp_pdn_cddn(k_,k,+j) * cp_pdn_cddn(k_,k,+l) * cp_i_cddn(k_,l,+stati)
                              if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                           enddo
                        enddo
                     enddo
                  enddo
                  if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                     xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                     xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);

                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi1+xi2)
                              enddo
                           enddo
                        enddo
                     enddo

                  enddo
             enddo;enddo;enddo
             !$OMP END DO NOWAIT

        endif

        if(ndn>1)then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_m2dn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
                  connect = .false.
                  do l__ = 1,norb
                     do l_ = 1, norb
                        do k__ = 1,norb
                           do k_ =1,norb
                              if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                              cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k_,stati,+l)* cp_m2dn_cddn(k_,l,+j) * cp_m2dn_cddn(k_,k,+j) * cp_mdn_cddn(k_,stati,+k)
                              if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                           enddo
                        enddo
                     enddo
                  enddo
                  if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                     yi1 =-PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                     yi2 = PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);

                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1+yi2)
                              enddo
                           enddo
                        enddo
                     enddo

                  enddo
             enddo;enddo;enddo
             !$OMP END DO NOWAIT

            diml = dim_E_mdn;
            dimk = dim_E_m2dn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
                  connect = .false.
                  do l__ = 1,norb
                     do l_ = 1, norb
                        do k__ = 1,norb
                           do k_ =1,norb
                              if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                              cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k_,stati,+l)* cp_m2dn_cddn(k_,l,+k) * cp_m2dn_cddn(k_,j,+k) * cp_mdn_cddn(k_,stati,+j)
                              if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                           enddo
                        enddo
                     enddo
                  enddo
                  if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                     yi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                     yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);

                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1+yi2)
                              enddo
                           enddo
                        enddo
                     enddo

                  enddo
             enddo;enddo;enddo
             !$OMP END DO NOWAIT
        endif

        if(ndn/=0)then

            diml = dim_E_i;
            dimk = dim_E_mdn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k_,stati,+k) * cp_mdn_cddn(k_,l,+k) * cp_mdn_cddn(k_,l,+j) * cp_mdn_cddn(k_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                        xi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);

                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                    chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi1+xi2)
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo
                  enddo;enddo;enddo
             !$OMP END DO NOWAIT
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_i;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k_,stati,+l) * cp_mdn_cddn(k_,j,+l) * cp_mdn_cddn(k_,j,+k) * cp_mdn_cddn(k_,stati,+k)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                        yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__==k_.and.l__==l_.and.k_==l__)) cycle
                                    chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1+yi2)
                                 enddo
                              enddo
                           enddo
                        enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
        endif





     endif

   if(op == 1)then


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

        if(ndn/=sites) then
            diml = dim_E_pdn ;
            dimk = dim_E_i;
            dimj = dim_E_pdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,j,+stati)* cp_i_cddn(k__,j,+k) * cp_i_cddn(l_,l,+k) * cp_i_cddn(l__,l,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle
                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                    chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi1)
                                 enddo
                              enddo
                           enddo
                        enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
               do j =1, dimj;
                  do k =1, dimk;
                     do l =1, diml;
                        connect = .false.
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                    cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,j,+stati)* cp_i_cddn(k__,j,+k) * cp_i_cddn(k_,l,+k) * cp_i_cddn(l__,l,+stati)
                                    if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                                 enddo
                              enddo
                           enddo
                        enddo
                        if(.not. connect) cycle


                        do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                           xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                           do l__ = 1,norb
                              do l_ = 1, norb
                                 do k__ = 1,norb
                                    do k_ =1,norb
                                       if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                       chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi2)
                                    enddo
                                 enddo
                              enddo
                           enddo

                        enddo
                     enddo;enddo;enddo
                     !$OMP END DO NOWAIT

               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,l,+stati)* cp_i_cddn(l__,l,+k) * cp_i_cddn(k_,j,+k) * cp_i_cddn(k__,j,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi1 = PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  ( yi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,l,+stati)* cp_i_cddn(l__,l,+k) * cp_i_cddn(l_,j,+k) * cp_i_cddn(k_,j,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_pdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (+yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT
        endif


        if((ndn/=sites .and. ndn/=0))then
            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;

            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k__,stati,+j)* cp_mdn_cddn(k_,k,+j) * cp_i_cddn(l_,l,+k) * cp_i_cddn(l__,l,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  xi1 =-PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1) + cccct(k_,k__,l_,l__)  * (xi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT

               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k__,stati,+j)* cp_mdn_cddn(l_,k,+j) * cp_i_cddn(k_,l,+k) * cp_i_cddn(l__,l,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1


                  xi2 = PhiM_ii(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1) + cccct(k_,k__,l_,l__)  * (xi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT

            diml = dim_E_mdn;
            dimk = dim_E_pdn;
            dimj = dim_E_i;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_mdn_cddn(l_,j,+l) * cp_i_cddn(k_,k,+j) * cp_i_cddn(k__,k,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_mdn_cddn(k_,j,+l) * cp_i_cddn(l_,k,+j) * cp_i_cddn(k__,k,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pdn_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (+yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,l,+stati)* cp_i_cddn(l__,l,+k) * cp_mdn_cddn(k__,k,+j) * cp_mdn_cddn(k_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  xi1 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT

               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,l,+stati)* cp_i_cddn(l__,l,+k) * cp_mdn_cddn(k__,k,+j) * cp_mdn_cddn(l_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  xi2 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (+xi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo

               enddo;enddo;enddo
               !$OMP END DO NOWAIT

            diml = dim_E_i;
            dimk = dim_E_pdn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,k,+stati)* cp_i_cddn(k__,k,+l) * cp_mdn_cddn(l__,l,+j) * cp_mdn_cddn(l_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,k,+stati)* cp_i_cddn(k__,k,+l) * cp_mdn_cddn(l__,l,+j) * cp_mdn_cddn(k_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pdn_E(k),cp_i_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
        endif

        if(ndn<sites-1)then

            diml = dim_E_p2dn;
            dimk = dim_E_pdn;
            dimj = dim_E_pdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,k,+stati)* cp_pdn_cddn(l_,l,+k) * cp_pdn_cddn(l__,l,+j) * cp_i_cddn(k__,j,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi1 =-PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  ( yi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,k,+stati)* cp_pdn_cddn(k_,l,+k) * cp_pdn_cddn(l__,l,+j) * cp_i_cddn(k__,j,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi2 = PhiM_ji(beta,cp_i_E(stati), cp_pdn_E(j),cp_pdn_E(k),cp_p2dn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (+yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
            diml = dim_E_pdn;
            dimk = dim_E_p2dn;
            dimj = dim_E_pdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,j,+stati)* cp_pdn_cddn(k_,k,+j) * cp_pdn_cddn(k__,k,+l) * cp_i_cddn(l__,l,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(k_,j,+stati)* cp_pdn_cddn(l_,k,+j) * cp_pdn_cddn(k__,k,+l) * cp_i_cddn(l__,l,+stati)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  xi2 =-PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_p2dn_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
        endif

        if(ndn>1)then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_m2dn;

            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_m2dn_cddn(k__,l,+j) * cp_m2dn_cddn(k_,k,+j) * cp_mdn_cddn(l_,stati,+k)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi1 =-PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_m2dn_cddn(k__,l,+j) * cp_m2dn_cddn(l_,k,+j) * cp_mdn_cddn(k_,stati,+k)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi2 = PhiM_li(beta,cp_i_E(stati), cp_m2dn_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT

            diml = dim_E_mdn;
            dimk = dim_E_m2dn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k__,stati,+l)* cp_m2dn_cddn(l__,l,+k) * cp_m2dn_cddn(l_,j,+k) * cp_mdn_cddn(k_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  ( yi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k__,stati,+l)* cp_m2dn_cddn(l__,l,+k) * cp_m2dn_cddn(k_,j,+k) * cp_mdn_cddn(l_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi2 =-PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_m2dn_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) = chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (+yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
        endif

        if(ndn/=0)then

            diml = dim_E_i;
            dimk = dim_E_mdn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k__,stati,+k) * cp_mdn_cddn(l_,l,+k) * cp_mdn_cddn(l__,l,+j) * cp_mdn_cddn(k_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(k__,stati,+k) * cp_mdn_cddn(k_,l,+k) * cp_mdn_cddn(l__,l,+j) * cp_mdn_cddn(l_,stati,+j)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  xi2 =-PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mdn_E(k),cp_i_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (xi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT

            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_i;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l) * cp_mdn_cddn(k_,j,+l) * cp_mdn_cddn(k__,j,+k) * cp_mdn_cddn(l_,stati,+k)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1

                  yi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  ( yi1)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
            do k =1, dimk;
               do l =1, diml;
	connect = .false.
	do l__ = 1,norb
            do l_ = 1, norb
                  do k__ = 1,norb
                     do k_ =1,norb
                        if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                        cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l) * cp_mdn_cddn(l_,j,+l) * cp_mdn_cddn(k__,j,+k) * cp_mdn_cddn(k_,stati,+k)
                         if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                      enddo
                   enddo
                enddo
             enddo
             if(.not. connect) cycle


                  do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                  yi2 =-PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                      do l__ = 1,norb
                         do l_ = 1, norb
                            do k__ = 1,norb
                               do k_ =1,norb
                                  if(.not.(k__/=k_.or.l__/=l_.or.k_/=l__)) cycle
                                  chi_loc(k_,k__,l_,l__,iomg,1) =chi_loc(k_,k__,l_,l__,iomg,1)+ cccct(k_,k__,l_,l__)  *  (+yi2)
                               enddo
                            enddo
                         enddo
                      enddo

                  enddo
               enddo;enddo;enddo
               !$OMP END DO NOWAIT
        endif


     endif


! **** UPDN *****************************************************************************

    if(op == 2)then
        if((ndn/=sites .and. nup /= sites))then

            dimj = dim_E_pup;
            dimk = dim_E_i;
            diml = dim_E_pdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,l,+stati)* cp_i_cddn(l__,l,+k) * cp_i_cdup(k_,j,+k) * cp_i_cdup(k__,j,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle

                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pup_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo
                  enddo;enddo;enddo
                  !$OMP END DO NOWAIT

            diml = dim_E_puppdn;
            dimk = dim_E_pup;
            dimj = dim_E_pup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cdup(k_,k,stati)* cp_pup_cddn(l_,l,+k) * cp_pup_cddn(l__,l,+j) * cp_i_cdup(k__,j,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle

                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ji(beta,cp_i_E(stati), cp_pup_E(j),cp_pup_E(k),cp_puppdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo
                  enddo;enddo;enddo
                  !$OMP END DO NOWAIT
                  !print*,3,chi_loc(1,1,1,1,1,2)/boltzZ,stati

            diml = dim_E_pdn;
            dimk = dim_E_puppdn;
            dimj = dim_E_pdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,j,+stati)* cp_pdn_cdup(k_,k,+j) * cp_pdn_cdup(k__,k,+l) * cp_i_cddn(l__,l,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle




                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1   = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_puppdn_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) = chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo
                  enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !print*,4,chi_loc(1,1,1,1,1,2)/boltzZ,stati

            diml = dim_E_pdn;
            dimk = dim_E_puppdn;
            dimj = dim_E_pup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cdup(k_,j,stati)* cp_pup_cddn(l_,k,+j) * cp_pdn_cdup(k__,k,+l) * cp_i_cddn(l__,l,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pup_E(j),cp_puppdn_E(k),cp_pdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - xi1 * cccct(k_,k__,l_,l__)
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo

                  enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !print*,5,chi_loc(1,1,1,1,1,2)/boltzZ,stati

            diml = dim_E_puppdn;
            dimk = dim_E_pdn;
            dimj = dim_E_pup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,k,+stati)* cp_pdn_cdup(k_,l,+k) * cp_pup_cddn(l__,l,+j) * cp_i_cdup(k__,j,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ji(beta,cp_i_E(stati), cp_pup_E(j),cp_pdn_E(k),cp_puppdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,6,chi_loc(1,1,1,1,1,2)/boltzZ,stati
        endif

        if((ndn/=sites .and. nup /= 0)) then
            diml = dim_E_muppdn;
            dimk = dim_E_mup;
            dimj = dim_E_mup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mup_cdup(k__,stati,+k)* cp_mup_cddn(l_,l,+k) * cp_mup_cddn(l__,l,+j) * cp_mup_cdup(k_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mup_E(j),cp_mup_E(k),cp_muppdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,7,chi_loc(1,1,1,1,1,2)/boltzZ,stati

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mup_cdup(k__,stati,+j)* cp_mup_cdup(k_,k,+j) * cp_i_cddn(l_,l,+k) * cp_i_cddn(l__,l,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ii(beta,cp_i_E(stati), cp_mup_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
            enddo;enddo;enddo
            !$OMP END DO NOWAIT
            !print*,8,chi_loc(1,1,1,1,1,2)/boltzZ,stati

            diml = dim_E_pdn;
            dimk = dim_E_i;
            dimj = dim_E_mup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,l,stati)* cp_i_cddn(l__,l,+k) * cp_mup_cdup(k__,k,+j) * cp_mup_cdup(k_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle

                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mup_E(j),cp_i_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


                  enddo;enddo;enddo
                  !$OMP END DO NOWAIT
                  !print*,9,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_muppdn;
            dimk = dim_E_pdn;
            dimj = dim_E_mup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ * cp_i_cddn(l_,k,+stati)* cp_muppdn_cdup(k__,k,+l) * cp_mup_cddn(l__,l,+j) * cp_mup_cdup(k_,stati,+j);
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mup_E(j),cp_pdn_E(k),cp_muppdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !print*,10,chi_loc(1,1,1,1,1,2)/boltzZ,stati

            diml = dim_E_pdn;
            dimk = dim_E_muppdn;
            dimj = dim_E_mup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mup_cdup(k__,stati,+j)* cp_mup_cddn(l_,k,+j) * cp_muppdn_cdup(k_,l,+k) * cp_i_cddn(l__,l,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ii(beta,cp_i_E(stati), cp_mup_E(j),cp_muppdn_E(k),cp_pdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


                  enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,11,chi_loc(1,1,1,1,1,2)/boltzZ,stati

            diml = dim_E_pdn;
            dimk = dim_E_muppdn;
            dimj = dim_E_pdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cddn(l_,j,+stati)* cp_muppdn_cdup(k__,j,+k) * cp_muppdn_cdup(k_,l,+k) * cp_i_cddn(l__,l,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ii(beta,cp_i_E(stati), cp_pdn_E(j),cp_muppdn_E(k),cp_pdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  *xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


                  enddo;enddo;enddo
                  !$OMP END DO NOWAIT
                  !print*,12,chi_loc(1,1,1,1,1,2)/boltzZ,stati
        endif

        if((ndn/=0 .and. nup /= sites))then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_pupmdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_mdn_cdup(k_,j,+l) * cp_mdn_cdup(k__,j,+k) * cp_mdn_cddn(l_,stati,+k)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_li(beta,cp_i_E(stati), cp_pupmdn_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,2),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,13,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_i;
            dimk = dim_E_pup;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cdup(k_,k,+stati)* cp_i_cdup(k__,k,+l) * cp_mdn_cddn(l__,l,+j) * cp_mdn_cddn(l_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_pup_E(k),cp_i_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,14,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_mdn;
            dimk = dim_E_pup;
            dimj = dim_E_i;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_mdn_cddn(l_,j,+l) * cp_i_cdup(k_,k,+j) * cp_i_cdup(k__,k,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_pup_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,15,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_mdn;
            dimk = dim_E_pup;
            dimj = dim_E_pupmdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,l)* cp_mdn_cdup(k_,j,+l) * cp_pupmdn_cddn(l_,k,+j) * cp_i_cdup(k__,k,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle


                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_li(beta,cp_i_E(stati), cp_pupmdn_E(j),cp_pup_E(k),cp_mdn_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) -  cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


                  enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !print*,16,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_pup;
            dimk = dim_E_pupmdn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cdup(k_,l,+stati)* cp_pupmdn_cddn(l__,l,+k) * cp_mdn_cdup(k__,k,+j) * cp_mdn_cddn(l_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle

                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_pupmdn_E(k),cp_pup_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo
            enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,17,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_pup;
            dimk = dim_E_pupmdn;
            dimj = dim_E_pup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_i_cdup(k_,l,stati)* cp_pupmdn_cddn(l__,l,+k) * cp_pupmdn_cddn(l_,j,+k) * cp_i_cdup(k__,j,+stati)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ki(beta,cp_i_E(stati), cp_pup_E(j),cp_pupmdn_E(k),cp_pup_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


                  enddo;enddo;enddo
                  !$OMP END DO NOWAIT
                  !print*,17,chi_loc(1,1,1,1,1,2)/boltzZ,stati
        endif

        if((ndn/=0 .and. nup /= 0))then
            diml = dim_E_mdn;
            dimk = dim_E_mdn;
            dimj = dim_E_mupmdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_mupmdn_cdup(k__,l,+j) * cp_mupmdn_cdup(k_,k,+j) * cp_mdn_cddn(l_,stati,+k)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle




                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_li(beta,cp_i_E(stati), cp_mupmdn_E(j),cp_mdn_E(k),cp_mdn_E(l),frequ_(iomg,2),frequ_(iomg,1),frequ_(iomg,3),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


                  enddo;enddo;enddo
                  !$OMP END DO NOWAIT
                  !print*,18,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_mup;
            dimk = dim_E_mupmdn;
            dimj = dim_E_mup;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mup_cdup(k__,stati,+l)* cp_mupmdn_cddn(l__,l,+k) * cp_mupmdn_cddn(l_,j,+k) * cp_mup_cdup(k_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mup_E(j),cp_mupmdn_E(k),cp_mup_E(l),frequ_(iomg,3),frequ_(iomg,1),frequ_(iomg,2),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,19,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_mup;
            dimk = dim_E_mupmdn;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mup_cdup(k__,stati,l)* cp_mupmdn_cddn(l__,l,+k) * cp_mupmdn_cdup(k_,j,+k) * cp_mdn_cddn(l_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ki(beta,cp_i_E(stati), cp_mdn_E(j),cp_mupmdn_E(k),cp_mup_E(l),frequ_(iomg,1),frequ_(iomg,3),frequ_(iomg,2),PHI_EPS);

                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


               enddo;enddo;enddo
               !$OMP END DO NOWAIT
               !print*,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_mdn;
            dimk = dim_E_mup;
            dimj = dim_E_mupmdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_mupmdn_cdup(k__,l,+j) * cp_mupmdn_cddn(l_,k,+j) * cp_mup_cdup(k_,stati,+k)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_li(beta,cp_i_E(stati), cp_mupmdn_E(j),cp_mup_E(k),cp_mdn_E(l),frequ_(iomg,2),frequ_(iomg,3),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) + cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_i;
            dimk = dim_E_mup;
            dimj = dim_E_mdn;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mup_cdup(k__,stati,+k)* cp_mup_cdup(k_,l,+k) * cp_mdn_cddn(l__,l,+j) * cp_mdn_cddn(l_,stati,+j)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle



                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_ji(beta,cp_i_E(stati), cp_mdn_E(j),cp_mup_E(k),cp_i_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) =chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo

                     enddo


         enddo;enddo;enddo
         !$OMP END DO NOWAIT
         !print*,chi_loc(1,1,1,1,1,2)/boltzZ,stati
            diml = dim_E_mdn;
            dimk = dim_E_mup;
            dimj = dim_E_i;
            !$OMP DO SCHEDULE(DYNAMIC)
            do j =1, dimj;
               do k =1, dimk;
                  do l =1, diml;
                     connect = .false.
                     do l__ = 1,norb
                        do l_ = 1, norb
                           do k__ = 1,norb
                              do k_ =1,norb
                                 cccct(k_,k__,l_,l__)  = boltzZ *  cp_mdn_cddn(l__,stati,+l)* cp_mdn_cddn(l_,j,+l) * cp_mup_cdup(k__,j,+k) * cp_mup_cdup(k_,stati,+k)
                                 if(abs(cccct(k_,k__,l_,l__)) >cccc_cutoff) connect = .true.
                              enddo
                           enddo
                        enddo
                     enddo
                     if(.not. connect) cycle
                     do iomg=iomg_mpi(rank),iomg_mpi(rank+1)-1
                        xi1 = PhiM_li(beta,cp_i_E(stati), cp_i_E(j),cp_mup_E(k),cp_mdn_E(l),frequ_(iomg,3),frequ_(iomg,2),frequ_(iomg,1),PHI_EPS);
                        do l__ = 1,norb
                           do l_ = 1, norb
                              do k__ = 1,norb
                                 do k_ =1,norb
                                    chi_loc(k_,k__,l_,l__,iomg,2) = chi_loc(k_,k__,l_,l__,iomg,2) - cccct(k_,k__,l_,l__)  * xi1
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo;enddo;enddo
            !$OMP END DO NOWAIT

        endif

     endif

 enddo
 !$OMP END PARALLEL

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
