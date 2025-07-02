! ====================================================================
subroutine read_restart
  use o_ARRAYS

  implicit none


  write(*,*) 'Reading restart'
  !
  ! Read restart
  open(25,file='restart.out', form='unformatted')
     read(25) eta_n, eta_n_1, eta_n_2, U_n_2D, U_n_1, U_n_2
     read(25) U_n, V_n, U_rhsAB, V_rhsAB, Wvel
     read(25) TF, SF, T_old, S_old
     if (ver_mix == 2) read(25) tke, Av_node, teps, Kv
     if (ver_mix == 3) read(25) bt, snu, Kv
     read(25) ini_time
	 !VF The nsteps changes (decreasing) if you use hot_start 
	 nsteps=nsteps-ini_time
     ini_time=ini_time+1
     read(25) time_jd, time_jd0, dt_restart, dt_2D_restart
     read(25) T_counter
  close(25)

   write(*,*) 'T_counter', T_counter
  write(*,*) 'ini_time', ini_time
  time_jd0=time_jd
  write(*,*) 'ini_time_jd', time_jd0
  lfirst=.false.
  
end subroutine read_restart

! ====================================================================

subroutine initial_state

  use o_MESH
  use o_ARRAYS
  use o_PARAM
  !aaa  use fv_ic

  use g_parsup
  use g_comm_auto

  implicit none

  integer       :: i,el, elnodes(4), n, m, j, jjj, nl, nz, c_riv, eledges(4), flag, elem
  real(kind=WP) :: x1, x2, y1, y2, amp, a, den_0
  real(kind=WP) :: dx, dy, period, hu, x, y
  real(kind=WP) :: pr, tt, ss, pp, d
  real(kind=WP) :: de, ga, ep, fe
  real(kind=WP), external :: theta
  real(kind=WP), allocatable, dimension(:,:)  ::   aux
  real, allocatable, dimension(:,:)  ::  S_t, T_t
  real(kind=WP) :: dmean

  integer :: node_size	
  real(kind=WP) :: mx_eta, mn_eta, mx_dep, mn_dep
  real(kind=WP) :: mx_t, mn_t, mx_s, mn_s
  integer :: ierror

  integer :: ttldim

  real(kind=WP) :: all_elem_area(elem2D)

  real(kind=WP), allocatable, dimension(:) :: Rj_buff,aj_buff, bj_buff
  real(kind=WP)                            :: Rval, aval, bval
  integer, allocatable, dimension(:) :: mapping
  integer                            :: iofs

  node_size=myDim_nod2D+eDim_nod2D

  !====================================
  ! Calculating of reference density	
  !====================================

#ifdef USE_MPI

  allocate(Rj_buff(nod2D), aj_buff(nod2D), bj_buff(nod2D))
  allocate(mapping(nod2D))
  Rj_buff(:)=0.0_WP; aj_buff(:)=0.0_WP; bj_buff(:)=0.0_WP
  mapping(:)=0

  do n=1, myDim_nod2D+eDim_nod2D
     iofs=myList_nod2D(n)
     mapping(iofs)=n
  end do


  if (mype==0) print *,'Gathering elements and calculating scale_area ...'
  
  call gather_elem(elem_area,all_elem_area)
  scale_area=sum(all_elem_area)/elem2D

  call MPI_BCast(scale_area, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)

#else
  scale_area=sum(elem_area)/elem2D
#endif

  if (mype==0) write(*,*) 'scale_area = ',scale_area
 
  !====================================
  ! Calculating coeff. to determine proper eta (see fv_average_dynamic) 
  ! in case of absence or presence tidal potential.
  ! Attention, in case of presence tidal potential, a_tp can be also calculating 
  ! for each time step to take into account the influence of internal waves presence 
  ! to tidal potential energy	
  !=====================================

  allocate(a_tp(node_size,2))
  if (T_potential) then
     a_tp(:,1)=0.918_WP
     a_tp(:,2)=0.69_WP
  else
     a_tp(:,1)=1.0_WP
     a_tp(:,2)=0.0_WP
  endif
  do i=1,node_size
     if (index_nod2D(i) == 2) then
        a_tp(i,1)=1.0_WP
        a_tp(i,2)=0.0_WP
     endif
  enddo

  ! =========
  ! initialize with zeros
  ! =========

  eta_n=0.0_WP
  U_n_2D=0.0_WP

  eta_n_1=0.0_WP
  eta_n_2=0.0_WP
  U_n_1=0.0_WP
  U_n_2=0.0_WP

!SHTEST
! Constant wind
!!$print *,'CONSTANT WIND FIELD'
!!$taux(:) = -0.05_WP
!!$tauy(:) = 0.0_WP

  do n=1,node_size
     d = depth(n) + eta_n(n)
     !print *, depth(n), eta_n(n), d
     if (d < 0.0_WP) then
        eta_n(n) = -depth(n)
        eta_n_1(n) = -depth(n)
        eta_n_2(n) = -depth(n)
     endif
     !  print*, eta_n(n)
  enddo

  !+++++++++++++++++++++++++++++
  ! for compute mask WaD
  !     we need a etaAB
  !+++++++++++++++++++++++++++++

  de=0.614_WP
  ga=0.088_WP
  ep=0.013_WP
  fe=1.0_WP-de-ga-ep
  ! Interpolate AM4
  ssh_rhs = eta_n
  etaAB=de*ssh_rhs+fe*eta_n+ga*eta_n_1+ep*eta_n_2  


#ifdef USE_MPI
  call MPI_REDUCE(maxval(eta_n), mx_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(minval(eta_n), mn_eta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(maxval(depth), mx_dep, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(minval(depth), mn_dep, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       0, MPI_COMM_FESOM_C, MPIerr)

  if (mype==0) then
     write(*,*) 'MAX_MIN_SSH_initial= ', mx_eta, mn_eta
     write(*,*) 'MAX_MIN_depth= ', mx_dep, mn_dep
  end if

#else

  write(*,*) 'MAX_MIN_SSH_initial= ', maxval(eta_n),minval(eta_n)
  write(*,*) 'MAX_MIN_depth= ', maxval(depth),minval(depth)

#endif

  ! AA initialzation for sediment
  do n=1,node_size
     dmean=max(Dmin,eta_n(n) + depth(n))
     con(n) = 0.0_WP/dmean
     con_bc(n) = con(n)
  enddo
  ! AA 

  !VF initialization of sediment module and taking into account
  !presence of sediments into the density calculation

!aa67  if (comp_sediment) call comp_sediment_ini

  density_0=1000.0_WP
  den_0=24.7630052361001_WP

  !VF, calculate density using info about T_const, S_const and Concentration of sus. sediments
!aaa   if (comp_sediment) then
!SH skipped for now max/min CF must be properly computed!!
!aaa      call densityJM(T_const, S_const, 0.0, den_0,(maxval(CF)+minval(CF))*0.5_WP)
!aaa   else
!aaa      call densityJM(T_const, S_const, 0.0, den_0)
!aaa   endif

 
  if (mype==0) write(*,*) 'den_0 = ', den_0
  density_0=density_0+den_0

  if (mype==0) write(*,*) 'density_0', density_0


  if (TF_presence) then  

     if (mype==0) then
        open(50,file=trim(meshpath)//'m2.out', status='old')
        read(50,*) n
     end if

#ifdef USE_MPI
     call MPI_BCast(n, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
#endif

     allocate(ampt(n,15),fazt(n,15))
     ampt=0.0_WP
     fazt=0.0_WP

     if (mype==0) then
        allocate(aux(n,n_th*2))
        do i=1,n
           read(50,*) jjj,aux(i,:)
           do j=2,n_th*2,2
              if (aux(i,j) < 0.0_WP) aux(i,j) = 360.0_WP + aux(i,j)
           end do
        enddo
        close(50)
        j=1
        do i=1,15
           if (a_th(i)==1) then
              ampt(:,i)=aux(:,j)
              fazt(:,i)=aux(:,j+1)*rad
              j=j+2
           endif
        end do
        write (*,*) 'max_amplitude', maxval(ampt)
        deallocate(aux)
     end if

#ifdef USE_MPI
     do i=1,n
        call MPI_BCast(ampt(i,1:15), 15, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        call MPI_BCast(fazt(i,1:15), 15, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     end do
#endif


  endif ! TF_presence
  
  !
  ! Irradiance in the upper ocean

  !
  if (type_task>2) then  

!SH OBACHT INITIALIZE JER
! NOT YET MPI



     if (.not.(jer_const)) then

        allocate(aj(node_size),bj(node_size),Rj(node_size))

#ifdef USE_MPI

        if (mype==0) then

           open(50,file=trim(meshpath)//'jer.out', status='old')

           do j=1,nod2D
              read(50,*) Rj_buff(j),aj_buff(j),bj_buff(j)
           enddo

        end if

        call MPI_BCast(Rj_buff(:), nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        call MPI_BCast(aj_buff(:), nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        call MPI_BCast(bj_buff(:), nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
          
        do j=1,nod2D

           Rval=Rj_buff(n)
           aval=aj_buff(n)
           bval=bj_buff(n)

           if (mapping(j)>0) then
              Rj(mapping(j))=Rval
              aj(mapping(j))=aval
              bj(mapping(j))=bval
           end if
        end do
        
        if (mype==0) close(50)

#else
        open(50,file=trim(meshpath)//'jer.out', status='old')

        do j=1,nod2D
           read(50,*) Rj(j),aj(j),bj(j)
        enddo

        close(50)	

#endif

     endif

     !aaa  if (key_ic) then
     !aaa     call ic_do  ! from fv_ic.f90
     !aaa  else
     !aaa     if (TEMP_ON) then
     !aaa      open(50,file=trim(meshpath)//trim(TITLE)//'_temp.out', status='old')
     !aaa      do n=1,nod2D
     !aaa       read(50,*) (TF(nz,n),nz=1,nsigma-1)
     !aaa      enddo
     !aaa      close(50)
     !aaa      else
     !aaa       TF = T_const
     !aaa      endif
     !aaa      if (SALINITY_ON) then
     !aaa      open(50,file=trim(meshpath)//trim(TITLE)//'_sal.out', status='old')
     !aaa      do n=1,nod2D
     !aaa       read(50,*) (SF(nz,n),nz=1,nsigma-1)
     !aaa      enddo
     !aaa      close(50)
     !aaa      else
     !aaa       SF = S_const
     !aaa     endif
     !aaa  endif
!!!!!!!!!!!!!!!!!!!!!
     ! LE
     !aaa  if (key_ic) then
     !aaa     call ic_do  ! from fv_ic.f90
     !aaa  else

     if (TEMP_ON) then
        !aa experiment lock exchange ONLY
        do n=1,node_size
           a = coord_nod2D(1,n)
           do nl=1,nsigma-1
              if (a < 0.0_WP) then
                 TF(nl,n) = 5.0_WP/0.28_WP
              else
                 !  if (a > 1.E-4*rad) then
                 TF(nl,n) = 10.0_WP/0.28_WP
                 !  else
                 !                        TF(nl,n) = 0.5_WP*(5.0_WP + 10.0_WP)/0.28_WP
              endif
           enddo
        enddo
     endif
     !aa END for lock exchange experiment

     !
     ! Climatology
     !
     !   convert in situ temperature into potential temperature
     !
     pr=0.0_WP

     do n=1,node_size
        do nz=1,nsigma-1
           tt = TF(nz,n)
           ss = SF(nz,n)
           pp = Z(nz,n)
           !aaa     TF(nz,n)=theta(ss,tt,pp,pr)
        enddo
     enddo
     T_old= TF
     S_old= SF
     Tclim = TF
     Sclim = SF

#ifdef USE_MPI
     
     call MPI_REDUCE(maxval(TF), mx_t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          0, MPI_COMM_FESOM_C, MPIerr)
     call MPI_REDUCE(minval(TF), mn_t, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          0, MPI_COMM_FESOM_C, MPIerr)
     call MPI_REDUCE(maxval(SF), mx_s, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          0, MPI_COMM_FESOM_C, MPIerr)
     call MPI_REDUCE(minval(SF), mn_s, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          0, MPI_COMM_FESOM_C, MPIerr)

     if (mype==0) then
        write(*,*)'T_max, T_min', mx_t, mn_t
        write(*,*)'S_max, S_min', mx_s, mn_s
     end if

#else

     write(*,*)'T_max, T_min', maxval(TF), minval(TF)
     write(*,*)'S_max, S_min', maxval(SF), minval(SF)      

#endif
  endif
  !VF: Initialisation of C_d depends on user choice
  if (BFC==2) then
     DO elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        dmean = max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))
        C_d_el(elem)=1.0_WP/(log(dmean/(z0b_min+za))/cka)**2
     enddo
#ifdef USE_MPI
     call exchange_elem(C_d_el)
#endif
  else
     C_d_el=C_d
  endif

  if (mype==0) write(*,*)'finish initial_state'

  !a call solve_tracer_second_order_rec(SF, 's')
  !a end if

end subroutine initial_state

! ====================================================================

subroutine initial_riv

  !VF,initial_riv
  use o_MESH
  use o_ARRAYS
  use o_PARAM

  use g_parsup

  implicit none

  integer       :: i,el, n,m,j, nz, eledges(4), flag, k, ed(2), elem
  real(kind=WP) :: xv, yv, aux, auy, xc(2), n1(2), n2(2), A1, A2, B1, B2, C1, C2, amin
  real(kind=WP) :: ss,tt,pp, pr
  real(kind=WP), external :: theta
!aa67  real(kind=WP), allocatable, dimension(:,:)  ::  S_t, T_t

  integer, allocatable, dimension(:) :: mapping
  integer :: ipos
  integer :: ierror

!aa67  if (riv_OB) then
!aa67     if (mype==0) then
!aa67        open(50,file=trim(meshpath)//trim(TITLE)//'_riv_ob.out', status='old')
!aa67        read(50,*) n
!aa67        read(50,*) m
!aa67        close(50)
!aa67     end if

!aa67#ifdef USE_MPI
!aa67     call MPI_BCast(n, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
!aa67     call MPI_BCast(m, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
!aa67#endif
     
!aa67     allocate(riv_node_ob(n), riv_elev(n),riv_vt2(m),riv_vt_jd(m+1)) 
!aa67     allocate(Tr_distr2(nsigma-1,n),Sr_distr2(nsigma-1,n))
!aa67     allocate(Tr_distr_t2(nsigma-1,n,m),Sr_distr_t2(nsigma-1,n,m), riv_elev_t(n,m))


!aa67     if (mype==0) then
!aa67        open(50,file=trim(meshpath)//trim(TITLE)//'_riv_ob.out', status='old')
!aa67        read(50,*) n
!aa67        read(50,*) m

!aa67        do j=1,m
!aa67           read(50,*) riv_vt2(j)
!aa67           do i=1,n
!aa67              read(50,*) riv_node_ob(i),riv_elev_t(i,j),Sr_distr_t2(1:nsigma-1,i,j),Tr_distr_t2(1:nsigma-1,i,j)
!aa67           enddo
!aa67        enddo
!aa67        close(50)

!aa67        Tr_distr2(:,:)=Tr_distr_t2(:,:,1)
!aa67        Sr_distr2(:,:)=Sr_distr_t2(:,:,1) 
!aa67        riv_elev = riv_elev_t(:,1)
!aa67        riv_amount_ob=n
!aa67        pr=0.0_WP

!aa67        do k=1,n
!aa67           do nz=1,nsigma-1
!aa67              tt = Tr_distr2(nz,k)
!aa67              ss = Sr_distr2(nz,k)
!aa67              pp = Z(nz,riv_node_ob(k))
!aa67              Tr_distr2(nz,k)=theta(ss,tt,pp,pr)
!aa67           enddo
!aa67        enddo
!aa67     end if

!aa67#ifdef USE_MPI
!aa67
!aa67     call MPI_BCast(riv_node_ob(:), n, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
!aa67     call MPI_BCast(riv_elev(:), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67     call MPI_BCast(riv_vt2(:), m, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67     call MPI_BCast(riv_vt_jd(:), m+1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67     do i=1,nsigma-1
!aa67        call MPI_BCast(Tr_distr2(i,:), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67        call MPI_BCast(Sr_distr2(i,:), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67     end do
!aa67     do i=1,nsigma-1
!aa67        do j=1,n
!aa67           call MPI_BCast(Tr_distr_t2(i,j,:), m, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67           call MPI_BCast(Sr_distr_t2(i,j,:), m, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67        end do
!aa67     end do
!aa67     do j=1,n
!aa67        call MPI_BCast(riv_elev_t(i,:), m, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
!aa67     end do

!aa67     !SH The global indices in these fields must be converted to local ones!!!!!

!aa67     allocate(mapping(nod2D))
!aa67     mapping(:)=0
!aa67     do i=1,myDim_nod2D+eDim_nod2D
!aa67        ipos=myList_nod2D(i)
!aa67        mapping(ipos)=i
!aa67     end do

!aa67     !SH DAS MUSS NOCH FERTIG!!

!aa67#endif 

!aa67     index_nod2D(riv_node_ob)=3
!aa67     write(*,*) 'Rivers set via OB, nodes are ', riv_node_ob
!aa67     riv_vt_jd(1)=Time_riv_begin
!aa67     DO i=2,m+1 
!aa67        riv_vt_jd(i)=jdh*riv_vt2(i-1)+riv_vt_jd(i-1)
!aa67     enddo
!aa67     write(*,*) 'riv_vt_jd', riv_vt_jd
!aa67  endif

  if (riv) then
     open(50,file=trim(meshpath)//'riv.out', status='old')
     read(50,*) n
     read(50,*) m
     !write(*,*) 'n,m'
     allocate(riv_ind_el(n), riv_ind_eg(n), Qr(edge2D - edge2D_in),riv_vt(m),riv_vt_jd(m+1)) 
!aa67     allocate(Tr_distr(nsigma-1,edge2D - edge2D_in),Sr_distr(nsigma-1,edge2D - edge2D_in),Qr_sig(nsigma-1,edge2D - edge2D_in))
!aa67     allocate(Tr_distr_t(nsigma-1,n,m),Sr_distr_t(nsigma-1,n,m),Qr_sig_t(nsigma-1,n,m), Qr_t(n,m))
     allocate(Qr_t(n,m))
     !write(*,*) 'alloc'
!aa67     if ((Q_sigma).and.(Stratif_sigma)) then 
!aa67        do j=1,m
!aa67           read(50,*) riv_vt(j)
!aa67           do i=1,n
!aa67              read(50,*) riv_ind_el(i),Qr_t(i,j),Qr_sig_t(1:nsigma-1,i,j),Sr_distr_t(1:nsigma-1,i,j),Tr_distr_t(1:nsigma-1,i,j)
!aa67           enddo
!aa67        enddo
!aa67        close(50)	

!aa67     else
!aa67	if (Q_sigma)then 
!aa67    !write(*,*) 'Q_sigma'
!aa67           allocate(S_t(n,m),T_t(n,m))
!aa67           do j=1,m
!aa67              read(50,*) riv_vt(j)
!aa67              !write(*,*) 'riv_vt'
!aa67              do i=1,n
!aa67                 read(50,*) riv_ind_el(i),Qr_t(i,j),Qr_sig_t(1:nsigma-1,i,j),S_t(i,j),T_t(i,j)
!aa67              enddo
!aa67           enddo
!aa67           close(50)	
!aa67           !write(*,*) 'river file is read'
!aa67           do j=1,nsigma-1
!aa67              Sr_distr_t(j,:,:)=S_t
!aa67              Tr_distr_t(j,:,:)=T_t
!aa67           enddo
!aa67           !write(*,*) 'river file is read 2'
!aa67           deallocate(S_t,T_t)

!aa67	else
!aa67           if (Stratif_sigma)then
!aa67              do j=1,m
!aa67                 read(50,*) riv_vt(j)
!aa67                 do i=1,n
!aa67                    read(50,*) riv_ind_el(i),Qr_t(i,j),Sr_distr_t(1:nsigma-1,i,j),Tr_distr_t(1:nsigma-1,i,j)
!aa67                 enddo
!aa67              enddo
!aa67              close(50)	
!aa67              do j=1,nsigma-1
!aa67                 Qr_sig_t(j,:,:)=(sigma(j) - sigma(j+1))
!aa67              enddo

!aa67           else
!aa67              allocate(S_t(n,m),T_t(n,m))
              do j=1,m
                 read(50,*) riv_vt(j)
                 do i=1,n
                    read(50,*) riv_ind_el(i),Qr_t(i,j)   !aa67,S_t(i,j),T_t(i,j)
                 enddo
              enddo
              close(50)
!aa67              do j=1,nsigma-1
!aa67                 Sr_distr_t(j,:,:)=S_t
!aa67                 Tr_distr_t(j,:,:)=T_t
!aa67                 Qr_sig_t(j,:,:)=(sigma(j) - sigma(j+1))
!aa67              enddo
!aa67              deallocate(S_t,T_t)
!aa67           endif
!aa67        endif
     endif
     !write(*,*) 'river1'
     do i=1,n
        eledges=elem_edges(:,riv_ind_el(i))
        flag=0
        do j=1,4
           if (eledges(j)>edge2D_in) then
              flag=1
              riv_ind_eg(i)=eledges(j)
           endif
        enddo
        if (flag==0) then
           write(*,*) 'Stop, wrong river input: element does not contain edge at the boundary', riv_ind_el(i)
           stop
        endif
     enddo

     Qr=0.0_WP
     Qr(riv_ind_eg-edge2D_in)=Qr_t(:,1)
!aa67     Qr_sig(:,riv_ind_eg-edge2D_in)=Qr_sig_t(:,:,1)
!aa67     Tr_distr(:,riv_ind_eg-edge2D_in)=Tr_distr_t(:,:,1)
!aa67     Sr_distr(:,riv_ind_eg-edge2D_in)=Sr_distr_t(:,:,1)

     allocate(ssin(edge2D-edge2D_in), scos(edge2D-edge2D_in), riv_w(2,edge2D-edge2D_in))
     riv_w=0.0_WP
     ssin=0.0_WP
     scos=0.0_WP
     DO i=1,n
        ed=edge_nodes(:,riv_ind_eg(i))
        elem=edge_tri(1,riv_ind_eg(i))
        n1=coord_nod2D(:,ed(1))
        n2=coord_nod2D(:,ed(2))
        call elem_center(elem, xc(1), xc(2), amin)
        A1=n2(1)-n1(1)
        B1=n2(2)-n1(2)
        C1=A1*xc(1)+B1*xc(2)
        C2=A1*n1(2)-B1*n1(1)
        aux=(C1*A1-B1*C2)/(A1**2+B1**2)
        auy=(B1*C1+A1*C2)/(A1**2+B1**2)
        !aux, auy - coordinates of the intersection point of normal from the triangle center and boundary edge
        xv=xc(1)-aux
        yv=xc(2)-auy
        ! xv, zv - vector coordinates, directed to the center of triangle from the intersection point    
        scos(riv_ind_eg(i)-edge2D_in)=xv/sqrt(xv**2+yv**2)
        ssin(riv_ind_eg(i)-edge2D_in)=yv/sqrt(xv**2+yv**2)
        !calculate weights: normal from the triangle center divides edge into two,
        !not necessarily equivalent pieces. The weights are ratio between these pieces and length of edge
        riv_w(1,riv_ind_eg(i)-edge2D_in)=sqrt((n1(1)-aux)**2+(n1(2)-auy)**2)/sqrt(A1**2+B1**2)
        riv_w(2,riv_ind_eg(i)-edge2D_in)=sqrt((n2(1)-aux)**2+(n2(2)-auy)**2)/sqrt(A1**2+B1**2)
        if ((riv_w(1,riv_ind_eg(i)-edge2D_in)>1).or.(riv_w(2,riv_ind_eg(i)-edge2D_in)>1)) then
           write(*,*) 'Normal from the center of element to the boundary edge line does not intersect boundary edge, ', &
                      i, riv_ind_el(i)
           riv_w(1,riv_ind_eg(i)-edge2D_in)=0.5_WP
           riv_w(2,riv_ind_eg(i)-edge2D_in)=0.5_WP
        endif
     enddo
     write(*,*) 'river4'
     riv_amount=n
     allocate(riv_vel(2,n))
!aa67     allocate(riv_vel_u(1:nsigma-1,n),riv_vel_v(1:nsigma-1,n))
     call calc_disch_node_area(1)
     !write(*,*) 'river9'
     riv_vt_jd(1)=Time_riv_begin
     DO i=2,m+1 
        riv_vt_jd(i)=jdh*riv_vt(i-1)+riv_vt_jd(i-1)
     enddo
     write(*,*) 'riv_vt_jd', riv_vt_jd
!aa67  endif

end subroutine initial_riv


subroutine calc_disch_node_area(rk)
!VF,calc_disch_node_area
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  use o_UTILIT

  implicit none
  type calc_ir
  integer, allocatable, dimension(:)  :: ci
  real(kind=WP),allocatable, dimension(:)  :: cr
  end type calc_ir
 
  integer       :: i,k, indc, ind0, irn, irn2, iln, inds(riv_amount), nz,n, rk
  real(kind=WP) :: ai(2), ai2(2), tt, ff, pp, pr, ss
  TYPE(calc_ir) :: Calc
  real(kind=WP), external :: theta
  
! The task of this subroutine is to distribute the total discharge on edges among nodes,
! additionally, distribute salinity and temperature fluxes. This subroutine is working based on
! weights riv_w, which were calculated in initial_riv.

  riv_num_nodes=riv_amount
  allocate(Calc%ci(2*riv_amount))
  allocate(Calc%cr(2*riv_amount))
  Calc%cr=0.0_WP
  Calc%ci=0
  Calc%cr(1)=1
  indc=2
  ind0=1
  k=1
  iln=0
  irn=0
  inds=riv_ind_eg
  !write(*,*) 'river6'
  do while(any(inds>0))

     do while (inds(k)==0)
        k=k+1
        if (k>size(inds)) k=2
     enddo

     call  find_left_ne(riv_ind_eg, inds(k),iln,riv_amount)
     call  find_right_ne(riv_ind_eg, inds(k),irn,riv_amount)

     if ((iln==0).and.(irn==0)) then
        Calc%ci(indc)=inds(k)
        Calc%cr(indc)=Qr(inds(k)-edge2D_in)
        indc=indc+1
        Calc%cr(indc-ind0-1)=ind0
        Calc%cr(indc)=0.0_WP
        indc=indc+1
        inds(k)=0
        riv_num_nodes=riv_num_nodes+1
     else
        if (iln==0) then
           Calc%ci(indc)=inds(k)
           Calc%cr(indc)=Qr(inds(k)-edge2D_in)
           indc=indc+1
           riv_num_nodes=riv_num_nodes+1
           inds(k)=0
           do while (irn>0) 
              ind0=ind0+1
              Calc%ci(indc)=inds(irn)
              Calc%cr(indc)=Qr(inds(irn)-edge2D_in)
              indc=indc+1
              call  find_right_ne(riv_ind_eg, inds(irn),irn2,riv_amount)
              inds(irn)=0
              irn=irn2
           enddo
           Calc%cr(indc-ind0-1)=ind0
           ind0=1
           Calc%cr(indc)=0.0_WP
           indc=indc+1
        else
           k=k+1
           if (k>size(inds)) k=2
        endif
     endif
  enddo
  if (rk<2) then
!aa67        allocate (Tr_node_sig(nsigma-1,riv_num_nodes),Sr_node_sig(nsigma-1,riv_num_nodes), riv_node(riv_num_nodes))
!aa67       allocate(Qr_node(riv_num_nodes),Qr_node_sig(nsigma-1,riv_num_nodes))
     allocate (riv_node(riv_num_nodes))
     allocate(Qr_node(riv_num_nodes))
  endif
  k=2
  ind0=1

  do while (k<=riv_num_nodes)
     riv_node(ind0)=edge_nodes(1,Calc%ci(k))
     !write(*,*) 'ed_nod1', edge_nodes(1,int(Calc%ci(k))),k,ind0, riv_node(ind0)
     Qr_node(ind0)=Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)
!aa67       Qr_node_sig(:,ind0)=Qr_node(ind0)*Qr_sig(:,Calc%ci(k)-edge2D_in)
!aa67       Tr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Tr_distr(:,Calc%ci(k)-edge2D_in)
!aa67       Sr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Sr_distr(:,Calc%ci(k)-edge2D_in)
     ind0=ind0+1
     k=k+1
     do while (Calc%ci(k)>0.0_WP)
        riv_node(ind0)=edge_nodes(1,int(Calc%ci(k)))
        !write(*,*) 'ed_nod2', edge_nodes(1,int(Calc%ci(k))),k, ind0, riv_node(ind0)
        Qr_node(ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)+Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)
!aa67          Qr_node_sig(:,ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)*Qr_sig(:,Calc%ci(k-1)-edge2D_in)&
!aa67          +Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)*Qr_sig(:,Calc%ci(k)-edge2D_in)
!aa67          Tr_node_sig(:,ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)*Qr_sig(:,Calc%ci(k-1)-edge2D_in)&
!aa67         *Tr_distr(:,Calc%ci(k-1)-edge2D_in)+Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)&
!aa67         *Qr_sig(:,Calc%ci(k)-edge2D_in)*Tr_distr(:,Calc%ci(k)-edge2D_in)
!aa67          Sr_node_sig(:,ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)&
!aa67          *Qr_sig(:,Calc%ci(k-1)-edge2D_in)*Sr_distr(:,Calc%ci(k-1)-edge2D_in) &
!aa67          +Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)*Qr_sig(:,Calc%ci(k)-edge2D_in)*Sr_distr(:,Calc%ci(k)-edge2D_in)
        ind0=ind0+1
        k=k+1
     enddo
     riv_node(ind0)=edge_nodes(2,int(Calc%ci(k-1)))
     !write(*,*) 'ed_nod3', edge_nodes(2,int(Calc%ci(k-1))),k, ind0, riv_node(ind0)
     Qr_node(ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)
!aa67       Qr_node_sig(:,ind0)=Qr_node(ind0)*Qr_sig(:,Calc%ci(k-1)-edge2D_in)
!aa67       Tr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Tr_distr(:,Calc%ci(k-1)-edge2D_in)
!aa67       Sr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Sr_distr(:,Calc%ci(k-1)-edge2D_in)
     ind0=ind0+1
     k=k+1
  enddo
  write(*,*)'Rivers set via SB, elements are', riv_ind_eg
  !write(*,*) 'Qr_node', riv_num_nodes,Qr_node
  write(*,*) 'riv_node',riv_node
  ! write(*,*) 'discharge', sum (Qr_sig(:,riv_ind_eg(1)-edge2D_in)),sum (Qr_sig(:,riv_ind_eg(2)-edge2D_in)),&
  !sum (Qr_sig(:,riv_ind_eg(3)-edge2D_in)),sum (Qr_sig(:,riv_ind_eg(4)-edge2D_in)),&
  !sum (Qr_sig(:,riv_ind_eg(5)-edge2D_in)),sum (Qr_sig(:,riv_ind_eg(6)-edge2D_in)),Qr_node_sig
!aa67   write(*,*) 'discharge',sum(Qr_node_sig(:,:)),sum(Qr_node), Qr_node !, Qr(riv_ind_eg-edge2D_in)
  !write(*,*) 'river8', Sr_distr(:,riv_ind_eg-edge2D_in)

  !  write(*,*) 'Qr_sig', Qr_sig
  deallocate(Calc%ci,Calc%cr)
  pr=0.0_WP
!aa67       do n=1,riv_num_nodes
!aa67      do nz=1,nsigma-1
!aa67       tt = Tr_node_sig(nz,n)/Qr_node_sig(nz,n)
!aa67       ss = Sr_node_sig(nz,n)/Qr_node_sig(nz,n)
!aa67       pp = Z(nz,riv_node(n))
!aa67       Tr_node_sig(nz,n)=theta(ss,tt,pp,pr)*Qr_node_sig(nz,n)
!aa67      enddo
!aa67     enddo

end subroutine calc_disch_node_area

subroutine initial_state_wind

  use o_MESH
  use o_ARRAYS
  use o_PARAM
  !
  implicit none

  integer             :: el, q, elnodes(4), wind, mass
  real(kind=WP) :: y1, y, dst
  
  
   
  wind=0
  mass=0 
  
  ! =========
  ! initialize with zeros
  ! =========
  eta_n=0.0_WP
  U_n_2D=0.0_WP
  if(Tstepping==1) then
  eta_n_1=0.0_WP
  eta_n_1=0.0_WP
  U_n_1=0.0_WP
  U_n_2=0.0_WP
  end if  
  
  
  IF (wind==1) THEN
  ! Stress divided by density_0
  
  ! initialize wind
  y=minval(coord_nod2D(2,:))
  y1=maxval(coord_nod2D(2,:))
  DO el=1,elem2D
  elnodes=elem2D_nodes(:,el)
  q=4
  if(elnodes(1)==elnodes(4)) q=3
     taux(el)=0.0001_WP*cos(pi*(sum(coord_nod2d(2,elnodes(1:q)))/dble(q)-y)/(y1-y))
     tauy(el)=0.0_WP
  END DO
    ! Southern ocean forcing
  !DO el=1, elem2D
  !   elnodes=elem2D_nodes(:, el)
  !   q=4
  ! if(elnodes(1)==elnodes(4)) q=3
  !   y=sum(coord_nod2D(2,elnodes(1:q)))/dble(q)
  !   if(y<-32.*rad) then
  !  dst=1.0
  !   if(y<-64.*rad) dst=(y/rad+85)/21.0
  !   taux(el)=0.2*dst*dst*sin(pi*abs(y+32.0*rad)/(32.0*rad))
  !   end if
  !  tauy(el)=0.
  !END DO
  END IF
  if(mass==1) THEN
  ! provide a pattern with relaxation to prescribed elevation
  DO el=1,nod2D
  y=coord_nod2D(2,el)
  y1=coord_nod2D(1,el)
  dst=((y-50.0_WP*rad)**2+y1**2)/rad**2/4.0_WP   ! radius of 2 degrees
  if(dst<1) then  
  relax_coef(el)=cos(pi*sqrt(dst)/2.0_WP)/3600.0_WP/24.0_WP/3.0_WP
  end if
  ! Set signal frequency
  !O_ext=2*pi/24.0_WP/3600.0_WP/365.0_WP/2.0_WP
  
  ENDDO
  END IF
  
end subroutine initial_state_wind

!===========================================================================

subroutine comp_sediment_ini

  !VF,initial_riv
  use o_MESH
  use o_ARRAYS
  use o_PARAM

  use g_parsup

  implicit none
  integer       :: n,m,j,i
  real(kind=WP) :: b,s,dst, teta,w_s1,w_s2,w_s3,c

  integer  :: ierror

  CF=0.0_WP
  c_old=CF
  Cclim = CF

  if (sed_boun_flux) then
     if (mype==0) then

        open(50,file=trim(meshpath)//trim(TITLE)//'_sed_ob.out', status='old')
        read(50,*) n
        read(50,*) m

        close(50)
     end if

#ifdef USE_MPI
     call MPI_BCast(n, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     call MPI_BCast(m, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
#endif

     allocate(cr_distr(nsigma-1,nobn,m),cr_distr2(nsigma-1,nobn),sed_vt2(m)) 

     if (mype==0) then
        open(50,file=trim(meshpath)//trim(TITLE)//'_sed_ob.out', status='old')
        read(50,*) n
        read(50,*) m

        do j=1,m
           read(50,*) sed_vt2(j)
           do i=1,n
              read(50,*)cr_distr(1:nsigma-1,i,j)
           enddo
        enddo
        close(50)

        cr_distr2(:,:)=cr_distr(:,:,1)
     end if

#ifdef USE_MPI     
     call MPI_BCast(sed_vt2(:), m, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     do i=1,nsigma-1
        call MPI_BCast(cr_distr2(i,1:nobn), nobn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     end do
     do i=1,nsigma-1
        do j=1,nobn
           call MPI_BCast(cr_distr(i,j,1:m), nobn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        end do
     end do
#endif

  end if


  !************************************************************************
  ! average settling velocity of the material available for transport (w_s)
  ! Simpson and Castelltort 2006, Coupled model of surface water flow, 
  ! sediment transport and morphological evolution. Computers & Geosciences,
  ! 32, pp. 1600-1614.
  ! ************************************************************************
  b = 1.0_WP/3.0_WP
  s = plop_s/density_0 - 1.0_WP
  dst = d_m*(g*s/snu_kin**2)**b


  ! Criteria is from Miller, Cave and Komar, 1977, Threshold of sediment motion under unidirectional current + Shields criterium (dst>=4)

  if (dst<4) then
     teta=0.115_WP*(dst)**(-0.5_WP)
  else
     if (dst<=10) then
        teta=0.14_WP*(dst)**(-0.64_WP)
     else
        if (dst<=20) then
           teta=0.04_WP*(dst)**(-0.1_WP)
        else
           if (dst<=150) then
              teta=0.013_WP*(dst)**(0.29_WP)
           else
              teta= 0.055_WP
           endif
        endif
     endif
  endif

  c = 13.95_WP*snu_kin/d_m
  w_s1 = sqrt(c*c + 1.09_WP*d_m*g*s) - c
  w_s2 =snu_kin/d_m*((10.36_WP**2+1.049*teta**3)**(0.5)-1-.36_WP)
  w_s3=g*(plop_s-density_0)*d_m**2/(18.0_WP*snu_kin*density_0)

  if (mype==0) write(*,*) 'w_s1, w_s2, w_s3,...', w_s1,w_s2,w_s3

end subroutine comp_sediment_ini

