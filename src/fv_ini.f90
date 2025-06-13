! ====================================================================
subroutine read_restart
!VF
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
!VF revised 27.11.18
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  use fv_ic

  !
  implicit none
  integer       :: i,el, elnodes(4), n,m,j,jjj, nl, nz, c_riv, eledges(4), flag, elem
  real(kind=WP) :: x1, x2, y1, y2, amp, a, den_0
  real(kind=WP) :: dx, dy, period, hu, x, y
  real(kind=WP) :: pr, tt, ss, pp, d
  real(kind=WP) :: de, ga, ep, fe
  real(kind=WP), external :: theta
  real(kind=WP), allocatable, dimension(:,:)  ::   aux
  real, allocatable, dimension(:,:)  ::  S_t, T_t
! AA for compute baroclinic pressure for 2D barotropic task
  real, allocatable, dimension(:)    :: tmp_plop

  real(kind=WP) :: dmean
	
!====================================
! Calculating of reference density	
!=====================================

scale_area=sum(elem_area)/elem2D
write(*,*) 'scale_area = ',scale_area

write(*,*) '**************************************************'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) '**************************************************'
write(*,*) 'FOR ATMOSPHERIC USED ONLY WIND AND PRESSURE'
write(*,*) 'SEE fv_sbc.f90 lines: 856 and 1006'
write(*,*) 'VERSION WITH SOME OFFSET IN SBC FOR WIND AND ATMOSPHERIC PRESSURE'
write(*,*) 'SEE fv_sbc.f90 line: 884'
write(*,*) '**************************************************'
write(*,*) '**************************************************'

 
!====================================
! Calculating coeff. to determine proper eta (see fv_average_dynamic) 
! in case of absence or presence tidal potential.
! Attention, in case of presence tidal potential, a_tp can be also calculating 
! for each time step to take into account the influence of internal waves presence 
! to tidal potential energy	
!=====================================
allocate(a_tp(nod2D,2))
 if (T_potential) then
 a_tp(:,1)=0.918_WP
 a_tp(:,2)=0.69_WP
 else
 a_tp(:,1)=1.0_WP
 a_tp(:,2)=0.0_WP
 endif
 do i=1,nod2D
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


 do n=1,nod2d
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
    
write(*,*) 'MAX_MIN_SSH_initial= ', maxval(eta_n),minval(eta_n)
write(*,*) 'MAX_MIN_depth= ', maxval(depth),minval(depth)

if (comp_sediment) write(*,*) 'comp_sed'
write(*,*) 'hu hu', max(Dmin,eta_n(n) + depth(n))
write(*,*) 'hu hu', con(1),con_bc(2)

! AA initialzation for sediment for 2d case
        do n=1,nod2D
          dmean=max(Dmin,eta_n(n) + depth(n))
          con(n) = 0.0_WP/dmean
	  con_bc(n) = con(n)
        enddo
! AA 
!VF initialization of sediment module and taking into account
!presence of sediments into the density calculation - 3d case

!if (comp_sediment) call comp_sediment_ini

density_0=1023.66_WP

!VF, calculate density using info about T_const, S_const and Concentration of sus. sediments
! if (comp_sediment) then
! call densityJM(T_const, S_const, 0.0, den_0,(maxval(CF)+minval(CF))*0.5_WP)
! else
! call densityJM(T_const, S_const, 0.0, den_0)
! endif

!write(*,*) 'den_0 = ', den_0
! density_0=density_0+den_0

write(*,*) 'density_0', density_0


if (TF_presence) then  
  open(50,file=trim(meshpath)//trim(TITLE)//'_m2.out', status='old')
  read(50,*) n
  
  allocate(ampt(n,15),fazt(n,15))
  ampt=0.0_WP
  fazt=0.0_WP
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
   enddo
   write (*,*) 'max_amplitude', maxval(ampt)
   deallocate(aux)
endif   

!
! Irradiance in the upper ocean

!
if (type_task>2) then  

if (.not.(jer_const)) then
allocate(aj(nod2D),bj(nod2D),Rj(nod2D))

open(50,file=trim(meshpath)//trim(TITLE)//'_jer.out', status='old')

   do j=1,nod2D
   read(50,*) Rj(j),aj(j),bj(j)
    enddo
 close(50)	
endif

  if (key_ic) then
    call ic_do  ! from fv_ic.f90
 else
     if (TEMP_ON) then
      open(50,file=trim(meshpath)//trim(TITLE)//'_temp.out', status='old')
     do n=1,nod2D
       read(50,*) (TF(nz,n),nz=1,nsigma-1)
      enddo
      close(50)
      else
       TF = T_const
       TF2 = T_const
     endif
     if (SALINITY_ON) then
      open(50,file=trim(meshpath)//trim(TITLE)//'_sal.out', status='old')
     do n=1,nod2D
      read(50,*) (SF(nz,n),nz=1,nsigma-1)
      enddo
      close(50)
      else
       SF = S_const
       SF2 = S_const
       TF = T_const
       TF2 = T_const
     endif
  endif
!
!   convert in situ temperature into potential temperature
!
   pr=0.0_WP
   do n=1,nod2D
    do nz=1,nsigma-1
     tt = TF(nz,n)
     ss = SF(nz,n)
     pp = Z(nz,n)
     TF(nz,n)=theta(ss,tt,pp,pr)
    enddo
if (coord_nod2D(2,n)<0) then
TF(:,n)=15.0_WP
SF(:,n)=0.0_WP
endif
   enddo

      T_old= TF
      S_old= SF
      T_old2= TF
      S_old2= SF
      Tclim = TF
      Sclim = SF
   write(*,*)'T_max, T_min', maxval(TF), minval(TF)
   write(*,*)'S_max, S_min', maxval(SF), minval(SF)      
endif

!VF: Initialisation of C_d depends on user choice
if (BFC==2) then
DO elem=1,elem2D
   elnodes=elem2D_nodes(:,elem)
    dmean = max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))
    C_d_el(elem)=(cka/(log(dmean/(z0b_min+za))-1))**2
enddo
else
    C_d_el=C_d
endif
if ((tracer_adv==3).or. (tracer_adv==5)) then
call fct_init
endif

  !   open(50,file=trim(meshpath)//trim(TITLE)//'_Bar2D.out', status='old')
  !   do elem=1,elem2D
  !    read(50,*) Bar_pr_2D(1,elem), Bar_pr_2D(2,elem)
  !   enddo

  !  write(*,*) 'Maxmin baroclinic pressure gradient for 2D: ',maxval(Bar_pr_2D(1,:)),minval(Bar_pr_2D(1,:))
  !  write(*,*) 'Maxmin baroclinic pressure gradient for 2D: ',maxval(Bar_pr_2D(2,:)),minval(Bar_pr_2D(2,:))


   write(*,*)'finish initial_state'

end subroutine initial_state
! ====================================================================
subroutine initial_riv
!VF,initial_riv
  use o_MESH
  use o_ARRAYS
  use o_PARAM

  implicit none
  integer       :: i,el, n,m,j, nz, eledges(4), flag, k, ed(2), elem
  real(kind=WP) :: xv, yv, aux, auy, xc(2), n1(2), n2(2), A1, A2, B1, B2, C1, C2, amin
  real(kind=WP) :: ss,tt,pp, pr
  real(kind=WP), external :: theta
  real(kind=WP), allocatable, dimension(:,:)  ::  S_t, T_t

if (riv_OB) then
   open(50,file=trim(meshpath)//trim(TITLE)//'_riv_ob.out', status='old')
   read(50,*) n
   read(50,*) m

   allocate(riv_node_ob(n), riv_elev(n),riv_vt2(m),riv_vt_jd(m+1)) 
   allocate(Tr_distr2(nsigma-1,n),Sr_distr2(nsigma-1,n))
   allocate(Tr_distr_t2(nsigma-1,n,m),Sr_distr_t2(nsigma-1,n,m), riv_elev_t(n,m))

   do j=1,m
   read(50,*) riv_vt2(j)
   do i=1,n
    read(50,*) riv_node_ob(i),riv_elev_t(i,j),Sr_distr_t2(1:nsigma-1,i,j),Tr_distr_t2(1:nsigma-1,i,j)
    enddo
    enddo
    close(50)

   Tr_distr2(:,:)=Tr_distr_t2(:,:,1)
   Sr_distr2(:,:)=Sr_distr_t2(:,:,1) 
   riv_elev = riv_elev_t(:,1)
   riv_amount_ob=n
   pr=0.0_WP

   do k=1,n
    do nz=1,nsigma-1
     tt = Tr_distr2(nz,k)
     ss = Sr_distr2(nz,k)
     pp = Z(nz,riv_node_ob(k))
     Tr_distr2(nz,k)=theta(ss,tt,pp,pr)
    enddo
   enddo

index_nod2D(riv_node_ob)=3
write(*,*) 'Rivers set via OB, nodes are ', riv_node_ob
riv_vt_jd(1)=Time_riv_begin
DO i=2,m+1 
riv_vt_jd(i)=jdh*riv_vt2(i-1)+riv_vt_jd(i-1)
enddo
write(*,*) 'riv_vt_jd', riv_vt_jd
endif		
if (riv) then
   open(50,file=trim(meshpath)//trim(TITLE)//'_riv.out', status='old')
   read(50,*) n
   read(50,*) m
!write(*,*) 'n,m'
   allocate(riv_ind_el(n), riv_ind_eg(n), Qr(edge2D - edge2D_in),riv_vt(m),riv_vt_jd(m+1)) 
   allocate(Tr_distr(nsigma-1,edge2D - edge2D_in),Sr_distr(nsigma-1,edge2D - edge2D_in),Qr_sig(nsigma-1,edge2D - edge2D_in))
   allocate(Tr_distr_t(nsigma-1,n,m),Sr_distr_t(nsigma-1,n,m),Qr_sig_t(nsigma-1,n,m), Qr_t(n,m))
!write(*,*) 'alloc'
   if ((Q_sigma).and.(Stratif_sigma)) then 
   do j=1,m
   read(50,*) riv_vt(j)
   do i=1,n
   read(50,*) riv_ind_el(i),Qr_t(i,j),Qr_sig_t(1:nsigma-1,i,j),Sr_distr_t(1:nsigma-1,i,j),Tr_distr_t(1:nsigma-1,i,j)
   enddo
   enddo
   close(50)	

	else
	if (Q_sigma)then 
!write(*,*) 'Q_sigma'
	allocate(S_t(n,m),T_t(n,m))
   do j=1,m
   read(50,*) riv_vt(j)
!write(*,*) 'riv_vt'
   do i=1,n
    read(50,*) riv_ind_el(i),Qr_t(i,j),Qr_sig_t(1:nsigma-1,i,j),S_t(i,j),T_t(i,j)
   enddo
    enddo
    close(50)	
!write(*,*) 'river file is read'
   do j=1,nsigma-1
	Sr_distr_t(j,:,:)=S_t
	Tr_distr_t(j,:,:)=T_t
   enddo
!write(*,*) 'river file is read 2'
	deallocate(S_t,T_t)
	
	else
    if (Stratif_sigma)then
   do j=1,m
   read(50,*) riv_vt(j)
   do i=1,n
    read(50,*) riv_ind_el(i),Qr_t(i,j),Sr_distr_t(1:nsigma-1,i,j),Tr_distr_t(1:nsigma-1,i,j)
    enddo
    enddo
    close(50)	
	 do j=1,nsigma-1
	Qr_sig_t(j,:,:)=(sigma(j) - sigma(j+1))
         enddo
	
	else
	allocate(S_t(n,m),T_t(n,m))
   do j=1,m
   read(50,*) riv_vt(j)
   do i=1,n
   read(50,*) riv_ind_el(i),Qr_t(i,j),S_t(i,j),T_t(i,j)
   enddo
   enddo
   close(50)
 do j=1,nsigma-1
	Sr_distr_t(j,:,:)=S_t
	Tr_distr_t(j,:,:)=T_t
	Qr_sig_t(j,:,:)=(sigma(j) - sigma(j+1))
enddo
	deallocate(S_t,T_t)
   endif
   endif
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
 	Qr_sig(:,riv_ind_eg-edge2D_in)=Qr_sig_t(:,:,1)
	Tr_distr(:,riv_ind_eg-edge2D_in)=Tr_distr_t(:,:,1)
	Sr_distr(:,riv_ind_eg-edge2D_in)=Sr_distr_t(:,:,1) 

   allocate(ssin(n), scos(n), riv_w(2,edge2D-edge2D_in))
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
   xv=(xc(1)-aux)*elem_cos(elem)
   yv=xc(2)-auy
! xv, zv - vector coordinates, directed to the center of triangle from the intersection point    
   scos(i)=xv/sqrt(xv**2+yv**2)
   ssin(i)=yv/sqrt(xv**2+yv**2)
!calculate weights: normal from the triangle center divides edge into two,
!not necessarily equivalent pieces. The weights are ratio between these pieces and length of edge
 !  riv_w(1,riv_ind_eg(i)-edge2D_in)=sqrt((n1(1)-aux)**2+(n1(2)-auy)**2)/sqrt(A1**2+B1**2)
 !  riv_w(2,riv_ind_eg(i)-edge2D_in)=sqrt((n2(1)-aux)**2+(n2(2)-auy)**2)/sqrt(A1**2+B1**2)
 !  if ((riv_w(1,riv_ind_eg(i)-edge2D_in)>1).or.(riv_w(2,riv_ind_eg(i)-edge2D_in)>1)) then
 !  write(*,*) 'Normal from the center of element to the boundary edge line does not intersect boundary edge, ', i, riv_ind_el(i)
 !  riv_w(1,riv_ind_eg(i)-edge2D_in)=0.5_WP
 !  riv_w(2,riv_ind_eg(i)-edge2D_in)=0.5_WP
 ! endif
 
! Attention: at the moment the node control volume is based on medians, not  heights, 
! therefore there is a simplification with the weights

 riv_w(1,riv_ind_eg(i)-edge2D_in)=0.5_WP
 riv_w(2,riv_ind_eg(i)-edge2D_in)=0.5_WP
   enddo
!write(*,*) 'scos, ssin', scos, ssin, xc(1), aux, xc(2), auy
   riv_amount=n
    allocate(riv_vel(2,n))
    allocate(riv_vel_u(1:nsigma-1,n),riv_vel_v(1:nsigma-1,n))
   call calc_disch_node_area(1)
!write(*,*) 'river9'
riv_vt_jd(1)=Time_riv_begin
DO i=2,m+1 
riv_vt_jd(i)=jdh*riv_vt(i-1)+riv_vt_jd(i-1)
enddo
write(*,*) 'riv_vt_jd', riv_vt_jd
endif

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
   allocate (Tr_node_sig(nsigma-1,riv_num_nodes),Sr_node_sig(nsigma-1,riv_num_nodes), riv_node(riv_num_nodes))
  allocate(Qr_node(riv_num_nodes),Qr_node_sig(nsigma-1,riv_num_nodes))
endif
  k=2
  ind0=1

  do while (k<=riv_num_nodes)
  riv_node(ind0)=edge_nodes(1,Calc%ci(k))
!write(*,*) 'ed_nod1', edge_nodes(1,int(Calc%ci(k))),k,ind0, riv_node(ind0)
  Qr_node(ind0)=Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)
  Qr_node_sig(:,ind0)=Qr_node(ind0)*Qr_sig(:,Calc%ci(k)-edge2D_in)
  Tr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Tr_distr(:,Calc%ci(k)-edge2D_in)
  Sr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Sr_distr(:,Calc%ci(k)-edge2D_in)
  ind0=ind0+1
  k=k+1
  do while (Calc%ci(k)>0.0_WP)
  riv_node(ind0)=edge_nodes(1,int(Calc%ci(k)))
!write(*,*) 'ed_nod2', edge_nodes(1,int(Calc%ci(k))),k, ind0, riv_node(ind0)
  Qr_node(ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)+Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)
  Qr_node_sig(:,ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)*Qr_sig(:,Calc%ci(k-1)-edge2D_in)&
  +Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)*Qr_sig(:,Calc%ci(k)-edge2D_in)
  Tr_node_sig(:,ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)*Qr_sig(:,Calc%ci(k-1)-edge2D_in)&
 *Tr_distr(:,Calc%ci(k-1)-edge2D_in)+Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)&
 *Qr_sig(:,Calc%ci(k)-edge2D_in)*Tr_distr(:,Calc%ci(k)-edge2D_in)
  Sr_node_sig(:,ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)&
  *Qr_sig(:,Calc%ci(k-1)-edge2D_in)*Sr_distr(:,Calc%ci(k-1)-edge2D_in) &
  +Calc%cr(k)*riv_w(1,Calc%ci(k)-edge2D_in)*Qr_sig(:,Calc%ci(k)-edge2D_in)*Sr_distr(:,Calc%ci(k)-edge2D_in)
  ind0=ind0+1
  k=k+1
  enddo
  riv_node(ind0)=edge_nodes(2,int(Calc%ci(k-1)))
!write(*,*) 'ed_nod3', edge_nodes(2,int(Calc%ci(k-1))),k, ind0, riv_node(ind0)
  Qr_node(ind0)=Calc%cr(k-1)*riv_w(2,Calc%ci(k-1)-edge2D_in)
  Qr_node_sig(:,ind0)=Qr_node(ind0)*Qr_sig(:,Calc%ci(k-1)-edge2D_in)
  Tr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Tr_distr(:,Calc%ci(k-1)-edge2D_in)
  Sr_node_sig(:,ind0)=Qr_node_sig(:,ind0)*Sr_distr(:,Calc%ci(k-1)-edge2D_in)
  ind0=ind0+1
  k=k+1
  enddo
  write(*,*)'Rivers set via SB, elements are', riv_ind_eg
  write(*,*) 'riv_node, Qr_node',riv_node, Qr_node
  write(*,*) 'discharge sum checking', sum(Qr_node_sig(:,:)),sum(Qr_node)
!  write(*,*) 'Qr_sig', Qr_sig

  deallocate(Calc%ci,Calc%cr)
pr=0.0_WP
     do n=1,riv_num_nodes
    do nz=1,nsigma-1
     tt = Tr_node_sig(nz,n)/Qr_node_sig(nz,n)
     ss = Sr_node_sig(nz,n)/Qr_node_sig(nz,n)
     pp = Z(nz,riv_node(n))
     Tr_node_sig(nz,n)=theta(ss,tt,pp,pr)*Qr_node_sig(nz,n)
    enddo
   enddo
end subroutine calc_disch_node_area



subroutine initial_state_wind
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  !
  implicit none
  integer             :: el, q, elnodes(4), wind, mass
  real(kind=8) :: y1, y, dst
  
  
   
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

  implicit none
  integer       :: n,m,j,i
  real(kind=WP) :: b,s,dst, teta,w_s1,w_s2,w_s3,c

CF=0.0_WP
c_old=CF
Cclim = CF
if (sed_boun_flux) then
   open(50,file=trim(meshpath)//trim(TITLE)//'_sed_ob.out', status='old')
   read(50,*) n
   read(50,*) m

   allocate(cr_distr(nsigma-1,nobn,m),cr_distr2(nsigma-1,nobn),sed_vt2(m)) 

   do j=1,m
   read(50,*) sed_vt2(j)
   do i=1,n
    read(50,*)cr_distr(1:nsigma-1,i,j)
    enddo
    enddo
    close(50)

   cr_distr2(:,:)=cr_distr(:,:,1)
endif
       

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
write(*,*) 'w_s1, w_s2, w_s3,...', w_s1,w_s2,w_s3

end subroutine comp_sediment_ini

subroutine comp_sediment_ini_2D
!VF,initial_sed
  use o_MESH
  use o_ARRAYS
  use o_PARAM

  implicit none
  integer       :: n,m,j,i
  real(kind=WP) :: b,s,dst, teta,w_s1,w_s2,w_s3,c


if (sed_boun_flux) then
   open(50,file=trim(meshpath)//trim(TITLE)//'_sed_ob_2d.out', status='old')
   read(50,*) n
   read(50,*) m

   allocate(cr_distr_2D(n,m),cr_distr2_2D(nod2d),sed_vt2(m), sed_ind(n), Qr_sed(nod2d)) 
   cr_distr2_2D=0.0_WP
   Qr_sed=0.0_WP

   do j=1,m
   read(50,*) sed_vt2(j)
   do i=1,n
    read(50,*)   sed_ind(i), cr_distr_2D(i,j)
    enddo
    enddo
    close(50)

   cr_distr2_2D(sed_ind)=cr_distr_2D(:,1)
   Qr_sed(sed_ind)=1.0_WP
   
   write(*,*) 'SPM', cr_distr2_2D (sed_ind)
endif

end subroutine comp_sediment_ini_2D
