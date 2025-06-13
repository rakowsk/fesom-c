! Shallow water code based cell-vertex finite volume discretization
! AB3-AM4 time stepping as suggested in Shchepetkin and McWilliams
! (2005)
! Serial version
! July 2012
! sergey.danilov@awi.de



!===========================================================================
subroutine read_mesh2D
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE o_UTILIT
!
IMPLICIT NONE
!
INTEGER               :: nq, nt, n_quad, n_tri
INTEGER               :: n1,n2,n3, n4, nod(4)
INTEGER               :: n,ind
REAL(kind=WP)         :: x1, x2
INTEGER, allocatable  :: elem_data(:)
INTEGER               :: i_error

! Requires new mesh format: elem2D lists 4 nodes in each row.
! For triangles the forth node coincides with the first one or just list the numbers of 3 nodes

  open (20,file=trim(meshpath)//trim(TITLE)//'_nod2d.out', status='old')
  open (21,file=trim(meshpath)//trim(TITLE)//'_elem2d.out', status='old')
  open (22,file=trim(meshpath)//trim(TITLE)//'_depth.out', status='old')

READ(20,*) nod2D
  ALLOCATE(coord_nod2D(2,nod2D),index_nod2D(nod2D))
  nobn=0

  if (cartesian) then
    do n=1,nod2D
     read(20,*) nq, x1, x2, ind
     index_nod2D(nq) = ind
	 coord_nod2D(1,nq)=x1/r_earth
     coord_nod2D(2,nq)=x2/r_earth
	 if (ind==2) then
	 nobn=nobn+1
	 endif
  end do
  else
     do n=1,nod2D
     read(20,*) nq, x1, x2, ind
     index_nod2D(nq) = ind
     coord_nod2D(1,nq)=x1*rad
     coord_nod2D(2,nq)=x2*rad
	 if (ind==2) then
	 nobn=nobn+1
	 endif
  end do
  endif

  close(20)
  !VF, create in_obn array with OBN indexes
allocate(in_obn(nobn))
ind=1
 do n=1,nod2D
	 if (index_nod2D(n)==2) then
	 in_obn(ind)=n
         ind=ind+1
	 endif
  end do

  read(21,*)  elem2D
  ALLOCATE(elem2D_nodes(4,elem2D))
  ALLOCATE(elem_data(4*elem2D))
  elem_data(:)=-1

  ! meshes with quads have 4 columns, but TsunAWI grids may be
  ! purely triangular, with 3 columns each. Test, how many
  ! columns there are!

  read(21,*,iostat=i_error) elem_data(1:4*elem2D)
write(*,*) i_error
  if (i_error == 0) then
     ! There is a fourth column => quad or mixed mesh
     n_quad = 0
     n_tri = 0
     do n=1,elem2D
        nod(1:4) = elem_data((n-1)*4+1:n*4)
        ! It might be important to have the elements sorted:
        ! triangles first.
        if (nod(1) == nod(4)) then   ! triangle
           n_tri = n_tri+1
           elem2D_nodes(1:4,n_tri) = nod(1:4)
        endif
     enddo
     do n=1,elem2D
        nod(1:4) = elem_data((n-1)*4+1:n*4)
        if (nod(1) /= nod(4)) then   ! quad
           n_quad = n_quad+1
           elem2D_nodes(1:4,n_tri + n_quad) = nod(1:4)
        endif
     enddo

  else
     ! No third column => triangles only
     n_quad=0
     n_tri=elem2D
     do n=1,elem2D
        elem2D_nodes(1:3,n) = elem_data((n-1)*3+1:n*3)
        elem2D_nodes(4,n) = elem_data((n-1)*3+1)
     enddo
  end if

  deallocate(elem_data)

  elem2D_tri = n_tri


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION: The list of elements may now differ from that in elem2d.out.
!  In order to work with the same lists put triangles first.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 allocate(depth(nod2D))
  read(22,*) depth
  CLOSE(22)
  do n=1,nod2D
    if(index_nod2D(n) == 1 ) depth(n) = -15.0_WP ! Wind Farms
 enddo

  ALLOCATE(ac(nod2D))

  ac=1.0_WP
  if (SL_obn_user_def) then
  open (23,file=trim(meshpath)//trim(TITLE)//'_ac.out', status='old')
  read(23,*) ac
  CLOSE(23)
  endif

  if ((SL_obn).and.(.not.(SL_obn_user_def))) then
     call ac_create(ac)
  endif

if (type_task>1) then
  call SET_SIGMA
endif

  write(*,*) 'Mesh is read     ', 'nod2D=', nod2D,' elem2D=', elem2D
  write(*,*) 'Mesh includes    ', elem2D_tri, 'triangles'
  write(*,*) 'Mesh includes    ', elem2D-elem2D_tri, 'quads'
  write(*,*) 'Amount of open boundary nodes ', nobn
  write(*,*) 'The list of elements starts from triangles !!!'
END SUBROUTINE read_mesh2D

SUBROUTINE array_setup
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

  integer            :: sbc_alloc

  allocate(mslp(nod2D), STAT=sbc_alloc )
  allocate(emp(nod2D), STAT=sbc_alloc )
  mslp = P_ref
  emp = 0.0_WP
  if( sbc_alloc /= 0 )   STOP 'array_setup: failed to allocate arrays'

  allocate(U_n_2D(2, elem2D), U_n_1(2, elem2D), U_n_2(2, elem2D))
  allocate(U_n_2Dav(2, elem2D))
  Allocate(U_n_2D_old(2,elem2D),U_r(2,elem2D)) 
  U_n_2D=0.0_WP
  U_n_1=0.0_WP
  U_n_2=0.0_WP
  U_n_2Dav=0.0_WP
  U_n_2D_old=0.0_WP
  U_r=0.0_WP

 allocate(dmean_n(elem2D))
 allocate(eta_p(nod2D))
 eta_p=0.0_WP

 allocate(UAB(2, elem2D), U_rhs_2D(2, elem2D))
 UAB=0.0_WP
 U_rhs_2D=0.0_WP
 allocate(U_rhs_2D_3D(2,elem2D))
 U_rhs_2D_3D=0.0_WP
 allocate( UV2_rhs(2,elem2D))
 UV2_rhs = 0.0_WP

 allocate(eta_n(nod2D), eta_n_1(nod2D), eta_n_2(nod2d))
 eta_n=0.0_WP
 eta_n_1=0.0_WP
 eta_n_2=0.0_WP

 allocate(etaAB(nod2d),ssh_rhs(nod2D), ssh_gp(nod2D))
 allocate(vel_grad(4,elem2D))
 U_rhs_2D=0.0_WP
 ssh_rhs=0.0_WP
 vel_grad=0.0_WP
 ssh_gp=0.0_WP

 allocate(taux(elem2D), tauy(elem2D))
 allocate(taux_node(nod2D), tauy_node(nod2D))
 taux=0.0_WP
 tauy=0.0_WP
 taux_node=0.0_WP
 tauy_node=0.0_WP

 allocate(windx(nod2D), windy(nod2D))
 windx=0.0_WP
 windy=0.0_WP

 allocate(relax_coef(nod2D))
 relax_coef=0.0_WP

 ! add 2D part
 allocate(Bar_pr_2D(2,elem2D),hpre_2D(nod2D))
 Bar_pr_2D=0.0_WP
 hpre_2D=0.0_WP

! 3D part
!VF, allocate depends on the type_task
if (type_task>1) then
allocate(TF(nsigma-1,nod2D), SF(nsigma-1,nod2D),TF2(nsigma-1,nod2D), SF2(nsigma-1,nod2D))
allocate(T_old(nsigma-1,nod2D), S_old(nsigma-1,nod2D),T_old2(nsigma-1,nod2D), S_old2(nsigma-1,nod2D))
allocate(Tclim(nsigma-1,nod2D), Sclim(nsigma-1,nod2D))
endif
!if (comp_sediment) allocate(CF(nsigma-1,nod2D),c_old(nsigma-1,nod2D),w_s(nsigma-1,nod2D),Cclim(nsigma-1,nod2D))
allocate(rho_c(nsigma-1,nod2D))
rho_c=0.0_WP

Allocate(C_d_el(elem2D))

if (type_task>1) then
allocate(z0b_gotm_el(elem2D))
z0b_gotm_el=z0b_min
allocate(U_n(nsigma-1,elem2D), V_n(nsigma-1,elem2D))
U_n=0.0_WP
V_n=0.0_WP
allocate(hpressure(nsigma,nod2D))
hpressure=0.0_WP
allocate(Bar_pru_3D(nsigma-1,elem2D), Bar_prv_3D(nsigma-1,elem2D))
Bar_pru_3D=0.0_WP
Bar_prv_3D=0.0_WP
allocate(Bar_pru_3D_clim(nsigma-1,elem2D), Bar_prv_3D_clim(nsigma-1,elem2D))
Bar_pru_3D_clim=0.0_WP
Bar_prv_3D_clim=0.0_WP
 allocate(snu(nsigma,nod2D), bt(nsigma,nod2D))
 snu = snul
 bt = 1.0d-8
 allocate(UV_rhs(2,nsigma-1,elem2D))
 UV_rhs = 0.0_WP

Allocate(Jc(nsigma-1,nod2D),Jc_old(nsigma-1,nod2D))
Allocate(Je(nsigma-1,elem2D),Jd(nsigma-1,edge2D))
Jc=0.0_WP
Jc_old=0.0_WP
Je=0.0_WP
Jd=0.0_WP
Allocate(Unode(nsigma-1,nod2D), Vnode(nsigma-1,nod2D))
Unode=0.0_WP
Vnode=0.0_WP
allocate(Kv(nsigma,nod2D), Kv2(nsigma,nod2D),Av_node(nsigma,nod2D),Av(nsigma,elem2D))
allocate(L(nsigma,nod2D), teps(nsigma,nod2D), tepsb(nsigma,nod2D),tke(nsigma,nod2D))
tke=1.0e-8
teps=1.0e-12
Kv=snul
Kv2=snul
Av=snul
Av_node=snul
allocate(U_rhs(nsigma-1,elem2D), V_rhs(nsigma-1,elem2D))
Allocate(U_rhsAB(nsigma-1,elem2D), V_rhsAB(nsigma-1,elem2D))
allocate(vel_grad_ux(nsigma-1,elem2D), vel_grad_uy(nsigma-1,elem2D))
allocate(vel_grad_vx(nsigma-1,elem2D), vel_grad_vy(nsigma-1,elem2D))
vel_grad_ux=0.0_WP
vel_grad_uy=0.0_WP
vel_grad_vx=0.0_WP
vel_grad_vy=0.0_WP
allocate(Visc(nsigma-1,elem2D))
Visc=0.0_WP

Allocate(U_puls(nsigma-1,elem2D), V_puls(nsigma-1,elem2D))
Allocate(Unode_p(nsigma-1,nod2D), Vnode_p(nsigma-1,nod2D))
Allocate(U_n_filt(nsigma-1,elem2D), V_n_filt(nsigma-1,elem2D))

Allocate(W_n(nsigma,nod2D))
Unode_p=0.0_WP
Vnode_p=0.0_WP
U_puls=0.0_WP
V_puls=0.0_WP
U_n_filt=0.0_WP
V_n_filt=0.0_WP
U_rhsAB=0.0_WP
V_rhsAB=0.0_WP
U_rhs=0.0_WP
V_rhs=0.0_WP
W_n=0.0_WP

Allocate(vorticity_2D(nod2D), Vorticity_3D(nsigma-1,nod2D))

allocate(Z(nsigma-1,nod2D),zbar(nsigma,nod2D))

Allocate(Wvel(nsigma,nod2D))
Wvel=0.0_WP

allocate(edge_up_dn_grad(4,nsigma-1,edge2D))

edge_up_dn_grad=0.0_WP

endif

allocate(Visc2D(elem2D))
Visc2D=0.0_WP
Allocate(edge_up_dn_tri(2,edge2D))
edge_up_dn_tri=0
allocate(relax2clim(nod2D))
relax2clim = 0.0_WP          ! relaxation to climatology
Allocate(U_filt_2D(elem2D), V_filt_2D(elem2D))
U_filt_2D=0.0_WP
V_filt_2D=0.0_WP
allocate(mask_wd(elem2D))
mask_wd = 1.0_WP          ! wet elements
allocate(mask_wd_node(nod2D),mask_wd_node_old(nod2D))
mask_wd_node = 1.0_WP          ! wet elements
mask_wd_node_old = 1.0_WP          ! wet elements
Allocate(mask_ad(nod2D))
mask_ad=1.0_WP            ! mask for order advection scheme for tracer
                                    ! 0.0  ---> upwind
			           ! 1.0 ---> 2 order (MIURA....)
				   ! 0.0 ---> 1.0 reduce accuracy near critical depth
				   ! H_advtr_crit  ---> if full depth > H_adv_crit ---> high order scheme for tracer
				   ! H_advtr_min ---> if full depth < H_advtr_min ---> UPWIND
ALLOCATE(mask_bpr(elem2D))
mask_bpr=1.0_WP


!++++++++++++++++++++++++++++++++++++++++++++++++++++
! SEDIMENT
!++++++++++++++++++++++++++++++++++++++++++++++++++++
Allocate(qbu(elem2D), qbv(elem2D))
qbu = 0.0_WP
qbv = 0.0_WP
Allocate(hama_v(nod2D), E_sed(nod2D), h_var(nod2D), h_var_old(nod2D))
Allocate(Er_Dep(nod2D), h_var2(nod2D), h_var_old2(nod2D))
Allocate(qb(2,nod2D))
qb=0.0_WP
hama_v = 0.0_WP
E_sed = 0.0_WP
Er_Dep = 0.0_WP
h_var = 0.0_WP
h_var_old = 0.0_WP
h_var2 = 0.0_WP
h_var_old2 = 0.0_WP
allocate(con(nod2D), con_bc(nod2D))
con = 0.0_WP
con_bc = 0.0_WP

END SUBROUTINE array_setup

!====================================================================================
subroutine jacobian
!++++++++++++++++++++++++++++++++++
! Compute vertical Jacobian
! 19.09.2014
! Androsov Alexey
!++++++++++++++++++++++++++++++++++
USE o_MESH
USE o_ARRAYS
USE o_PARAM
!a USE g_PARSUP
IMPLICIT NONE

INTEGER        :: n, nz, el, ed
real(kind=WP)  :: a


!$OMP DO
do n=1,nod2D
   Jc_old(1:nsigma-1,n) = Jc(1:nsigma-1,n)

   a = max(depth(n) + eta_n(n), Dmin)

   zbar(1:nsigma,n) = (1.0_WP - sigma(1:nsigma))*a

   Jc(1:nsigma-1,n) = a*(    sigma(1:nsigma-1) -   sigma(2:nsigma))

   Z(1:nsigma-1,n)  = 0.5_WP*(zbar(1:nsigma-1,n) + zbar(2:nsigma,n))

end do
!$OMP END DO
!$OMP DO
do el=1,elem2D

   Je(1:nsigma-1,el) = w_cv(1,el)*Jc(1:nsigma-1,elem2D_nodes(1,el)) &
                     + w_cv(2,el)*Jc(1:nsigma-1,elem2D_nodes(2,el)) &
                     + w_cv(3,el)*Jc(1:nsigma-1,elem2D_nodes(3,el)) &
                     + w_cv(4,el)*Jc(1:nsigma-1,elem2D_nodes(4,el))
enddo
!$OMP END DO NOWAIT
!$OMP DO
DO ed=1, edge2D
   Jd(1:nsigma-1,ed) = 0.5_WP*(Jc(1:nsigma-1,edge_nodes(1,ed)) + Jc(1:nsigma-1,edge_nodes(2,ed)))
enddo
!$OMP END DO


End subroutine jacobian
!===========================================================================
subroutine wad_mask
! IK: Calculate masks (1/0) for wetting and drying on elemens and on nodes.
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

integer         :: el,nod
mask_wd_node_old=mask_wd_node
!$OMP DO
do el=1,elem2D

  !NR  .false. = dry element, .true. = wet element

  mask_wd(el) = merge(1.0_WP,0.0_WP,( min(depth(elem2D_nodes(1,el)),depth(elem2D_nodes(2,el)),   &
                      depth(elem2D_nodes(3,el)),depth(elem2D_nodes(4,el)))   &
                + max(etaAB(elem2D_nodes(1,el)),etaAB(elem2D_nodes(2,el)),   &
                      etaAB(elem2D_nodes(3,el)),etaAB(elem2D_nodes(4,el))) > Dmin))

 enddo
!$OMP END DO NOWAIT
!$OMP DO
 do nod=1,nod2D
    mask_wd_node(nod) = merge(1.0_WP,0.0_WP,( depth(nod) + etaAB(nod) > Dmin))
 enddo
!$OMP END DO
end subroutine wad_mask


subroutine wad_mask_part
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

real(kind=WP)   :: dmean
integer         :: el,nod

!$OMP DO
do nod=1,nod2D
   dmean = depth(nod) + eta_n(nod)
    if (dmean > Dmin) then
      mask_wd_node(nod) = 1.0_WP  ! high order advection for tracer
    elseif (dmean < Dcr) then
      mask_wd_node(nod) = 0.0_WP
    else
      mask_wd_node(nod) = (dmean - Dcr)/(Dmin - Dcr)
    endif
 enddo
!$OMP END DO NOWAIT
!$OMP DO
do el=1,elem2D
 mask_wd(el) = w_cv(1,el)*mask_wd_node(elem2D_nodes(1,el)) &
             + w_cv(2,el)*mask_wd_node(elem2D_nodes(2,el)) &
             + w_cv(3,el)*mask_wd_node(elem2D_nodes(3,el)) &
             + w_cv(4,el)*mask_wd_node(elem2D_nodes(4,el))
enddo
!$OMP END DO

 end subroutine wad_mask_part

!===========================================================================
subroutine adv_tracer_mask
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
real(kind=WP)   :: dmean,x,y
integer              :: n

mask_ad = 1.0_WP
do n=1,nod2D
   dmean = max(Dmin,(depth(n) + eta_n(n)))
    if (dmean > H_advtr_crit) then
      mask_ad(n) = 1.0_WP  ! high order advection for tracer
    elseif (dmean < H_advtr_min) then
      mask_ad(n) = 0.0_WP  ! first order advection for tracer
    else
      mask_ad(n) = (dmean - H_advtr_min)/(H_advtr_crit - H_advtr_min)
    endif
 enddo

 end subroutine adv_tracer_mask
!===========================================================================
subroutine bpr_mask
   USE o_MESH
   USE o_ARRAYS
   USE o_PARAM

   IMPLICIT NONE

   real(kind=WP)   :: dmean
   integer         :: elem

   mask_bpr = 1.0_WP

   do elem=1,elem2D
      dmean_n(elem) = max(Dmin, w_cv(1,elem)*(eta_n(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
                             +  w_cv(2,elem)*(eta_n(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
                             +  w_cv(3,elem)*(eta_n(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
                             +  w_cv(4,elem)*(eta_n(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem))))

      if (dmean_n(elem) > H_bpr_crit) then
         mask_bpr(elem) = 1.0_WP  ! compute Baroclinic pressure
      elseif (dmean_n(elem) < H_bpr_min) then
         mask_bpr(elem) = 0.0_WP  ! without Baroclinic pressure
      else
         mask_bpr(elem) = (dmean_n(elem) - H_bpr_min)/(H_bpr_crit - H_bpr_min)
      endif
   enddo

 end subroutine bpr_mask
!===========================================================================
!
SUBROUTINE timestep_AB_2D(step)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE fv_sbc

IMPLICIT NONE
real(kind=WP)   :: dmean, fD(4), fDD(4), ibv, rho_inv
integer         :: step, elem, i, elnodes(4),k,ed
integer         :: n, el
! No flux form, only U\nablau

  ! =============
  !  AB3 interpolate
  ! =============
!write(*,*) 'kkp'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,el)
!$OMP DO
do n=1,nod2D
  etaAB(n)  = (1.5_WP+beta)*eta_n(n) - (0.5_WP+2.0_WP*beta)*eta_n_1(n)  + beta*eta_n_2(n)
enddo
!$OMP END DO NOWAIT
!$OMP DO
do el=1,elem2D
  UAB(1,el) = (1.5_WP+beta)*U_n_2D(1,el) - (0.5_WP+2.0_WP*beta)*U_n_1(1,el) + beta*U_n_2(1,el)
  UAB(2,el) = (1.5_WP+beta)*U_n_2D(2,el) - (0.5_WP+2.0_WP*beta)*U_n_1(2,el) + beta*U_n_2(2,el)
enddo
!$OMP END DO
!$OMP END PARALLEL

!write(*,*) 'kkp2'
!==================
!  compute ssh_rhs
!==================
  call compute_ssh_rhs_elem
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO
 do i=1,nod2D
    !Calculate eta_n+1
    if (index_nod2D(i) < 2) ssh_rhs(i) = eta_n(i) + dt_2D*ssh_rhs(i)/area(i)
 enddo

!$OMP END DO
!$OMP END PARALLEL
!NROMP Not nice for parallelization, and not necessary in all 3D and 2D time steps.
!NROMP OpenMP skipped for the first round, but check added to compute this only
!NROMP in time steps with diagnostic output
!write(*,*) 'kkp3p'
  if (mod(n_dt,IREP)==0 .and. step==Mt ) then
     ib_mean=0._WP
     k=0
     do i=1,nod2D
        if (index_nod2D(i) < 2) then
           k = k+1
           !VF, calculate imbalance of the mass conservation law using implicit scheme
           ibv = -eta_n(i) +eta_n_1(i) +dt_2D*ssh_rhs(i)/area(i)
           ib_mean = ib_mean+abs(ibv)
           if (abs(ibv) > ib) then
              ib = abs(ibv)
              im_index = i
           endif
        endif
     ENDDO
     ib_mean = ib_mean / real(k,WP)
  endif

!write(*,*) 'ssh_rhs1',maxval(ssh_rhs)


!write(*,*) 'ssh_rhs,eta_n',maxval(ssh_rhs),maxval(eta_n)
  !write(*,*) 'SSH=',maxval(ssh_rhs), minval(ssh_rhs)
    ! ssh_rhs contains eta_n+1. It will be copied on eta_n, but
    ! after computations of the velocity are made

! calculation of the fluctuation velocity
 if ( (step==1) .and. (type_task>1) ) then
!================================================
! right hand from 3D advection and diffusion to 2D
!================================================
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem,i)
!$OMP DO
    DO elem=1,elem2D
       U_rhs_2D_3D(1,elem) =0.0_WP
       U_rhs_2D_3D(2,elem) =0.0_WP
    ENDDO
!$OMP END DO NOWAIT
!++++++++++++++++++++++++++++++++++
! compute depth and ssh
! on predictor step
!++++++++++++++++++++++++++++++++++
!$OMP DO
   DO elem=1,elem2D

      dmean_n(elem) = max(Dmin, w_cv(1,elem)*(eta_n(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
                             +  w_cv(2,elem)*(eta_n(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
                             +  w_cv(3,elem)*(eta_n(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
                             +  w_cv(4,elem)*(eta_n(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem))))
   enddo

!$OMP END DO NOWAIT
!$OMP DO
   DO i=1,nod2D
      eta_p(i) = eta_n(i)
   END DO
!$OMP END DO
   call compute_puls_vel

   call compute_el2nodes_3D(U_puls,V_puls,Unode_p,Vnode_p)

!$OMP END PARALLEL
!write(*,*) 'dmean_n',  maxval(dmean_n), minval(dmean_n)
!write(*,*) 'Unode_p',  maxval(Unode_p), minval(Unode_p)

   call momentum_adv_P1_3D_to_2D

   if (filt_3D)  call viscosity_filt_3D_to_2D

   if (bih_3D)  call biharmonic_viscosity_3D_to_2D

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem)
   DO elem=1,elem2D
 !NR     elnodes=elem2D_nodes(:,elem)
 !NR     fD=depth(elnodes) + eta_n(elnodes)
 !NR     dmean = max(Dmin,sum(w_cv(1:4,elem)*fD))
 !NR dmean ist already calculated as dmean_n(elem), see loop above
      U_rhs_2D_3D(1,elem) = U_rhs_2D_3D(1,elem)/dmean_n(elem)
      U_rhs_2D_3D(2,elem) = U_rhs_2D_3D(2,elem)/dmean_n(elem)
   END DO
!$OMP END PARALLEL DO
endif
!write(*,*) ' U_rhs_2D_3D(1,elem)',  maxval (U_rhs_2D_3D(1,:))
 call update_2D_vel(step)

  ! =============
  ! Update elevation and velocity
  ! =============

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,n,elem,dmean)
!$OMP DO
do n=1,nod2D
  eta_n_2(n) = eta_n_1(n)
  eta_n_1(n) = eta_n(n)
  eta_n(n)   = ssh_rhs(n)
end do
!$OMP END DO NOWAIT
!$OMP DO
do elem=1,elem2D
  U_n_2(1,elem)      = U_n_1(1,elem)
  U_n_2(2,elem)      = U_n_1(2,elem)
  U_n_1(1,elem)      = U_n_2D(1,elem)
  U_n_1(2,elem)      = U_n_2D(2,elem)
  U_n_2D_old(1,elem) = U_n_2D(1,elem)  
  U_n_2D_old(2,elem) = U_n_2D(2,elem)
  U_n_2D(1,elem)     = U_rhs_2D(1,elem)
  U_n_2D(2,elem)     = U_rhs_2D(2,elem)
enddo
!$OMP END DO NOWAIT

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  preparation for compute Filtering 2D velocity
!
!$OMP DO
  do elem=1,elem2D

    dmean = max(Dmin, w_cv(1,elem)*(eta_n(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
                   +  w_cv(2,elem)*(eta_n(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
                   +  w_cv(3,elem)*(eta_n(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
                   +  w_cv(4,elem)*(eta_n(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem))))

!NR U_rhs_2D allows the loop above to "NOWAIT"
!NR    U_filt_2D(elem) = U_filt_2D(elem) + U_n_2D(1,elem)*dmean
!NR    V_filt_2D(elem) = V_filt_2D(elem) + U_n_2D(2,elem)*dmean
    U_filt_2D(elem) = U_filt_2D(elem) + U_rhs_2D(1,elem)*dmean
    V_filt_2D(elem) = V_filt_2D(elem) + U_rhs_2D(2,elem)*dmean
  enddo
!$OMP END DO
!$OMP END PARALLEL
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! VF calling of sediment procedure is shifted here,
! so do not spoil parallelization
! AA compute sediment model
   If (comp_sediment) call sediment
! AA

END SUBROUTINE timestep_AB_2D
!===========================================================================
SUBROUTINE oce_barotropic_timestep
!
!a modification by Alexey Androsov
!a (for sigma coordinates)
!a 14.10.14
!
USE o_MESH
USE o_ARRAYS
USE o_PARAM
!a USE g_PARSUP
IMPLICIT NONE

integer       :: n, el
real(kind=WP) :: Mt_inv

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el)
do el=1,elem2D
   U_filt_2D(el) = 0.0_WP
   V_filt_2D(el) = 0.0_WP
enddo
!$OMP END PARALLEL DO

do n=1,Mt
   time_2D=time + dt_2D*n
   call timestep_AB_2D(n)
enddo

! Filt. velocity

Mt_inv = 1._WP/real(Mt,WP)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el)
do el=1,elem2D
   U_filt_2D(el) = U_filt_2D(el) * Mt_inv
   V_filt_2D(el) = V_filt_2D(el) * Mt_inv
enddo
!$OMP END PARALLEL DO

end subroutine oce_barotropic_timestep
!===================================================================================
SUBROUTINE solve_tracers
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
real(kind=WP) :: T_min, T_max, S_min, S_max
integer       :: n

T_min = TF(1,1)
T_max = TF(1,1)
S_min = SF(1,1)
S_max = SF(1,1)


!$OMP PARALLEL REDUCTION(min:T_min, S_min) &
!$OMP&         REDUCTION(max:T_max, S_max)
!$OMP DO
do n=1,nod2D
   T_min = min(T_min,minval(TF(:,n)))
   T_max = max(T_max,maxval(TF(:,n)))
enddo
!$OMP END DO NOWAIT
!$OMP DO
do n=1,nod2D
   S_min = min(S_min,minval(SF(:,n)))
   S_max = max(S_max,maxval(SF(:,n)))
enddo
!$OMP END DO
!$OMP END PARALLEL
T_aver=0.5_WP*(T_min + T_max)
S_aver=0.5_WP*(S_min + S_max)

call adv_tracer_mask

 !write(*,*) 'update T_aver, S_aver'


if(tracer_adv==2) then
TF2=TF 
SF2=SF
T_old2=T_old
S_old2=S_old

   call solve_tracer_muscl(TF, T_old, 't')
   call solve_tracer_muscl(SF, S_old, 's')
!   call fct_muscl_2(TF2, T_old2, mu_w) 
!   call fct_muscl_2(SF2, S_old2, mu_w)
!write(*,*) 'maxval, muscl, muscl2', maxval(TF-TF2), minval(TF-TF2),maxval(SF-SF2), minval(SF-SF2)

elseif (tracer_adv==3) then


 call fct_muscl_LH(SF, S_old, mu_w)    
 call fct(SF, S_old, 's')

 call fct_muscl_LH(TF, T_old, mu_w)    

 call fct(TF, T_old, 't')



elseif(tracer_adv==1) then

   call solve_tracer_upwind(TF, T_old, SF, S_old)

elseif (tracer_adv==4) then

   call solve_tracer_miura(TF, T_old, 't')
   call solve_tracer_miura(SF, S_old, 's')

elseif (tracer_adv==5) then
 call fct_miura_LH(TF, T_old)    
 call fct(TF, T_old, 't')
 call fct_miura_LH(SF, S_old)    
 call fct(SF, S_old, 's')

end if



  call tracer_riv(TF, T_old, SF, S_old)
  call tracer_impl_vert_diff

end subroutine solve_tracers
!==========================================================================
!

!===================================================================================
SUBROUTINE solve_tracers_sed
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
real(kind=WP) :: c_min, c_max
integer       :: n

c_min = CF(1,1)
c_max = CF(1,1)

!$OMP PARALLEL REDUCTION(min:c_min) &
!$OMP&         REDUCTION(max:c_max)
!$OMP DO
do n=1,nod2D
   c_min = min(c_min,minval(CF(:,n)))
   c_max = max(c_max,maxval(CF(:,n)))
enddo
!$OMP END DO
!$OMP END PARALLEL
c_aver=0.5_WP*(c_min + c_max)


call adv_tracer_mask
! IK adv_tracer_mask does not work after step back (cover 429) with old/new fv_tracer

STOP "ERROR: check here,fv_main,solve_tracers_sed,call adv_tracer_mask"

 call solve_tracer_upwind_sed(CF, c_old)

 call tracer_impl_vert_diff_sed
 call bottom_evolution2

end subroutine solve_tracers_sed
!==========================================================================
!
SUBROUTINE oce_timestep
!
!a modification by Alexey Androsov
!a (for sigma coordinates)
!a 13.10.14
!

USE o_MESH
USE o_ARRAYS
USE o_PARAM
!a USE g_PARSUP

IMPLICIT NONE
!$ real(kind=WP)      :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
integer                :: tracer_sch, n, nz
if (T_potential) call potential           ! AA compute tidal geopotential here

select case (type_task)

case(1)
!$ if (iverbosity >= 3) t1=omp_get_wtime()
   call oce_barotropic_timestep
!$ if (iverbosity >= 3) then
!$    t2=omp_get_wtime()
!$    write(*,'("oce_barotropic_timestep took ",f10.4,"s")') t2-t1
!$ endif

case(2)

!$ if (iverbosity >= 3) t1=omp_get_wtime()
!$OMP PARALLEL

   call compute_el2nodes_3D(U_n,V_n,Unode,Vnode)
!$OMP END PARALLEL
!+++++++++++++++++++++++++
! vertical mixing scheme
!+++++++++++++++++++++++++
!$ if (iverbosity >= 3) t2=omp_get_wtime()

   if (ver_mix == 1) then
   call oce_mixing_PP
   Kv2=Kv
   endif
   if (ver_mix == 2)then
 call GOTM
!     Av_node = snul !eddy viscosity
!     Kv =  snul    !eddy diffusivity
!     Kv2=   snul
endif
   if (ver_mix == 3) then
   call d3_end
   Kv2=Kv
   endif
!$ if (iverbosity >= 3) t3=omp_get_wtime()

   call compute_vel_rhs

!$ if (iverbosity >= 3) t4=omp_get_wtime()
   call oce_barotropic_timestep
!$ if (iverbosity >= 3) t5=omp_get_wtime()

!$OMP PARALLEL
   call jacobian

!$ if (iverbosity >= 3) t6=omp_get_wtime()
   call update_3D_vel

!$ if (iverbosity >= 3) t7=omp_get_wtime()
   call vert_vel_sigma

!$OMP END PARALLEL

!$ if (iverbosity >= 3) then
!$   t8=omp_get_wtime()
!$    write(*,'("compute_vel_nodes took       ",f10.4," s")') t2-t1
!$    write(*,'("mixing took                  ",f10.4," s")') t3-t2
!$    write(*,'("compute_vel_rhs took         ",f10.4," s")') t4-t3
!$    write(*,'("oce_barotropic_timestep took ",f10.4," s")') t5-t4
!$    write(*,'("jacobian took                ",f10.4," s")') t6-t5
!$    write(*,'("update_3D_vel took           ",f10.4," s")') t7-t6
!$    write(*,'("vert_vel_sigma took          ",f10.4," s")') t8-t7
!$ endif

 !if (comp_sediment) then
!!VF If comp_sediment and no other tracers, the density still changes due to varying concentration
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz)
 !  DO n=1,nod2D
!! VF, compute density in the vertical column and add it to rho_c - current rho, needed in fv_tracer
!! VF, add effect of the susp. sediments
 !     DO nz=1,nsigma-1
 !        call densityJM(T_const, S_const, -Z(nz,n), rho_c(nz,n),CF(nz,n))
 !     END DO
 !  END DO
!!$OMP END PARALLEL DO

!call solve_tracers_sed

!endif

case(3)
if (Mask_Bar_pr) call bpr_mask

!$ if (iverbosity >= 3) t1=omp_get_wtime()
call pressure
! if ((time/3600.0_WP/24.0_WP)>4) call pressure

!$ if (iverbosity >= 3) t2=omp_get_wtime()

!$OMP PARALLEL
call compute_el2nodes_3D(U_n,V_n,Unode,Vnode)
!$OMP END PARALLEL
!$ if (iverbosity >= 3) t3=omp_get_wtime()
!+++++++++++++++++++++++++
! vertical mixing scheme
!+++++++++++++++++++++++++
   if (ver_mix == 1) then
   call oce_mixing_PP
   Kv2=Kv
   endif
   if (ver_mix == 2) then
   call GOTM
  !   Av = snul
  !   Av_node = snul !eddy viscosity
  !   Kv =  snul     !eddy diffusivity
  !   Kv2=   snul
endif
   if (ver_mix == 3) then
   call d3_end
   Kv2=Kv
   endif

!$ if (iverbosity >= 3) t4=omp_get_wtime()
  call compute_vel_rhs

!$ if (iverbosity >= 3) t5=omp_get_wtime()
   call oce_barotropic_timestep

!$ if (iverbosity >= 3) t6=omp_get_wtime()

!$OMP PARALLEL
   call jacobian

!$ if (iverbosity >= 3) t7=omp_get_wtime()
   call update_3D_vel

!$ if (iverbosity >= 3) t8=omp_get_wtime()
   call vert_vel_sigma

!$OMP END PARALLEL

!$ if (iverbosity >= 3) t9=omp_get_wtime()

!   print *,'================================ CALL TRACERS ===================='
   call solve_tracers


  ! if (comp_sediment) call solve_tracers_sed

!$ if (iverbosity >= 3) then
!$   t10=omp_get_wtime()
!$    write(*,'("pressure took                ",f10.4," s")') t2-t1
!$    write(*,'("compute_vel_nodes took       ",f10.4," s")') t3-t2
!$    write(*,'("mixing took                  ",f10.4," s")') t4-t3
!$    write(*,'("compute_vel_rhs took         ",f10.4," s")') t5-t4
!$    write(*,'("oce_barotropic_timestep took ",f10.4," s")') t6-t5
!$    write(*,'("jacobian took                ",f10.4," s")') t7-t6
!$    write(*,'("update_3D_vel took           ",f10.4," s")') t8-t7
!$    write(*,'("vert_vel_sigma took          ",f10.4," s")') t9-t8
!$    write(*,'("solve_tracers took           ",f10.4," s")') t10-t9
!$ endif


end select

END SUBROUTINE oce_timestep
!==========================================================================
program main
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE o_UTILIT
USE turbulence, only: init_turbulence
USE mtridiagonal, only: init_tridiagonal
USE fv_ncoutput
USE fv_sbc
USE fv_obc
USE fv_obc_2d

IMPLICIT NONE

integer      :: n, period_m2, out_fft, istep, ist,elnodes(4), nz, rk, k,rk2, turn_on_riv, n_dt2,i, flag_riv,vsp(11)
integer      :: ssh_fft=111, vel2D_fft=100, vel3Du_fft=3, vel3Dv_fft=4, Tfile=5, Sfile=6, tracer_sch, elem,  nn, sk

integer      :: ed, nodes(2), el(2)

real(kind=WP) :: x, y, dt_2D_old, dt_old, tt, ss, pp, pr
real(kind=WP) :: eout,Efric_b,Efric_s,E_obn_fl,eout_nwd,Efric_b_nwd,Efric_s_nwd,E_obn_fl_nwd
real(kind=WP), external :: theta
real(kind=WP)    :: hour_out, day_out
double precision :: L_min

!character*22                      :: output_vel='Output/UV-NS000-00.dat' 



!$ real(kind=8) :: t0, t1, t2, t3, t4, t5


!open(51,file='max_min_ssh.dat')
!open(52,file='max_min_vel.dat')
open(53,file='energy.dat')
!open(54,file='cpu_time.dat')


!++++++++++++++++++++++++++
! Set model parameters
!++++++++++++++++++++++++++
 call READ_DATA_RUN
 dt_2D_old=dt_2D
 dt_old=dt
!FV The initial time(Julian days) is taken from run file
! and used for the calculation of the total amount of simulated days
! from the beginning of calculation
 time_jd0_orig=time_jd0

!
!++++++++++++++++++++++++++
! Load the mesh
!++++++++++++++++++++++++++

 call read_mesh2D

 !
!++++++++++++++++++++++++++
! Assemble mesh arrays
!++++++++++++++++++++++++++
!
call test_elem
write(*,*) 'test_elem'
call find_edges

write(*,*) 'find_edges'
call find_elem_neighbors
write(*,*) 'find_neigh'
call mesh_arrays1
write(*,*) 'find_arrays'
call mesh_arrays2
write(*,*) 'find_arrays2'
call array_setup
write(*,*) 'array_setup'
call find_up_downwind_triangles
write(*,*) 'find_up_down'
if (type_task>1) call jacobian
    !=====================
	! Initialize fields
	!=====================
     call initial_state
	!=================================
    ! VF, inizialization of river forcing
    !=================================

     if ((riv).or.(riv_ob)) call initial_riv

    write(*,*) 'initialization done'
	if (type_task>2) then
        T_aver=0.5_WP*(minval(TF)+maxval(TF))
        S_aver=0.5_WP*(minval(SF)+maxval(SF))
	T_maxini = maxval(TF)
	write(*,*) 'aver TS= ', T_aver, S_aver
	endif
	!=====================
	! Initialize turbulence related stuff
	!=====================

         if ((type_task>1).and.(ver_mix == 2)) then
	 call init_turbulence(namlst,'gotmturb.nml',nsigma-1)
         call init_tridiagonal(nsigma-1)
         endif

!Av_node=0.0_WP
!Av=0.0_WP
!Kv=0.0_WP
!Kv2=0.0_WP
!L=0.0_WP


	!=====================
	! Time stepping
	!=====================
write(*,*) 'ini_timestep'
!VF, files for fft analysis
	if (T_out .and. TF_presence) then
           if (T_out_el_va) then
              open(ssh_fft,file=trim(OUTDIR)//'ssh_fft.out')
              open(vel2D_fft,file=trim(OUTDIR)//'vel2D_fft.out')
           endif
           if (T_out_3Dvel) then
              open(vel3Du_fft,file=trim(OUTDIR)//'vel3Du_fft.out')
              open(vel3Dv_fft,file=trim(OUTDIR)//'vel3Dv_fft.out')
           endif
	endif


	if(restart) call read_restart
	write(*,*) ini_time
!++++++++++++++++++++++++++++++++++
! compute mask for wetting/drying
!++++++++++++++++++++++++++++++++++
 if (WET_DRY_ON) then
!$OMP PARALLEL
    call wad_mask
!$OMP END PARALLEL
 endif



    time_jd = time_jd0
    turn_on_riv=0

!VF Here is a control when we should turn on rivers
    if (riv_control_ob) then
    riv_OB=.FALSE.
    flag_riv=0
    if (time_jd>=riv_vt_jd(1)) then

    do i=2, size(riv_vt_jd)
    if (time_jd<riv_vt_jd(i)) then
!VF We should turn on rivers immediately and calculate how long the current river record is valid
    rk2=i-1
    riv_vt2(rk2)=(riv_vt_jd(i)-time_jd)/jdh
    flag_riv=1
    riv_OB=.TRUE.
     write(*,*)'Rivers are activated, start from record :', rk2, 'the record valid (hours): ', riv_vt2(rk2)
    endif
    if (flag_riv==1) exit
    enddo
    if (flag_riv==0) then
    riv_control_OB=.FALSE.
    write(*,*) 'Rivers will be not activated, initial Julian date is large than the date of last river record!'
    endif
    else
!VF Here we calculate how many steps should be passed before enabling rivers
!VF (procedure does not take into account possible time step changing!))
    turn_on_riv=nint(((riv_vt_jd(1)-time_jd)/jdh)*3600.0_WP/dt)+1
    write(*,*)'Rivers will be activated after ', turn_on_riv, 'steps'
    rk2=1
    endif
    endif

    if (riv_control) then
    riv=.FALSE.
    flag_riv=0
    if (time_jd>=riv_vt_jd(1)) then

    do i=2, size(riv_vt_jd)
    if (time_jd<riv_vt_jd(i)) then
    rk=i-1
    riv_vt(rk)=(riv_vt_jd(i)-time_jd)/jdh
    flag_riv=1
    riv=.TRUE.
    write(*,*)'Rivers are activated, start from record :', rk, 'the record valid (hours): ', riv_vt(rk)
    endif
    if (flag_riv==1) exit
    enddo
    if (flag_riv==0) then
    riv_control=.FALSE.
    write(*,*) 'Rivers will be not activated, initial julian date large than date of last river record'
    endif
    else
    turn_on_riv=nint(((riv_vt_jd(1)-time_jd)/jdh)*3600.0_WP/dt)+1
    write(*,*)'Rivers will be activated after ', turn_on_riv, 'steps'
    rk=1
    endif
    endif

    !=================================
    ! inizialization of ocean forcing
    !=================================

    if (key_atm) call sbc_ini

    !=================================
    ! inizialization of open boundary
    !=================================

    if (key_obc)    call obc_ini
    if (key_obc_2d) call obc_2d_ini

    !=================================
    ! Initialization of NC output
    !=================================
    if (key_nc_output) call output_ini
   !VF, call Jacobian one more time to initialize Jc_old to non-zero value
   if (type_task>1) call jacobian
   If (comp_sediment) then
   h_var_old = h_var
   h_var_old2 = h_var2
   endif

open(75,file='area-WFALL.dat')
do n=1,nod2D
write(75,*) area(n)
enddo
close(75)
open(75,file='elem_area-WFALL.dat')
do n=1,elem2D
write(75,*) elem_area(n)
enddo
close(75)

DO ed=1, edge2D
   nodes = edge_nodes(:,ed)   
   el = edge_tri(:,ed)
   if (el(1) > 0 .and. el(2) <= 0) then
    write(70,*) nodes
   elseif (el(1) <= 0 .and. el(2) > 0) then
    write(71,*) nodes
   endif
enddo


!===============================
!
!        MAIN LOOP
!
!===============================
   n_dt2=0
   sk=1;
!VF The loop starts from 1 now, nsteps changes (decreasing) if you use hot_start
   do n_dt = 1, nsteps
!write(*,*) "time step: ", n_dt
   n_dt2=n_dt2+1

    if (n_dt <= 350) amA = 0.1_WP
    if (n_dt > 350 .and. n_dt <= 550) amA = 0.2_WP
    if (n_dt > 550 .and. n_dt <= 750) amA = 0.3_WP
    if (n_dt > 750 .and. n_dt <= 1000) amA = 0.4_WP
    if (n_dt > 1000 .and. n_dt <= 1250) amA = 0.5_WP
    if (n_dt > 1250 .and. n_dt <= 1300) amA = 0.6_WP
    if (n_dt > 1300 .and. n_dt <= 1350) amA = 0.7_WP
    if (n_dt > 1350 .and. n_dt <= 1400) amA = 0.8_WP
    if (n_dt > 1400 .and. n_dt <= 1800) amA = 0.9_WP
    if (n_dt > 1800) amA = 1.0_WP


!$ if (iverbosity >= 2) t1=omp_get_wtime()

   time_jd = time_jd+dt/86400.0_WP
   time = time_jd*86400.0_WP

!write(*,*) ' I before sbc'
! surface boundary conditions
!$ if (iverbosity >= 3) t2=omp_get_wtime()


    if (key_atm)   call sbc_do

!$ if (iverbosity >= 3) t3=omp_get_wtime()

    if (key_obc)   call obc_do
!write(*,*) 'here am I before oce_timestep'
    call oce_timestep
!*****************************************
! part of sediment model (every tidal cycle change the bottom due to sediment)
! AA
!write(*,*) 'here am I before h updating'
       If (comp_sediment) then

! AA compute sediment model
!         call sediment
! AA


        if ( mod(n_dt,IREP*6)==0) then
	 do nn=1,nod2D
	  depth(nn) = depth(nn) + (h_var(nn) - h_var_old(nn))
	  h_var_old(nn) = h_var(nn)
!VF: alternative calculation
!aa67          h_var_old2(nn) = h_var2(nn)
         enddo
        endif

        endif

!AA
!*************************************

        !if (key_obc_2d) call obc_2d_do
        !call compute_vortex_2D
        !write(*,*)  ' vortex= ',maxval(vorticity_2D), minval(vorticity_2D)
        !make NetCDF output
!write(*,*) 'here am I before output'

        if (key_nc_output) call output_do

!$ if (iverbosity >= 3) t4=omp_get_wtime()
    !    write(51,'(3e13.5)') (time- time_jd0*86400.0_WP)/3600.0_WP,maxval(eta_n),minval(eta_n)

        dt=dt_old
        dt_2D=dt_2D_old
        If ((T_out).and.(TF_presence)) then
        write(*,*) 'fft writing starts'
           if ((time - time_jd0_orig*86400.0_WP>= T_period*3600.0_WP*T_num).and.(T_counter<T_aval)) then
              dt=T_step
              dt_2D=dt/Mt
              if (T_out_el_va) then
                 write(ssh_fft,'(e14.6)') eta_n
                 write(vel2D_fft,'(2f10.5)') U_n_2D
              endif
              if ((T_out_3Dvel).and.(type_task>1)) then
                 write(vel3Du_fft,'(2f10.5)') U_n
                 write(vel3Dv_fft,'(2f10.5)') V_n
              endif
              T_counter=T_counter+1
           endif
        endif
		
        
!==================================================
!VF, update river input characteristics, if necessary
!==================================================

   if (riv_control)     call update_info_riv(rk,n_dt2, turn_on_riv)
   if (riv_control_ob)  call update_info_riv_ob(rk2,n_dt2, turn_on_riv)

!==================================================
!VF, update sediments concentration at the ob, if necessary
!================================================== 

   if (sed_boun_flux)   call update_info_sed_2D(sk)

!=======================
!  ouput only for control
!=======================
!write(*,*) 'here am I output'
        if ( mod(n_dt,IREP)==0 ) then

        !   call energy(eout,E_obn_fl,Efric_b,Efric_s)
        !   call energy_nwd(eout_nwd,E_obn_fl_nwd,Efric_b_nwd,Efric_s_nwd)
          !write(52,'(4e16.7)') time/86400.0_WP-time_jd0,eout_nwd,E_obn_fl_nwd,Efric_b_nwd
        !   write(53,'(4e16.7)') time/86400.0_WP-time_jd0,eout,E_obn_fl,Efric_b

           if (iverbosity >= 1) then
              write(*,*) 'time= ',time/86400.0_WP-time_jd0,'time_all= ',time_2D/86400.0_WP-time_jd0_orig,&
              time/86400.0_WP-time_jd0_orig
           !   write(*,*)  ' energy= ',eout
          !   write(*,*) 'Mass Conserv. imbalance: max, index_max, mean', ib, im_index, ib_mean
              write(*,*) 'max_min_ssh:' , maxval(eta_n), minval(eta_n)
              write(*,*) 'max_vel 2D:', maxval(U_n_2D),minval(U_n_2D)
           !   write(*,*) 'max-wind:',maxval(taux),maxval(tauy)
              write(*,*) 'min/max-slp',maxval(mslp), minval(mslp)
           !   write(*,*) 'mask', maxval(mask_wd), minval(mask_wd)
              if (type_task>1) then
                 write(*,*) 'max_vel 3D:', maxval(U_n),maxval(V_n)
                 write(*,*) 'vert_vel= ', maxval(Wvel),minval(Wvel)
                 write(*,*) 'maxmin_2Dpr=',maxval(Bar_pr_2D), minval(Bar_pr_2D)
                 write(*,*) 'max_3Dpr=',maxval(Bar_pru_3D), maxval(Bar_prv_3D)
                 write(*,*) 'C_d_el=',maxval(C_d_el), minval(C_d_el)
                 write(*,*) 'max_vert_vis_dif= ', maxval(Av), maxval(Kv), maxval(Kv2)
                 if (type_task>2) then
                    write(*,*) 'T_maxmin= ', maxval(TF),minval(TF)
                    write(*,*) 'rho= ', maxval(rho_c),minval(rho_c)
                    write(*,*) 'S_maxmin= ', maxval(SF),minval(SF)
                    write(*,*) 'dt', dt
                   !write(*,*) 'concentration of sedimentation 3d = ' , maxval(CF),minval(CF)
                 endif
              endif
          !    write(*,*) 'dt_2D', dt_2D
   If (comp_sediment) then
    write(*,*) 'concentration of sedimentation max, min: con = ' , maxval(con), maxval(con)*plop_s*10,minval(con)*plop_s*10
  !  write(*,*) 'concentration of sedimentation 3d = ' , maxval(CF),minval(CF)
   ! write(*,*) 'concentration riv = ' , maxval(CF(:,riv_node(:))),minval(CF(:,riv_node(:)))
   ! write(*,*) 'w_s = ' ,maxval(w_s),minval(w_s)
   ! write(*,*) 'Visc = ' ,maxval(Visc),minval(Visc)
    write(*,*) 'bottom variation max, min: h_var = ' , maxval(h_var),minval(h_var)
   ! write(*,*) 'bottom variation max, min: h_var2 = ' , maxval(h_var2),minval(h_var2)
   ! write(*,*) 'Er_Dep', maxval(Er_Dep),minval(Er_Dep), sum(Er_Dep)/nod2D
   ! write(*,*) 'riv er_dep', maxval(Er_Dep(riv_node(:))),minval(Er_Dep(riv_node(:)))
    write(*,*) 'C_d = ', maxval(C_d_el), minval(C_d_el)
    !write(*,*) 'riv node = ', Unode(:,riv_node(:)), maxval(Vnode(:,riv_node(:)))
    !write(*,*) 'riv node = ', Er_Dep(riv_node(:)), maxval(w_s(:,riv_node(:))),minval(w_s(:,riv_node(:)))
   endif
           endif
        endif

!AA output for control	
  !      if ( mod(n_dt, IREC)==0 ) then
!	   if (type_task>1) then
 !          call cross_sec_LE
!	   write(*,*) 'output cross section LE experiment', time/60.0_WP
           !  call compute_vortex_2D
           !  call compute_vortex_3D
           !  call output_vert_vel
           !  write(ssh_fft,'(e14.6)') eta_n
           !  write(vel2D_fft,'(2f10.5)') U_n_2D
 !         endif
  !       endif

     !        call cross_sec
     !        write(*,*) 'output cross section TF'
     !  endif
       if ( mod(n_dt, IRESTART)==0 ) call write_restart(n_dt)

!$ if (iverbosity >= 2) then
!$    t5=omp_get_wtime()
!$    if (iverbosity >= 3) then
!$       write(*,'("-- sbc_do took           ",f10.4," s --")') t3-t2
!$       write(*,'("-- oce_timestep took     ",f10.4," s --")') t4-t3
!$       write(*,'("-- output+diagnosis took ",f10.4," s --")') t5-t4
!$    endif
!$       write(*, '(" =========== STEP took",f10.4," s =========== " )')  t5-t1
!$       write(*,*)  " "
!$ end if
    end do
    close(52,status='KEEP')  ! close file ---> energy_nwd.dat
    close(53,status='KEEP')  ! close file ---> energy.dat

end program main

