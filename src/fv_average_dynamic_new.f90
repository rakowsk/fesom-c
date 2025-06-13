SUBROUTINE compute_ssh_rhs_elem
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE fv_obc_2d

IMPLICIT NONE
!
! -div(Hu) is computed in cycle over edges
!
! The depth is estimated at elements
! c1 is -u_n*L_left*d1, c2  -u_n*L_right*d1

integer        :: ed, el(2), elem, j, n, q, i
real(kind=WP)  :: c1(edge2D)
real(Kind=WP)  :: aux_elem(elem2D)


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,ed,el,elem)
!$OMP DO
DO n=1,nod2D
   ssh_rhs(n)=0.0_WP
ENDDO
!$OMP END DO NOWAIT

 ! ==============
 ! fill mean depth
 ! ==============
!$OMP DO
  DO elem=1,elem2D

    aux_elem(elem) = w_cv(1,elem) *max(Dmin, etaAB(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
                   + w_cv(2,elem) *max(Dmin, etaAB(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
                   + w_cv(3,elem) *max(Dmin, etaAB(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
                   + w_cv(4,elem) *max(Dmin, etaAB(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem)))
 END DO
!$OMP END DO

 ! ==============
 ! internal edges
 ! ==============
!$OMP DO
 DO ed=1,edge2D_in

    el     = edge_tri(:,ed)

    c1(ed) = ( UAB(2,el(1))*edge_cross_dxdy(1,ed) -UAB(1,el(1))*edge_cross_dxdy(2,ed)) *aux_elem(el(1))  &
            +(-UAB(2,el(2))*edge_cross_dxdy(3,ed) +UAB(1,el(2))*edge_cross_dxdy(4,ed)) *aux_elem(el(2))
 END DO
!$OMP END DO NOWAIT
 ! =============
 ! boundary edges
 ! only the left element (1) is available
 ! =============

!$OMP DO
  DO ed=1+edge2D_in, edge2D
     c1(ed) = ( UAB(2,edge_tri(1,ed))*edge_cross_dxdy(1,ed) &
               -UAB(1,edge_tri(1,ed))*edge_cross_dxdy(2,ed) ) * aux_elem(edge_tri(1,ed))
  END DO
!$OMP END DO
!$OMP END PARALLEL

!NROMP  To be parallelized:
  DO ed=1,edge2D
     ssh_rhs(edge_nodes(1,ed)) = ssh_rhs(edge_nodes(1,ed))+c1(ed)
     ssh_rhs(edge_nodes(2,ed)) = ssh_rhs(edge_nodes(2,ed))-c1(ed)
  END DO
!VF, riv
  if (riv) then
 DO n=1,riv_num_nodes
     ssh_rhs(riv_node(n)) = ssh_rhs(riv_node(n)) + Qr_node(n)
  END DO
endif

!VF, sediment
  if (comp_sediment) then
ssh_rhs = ssh_rhs + Qr_sed
endif



  !    open boundary nodes

!VF, TF_presence....., ssh_rhs
if (TF_presence)  then
ssh_rhs(in_obn)=0.0_WP
     do i=1,nobn
      n=in_obn(i)
           q=i
           do j=1,15
              ssh_rhs(n) = ssh_rhs(n)+ampt(q,j)*cos((time_2D-time_jd0_orig*86400.0_WP)&
              *2.0_WP*pi/(Harmonics_tf_period(j)*3600._WP) - fazt(q,j))
           end do
     enddo
  endif
!VF: Attention: here we assume that the phases and amplitudes are already corrected 
! by user according to the initial time of the calculation 

!VF, riv_OB
  if (riv_OB)  then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO
     do i=1,riv_amount_ob
              ssh_rhs(riv_node_ob(i)) = riv_elev(i)
     end do
!$OMP END DO
!$OMP END PARALLEL
  endif

END subroutine compute_ssh_rhs_elem
!===========================================================================
SUBROUTINE update_2D_vel(step)
! It is intended for AB3
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE fv_sbc

!
IMPLICIT NONE
integer, intent(in)   :: step
integer               :: el, elem, elnodes(4), n, rie, kk
real(kind=WP)    :: dmean, dmean2, ff, friction
real(kind=WP)    :: InBr           ! inverse barometer
real(kind=WP)    :: fD(4), fDD(4), eta(nod2D), dmeanv(riv_amount)
real(kind=WP), parameter ::  de=0.614_WP
real(kind=WP), parameter ::  ga=0.088_WP
real(kind=WP), parameter ::  ep=0.013_WP
real(kind=WP), parameter ::  fe=1.0_WP-de-ga-ep
real(kind=WP)    :: acc, a_Cd, density0_inv, C_d2, C_d3, sqfv

 call vel_gradients_2D(UAB)   !AA ! We need them only for biharmonic viscosity
                               !NROMP and for something else? It may not be skipped,
                               !NROMP even if filtering viscosity is used.
 call bottom_friction_C_d_el(step)
! ====================
! Sea level contribution   -D\nabla\eta
! and  the Coriolis force (elemental part)
! ====================
    ! Below are the interpolation coefficients for the elevation
    ! as recommended in Shchepetkin and McWilliams (2005).
    ! CFL stabily limit is essentially improved with them.
    ! Interpolate AM4

density0_inv = 1._WP/(density_0)

!if (T_potential) call potential           ! AA geopotential now computed in:   oce_timestep


!$OMP PARALLEL PRIVATE(n,density0_inv,el,dmean)
!$OMP DO
   DO n=1,nod2D
      eta_n_2(n) = de*ssh_rhs(n) + fe*eta_n(n) + ga*eta_n_1(n) + ep*eta_n_2(n)
      etaAB(n)   = eta_n_2(n)
       eta(n)     = g*(a_tp(n,1)*eta_n_2(n) - emp(n)*density0_inv - a_tp(n,2)*ssh_gp(n)) &
                  + (P_ref - mslp(n))*density0_inv ! Use AM4 interpolation


   END DO

!$OMP END DO


!++++++++++++++++++++++++++++++++++
! compute mask for wetting/drying
!++++++++++++++++++++++++++++++++++

if (WET_DRY_ON) call wad_mask

!$OMP DO
  DO el=1,elem2D

    dmean =  max(Dmin, w_cv(1,el)*(etaAB(elem2D_nodes(1,el)) + depth(elem2D_nodes(1,el))) &
                    +  w_cv(2,el)*(etaAB(elem2D_nodes(2,el)) + depth(elem2D_nodes(2,el))) &
                    +  w_cv(3,el)*(etaAB(elem2D_nodes(3,el)) + depth(elem2D_nodes(3,el))) &
                    +  w_cv(4,el)*(etaAB(elem2D_nodes(4,el)) + depth(elem2D_nodes(4,el))))


    U_rhs_2D(1,el) = elem_area(el) * (-dmean*( gradient_sca(1,el)*eta(elem2D_nodes(1,el)) &
                                             + gradient_sca(2,el)*eta(elem2D_nodes(2,el)) &
                                             + gradient_sca(3,el)*eta(elem2D_nodes(3,el)) &
                                             + gradient_sca(4,el)*eta(elem2D_nodes(4,el)) ) &
             + mask_wd(el)*UAB(2,el)*dmean*coriolis(el)                                     &
             + taux(el)*density0_inv  )

    U_rhs_2D(2,el) = elem_area(el) * (-dmean*( gradient_sca(5,el)*eta(elem2D_nodes(1,el)) &
                                             + gradient_sca(6,el)*eta(elem2D_nodes(2,el)) &
                                             + gradient_sca(7,el)*eta(elem2D_nodes(3,el)) &
                                             + gradient_sca(8,el)*eta(elem2D_nodes(4,el)) ) &
             - mask_wd(el)*UAB(1,el)*dmean*coriolis(el)                                      &
             + tauy(el)*density0_inv  )

 END DO
!$OMP END DO
!$OMP END PARALLEL

 ! ======================
 ! Momentum advection   -\int div(uDu)dS=-\sum uD(un)l
 ! and viscosity: (assembly over edges)
 ! ======================
 if (mom_adv_2D == 1) call momentum_adv_scalar_2D
 if (mom_adv_2D == 2) call momentum_adv_upwind_2D


 If (filt_2D) call viscosity_filt_2D
 If (filt_bi_2D) call viscosity_filt2x_2D
 IF (bih_2D) call biharmonic_2D

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem,ff,dmean, friction,a_Cd,elnodes)
!AA connection with 3D part of the model. If we compute only barotropic task 2D or 3D this is not needed.
if (type_task>2) then
!$OMP DO SCHEDULE(STATIC)
   do elem=1,elem2D
      U_rhs_2D(1,elem) = U_rhs_2D(1,elem) + U_rhs_2D_3D(1,elem)
      U_rhs_2D(2,elem) = U_rhs_2D(2,elem) + U_rhs_2D_3D(2,elem)
   enddo
!$OMP END DO NOWAIT
endif

 ! ==============
 ! Final update and bottom drag
 ! The contribution from the bottom drag is taken implicitly.
 ! ==============
!VF, task, friction
!$OMP DO SCHEDULE(STATIC)
 DO elem=1,elem2D
    elnodes = elem2D_nodes(:,elem)

    dmean = max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))
       a_Cd = C_d_el(elem)*((1.0_WP - sum(w_cv(1:4,elem)*ac(elnodes)))*Cd_bz + 1.0_WP)
    if (type_task == 1) then

       friction = mask_wd(elem)*a_Cd*sqrt(U_n_2D(1,elem)**2+U_n_2D(2,elem)**2)

       ff = 1._WP/(max(Dmin,sum(w_cv(1:4,elem)*(ssh_rhs(elnodes) + depth(elnodes)))) + friction*dt_2D)

!       U_r(1,elem)= (dmean*U_n_2D(1,elem) &
!            + dt_2d*Bar_pr_2D(1,elem) +&
!            dt_2D*U_rhs_2D(1,elem)/elem_area(elem))*ff

!      U_r(2,elem)= (dmean*U_n_2D(2,elem) &
!           + dt_2d*Bar_pr_2D(2,elem) +&
!           dt_2D*U_rhs_2D(2,elem)/elem_area(elem))*ff 

       U_rhs_2D(1,elem) = ((dmean*U_n_2D(1,elem) &
           + dt_2d*Bar_pr_2D(1,elem) +&
            dt_2D*U_rhs_2D(1,elem)/elem_area(elem))*ff + &
            UV2_rhs(1,elem))

       U_rhs_2D(2,elem) = ((dmean*U_n_2D(2,elem) &
            + dt_2d*Bar_pr_2D(2,elem) +&
            dt_2D*U_rhs_2D(2,elem)/elem_area(elem))*ff + &
            UV2_rhs(2,elem))

    else

       friction = mask_wd(elem)*a_Cd*sqrt(U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2)

       ff = 1._WP/(max(Dmin,sum(w_cv(1:4,elem)*(ssh_rhs(elnodes) + depth(elnodes)))) + friction*dt_2D)

       U_rhs_2D(1,elem) = ((dmean*U_n_2D(1,elem) - &
            dt_2D*friction*(U_n(nsigma-1,elem) - U_n_2D(1,elem)) + dt_2d*Bar_pr_2D(1,elem) +&
            dt_2D*U_rhs_2D(1,elem)/elem_area(elem))*ff + &
            UV2_rhs(1,elem))
       U_rhs_2D(2,elem) = ((dmean*U_n_2D(2,elem) - &
            dt_2D*friction*(V_n(nsigma-1,elem) - U_n_2D(2,elem)) + dt_2d*Bar_pr_2D(2,elem) +&
            dt_2D*U_rhs_2D(2,elem)/elem_area(elem))*ff + &
            UV2_rhs(2,elem))
    endif

    ! Now it contains the updated velocity
 END DO
!$OMP END DO
!$OMP END PARALLEL
!write(*,*) 'U_rhs_2D',  maxval(U_rhs_2D), minval(U_rhs_2D)
!write(*,*) 'elem_area',  maxval(elem_area), minval(elem_area), maxval(Bar_pr_2D(1,:)),  maxval(Bar_pr_2D(2,:))
if (riv) then
 call riv_mom_adv_2D

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(kk,elnodes,ff,el)
!$OMP DO
DO el=1,riv_amount

kk = edge_tri(1,riv_ind_eg(el))

elnodes = elem2D_nodes(:,kk)

ff = 1._WP/max(Dmin,sum(w_cv(1:4,kk)*(ssh_rhs(elnodes) + depth(elnodes))))

U_rhs_2D(:,kk) = U_rhs_2D(:,kk) + dt_2D*ff*riv_vel(:,el)/elem_area(kk)

enddo
!$OMP END DO
!$OMP END PARALLEL
endif

END SUBROUTINE update_2D_vel
!==========================================================================
subroutine compute_vortex_2D
  use o_MESH
  use o_ARRAYS
  use o_PARAM

  implicit none
  integer               :: ed, el(2), enodes(2) , n
  real(kind=WP)     :: deltaX1, deltaX2, deltaY1, deltaY2, c1

 vorticity_2D=0.0_WP

do ed=1,edge2D
   enodes=edge_nodes(:,ed)
   el=edge_tri(:,ed)
   deltaX1=edge_cross_dxdy(1,ed)
   deltaY1=edge_cross_dxdy(2,ed)
   c1=deltaX1*U_n_2D(1,el(1))+deltaY1*U_n_2D(2,el(1))
   if(el(2)>0) then
    deltaX2=edge_cross_dxdy(3,ed)
    deltaY2=edge_cross_dxdy(4,ed)
    c1=c1-deltaX2*U_n_2D(1,el(2))-deltaY2*U_n_2D(2,el(2))
   endif
   vorticity_2D(enodes(1))=vorticity_2D(enodes(1))+c1
   vorticity_2D(enodes(2))=vorticity_2D(enodes(2))-c1
end do
do n=1,nod2D
 vorticity_2D(n)=vorticity_2D(n)/area(n)
enddo

end subroutine compute_vortex_2D
!===============================================================================================
subroutine energy(eout,E_obn_fl,Efric_b,Efric_s)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF 08.02.2019
!this procedure calculates the energy balance
!in 2D barotropic case
  implicit none
  integer          :: n, el, elnodes(4),elem,i,enum,ed(2),kk
  real(kind=8)  :: eout, Efric_b,Efric_s,E_obn_fl, fD(4),D,felev


   E_obn_fl=0.0_WP
   Efric_b=0.0_WP
   Efric_s=0.0_WP
   eout=0.0_WP


DO i=1,nobn-1

enum=obn_edge(i)
elem=edge_tri(1,enum)
ed=edge_nodes(:,enum)
D=sum(max(Dmin,depth(ed) + eta_n_1(ed)))*0.5_WP
felev=sum(eta_n_1(ed))*0.5_WP

E_obn_fl=E_obn_fl-edge_leng(enum)*D*(obn_norm(1,i)*U_n_2D_old(1,elem)+obn_norm(2,i)*U_n_2D_old(2,elem))*&
(0.5_WP*(U_n_2D_old(1,elem)**2+U_n_2D_old(2,elem)**2)+g*felev)

end do

E_obn_fl=E_obn_fl*density_0


   DO n=1,nod2D
   eout=eout + 0.5_WP*mask_wd_node_old(n)*g*area(n)*eta_n_1(n)**2
   END DO

   do el=1,elem2D
      elnodes=elem2D_nodes(:,el)
      fD=max(Dmin,depth(elnodes) + eta_n_1(elnodes))
      kk=sum(w_cv(1:4,el)*fD)*elem_area(el)
      eout=eout+mask_wd(el)*0.5_WP*kk*&
                (U_n_2D(1,el)**2+U_n_2D(2,el)**2)

      Efric_b=Efric_b -mask_wd(el)*elem_area(el)*C_d_el(el)*((1.0_WP - sum(w_cv(1:4,el)*ac(elnodes)))*Cd_bz + 1.0_WP)*&
                       sqrt(U_n_2D(1,el)**2+U_n_2D(2,el)**2)*&
                    sqrt(U_n_2D_old(1,el)**2+U_n_2D_old(2,el)**2)*&
                      sqrt(U_r(1,el)**2+U_r(2,el)**2)
    
      Efric_s=Efric_s +elem_area(el)*(U_n_2D_old(1,el)*taux(el)+U_n_2D_old(2,el)*tauy(el))
   end do
   Efric_b=Efric_b*density_0
   eout=eout*density_0   ! Energy for area
!write(*,*) 'E_obn_fl, Efric_b, U_n_2D, U_n_2D_old, U_r, koeff', E_obn_fl, Efric_b, &
!sqrt(U_n_2D(1,el)**2+U_n_2D(2,el)**2),&
!sqrt(U_n_2D_old(1,el-1)**2+U_n_2D_old(2,el-1)**2),&
!sqrt(U_r(1,el)**2+U_r(2,el)**2),&
!((1.0_WP - sum(w_cv(1:4,el)*ac(elnodes)))*Cd_bz + 1.0_WP)

end subroutine energy

subroutine energy_nwd(eout,E_obn_fl,Efric_b,Efric_s)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF 08.02.2019
  implicit none
  integer          :: n, el, elnodes(4),elem,i,enum,ed(2),kk
  real(kind=8)  :: eout, Efric_b,Efric_s,E_obn_fl, fD(4),D,felev

!this procedure calculates the energy balance as above proc.,
!but do not kill potential energy in the dry zone,
!prescribing here the minimum depth
   E_obn_fl=0.0_WP
   Efric_b=0.0_WP
   Efric_s=0.0_WP
   eout=0.0_WP


DO i=1,nobn-1

enum=obn_edge(i)
elem=edge_tri(1,enum)
ed=edge_nodes(:,enum)
D=sum(max(Dmin,depth(ed) + eta_n_1(ed)))*0.5_WP
felev=sum(eta_n_1(ed))*0.5_WP

E_obn_fl=E_obn_fl-edge_leng(enum)*D*(obn_norm(1,i)*U_n_2D_old(1,elem)+obn_norm(2,i)*U_n_2D_old(2,elem))*&
(0.5_WP*(U_n_2D_old(1,elem)**2+U_n_2D_old(2,elem)**2)+g*felev)

end do

E_obn_fl=E_obn_fl*density_0


   DO n=1,nod2D
   eout=eout + 0.5_WP*g*area(n)*eta_n_1(n)**2
   END DO

   do el=1,elem2D
      elnodes=elem2D_nodes(:,el)
      fD=max(Dmin,depth(elnodes) + eta_n_1(elnodes))
      kk=sum(w_cv(1:4,el)*fD)*elem_area(el)
      eout=eout+0.5_WP*kk*&
                (U_n_2D(1,el)**2+U_n_2D(2,el)**2)

      Efric_b=Efric_b -elem_area(el)*C_d_el(el)*((1.0_WP - sum(w_cv(1:4,el)*ac(elnodes)))*Cd_bz + 1.0_WP)*&
                       (U_n_2D(1,el)**2+U_n_2D(2,el)**2)**0.5*&
                    (U_n_2D_old(1,el)**2+U_n_2D_old(2,el)**2)**0.5*&
                      (U_r(1,el)**2+U_r(2,el)**2)**0.5
    
      Efric_s=Efric_s +elem_area(el)*(U_n_2D_old(1,el)*taux(el)+U_n_2D_old(2,el)*tauy(el))
   end do
   Efric_b=Efric_b*density_0
   eout=eout*density_0   ! Energy for area
end subroutine energy_nwd

!=============================================================================================
subroutine C_d_el_calc_2d(elem)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF: please, look on Bi and Toorman, 2015 "Mixed sediment transport modelling in Scheldt estuary
!with a physics-based bottom friction law"
  implicit none

  integer,INTENT(IN)            :: elem
  integer                       :: elnodes(4),n
  real(kind=WP)  :: dmean, vel,sqfv, con_el, h_p, A_p, phi, aux, fA, C_ref, h_ref,snu_aux

 if (comp_sediment) then
    za=2.85_WP*d_m/30.0_WP
 else
    za=0.0_WP
 endif
 
elnodes = elem2D_nodes(:,elem) 
dmean=max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))

!  vel=dsqrt(U_n_2D(1,elem)**2+U_n_2D(2,elem)**2)


!  If (comp_sediment) then
!  con_el=sum(w_cv(1:4,elem)*con(elnodes))*plop_s*dmean
!  else
!  con_el=0.0_WP
!  endif
!  phi=con_el/plop_s !concentr./plop_s
!  !version 1: Suspension viscosity calculation, C_ref, h_ref can be calibrated!!!
!  !C_ref=0.222_WP  !g/l (kg/m**3)
!  !h_ref=0.12_WP   !m
!  !snu_aux=snu_kin*(1+con_el*dmean/(C_ref*h_ref))
 
!  !version 2: approach by Fei and Xiangjun, 1982
!  snu_aux=snu_kin*(1.0_WP-1.35_WP*phi)**(-2.5_WP)
!  !Non-dimensionalised water depth
!  h_p=dsqrt(3*vel*dmean/snu_aux)

!  !Turbulence dampig factor,A_+(A_p) is an empirical value
!   A_p=17.0_WP
!   fA=(1-exp(-h_p/A_p))**2

!  !Suspension friction, empirical value
!   beta=0.045_WP

!  ! Volumetric suspended particle concentration

!   aux=(z0b_min+beta*phi*dmean)*1.5_WP*vel/dmean
!  !Calculation of the squared friction velocity
!  sqfv=fA*(cka*vel/(log(dmean/z0b_min)-1.0_WP+z0b_min/dmean))**2 &
!  + (dsqrt(aux**2 + 3.0_WP*snu_aux*vel/dmean)+aux)**2
!  
!  C_d_el(elem)=min(sqfv/vel**2,0.05)
 
!  Attention the written above approximation works not smoothly everywhere, should be tested more.
!  So, it is commented. Please, uncomment, when it is necessary.

 C_d_el(elem)=(cka/(log(dmean/(z0b_min+za))-1.0_WP))**2
end subroutine C_d_el_calc_2d

!=============================================================================================
subroutine C_d_el_calc_3d(elem)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF
  implicit none

  integer,INTENT(IN)            :: elem
  integer                       :: elnodes(4),j
  real(kind=WP)  :: sqfv, vel, con_el, phi, snu_aux,nu, z0, u_taub
  real(kind=WP)  :: v_sqrt, z0b_gotm, rr



! elnodes = elem2D_nodes(:,elem)

 !vel=dsqrt(U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2)

!If (comp_sediment) then
!con_el=sum(w_cv(1:4,elem)*CF(nsigma-1,elnodes))
!else
!con_el=0.0_WP
!endif
!nu=sum(w_cv(1:4,elem)*Z(nsigma-1,elnodes))/dmean
!nu=1-0.5_WP*(sigma(nsigma-1)-sigma(nsigma))
!phi=con_el/plop_s
!snu_aux=snu_kin*(1.0_WP-1.35_WP*phi)**(-2.5_WP)
!sqfv=vel*(Av(nsigma-1,elem)+snu_aux)/((1.0_WP-nu)*0.5_WP*(sigma(nsigma-1)-sigma(nsigma))*dmean)
!C_d_el(elem)= sqfv/(vel**2)
!endif
! please, look on Bi and Toorman, 2015 "Mixed sediment transport modelling in Scheldt estuary
! with a physics-based botttom friction law" for the above approximation details
 
 nu= max(Je(nsigma-1,elem),dmin)*0.5_WP

 if (comp_sediment) then
    za=2.85_WP*d_m/30.0_WP
 else
    za=0.0_WP
 endif
!  !Iterate bottom roughness length 10 times, see Burchard, 2001; Kagan, 1995
! vel = U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2

! if (vel > 0.0_WP) then
!    v_sqrt = sqrt(vel)*cka
!    z0b_gotm=z0b_min+za
!    u_taub = v_sqrt / log(1._WP + nu/z0b_gotm)

!     DO j=1,2
!     z0b_gotm = 0.1_WP*min(snul/u_taub,1.0_WP)+z0b_min+za
!     u_taub = v_sqrt / log(1._WP + nu/(z0b_gotm))
!     END DO
!     z0b_gotm_el(elem)=z0b_gotm!0.1*snu_kin/max(snu_kin,u_taub)+0.03_WP*(z0b_min+za)
!     C_d_el(elem) = min(u_taub**2/vel,0.05_WP)
!  else
!     z0b_gotm_el(elem)=z0b_min+za
!     C_d_el(elem) = (cka/(log(1._WP + nu/(z0b_min+za))))**2
! endif
  C_d_el(elem) = (cka/(log(1._WP + nu/(z0b_min+za))))**2
  
end subroutine C_d_el_calc_3d

subroutine bottom_friction_C_d_el(step)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF: calculation of C_d_el if it is needed (C_d is not fixed constant)
  implicit none
  integer, intent(in)           :: step
  integer                       :: elem, elnodes(4)
  real(kind=WP)                 :: dmean, Cd_crit_depth, Cd_crit_depth_coef

  if (BFC==1) then
     C_d_el=C_d

  elseif (BFC==3) then
    !IK : increase Cd in the shallow area
     C_d_el = C_d
     Cd_crit_depth      = 3.0  ! critical depth, below this depth Cd will be increased
     Cd_crit_depth_coef = 10.0  !  multiplication coefficient for Cd at depth = 0.0
!$OMP PARALLEL DO private(elem,elnodes,dmean)
     do elem = 1, elem2D
        elnodes = elem2D_nodes(:,elem)
        dmean = sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes)))
        if ( dmean <= Cd_crit_depth ) then
!            C_d_el(elem) = C_d*Cd_crit_depth_coef
            C_d_el(elem) = C_d*( dmean * (1-Cd_crit_depth_coef)/Cd_crit_depth+Cd_crit_depth_coef )
        end if
     enddo
!$OMP END PARALLEL DO

  elseif (BFC==2) then
     if (type_task==1) then
!$OMP PARALLEL DO private(elem)
        do elem=1,elem2D
              call C_d_el_calc_2d(elem)
        enddo
!$OMP END PARALLEL DO
     elseif (step==1) then ! because jacobian calls only every baroclinic step
	 
!$OMP PARALLEL DO private(elem)
        do elem=1,elem2D
           call C_d_el_calc_3d(elem)
        enddo
!$OMP END PARALLEL DO
     endif
  endif

end subroutine bottom_friction_C_d_el
