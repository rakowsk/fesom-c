! ============================================================================
SUBROUTINE momentum_adv_p1
!
!a modification by Alexey Androsov
!a (for sigma coordinates)
!a 09.10.14
! 
USE o_MESH
USE o_ARRAYS
USE o_PARAM
!a USE g_PARSUP
IMPLICIT NONE
real(kind=WP)  :: tx, ty, uu, vv
real(kind=WP)  :: un, u1, v1, xe, ye, acc
integer       :: n, elem, nz, k,  m
integer       :: nodes(2), el(2), ed
 
 DO ed=1, edge2D_in
   nodes=edge_nodes(:,ed)   
   el=edge_tri(:,ed)

   xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye=edge_dxdy(2,ed)*r_earth
   acc=sum(ac(nodes))/2.0_WP
! since velocities are obtained by averaging, 
! one cannot introduce clean metrics here!   
! Only faces that do not belong to
! vertical walls can contribute to
! the momentum advection 
   DO nz=1,nsigma-1
   !====== 
   ! The piece below gives second order spatial accuracy for
   ! the momentum fluxes. 
   !======
   
   u1=0.5_WP*(Unode(nz,nodes(1))+Unode(nz,nodes(2)))
   v1=0.5_WP*(Vnode(nz,nodes(1))+Vnode(nz,nodes(2)))
   
   !======
   ! Normal velocity at edge ed directed to el(2)
   ! (outer to el(1)) multiplied with the length of the edge
   !======
   un=u1*ye-v1*xe
   un = un*Jd(nz,ed)
   U_rhsAB(nz,el(1))=U_rhsAB(nz,el(1))-u1*un*acc
   V_rhsAB(nz,el(1))=V_rhsAB(nz,el(1))-v1*un*acc
   U_rhsAB(nz,el(2))=U_rhsAB(nz,el(2))+u1*un*acc
   V_rhsAB(nz,el(2))=V_rhsAB(nz,el(2))+v1*un*acc
   END DO
 END DO  
END subroutine momentum_adv_p1


!==========================================================================================
SUBROUTINE momentum_adv_upwind
! 
! Linear reconstruction upwind horizontal momentum advection.
! It is typically too damping or too noisy, so other options are
! recommended instead.
!
!a modification by Alexey Androsov
!a (for sigma coordinates)
!a 09.10.14
! 
USE o_PARAM
USE o_MESH
USE o_ARRAYS
!a USE g_PARSUP

IMPLICIT NONE
real(kind=WP)    :: u1, u2, v1, v2
real(kind=WP)    :: un, acc
integer          :: ed, el(2), nz
integer          :: nthreads, chunk, me_el1, me_el2, me
real(kind=WP)    :: uu(nsigma-1,edge2D_in), vv(nsigma-1,edge2D_in)


! ======================
! Horizontal momentum advection   -\int div(uu)dS=-\sum u(un)l
! (assembly over edges)        
! ======================
!$OMP PARALLEL PRIVATE(u1,u2,v1,v2,un,acc,ed,el,nz,me_el1, me_el2,chunk,nthreads)
 
!$  nthreads = omp_get_num_threads()
!$  me = omp_get_thread_num()
!$  chunk = elem2D / nthreads + 1 
!$  me_el1 = me*chunk +1
!$  me_el2 = min((me+1)*chunk, elem2D)

!$OMP DO
 DO ed=1, edge2D_in

   el=edge_tri(:,ed)

   acc = 0.5_WP * sum(ac(edge_nodes(:,ed)))
   
! Only faces that do not belong to
! vertical walls can contribute to
! the momentum advection 
   do nz=1,nsigma-1
   !====== 
   ! The piece below gives second order spatial accuracy for
   ! the momentum fluxes. 
   !======
      u1 = U_n(nz,el(1))- vel_grad_ux(nz,el(1))*edge_cross_dxdy(1,ed) &
                         -vel_grad_uy(nz,el(1))*edge_cross_dxdy(2,ed)
      v1 = V_n(nz,el(1))- vel_grad_vx(nz,el(1))*edge_cross_dxdy(1,ed) &
                         -vel_grad_vy(nz,el(1))*edge_cross_dxdy(2,ed)
      u2 = U_n(nz,el(2))- vel_grad_ux(nz,el(2))*edge_cross_dxdy(3,ed) &
                         -vel_grad_uy(nz,el(2))*edge_cross_dxdy(4,ed)
      v2 = V_n(nz,el(2))- vel_grad_vx(nz,el(2))*edge_cross_dxdy(3,ed) &
                         -vel_grad_vy(nz,el(2))*edge_cross_dxdy(4,ed)
   
   !======
   ! Normal velocity at edge ed directed to el(2)
   ! (outer to el(1)) multiplied with the length of the edge
   ! and mean depth dmean
   !======
      un = 0.5_WP* r_earth* Jd(nz,ed)* &
           ((u1+u2)*edge_dxdy(2,ed) - (v1*elem_cos(el(1))+v2*elem_cos(el(2)))*edge_dxdy(1,ed))
      
   !======
   ! If it is positive, take velocity in the left element (el(1)),
   ! and use the velocity at el(2) otherwise.
   !======  
   
      if(un>=0.0_WP) then    
         uu(nz,ed) = u1*un*acc    
         vv(nz,ed) = v1*un*acc
      else    
         uu(nz,ed) = u2*un*acc    
         vv(nz,ed) = v2*un*acc
      end if
   
 
   END DO
 END DO  
!$OMP END DO

!NROMP   To be better parallelized! This workarround is not optimal:
 DO ed=1, edge2D_in

    el=edge_tri(:,ed)

!$ if (el(1) >= me_el1 .and. el(1)<= me_el2 ) then
    U_rhsAB(1:nsigma-1,el(1)) = U_rhsAB(1:nsigma-1,el(1)) -uu(1:nsigma-1,ed)
    V_rhsAB(1:nsigma-1,el(1)) = V_rhsAB(1:nsigma-1,el(1)) -vv(1:nsigma-1,ed)
!$ end if
!$ if (el(2) >= me_el1 .and. el(2)<= me_el2 ) then
    U_rhsAB(1:nsigma-1,el(2)) = U_rhsAB(1:nsigma-1,el(2)) +uu(1:nsigma-1,ed)
    V_rhsAB(1:nsigma-1,el(2)) = V_rhsAB(1:nsigma-1,el(2)) +vv(1:nsigma-1,ed)
!$ end if
 ENDDO

!$OMP END PARALLEL

END SUBROUTINE momentum_adv_upwind

! ============================================================================
SUBROUTINE momentum_vert_adv
! 
! Vertical momentum advection and viscosity
! For advection, quadratic upwind reconstruction is used.
! 
 
  USE o_PARAM
  USE o_MESH
  USE o_ARRAYS
  !USE g_PARSUP

  IMPLICIT NONE

  integer          :: el,  nz, is, ie 
  real(kind=WP)    :: w(nsigma-1), acc
  real(kind=WP)    :: a, b, c, d, dg, db, da, f
  real(kind=WP)    :: Z1(2:nsigma-1), Ze(nsigma-1)

  real(kind=WP)    :: u_mean(nsigma), v_mean(nsigma)

  ! =======================
  ! Vertical momentum advection 
  ! and vertical viscosity
  ! =======================

!$OMP DO 

!NR  I have split the inner loop, nz, to allow for at least
!NR  some vectorization.

  DO el=1, elem2D
  
     If (mask_wd(el) /= 0.0_WP) then

     ! Nodal values projected to the element

     acc=sum(w_cv(1:4,el)*ac(elem2D_nodes(1:4,el)))

     ! Wind stress and bottom drag
     w(1:nsigma-1) = elem_area(el)*(w_cv(1,el)*Wvel(1:nsigma-1,elem2D_nodes(1,el)) &
                                  + w_cv(2,el)*Wvel(1:nsigma-1,elem2D_nodes(2,el)) &
                                  + w_cv(3,el)*Wvel(1:nsigma-1,elem2D_nodes(3,el)) &
                                  + w_cv(4,el)*Wvel(1:nsigma-1,elem2D_nodes(4,el)))

     Z1(2:nsigma-1) = w_cv(1,el)*zbar(2:nsigma-1,elem2D_nodes(1,el)) &
                    + w_cv(2,el)*zbar(2:nsigma-1,elem2D_nodes(2,el)) &
                    + w_cv(3,el)*zbar(2:nsigma-1,elem2D_nodes(3,el)) &
                    + w_cv(4,el)*zbar(2:nsigma-1,elem2D_nodes(4,el))

     Ze(1:nsigma-1) = w_cv(1,el)*Z(1:nsigma-1,elem2D_nodes(1,el)) &
                    + w_cv(2,el)*Z(1:nsigma-1,elem2D_nodes(2,el)) &
                    + w_cv(3,el)*Z(1:nsigma-1,elem2D_nodes(3,el)) &
                    + w_cv(4,el)*Z(1:nsigma-1,elem2D_nodes(4,el))
     

     u_mean(1)= -w(1)*U_n(1,el)
     v_mean(1)= -w(1)*V_n(1,el)
     
     is = 2
     if(w(2) <= 0.0_WP) then
        is=3
        u_mean(2) = -0.5_WP* w(2) *(U_n(1,el)+U_n(2,el))
        v_mean(2) = -0.5_WP* w(2) *(V_n(1,el)+V_n(2,el))
     endif

     ie = nsigma-1   
     if(w(nsigma-1) >=0.0_WP) then 
        ie=nsigma-2
        u_mean(nsigma-1) = -0.5_WP* w(nsigma-1) *(U_n(nsigma-2,el)+U_n(nsigma-1,el))
        v_mean(nsigma-1) = -0.5_WP* w(nsigma-1) *(V_n(nsigma-2,el)+V_n(nsigma-1,el))
        ! or replace this with first
        ! order upwind 
     endif

     u_mean(nsigma)=0.0_WP 
     v_mean(nsigma)=0.0_WP 

     DO nz=is, ie

        a =   Z1(nz) - Ze(nz-1)
        b = -(Z1(nz) - Ze(nz))

        if(w(nz) >=0.0_WP) then

           c = -(Z1(nz) - Ze(nz+1)) 

           d = 1._WP/(( Ze(nz+1)-Ze(nz-1))*(b*b-a*a) - (c*c-a*a)*(Ze(nz)- Ze(nz-1)))

           dg= a*b*(Ze(nz)   - Ze(nz-1))*d
           db=-a*c*(Ze(nz+1) - Ze(nz-1))*d
           !da= 1.0_WP-dg-db

           u_mean(nz) = -w(nz)*( U_n(nz-1,el)*(1.0_WP-dg-db) + U_n(nz,el)*db + U_n(nz+1,el)*dg)
           v_mean(nz) = -w(nz)*( V_n(nz-1,el)*(1.0_WP-dg-db) + V_n(nz,el)*db + V_n(nz+1,el)*dg)
        

        else !            if(w<0.0_WP) then

           f  = Z1(nz) - Ze(nz-2)
           d  = 1._WP/((Ze(nz) - Ze(nz-2))*(a*a-b*b) - (f*f-b*b)*(Ze(nz) - Ze(nz-1)))

           dg = a*b*(Ze(nz) - Ze(nz-1))*d
           db =-b*f*(Ze(nz) - Ze(nz-2))*d
           !da= 1.0_WP-dg-db
           
           u_mean(nz) = -w(nz) *( U_n(nz,el)*(1.0_WP-dg-db) + U_n(nz-1,el)*db + U_n(nz-2,el)*dg)
           v_mean(nz) = -w(nz) *( V_n(nz,el)*(1.0_WP-dg-db) + V_n(nz-1,el)*db + V_n(nz-2,el)*dg)
           
        end if

     END DO
 
	!  if(i_mean_visc) then
	!  uvert(1:2,1)=0.0_WP
	!  uvert(1:2,nsigma)=0.0_WP
	!  else
	!  uvert(1,1)= stress_surf(1,el)*elem_area(el)/density_0
	!  uvert(2,1)= stress_surf(2,el)*elem_area(el)/density_0
        !	!  end if


     DO nz=1,nsigma-1
        U_rhsAB(nz,el) = U_rhsAB(nz,el) +(u_mean(nz)-u_mean(nz+1))*acc 
        !/(zbar(nz,elem)-zbar(nz+1,elem))
        
        V_rhsAB(nz,el) = V_rhsAB(nz,el) +(v_mean(nz)-v_mean(nz+1))*acc 
        !/(zbar(nz,elem)-zbar(nz+1,elem))
        
     END DO
     
  else  ! dry elements
   DO nz=1,nsigma-1
       U_rhsAB(nz,el)=0.0_WP
       V_rhsAB(nz,el)=0.0_WP
   END DO
  endif  
  END DO
!$OMP END DO


END SUBROUTINE momentum_vert_adv

 ! ===================================================================
SUBROUTINE momentum_vert_adv_upwind
!
! Vertical momentum advection
! For advection, quadratic upwind reconstruction is used.
!
 
USE o_PARAM
USE o_MESH
USE o_ARRAYS
!USE g_PARSUP

integer          :: elem, elnodes(4), nz
real(kind=WP)    :: friction
real(kind=WP)    :: w, acc
real(kind=WP)    :: umean, vmean, a, b, c, d, dg, da, db
real(kind=WP)    :: Z1
real(kind=WP), allocatable  :: uvertAB(:,:)

allocate(uvertAB(2,nsigma))

! =======================
! Vertical momentum advection
! and vertical viscosity
! =======================
 uvertAB=0.0_WP
 !!$OMP DO
 DO elem=1, elem2D
 
   If (mask_wd(elem) /= 0.0_WP) then
   
   elnodes=elem2D_nodes(:,elem)
  
      acc=sum(w_cv(1:4,elem)*ac(elnodes))

   DO nz=2, nsigma-1
          w=sum(w_cv(1:4,elem)*Wvel(nz,elnodes))*elem_area(elem)
      !Aumean=0.5_WP*(U_n(nz-1,elem)+U_n(nz,elem))
        !Avmean=0.5_WP*(V_n(nz-1,elem)+V_n(nz,elem))

      if(w>=0.0_WP) then
        if(nz==nsigma-1) then
        umean=0.5_WP*(U_n(nz-1,elem)+U_n(nz,elem))
          vmean=0.5_WP*(V_n(nz-1,elem)+V_n(nz,elem))
                                           ! or replace this with first
                                           ! order upwind 
    else
        Z1 = sum(w_cv(1:4,elem)*zbar(nz,elnodes))
    a=Z1 - sum(w_cv(1:4,elem)*Z(nz-1,elnodes))
    b=sum(w_cv(1:4,elem)*Z(nz,elnodes)) - Z1
    c=sum(w_cv(1:4,elem)*Z(nz+1,elnodes)) - Z1
    d=(c+a)*(b**2-a**2)-(c**2-a**2)*(b+a)
        dg=a*b*(a+b)/d
    db=-a*c*(a+c)/d
    da=1.0_WP-dg-db
    umean=U_n(nz-1,elem)*da+U_n(nz,elem)*db+U_n(nz+1,elem)*dg
    vmean=V_n(nz-1,elem)*da+V_n(nz,elem)*db+V_n(nz+1,elem)*dg
    end if
      end if

      if(w<0.0_WP) then
        if(nz==2) then
        umean=0.5_WP*(U_n(nz-1,elem)+U_n(nz,elem))
          vmean=0.5_WP*(V_n(nz-1,elem)+V_n(nz,elem))
    else
        Z1 = sum(w_cv(1:4,elem)*zbar(nz,elnodes))
    a=sum(w_cv(1:4,elem)*Z(nz,elnodes)) - Z1
    b=Z1 - sum(w_cv(1:4,elem)*Z(nz-1,elnodes))
    c=Z1 - sum(w_cv(1:4,elem)*Z(nz-2,elnodes))
    d=(c+a)*(b**2-a**2)-(c**2-a**2)*(b+a)
        dg=a*b*(a+b)/d
    db=-a*c*(a+c)/d
    da=1.0_WP-dg-db
    umean=U_n(nz,elem)*da+U_n(nz-1,elem)*db+U_n(nz-2,elem)*dg
    vmean=V_n(nz,elem)*da+V_n(nz-1,elem)*db+V_n(nz-2,elem)*dg
    end if
      end if

     uvertAB(1,nz)= -umean*w
     uvertAB(2,nz)= -vmean*w

   END DO

          w=sum(w_cv(1:4,elem)*Wvel(1,elnodes))*elem_area(elem)
     
      uvertAB(1,1)= -w*U_n(1,elem)
      uvertAB(2,1)= -w*V_n(1,elem)

      uvertAB(1,nsigma)=0.0_WP
      uvertAB(2,nsigma)=0.0_WP

   DO nz=1,nsigma-1
       U_rhsAB(nz,elem)=U_rhsAB(nz,elem) +(uvertAB(1,nz)-uvertAB(1,nz+1))*acc
       V_rhsAB(nz,elem)=V_rhsAB(nz,elem) +(uvertAB(2,nz)-uvertAB(2,nz+1))*acc

   END DO
   
  else  ! dry elements
   DO nz=1,nsigma-1
       U_rhsAB(nz,elem)=0.0_WP
       V_rhsAB(nz,elem)=0.0_WP
   END DO
  endif  
END DO
!!$OMP END DO
deallocate(uvertAB) 


END SUBROUTINE momentum_vert_adv_upwind

!========================================================================================
SUBROUTINE momentum_adv_scalar_2D
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

real(kind=WP)      :: un1, un2, ux, vy, x1, y1, x2, y2, dmean ,acc , fU(4), fV(4)
real(kind=WP), allocatable  ::  Unode_2D(:), Vnode_2D(:) 
integer              :: n, elem, elnodes(4)
integer              :: nodes(2), el(2), ed

! ======
! Momentum advection is assembled on on scalar control volumes
! and then averaged to triangles. This removes noise
! ======
 allocate(Unode_2D(nod2D),Vnode_2D(nod2D))
 Unode_2D=0.0_WP
 Vnode_2D=0.0_WP
 ! ==============
 ! internal edges
 ! ==============
DO ed=1, edge2D_in 
   nodes=edge_nodes(:,ed)   
   el=edge_tri(:,ed)
   x1=edge_cross_dxdy(1,ed)
   y1=edge_cross_dxdy(2,ed)
   x2=edge_cross_dxdy(3,ed)
   y2=edge_cross_dxdy(4,ed)
   dmean=sum(etaAB(nodes)+depth(nodes))/2.0_WP
   dmean=max(Dmin,dmean)
   acc=sum(ac(nodes))/2.0_WP
   un1=(UAB(2,el(1))*x1- UAB(1,el(1))*y1)*dmean
   un2=-(UAB(2,el(2))*x2- UAB(1,el(2))*y2)*dmean
   ux=un1*UAB(1,el(1))+un2*UAB(1,el(2))
   vy=un1*UAB(2,el(1))+un2*UAB(2,el(2))
   
   Unode_2D(nodes(1))=Unode_2D(nodes(1))+acc*ux/area(nodes(1))
   Unode_2D(nodes(2))=Unode_2D(nodes(2))-acc*ux/area(nodes(2))
   Vnode_2D(nodes(1))=Vnode_2D(nodes(1))+acc*vy/area(nodes(1))
   Vnode_2D(nodes(2))=Vnode_2D(nodes(2))-acc*vy/area(nodes(2))
 END DO
 ! =============
 ! boundary edges
 ! =============

 DO ed=1+edge2D_in, edge2D 
   nodes=edge_nodes(:,ed)   
   el=edge_tri(:,ed)
   x1=edge_cross_dxdy(1,ed)
   y1=edge_cross_dxdy(2,ed)
   dmean=sum(etaAB(nodes)+depth(nodes))/2.0_WP
   dmean=max(Dmin,dmean)
   acc=sum(ac(nodes))/2.0_WP
   un1=(UAB(2,el(1))*x1- UAB(1,el(1))*y1)*dmean
   ux=un1*UAB(1,el(1))
   vy=un1*UAB(2,el(1))
 
   Unode_2D(nodes(1))=Unode_2D(nodes(1))+acc*ux/area(nodes(1))
   Unode_2D(nodes(2))=Unode_2D(nodes(2))-acc*ux/area(nodes(2))
   Vnode_2D(nodes(1))=Vnode_2D(nodes(1))+acc*vy/area(nodes(1))
   Vnode_2D(nodes(2))=Vnode_2D(nodes(2))-acc*vy/area(nodes(2))
 END DO
 DO elem=1,  elem2D
   elnodes=elem2D_nodes(:,elem)
   fU=Unode_2D(elnodes)
   fV=Vnode_2D(elnodes)
   U_rhs_2D(1,elem)=U_rhs_2D(1,elem)+sum(w_cv(1:4,elem)*fU)*elem_area(elem)
   U_rhs_2D(2,elem)=U_rhs_2D(2,elem)+sum(w_cv(1:4,elem)*fV)*elem_area(elem)
 END DO
 deallocate(Unode_2D,Vnode_2D)
end subroutine momentum_adv_scalar_2D
!========================================================================================
SUBROUTINE riv_mom_adv_2D
!VF, riv_mom_adv_2D
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

real(kind=WP)        ::  dd, dmean
integer              :: n, ne, nel
integer              :: nodes(2)
 ! =============
 ! boundary edges
 ! =============
 
 riv_vel=0.0_WP
 DO n=1, riv_amount
   ne=riv_ind_eg(n)
   dd=edge_leng(ne)
   nodes=edge_nodes(:,ne)   
   ne=ne-edge2D_in

!if (type_task>1) then
dmean=max(sum(eta_n(nodes)+depth(nodes))*0.5_WP,Dmin)
!else
!   dmean=sum(eta_n(nodes)+depth(nodes))*0.5_WP
!   dmean=max(Dmin,dmean)
!endif
! Attention: if type_task = 1, rivers cannot be used!

   riv_vel(1,n)=Qr(ne)**2*scos(n)/(dd*dmean)
   riv_vel(2,n)=Qr(ne)**2*ssin(n)/(dd*dmean)
  ! riv_vel(1,n)=UAB(1,nel)*Qr(ne)*scos(ne)/(dd*dmean)
  ! riv_vel(2,n)=UAB(2,nel)*Qr(ne)*ssin(ne)/(dd*dmean)
 END DO

end SUBROUTINE riv_mom_adv_2D

!========================================================================================
SUBROUTINE riv_mom_adv_3D
!VF, riv_mom_adv_3D
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

real(kind=WP)        ::  dd, dmean, dmean_sq(nsigma-1)
integer              :: n, ne, nel
integer              :: nodes(2)
 ! =============
 ! boundary edges
 ! =============
 
 riv_vel=0.0_WP
 DO n=1, riv_amount
   ne=riv_ind_eg(n)
   nodes=edge_nodes(:,ne)   
   dd=edge_leng(ne)
   ne=ne-edge2D_in
   dmean=sum(eta_n(nodes)+depth(nodes))*0.5_WP
   dmean_sq=max(Dmin,dmean)*(sigma(1:nsigma-1) - sigma(2:nsigma))
   riv_vel_u(1:nsigma-1,n)=scos(n)*(Qr(ne)*Qr_sig(1:nsigma-1,ne))**2/(dd*dmean_sq)
   riv_vel_v(1:nsigma-1,n)=ssin(n)*(Qr(ne)*Qr_sig(1:nsigma-1,ne))**2/(dd*dmean_sq)
  ! riv_vel_u(1:nsigma-1,n)=U_n(1:nsigma-1,nel)*Qr(ne)*scos(ne)*Qr_sig(1:nsigma-1,ne)&
  ! /(dd*Jd(1:nsigma-1,ne+edge2D_in))
  ! riv_vel_v(1:nsigma-1,n)=V_n(1:nsigma-1,nel)*Qr(ne)*ssin(ne)*Qr_sig(1:nsigma-1,ne)&
  ! /(dd*Jd(1:nsigma-1,ne+edge2D_in))

 END DO

end SUBROUTINE riv_mom_adv_3D

! ============================================================
SUBROUTINE momentum_adv_upwind_2D
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

integer       :: el(2), ed
real(kind=WP) :: un, u1, u2, v1, v2,  dmean  ,acc
real(kind=WP) :: uu(edge2D_in), vv(edge2D_in)

! ======
! Upwind momentum advection with linear velocity reconstruction
! ======

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el,ed,un,u1,u2,v1,v2,dmean, acc)
   DO ed=1,edge2D_in
 
   dmean = max(Dmin,0.5_WP * sum(etaAB(edge_nodes(:,ed))+depth(edge_nodes(:,ed))))
   acc   = sum(ac(edge_nodes(:,ed)))*0.5_WP      

   el=edge_tri(:,ed) 
   ! edge_cross_dxdy(1:4,ed)  ! They are in physical measure
   ! edge_dxdy(1:2,ed)        ! They are in radian measure        
   
   !====== 
   ! Linear reconstruction upwind
   !====== 
   u1 = UAB(1,el(1)) - vel_grad(1,el(1))*edge_cross_dxdy(1,ed) -vel_grad(2,el(1))*edge_cross_dxdy(2,ed)
   v1 = UAB(2,el(1)) - vel_grad(3,el(1))*edge_cross_dxdy(1,ed) -vel_grad(4,el(1))*edge_cross_dxdy(2,ed)
   u2 = UAB(1,el(2)) - vel_grad(1,el(2))*edge_cross_dxdy(3,ed) -vel_grad(2,el(2))*edge_cross_dxdy(4,ed)
   v2 = UAB(2,el(2)) - vel_grad(3,el(2))*edge_cross_dxdy(3,ed) -vel_grad(4,el(2))*edge_cross_dxdy(4,ed)
   !======
   ! Normal velocity at edge ed in the direction to el(2)
   ! (outer to el(1)), multiplied with the length of the edge
   ! and mean depth dmean
   !======
   un = 0.5_WP*r_earth*dmean* &
        ((u1+u2)*edge_dxdy(2,ed)-(v1*elem_cos(el(1))+v2*elem_cos(el(2)))*edge_dxdy(1,ed))
   
   !======
   ! If it is positive, take velocity in the left element (el(1)),
   ! and use the velocity at el(2) otherwise.
   !======  
   
   ! According to FVCOM paper 2003, they do not use linear
   ! reconstruction other than to estimate the normal velocity.
   ! To get their scheme, uncomment 4 lines below:
    u1=UAB(1,el(1))
    v1=UAB(2,el(1))
    u2=UAB(1,el(2))
    v2=UAB(2,el(2))
   

   !NR u2=0.5_WP*((un+abs(un))*u1+(un-abs(un))*u2)
   !NR v2=0.5_WP*((un+abs(un))*v1+(un-abs(un))*v2)
   !NR On modern (Xeon) architecture, an "if" is rather cheap
    if (un >= 0._WP) then
       uu(ed) = un*u1*acc
       vv(ed) = un*v1*acc
    else
       uu(ed) = un*u2*acc
       vv(ed) = un*v2*acc
    endif
   ENDDO
!$OMP END PARALLEL DO

!NROMP To be parallelized!
   DO ed=1,edge2D_in
      el=edge_tri(:,ed) 
      U_rhs_2D(1,el(1)) = mask_wd(el(1))*(U_rhs_2D(1,el(1)) -uu(ed)) !*un  ! Changed
      U_rhs_2D(2,el(1)) = mask_wd(el(1))*(U_rhs_2D(2,el(1)) -vv(ed)) !*un
      U_rhs_2D(1,el(2)) = mask_wd(el(2))*(U_rhs_2D(1,el(2)) +uu(ed)) !*un
      U_rhs_2D(2,el(2)) = mask_wd(el(2))*(U_rhs_2D(2,el(2)) +vv(ed)) !*un
   END DO

end subroutine momentum_adv_upwind_2D
! AB3 time stepping, with the momentm advection in the 
! form u\nabla u. It allows to make the code faster, but needs for 
! that special versions of other routines included in this file.  

!=====================================================================================
SUBROUTINE momentum_adv_P1_3D_to_2D
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

integer       :: elem
integer       :: nodes(2), el(2), ed

real(kind=WP) :: xe, ye, un(nsigma-1), acc, uu(nsigma-1), vv(nsigma-1)

real(kind=WP) :: u_e1(edge2D_in), u_e2(edge2D_in)
real(kind=WP) :: v_e1(edge2D_in), v_e2(edge2D_in)

! ======
! Momentum advection is assembled on scalar control volumes
! and then averaged to triangles. This removes noise
! ======
 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem,nodes,el,ed,xe,ye,un,acc,uu,vv)

!=================================== 
!   Only faces that do not belong to
!    vertical walls can contribute to
!         the momentum advection
!===================================
!$OMP DO
DO ed=1, edge2D_in
   nodes = edge_nodes(:,ed)   
   el = edge_tri(:,ed)
   acc= sum(ac(nodes))*.5_WP 

   xe = edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye = edge_dxdy(2,ed)*r_earth
                           ! since velocities are obtained by averaging, 
			   ! one cannot introduce clean metrics here!   



   !===================================
   !        The piece below gives
   !  second order spatial accuracy for
   !        the momentum fluxes. 
   !===================================
   uu(1:nsigma-1) = 0.5_WP*(Unode_p(1:nsigma-1,nodes(1)) + Unode_p(1:nsigma-1,nodes(2)))
   vv(1:nsigma-1) = 0.5_WP*(Vnode_p(1:nsigma-1,nodes(1)) + Vnode_p(1:nsigma-1,nodes(2)))

   !======
   ! Normal velocity at edge ed directed to el(2)
   ! (outer to el(1)) multiplied with the length of the edge)
   !======
   un(1:nsigma-1) =  (uu(1:nsigma-1)*ye -vv(1:nsigma-1)*xe) *Jd(1:nsigma-1,ed)


   !NR Eventually, we are interested in the vertical average. Do it already here, 
   !NR as soon as possible.
   u_e1(ed) = acc* sum(uu(1:nsigma-1) *un(1:nsigma-1) *Je(1:nsigma-1,el(1)))
   v_e1(ed) = acc* sum(vv(1:nsigma-1) *un(1:nsigma-1) *Je(1:nsigma-1,el(1)))
   u_e2(ed) = acc* sum(uu(1:nsigma-1) *un(1:nsigma-1) *Je(1:nsigma-1,el(2)))
   v_e2(ed) = acc* sum(vv(1:nsigma-1) *un(1:nsigma-1) *Je(1:nsigma-1,el(2)))
   

 END DO  
!$OMP END DO
!$OMP END PARALLEL 

!NROMP  Loop to be parallelized!
 DO  ed=1, edge2D_in
    el = edge_tri(:,ed)
    U_rhs_2D_3D(1,el(1)) = U_rhs_2D_3D(1,el(1)) - u_e1(ed)
    U_rhs_2D_3D(2,el(1)) = U_rhs_2D_3D(2,el(1)) - v_e1(ed)
    U_rhs_2D_3D(1,el(2)) = U_rhs_2D_3D(1,el(2)) + u_e2(ed)
    U_rhs_2D_3D(2,el(2)) = U_rhs_2D_3D(2,el(2)) + v_e2(ed)
 ENDDO

! average on vertical   !NR Sum is moved to first loop.
!
!NR DO elem=1,elem2D
!NR    U_rhs_2D_3D(1,elem) = U_rhs_2D_3D(1,elem) + sum(U_3Dadv_puls(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
!NR    U_rhs_2D_3D(2,elem) = U_rhs_2D_3D(2,elem) + sum(V_3Dadv_puls(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
!NR END DO


end subroutine momentum_adv_P1_3D_to_2D
! ============================================================================
SUBROUTINE momentum_adv_scalar_3D
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

real(kind=WP)     :: tvol, tx, ty, uu, vv
real(kind=WP)      :: un, x1, y1, x2, y2, wu(1:nsigma),wv(1:nsigma)
REAL(KIND=WP)      :: wuw(1:nsigma),wvw(1:nsigma), ff1, ff2, ff, acc
real(kind=WP), allocatable  ::  Unode_rhs(:,:), Vnode_rhs(:,:)
integer                :: n, elem, nz, k, elnodes(4), m
integer                :: nl1, nl2, nodes(2), el(2), ed
integer                :: elnodes_1(4), elnodes_2(4)

! ======
! Momentum advection on scalar control volumes
! ======
 allocate(Unode_rhs(nsigma-1,nod2d), Vnode_rhs(nsigma-1, nod2d))
 Unode_rhs=0.0_WP
 Vnode_rhs=0.0_WP
! ============
! The horizontal part
! ============

 DO ed=1, edge2D
 
   nodes=edge_nodes(:,ed)  
   el=edge_tri(:,ed)
   x1=edge_cross_dxdy(1,ed)
   y1=edge_cross_dxdy(2,ed)
   if(el(2)>0) then
    x2=edge_cross_dxdy(3,ed)
    y2=edge_cross_dxdy(4,ed)
   else
    x2=0.0_WP
    y2=0.0_WP
   end if
   DO nz=1, nsigma-1
      un=(V_n(nz,el(1))*x1- U_n(nz,el(1))*y1)
      Unode_rhs(nz,nodes(1))=Unode_rhs(nz,nodes(1))+un*U_n(nz,el(1))/area(nodes(1))
      Unode_rhs(nz,nodes(2))=Unode_rhs(nz,nodes(2))-un*U_n(nz,el(1))/area(nodes(2))
      Vnode_rhs(nz,nodes(1))=Vnode_rhs(nz,nodes(1))+un*V_n(nz,el(1))/area(nodes(1))
      Vnode_rhs(nz,nodes(2))=Vnode_rhs(nz,nodes(2))-un*V_n(nz,el(1))/area(nodes(2))
   END DO
   if(el(2)>0)then
   DO nz=1, nsigma-1
      un=-(V_n(nz,el(2))*x2- U_n(nz,el(2))*y2)
      Unode_rhs(nz,nodes(1))=Unode_rhs(nz,nodes(1))+un*U_n(nz,el(2))/area(nodes(1))
      Unode_rhs(nz,nodes(2))=Unode_rhs(nz,nodes(2))-un*U_n(nz,el(2))/area(nodes(2))
      Vnode_rhs(nz,nodes(1))=Vnode_rhs(nz,nodes(1))+un*V_n(nz,el(2))/area(nodes(1))
      Vnode_rhs(nz,nodes(2))=Vnode_rhs(nz,nodes(2))-un*V_n(nz,el(2))/area(nodes(2))
   END DO
   END IF
 END DO
 
 if (vert_adv == 2) then

! =============
! The vertical part
! =============
 DO ed=1, edge2D
 
   nodes=edge_nodes(:,ed)  
   el=edge_tri(:,ed)
   elnodes_1=elem2D_nodes(:,el(1))
   elnodes_2=elem2D_nodes(:,el(2))
   ff1 = 6.0_WP
   ff2 = 6.0_WP
   if(elnodes_1(1)/=elnodes_1(4)) ff1 = 8.0_WP
   if(elnodes_2(1)/=elnodes_2(4)) ff2 = 8.0_WP
   nl1=nsigma-1
   ! Here 1/6 because the 1/6 of area is related to the edge (for triangle)
   !                               1/8 of area is related to the edge (for quad)

   wu(1)=U_n(1,el(1))*elem_area(el(1))/ff1
   wu(2:nl1)=0.5_WP*(U_n(2:nl1,el(1))+U_n(1:nl1-1,el(1)))*elem_area(el(1))/ff1
   wu(nl1+1)=0.0_WP
  
   wv(1)=V_n(1,el(1))*elem_area(el(1))/ff1
   wv(2:nl1)=0.5_WP*(V_n(2:nl1,el(1))+V_n(1:nl1-1,el(1)))*elem_area(el(1))/ff1
   wv(nl1+1)=0.0_WP
  
   wuw(1:nl1+1)=wu*Wvel(1:nl1+1,nodes(1))
   wvw(1:nl1+1)=wv*Wvel(1:nl1+1,nodes(1))
   DO nz=1,nl1
   Unode_rhs(nz,nodes(1))=Unode_rhs(nz,nodes(1))-(wuw(nz)-wuw(nz+1))/ &
    (zbar(nz+1,nodes(1))-zbar(nz,nodes(1)))/area(nodes(1))
   Vnode_rhs(nz,nodes(1))=Vnode_rhs(nz,nodes(1))-(wvw(nz)-wvw(nz+1))/ &
    (zbar(nz+1,nodes(1))-zbar(nz,nodes(1)))/area(nodes(1))
   END DO
  
   wuw(1:nl1+1)=wu*Wvel(1:nl1+1,nodes(2))
   wvw(1:nl1+1)=wv*Wvel(1:nl1+1,nodes(2))
   DO nz=1,nl1
   Unode_rhs(nz,nodes(2))=Unode_rhs(nz,nodes(2))-(wuw(nz)-wuw(nz+1))/ &
    (zbar(nz+1,nodes(2))-zbar(nz,nodes(2)))/area(nodes(2))
   Vnode_rhs(nz,nodes(2))=Vnode_rhs(nz,nodes(2))-(wvw(nz)-wvw(nz+1))/ &
    (zbar(nz+1,nodes(2))-zbar(nz,nodes(2)))/area(nodes(2))
   END DO
  
   if(el(2)>0) then
   wu(1)=U_n(1,el(2))*elem_area(el(2))/ff2
   wu(2:nl1)=0.5_WP*(U_n(2:nl1,el(2))+U_n(1:nl1-1,el(2)))*elem_area(el(2))/ff2
   wu(nl1+1)=0.0_WP
  
   wv(1)=V_n(1,el(2))*elem_area(el(2))/ff2
   wv(2:nl1)=0.5_WP*(V_n(2:nl1,el(2))+V_n(1:nl1-1,el(2)))*elem_area(el(2))/ff2
   wv(nl1+1)=0.0_WP
  
   wuw(1:nl1+1)=wu*Wvel(1:nl1+1,nodes(1))
   wvw(1:nl1+1)=wv*Wvel(1:nl1+1,nodes(1))
   DO nz=1,nl1
   Unode_rhs(nz,nodes(1))=Unode_rhs(nz,nodes(1))-(wuw(nz)-wuw(nz+1))/ &
    (zbar(nz+1,nodes(1))-zbar(nz,nodes(1)))/area(nodes(1))
   Vnode_rhs(nz,nodes(1))=Vnode_rhs(nz,nodes(1))-(wvw(nz)-wvw(nz+1))/ &
    (zbar(nz+1,nodes(1))-zbar(nz,nodes(1)))/area(nodes(1))
   END DO
  
   wuw(1:nl1+1)=wu*Wvel(1:nl1+1,nodes(2))
   wvw(1:nl1+1)=wv*Wvel(1:nl1+1,nodes(2))
   DO nz=1,nl1
   Unode_rhs(nz,nodes(2))=Unode_rhs(nz,nodes(2))-(wuw(nz)-wuw(nz+1))/ &
    (zbar(nz+1,nodes(2))-zbar(nz,nodes(2)))/area(nodes(2))
   Vnode_rhs(nz,nodes(2))=Vnode_rhs(nz,nodes(2))-(wvw(nz)-wvw(nz+1))/ &
    (zbar(nz+1,nodes(2))-zbar(nz,nodes(2)))/area(nodes(2))
   END DO
   end if
 END DO
 
 endif

DO elem=1, elem2D
   elnodes=elem2D_nodes(:,elem)
   acc = sum(w_cv(1:4,elem)*ac(elnodes))
   DO nz=1,nsigma-1
   ff = elem_area(elem)*Je(nz,elem)*acc
   U_rhsAB(nz,elem)=U_rhsAB(nz,elem)+ff*sum(w_cv(1:4,elem)*Unode_rhs(nz,elnodes))
   V_rhsAB(nz,elem)=V_rhsAB(nz,elem)+ff*sum(w_cv(1:4,elem)*Vnode_rhs(nz,elnodes))
   END DO
 END DO

   deallocate(Vnode_rhs,Unode_rhs)
END subroutine momentum_adv_scalar_3D
 ! ===================================================================
