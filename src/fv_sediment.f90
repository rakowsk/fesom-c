 Subroutine sediment
!++++++++++++++++++++++++++
! sediment transport model: Ver. 01
! 05.04.2016
! Androsov Alexey
! alexey.androsov@awi.de
!++++++++++++++++++++++++++
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  implicit none
  
  integer             :: elnodes(4), elem, n, k

  real(kind=WP) :: Zw, a, b, c, uv, a_min, c_e, s
  real(kind=WP) :: Q_cr, u_st_cr, u_st, Chez, dmean
  real(kind=WP) :: tvol, tx, ty, hama_c, D_st, dmean_MOST
  real(kind=WP) :: Tsp_el, D_st_el, ssk,w_s1
  !real(Kind=WP), allocatable   :: U2D_node(:), V2D_node(:)
  real(Kind=WP), allocatable   ::  Tsp(:)
  real(kind=WP), Allocatable   :: qsu(:), qsv(:)
  
  Allocate(U2D_node(nod2D),V2D_node(nod2D))
U2D_node=0.0_WP
V2D_node=0.0_WP
  Allocate(Tsp(nod2D))
  allocate(qsu(nod2D),qsv(nod2D))
  
!         Kh_c - horizontal diffusion for concentration [1/s]
!         w_s1    --->  settling velocity of a single particle in tranquil water [m/s]
!         D_sed --->  deposition of non-cohesive sediment (k1*l)
!         C_a    --->  near-bed volumetric sediment concentration
!         f_DW ---> Darcy-Weisbach friction factor (const in this simulation = 0.1)
!         u_st   ---> friction velocity
!         u_sur  --->  surface velocity
!         Ren     ---> particle Reynolds number
!         Q_c   ----> critical Shields parameter
!         Q_p   ----> Shields parameter
!
!-------------------------------------------------------------------------------------
!   compute particle parameter D_st
!   Leo and van Rijn "Sediment Transport, Part II:
!   Suspended load Transport". p 1613.
!

    b = 1.0_WP/3.0_WP
    s = plop_s/density_0 - 1.0_WP
    D_st = d_m*(g*s/snu_kin**2)**b
!    write(*,*) 'particle parameter max, min: D_st = ' , D_st
!-------------------------------------------------------------------------------------
!   compute 2D velocity in the node points
!
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,k,tvol,tx,ty)
!$OMP DO
    DO n=1, nod2D 
       tvol=0.0_WP
       tx=0.0_WP
       ty=0.0_WP
       DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
           tvol=tvol+elem_area(elem)
           tx=tx+U_n_2D(1,elem)*elem_area(elem)
	   ty=ty+U_n_2D(2,elem)*elem_area(elem)
        END DO
	U2D_node(n)=tx/tvol
	V2D_node(n)=ty/tvol
     END DO
!$OMP END DO
!$OMP END PARALLEL
   !write(*,*) 'sediment module  = 1 '
          za=2.85_WP*d_m/30.0_WP

         if (D_st <= 4.0_WP) Q_cr = 0.24_WP/D_st
         If (D_st > 4.0_WP .and. D_st <= 10.0_WP) Q_cr = 0.14_WP/D_st**0.64_WP
         If (D_st > 10.0_WP .and. D_st <= 20.0_WP) Q_cr = 0.04_WP/D_st**0.1_WP
         If (D_st > 20.0_WP .and. D_st <= 150.0_WP) Q_cr = 0.013_WP*D_st**0.29_WP
         If (D_st > 150.0_WP) Q_cr = 0.055_WP

	  b = d_m*g*s
          u_st_cr = sqrt(b*Q_cr)
          ssk = (z0b_min+za)*30.0_WP
          c = 13.95_WP*snu_kin/d_m
          w_s1 = sqrt(c*c + 1.09_WP*b ) - c

 !write(*,*) 'critical velocity:', u_st_cr
! write(*,*) u_st_cr 
 !  write(*,*) 'sediment module  = 2 '
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,elem,a,dmean,c,Chez,u_st,Zw,hama_c,b,a_min,c_e,elnodes,Tsp_el,uv)
!$OMP DO

    do n=1,nod2D

!******************************************************************
!   compute critical velocity (Shields criterium) u_st_cr=f(D_st,s,d_m)
!******************************************************************
!************************************
!   compute dynamical velocity  (u_st)
!************************************
          a = sqrt(U2D_node(n)**2 + V2D_node(n)**2)
         !VF z0r is replaced by z0b_min+za, where
         !za is roughness caused by suspended sediment
         !z0b_min is bottom roughness height, minimum value

          dmean=eta_n(n) + depth(n)
          dmean=max(Dmin,dmean)
	  c=max(50.0,12.0_WP*dmean/(ssk))
          Chez = 18.0_WP*dlog10(c)
          u_st=sqrt(g)*a/(Chez)
!          u_st = u_st_BOOK(n)
!	 if (dmean_MOST .gt. 1.0_WP) then
!	   u_st = u_st_MOST(n)
!	 endif


!************************************************************************
! average settling velocity of the material available for transport (w_s1)
! Simpson and Castelltort 2006, Coupled model of surface water flow, 
! sediment transport and morphological evolution. Computers & Geosciences,
! 32, pp. 1600-1614.
! ************************************************************************
 
!w_s1=0.001
!************************
! compute transport stage
!************************
          a = u_st**2/(u_st_cr)**2 - 1.0_WP
	  ! a = u_st_MOST(n)**2/(u_st_cr)**2 - 1.0_WP
	  Tsp(n) = max(0.0_WP,a)

!******************************
! compute weight parameter Zw
!******************************
          Zw = w_s1/(cka*(u_st + 2.0_WP*w_s1))
	  ! Zw = w_s1/(cka*(u_st_MOST(n) + 2.0_WP*w_s1))
	  Zw = min(0.9999_WP,Zw)

          hama_v(n) = 0.980_WP - 0.198_WP*Zw + 0.032_WP*Zw**2
          hama_c = 0.434_WP + 5.975_WP*Zw + 2.888_WP*Zw**2
	  
!*********************************
! compute weighted sediment transport qs
!*********************************

          b = hama_v(n)*dmean*con(n)
          qsu(n) = b*U2D_node(n)
          qsv(n) = b*V2D_node(n)
	  
!*********************************
! akkumulation sediment particles
!*********************************

          a_min=0.8_WP*dmean  ! for test *10
          ! a_min=0.005_WP*dmean
          c_e = 0.015_WP*d_m*Tsp(n)**1.5_WP/(D_st**0.3_WP*a_min*dmean)

          E_sed(n) = hama_c*w_s1*(c_e - con(n))
    !      if (E_sed(n) > 0.0_WP) E_sed(n)=min(0.04_WP,E_sed(n))
    !      if (E_sed(n) < 0.0_WP) E_sed(n)=max(-0.04_WP,E_sed(n))

   enddo

!$OMP END DO
  ! write(*,*) 'sediment module  = 5 '

 !    write(*,*) 'akkumulation of sediment part max, min: E_sed = ' , maxval(E_sed),minval(E_sed)
 !    write(*,*) 'sediment transport max: qsu,qsv = ' , maxval(qsu),minval(qsv)
 !    write(*,*) 'hama =                                         ', maxval(hama_v),maxval(test1)
 !   write(*,*) 'c_e, max, min = ', maxval(test4),minval(test4)
 !   write(*,*) 'transport stage, max, min =', maxval(Tsp),minval(Tsp)
 !    write(*,*) 'dynamical velocity, max, min = ', maxval(test2),minval(test2)
 !    write(*,*) 'critical dynamical velocity, max, min = ', maxval(test3),minval(test3)
 !    write(*,*) '.............................'
!*********************************
! compute bed-load transport qb
!*********************************
          a = 0.053_WP*dsqrt(g*(plop_s/density_0 - 1.0_WP))*d_m**1.5_WP/D_st**0.3_WP

!write(*,*) 'sediment module  = 6 '
!$OMP DO
    do elem=1,elem2D
         elnodes=elem2D_nodes(:,elem)
          Tsp_el=sum(w_cv(1:4,elem)*Tsp(elnodes))
!          D_st_el=sum(w_cv(1:4,elem)*D_st(elnodes))
!          a = 0.053_WP*dsqrt(g*(plop_s/density_0 - 1.0_WP))*d_m**1.5_WP/D_st_el**0.3_WP
          uv = dsqrt(U_n_2D(1,elem)**2 + U_n_2D(2,elem)**2)
!aa050719	  if (uv < 1.e-12) then
	  if (uv < 1.e-8) then
	   qbu(elem) = 0.0_WP
	   qbv(elem) = 0.0_WP
         else
          qbu(elem) = a*Tsp_el**2.1_WP*U_n_2D(1,elem)/uv
          qbv(elem) = a*Tsp_el**2.1_WP*U_n_2D(2,elem)/uv
         endif
     enddo
!$OMP END DO
!$OMP END PARALLEL  

!*********************************
! compute equation for concentration
!*********************************
!write(*,*) 'sediment module  = '
       call concentration_sed
 !  write(*,*) 'concentration of sedimentation max, min: con = ' , maxval(con),minval(con)
 ! if (n_dt>3010) then
!write(*,*) 'sediment module  = 7 '
         call bottom_evolution
 ! endif 
!  write(*,*) 'bottom variation max, min: h_var = ' , maxval(h_var),minval(h_var)
 
     deallocate(U2D_node,V2D_node)
     deallocate(Tsp)
     deallocate(qsu,qsv)

END subroutine sediment
!============================================================================
SUBROUTINE concentration_sed
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

 integer      :: el(2), enodes(2), n, edge, ed, nodes(2)
 
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2
 real(kind=WP) :: Cmean, fD, fDold, x1, y1, dmean, un1 !, relax2clim_ac
 
 real(kind=WP), allocatable :: cHrhs(:), Vel_nor(:)

 allocate(cHrhs(nod2D))
 allocate(Vel_nor(nod2D))

! =========================
! momentum equation for
! concentration in conservative form
!!!!!!!!!!!!!!! cr_distr2_2D -  for boundary with sediment transport!!!!!!!!!!
! =========================
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,edge,ed,fDold,c1,c2,deltaX1,deltaX2,deltaY1,deltaY2,Cmean,nodes,el,x1,y1,dmean,un1,fD)
!$OMP DO

!!!!!!!!!!!  DO n=1, nod2D
!!!!!!!!!!!     fDold=max(1.0_WP,depth(n) + eta_n_1(n))
!aa67    if (cr_distr2_2D(n) .eq. 0.0_WP) then
!!!!!!!!!!!!!     con(n) = con(n)*fDold
!aa67    else 
!aa67     con(n) = cr_distr2_2D(n)*fDold
!aa67    endif
!!!!!!!!!  END DO
!!!!!!!!!!!!enddo



  DO n=1, nod2D
!aa      if (cr_distr2_2D(n) .eq. 0.0_WP) then
     fDold=max(Dmin,depth(n) + eta_n_1(n))
     con(n) = con(n)*fDold
!aa     endif
  END DO


!$OMP END DO
! =================
! Clean the rhs
! =================          
 cHrhs=0.0_WP
! =================
! Horizontal advection
! =================
!   c1 = 0.0_WP
!   c2 = 0.0_WP

!$OMP DO

  DO edge=1, edge2D 
   enodes=edge_nodes(:,edge)   
   el=edge_tri(:,edge)
   c1 = 0.0_WP
   c2 = 0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
!========== ============
! First segment
! Linear upwind reconstruction
! ======================
   if(U_n_2D(2,el(1))*deltaX1- U_n_2D(1,el(1))*deltaY1>0.0_WP) then   
      Cmean=con(enodes(2))*hama_v(enodes(2))
   else
      Cmean=con(enodes(1))*hama_v(enodes(1))
   endif
   c1=U_n_2D(2,el(1))*Cmean*deltaX1 - U_n_2D(1,el(1))*Cmean*deltaY1
   cHrhs(enodes(1)) = cHrhs(enodes(1)) + c1
   cHrhs(enodes(2)) = cHrhs(enodes(2)) - c1

! ======================
! Second segment
! Linear upwind reconstruction
!========== ============
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   if(U_n_2D(2,el(2))*deltaX2- U_n_2D(1,el(2))*deltaY2<0.0_WP) then   
      Cmean= con(enodes(2))*hama_v(enodes(2))
   else
      Cmean= con(enodes(1))*hama_v(enodes(1))
   endif
   c2=-U_n_2D(2,el(2))*Cmean*deltaX2 + U_n_2D(1,el(2))*Cmean*deltaY2
   cHrhs(enodes(1))=cHrhs(enodes(1)) + c2
   cHrhs(enodes(2))=cHrhs(enodes(2)) - c2

  end if
   
 END DO
!$OMP END DO 
!$OMP DO
  DO n=1,nod2D
     cHrhs(n) = cHrhs(n)/area(n) !+ E_sed(n)
  END DO
!$OMP END DO

!=================
! boundary conditions
! =================
Vel_nor = 0.0_WP

!$OMP DO
DO ed=1+edge2D_in, edge2D 
   nodes=edge_nodes(:,ed)   
   el=edge_tri(:,ed)
   x1=edge_cross_dxdy(1,ed)
   y1=edge_cross_dxdy(2,ed)
   dmean=sum(eta_n(nodes)+depth(nodes))/2.0_WP
   dmean=max(Dmin,dmean)
   un1=(U_n_2D(2,el(1))*x1- U_n_2D(1,el(1))*y1)*dmean
 
   Vel_nor(nodes(1))=Vel_nor(nodes(1))+un1/area(nodes(1))
   Vel_nor(nodes(2))=Vel_nor(nodes(2))-un1/area(nodes(2))
 END DO
!$OMP END DO

     relax2clim_ac = 0.0 !000000000001
     clim_relax = 0.0 ! 0000000001
!$OMP DO
  DO n=1, nod2D
     if (index_nod2D(n) == 2) then
       fD=max(Dmin,depth(n) + eta_n_1(n))
       if (Vel_nor(n) > 0.0_WP) then
        cHrhs(n)=cHrhs(n) + relax2clim_ac*(con_bc(n)*fD - con(n))
       else
       cHrhs(n)=clim_relax*(con_bc(n)*fD - con(n))
       endif
      endif
 
  END DO
!$OMP END DO

! =================
! Update ttfold (to be used on the next time level)
! and compute new ttf
! =================          
!$OMP DO

  DO n=1,nod2D
      fD=max(Dmin,depth(n) + eta_n(n))

!aa67      if (cr_distr2_2D(n) .eq. 0.0_WP) then
!     con(n) = ac(n)**2*(con(n) + (cHrhs(n)+E_sed(n))*dt_2D)/fD
     if (ac(n) .lt. 1._WP) then
     con(n) = (1._WP - ac(n))**3.*ac(n)**10.*(con(n) + (cHrhs(n)+E_sed(n))*dt_2D)/fD
     else
     con(n) = (con(n) + (cHrhs(n)+E_sed(n))*dt_2D)/fD
     endif

!aa67      else
!aa67 011220      con(n) = 1000.0_WP*cr_distr2_2D(n)
!aa67      con(n) = cr_distr2_2D(n)
!aa67      endif 
  END DO

!$OMP END DO
!$OMP END PARALLEL

deallocate(cHrhs, Vel_nor)
end subroutine concentration_sed
!=============================================
   subroutine bottom_evolution

  use o_MESH
  use o_ARRAYS
  use o_PARAM
  implicit none

  real(kind=WP) :: c1, c2, a
  integer             :: ed, el(2), enodes(2), n

  real(kind=WP), allocatable :: cHrhs(:)

  allocate(cHrhs(nod2D))
  cHrhs = 0.0_WP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ed,enodes,el,c1,c2,n)
!$OMP DO

 ! ==============
 ! internal edges
 ! ==============
 DO ed=1,edge2D_in
    enodes=edge_nodes(:,ed)
    el=edge_tri(:,ed)
    c1=-qbv(el(1))*edge_cross_dxdy(1,ed)+qbu(el(1))*edge_cross_dxdy(2,ed)
    c2=qbv(el(2))*edge_cross_dxdy(3,ed)-qbu(el(2))*edge_cross_dxdy(4,ed)
    c1 = c1 + c2
    cHrhs(enodes(1)) = cHrhs(enodes(1)) + c1
    cHrhs(enodes(2)) = cHrhs(enodes(2)) - c1
 END DO
!$OMP END DO
 ! =============
 ! boundary edges
 ! only the left element (1) is available
 ! =============
!$OMP DO
  DO ed=1+edge2D_in, edge2D
    enodes=edge_nodes(:,ed)
    el=edge_tri(:,ed)
    c1=-qbv(el(1))*edge_cross_dxdy(1,ed)+qbu(el(1))*edge_cross_dxdy(2,ed)
    cHrhs(enodes(1)) = cHrhs(enodes(1)) + c1
    cHrhs(enodes(2)) = cHrhs(enodes(2)) - c1
 END DO
!$OMP END DO
   a = dt_2d/(1.0_WP - e_p)
!$OMP DO
   do n=1,nod2D
     if (ac(n) .lt. 1._WP) then
    h_var(n)=h_var(n) + (1.0-ac(n))**3.*ac(n)**10*a*(cHrhs(n)/area(n) + E_sed(n))
!     else
!      if (ac_sed(n) .lt. 1.0_WP) then
!    h_var(n)=h_var(n) + (1.0-ac_sed(n))**3.*ac_sed(n)**10*a*(cHrhs(n)/area(n) + E_sed(n))
      else
    h_var(n)=h_var(n) + a*(cHrhs(n)/area(n) + E_sed(n))
     endif
!     endif

!    h_var(n)=h_var(n) + ac(n)*a*(cHrhs(n)/area(n) + E_sed(n))
!    h_var(n)=h_var(n) + ac(n)**10*a*(cHrhs(n)/area(n) + E_sed(n))
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  deallocate(cHrhs)
  end subroutine bottom_evolution
