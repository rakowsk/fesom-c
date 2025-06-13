
!===================================================================================
! It is assumed that velocity is at n+1/2, hence only tracer field
! is AB2 interpolated to n+1/2.
!
! Alexey Androsov
! 12.01.15
!
SUBROUTINE solve_tracer_muscl(ttf, ttfold, tracer_type)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
!a USE g_PARSUP
IMPLICIT NONE

 integer      :: el(2), enodes(2), n, nz, edge, nodes(2), ed, i
 integer      :: nl1, nl2, nn, nod1, nod2, el1
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, flux=0.0 , delta_sigma
 real(kind=WP) :: tvert(nsigma), a, b, c, da, db, dg, aux
 real(kind=WP) :: x1,y1, un1, dmean,un2,num_ord
 real(kind=WP) :: Kh, Kvv, Tx, Ty, Tmean, Tmean1, Tmean2, rdata=0.0
 real(kind=WP) :: Vel_nor(nod2D)

 real(kind=WP) :: ttf(nsigma-1, nod2D)
 real(kind=WP) :: ttfold(nsigma-1, nod2D)
 real(kind=WP) :: ttrhs(nsigma-1, nod2D), ttrhs2(nsigma-1, nod2D),ttrhs_c(nsigma-1,nod2D)
 real(kind=WP) :: ttx(nsigma-1,elem2D), tty(nsigma-1,elem2D)

character*1  :: tracer_type


num_ord=0.85_WP
! =================
! AB interpolation
! =================

 ttfold=-(0.5_WP+epsilon)*ttfold+(1.5_WP+epsilon)*ttf

  if(tracer_type=='t') then
   ttfold = ttfold - T_aver
   ttf = ttf - T_aver
  else if (tracer_type=='s') then
   ttfold = ttfold - S_aver
   ttf = ttf - S_aver
  endif

  ! ttfold now contains AB-interpolated tracer! (n+1/2)

! =================
! Compute gradients
! =================
  call tracer_gradient_elements(ttfold, ttx, tty)
  call fill_up_dn_grad(ttx,tty)

  ! elemental gradients are in ttx, tty
  ! They are based on AB-interpolated values.
! =================
! Clean the rhs
! =================
 ttrhs=0.0_WP
 ttrhs2=0.0_WP
 ttrhs_c=0.0_WP
 tvert=0.0_WP

! =================
! Horizontal advection and diffusion
! =================
  DO edge=1, edge2D
   enodes=edge_nodes(:,edge)
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   Kh=elem_area(el(1))

   a=r_earth*elem_cos(el(1))
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)

   b=r_earth*elem_cos(el(2))
   Kh=0.5_WP*(Kh+elem_area(el(2)))
   end if

   Kh=K_hor*Kh/scale_area
   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! Average diffusive flux
   ! ============
   if (el(2)>0) then
   Tx=0.5_WP*Kh*(ttx(nz,el(1))+ttx(nz,el(2)))
   Ty=0.5_WP*Kh*(tty(nz,el(1))+tty(nz,el(2)))
   else
   Tx=Kh*ttx(nz,el(1))
   Ty=Kh*tty(nz,el(1))
   end if
   !Tx=0.0_WP
   !Ty=0.0_WP
   ! ============
   ! MUSCL type reconstruction
   ! ============
      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8

      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8

   if((V_n_filt(nz,el(1))*deltaX1- U_n_filt(nz,el(1))*deltaY1)>0) then
      Tmean=0.15_WP*Tmean2+0.85_WP*(Tmean1+Tmean2)/2.0_WP
   else
      Tmean=0.15_WP*Tmean1+0.85_WP*(Tmean1+Tmean2)/2.0_WP
   end if
   c1=(Je(nz,el(1))*V_n_filt(nz,el(1))*Tmean-Ty)*deltaX1- &
        (Je(nz,el(1))*U_n_filt(nz,el(1))*Tmean-Tx)*deltaY1
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1)) + c1
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2)) - c1

   end do

   ! ============
   ! Second segment
   ! ============
   if(el(2)>0) then
   DO nz=1, nsigma-1
   Tx=0.5_WP*Kh*(ttx(nz,el(1))+ttx(nz,el(2)))
   Ty=0.5_WP*Kh*(tty(nz,el(1))+tty(nz,el(2)))
   Tx=0.0_WP
   Ty=0.0_WP
   ! ============
   ! Linear upwind reconstruction
   ! ============

      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*b*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8

      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*b*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8

   if(V_n_filt(nz,el(2))*deltaX2- U_n_filt(nz,el(2))*deltaY2<0) then
      Tmean=0.15_WP*Tmean2+0.85_WP*(Tmean1+Tmean2)/2.0_WP
   else
      Tmean=0.15_WP*Tmean1+0.85_WP*(Tmean1+Tmean2)/2.0_WP
   end if

   c2=-(Je(nz,el(2))*V_n(nz,el(2))*Tmean-Ty)*deltaX2+ &
          (Je(nz,el(2))*U_n(nz,el(2))*Tmean-Tx)*deltaY2
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1))+c2
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2))-c2

 END DO
  end if
enddo




 DO edge=1, edge2D
   enodes=edge_nodes(:,edge)   
   el=edge_tri(:,edge)
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   a=r_earth*elem_cos(el(1))

   if(el(2)>0) then
     deltaX2=edge_cross_dxdy(3,edge)
     deltaY2=edge_cross_dxdy(4,edge)
     a=0.5_WP*(a+r_earth*elem_cos(el(2)))
   end if     

   DO nz=1, nsigma-1

 ! ============
 ! MUSCL-type reconstruction (high order) and upwind (low order)
 ! ============

      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP   
   
      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP    
   
   un1= Je(nz,el(1))*(-V_n_filt(nz,el(1))*deltaX1+ U_n_filt(nz,el(1))*deltaY1)   ! normal outer velocity
   un2= Je(nz,el(2))*( V_n_filt(nz,el(2))*deltaX2- U_n_filt(nz,el(2))*deltaY2)   ! from 1 to 2
   un1=un1+un2

   c1=(un1+abs(un1))*Tmean1+(un1-abs(un1))*Tmean2
   c1=-0.5_WP*((1.0_WP-num_ord)*c1+un1*num_ord*(Tmean1+Tmean2))


   ttrhs2(nz,enodes(1))=ttrhs2(nz,enodes(1)) + c1
   ttrhs2(nz,enodes(2))=ttrhs2(nz,enodes(2)) - c1

   END DO
 END DO





write(*,*) 'ttrhs2, ttrhs1 diff', maxval(ttrhs2-ttrhs), minval(ttrhs2-ttrhs)















 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO edge=1, edge2D
   enodes=edge_nodes(:,edge)
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)

   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! Average diffusive flux
   ! ============
  c1=Je(nz,el(1))*V_n_filt(nz,el(1))*deltaX1- &
       Je(nz,el(1))*U_n_filt(nz,el(1))*deltaY1
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1)) + c1
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2)) - c1
   END DO
   ! ============
   ! Second segment
   ! ============
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   DO nz=1, nsigma-1
   c2=-Je(nz,el(2))*V_n_filt(nz,el(2))*deltaX2+ &
          Je(nz,el(2))*U_n_filt(nz,el(2))*deltaY2
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1))+c2
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2))-c2
   END DO
   end if

 END DO



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ===================
! Vertical advection and diffusion
! ===================
 DO n=1, nod2D
!a  if(tracer_type=='t') then
!a  flux=-heat_flux(n)/density_0/4200.0_WP
!a  rdata=surf_relax_T*(Tsurf(n)-ttf(1,n))
!a  else if (tracer_type=='s') then
!a  flux=water_flux(n)*SF(1,n)
!a  rdata=surf_relax_S*(Ssurf(n)-ttf(1,n))
!a  end if
  flux=0.0_WP
  rdata=0.0_WP
  ! ===========
  ! Fluxes in the column
  ! ===========
  ! Surface forcing
          if (i_vert_diff) then
	  tvert(1)= -Wvel(1,n)*ttfold(1,n)
	  Kvv=0.0_WP
	  else
	  tvert(1)= (flux + rdata-Wvel(1,n)*ttfold(1,n))
	  Kvv=K_ver
          end if
  ! Bottom conditions
	  tvert(nsigma)=0.0_WP

DO nz=2, nsigma-1

      ! ============
      ! QUICK upwind (3rd order)
      ! ============
      if(Wvel(nz,n)>0.0_WP) then
        if(nz==nsigma-1) then
	  Tmean= ttfold(nz,n)  ! or replace this with
	                                             ! the first order
						     ! upwind  tttfold(nz,n)
	else
	a=zbar(nz,n) - Z(nz-1,n)
	b=Z(nz,n) - zbar(nz,n)
	c=Z(nz+1,n) - zbar(nz,n)
	dg=a*b/(c+a)/(b-c)
	db=-a*c/(b+a)/(b-c)
	da=1.0_WP-dg-db
	Tmean=ttfold(nz-1,n)*da+ttfold(nz,n)*db+ttfold(nz+1,n)*dg
	end if
      end if

      if(Wvel(nz,n)<0.0_WP) then
        if(nz==2) then
	  Tmean=ttfold(nz-1,n)        ! or ttfold(nz-1,n)
	else
	a=Z(nz,n) - zbar(nz,n)
	b=zbar(nz,n) - Z(nz-1,n)
	c=zbar(nz,n) - Z(nz-2,n)
	dg=a*b/(c+a)/(b-c)
	db=-a*c/(b+a)/(b-c)
	da=1.0_WP-dg-db
	Tmean=ttfold(nz,n)*da+ttfold(nz-1,n)*db+ttfold(nz-2,n)*dg
	end if
      end if
         tvert(nz)= -Tmean*Wvel(nz,n)! + &
!	      Kvv*(ttf(nz-1,n)-ttf(nz,n))/(Z(nz-1,n)-Z(nz,n)))
   END DO

   DO nz=1,nsigma-1
       a = ttfold(nz,n)
       ttrhs(nz,n)=( ttrhs(nz,n)/(area(n)*Jc(nz,n))- &
                       a*ttrhs_c(nz,n)/(area(n)*Jc(nz,n)) +&
		       a*(Wvel(nz,n) - Wvel(nz+1,n))/Jc(nz,n) + &
                    (tvert(nz)-tvert(nz+1))/Jc(nz,n) )
   END DO
 END DO

if (cl_relax) then

  do n=1,nod2D
     Vel_nor(n) = 0.0_WP
  enddo

  DO ed = 1+edge2D_in, edge2D
     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1  = edge_tri(1,ed)

     dmean = max(Dmin, 0.5_WP*(eta_n(nod1)+depth(nod1) + eta_n(nod2)+depth(nod2)))

     un1 = (V_filt_2D(el1)*edge_cross_dxdy(1,ed)- U_filt_2D(el1)*edge_cross_dxdy(2,ed))*dmean

     Vel_nor(nod1) = Vel_nor(nod1) + un1/area(nod1)
     Vel_nor(nod2) = Vel_nor(nod2) - un1/area(nod2)
  END DO

 if(tracer_type=='t') then
  DO n=1, nod2D
  !   relax2clim_ac = clim_relax*2.0_WP   !/24.0_WP/10.0_WP
     if (Vel_nor(n) > 0.0_WP) then
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=ttrhs(nz,n) +relax2clim_ac*(Tclim(nz,n)-(ttf(nz,n)+T_aver))
     END DO
     else
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=clim_relax*(Tclim(nz,n)-(ttf(nz,n)+T_aver))
    END DO
     endif
  END DO

 end if ! tracer_type = 't'

  if(tracer_type=='s') then

  DO n=1, nod2D
!     relax2clim_ac = clim_relax*2.0_WP     !/24.0_WP/10._WP
    if (Vel_nor(n) > 0.0_WP) then
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=ttrhs(nz,n) +relax2clim_ac*(Sclim(nz,n)-(ttf(nz,n)+S_aver))
     END DO
     else
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=clim_relax*(Sclim(nz,n)-(ttf(nz,n)+S_aver))
     END DO
     endif
  END DO

  end if ! tracer_type = 's'

endif ! climatology

! =================
! Update ttfold (to be used on the next time level)
! and compute new ttf
! =================


      ttfold=ttf
      !if (index_nod2D(n) /= 2)
      ttf=ttf+ttrhs*dt


   if(tracer_type=='t') then
   ttfold = ttfold + T_aver
   ttf = ttf + T_aver
  endif

   if (tracer_type=='s') then
   ttfold = ttfold + S_aver
   ttf = ttf + S_aver
   endif

!  ===============================================
! relaxation to climat near solid boundary
!===============================================
  if (cl_solid_boder_relax) then

   if (tracer_type=='t') then
  DO n=1, nod2D
     DO nz=1,nsigma-1
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Tclim(nz,n)-ttf(nz,n))
     END DO
  END DO
  end if
  if (tracer_type=='s') then
  DO n=1, nod2D
     DO nz=1,nsigma-1
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Sclim(nz,n)-ttf(nz,n))
     END DO
  END DO
  end if

  end if

end subroutine solve_tracer_muscl

! ==========================================================================
SUBROUTINE fct_init
USE o_MESH
USE o_ARRAYS
USE o_PARAM
!USE g_PARSUP
IMPLICIT NONE
integer      :: my_size

my_size=  nod2D  !myDim_nod2D+eDim_nod2D
allocate(fct_aec(nsigma-1,edge2D))   ! Antidiffusive edge contributions
allocate(fct_LO(nsigma-1, my_size))  ! Low-order solution 
allocate(fct_aec_ver(nsigma, nod2D)) ! antidiffusive 
                                     ! vertical fluxes
allocate(fct_ttf_max(nsigma-1, my_size),fct_ttf_min(nsigma-1, my_size))
allocate(fct_plus(nsigma-1, my_size),fct_minus(nsigma-1, my_size))
! Initialize with zeros: 
 fct_aec=0.0_WP
 fct_LO=0.0_WP
 
 fct_aec_ver=0.0_WP
 fct_ttf_max=0.0_WP
 fct_ttf_min=0.0_WP
 fct_plus=0.0_WP
 fct_minus=0.0_WP
 
 write(*,*) 'FCT is initialized'
 
END SUBROUTINE fct_init
!===========================================================================

SUBROUTINE fct(ttf, ttfold, tracer_type)
!
! 3D Flux Corrected Transport scheme
! Limits antidiffusive fluxes==the difference in flux HO-LO
! LO ==Low-order  (first-order upwind)
! HO ==High-order (3rd/4th order gradient reconstruction method)
! Adds limited fluxes to the LO solution   
USE o_MESH
USE o_ARRAYS
USE o_PARAM
!USE g_PARSUP

IMPLICIT NONE
integer       :: n, nz, k, elem, enodes(3), num, el(2), nl1, nl2, edge
integer       :: nn, nod1, nod2, el1, un1
real(kind=WP) :: flux, ae,tvert_max(nsigma-1),tvert_min(nsigma-1), dmean 
real(kind=WP) :: ttf(nsigma-1, nod2D),ttfold(nsigma-1, nod2D),urhs(nsigma-1, elem2D),vrhs(nsigma-1, elem2D)
real(kind=WP) :: dttf(nsigma-1, nod2D)
real*8        :: flux_eps=1e-16
integer       :: vlimit=1
character*1   :: tracer_type
real(kind=WP) :: Vel_nor(nod2D)

! ----------------------------------
! ttf is the tracer field on step n
! dttf is the increment 
! vlimit sets the version of limiting, see below
! ----------------------------------

dttf=0.0_WP
fct_ttf_max=0.0_WP
fct_ttf_min=0.0_WP
fct_plus=0.0_WP
fct_minus=0.0_WP

! ==========
! a1. max, min between old solution and updated low-order solution 
! ==========
DO n=1,nod2D 
   DO nz=1, nsigma-1 
      fct_ttf_max(nz,n)=max(fct_LO(nz,n), ttf(nz,n))
      fct_ttf_min(nz,n)=min(fct_LO(nz,n), ttf(nz,n))
   END DO
END DO       

!write(*,*) 'max,min fct_ttf_max,min',  maxval(fct_ttf_max), minval(fct_ttf_max), maxval(fct_ttf_min), minval(fct_ttf_min)
! ==========
! a2. Admissible increments on elements
!     (only layers below the first and above the last layer)
! ==========

DO elem=1, elem2D
   enodes=elem2D_nodes(:,elem)
   DO nz=1, nsigma-1
      urhs(nz,elem)=maxval(fct_ttf_max(nz,enodes))
      vrhs(nz,elem)=minval(fct_ttf_min(nz,enodes))
   END DO
END DO
!write(*,*) 'max,min urhs,vrhs',  maxval(urhs), minval(urhs), maxval(vrhs), minval(vrhs)
! ==========
! a3. Bounds on clusters and admissible increments
! ==========

if(vlimit==1) then
!Horizontal
DO n=1, nod2D
   DO nz=1,nsigma-1
      tvert_max(nz)= maxval(urhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
      tvert_min(nz)= minval(vrhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
   END DO
   !fct_ttf_max(:,n)=tvert_max
   !fct_ttf_min(:,n)=tvert_min 
   fct_ttf_max(1,n)=tvert_max(1)-fct_LO(1,n)
   fct_ttf_min(1,n)=tvert_min(1)-fct_LO(1,n)
   DO nz=2,nsigma-2  
   fct_ttf_max(nz,n)=maxval(tvert_max(nz-1:nz+1))-fct_LO(nz,n)
   fct_ttf_min(nz,n)=minval(tvert_min(nz-1:nz+1))-fct_LO(nz,n)
   END DO
   nz=nsigma-1
   fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
   fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)
END DO
!Vertical1: In this version we look at the bounds on the clusters
!           above and below, which leaves wide bounds because typically 
!           vertical gradients are larger.  
end if

if(vlimit==2) then
DO n=1, nod2D
   DO nz=1,nsigma-1
      tvert_max(nz)= maxval(urhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
      tvert_min(nz)= minval(vrhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
   END DO
   DO nz=2, nsigma-2
      tvert_max(nz)=max(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
      tvert_min(nz)=min(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
   END DO
   DO nz=1,nsigma-1
       fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
       fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
   END DO
END DO
!Vertical2: Similar to the version above, but the vertical bounds are more local  
END IF

if(vlimit==3) then
DO n=1, nod2D
   DO nz=1,nsigma-1
      tvert_max(nz)= maxval(urhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
      tvert_min(nz)= minval(vrhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
   END DO
   DO nz=2, nsigma-2
      tvert_max(nz)=min(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
      tvert_min(nz)=max(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
   END DO
   DO nz=1,nsigma-1
       fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
       fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
   END DO
END DO
!Vertical3: Vertical bounds are taken into account only if they are narrower than the
!           horizontal ones  
END IF

!write(*,*) 'max,min fct_ttf_max,min',  maxval(fct_ttf_max), minval(fct_ttf_max), maxval(fct_ttf_min), minval(fct_ttf_min)
! ==========
! b1. Split positive and negative antidiffusive contributions
! ==========


!Vertical
DO n=1, nod2D
   DO nz=1,nsigma-1
      fct_plus(nz,n)=fct_plus(nz,n)+ &
                    (max(0.0_WP,fct_aec_ver(nz,n))+max(0.0_WP,-fct_aec_ver(nz+1,n)))
      fct_minus(nz,n)=fct_minus(nz,n)+ &
                     (min(0.0_WP,fct_aec_ver(nz,n))+min(0.0_WP,-fct_aec_ver(nz+1,n)))  
   END DO
END DO

!write(*,*) 'max,min fct_minus, fct_plus -0000',  maxval(fct_plus), minval(fct_plus),  maxval(fct_minus), minval(fct_minus)

!Horizontal
DO edge=1, edge2D
   enodes(1:2)=edge_nodes(:,edge)   
   el=edge_tri(:,edge)
 
   DO nz=1, nsigma-1
      fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0, fct_aec(nz,edge))
      fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0, fct_aec(nz,edge))
      fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0, -fct_aec(nz,edge))
      fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0, -fct_aec(nz,edge))
   END DO   
END DO 
! ==========
! b2. Limiting factors
! ==========

!write(*,*) 'max,min fct_minus, fct_plus 0000',  maxval(fct_plus), minval(fct_plus),  maxval(fct_minus), minval(fct_minus)
!write(*,*) 'max,min fct_minus, fct_plus 0000',  maxval(Jc), minval(Jc),  maxval(fct_aec), minval(fct_aec)


DO n=1,nod2D
   DO nz=1,nsigma-1
     flux=fct_plus(nz,n)*dt/(area(n)*Jc(nz,n))+flux_eps
     fct_plus(nz,n)=min(1.0_WP,fct_ttf_max(nz,n)/flux)
     flux=fct_minus(nz,n)*dt/(area(n)*Jc(nz,n))-flux_eps
     fct_minus(nz,n)=min(1.0_WP,fct_ttf_min(nz,n)/flux)
   END DO
END DO 

!write(*,*) 'max,min fct_minus, fct_plus',  maxval(fct_plus), minval(fct_plus),  maxval(fct_minus), minval(fct_minus)
	   
! 
!========================	 
! b3. Limiting   
!========================	 
!Vertical
DO n=1, nod2D
   nz=1
   ae=1.0_WP
   flux=fct_aec_ver(nz,n)
   if(flux>=0.0) then 
      ae=min(ae,fct_plus(nz,n))
   else
      ae=min(ae,fct_minus(nz,n))
   end if
   fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)    
   DO nz=2,nsigma-1
      ae=1.0
      flux=fct_aec_ver(nz,n)
        if(flux>=0.0_WP) then 
	   ae=min(ae,fct_minus(nz-1,n))
	   ae=min(ae,fct_plus(nz,n))
	else
	   ae=min(ae,fct_plus(nz-1,n))
	   ae=min(ae,fct_minus(nz,n))
	end if   
      fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)
   END DO
   ! the bottom flux is always zero 
END DO

!write(*,*) 'max,min vert  ae 1',  maxval(fct_aec_ver), minval(fct_aec_ver)

!Horizontal
DO edge=1, edge2D
   enodes(1:2)=edge_nodes(:,edge)
   el=edge_tri(:,edge)
  
   DO nz=1, nsigma-1
      ae=1.0
      flux=fct_aec(nz,edge)
      if(flux>=0.) then
          ae=min(ae,fct_plus(nz,enodes(1)))
          ae=min(ae,fct_minus(nz,enodes(2)))
      else
          ae=min(ae,fct_minus(nz,enodes(1)))
          ae=min(ae,fct_plus(nz,enodes(2)))
      endif
      fct_aec(nz,edge)=ae*fct_aec(nz,edge)
   END DO
END DO

!write(*,*) 'max,min vert  ae2',  maxval(fct_aec), minval(fct_aec)
!==========================
! c. Update the solution
!==========================
      ttfold=ttf
!write(*,*) 'max,min ttf', maxval(ttf), minval(ttf)

   ! Vertical part
DO n=1, nod2D
   DO nz=1,nsigma-1  
   ttf(nz,n)=fct_LO(nz,n)+ dt* &
             (fct_aec_ver(nz,n)-fct_aec_ver(nz+1,n))/(Jc(nz,n)*area(n))    
   END DO
END DO
!write(*,*) 'max,min vert  fct_aec_ver, ttf',  maxval(fct_aec_ver), minval(fct_aec_ver), maxval(ttf), minval(ttf)

   ! Horizontal part
DO edge=1, edge2D
   enodes(1:2)=edge_nodes(:,edge)
   el=edge_tri(:,edge)
  DO nz=1, nsigma-1
      ttf(nz,enodes(1))=ttf(nz,enodes(1))+fct_aec(nz,edge)*dt/(area(enodes(1))*Jc(nz,enodes(1)))
      ttf(nz,enodes(2))=ttf(nz,enodes(2))-fct_aec(nz,edge)*dt/(area(enodes(2))*Jc(nz,enodes(2)))
   END DO
END DO  

!write(*,*) 'max,min vert  fct_aec_ver',  maxval(fct_aec), minval(fct_aec), maxval(ttf), minval(ttf)


!===============================================
! relaxation to climat near solid boundary
!===============================================
  if (cl_solid_boder_relax) then

   if (tracer_type=='t') then
  DO n=1, nod2D
     DO nz=1,nsigma-1
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Tclim(nz,n)-ttf(nz,n))
     END DO
  END DO
  end if
  if (tracer_type=='s') then
  DO n=1, nod2D
     DO nz=1,nsigma-1
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Sclim(nz,n)-ttf(nz,n))
     END DO
  END DO
  end if

  end if


END SUBROUTINE fct


!===========================================================================
SUBROUTINE fct_muscl_2(ttf, ttfold, num_ord)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE
 integer       :: el(2), enodes(2), n, nz, edge
 integer       :: n2
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2,Kvv, rdata
 real(kind=WP) :: tvert(nsigma),tvert2(nsigma), a, qc, qu, qd, un1, un2, Tupw1,dg,db,da,flux,b,c
 real(kind=WP) :: Tmean1, Tmean2, Tmean, num_ord
 real(kind=WP) :: ttf(nsigma-1,nod2D), ttfold(nsigma-1,nod2D)
 real(kind=WP) :: ttrhs_c(nsigma-1,nod2D),ttrhs(nsigma-1,nod2D), fct_HO(nsigma-1,nod2D)

character*1  :: tracer_type

 ! ----------------------------------------
 ! It is assumed that velocity is at n+1/2, hence only tracer field 
 ! is AB2 interpolated to n+1/2. 
 ! ttf contains tracer at step n
 ! ttfold is AB interpolated (n+1/2)
 ! num_ord is the fraction of fourth-order contribution in the HO solution
 ! The result is the low-order solution
 ! and vertical and horizontal antidiffusive fluxes
 ! They are put into fct_LO
 !                   fct_aec
 !                   fct_aec_ver
 ! ----------------------------------------


! =================
! AB interpolation
! =================

 ttfold=-(0.5_WP+epsilon)*ttfold+(1.5_WP+epsilon)*ttf


  ! ttfold now contains AB-interpolated tracer! (n+1/2)

! =================
! Clean low-order solution
! =================

 tvert=0.0_WP
 tvert2=0.0_WP
 fct_LO=0.0_WP 
 fct_aec=0.0_WP  
 fct_aec_ver=0.0_WP 
 ttrhs_c=0.0_WP
 ttrhs=0.0_WP
 fct_HO=0.0_WP
! =================
 ! Horizontal advection
 ! =================
  DO edge=1, edge2D
   enodes=edge_nodes(:,edge)   
   el=edge_tri(:,edge)
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   a=r_earth*elem_cos(el(1))

   if(el(2)>0) then
     deltaX2=edge_cross_dxdy(3,edge)
     deltaY2=edge_cross_dxdy(4,edge)
     a=0.5_WP*(a+r_earth*elem_cos(el(2)))
   end if     

   DO nz=1, nsigma-1

 ! ============
 ! MUSCL-type reconstruction (high order) and upwind (low order)
 ! ============

      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP   
   
      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP    
   
   un1= Je(nz,el(1))*(-V_n_filt(nz,el(1))*deltaX1+ U_n_filt(nz,el(1))*deltaY1)   ! normal outer velocity
   un2= Je(nz,el(2))*( V_n_filt(nz,el(2))*deltaX2- U_n_filt(nz,el(2))*deltaY2)   ! from 1 to 2
   un1=un1+un2
   Tupw1=0.5_WP*(ttf(nz, enodes(1))*(un1+abs(un1))+ttf(nz, enodes(2))*(un1-abs(un1)))
   c1=-Tupw1 
   
   fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+c1            ! Low-order upwind 
   fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-c1            !
   c1=(un1+abs(un1))*Tmean1+(un1-abs(un1))*Tmean2
   c1=-0.5_WP*((1.0_WP-num_ord)*c1+un1*num_ord*(Tmean1+Tmean2))


   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1)) + c1
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2)) - c1


   fct_aec(nz,edge)=c1+Tupw1                               ! Antidiffusive edge flux
   END DO
 END DO


  DO edge=1, edge2D
   enodes=edge_nodes(:,edge)
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)

   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! Average diffusive flux
   ! ============
  c1=Je(nz,el(1))*V_n_filt(nz,el(1))*deltaX1- &
       Je(nz,el(1))*U_n_filt(nz,el(1))*deltaY1
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1)) + c1
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2)) - c1
   END DO
   ! ============
   ! Second segment
   ! ============
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   DO nz=1, nsigma-1
   c2=-Je(nz,el(2))*V_n_filt(nz,el(2))*deltaX2+ &
          Je(nz,el(2))*U_n_filt(nz,el(2))*deltaY2
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1))+c2
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2))-c2
   END DO
   end if

 END DO

! ===================
! Vertical advection
! ===================

DO n=1, nod2D

  ! Surface conditions

	  tvert(1)= -Wvel(1,n)*ttfold(1,n)*area(n)
          tvert2(1)= -Wvel(1,n)*ttfold(1,n)*area(n)
  
! Bottom conditions

	  tvert(nsigma)=0.0_WP
          tvert2(nsigma)=0.0_WP

    nz=2  
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1*area(n)
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    nz=nsigma-1
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1*area(n)
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    tvert2(nz)= (-Tmean)*area(n)
    ! the remaining levels    
    DO nz=3,nsigma-2
    ! ------
    ! First-order upwind estimate
    ! ------
    Tupw1=0.5*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
              ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    ! -----
    ! centered (4th order)
    ! -----
    qc=(ttfold(nz-1,n)-ttfold(nz,n))/(Z(nz-1,n)-Z(nz,n))
    qu=(ttfold(nz,n)-ttfold(nz+1,n))/(Z(nz,n)-Z(nz+1,n))    
    qd=(ttfold(nz-2,n)-ttfold(nz-1,n))/(Z(nz-2,n)-Z(nz-1,n))
    Tmean1=ttfold(nz,n)+(2*qc+qu)*(zbar(nz,n)-Z(nz,n))/3.0_WP
    Tmean2=ttfold(nz-1,n)+(2*qc+qd)*(zbar(nz,n)-Z(nz-1,n))/3.0_WP
    Tmean=(Wvel(nz,n)+abs(Wvel(nz,n)))*Tmean1+(Wvel(nz,n)-abs(Wvel(nz,n)))*Tmean2
    tvert(nz)=-Tupw1*area(n)
    fct_aec_ver(nz,n)=(-0.5_WP*(num_ord*(Tmean1+Tmean2)*Wvel(nz,n)+(1.0_WP-num_ord)*Tmean)+Tupw1)*area(n)
    tvert2(nz)=(-0.5_WP*(num_ord*(Tmean1+Tmean2)*Wvel(nz,n)+(1.0_WP-num_ord)*Tmean))*area(n)

    end do
    fct_aec_ver(1,n)=0.0_WP
    fct_aec_ver(nsigma,n)=0.0_WP
   DO nz=1,nsigma-1
       a = ttfold(nz,n)
      fct_LO(nz,n)=ttf(nz,n)+dt*(fct_LO(nz,n) + &
                      (tvert(nz)-tvert(nz+1)))/(Jc(nz,n)*area(n))- &
                       a*ttrhs_c(nz,n)/(area(n)*Jc(nz,n)) +&
		       a*(Wvel(nz,n) - Wvel(nz+1,n))/Jc(nz,n)
      fct_HO(nz,n)=ttf(nz,n)+dt*(ttrhs(nz,n) + &
                      (tvert2(nz)-tvert2(nz+1)))/(Jc(nz,n)*area(n))- &
                       a*ttrhs_c(nz,n)/(area(n)*Jc(nz,n)) +&
		       a*(Wvel(nz,n) - Wvel(nz+1,n))/Jc(nz,n)

   END DO

 ENDDO

! ===================
! Vertical advection and diffusion
! ===================
 DO n=1, nod2D
!a  if(tracer_type=='t') then
!a  flux=-heat_flux(n)/density_0/4200.0_WP
!a  rdata=surf_relax_T*(Tsurf(n)-ttf(1,n))
!a  else if (tracer_type=='s') then
!a  flux=water_flux(n)*SF(1,n)
!a  rdata=surf_relax_S*(Ssurf(n)-ttf(1,n))
!a  end if
  flux=0.0_WP
  rdata=0.0_WP
  ! ===========
  ! Fluxes in the column
  ! ===========
  ! Surface forcing
          if (i_vert_diff) then
	  tvert(1)= -Wvel(1,n)*ttfold(1,n)*area(n)
	  Kvv=0.0_WP
	  else
	  tvert(1)= (flux + rdata-Wvel(1,n)*ttfold(1,n))*area(n)
	  Kvv=K_ver
          end if
  ! Bottom conditions
	  tvert(nsigma)=0.0_WP

DO nz=2, nsigma-1

      ! ============
      ! QUICK upwind (3rd order)
      ! ============
      if(Wvel(nz,n)>0.0_WP) then
        if(nz==nsigma-1) then
	  Tmean= ttfold(nz,n)  ! or replace this with
	                                             ! the first order
						     ! upwind  tttfold(nz,n)
	else
	a=zbar(nz,n) - Z(nz-1,n)
	b=Z(nz,n) - zbar(nz,n)
	c=Z(nz+1,n) - zbar(nz,n)
	dg=a*b/(c+a)/(b-c)
	db=-a*c/(b+a)/(b-c)
	da=1.0_WP-dg-db
	Tmean=ttfold(nz-1,n)*da+ttfold(nz,n)*db+ttfold(nz+1,n)*dg
	end if
      end if

      if(Wvel(nz,n)<0.0_WP) then
        if(nz==2) then
	  Tmean=ttfold(nz-1,n)        ! or ttfold(nz-1,n)
	else
	a=Z(nz,n) - zbar(nz,n)
	b=zbar(nz,n) - Z(nz-1,n)
	c=zbar(nz,n) - Z(nz-2,n)
	dg=a*b/(c+a)/(b-c)
	db=-a*c/(b+a)/(b-c)
	da=1.0_WP-dg-db
	Tmean=ttfold(nz,n)*da+ttfold(nz-1,n)*db+ttfold(nz-2,n)*dg
	end if
      end if
         tvert(nz)= -Tmean*Wvel(nz,n)*area(n)! + &
!	      Kvv*(ttf(nz-1,n)-ttf(nz,n))/(Z(nz-1,n)-Z(nz,n)))
   END DO

   DO nz=1,nsigma-1
       a = ttfold(nz,n)
      ttrhs(nz,n)=ttf(nz,n)+dt*(ttrhs(nz,n) + &
                      (tvert(nz)-tvert(nz+1)))/(Jc(nz,n)*area(n))- &
                       a*ttrhs_c(nz,n)/(area(n)*Jc(nz,n)) +&
		       a*(Wvel(nz,n) - Wvel(nz+1,n))/Jc(nz,n)


   END DO
 END DO



write(*,*) 'max,min tvert', maxval(tvert), minval(tvert)
ttfold=ttf
ttf=ttrhs

 ! Summary:   
 ! fct_LO contains full low-order solution
 ! fct_aec contains antidiffusive component of horizontal flux 
 ! fct_aec_ver contains antidiffusive component of vertical fluxes
end subroutine fct_muscl_2
!===========================================================================



!===========================================================================
SUBROUTINE fct_muscl_LH(ttf, ttfold, num_ord)
!VF 09.01.20  
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE
 integer       :: el(2), enodes(2), n, nz, edge
 integer       :: n2
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2
 real(kind=WP) :: tvert(nsigma), a, qc, qu, qd, un1, un2, Tupw1
 real(kind=WP) :: Tmean1, Tmean2, Tmean, num_ord
 real(kind=WP) :: ttf(nsigma-1,nod2D), ttfold(nsigma-1,nod2D)
 real(kind=WP) :: ttrhs_c(nsigma-1,nod2D), ttrhs_c2(nsigma-1,nod2D)

character*1  :: tracer_type

 ! ----------------------------------------
 ! It is assumed that velocity is at n+1/2, hence only tracer field 
 ! is AB2 interpolated to n+1/2. 
 ! ttf contains tracer at step n
 ! ttfold is AB interpolated (n+1/2)
 ! num_ord is the fraction of fourth-order contribution in the HO solution
 ! The result is the low-order solution
 ! and vertical and horizontal antidiffusive fluxes
 ! They are put into fct_LO
 !                   fct_aec
 !                   fct_aec_ver
 ! ----------------------------------------
  
 ! fct_LO is full low-order solution
 ! fct_aec is antidiffusive horizontal fluxes without limitation
 ! fct_aec_ver is antidiffusive vertical fluxes without limitation


! =================
! AB interpolation
! =================

 ttfold=-(0.5_WP+epsilon)*ttfold+(1.5_WP+epsilon)*ttf


  ! ttfold now contains AB-interpolated tracer! (n+1/2)

! =================
! Clean low-order solution
! =================

 tvert=0.0_WP
 fct_LO=0.0_WP 
 fct_aec=0.0_WP  
 fct_aec_ver=0.0_WP 
 ttrhs_c=0.0_WP

! =================
 ! Horizontal advection
 ! =================
  DO edge=1, edge2D
   enodes=edge_nodes(:,edge)   
   el=edge_tri(:,edge)
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   a=r_earth*elem_cos(el(1))

   if(el(2)>0.0_WP) then
     deltaX2=edge_cross_dxdy(3,edge)
     deltaY2=edge_cross_dxdy(4,edge)
     a=0.5_WP*(a+r_earth*elem_cos(el(2)))
   else
     deltaX2=0.0_WP
     deltaY2=0.0_WP
   end if     

   DO nz=1, nsigma-1

 ! ============
 ! MUSCL-type reconstruction (high order) and upwind (low order)
 ! ============

      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP   
   
      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP    

   
   un1= Je(nz,el(1))*(-V_n_filt(nz,el(1))*deltaX1+ U_n_filt(nz,el(1))*deltaY1)   ! normal outer velocity
   un2= Je(nz,el(2))*( V_n_filt(nz,el(2))*deltaX2- U_n_filt(nz,el(2))*deltaY2)   ! from 1 to 2
   un1=un1+un2

   Tupw1=0.5_WP*(ttf(nz, enodes(1))*(un1+abs(un1))+ttf(nz, enodes(2))*(un1-abs(un1)))
   c1=-Tupw1 
   
   fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+c1            ! Low-order upwind 
   fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-c1            !
   c1=(un1+abs(un1))*Tmean1+(un1-abs(un1))*Tmean2
   c1=-0.5_WP*((1.0_WP-num_ord)*c1+un1*num_ord*(Tmean1+Tmean2))
   fct_aec(nz,edge)=c1+Tupw1                               ! Antidiffusive edge flux
   END DO
 END DO

   ! ============
   ! Average diffusive flux
   ! ============

  DO edge=1, edge2D
   enodes=edge_nodes(:,edge)
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
  c1=Je(nz,el(1))*V_n_filt(nz,el(1))*deltaX1- &
       Je(nz,el(1))*U_n_filt(nz,el(1))*deltaY1
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1)) + c1
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2)) - c1
   END DO
   ! ============
   ! Second segment
   ! ============
   if(el(2)>0.0_WP) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   DO nz=1, nsigma-1
   c2=-Je(nz,el(2))*V_n_filt(nz,el(2))*deltaX2+ &
          Je(nz,el(2))*U_n_filt(nz,el(2))*deltaY2
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1))+c2
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2))-c2
   END DO
   end if

 END DO

! ===================
! Vertical advection
! ===================


DO n=1, nod2D

  ! Surface conditions

	  tvert(1)= -Wvel(1,n)*ttf(1,n)

  ! Bottom conditions
	  tvert(nsigma)=0.0_WP


    nz=2  
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    nz=nsigma-1
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    ! the remaining levels    
    DO nz=3,nsigma-2
    ! ------
    ! First-order upwind estimate
    ! ------
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
              ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    ! -----
    ! centered (4th order)
    ! -----
    qc=(ttfold(nz-1,n)-ttfold(nz,n))/(Z(nz-1,n)-Z(nz,n))
    qu=(ttfold(nz,n)-ttfold(nz+1,n))/(Z(nz,n)-Z(nz+1,n))    
    qd=(ttfold(nz-2,n)-ttfold(nz-1,n))/(Z(nz-2,n)-Z(nz-1,n))
    Tmean1=ttfold(nz,n)+(2.0_WP*qc+qu)*(zbar(nz,n)-Z(nz,n))/3.0_WP
    Tmean2=ttfold(nz-1,n)+(2.0_WP*qc+qd)*(zbar(nz,n)-Z(nz-1,n))/3.0_WP
    Tmean=(Wvel(nz,n)+abs(Wvel(nz,n)))*Tmean1+(Wvel(nz,n)-abs(Wvel(nz,n)))*Tmean2
    tvert(nz)=-Tupw1
    fct_aec_ver(nz,n)=(-0.5_WP*(num_ord*(Tmean1+Tmean2)*Wvel(nz,n)+(1.0_WP-num_ord)*Tmean)+Tupw1)*area(n)
    end do

    fct_aec_ver(1,n)=0.0_WP
    fct_aec_ver(nsigma,n)=0.0_WP

   DO nz=1,nsigma-1
       a = ttf(nz,n)
      fct_LO(nz,n)=ttf(nz,n)+dt*((fct_LO(nz,n)/area(n) + &
                      (tvert(nz)-tvert(nz+1)))/Jc(nz,n) - &
                       a*(ttrhs_c(nz,n)/area(n)- &
		       (Wvel(nz,n) - Wvel(nz+1,n)))/Jc(nz,n))

   END DO

 ENDDO

!write(*,*) 'max,min fct_LO', maxval(fct_LO), minval(fct_LO)
!write(*,*) 'max,min fct_aec', maxval(fct_aec), minval(fct_aec)
!write(*,*) 'max,min fct_aec_ver', maxval(fct_aec_ver), minval(fct_aec_ver)

end subroutine fct_muscl_LH
SUBROUTINE fct_muscl_LH2(ttf, ttfold, num_ord)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE
 integer       :: el(2), enodes(2), n, nz, edge
 integer       :: n2
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, Tx, Ty, b
 real(kind=WP) :: tvert(nsigma), a, qc, qu, qd, un1, un2, Tupw1
 real(kind=WP) :: Tmean1, Tmean2, Tmean, num_ord
 real(kind=WP) :: ttf(nsigma-1,nod2D), ttfold(nsigma-1,nod2D)
 real(kind=WP) :: ttrhs_c(nsigma-1,nod2D), ttrhs_c2(nsigma-1,nod2D)


! =================
! AB interpolation
! =================

 ttfold=-(0.5_WP+epsilon)*ttfold+(1.5_WP+epsilon)*ttf


  ! ttfold now contains AB-interpolated tracer! (n+1/2)

! =================
! Clean low-order solution
! =================

 tvert=0.0_WP
 fct_LO=0.0_WP 
 fct_aec=0.0_WP  
 fct_aec_ver=0.0_WP 
 ttrhs_c=0.0_WP

! =================
! Horizontal advection and diffusion
! =================
  DO edge=1, edge2D
   enodes=edge_nodes(:,edge)
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)

   a=r_earth*elem_cos(el(1))
   if(el(2)>0.0_WP) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   b=r_earth*elem_cos(el(2))
   end if

   Tx=0.0_WP
   Ty=0.0_WP
   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! MUSCL type reconstruction
   ! ============
      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP


   if((V_n_filt(nz,el(1))*deltaX1- U_n_filt(nz,el(1))*deltaY1)>0.0_WP) then
      Tmean=(1.0_WP-num_ord)*Tmean2+num_ord*(Tmean1+Tmean2)/2.0_WP
   else
      Tmean=(1.0_WP-num_ord)*Tmean1+num_ord*(Tmean1+Tmean2)/2.0_WP
   end if
 
   c1=(Je(nz,el(1))*V_n_filt(nz,el(1))*Tmean-Ty)*deltaX1- &
        (Je(nz,el(1))*U_n_filt(nz,el(1))*Tmean-Tx)*deltaY1
		
   un1= Je(nz,el(1))*(-V_n_filt(nz,el(1))*deltaX1+ U_n_filt(nz,el(1))*deltaY1)   ! normal outer velocity
   Tupw1=0.5_WP*(ttf(nz, enodes(1))*(un1+abs(un1))+ttf(nz, enodes(2))*(un1-abs(un1))) 
   
   fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))-Tupw1           ! Low-order upwind 
   fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))+Tupw1   
   fct_aec(nz,edge)=c1+Tupw1

   end do

   ! ============
   ! Second segment
   ! ============
   if(el(2)>0.0_WP) then
   DO nz=1, nsigma-1


      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*b*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_WP*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*b*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP

   if(V_n_filt(nz,el(2))*deltaX2- U_n_filt(nz,el(2))*deltaY2<0.0_WP) then
      Tmean=(1.0_WP-num_ord)*Tmean2+num_ord*(Tmean1+Tmean2)/2.0_WP
   else
      Tmean=(1.0_WP-num_ord)*Tmean1+num_ord*(Tmean1+Tmean2)/2.0_WP
   end if

   c2=-(Je(nz,el(2))*V_n(nz,el(2))*Tmean-Ty)*deltaX2+ &
          (Je(nz,el(2))*U_n(nz,el(2))*Tmean-Tx)*deltaY2
		  
   un1= Je(nz,el(2))*( V_n_filt(nz,el(2))*deltaX2- U_n_filt(nz,el(2))*deltaY2)   ! normal outer velocity
   Tupw1=0.5_WP*(ttf(nz, enodes(1))*(un1+abs(un1))+ttf(nz, enodes(2))*(un1-abs(un1)))
		  fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))-Tupw1           ! Low-order upwind 
          fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))+Tupw1 
		  
			  
		 fct_aec(nz,edge)=fct_aec(nz,edge)+c2+Tupw1

 END DO
  end if
enddo


   DO edge=1, edge2D
   enodes=edge_nodes(:,edge)
   el=edge_tri(:,edge)
   c1=0.0_WP
   c2=0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)

   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! Average diffusive flux
   ! ============
  c1=Je(nz,el(1))*V_n_filt(nz,el(1))*deltaX1- &
       Je(nz,el(1))*U_n_filt(nz,el(1))*deltaY1
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1)) + c1
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2)) - c1
   END DO
   ! ============
   ! Second segment
   ! ============
   if(el(2)>0.0_WP) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   DO nz=1, nsigma-1
   c2=-Je(nz,el(2))*V_n_filt(nz,el(2))*deltaX2+ &
          Je(nz,el(2))*U_n_filt(nz,el(2))*deltaY2
   ttrhs_c(nz,enodes(1))=ttrhs_c(nz,enodes(1))+c2
   ttrhs_c(nz,enodes(2))=ttrhs_c(nz,enodes(2))-c2
   END DO
   end if

 END DO

! ===================
! Vertical advection
! ===================

DO n=1, nod2D

  ! Surface conditions

	  tvert(1)= -Wvel(1,n)*ttf(1,n)

  ! Bottom conditions
	  tvert(nsigma)=0.0_WP


    nz=2  
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    nz=nsigma-1
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    ! the remaining levels    
    DO nz=3,nsigma-2
    ! ------
    ! First-order upwind estimate
    ! ------
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
              ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    ! -----
    ! centered (4th order)
    ! -----
    qc=(ttfold(nz-1,n)-ttfold(nz,n))/(Z(nz-1,n)-Z(nz,n))
    qu=(ttfold(nz,n)-ttfold(nz+1,n))/(Z(nz,n)-Z(nz+1,n))    
    qd=(ttfold(nz-2,n)-ttfold(nz-1,n))/(Z(nz-2,n)-Z(nz-1,n))
    Tmean1=ttfold(nz,n)+(2*qc+qu)*(zbar(nz,n)-Z(nz,n))/3.0_WP
    Tmean2=ttfold(nz-1,n)+(2*qc+qd)*(zbar(nz,n)-Z(nz-1,n))/3.0_WP
    Tmean=(Wvel(nz,n)+abs(Wvel(nz,n)))*Tmean1+(Wvel(nz,n)-abs(Wvel(nz,n)))*Tmean2
    tvert(nz)=-Tupw1
    fct_aec_ver(nz,n)=(-0.5_WP*(num_ord*(Tmean1+Tmean2)*Wvel(nz,n)+(1.0_WP-num_ord)*Tmean)+Tupw1)*area(n)
    end do

    fct_aec_ver(1,n)=0.0_WP
    fct_aec_ver(nsigma,n)=0.0_WP

   DO nz=1,nsigma-1
       a = ttf(nz,n)
      fct_LO(nz,n)=ttf(nz,n)+dt*((fct_LO(nz,n)/area(n) + &
                      (tvert(nz)-tvert(nz+1)))/Jc(nz,n) - &
                       a*(ttrhs_c(nz,n)/area(n)- &
		       (Wvel(nz,n) - Wvel(nz+1,n)))/Jc(nz,n))

   END DO

 ENDDO

!write(*,*) 'max,min fct_LO', maxval(fct_LO), minval(fct_LO)
!write(*,*) 'max,min fct_aec', maxval(fct_aec), minval(fct_aec)
!write(*,*) 'max,min fct_aec_ver', maxval(fct_aec_ver), minval(fct_aec_ver)

end subroutine fct_muscl_LH2

!===========================================================================
SUBROUTINE fct_miura_LH(ttf, ttfold)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
 integer      :: n, nz, edge, ed, nod(2)
 integer      :: nl1, nl2, el1, el2, nod1, nod2
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, delta_sigma
 real(kind=WP) :: a, b, c, d, da, db, dg, fdepth, un1, Tupw1,qc, qu, qd
 real(kind=WP) :: Tmean1, Tmean2, Tmean, num_ord,tvert(nsigma),tvert2
 real(kind=WP) :: Tx, Ty, crt1
 real(kind=WP) :: ttxnod(nsigma-1,nod2D), ttynod(nsigma-1,nod2D)
 real(kind=WP) :: ttf(nsigma-1, nod2D)
 real(kind=WP) :: ttfold(nsigma-1,nod2D)
 real(kind=WP) :: ttrhs_c(nsigma-1,nod2D)
 real(kind=WP) :: ttx(nsigma-1,elem2D),tty(nsigma-1,elem2D)
 
! =================
! AB interpolation
! =================

 ttfold=-(0.5_WP+epsilon)*ttfold+(1.5_WP+epsilon)*ttf
! =================
! Clean the rhs
! =================
	 
     tvert=0.0_WP
     fct_LO=0.0_WP 
     fct_aec=0.0_WP  
     fct_aec_ver=0.0_WP 
     ttrhs_c=0.0_WP

  call tracer_gradient_elements(ttf, ttx, tty)
  call tracer_gradient_nodes(ttx, tty, ttxnod, ttynod)



! =================
! Horizontal advection and diffusion
! =================
  DO edge=1, edge2D
   nod=edge_nodes(:,edge)

   crt1 = 0.5_WP*(mask_ad(nod(1))+mask_ad(nod(2)))

  ! dmean = sum(depth(nod) + eta_n(nod))/2.0_WP
  ! dmean = max(Dmin,dmean)

  ! crt1 = (dmean/(5.0_WP + Dmin) - 1.0_WP)/20.0_WP
  ! crt1=min(1.0_WP,abs(crt1))

   el1 = edge_tri(1,edge)
   el2 = edge_tri(2,edge)

   deltaX1 = edge_cross_dxdy(1,edge)
   deltaY1 = edge_cross_dxdy(2,edge)


   a=r_earth*elem_cos(el1)
   if(el2>0) then
      deltaX2=edge_cross_dxdy(3,edge)
      deltaY2=edge_cross_dxdy(4,edge)

      b=r_earth*elem_cos(el2)
 
   end if

   Tx=0.0_WP
   Ty=0.0_WP
  
   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1


! ================================================================
! Miura upwind implementation (second order, linear reconstruction)
! ================================================================
      if (V_n_filt(nz,el1)*deltaX1- U_n_filt(nz,el1)*deltaY1 > 0) then

         Tmean=ttf(nz, nod(2))+( 0.5_WP*(-edge_dxdy(1,edge)*a+ deltaX1- &
                               U_n_filt(nz,el1)*dt)*ttxnod(nz,nod(2))+&
                               0.5_WP*(-edge_dxdy(2,edge)*r_earth+ deltaY1- &
			      V_n_filt(nz,el1))*ttynod(nz,nod(2)) )*crt1
      else
         Tmean=ttf(nz, nod(1))+( 0.5_WP*(edge_dxdy(1,edge)*a+ deltaX1 - &
                            U_n_filt(nz,el1)*dt)*ttxnod(nz,nod(1))+&
                               0.5_WP*(edge_dxdy(2,edge)*r_earth+ deltaY1 - &
			      V_n_filt(nz,el1)*dt)*ttynod(nz,nod(1))  )*crt1
      end if
	  
   un1= Je(nz,el1)*(-V_n_filt(nz,el1)*deltaX1+ U_n_filt(nz,el1)*deltaY1)   
   Tupw1=0.5_WP*(ttf(nz, nod(1))*(un1+abs(un1))+ttf(nz, nod(2))*(un1-abs(un1))) 
   
   fct_LO(nz,nod(1))=fct_LO(nz,nod(1))-Tupw1           ! Low-order upwind 
   fct_LO(nz,nod(2))=fct_LO(nz,nod(2))+Tupw1   
    c1 = (Je(nz,el1)*V_n_filt(nz,el1)*Tmean-Ty)*deltaX1- &
           (Je(nz,el1)*U_n_filt(nz,el1)*Tmean-Tx)*deltaY1
   fct_aec(nz,edge)=c1+Tupw1
   END DO
   ! ============
   ! the second column
   ! ============
   if(el2>0)then
      DO nz=1, nsigma-1

   ! ============
   ! Miura upwind
   ! ============
         if(V_n_filt(nz,el2)*deltaX2- U_n_filt(nz,el2)*deltaY2<0) then
            Tmean=ttf(nz, nod(2))+(0.5_WP*(-edge_dxdy(1,edge)*b+ deltaX2- &
                                 U_n_filt(nz,el2)*dt)*ttxnod(nz,nod(2))+&
                               0.5_WP*(-edge_dxdy(2,edge)*r_earth+ deltaY2- &
			         V_n_filt(nz,el2)*dt)*ttynod(nz,nod(2)) )*crt1
         else
            Tmean=ttf(nz, nod(1))+(0.5_WP*(edge_dxdy(1,edge)*b+ deltaX2- &
                                 U_n_filt(nz,el2)*dt)*ttxnod(nz,nod(1))+&
                               0.5_WP*(edge_dxdy(2,edge)*r_earth+ deltaY2- &
			         V_n_filt(nz,el2)*dt)*ttynod(nz,nod(1)) )*crt1
         end if
		    un1= Je(nz,el2)*( V_n_filt(nz,el2)*deltaX2- U_n_filt(nz,el2)*deltaY2)  
            Tupw1=0.5_WP*(ttf(nz, nod(1))*(un1+abs(un1))+ttf(nz, nod(2))*(un1-abs(un1)))
		  fct_LO(nz,nod(1))=fct_LO(nz,nod(1))-Tupw1           ! Low-order upwind 
          fct_LO(nz,nod(2))=fct_LO(nz,nod(2))+Tupw1 
		  
         c2 =-(Je(nz,el2)*V_n_filt(nz,el2)*Tmean-Ty)*deltaX2 + &
              (Je(nz,el2)*U_n_filt(nz,el2)*Tmean-Tx)*deltaY2
			  
		 fct_aec(nz,edge)=fct_aec(nz,edge)+c2+Tupw1
      END DO
   end if
   
   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! Average diffusive flux
   ! ============
      c1 = Je(nz,el1)*V_n_filt(nz,el1)*deltaX1- &
           Je(nz,el1)*U_n_filt(nz,el1)*deltaY1
      ttrhs_c(nz,nod(1)) = ttrhs_c(nz,nod(1)) + c1
      ttrhs_c(nz,nod(2)) = ttrhs_c(nz,nod(2)) - c1
   END DO
   ! ============
   ! Second segment
   ! ============
   if(el2>0) then

      DO nz=1, nsigma-1
         c2 = -Je(nz,el2)*V_n_filt(nz,el2)*deltaX2+ &
               Je(nz,el2)*U_n_filt(nz,el2)*deltaY2
         ttrhs_c(nz,nod(1)) = ttrhs_c(nz,nod(1))+c2
         ttrhs_c(nz,nod(2)) = ttrhs_c(nz,nod(2))-c2
      END DO
   end if

 END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO n=1, nod2D

  ! Surface conditions

	  tvert(1)= -Wvel(1,n)*ttf(1,n)

  ! Bottom conditions
	  tvert(nsigma)=0.0_WP


    nz=2  
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    nz=nsigma-1
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_WP*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(n)
    ! the remaining levels    
    DO nz=3,nsigma-2
    ! ------
    ! First-order upwind estimate
    ! ------
    Tupw1=0.5_WP*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
              ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    ! -----
    ! centered (4th order)
    ! -----
    qc=(ttfold(nz-1,n)-ttfold(nz,n))/(Z(nz-1,n)-Z(nz,n))
    qu=(ttfold(nz,n)-ttfold(nz+1,n))/(Z(nz,n)-Z(nz+1,n))    
    qd=(ttfold(nz-2,n)-ttfold(nz-1,n))/(Z(nz-2,n)-Z(nz-1,n))
    Tmean1=ttfold(nz,n)+(2.0_WP*qc+qu)*(zbar(nz,n)-Z(nz,n))/3.0_WP
    Tmean2=ttfold(nz-1,n)+(2.0_WP*qc+qd)*(zbar(nz,n)-Z(nz-1,n))/3.0_WP
    Tmean=(Wvel(nz,n)+abs(Wvel(nz,n)))*Tmean1+(Wvel(nz,n)-abs(Wvel(nz,n)))*Tmean2
    tvert(nz)=-Tupw1
    fct_aec_ver(nz,n)=(-0.5_WP*(num_ord*(Tmean1+Tmean2)*Wvel(nz,n)+(1.0_WP-num_ord)*Tmean)+Tupw1)*area(n)
    end do

    fct_aec_ver(1,n)=0.0_WP
    fct_aec_ver(nsigma,n)=0.0_WP

   DO nz=1,nsigma-1
       a = ttf(nz,n)
      fct_LO(nz,n)=ttf(nz,n)+dt*((fct_LO(nz,n)/area(n) + &
                      (tvert(nz)-tvert(nz+1)))/Jc(nz,n) - &
                       a*(ttrhs_c(nz,n)/area(n)- &
		       (Wvel(nz,n) - Wvel(nz+1,n)))/Jc(nz,n))

   END DO

 ENDDO

end subroutine fct_miura_LH

!=====================================================================================
! It is assumed that velocity is at n+1/2, hence only tracer field
! is AB2 interpolated to n+1/2.
!
! Alexey Androsov
! 12.01.15
!

SUBROUTINE solve_tracer_upwind(ttf, ttfold, stf, stfold)
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  !a USE g_PARSUP
  IMPLICIT NONE
  real(kind=WP), intent(inout) :: ttf(nsigma-1, nod2D)
  real(kind=WP), intent(inout) :: ttfold(nsigma-1, nod2D)
  real(kind=WP), intent(inout) :: stf(nsigma-1, nod2D)
  real(kind=WP), intent(inout) :: stfold(nsigma-1, nod2D)
 real(kind=WP) :: Vel_nor(nod2D)

  integer      :: nod1, nod2, n, nz, ed, el1, el2, me_nod1, me_nod2
  real(kind=WP) :: c1, c2,  flux=0.0
  real(kind=WP) :: tvert(nsigma), svert(nsigma), db, dg
  real(kind=WP) :: Kh, Kvv,area_inv, dmean, un1
  real(kind=WP) :: t_aux(nsigma-1),s_aux(nsigma-1), aux_c(nsigma-1)
  real(kind=WP) :: ttrhs(nsigma-1, nod2D),strhs(nsigma-1, nod2D), trhs_c(nsigma-1,nod2D)



!NR OpenMP domain decomposition of nodes
me_nod1=1
me_nod2=nod2D

  ! =================
  ! AB interpolation
  ! =================

DO n=me_nod1, me_nod2
   ttfold(:,n) = -(0.5_WP+epsilon)*ttfold(:,n)+(1.5_WP+epsilon)*ttf(:,n)  - T_aver
   ttf(:,n)    = ttf(:,n) - T_aver
   stfold(:,n) = -(0.5_WP+epsilon)*stfold(:,n)+(1.5_WP+epsilon)*stf(:,n)  - S_aver
   stf(:,n)    = stf(:,n) - S_aver

   ! =================
   ! Clean the rhs
   ! =================
   ttrhs(:,n)   = 0.0_WP
   strhs(:,n)   = 0.0_WP
   trhs_c(:,n) = 0.0_WP
  END DO


  ! ttfold now contains AB-interpolated tracer! (n+1/2)

  ! =================
  ! Horizontal advection and diffusion and
  !  Average diffusive flux
  ! =================

  DO ed=1, edge2D_in

     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1  = edge_tri(1,ed)
     el2  = edge_tri(2,ed)

     !a  = r_earth*elem_cos(el1)   !NR not used
     !b  = r_earth*elem_cos(el2)

     ! ============
     ! Average diffusive flux
     ! ============


     DO nz=1, nsigma-1

        c1 =   V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) - U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)
        c2 = - V_n_filt(nz,el2)*edge_cross_dxdy(3,ed) + U_n_filt(nz,el2)*edge_cross_dxdy(4,ed)

        aux_c(nz) = Je(nz,el1)*c1 + Je(nz,el2)*c2

        if (c1 > 0._WP .and. c2 > 0._WP ) then
           t_aux(nz) = Je(nz,el1) *c1 * ttfold(nz,nod2) + Je(nz,el2) *c2 * ttfold(nz,nod2)
           s_aux(nz) = Je(nz,el1) *c1 * stfold(nz,nod2) + Je(nz,el2) *c2 * stfold(nz,nod2)
        elseif (c1 > 0._WP) then
           t_aux(nz) = Je(nz,el1) *c1 * ttfold(nz,nod2) + Je(nz,el2) *c2 * ttfold(nz,nod1)
           s_aux(nz) = Je(nz,el1) *c1 * stfold(nz,nod2) + Je(nz,el2) *c2 * stfold(nz,nod1)
         elseif (c2 > 0._WP) then
           t_aux(nz) = Je(nz,el1) *c1 * ttfold(nz,nod1) + Je(nz,el2) *c2 * ttfold(nz,nod2)
           s_aux(nz) = Je(nz,el1) *c1 * stfold(nz,nod1) + Je(nz,el2) *c2 * stfold(nz,nod2)
        else
           t_aux(nz) = Je(nz,el1) *c1 * ttfold(nz,nod1) + Je(nz,el2) *c2 * ttfold(nz,nod1)
           s_aux(nz) = Je(nz,el1) *c1 * stfold(nz,nod1) + Je(nz,el2) *c2 * stfold(nz,nod1)
        end if

     END DO


     trhs_c(1:nsigma-1,nod1) = trhs_c(1:nsigma-1,nod1) + aux_c(1:nsigma-1)
     ttrhs( 1:nsigma-1,nod1) = ttrhs( 1:nsigma-1,nod1) + t_aux(1:nsigma-1)
     strhs( 1:nsigma-1,nod1) = strhs( 1:nsigma-1,nod1) + s_aux(1:nsigma-1)

     trhs_c(1:nsigma-1,nod2) = trhs_c(1:nsigma-1,nod2) - aux_c(1:nsigma-1)
     ttrhs( 1:nsigma-1,nod2) = ttrhs( 1:nsigma-1,nod2) - t_aux(1:nsigma-1)
     strhs( 1:nsigma-1,nod2) = strhs( 1:nsigma-1,nod2) - s_aux(1:nsigma-1)

  END DO

  DO ed = edge2D_in+1, edge2D

     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1    = edge_tri(1,ed)

     ! a = r_earth*elem_cos(el1)   !NR not used

     ! ============
     ! First segment only
     ! ============
     DO nz=1, nsigma-1
        ! ============
        ! Average diffusive flux
        ! ============
        c1 =   V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) - U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)
        aux_c(nz) = Je(nz,el1)*c1

        ! ============
        ! MUSCL type reconstruction
        ! ============
        if(c1 > 0._WP) then
           t_aux(nz) = aux_c(nz) * ttfold(nz, nod2)    ! Tmean
           s_aux(nz) = aux_c(nz) * stfold(nz, nod2)
        else
           t_aux(nz) = aux_c(nz) * ttfold(nz, nod1)    ! Tmean
           s_aux(nz) = aux_c(nz) * stfold(nz, nod1)
        end if
     END DO


     trhs_c(1:nsigma-1,nod1) = trhs_c(1:nsigma-1,nod1) + aux_c(1:nsigma-1)
     ttrhs( 1:nsigma-1,nod1) = ttrhs( 1:nsigma-1,nod1) + t_aux(1:nsigma-1)
     strhs( 1:nsigma-1,nod1) = strhs( 1:nsigma-1,nod1) + s_aux(1:nsigma-1)

     trhs_c(1:nsigma-1,nod2) = trhs_c(1:nsigma-1,nod2) - aux_c(1:nsigma-1)
     ttrhs( 1:nsigma-1,nod2) = ttrhs( 1:nsigma-1,nod2) - t_aux(1:nsigma-1)
     strhs( 1:nsigma-1,nod2) = strhs( 1:nsigma-1,nod2) - s_aux(1:nsigma-1)

  END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Alternative Option with Tx, Ty
! The case above is the special case with Tx=Ty=0, which allows to save some
! computations.
!!$  DO ed=1, edge2D_in
!!$     nod1 = edge_nodes(1,ed)
!!$     nod2 = edge_nodes(2,ed)
!!$     el1  = edge_tri(1,ed)
!!$     el2  = edge_tri(2,ed)
!!$
!!$     !a  = r_earth*elem_cos(el1)   !NR not used
!!$     !b  = r_earth*elem_cos(el2)
!!$
!!$     Kh = 0.5_WP* K_hor * (elem_area(el1)+elem_area(el2)) * scalearea_inv
!!$
!!$     DO nz=1, nsigma-1
!!$
!!$        ! First segment
!!$
!!$        ! ============
!!$        ! Average diffusive flux
!!$        ! ============
  ! elemental gradients are in ttx, tty
  ! They are based on AB-interpolated values.
  !NR No, ttx, tty are not calculated in this subroutine!

!!$        Tx=0.5_WP*Kh*(ttx(nz,el1)+ttx(nz,el2))
!!$        Ty=0.5_WP*Kh*(tty(nz,el1)+tty(nz,el2))
!!$
!!$        ! ============
!!$        ! MUSCL type reconstruction
!!$        ! ============
!!$        if (V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) > U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)) then
!!$           Tmean = ttfold(nz, nod2)
!!$        else
!!$           Tmean = ttfold(nz, nod1)
!!$        end if
!!$        c1= (Je(nz,el1)*V_n_filt(nz,el1)*Tmean - Ty)*edge_cross_dxdy(1,ed)  &
!!$          - (Je(nz,el1)*U_n_filt(nz,el1)*Tmean - Tx)*edge_cross_dxdy(2,ed)
!!$
!!$     ! Second segment
!!$
!!$           ! ============
!!$           ! Linear upwind reconstruction
!!$           ! ============
!!$        if (V_n_filt(nz,el2)*edge_cross_dxdy(3,ed) < U_n_filt(nz,el2)*edge_cross_dxdy(4,ed) ) then
!!$           Tmean= ttfold(nz, nod2)
!!$        else
!!$           Tmean= ttfold(nz, nod1)
!!$        end if
!!$        c2 = - (Je(nz,el2)*V_n_filt(nz,el2)*Tmean -Ty)*edge_cross_dxdy(3,ed) &
!!$             + (Je(nz,el2)*U_n_filt(nz,el2)*Tmean -Tx)*edge_cross_dxdy(4,ed)
!!$
!!$        aux(nz,ed) = c1 + c2
!!$
!!$        ! ============
!!$        ! Average diffusive flux
!!$        ! ============
!!$        ! First segment
!!$        c1_c = Je(nz,el1)*V_n_filt(nz,el1)*edge_cross_dxdy(1,ed)  &
!!$             - Je(nz,el1)*U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)
!!$
!!$        ! Second segment
!!$        c2_c = - Je(nz,el2)*V_n_filt(nz,el2)*edge_cross_dxdy(3,ed) &
!!$               + Je(nz,el2)*U_n_filt(nz,el2)*edge_cross_dxdy(4,ed)
!!$
!!$        aux_c(nz,ed) = c1_c + c2_c
!!$     END DO
!!$
!!$  END DO
!!$
!!$  DO ed = edge2D_in+1, edge2D
!!$     enodes = edge_nodes(:,ed)
!!$     el1    = edge_tri(1,ed)
!!$     Kh     = elem_area(el1)
!!$
!!$     ! a = r_earth*elem_cos(el1)   !NR not used
!!$
!!$     Kh = K_hor * Kh * scalearea_inv
!!$     ! ============
!!$     ! First segment only
!!$     ! ============
!!$     DO nz=1, nsigma-1
!!$        ! ============
!!$        ! Average diffusive flux
!!$        ! ============
!!$        Tx = Kh*ttx(nz,el1)
!!$        Ty = Kh*tty(nz,el1)
!!$
!!$        ! ============
!!$        ! MUSCL type reconstruction
!!$        ! ============
!!$        if(V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) > U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)) then
!!$           Tmean = ttfold(nz, nod2)
!!$        else
!!$           Tmean = ttfold(nz, nod1)
!!$        end if
!!$        aux(nz,ed) = (Je(nz,el1)*V_n_filt(nz,el1)*Tmean - Ty)*edge_cross_dxdy(1,ed)  &
!!$                         - (Je(nz,el1)*U_n_filt(nz,el1)*Tmean - Tx)*edge_cross_dxdy(2,ed)
!!$
!!$        ! ============
!!$        ! Average diffusive flux
!!$        ! ============
!!$        aux_c(nz,ed) = Je(nz,el1)*(V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) &
!!$                                       - U_n_filt(nz,el1)*edge_cross_dxdy(2,ed))
!!$
!!$     END DO
!!$
!!$  END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ===================
  ! Vertical advection and diffusion
  ! ===================

  tvert(1:nsigma)   = 0.0_WP
  svert(1:nsigma)   = 0.0_WP

  DO n=me_nod1, me_nod2
     ! ===========
     ! Fluxes in the column
     ! ===========
     ! Surface forcing
     if (i_vert_diff) then
        tvert(1)= -Wvel(1,n)*ttfold(1,n)
        svert(1)= -Wvel(1,n)*stfold(1,n)
        Kvv=0.0_WP
     else

        tvert(1)= -Wvel(1,n)*ttfold(1,n)
        svert(1)= -Wvel(1,n)*stfold(1,n)
        Kvv=K_ver
     end if
     ! Bottom conditions
     tvert(nsigma)=0.0_WP

     ! Extract first special case nz = 2
     if(Wvel(2,n) > 0.0_WP) then
        dg =  (zbar(2,n) - Z(1,n))*(Z(2,n) - zbar(2,n)) /((Z(3,n) - Z(1,n))*(Z(2,n) - Z(3,n)))
        db = -(zbar(2,n) - Z(1,n))*(Z(3,n) - zbar(2,n)) /((Z(2,n) - Z(1,n))*(Z(2,n) - Z(3,n)))

        tvert(2) = - (ttfold(1,n)*(1.0_WP-dg-db) +ttfold(2,n)*db +ttfold(3,n)*dg)*Wvel(2,n)
        svert(2) = - (stfold(1,n)*(1.0_WP-dg-db) +stfold(2,n)*db +stfold(3,n)*dg)*Wvel(2,n)
     else
        tvert(2) = - ttfold(1,n) *Wvel(2,n)
        svert(2) = - stfold(1,n) *Wvel(2,n)
     end if
     ! tvert(2)= -Tmean*Wvel(2,n)! + &
     !	      Kvv*(ttf(2-1,n)-ttf(2,n))/(Z(2-1,n)-Z(2,n)))

     DO nz=3, nsigma-2

        ! ============
        ! QUICK upwind (3rd order)
        ! ============
        if(Wvel(nz,n) >= 0.0_WP) then
           !a = zbar(nz,n) - Z(nz-1,n)
           !b = Z(nz,n) - zbar(nz,n)
           !c = Z(nz+1,n) - zbar(nz,n)
           dg =  (zbar(nz,n) - Z(nz-1,n))*(Z(nz,n) - zbar(nz,n)) &
                /((Z(nz+1,n) - Z(nz-1,n))*(Z(nz,n) - Z(nz+1,n)))    !  a*b/((c+a)*(b-c))
           db = -(zbar(nz,n) - Z(nz-1,n))*(Z(nz+1,n) - zbar(nz,n)) &
                /((Z(nz,n)   - Z(nz-1,n))*(Z(nz,n) - Z(nz+1,n)))    ! -a*c/((b+a)*(b-c))

           tvert(nz) = -(ttfold(nz-1,n)*(1.0_WP-dg-db) +ttfold(nz,n)*db +ttfold(nz+1,n)*dg)*Wvel(nz,n)
           svert(nz) = -(stfold(nz-1,n)*(1.0_WP-dg-db) +stfold(nz,n)*db +stfold(nz+1,n)*dg)*Wvel(nz,n)

        else
           !a = Z(nz,n) - zbar(nz,n)
           !b = zbar(nz,n) - Z(nz-1,n)
           !c = zbar(nz,n) - Z(nz-2,n)
           dg =  (Z(nz,n) - zbar(nz,n))*(zbar(nz,n) - Z(nz-1,n)) &
                /((Z(nz,n) - Z(nz-2,n) )*(Z(nz-2,n) - Z(nz-1,n)))    !  a*b/((c+a)*(b-c))
           db = -(Z(nz,n) - zbar(nz,n))*(zbar(nz,n) - Z(nz-2,n)) &
                /((Z(nz,n) - Z(nz-1,n) )*(Z(nz-2,n) - Z(nz-1,n)))    ! -a*c/((b+a)*(b-c))

           tvert(nz) = -(ttfold(nz,n)*(1.0_WP-dg-db) +ttfold(nz-1,n)*db +ttfold(nz-2,n)*dg)*Wvel(nz,n)
           svert(nz) = -(stfold(nz,n)*(1.0_WP-dg-db) +stfold(nz-1,n)*db +stfold(nz-2,n)*dg)*Wvel(nz,n)
        end if

        ! tvert(nz)= -Tmean*Wvel(nz,n)! + &
        !	      Kvv*(ttf(nz-1,n)-ttf(nz,n))/(Z(nz-1,n)-Z(nz,n)))
     END DO

     ! And the other special case: nz = nsigma-1
     if(Wvel(nsigma-1,n) >= 0.0_WP) then
        tvert(nsigma-1) = - ttfold(nsigma-1,n) * Wvel(nsigma-1,n)
        svert(nsigma-1) = - stfold(nsigma-1,n) * Wvel(nsigma-1,n) ! or replace this with the first order upwind
                             ! tttfold(nsigma-1,n)
     else
        dg =  (Z(nsigma-1,n) - zbar(nsigma-1,n))*(zbar(nsigma-1,n) - Z(nsigma-2,n)) &
            /((Z(nsigma-1,n) - Z(nsigma-3,n) )  *(Z(   nsigma-3,n) - Z(nsigma-2,n)))
        db = -(Z(nsigma-1,n) - zbar(nsigma-1,n))*(zbar(nsigma-1,n) - Z(nsigma-3,n)) &
            /((Z(nsigma-1,n) - Z(nsigma-2,n) )  *(Z(   nsigma-3,n) - Z(nsigma-2,n)))

        tvert(nsigma-1) = -(ttfold(nsigma-1,n)*(1.0_WP-dg-db) +ttfold(nsigma-2,n)*db &
                          + ttfold(nsigma-3,n)*dg ) *  Wvel(nsigma-1,n)
        svert(nsigma-1) = -(stfold(nsigma-1,n)*(1.0_WP-dg-db) +stfold(nsigma-2,n)*db &
                          + stfold(nsigma-3,n)*dg ) *  Wvel(nsigma-1,n)
     end if

     ! tvert(nsigma-1)= -Tmean*Wvel(nsigma-1,n)! + &
     !	      Kvv*(ttf(nsigma-2,n)-ttf(nsigma-1,n))/(Z(nsigma-2,n)-Z(nsigma-1,n)))


     area_inv = 1._WP / area(n)
     DO nz=1,nsigma-1
        ttrhs(nz,n) = ( (ttrhs(nz,n) - ttfold(nz,n)*trhs_c(nz,n))*area_inv   &
                     + ttfold(nz,n)*(Wvel(nz,n) - Wvel(nz+1,n)) + (tvert(nz)-tvert(nz+1))) &
                   / Jc(nz,n)
        strhs(nz,n) = ( (strhs(nz,n) - stfold(nz,n)*trhs_c(nz,n))*area_inv   &
                     + stfold(nz,n)*(Wvel(nz,n) - Wvel(nz+1,n)) + (svert(nz)-svert(nz+1))) &
                   / Jc(nz,n)
     END DO
  END DO

if (cl_relax) then

  do n=1,nod2D
     Vel_nor(n) = 0.0_WP
  enddo

  DO ed = 1+edge2D_in, edge2D
     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1  = edge_tri(1,ed)

     dmean = max(Dmin, 0.5_WP*(eta_n(nod1)+depth(nod1) + eta_n(nod2)+depth(nod2)))

     un1 = (V_filt_2D(el1)*edge_cross_dxdy(1,ed)- U_filt_2D(el1)*edge_cross_dxdy(2,ed))*dmean

     Vel_nor(nod1) = Vel_nor(nod1) + un1/area(nod1)
     Vel_nor(nod2) = Vel_nor(nod2) - un1/area(nod2)
  END DO

  DO n=1, nod2D
     if (Vel_nor(n) > 0.0_WP) then
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) then
       ttrhs(nz,n)=ttrhs(nz,n) +relax2clim_ac*(Tclim(nz,n)-(ttfold(nz,n)+T_aver))
       strhs(nz,n)=strhs(nz,n) +relax2clim_ac*(Sclim(nz,n)-(stfold(nz,n)+S_aver))
    endif
     END DO
     else
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) then
       ttrhs(nz,n)=clim_relax*(Tclim(nz,n)-(ttfold(nz,n)+T_aver))
       strhs(nz,n)=clim_relax*(Sclim(nz,n)-(stfold(nz,n)+S_aver))
    endif
    END DO
     endif
  END DO

endif ! climatology

  DO n=1,nod2D
     DO nz=1,nsigma-1
        ttfold(nz,n) = ttf(nz,n) + T_aver
        ttf(nz,n)    = ttf(nz,n) + ttrhs(nz,n)*dt + T_aver
        stfold(nz,n) = stf(nz,n) + S_aver
        stf(nz,n)    = stf(nz,n) + strhs(nz,n)*dt + S_aver
     END DO
  END DO

 !  ===============================================
! relaxation to climat near solid boundary
!===============================================
  if (cl_solid_boder_relax) then

  DO n=1, nod2D
     DO nz=1,nsigma-1
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Tclim(nz,n)-ttf(nz,n))
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Sclim(nz,n)-ttf(nz,n))
     END DO
  END DO

  end if
   !  DO n=1, nod2D
   !     DO nz=1,nsigma-1
   !        if (abs( SF(nz,n))< 1.d-12) SF(nz,n) = abs(SF(nz,n))
   !         if (SF(nz,n)< 0.0_WP) then
   !         write(*,*) 'Negative salinity, time resolution is insufficient, stop', SF(nz,n), nz, n
   !        SF(nz,n)=0.0_WP
   !         endif
   !          if (TF(nz,n)< -1.81_WP) then
   !          TF(nz,n)= -1.81_WP
   !         endif
   !     enddo
   !  END DO

end subroutine solve_tracer_upwind

subroutine tracer_riv(ttf, ttfold, stf, stfold)
!VF, 31.08.16
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  !a USE g_PARSUP
  IMPLICIT NONE
  real(kind=WP), intent(inout) :: ttf(nsigma-1, nod2D)
  real(kind=WP), intent(in) :: ttfold(nsigma-1, nod2D)
  real(kind=WP), intent(inout) :: stf(nsigma-1, nod2D)
  real(kind=WP), intent(in) :: stfold(nsigma-1, nod2D)

  integer      :: nod1, nod2, n, nz, ed, el1, i, nn
  real(kind=WP) :: dmean, un1
  real(kind=WP) :: Vel_nor(nod2D)
   

 if (riv_OB) then

  do n=1,nod2D
     Vel_nor(n) = 0.0_WP
  enddo

  DO ed = 1+edge2D_in, edge2D
     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1  = edge_tri(1,ed)

     dmean = max(Dmin, 0.5_WP*(eta_n(nod1)+depth(nod1) + eta_n(nod2)+depth(nod2)))

     un1 = (V_filt_2D(el1)*edge_cross_dxdy(1,ed)- U_filt_2D(el1)*edge_cross_dxdy(2,ed))*dmean

     Vel_nor(nod1) = Vel_nor(nod1) + un1/area(nod1)
     Vel_nor(nod2) = Vel_nor(nod2) - un1/area(nod2)
  END DO

  !VF, riv_OB

     DO n=1, riv_amount_ob
        nn=riv_node_ob(n)
        if (Vel_nor(nn) > 0.0_WP) then
           DO nz=1,nsigma-1
              ttf(nz,nn) = ttf(nz,nn) + dt*relax2riv*(Tr_distr2(nz,n)-ttfold(nz,nn))
              stf(nz,nn) = stf(nz,nn) + dt*relax2riv*(Sr_distr2(nz,n)-stfold(nz,nn))
           END DO
        else
           DO nz=1,nsigma-1
              ttf(nz,nn) = ttfold(nz,nn) + dt*riv_relax*(Tr_distr2(nz,n)-ttfold(nz,nn))
              stf(nz,nn) = stfold(nz,nn) + dt*riv_relax*(Sr_distr2(nz,n)-stfold(nz,nn))
           END DO
        endif
     END DO

 endif

 !VF, riv
  if (riv) then


     do i=1, riv_num_nodes
        DO nz=1,nsigma-1
              
           ttf(nz,riv_node(i)) = ttf(nz,riv_node(i)) &
                + dt*(Tr_node_sig(nz,i)-ttfold(nz,riv_node(i))*Qr_node_sig(nz,i)) &
                   /    (area(riv_node(i))*Jc(nz,riv_node(i)))
              
           stf(nz,riv_node(i)) = stf(nz,riv_node(i)) &
                + dt*(Sr_node_sig(nz,i) - stfold(nz,riv_node(i))*Qr_node_sig(nz,i)) &
                /    (area(riv_node(i))*Jc(nz,riv_node(i)))
        enddo
     enddo

  endif

end subroutine tracer_riv
!=====================================================================================
!
! Alexey Androsov
! 03.06.16
!
SUBROUTINE solve_tracer_miura(ttf, ttfold, tracer_type)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
 integer      :: n, nz, edge, ed, nod(2)
 integer      :: nl1, nl2, el1, el2, nod1, nod2
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, flux=0.0 , delta_sigma
 real(kind=WP) :: tvert(nsigma), a, b, c, d, da, db, dg, fdepth
 real(kind=WP) :: Z1, Z2, Z3, Z4, Z5, dmean, x1,y1, un1
 real(kind=WP) :: Kh, Kvv, Tx, Ty, Tmean, rdata=0.0, crt1, TS_aver

 real(kind=WP) :: ttf(nsigma-1, nod2D)
 real(kind=WP) :: ttfold(nsigma-1,nod2D)
 real(kind=WP) :: ttrhs( nsigma-1,nod2D), ttrhs_c(nsigma-1,nod2D)
 real(kind=WP) :: ttx(nsigma-1,elem2D),tty(nsigma-1,elem2D)
 real(kind=WP) :: Vel_nor(nod2D)
 real(kind=WP) :: ttxnod(nsigma-1,nod2D), ttynod(nsigma-1,nod2D)
character*1  :: tracer_type

! =================
! AB interpolation
! =================

  if(tracer_type=='t') then
     TS_aver = T_aver
  elseif (tracer_type=='s') then
     TS_aver = S_aver
  else
     TS_aver = 0._WP
  endif

!$OMP PARALLEL private(n,nz)
!$OMP DO
  DO n=1, nod2D
     DO nz=1,nsigma-1
        ttfold(nz,n) = -(0.5_WP+epsilon)*ttfold(nz,n) + (1.5_WP+epsilon)*ttf(nz,n) -TS_aver
        ttf(nz,n)    = ttf(nz,n) - TS_aver
     END DO
! =================
! Clean the rhs
! =================
     ttrhs(:,n)   = 0.0_WP
     ttrhs_c(:,n) = 0.0_WP
  END DO

!   =================
! Compute gradients
! =================
!  call tracer_gradient_elements(ttfold, ttx, tty)
  call tracer_gradient_elements(ttf, ttx, tty)
  call tracer_gradient_nodes(ttx, tty, ttxnod, ttynod)

  ! ttfold now contains AB-interpolated tracer! (n+1/2)

  ! elemental gradients are in ttx, tty
  ! They are based on AB-interpolated values.
!$OMP END PARALLEL

! =================
! Horizontal advection and diffusion
! =================
  DO edge=1, edge2D
   nod=edge_nodes(:,edge)

   crt1 = 0.5_WP*(mask_ad(nod(1))+mask_ad(nod(2)))

  ! dmean = sum(depth(nod) + eta_n(nod))/2.0_WP
  ! dmean = max(Dmin,dmean)

  ! crt1 = (dmean/(5.0_WP + Dmin) - 1.0_WP)/20.0_WP
  ! crt1=min(1.0_WP,abs(crt1))

   el1 = edge_tri(1,edge)
   el2 = edge_tri(2,edge)

   deltaX1 = edge_cross_dxdy(1,edge)
   deltaY1 = edge_cross_dxdy(2,edge)
   Kh = elem_area(el1)

   a=r_earth*elem_cos(el1)
   if(el2>0) then
      deltaX2=edge_cross_dxdy(3,edge)
      deltaY2=edge_cross_dxdy(4,edge)

      b=r_earth*elem_cos(el2)
      Kh=0.5_WP*(Kh+elem_area(el2))
   end if

   Kh = K_hor*Kh/scale_area
   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! Average diffusive flux
   ! ============
      if (el2>0) then
         Tx=0.5_WP*Kh*(ttx(nz,el1)+ttx(nz,el2))
         Ty=0.5_WP*Kh*(tty(nz,el1)+tty(nz,el2))
      else
         Tx=Kh*ttx(nz,el1)
         Ty=Kh*tty(nz,el1)
      end if
   Tx=0.0_WP
   Ty=0.0_WP
! ================================================================
! Miura upwind implementation (second order, linear reconstruction)
! ================================================================
      if (V_n_filt(nz,el1)*deltaX1- U_n_filt(nz,el1)*deltaY1 > 0) then

         Tmean=ttf(nz, nod(2))+( 0.5_WP*(-edge_dxdy(1,edge)*a+ deltaX1- &
                               U_n_filt(nz,el1)*dt)*ttxnod(nz,nod(2))+&
                               0.5_WP*(-edge_dxdy(2,edge)*r_earth+ deltaY1- &
			      V_n_filt(nz,el1))*ttynod(nz,nod(2)) )*crt1
      else
         Tmean=ttf(nz, nod(1))+( 0.5_WP*(edge_dxdy(1,edge)*a+ deltaX1 - &
                            U_n_filt(nz,el1)*dt)*ttxnod(nz,nod(1))+&
                               0.5_WP*(edge_dxdy(2,edge)*r_earth+ deltaY1 - &
			      V_n_filt(nz,el1)*dt)*ttynod(nz,nod(1))  )*crt1
      end if
      c1 = (Je(nz,el1)*V_n_filt(nz,el1)*Tmean-Ty)*deltaX1- &
           (Je(nz,el1)*U_n_filt(nz,el1)*Tmean-Tx)*deltaY1
      ttrhs(nz,nod(1))=ttrhs(nz,nod(1))+c1
      ttrhs(nz,nod(2))=ttrhs(nz,nod(2))-c1
   END DO
   ! ============
   ! the second column
   ! ============
   if(el2>0)then
      DO nz=1, nsigma-1
         Tx=0.5*Kh*(ttx(nz,el1)+ttx(nz,el2))
         Ty=0.5*Kh*(tty(nz,el1)+tty(nz,el2))

   ! ============
   ! Miura upwind
   ! ============
         if(V_n_filt(nz,el2)*deltaX2- U_n_filt(nz,el2)*deltaY2<0) then
            Tmean=ttf(nz, nod(2))+(0.5_WP*(-edge_dxdy(1,edge)*b+ deltaX2- &
                                 U_n_filt(nz,el2)*dt)*ttxnod(nz,nod(2))+&
                               0.5_WP*(-edge_dxdy(2,edge)*r_earth+ deltaY2- &
			         V_n_filt(nz,el2)*dt)*ttynod(nz,nod(2)) )*crt1
         else
            Tmean=ttf(nz, nod(1))+(0.5_WP*(edge_dxdy(1,edge)*b+ deltaX2- &
                                 U_n_filt(nz,el2)*dt)*ttxnod(nz,nod(1))+&
                               0.5_WP*(edge_dxdy(2,edge)*r_earth+ deltaY2- &
			         V_n_filt(nz,el2)*dt)*ttynod(nz,nod(1)) )*crt1
         end if
         c2 =-(Je(nz,el2)*V_n_filt(nz,el2)*Tmean-Ty)*deltaX2 + &
              (Je(nz,el2)*U_n_filt(nz,el2)*Tmean-Tx)*deltaY2
         ttrhs(nz,nod(1))=ttrhs(nz,nod(1))+c2
         ttrhs(nz,nod(2))=ttrhs(nz,nod(2))-c2
      END DO
   end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   ! ============
   ! First segment
   ! ============
   DO nz=1, nsigma-1
   ! ============
   ! Average diffusive flux
   ! ============
      c1 = Je(nz,el1)*V_n_filt(nz,el1)*deltaX1- &
           Je(nz,el1)*U_n_filt(nz,el1)*deltaY1
      ttrhs_c(nz,nod(1)) = ttrhs_c(nz,nod(1)) + c1
      ttrhs_c(nz,nod(2)) = ttrhs_c(nz,nod(2)) - c1
   END DO
   ! ============
   ! Second segment
   ! ============
   if(el2>0) then

      DO nz=1, nsigma-1
         c2 = -Je(nz,el2)*V_n_filt(nz,el2)*deltaX2+ &
               Je(nz,el2)*U_n_filt(nz,el2)*deltaY2
         ttrhs_c(nz,nod(1)) = ttrhs_c(nz,nod(1))+c2
         ttrhs_c(nz,nod(2)) = ttrhs_c(nz,nod(2))-c2
      END DO
   end if

 END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ===================
! Vertical advection and diffusion
! ===================

!aa!$OMP PARALLEL DO private(n,flux,rdata,tvert,Kvv,nz,a,b,c,dg,db,da,Tmean)
 DO n=1, nod2D
!a  if(tracer_type=='t') then
!a  flux=-heat_flux(n)/density_0/4200.0_WP
!a  rdata=surf_relax_T*(Tsurf(n)-ttf(1,n))
!a  else if (tracer_type=='s') then
!a  flux=water_flux(n)*SF(1,n)
!a  rdata=surf_relax_S*(Ssurf(n)-ttf(1,n))
!a  end if
  flux=0.0_WP
  rdata=0.0_WP
  ! ===========
  ! Fluxes in the column
  ! ===========
  ! Surface forcing
  if (i_vert_diff) then
     tvert(1)= -Wvel(1,n)*ttf(1,n)
     Kvv=0.0_WP
  else
     tvert(1)= (flux + rdata-Wvel(1,n)*ttf(1,n))
     Kvv=K_ver
  end if
  ! Bottom conditions
  tvert(nsigma)=0.0_WP

  DO nz=2, nsigma-1

      ! ============
      ! QUICK upwind (3rd order)
      ! ============
   if(Wvel(nz,n)>0.0_WP) then
      if(nz==nsigma-1) then
         Tmean= ttf(nz,n)  ! or replace this with
         ! the first order
						     ! upwind  tttfold(nz,n)
      else
         a=zbar(nz,n) - Z(nz-1,n)
         b=Z(nz,n) - zbar(nz,n)
         c=Z(nz+1,n) - zbar(nz,n)
         dg=a*b/((c+a)*(b-c))
         db=-a*c/((b+a)*(b-c))
         da=1.0_WP-dg-db
         Tmean=ttf(nz-1,n)*da+ttf(nz,n)*db+ttf(nz+1,n)*dg
      end if
   end if

      if(Wvel(nz,n)<0.0_WP) then
         if(nz==2) then
            Tmean=ttf(nz-1,n)        ! or ttfold(nz-1,n)
         else
            a=Z(nz,n) - zbar(nz,n)
            b=zbar(nz,n) - Z(nz-1,n)
            c=zbar(nz,n) - Z(nz-2,n)
            dg=a*b/((c+a)*(b-c))
            db=-a*c/((b+a)*(b-c))
            da=1.0_WP-dg-db
            Tmean=ttf(nz,n)*da+ttf(nz-1,n)*db+ttf(nz-2,n)*dg
         end if
      end if
         tvert(nz)= -Tmean*Wvel(nz,n)! + &
!	      Kvv*(ttf(nz-1,n)-ttf(nz,n))/(Z(nz-1,n)-Z(nz,n)))
   END DO

   DO nz=1,nsigma-1
       ttrhs(nz,n) =  ((ttrhs(nz,n) - ttfold(nz,n)*ttrhs_c(nz,n))/area(n) &
		     + ttfold(nz,n)*(Wvel(nz,n) - Wvel(nz+1,n)) &
                      + tvert(nz)-tvert(nz+1))  / Jc(nz,n)
   END DO

 enddo

  !=================
  ! Relax to climatology if needed
  ! =================

 if (cl_relax) then

  do n=1,nod2D
     Vel_nor(n) = 0.0_WP
  enddo

  DO ed = 1+edge2D_in, edge2D
     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1  = edge_tri(1,ed)

     dmean = max(Dmin, 0.5_WP*(eta_n(nod1)+depth(nod1) + eta_n(nod2)+depth(nod2)))

     un1 = (V_filt_2D(el1)*edge_cross_dxdy(1,ed)- U_filt_2D(el1)*edge_cross_dxdy(2,ed))*dmean

     Vel_nor(nod1) = Vel_nor(nod1) + un1/area(nod1)
     Vel_nor(nod2) = Vel_nor(nod2) - un1/area(nod2)
  END DO

  if (cl_relax) then

 if(tracer_type=='t') then
  DO n=1, nod2D
  !   relax2clim_ac = clim_relax*2.0_WP   !/24.0_WP/10.0_WP
     if (Vel_nor(n) > 0.0_WP) then
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=ttrhs(nz,n) +relax2clim_ac*(Tclim(nz,n)-(ttf(nz,n)+T_aver))
     END DO
     else
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=clim_relax*(Tclim(nz,n)-(ttf(nz,n)+T_aver))
    END DO
     endif
  END DO

 end if ! tracer_type = 't'

  if(tracer_type=='s') then

  DO n=1, nod2D
!     relax2clim_ac = clim_relax*2.0_WP     !/24.0_WP/10._WP
    if (Vel_nor(n) > 0.0_WP) then
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=ttrhs(nz,n) +relax2clim_ac*(Sclim(nz,n)-(ttf(nz,n)+S_aver))
     END DO
     else
     DO nz=1,nsigma-1
    if (index_nod2D(n) == 2) ttrhs(nz,n)=clim_relax*(Sclim(nz,n)-(ttf(nz,n)+S_aver))
     END DO
     endif
  END DO

  end if ! tracer_type = 's'
 endif ! cl_relax

endif ! climatology

! =================
! Update ttfold (to be used on the next time level)
! and compute new ttf
! =================
  DO n=1, nod2D
   DO nz=1,nsigma-1
      ttfold(nz,n) = ttf(nz,n) + TS_aver
      !if (index_nod2D(n) /= 2)
      ttf(nz,n) = ttf(nz,n)+ttrhs(nz,n)*dt + TS_aver
     END DO
  END DO
!aa!$OMP END PARALLEL DO

!  ===============================================
! relaxation to climat near solid boundary
!===============================================
  if (cl_solid_boder_relax) then

   if (tracer_type=='t') then
  DO n=1, nod2D
     DO nz=1,nsigma-1
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Tclim(nz,n)-ttf(nz,n))
     END DO
  END DO
  end if
  if (tracer_type=='s') then
  DO n=1, nod2D
     DO nz=1,nsigma-1
      ttf(nz,n)=ttf(nz,n)+relax2clim(n)*dt*(Sclim(nz,n)-ttf(nz,n))
     END DO
  END DO
  end if

  end if


end subroutine solve_tracer_miura
!=====================================================================================

subroutine tracer_impl_vert_diff
USE o_MESH
USE o_PARAM
USE o_ARRAYS
!USE g_PARSUP
IMPLICIT NONE

real(kind=WP)              ::  a(nsigma), b(nsigma), c(nsigma), tr(nsigma)
real(kind=WP)              ::  sr(nsigma)
real(kind=WP)              ::  cp(nsigma), tp(nsigma), sp(nsigma)
real(kind=WP)              ::  t1, s1, dt_inv, m, swr_source,frw_source(nsigma-1),frw
integer                    ::  nz, n

   a = 0.0_WP
   b = 0.0_WP
   c = 0.0_WP
   tr = 0.0_WP
   sr = 0.0_WP
   cp  = 0.0_WP
   tp = 0.0_WP
   sp = 0.0_WP

corr_TS = TF(1,1)

do n=1,nod2D
 if (mask_wd_node(n) /= 0.0_WP) then
   corr_TS = max(corr_TS,maxval(TF(:,n)))
   TF(:,n) = TF(:,n) - T_aver
   SF(:,n) = SF(:,n) - S_aver
endif
enddo

corr_TS = corr_TS/T_maxini

! if (corr_TS .gt. 1.0_WP) corr_TS=corr_TS*10.0_WP
 !VF, create freshwater source(precipitation) and heat source (short wave radiation)
frw_source=0.0_WP
frw=0.0_WP
swr_source=0.0_WP

dt_inv=1.0_WP/dt


DO n=1,nod2D
 if (mask_wd_node(n) /= 0.0_WP) then
   if ((key_atm).and.(mask_wd_node(n) /= 0.0_WP)) call calc_evap_precip(n,frw)
   frw_source(1)=frw
   ! Regular part of coefficients:
   DO nz=2, nsigma-2
      a(nz) =-dt*Kv(nz,n)/((Z(nz-1,n)-Z(nz,n))*(zbar(nz,n)-zbar(nz+1,n)))
      c(nz) =-dt*Kv(nz,n)/((Z(nz,n)-Z(nz+1,n))*(zbar(nz,n)-zbar(nz+1,n)))
      b(nz) = -a(nz)-c(nz)+1.0_WP   !dt_inv
   END DO
   ! The first row
   c(1) = -dt*Kv(2,n)/((Z(1,n)-Z(2,n))*(zbar(1,n)-zbar(2,n)))
   a(1) = 0.0_WP
   b(1) = -c(1)+1.0_WP !dt_inv
   ! The last row
   a(nsigma-1) = -dt*Kv(nsigma-1,n)/((Z(nsigma-2,n)-Z(nsigma-1,n))*(zbar(nsigma-1,n)-zbar(nsigma,n)))
   b(nsigma-1) = -a(nsigma-1)+1.0_WP !dt_inv
   c(nsigma-1) = 0.0_WP
   ! ===========================================
 !VF, check do we have absorption of swr by sea bed
!   if ((swr_bot_refl_part>0.0_WP).and.(key_atm)) call sw_bot_corr(n)
     if ((key_atm).and. (mask_wd_node(n) /= 0.0_WP)) then

 !VF, add heat and fresh sources
      nz=1
      call calc_heat_flux_sw(nz,n,swr_source)
      tr(nz)=c(nz)*(-TF(nz+1,n)+TF(nz,n))+dt*swr_source
      sr(nz)=c(nz)*(-SF(nz+1,n)+SF(nz,n))+dt*frw_source(nz)

      DO nz=2,nsigma-2
         call calc_heat_flux_sw(nz,n,swr_source)
         tr(nz)=-a(nz)*TF(nz-1,n)-c(nz)*TF(nz+1,n)+(a(nz)+c(nz))*TF(nz,n)+dt*swr_source
         sr(nz)=-a(nz)*SF(nz-1,n)-c(nz)*SF(nz+1,n)+(a(nz)+c(nz))*SF(nz,n)+dt*frw_source(nz)
      END DO
      nz=nsigma-1
      call calc_heat_flux_sw(nz,n,swr_source)
      tr(nz)=a(nz)*(TF(nz,n)-TF(nz-1,n))+dt*swr_source
      sr(nz)=a(nz)*(SF(nz,n)-SF(nz-1,n))+dt*frw_source(nz)
   else
      tr(1) = c(1)*(-TF(2,n)+TF(1,n))
      sr(1) = c(1)*(-SF(2,n)+SF(1,n))
      DO nz=2,nsigma-2
         tr(nz) = -a(nz)*TF(nz-1,n)-c(nz)*TF(nz+1,n)+(a(nz)+c(nz))*TF(nz,n)
         sr(nz) = -a(nz)*SF(nz-1,n)-c(nz)*SF(nz+1,n)+(a(nz)+c(nz))*SF(nz,n)
      END DO
      tr(nsigma-1) = a(nz)*(TF(nsigma-1,n)-TF(nsigma-2,n))
      sr(nsigma-1) = a(nz)*(SF(nsigma-1,n)-SF(nsigma-2,n))
   endif

   !     nz=nsigma-1
   !     tr(nz)=-a(nz)*TF(nz-1,n)+a(nz)*TF(nz,n)
   !     sr(nz)=-a(nz)*SF(nz-1,n)+a(nz)*SF(nz,n)

   !  The first row contains also surface forcing
   !
   !   tr(1)=c(1)*(TF(1,n)- TF(2,n)) + dt*(-1.e-6*heat_flux(n)/(4.2_WP)+surf_relax_T*(Tsurf(n)-TF(1,n)))
   !   sr(1)=c(1)*(SF(1,n)-SF(2,n)) ! + &
   !      dt*(SF(1,n)*water_flux(n) +surf_relax_S*(Ssurf(n)-SF(1,n)))

   ! =============================================
   ! The sweep algorithm
   ! initialize c-prime and s,t-prime
   cp(1) = c(1)/b(1)
   tp(1) = tr(1)/b(1)
   sp(1) = sr(1)/b(1)
   ! solve for vectors c-prime and t, s-prime
   do nz = 2,nsigma-1
      m = 1._WP/(b(nz)-cp(nz-1)*a(nz))
      cp(nz) = c(nz)*m
      tp(nz) = (tr(nz)-tp(nz-1)*a(nz))*m
      sp(nz) = (sr(nz)-sp(nz-1)*a(nz))*m
   enddo
   ! initialize x
   tr(nsigma-1) = tp(nsigma-1)
   sr(nsigma-1) = sp(nsigma-1)
   ! solve for x from the vectors c-prime and d-prime
   do nz = nsigma-2, 1, -1
      tr(nz) = tp(nz)-cp(nz)*tr(nz+1)
      sr(nz) = sp(nz)-cp(nz)*sr(nz+1)
   end do

   DO nz=1,nsigma-1
      TF(nz,n) = TF(nz,n)+tr(nz) + T_aver
      SF(nz,n) = SF(nz,n)+sr(nz) + S_aver
   END DO
endif
  ! end for node mask
END DO   !!! cycle over nodes



     DO n=1, nod2D
        DO nz=1,nsigma-1
           ! if (SF(nz,n)< 0.0_WP) then
           ! write(*,*) 'Negative salinity, time resolution is insufficient, stop', SF(nz,n), nz, n
           !SF(nz,n)=0.0_WP
           ! endif
            ! if (TF(nz,n)< -1.81_WP) then
            ! TF(nz,n)= -1.81_WP
            !endif
call densityJM(TF(nz,n),   SF(nz,n),    -Z(nz,n), rho_c(nz,n))
        enddo
     END DO
density_0= sum(density_0+rho_c)/(nod2D*(nsigma-1))

            !write(*,*) '!!!!!!! ',density_0
 END subroutine tracer_impl_vert_diff
!=======================================================================
subroutine pp_diff(nz, node, diffcoeff)
! Computes the Richardson number dependent diffusivity as in Gent & Cane
! parameterization of mixing (see Large and Gent JPO 1999, 449-464).
! (slightly modified Pacanowski-Philander)
!
USE o_PARAM
USE o_MESH
USE o_ARRAYS
IMPLICIT NONE
integer, intent(in)         :: nz, node
real(kind=WP), intent(out)   :: diffcoeff

real(kind=WP)                :: dz_inv, bv, shear, a, rho_up, rho_dn, t, s
integer                     :: nn, elem, ne

  dz_inv=1./(Z(nz,node)-Z(nz-1,node))
  t=TF(nz-1, node)+T_aver
  s=SF(nz-1, node)+S_aver

  call densityJM(t, s, -zbar(nz,node), rho_up)

  t=TF(nz, node)+T_aver
  s=SF(nz, node)+S_aver

  call densityJM(t, s, -zbar(nz,node), rho_dn)


  ! both densities are referenced to the plane
  ! where flux is computed (zbar(nz)
  bv  = -g*dz_inv*(rho_up-rho_dn)/density_0

    if (bv<=0.0_WP) then
        a=1.0_WP
    else
        shear=0.0_WP
	nn=0
	DO ne=1, nod_in_elem2D_num(node)
	   elem=nod_in_elem2D(ne,node)
	   shear=shear+(U_n(nz-1,elem)-U_n(nz,elem))**2 +&
	               (V_n(nz-1,elem)-V_n(nz,elem))**2
           nn=nn+1
	END DO
        shear=shear*dz_inv*dz_inv/nn
	!a=shear/(shear+10*bv)
	 a=shear/(shear+10.0_WP*bv)
    end if


!    diffcoeff=K_ver+0.01_WP*a*a*a
    diffcoeff=K_ver+corr_TS*0.05_WP*a

end subroutine pp_diff

  !
  !======================
  !
 SUBROUTINE swr_soloviev(z,swr_z)
 !VF
    !
    ! Solar radiation absorbed by the ocean at the depth z
    ! following Soloviev (1982)
    !
    ! INPUT:
    ! z:       depth (m), z<0!
    !
    ! OUTPUT:
    ! swr_z: solar radiation absorbed by the ocean at depth z (W/m**2)
    !
    USE o_PARAM
    IMPLICIT NONE
    REAL(kind=WP),INTENT(in):: z
    REAL(kind=WP),INTENT(out):: swr_z
    REAL(kind=WP),DIMENSION(3), PARAMETER :: a= (/0.45,0.27,0.28/) & 
    , b=(/0.07,2.8,71.5/) 
    !REAL(kind=WP),DIMENSION(3), PARAMETER :: a= (/1.0,0.0,0.0/) &
    !, b=(/10.0,0.0,0.0/)
 
    swr_z=sum(a*exp(b*z))

  END SUBROUTINE swr_soloviev

 SUBROUTINE swr_paul_sim(z,R,a,b, swr_z)
 !VF
    !
    ! Solar radiation absorbed by the ocean at the depth z
    ! following Paulson and Simpson (1977) and Jerlov (1976) classification
    !
    ! INPUT:
    !
    ! z:       depth (m), z<0!
    !
    ! OUTPUT:
    ! swr_z: solar radiation absorbed by the ocean at depth z (W/m**2)
    !
    USE o_PARAM
    IMPLICIT NONE
    REAL(kind=WP),INTENT(in):: z,R,a,b
    REAL(kind=WP),INTENT(out):: swr_z

    swr_z=R*exp(z/a)+(1-R)*exp(z/b)

  END SUBROUTINE swr_paul_sim

!=======================================================================
SUBROUTINE sw_zeng_bel(swr_0,z,swr_z)
!VF
    !
    ! Solar radiation absorbed by the ocean at the depth z
    ! following Zeng and Beljaars (2005).
    !
    ! INPUT:
    ! swr_0: solar radiation at the ocean surface (W/m**2)
    ! z:       depth (m), z<0!
    !
    ! OUTPUT:
    ! swr_z: solar radiation absorbed by the ocean at depth z (W/m**2)
    !
    USE o_PARAM
    IMPLICIT NONE
    REAL(kind=WP),INTENT(in):: swr_0,z
    REAL(kind=WP),INTENT(out):: swr_z

     swr_z=(0.065_WP+11.0_WP*z-6.6e-5/z*(1.0_WP-EXP(-z/8.e-4)))

    !
END SUBROUTINE sw_zeng_bel

SUBROUTINE calc_heat_flux_sw(num_sig,i,swr_source)
!VF
    !
    ! Solar radiation heat source calculation swr_source=(qsr*int(d(qsr_z)/d(z))/(Jc*density_0*Cp))
    ! qsr: solar radiation at the ocean surface (W/m**2) with albedo correction from fv_sbc
    ! qsr_z: Solar radiation absorbed by the ocean at the depth z
    ! c_p: specific heat of seawater
    ! num: number of sigma layer
    ! i: number of the node
    ! OUTPUT:
    ! swr_source: radiation heat source (°C/sec) for the num sigma layer
    !
    USE o_PARAM
    USE o_MESH
    USE o_ARRAYS
    USE fv_sbc

    IMPLICIT NONE
    INTEGER, INTENT(in):: num_sig,i
    REAL(kind=WP),INTENT(out):: swr_source
    REAL(kind=WP):: d1,d2,swr1,swr2,c_p,P,T,S,aux,rho

d1=zbar(num_sig,i)
d2=zbar(num_sig+1,i)
T=T_old(num_sig,i)+T_aver
S=S_old(num_sig,i)+S_aver
aux=(d1+d2)*0.5_WP
rho=rho_c(num_sig,i)+density_0
if (mslp(i)==0.0_WP) then
 P=(rho*g*aux+ P_ref)*10.0d-5
else
 P=(rho*g*aux+mslp(i))*10.0d-5
endif

 call sw_cp(T, S, P,c_p)

if (type_swr_body==1) then
if (jer_const) then
  call swr_paul_sim(-d1,Rjc,ajc,bjc, swr1)
  call swr_paul_sim(-d2,Rjc,ajc,bjc,swr2)
else
  call swr_paul_sim(-d1, Rj(i),aj(i),bj(i), swr1)
  call swr_paul_sim(-d2, Rj(i),aj(i),bj(i), swr2)
 endif
else
if (type_swr_body==2) then
 call swr_soloviev(-d1,swr1)
 call swr_soloviev(-d2,swr2)
endif
endif

if (num_sig==1) then
swr_source=(qsr(i)*(swr1-swr2)+qns(i))/(Jc(num_sig,i)*c_p*rho)
else
swr_source= qsr(i)*(swr1-swr2)/(Jc(num_sig,i)*c_p*rho)
endif
!add rest of qsr to last layer: should be checked!
!if (num_sig==nsigma-1) then
!swr_source= swr_source + qsr(i)*(swr2)/(Jc(num_sig,i)*c_p*rho)
!endif

END SUBROUTINE calc_heat_flux_sw

SUBROUTINE calc_evap_precip(i,frw)
!VF
    ! INPUT:
    ! i: number of the node
    ! OUTPUT:
    ! frw_source: evaporation-precipitation freshwater source
    !
    USE o_PARAM
    USE o_MESH
    USE o_ARRAYS
    USE fv_sbc

    IMPLICIT NONE
    INTEGER, INTENT(in) :: i
    REAL(kind=WP),INTENT(out):: frw
    REAL(kind=WP)            :: T,S,rho

T=TF(1,i)+T_aver
S=SF(1,i)+S_aver


 call  densityJM(T,0.0_WP, -zbar(1,i), rho)


 if ((rho+density_0)*Jc(1,i)/=0) then
frw=S*emp(i)/((rho+density_0)*Jc(1,i))
  else
    frw = 0.0_WP
 endif


END SUBROUTINE calc_evap_precip


SUBROUTINE sw_bot_corr(i)
!VF
    USE o_PARAM
    USE fv_sbc

    IMPLICIT NONE
    INTEGER, INTENT(in):: i
    REAL(kind=WP):: swr_source_aux(2)
    REAL(kind=WP):: swr_bot_refl


           call calc_heat_flux_sw(1,i,swr_source_aux(1));
           call calc_heat_flux_sw(nsigma-1,i,swr_source_aux(2));
if (swr_source_aux(1)>0.0_WP) then
if ((swr_source_aux(2)/swr_source_aux(1))>swr_bot_min) then
               swr_bot_refl=swr_source_aux(2)*swr_bot_refl_part
               qsr(i)=qsr(i)-swr_bot_refl
endif
endif
END SUBROUTINE sw_bot_corr

SUBROUTINE sw_cp(S,Ti,Pr,cp0)
!VF

! SW_CP      Heat Capacity (Cp) of sea water
!=========================================================================
!
! USAGE: cp = sw_cp(S,T,P)
!
! DESCRIPTION:
!    Heat Capacity of Sea Water using UNESCO 1983 polynomial.
!
! INPUT:  (all must have same dimensions)
!   S = salinity    [psu ]
!   T = temperature [C]
!   P = pressure    [Bar]
!       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
!
! OUTPUT:
!   cp0 = Specific Heat Capacity  [J kg**-1 C**-1]
!
! REFERENCES:
!    Fofonff, P. and Millard, R.C. Jr
!    Unesco 1983. Algorithms for computation of fundamental properties of
!    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
! Range of validity - S: 0-40, T: 0-35
!=========================================================================
    USE o_PARAM
    USE o_MESH
    USE o_ARRAYS
    USE fv_sbc

    IMPLICIT NONE
    REAL(kind=WP),INTENT(in) :: S,Ti
!    REAL(kind=WP),INTENT(in), OPTIONAL :: Pr
    REAL(kind=WP),INTENT(in) :: Pr
    REAL(kind=WP),INTENT(out):: cp0
    REAL(kind=WP):: c0,c1,c2,c3,c4,a0,a1,a2,a3,a4,b0,del_Cp0t0,T68,P,b1,b2,b3,b4, S0, sS,T
    REAL(kind=WP):: d0,d1,d2,d3,d4,e0,e1,e2,f0,f1,f2,f3,g0,j1,h0,h1,h2,S3_2,Cpst0,del_Cpstp

    if (Ti<0) then
       T=0.0_WP
    else
       T=Ti
    endif

    if (S<0) then
       S0=0.0_WP
       sS=0.0_WP
    else
       S0=S
       sS=S*sqrt(S)
    endif
!if(.not.(PRESENT(Pr))) then
! P= 1.01325_WP
!else
P=Pr
!endif

T68 = T * 1.00024_WP


c0 = 4217.4_WP
c1 =   -3.720283_WP
c2 =    0.1412855_WP
c3 =   -2.654387e-3
c4 =    2.093236e-5

a0 = -7.64357_WP
a1 =  0.1072763_WP
a2 = -1.38385e-3

b0 =  0.1770383_WP
b1 = -4.07718e-3
b2 =  5.148e-5

Cpst0 = (((c4*T68 + c3)*T68 + c2)*T68 + c1)*T68 + c0 + &
        (a0 + a1*T68 + a2*T68**2)*S0 + &
    (b0 + b1*T68 + b2*T68**2)*sS

!------------
a0 = -4.9592e-1
a1 =  1.45747e-2
a2 = -3.13885e-4
a3 =  2.0357e-6
a4 =  1.7168e-8

b0 =  2.4931e-4
b1 = -1.08645e-5
b2 =  2.87533e-7
b3 = -4.0027e-9
b4 =  2.2956e-11

c0 = -5.422e-8
c1 =  2.6380e-9
c2 = -6.5637e-11
c3 =  6.136e-13

del_Cp0t0 =  (((((c3*T68 + c2)*T68 + c1)*T68 + c0)*P + &
             ((((b4*T68 + b3)*T68 + b2)*T68 + b1)*T68 + b0))*P + &
             ((((a4*T68 + a3)*T68 + a2)*T68 + a1)*T68 + a0))*P


d0 =  4.9247e-3
d1 = -1.28315e-4
d2 =  9.802e-7
d3 =  2.5941e-8
d4 = -2.9179e-10

e0 = -1.2331e-4
e1 = -1.517e-6
e2 =  3.122e-8

f0 = -2.9558e-6
f1 =  1.17054e-7
f2 = -2.3905e-9
f3 =  1.8448e-11

g0 =  9.971e-8

h0 =  5.540e-10
h1 = -1.7682e-11
h2 =  3.513e-13

j1 = -1.4300e-12
S3_2  = sS

del_Cpstp = (((((d4*T68 + d3)*T68 + d2)*T68 + d1)*T68 + d0)*S0 + &
             ((e2*T68 + e1)*T68 + e0)*S3_2)*P                + &
        ((((f3*T68 + f2)*T68 + f1)*T68 + f0)*S0            + &
         g0*S3_2)*P**2                                  + &
         (((h2*T68 + h1)*T68 + h0)*S0                      + &
         j1*T68*S3_2)*P**3


 cp0 = Cpst0 + del_Cp0t0 + del_Cpstp

END SUBROUTINE sw_cp

!=======================================================================

!SUBROUTINE tracer_gradient_nodes (ttx, tty, ttxnodes, ttynodes)

! ttx, tty  elemental gradient of tracer
! ttxnodes, ttynodes nodal gradients obtained by averaging
!USE o_PARAM
!USE o_MESH
!USE o_ARRAYS
!IMPLICIT NONE

!real(kind=WP)      :: ttx(nsigma-1,elem2D)
!real(kind=WP)      :: tty(nsigma-1,elem2D)
!real(kind=WP)      :: ttxnodes(nsigma-1,nod2D)
!real(kind=WP)      :: ttynodes(nsigma-1,nod2D)
!integer         :: n, nz, elem, k
!real(kind=WP)   :: tvol, tx, ty

!  DO n=1, nod2D
!     DO nz=1, nsigma-1
 !       tvol=0.0_WP
!	tx=0.0_WP
!	ty=0.0_WP
!	DO k=1, nod_in_elem2D_num(n)
!           elem=nod_in_elem2D(k,n)
!           tvol=tvol+elem_area(elem)
!          tx=tx+ttx(nz,elem)*elem_area(elem)
!	   ty=ty+tty(nz,elem)*elem_area(elem)
!        END DO
!	ttxnodes(nz,n)=tx/tvol
!	ttynodes(nz,n)=ty/tvol
!     END DO
!   END DO

!END SUBROUTINE tracer_gradient_nodes
!=======================================================================

