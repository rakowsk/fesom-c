MODULE o_PARAM
!VF, update, beginning of 2016, Cd_bz, Bp_bz, mixed_BP_grad_OB
!$ include "omp_lib.h"
!+++++++++++++++++++++++++++++++
!  basic parameters
!+++++++++++++++++++++++++++++++
integer, parameter            :: namlst=11  ! namlst used for init_turbulence
integer, parameter            :: WP=8        ! Working precision
real(kind=WP), parameter      :: pi=3.14159265358979_WP
real(kind=WP), parameter      :: rad=pi/180.0_WP
real(kind=WP), parameter      :: jdh=0.0416666666666666667_WP
real(kind=WP)                 :: time_jd0_orig
real(kind=WP)                 :: density_0 ! will be set in  fv_ini
real(kind=WP)                 :: density_0_inv ! 1./density_0, will be set in  fv_ini
real(kind=WP), parameter      :: g=9.8_WP
real(kind=WP), parameter      :: r_earth=6400000.0_WP
real(kind=WP), parameter      :: omega= 2.0_WP*pi/(3600.0_WP*24.0_WP)  !Coriolis parameter
real(kind=WP)                 :: cka = 0.4_WP !Karman constant
real(kind=WP)                 :: P_ref = 101325.0_WP      ! reference atmospheric pressure
!+++++++++++++++++++++++++++++++
!  roughness parameters
!+++++++++++++++++++++++++++++++
real(kind=WP)                 :: z0s_min      ! surface roughness height, minimum value(use in fv_mixing, gotm subroutine)
real(kind=WP)                 :: z0b_min      ! bottom roughness depth length scale, minimum value (use in fv_mixing, gotm subroutine)
real(kind=WP)                 :: za           ! roughness caused by suspended sediment
!+++++++++++++++++++++++++++++++
!  bottom friction
!+++++++++++++++++++++++++++++++
real(kind=WP)                 :: C_d
integer                       :: fic_point
!+++++++++++++++++++++++++++++++
!  diffusion and viscosity parameters
!+++++++++++++++++++++++++++++++
real(kind=WP)                 :: C_Smag ! Smagorinsky coefficient
real(kind=WP)                 :: snul ! minimum diffusion
real(kind=WP)                 :: A_hor
real(kind=WP)                 :: Abh0
real(kind=WP)                 :: scale_area
real(kind=WP)                 :: K_hor=.010_WP
real(kind=WP)                 :: A_ver=0.001_WP                   ! Vertical harm. visc.
real(kind=WP)                 :: K_ver=0.00001_WP
real(kind=WP)                 :: mix_coeff_PP=0.005_WP
real(kind=WP)                 :: PR_num=0.01_WP          ! Prandtl number ratio of momentum diffusivity
!logical, Parameter            :: laplacian=.true.       ! if .false. -----> biharmonic viscosity for 3D momentum eq.
 logical, parameter            :: i_vert_diff=.true.

 real(kind=WP)                 :: H_advtr_crit=5.0_WP   ! critical depth. For deeper part of the region
                                                                              ! horizontal advection for tracer will be higher order
 real(kind=WP)                 :: H_advtr_min=1.0_WP      ! befor this depth advection for tracer will be first order
!++++++++++++++++++++++++++++++++
! critical depth for reduce baroclinic pressure
!++++++++++++++++++++++++++++++++
 real(kind=WP)                 :: H_bpr_crit=3.0_WP      ! critical depth for compute BAROCLINIC PRESSURE GRADIENT = 1.0
 real(kind=WP)                 :: H_bpr_min=2.0_WP      ! befor this depth BAROCLINIC PRESSURE GRADIENT = 0.0 ! between this two value interpolation from 0 to 1
 logical, parameter          :: Mask_Bar_pr=.true.
!++++++++++++++++++++++++++++++
! diffusion for 2D
!++++++++++++++++++++++++++++++
Logical                      :: filt_2D
logical                      :: filt_bi_2D
Logical                      :: bih_2D
!++++++++++++++++++++++++++++++
! diffusion for3D
!++++++++++++++++++++++++++++++
Logical                      :: filt_3D
logical                      :: filt_bi_3D
Logical                      :: bih_3D
real(kind=WP)                :: tau_c
! (kinematic viscosity) to thermal diffusivity
!******************************
! parameter for buffering zone
! near open boundary
!******************************
Real(kind=WP)                 :: Cd_bz  ! Cd_bf -- increasement of the bottom friction coefficien
                                                                 ! in the buffering zone (ac), range from 0-10...
real(kind=WP)                 :: Bp_bz  ! Bp_bz -- increasement of the influence baroclinic pressure
                                                                 ! in the buffering zone (ac). If Bp_bz = 0.0 normal computation
							         ! of the Baroclinic pressure near OB. If Bp >= 1 reduce
                                                               ! baroclinicity near OB   ac**Bp_bz*Bar_pru(v)_3D
logical          :: mixed_BP_grad_OB  ! if true --> near open boundary mixed
                                                                                ! baroclinic pressure gradient with some
			                                                        ! part of modelling T and S and some part
									        ! with climatology data (Tclim, Sclim)
									        ! NOW linear connection between
										! climatology BP and modelling BP
										! BP ---> baroclinic pressure
										! (1-ac)*Bar_pru(v)_3D_cl + ac*Bar_pru(v)_3D
integer                              :: BFC
!+++++++++++++++++++++++++++++++
! parameters of the critical depth
!+++++++++++++++++++++++++++++++
real(kind=WP)                 :: Dcr, Dmin  ! Critical depth and Minimal depth
!+++++++++++++++++++++++++++++++
! time step parameters
!+++++++++++++++++++++++++++++++
real(kind=WP)                 :: dt, dt_2D   ! time step for baroclinic an barotropic modes
real(kind=WP)                 :: dt_restart, dt_2D_restart   ! time step for baroclinic an barotropic modes from restart file
integer                      :: Tstepping=1  !AB=1, RK=2    ABx=1
integer                      :: nsteps ! total amount of steps
integer                      :: n_dt   ! current step


!real(kind=WP)                 :: w_m2=2.0_WP*pi/(12.42_WP*3600._WP), w_m4=2*pi/(6.21*3600.)
real(kind=WP)                                 :: time, time0, time_2D, time_jd0, time_jd
real(kind=WP)                 :: beta=0.281105_WP !5.0_8/12.0_8 ! Coeff in AB3
real(kind=WP)                 :: alpha(4)=(/0.25_8,1._8/3.0_8, 0.5_8,1.0_8/) ! RK4
real(kind=WP)              :: epsilon=0.01_WP                 ! AB2 offset
real(kind=WP)              :: T_aver, S_aver, c_aver, T_maxini, corr_TS
!real(kind=WP)              :: relax2clim = 5787e-9	! relaxetion to climatology in 1/sec
real(kind=WP)              :: clim_relax, relax2clim_ac
real(kind=WP)              :: relax2riv , riv_relax


		                         ! 0 ---> uniform sigma layers
					 ! 1 --->  non-uniform (near surface)
!+++++++++++++++++++++++++++++++
!  advection scheme
!+++++++++++++++++++++++++++++++
integer                      :: mom_adv_2D
integer                      :: mom_adv_3D
integer                      :: tracer_adv, vert_adv
real(kind=WP)                :: mu_w, up_w
logical                      :: ex_vert_visc, im_vert_visc, im_vert_visc_fic
!++++++++++++++++++++++++++++++
!  vertical mixing scheme
!++++++++++++++++++++++++++++++
integer                     :: ver_mix
Real(kind=WP)               :: beta_scale=0.0_WP      ! 0<=beta_scale<=4  Cut off function for scale
                                                      ! of turbulence profile. IF beta_scale=0 eq. Montgomery
logical, parameter           :: charnock=.true.       ! used in the GOTM subroutine
!++++++++++++++++++++++++++++++
!  short wave radiation penetration
!++++++++++++++++++++++++++++++
real(kind=WP)                ::  swr_bot_refl_part, swr_bot_min !specific heat of seawater; part of swr, which is reflected from bottom ([0,1])
                                                                     !reflection is orginized if the ratio between the surface and

                                                                  !bottom values of swr is greater than swr_bot_min (default 0.01).


real(kind=WP), allocatable   :: aj(:),bj(:),Rj(:) !Jerlov_1976 classification
real(kind=WP)                :: ajc,bjc,Rjc
logical                      :: jer_const      !a,b - attenuation coefficients (blue-green)
logical                      :: key_atm        ! switch ON or OFF atmospheric forcing
logical                      :: key_obc        ! switch ON or OFF open boundary for 3D (TF,SF,..)
logical                      :: key_obc_2d     ! switch ON or OFF open boundary for 2D (eta_n)
logical                      :: key_ic         ! switch ON or OFF module for interpolation of initial conditions
integer                      :: type_swr_body               !R-percent of the flux associated with the longer wavelength irradiance
                                                           ! See also Simpson and Dickey (1981)
integer                      :: Atm_pr

logical                      :: key_nc_output  ! switch ON or OFF netCDF output

!++++++++++++++++++++++++++++++
! Tides
!++++++++++++++++++++++++++++++
logical                      :: TF_presence, T_out, T_out_el_va, T_out_3Dvel, T_potential
Integer                      :: T_aval, T_num, T_counter=0, a_th(15)=0, n_th
Real(kind=WP)                :: T_period, T_step
CHARACTER (LEN=80)           :: time_tf
character (len=3) :: Harmonics_tf_base(15) =(/'M2 ', 'S2 ', 'N2 ', 'K2 ','K1 ', 'O1 ', &
'P1 ', 'Q1 ', 'MM ', 'MF ', 'M4 ', 'MN ', 'MS ', '2N ', 'S1 '/)
real(kind=WP) :: Harmonics_tf_period(15) =(/12.4206_WP, 12.00_WP, 12.6583_WP, 11.9672_WP, 23.9345_WP, &
25.8193_WP, 24.0659_WP, 26.8684_WP, 661.3093_WP, 327.8590_WP, 6.2103_WP,6.2692_WP,6.1033_WP,12.9054_WP,24.0_WP/)
!++++++++++++++++++++++++++++++++++++++++++++++++
! Open boundary
!++++++++++++++++++++++++++++++++++++++++++++++++
real(kind=WP), allocatable    :: eta_n_obc(:) ! array of ssh at the open boundary, constructed in fv_obc_2d
real(kind=WP)                 :: amA   ! amplification of the amplitude on the open boundary (for example: from 0.1 to 1)

!++++++++++++++++++++++++++++++++++++++++++++++++++++
! SEDIMENT
!++++++++++++++++++++++++++++++++++++++++++++++++++++
logical                        :: comp_sediment, sed_boun_flux, cl_relax_sed
real(kind=WP)                  :: d_m            !         d_m  - grain -size [m], w_s - settling velocity [m/s]
real(kind=WP)                  :: plop_s                !         plop_s - density of individual grains [kg/m3]
Real(kind=WP)                  :: e_p                  !         e_p       - porous of the grain
real(kind=WP)                  :: snu_kin=1.0068e-6     !         kinematic viscosity coefficient [m2/s]
Real(kind=WP)                  :: z0r=0.3_WP            !         z0     - raf. parameter
Real(kind=WP)                  :: Er0                    !         empirical floc erosion rate
Real(kind=WP)                  :: sed_relax, relax2sed  !      relaxation parameters

!++++++++++++++++++++++++++++++
! Restart parameters
!++++++++++++++++++++++++++++++
logical                      :: restart
Integer                      :: irestart
!++++++++++++++++++++++++++++++
! Output
!++++++++++++++++++++++++++++++
Integer                      :: irep
Integer                      :: irec
Integer                      :: iverbosity

!++++++++++++++++++++++++++++++
! OBC sponge
!++++++++++++++++++++++++++++++
logical                      :: SL_obn
logical                      :: SL_obn_user_def
real(kind=WP)                :: SL_radius

!++++++++++++++++++++++++++++++
! OBC nudging
!++++++++++++++++++++++++++++++
logical                      :: OB_nudg
real(kind=WP)                :: radius
real(kind=WP)                :: ALPHA_nudg

!++++++++++++++++++++++++++++++
! sigma
!++++++++++++++++++++++++++++++
 integer              :: nsigma ! number of layers
integer              :: ind_sigma
real(kind=WP)        :: pow_sig
real(kind=WP)        :: lev_bot, lev_sur
real(kind=WP)        :: D_ref_max
real(kind=WP), allocatable, dimension(:)  ::  KSw, KBw
!++++++++++++++++++++++++++++++
! General parameters
!++++++++++++++++++++++++++++++
Integer                      :: ini_time = 1, type_task
Integer                      :: Mt  ! number of steps for barotropic mode
integer                      :: nobn, im_index ! amount of open boundary nodes; index of the node with maximum imbalance
real(kind=WP)                :: cyclic_length, ib, ib_mean, T_const, S_const
logical                      :: cartesian
logical                      :: barotropic, barotropic_3D
logical                      :: salinity_on
logical                      :: temp_on
logical                      :: wet_dry_on
logical                      :: riv,riv_control,riv_control_ob, Stratif_sigma, Q_sigma, riv_OB
integer                      :: riv_amount, riv_num_nodes, riv_amount_ob
CHARACTER (LEN=80)           :: title
logical                      :: Coriolis_TF
real(kind=WP)                :: lat_cartesian
logical, parameter           :: free_slip=.true.
logical, save                :: lfirst=.true.
logical                      :: cl_relax
logical                      :: cl_solid_boder_relax


END MODULE o_PARAM


!==========================================================
MODULE o_MESH
USE o_PARAM
integer                ::   nod2D      ! the number of 2D nodes
real(kind=WP), allocatable, dimension(:,:)  ::   coord_nod2D
integer                ::   edge2D     ! the number of 2D edges
integer                ::   edge2D_in  ! the number of internal 2D edges
integer                ::   elem2D     ! the number of 2D elements
integer                ::   elem2D_tri ! the number of triangles among them
integer, allocatable, dimension(:,:)       ::   elem2D_nodes
                                       ! elem2D_nodes(:,n) lists
				       ! 3 nodes of element n
integer, allocatable, dimension(:,:)       ::   edge_nodes
                                       ! edge_nodes(:,n) lists 2 nodes
				       ! edge n
integer, allocatable, dimension(:,:)       ::   edge_tri
                                       ! edge_tri(:,n) lists 2
				       ! triangles containing edge n
				       ! The first one is to left
				       ! of the line directed from the first
				       ! to the second node
integer, allocatable, dimension(:,:)       ::   elem_edges
                                       ! elem_edges(:,n) are edges of
                                       ! element n.
real(kind=WP), allocatable, dimension(:)    ::   elem_area, area
                                       ! areas of triangles and scalar
				       ! control volumes
real(kind=WP), allocatable, dimension(:,:)  ::   edge_dxdy
                                       ! x(2)-x(1) in radians

real(kind=WP), allocatable, dimension(:)  ::   edge_leng
                                       ! x(2)-x(1) in meters


real(kind=WP), allocatable, dimension(:,:)  ::   edge_cross_dxdy
                                       ! physical distances from the edge center to
				       ! left and right triangle
integer,allocatable,dimension(:,:)         ::   elem_neighbors
                                       ! Three neighb. triangles
Integer, allocatable, Dimension(:,:)    :: edge_up_dn_tri
integer,allocatable,dimension(:)         ::   nod_in_elem2D_num, in_obn
integer,allocatable,dimension(:,:)         ::   nod_in_elem2D
real(kind=WP),allocatable,dimension(:)      ::   depth, ac
                                       ! depth(n) is the depths at
				       ! node n
real(kind=WP),allocatable,dimension(:,:)    ::   gradient_sca, gradient_vec, obn_norm
real(kind=WP),allocatable,dimension(:)      ::   cos_elem, cos_edge
real(kind=WP),allocatable,dimension(:)      ::   coriolis
REal(kind=WP),ALlocatable, target, dimension(:,:)  :: ampt, fazt
real(kind=WP), allocatable, dimension(:)   ::   elem_cos, metric_factor
Real(kind=WP), Allocatable, dimension(:,:)   ::   w_cv,a_tp
Real(kind=WP), Allocatable, dimension(:)   ::C_d_el
real(kind=WP), Allocatable, dimension(:) :: sigma
Integer, allocatable, dimension(:)           :: index_nod2D
integer,allocatable,dimension(:,:)         ::   ne_pos, nn_pos
integer,allocatable,dimension(:)          ::   ne_num, nn_num, obn_edge
character(LEN=80)                           ::   indir
character(LEN=80)                           ::   meshpath
character(LEN=80)                           ::   outdir
END MODULE o_MESH

!===========================================================================
MODULE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
! Arrays are described in subroutine array_setup
! 2D velocity
real(kind=WP), allocatable    :: U_n_2D(:,:), U_n_1(:,:), U_n_2(:,:), UAB(:,:), U_n_2Dav(:,:)
real(kind=WP), allocatable    :: U_rhs_2D(:,:), U_n_2D_old(:,:), U_r(:,:)
real(kind=WP), allocatable    :: UV_rhs(:,:,:), UV2_rhs(:,:)
! ssh
real(kind=WP), allocatable    :: eta_n(:), eta_n_1(:), eta_n_2(:), etaAB(:)
real(kind=WP), allocatable    :: ssh_rhs(:), relax_coef(:), ssh_gp(:)
real(KIND=WP), allocatable    :: vel_grad(:,:), taux(:), tauy(:), taux_node(:), tauy_node(:)
real(KIND=WP), allocatable    :: windx(:), windy(:) ! wind at 10 m in X and Y directions

real(kind=WP), allocatable    :: Jc(:,:), Jc_old(:,:), Je(:,:), Jd(:,:)
real(kind=WP), allocatable    :: Z(:,:), zbar(:,:)     ! node depth

real(kind=WP), allocatable    :: dmean_n(:), eta_p(:), z0b_gotm_el(:)

Real(kind=WP), allocatable    :: hpressure(:,:), Bar_pru_3D(:,:), Bar_prv_3D(:,:), hpre_2D(:)
Real(kind=WP), allocatable    :: Bar_pru_3D_clim(:,:), Bar_prv_3D_clim(:,:)
real(kind=WP), allocatable    :: Bar_pr_2D(:,:), rho_c(:,:)
real(kind=WP), allocatable    :: mslp(:), emp(:)    ! mean sea level pressure [Pascals], evaporation minus precip. [kg/m2/s]

real(kind=WP), allocatable    :: Unode(:,:), Vnode(:,:)
real(kind=WP), allocatable    :: Kv(:,:), Kv2(:,:), Av(:,:), Av_node(:,:)
real(kind=WP), allocatable    :: teps (:,:), tepsb(:,:), tke(:,:), L(:,:)
real(kind=WP), allocatable    :: U_rhs(:,:), V_rhs(:,:)
Real(kind=WP), allocatable    :: U_rhsAB(:,:), V_rhsAB(:,:)
real(kind=WP), Allocatable    :: U_rhs_2D_3D(:,:)
real(kind=WP), allocatable    :: vel_grad_ux(:,:), vel_grad_uy(:,:)
real(kind=WP), allocatable    :: vel_grad_vx(:,:), vel_grad_vy(:,:)
real(kind=WP), allocatable    :: Visc(:,:), Visc2D(:)
real(kind=WP), Allocatable    :: U_puls(:,:), V_puls(:,:)
real(kind=WP), allocatable    :: Unode_p(:,:), Vnode_p(:,:)
real(kind=WP), allocatable    :: U_n(:,:), V_n(:,:), W_n(:,:)
real(kind=WP), Allocatable    :: U_n_filt(:,:), V_n_filt(:,:)
!real(kind=WP), allocatable    :: uvert(:,:)
real(kind=WP), allocatable    :: Wvel(:,:)
real(kind=WP), allocatable    :: vorticity_2D(:), vorticity_3D(:,:)
real(kind=WP), allocatable    :: bt(:,:), snu(:,:)

real(kind=WP), allocatable    :: U_filt_2D(:), V_filt_2D(:)

real(kind=WP), allocatable    :: edge_up_dn_grad(:,:,:)
real(kind=WP), allocatable    :: TF(:,:), T_old(:,:), TF2(:,:), T_old2(:,:)
real(kind=WP), allocatable    :: CF(:,:), w_s(:,:), c_old(:,:), Cclim(:,:), cr_distr2(:,:), cr_distr(:,:,:), sed_vt2(:) ! concentration
real(kind=WP), allocatable    :: cr_distr2_2D(:), cr_distr_2D(:,:)! concentration
integer, allocatable          :: sed_ind(:) 
real(kind=WP), Allocatable    :: SF(:,:), S_old(:,:), SF2(:,:), S_old2(:,:) !  temperature
real(kind=WP), allocatable    :: Tclim(:,:), Sclim(:,:) ! salinity
real(kind=WP), allocatable    :: relax2clim(:)
real(kind=WP), allocatable    :: Zl(:,:)
real(kind=WP), allocatable    :: mask_wd(:), mask_wd_node(:), mask_wd_node_old(:),mask_ad(:), mask_bpr(:)
character (len=3), allocatable    ::Harmonics_tf(:)


!++++++++++++++++++++++++++++++
!  FCT
!++++++++++++++++++++++++++++++
real(kind=WP), allocatable    ::  fct_aec_ver(:,:), fct_aec(:,:), fct_LO(:,:)
real(kind=WP), allocatable    ::  fct_ttf_max(:,:), fct_ttf_min(:,:), fct_plus(:,:), fct_minus(:,:)

!++++++++++++++++++++++++++++++
!  Rivers
!++++++++++++++++++++++++++++++
real(kind=WP), allocatable    :: Qr(:), Qr_t(:,:), ssin(:), scos(:), riv_vt(:), riv_vt2(:), riv_vt_jd(:)
real(kind=WP), allocatable    :: riv_vel(:,:),riv_vel_u(:,:), riv_vel_v(:,:),riv_w(:,:), riv_elev(:), riv_elev_t(:,:)
real(kind=WP), allocatable    :: Qr_sig(:,:), Qr_sig_t(:,:,:), Qr_node_sig(:,:), Qr_node(:), Tr_node_sig(:,:),Sr_node_sig(:,:)
real(kind=WP), allocatable    :: Tr_distr(:,:), Sr_distr(:,:), Tr_distr_t(:,:,:), Sr_distr_t(:,:,:)
real(kind=WP), allocatable    :: Tr_distr2(:,:), Sr_distr2(:,:), Tr_distr_t2(:,:,:), Sr_distr_t2(:,:,:)
integer, allocatable :: riv_ind_el(:), riv_ind_eg(:), riv_node(:), riv_node_ob(:)
real(kind=WP)                 :: Time_riv_begin

!++++++++++++++++++++++++++++++++++++++++++++++++++++
! SEDIMENT
!++++++++++++++++++++++++++++++++++++++++++++++++++++
real(kind=WP), allocatable    :: qbu(:), qbv(:), u_st_MOST(:), v_st_MOST(:), u_st_BOOK(:), v_st_BOOK(:)
real(kind=WP), allocatable    :: hama_v(:), E_sed(:), h_var(:), h_var_old(:), z0b_new_yy(:), z0b_new_xx(:)
real(kind=WP), allocatable    :: con(:), con_bc(:), Er_Dep(:), h_var2(:), h_var_old2(:), qb(:,:), Qr_sed(:)
real(kind=WP), allocatable    :: U2D_node(:), V2D_node(:) ,ZZ_1(:), ZZ_2(:), ZZ_3(:), zzx_MOST(:), zzy_MOST(:)
real(kind=WP), allocatable    :: u_st_stress(:), v_st_stress(:), u_st_stress_n(:), v_st_stress_n(:)
real(kind=WP), allocatable    :: zzx_MOST_n(:), zzy_MOST_n(:), z0b_new_xx_n(:), z0b_new_yy_n(:)
real(kind=WP), allocatable    :: ZZ_1_n(:), ZZ_2_n(:), ZZ_3_n(:)
real(kind=WP), allocatable    :: u_star(:,:), v_star(:,:)

END MODULE o_ARRAYS
!

! ==================================================================

module g_PARSUP

! Variables to organize parallel work
implicit none


character(LEN=150)                           ::   DistPath

  type com_struct
     integer    :: rPEnum                    ! the number of PE I receive info from
     integer, dimension(:), allocatable :: rPE   ! their list
     integer, dimension(:), allocatable :: rptr  ! allocatables to the list of nodes
     integer, dimension(:), allocatable :: rlist ! the list of nodes
     integer    :: sPEnum                    ! send part
     integer, dimension(:), allocatable :: sPE
     integer, dimension(:), allocatable :: sptr
     integer, dimension(:), allocatable :: slist
     integer, dimension(:), allocatable :: req  ! request for MPI_Wait
  end type com_struct

  type(com_struct)   :: com_nod2D
  type(com_struct)   :: com_edge2D
  type(com_struct), target :: com_elem2D
  type(com_struct), target :: com_elem2D_full

  ! MPI Datatypes for interface exchange

  ! Edge fields (2D)
  integer, allocatable       :: s_mpitype_edge2D(:),         r_mpitype_edge2D(:)

  ! Element fields (2D; 2D integer; 3D with nl-1 or nl levels, 1 - 4 values)
  !                 small halo and / or full halo
  integer, allocatable, target :: s_mpitype_elem2D(:,:),       r_mpitype_elem2D(:,:)
  integer, allocatable         :: s_mpitype_elem2D_full_i(:),  r_mpitype_elem2D_full_i(:)
  integer, allocatable, target :: s_mpitype_elem2D_full(:,:),  r_mpitype_elem2D_full(:,:)
  integer, allocatable, target :: s_mpitype_elem3D(:,:,:),     r_mpitype_elem3D(:,:,:)
  integer, allocatable, target :: s_mpitype_elem3D_full(:,:,:),r_mpitype_elem3D_full(:,:,:)

  ! Nodal fields (2D; 2D integer; 3D with nl-1 or nl levels, one, two, or three values)
  integer, allocatable       :: s_mpitype_nod2D(:),     r_mpitype_nod2D(:)
  integer, allocatable       :: s_mpitype_nod2D_i(:),   r_mpitype_nod2D_i(:)
  integer, allocatable       :: s_mpitype_nod3D(:,:,:), r_mpitype_nod3D(:,:,:)

  ! general MPI part
  integer            :: MPIERR
  integer            :: npes
  integer            :: mype
  integer            :: maxPEnum=100
  integer, allocatable, dimension(:)  :: part

  ! Mesh partition
  integer                             :: myDim_nod2D, eDim_nod2D
  integer, allocatable, dimension(:)  :: myList_nod2D
  integer                             :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
  integer, allocatable, dimension(:)  :: myList_elem2D
  integer                             :: myDim_edge2D, eDim_edge2D
  integer, allocatable, dimension(:)  :: myList_edge2D

  integer :: pe_status = 0 ! if /=0 then something is wrong

   integer, allocatable ::  remPtr_nod2D(:),  remList_nod2D(:)
   integer, allocatable ::  remPtr_elem2D(:), remList_elem2D(:)

   logical :: elem_full_flag
 end module g_PARSUP
