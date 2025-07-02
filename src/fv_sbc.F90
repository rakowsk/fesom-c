MODULE fv_sbc
   !!===========================================================================
   !! Ocean forcing:
   !!===========================================================================
   !! History: 0.1 ! 07/2015 I. Kuznetsov
   !!
   !! WARNING: in getcoeffld;nc_readTimeGrid
   !!   module will flip data in infiles for lat from -90 to 90 (NCEP-DOE Reanalysis 2 standart)
   !!   if you are going to use other coordinates in input files, please rewrite getcoeffld and nc_readTimeGrid functions.
   !!! WARNING : for now reading only files were time counts in hours from 1800 year; nc_readTimeGrid ; done FOR NCEP data
   !!! WARNING : move time forward for dt/2 , last is a last+dt from last -1 ; nc_readTimeGrid ; done FOR NCEP data

   !! Description:
   !!   read and interpolate atmpospheric forcing on model grid,
   !!     or use constants from namelist each time step
   !!
   !!   first initialization before first time step
   !!     model will read namelist, made some checks and prepare first interpolation coeficients
   !!   during model time steping, each time step sb_do is calling
   !!     during run model check if there is a time to read new data and construct new interpolation coefients
   !!     and interpolate data for new time step
   !!
   !!   taux and tuay defined and allocated outside of this module, but will be changed in this module
   !!   qns - Downward Non Solar heat flux over the ocean defined here and changed here, to use it outside:
   !!   USE fv_sbc
   !!   emp - evaporation minus precipitation defined here and changed here, to use it outside:
   !!   USE fv_sbc
   !!
   !!   Two formats are used or ASCII or NetCDF for input files
   !! NetCDF:
   !!   we assume that all NetCDF files have identical grid and time variable
   !!   nm_sbc=2  nm_sbc_ftype=2 nm_tauwind=2
   !!
   !! ASCII:
   !!   for now, model can read only stress interpolated on model mesh (nodes) from ASCII files
   !!   time interpolation will be done by model
   !!   nm_sbc=2  nm_sbc_ftype=1 nm_tauwind=1
   !!   file format:
   !!10
   !!2000/01/01 00:00:00
   !!0
   !!0
   !! ...
   !!2000/06/15 00:00:00
   !!0
   !!...
   !!were 10 is number of elem2D, 0 here are values of stress in x/y direction, 2000/01/01 00:00:00 - date string
   !!
   !! tips:
   !!       * if you want to have only one wind field and nothing else:
   !!           use ascii format and provide one field (two files for each x/y)
   !!           put one time string (one time step in files
   !!           model will recognize it and use it as a one stress field for all time steps
   !!       * if you are looking for a constant parameters use namelist to setup it and nm_sbc = 1 (check danger part)
   !! danger (or need to be done):
   !!       * if you use daily mean, then solar radiation will be interpolated in time and you will have lisght during nights
   !!       * with nm_sbc = 1 not everthing implemented, for now it is only wind stress will be transferd to ocean variables
   !!       * only SPHERICAL coordinates for now. it means:
   !!            in case of CARTESIAN sbc module will convert c. to SPHERICAL c. and try to interpolate on SPHERICAL coordinates
   !!
   !!
   !! public:
   !!   sbc_ini  -- inizialization atmpospheric forcing
   !!   sbc_do   -- provide a sbc (surface boundary conditions) each time step
   !!
   USE o_ARRAYS
   USE o_MESH
   USE o_PARAM

   USE g_parsup

!   USE ncar_ocean_fluxes_mode

   IMPLICIT NONE

   include 'netcdf.inc'

   public  sbc_ini  ! routine called before 1st time step (open files, read namelist,...)
   public  sbc_do   ! routine called each time step to provide a sbc fileds (wind,...)
   public  sbc_end  ! routine called after last time step
   public  julday   ! get julian day from date

   private

   real(wp), allocatable, save, dimension(:), public     :: qns   ! downward non solar heat over the ocean [W/m2]
   real(wp), allocatable, save, dimension(:), public     :: qsr   ! downward solar heat over the ocean [W/m2]
   real(wp), allocatable, save, dimension(:), public     :: emp   ! evaporation minus precipitation        [kg/m2/s]

!   real(wp), allocatable, save, dimension(:), public     :: qns_2   ! downward non solar heat over the ocean [W/m2]
!   real(wp), allocatable, save, dimension(:), public     :: qsr_2   ! downward solar heat over the ocean [W/m2]
!   real(wp), allocatable, save, dimension(:), public     :: emp_2   ! evaporation minus precipitation        [kg/m2/s]

!   real(wp), allocatable, save, dimension(:), public     :: taux_node_2 ! wind at 10m        [m/s]
!   real(wp), allocatable, save, dimension(:), public     :: tauy_node_2 ! wind at 10m        [m/s]

!   windx and windy are defined in fv_var and allocated in main
!   real(wp), allocatable, save, dimension(:), public     :: windx ! wind at 10m        [m/s]
!   real(wp), allocatable, save, dimension(:), public     :: windy ! wind at 10m        [m/s]
!   mslp now defined in fv_var and allocated in main module
!   real(wp), allocatable, save, dimension(:), public     :: mslp  ! mean sea level pressure [Pascals]
!


   ! namelists
   integer, save  :: nm_sbc_unit     = 101       ! unit to open namelist file
  !============== namelistatmdata variables ================
   integer, save  :: nm_sbc       = 1        ! data  1= constant, 2=from file
   integer, save  :: nm_sbc_ftype = 1        ! input file type if nm_sbc=2 : 1 = ASCII, 2 = netcdf
   integer, save  :: nm_tauwind   = 2        ! 1 = wind stress, 2 = wind 10 m  .
   character(len=256), save   :: nm_xwind_file = 'xwind.dat' ! name of file with winds/stress, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_ywind_file = 'ywind.dat' ! name of file with winds/stress, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_humi_file  = 'humidity.dat' ! name of file with humidity,  if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_qsr_file   = 'qsr.dat'   ! name of file with solar heat,   if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_qlw_file   = 'qlw.dat'   ! name of file with Long wave,    if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_tair_file  = 'tair.dat'  ! name of file with 2m air temperature, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_prec_file  = 'prec.dat'  ! name of file with total precipitation, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_mslp_file  = 'mslp.dat'  ! name of file with mean sea level pressure, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_cloud_file  = 'cloud.dat'  ! name of file with clouds, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model

   character(len=34), save   :: nm_xwind_var = 'uwnd' ! name of variable in file with wind
   character(len=34), save   :: nm_ywind_var = 'vwnd' ! name of variable in file with wind
   character(len=34), save   :: nm_humi_var  = 'shum' ! name of variable in file with humidity
   character(len=34), save   :: nm_qsr_var   = 'dswrf'! name of variable in file with solar heat
   character(len=34), save   :: nm_qlw_var   = 'dlwrf'! name of variable in file with Long wave
   character(len=34), save   :: nm_tair_var  = 'air'  ! name of variable in file with 2m air temperature
   character(len=34), save   :: nm_prec_var  = 'prate'! name of variable in file with total precipitation
   character(len=34), save   :: nm_mslp_var  = 'mslp' ! name of variable in file with mean sea level pressure
   character(len=34), save   :: nm_cloud_var = 'cloud'! name of variable in file with clouds

   real(wp), save :: nm_xwind0 = 0._wp    ! constant 10m. wind value in i-direction/stress, if nm_sbc=1
   real(wp), save :: nm_ywind0 = 0._wp    ! constant 10m. wind value in j-direction, if nm_sbc=1
   real(wp), save :: nm_humi0  = 0._wp    ! constant humidity
   real(wp), save :: nm_qsr0   = 0._wp    ! constant solar heat
   real(wp), save :: nm_qlw0   = 0._wp    ! constant Long wave
   real(wp), save :: nm_tair   = 0._wp    ! constant 2m air temperature
   real(wp), save :: nm_prec   = 0._wp    ! constant total precipitation
   real(wp), save :: nm_mslp   = 0._wp    ! constant mean sea level pressure
   real(wp), save :: nm_cloud  = 0._wp    ! constant clouds,

   real(wp),public, save :: depth_swr = 5._wp    ! depth of swr penetration
   real(wp), save :: nm_prec_coef = 1._wp ! precipitation will be devide by this constant (3600 for CoastDat, 1 for NCEP) (3600 - total precipitation per hour)
   integer , save :: nm_net_flux  = 0     ! constant for downward longwave heat over the ocean: 0 - downward, 1 - Net downward

   integer, save  :: nm_calc_flux = 0 ! calculate atm. flux =1 (based on GOTM subroutins), =0 use precalculated I0,... (based in NEMO subroutins)
   ! ========== netCDF time param
   integer, save :: nm_nc_iyear = 1948    ! initial year of time axis in netCDF (1948 like CoastDat,1800 NCEP)
   integer, save :: nm_nc_imm = 1         ! initial month of time axis in netCDF
   integer, save :: nm_nc_idd = 1         ! initial day of time axis in netCDF
   real, save :: nm_nc_secstep = 86400.0 ! time units coef (86400 CoastDat, 24 NCEP)

   integer,save            :: warn       ! warning switch node/element coordinate out of forcing bounds

   ! ========== interpolation coeficients
   integer,  allocatable, save, dimension(:)     :: bilin_indx_i ! indexs i for interpolation
   integer,  allocatable, save, dimension(:)     :: bilin_indx_j ! indexs j for interpolation

   real(wp), allocatable, save, dimension(:,:)   :: coef_b ! time inerp coef. b (x=a*t+b)
   real(wp), allocatable, save, dimension(:,:)   :: coef_a ! time inerp coef. a (x=a*t+b)

   real(wp), allocatable, save, dimension(:)   :: datawx_f ! wind X data from first time slice
   real(wp), allocatable, save, dimension(:)   :: datawx_s ! wind X data from second time slice
   real(wp), allocatable, save, dimension(:)   :: datawy_f ! wind Y data from first time slice
   real(wp), allocatable, save, dimension(:)   :: datawy_s ! wind Y data from second time slice

   logical, save :: one_field = .false.! only one field used for forcing

   real(wp), allocatable, save, dimension(:,:)   :: atmdata ! atmosperic data for current time step


  !=================ASCII =======================
   integer, save  :: xwind_file_unit = 102       ! unit to open X wind file
   integer, save  :: ywind_file_unit = 103       ! unit to open Y wind file

   integer, parameter        :: END_OF_FILE = -1
   integer, parameter        :: READ_ERROR  = -2

   integer, save :: datewx_f, datewx_s ! date of first and second profiles wind x
   integer, save :: secwx_f,  secwx_s  ! seconds of first and second profiles wind x
   integer, save :: datewy_f, datewy_s ! date of first and second profiles wind y
   integer, save :: secwy_f,  secwy_s  ! seconds of first and second profiles wind x

!   integer,save  :: unit_deb=110
!============== NETCDF ==========================================

!   character(len=256) :: tmp_str

   integer, parameter :: i_totfl = 9 ! total number of fluxes
   integer, parameter :: i_xwind = 1 ! index of 10m wind velocity (x-component) [m/s]
   integer, parameter :: i_ywind = 2 ! index of 10m wind velocity (y-component) [m/s]
   integer, parameter :: i_humi  = 3 ! index of specific humidity               [kg/kg]
   integer, parameter :: i_qsr   = 4 ! index of solar heat                      [W/m2]
   integer, parameter :: i_qlw   = 5 ! index of Long wave                       [W/m2]
   integer, parameter :: i_tair  = 6 ! index of 2m air temperature              [degK]
   integer, parameter :: i_prec  = 7 ! index of total precipitation (rain+snow) [Kg/m^2/s]
   integer, parameter :: i_mslp  = 8 ! index of mean sea level pressure         [Pascals]
   integer, parameter :: i_cloud = 9 ! index of mean sea level pressure         [0-1]


   type, public ::   flfi_type    !flux file informations
      character(len = 256) :: file_name ! file name
      character(len = 34)  :: var_name  ! variable name in the NetCDF file
   end type flfi_type

  type(flfi_type),save, dimension(i_totfl) :: sbc_flfi  !array for information about flux files

  ! arrays of time, lon and lat in INfiles
   real(wp), allocatable, save, dimension(:)  :: nc_lon
   real(wp), allocatable, save, dimension(:)  :: nc_lat
   real(wp), allocatable, save, dimension(:)  :: nc_time
  ! lenght of arrays in INfiles
   integer,save              :: nc_Nlon
   integer,save              :: nc_Nlat
   integer,save              :: nc_Ntime
   ! time index for NC time array
   integer,save              :: t_indx    ! now time index in nc_time array
   integer,save              :: t_indx_p1 ! now time index +1 in nc_time array

  ! flip latitude from infiles (for example  NCEP-DOE Reanalysis 2 standart)
  integer, save              :: flip_lat ! 1 if we need to flip
!============== NETCDF ==========================================


CONTAINS
   SUBROUTINE nc_readTimeGrid(flf)
   ! Read time array and grid from nc file
      IMPLICIT NONE

      type(flfi_type),intent(in) :: flf

      integer              :: iost !I/O status
      integer              :: ncid      ! netcdf file id
      integer              :: i
      ! ID dimensions and variables:
      integer              :: id_lon
      integer              :: id_lat
      integer              :: id_time
      integer              :: id_lond
      integer              :: id_latd
      integer              :: id_timed
!      integer              :: nf_dims(4) ! dimensions (temporal)
      integer              :: nf_start(4)
      integer              :: nf_edges(4)
      integer              :: zero_year,yyyy,mm,dd
      character(len = 256) :: att_string ! attribute
      integer              :: sbc_alloc                   !: allocation status


      !open file


      iost = nf_open(flf%file_name,NF_NOWRITE,ncid)
      call check_nferr(iost,flf%file_name)
      
      ! get dimensions
      iost = nf_inq_dimid(ncid, "latitude", id_latd)
      call check_nferr(iost,flf%file_name)

      iost = nf_inq_dimid(ncid, "longitude", id_lond)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimid(ncid, "time", id_timed)
      call check_nferr(iost,flf%file_name)
      
      ! get variable id
      ! iost = nf_inq_varid(ncid, "air", id_data)
      ! call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "longitude", id_lon)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "latitude", id_lat)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "time", id_time)
      call check_nferr(iost,flf%file_name)
      !  get dimensions size
      iost = nf_inq_dimlen(ncid, id_latd, nc_Nlat)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimlen(ncid, id_lond, nc_Nlon)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimlen(ncid, id_timed, nc_Ntime)
      call check_nferr(iost,flf%file_name)

      ALLOCATE( nc_lon(nc_Nlon), nc_lat(nc_Nlat), nc_time(nc_Ntime),&
           &      STAT=sbc_alloc )
      if( sbc_alloc /= 0 )   STOP 'read_sbc: failed to allocate arrays'

      !read variables from file
      ! coordinates
      nf_start(1)=1
      nf_edges(1)=nc_Nlat
      iost = nf_get_vara_double(ncid, id_lat, nf_start, nf_edges, nc_lat)
      call check_nferr(iost,flf%file_name)
      nf_start(1)=1
      nf_edges(1)=nc_Nlon
      iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, nc_lon)
      call check_nferr(iost,flf%file_name)
      nf_start(1)=1
      nf_edges(1)=nc_Ntime
      iost = nf_get_vara_double(ncid, id_time, nf_start, nf_edges, nc_time)
      call check_nferr(iost,flf%file_name)
      
      ! convert time to days
      !!NCEP     nc_time = nc_time / 24.0 + julday(1800,1,1)
      nc_time = nc_time / nm_nc_secstep + julday(nm_nc_iyear,nm_nc_imm,nm_nc_idd)
!!! WARNING : move time forward for dt/2 , last is a last+dt from last -1
      if (nc_Ntime > 1) then
         do i = 1, nc_Ntime-1
            nc_time(i) = (nc_time(i+1) + nc_time(i))/2.0
         end do
            
         nc_time(nc_Ntime) = nc_time(nc_Ntime) + (nc_time(nc_Ntime) - nc_time(nc_Ntime-1))/2.0
      end if
      !flip lat and data in case of lat from -90 to 90
!!!! WARNING this is temporal solution, needs some more checks
      flip_lat = 0
      if ( nc_Nlat > 1 ) then
         if ( nc_lat(1) > nc_lat(nc_Nlat) ) then
            flip_lat = 1
            nc_lat=nc_lat(nc_Nlat:1:-1)
            write(*,*) "fv_sbc: nc_readTimeGrid: FLIP lat and data while lat from -90 to 90"
         endif
      endif

      !      iost = nf_get_att(ncid, id_time, "units",att_string)
      !      call check_nferr(iost,flf%file_name)
      !      write(*,*) iost
      !      write(*,*) att_string
      iost = nf_close(ncid)
      call check_nferr(iost,flf%file_name)

   END SUBROUTINE nc_readTimeGrid

   SUBROUTINE nc_sbc_ini_fillnames(yyear)
      character(len=4),intent(in)   :: yyear

      !! ** Purpose : Fill names of sbc_flfi array (file names and variable names)

      !prepare proper nc file (add year and .nc to the end of the file name from namelist
      write(sbc_flfi(i_xwind)%file_name,*) trim(nm_xwind_file)  !,yyear,'.nc'
      write(sbc_flfi(i_ywind)%file_name,*) trim(nm_ywind_file)  !,yyear,'.nc'
      write(sbc_flfi(i_humi)%file_name, *) trim(nm_humi_file)  !,yyear,'.nc'
      write(sbc_flfi(i_qsr)%file_name, *) trim(nm_qsr_file)  !,yyear,'.nc'
      write(sbc_flfi(i_qlw)%file_name, *) trim(nm_qlw_file)  !,yyear,'.nc'
      write(sbc_flfi(i_tair)%file_name, *) trim(nm_tair_file)  !,yyear,'.nc'
      write(sbc_flfi(i_prec)%file_name, *) trim(nm_prec_file)  !,yyear,'.nc'
      write(sbc_flfi(i_mslp)%file_name, *) trim(nm_mslp_file)  !,yyear,'.nc'
      if (nm_calc_flux==1) then
         write(sbc_flfi(i_cloud)%file_name, *) trim(nm_cloud_file),yyear,'.nc'
      end if

      sbc_flfi(i_xwind)%file_name=ADJUSTL(trim(sbc_flfi(i_xwind)%file_name))
      sbc_flfi(i_ywind)%file_name=ADJUSTL(trim(sbc_flfi(i_ywind)%file_name))
      sbc_flfi(i_humi)%file_name=ADJUSTL(trim(sbc_flfi(i_humi)%file_name))
      sbc_flfi(i_qsr)%file_name=ADJUSTL(trim(sbc_flfi(i_qsr)%file_name))
      sbc_flfi(i_qlw)%file_name=ADJUSTL(trim(sbc_flfi(i_qlw)%file_name))
      sbc_flfi(i_tair)%file_name=ADJUSTL(trim(sbc_flfi(i_tair)%file_name))
      sbc_flfi(i_prec)%file_name=ADJUSTL(trim(sbc_flfi(i_prec)%file_name))
      sbc_flfi(i_mslp)%file_name=ADJUSTL(trim(sbc_flfi(i_mslp)%file_name))
      if (nm_calc_flux==1) then
         sbc_flfi(i_cloud)%file_name=ADJUSTL(trim(sbc_flfi(i_cloud)%file_name))
      end if

      sbc_flfi(i_xwind)%var_name=ADJUSTL(trim(nm_xwind_var))
      sbc_flfi(i_ywind)%var_name=ADJUSTL(trim(nm_ywind_var))
      sbc_flfi(i_humi)%var_name=ADJUSTL(trim(nm_humi_var))
      sbc_flfi(i_qsr)%var_name=ADJUSTL(trim(nm_qsr_var))
      sbc_flfi(i_qlw)%var_name=ADJUSTL(trim(nm_qlw_var))
      sbc_flfi(i_tair)%var_name=ADJUSTL(trim(nm_tair_var))
      sbc_flfi(i_prec)%var_name=ADJUSTL(trim(nm_prec_var))
      sbc_flfi(i_mslp)%var_name=ADJUSTL(trim(nm_mslp_var))
      if (nm_calc_flux==1) then
         sbc_flfi(i_cloud)%var_name=ADJUSTL(trim(nm_cloud_var))
      end if

   END SUBROUTINE nc_sbc_ini_fillnames

   SUBROUTINE nc_sbc_ini(idate,isec)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ascii_ini ***
      !!
      !! ** Purpose : initialization of ocean forcing from NETCDF file
      !!----------------------------------------------------------------------

      IMPLICIT NONE
      integer,intent(in) :: idate ! initialization date
      integer,intent(in) :: isec ! initialization seconds

      character(len=4)   :: yyear
      integer            :: yyyy,mm,dd

      integer            :: i
      integer            :: sbc_alloc

      integer            :: elnodes(4) !4 nodes from one element
      integer            :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )
      real(wp)           :: x, y       ! coordinates of elements


! used for interpolate on elements
!      ALLOCATE( bilin_indx_i(elem2D),bilin_indx_j(elem2D), &
!              & qns(elem2D), emp(elem2D), qsr(elem2D),     &
!                   &      STAT=sbc_alloc )
! used to inerpolate on nodes
      warn = 0

      ! get ini year; Fill names of sbc_flfi
      call calendar_date(idate,yyyy,dd,mm)
      write(yyear,"(I4)") yyyy
      call nc_sbc_ini_fillnames(yyear)

      ! we assume that all NetCDF files have identical grid and time variable
      call nc_readTimeGrid(sbc_flfi(i_xwind))

      ! prepare nearest coordinates in INfile , save to bilin_indx_i/j
      do i = 1, myDim_nod2D+eDim_nod2D
!         ! get coordinates of elements
!         elnodes = elem2D_nodes(:,i) !! 4 nodes in element
!         numnodes = 4
!         if(elnodes(1)==elnodes(4)) numnodes = 3  !! set to 3 if we have triangle
!         ! coord_nod2d is array(:,:) for nodes coordinates
!         ! rad is pi/180
!         y = sum(coord_nod2d(2,elnodes(1:numnodes)))/dble(numnodes) /rad
!         x = sum(coord_nod2d(1,elnodes(1:numnodes)))/dble(numnodes) /rad
!  use these lines in case if we use sbc on nodes
         x  = coord_nod2d(1,i)/rad
         y  = coord_nod2d(2,i)/rad

         ! find nearest
         if ( x < nc_lon(nc_Nlon) .and. x >= nc_lon(1) ) then
            call binarysearch(nc_Nlon, nc_lon, x, bilin_indx_i(i))
         else ! NO extrapolation in space
            if ( x < nc_lon(1) ) then
               bilin_indx_i(i)=-1
            else
               bilin_indx_i(i)=0
            end if
         end if
         if ( y < nc_lat(nc_Nlat) .and. y >= nc_lat(1) ) then
            call binarysearch(nc_Nlat, nc_lat, y, bilin_indx_j(i))
         else ! NO extrapolation in space
            if ( y < nc_lat(1) ) then
               bilin_indx_j(i)=-1
            else
               bilin_indx_j(i)=0
            end if
         end if
         if (warn == 0) then
            if (bilin_indx_i(i) < 1 .or. bilin_indx_j(i) < 1) then
               WRITE(*,*) '     WARNING:  node/element coordinate out of forcing bounds,'
               WRITE(*,*) '        nearest value will be used as a constant field'
               warn = 1
            end if
         end if




      end do
      ! get first coefficients for time interpolation on model grid for all data
      call getcoeffld(time_jd)

      ! interpolate in time

      call data_timeinterp(time_jd)

      if (nm_calc_flux==1) then
         call apply_atm_fluxes_gotm(time_jd)
      else
         call apply_atm_fluxes
      end if


   END SUBROUTINE nc_sbc_ini

   SUBROUTINE apply_atm_fluxes
      !!----------------------------------------------------------------------
      !! ** Purpose : Change model variables according to atm fluxes
      !! source of original code: NEMO 3.1.1 + NCAR
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer             :: i
      integer             :: num, n
      real(wp)            :: rtmp    ! temporal real
      real(wp)            :: wndm    ! delta of wind module and ocean curent module
      real(wp)            :: wdx,wdy ! delta of wind x/y and ocean curent x/y
      real(wp)            :: q_sat   ! sea surface specific humidity         [kg/kg]
      real(wp), parameter :: rhoa = 1.22 ! air density
      real(wp), parameter :: cpa  = 1000.5         ! specific heat of air
      real(wp), parameter :: Lv   =    2.5e6       ! latent heat of vaporization
      real(wp), parameter :: Stef =    5.67e-8     ! Stefan Boltzmann constant
      real(wp), parameter :: albo =    0.066       ! ocean albedo assumed to be contant
      real(wp)            :: zst     ! surface temperature in Kelvin


      real(wp)           ::  &
         Cd,       &     ! transfer coefficient for momentum         (tau)
         Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
         Ce,       &     ! transfert coefficient for evaporation   (Q_lat)
         t_zu,     &     ! air temp. shifted at zu                     [K]
         q_zu            ! spec. hum.  shifted at zu               [kg/kg]

      real(wp)           :: zevap, zqsb, zqla, zqlw
!!$OMP PARALLEL
!!$OMP DO
      do i = 1,nod2D
         windx(i) = atmdata(i_xwind,i)
         windy(i) = atmdata(i_ywind,i)
         ! IF you going to use winds on nodes (nod2d) U_n and V_n should be changed to Unode and Vnode
         if (type_task>1) then
            wdx = atmdata(i_xwind,i) - Unode(1,i) ! wind from data - ocean current ( x direction)
            wdy = atmdata(i_ywind,i) - Vnode(1,i) ! wind from data - ocean current ( y direction)
         else
! WARNING:: surface curent U and V should be used here, while it should be interpolate by compute_el2nodes_2D_2fields from fv_utilit
!           so if u use type_task = 1 and winds you are welcome to add some code
            ! convert 2D velocity on element to nodes (wdx), than modify wdx by
            ! this value
            n = i

            num = nod_in_elem2D_num(n)
            wdx = sum(U_n_2D(1,nod_in_elem2D(1:num,n))*elem_area(nod_in_elem2D(1:num,n))) &
                                                  /sum(elem_area(nod_in_elem2D(1:num,n)))
            wdy = sum(U_n_2D(2,nod_in_elem2D(1:num,n))*elem_area(nod_in_elem2D(1:num,n))) &
                                                  /sum(elem_area(nod_in_elem2D(1:num,n)))
           ! if (i == 10000) then
           !  write(*,*) "wdx = ", wdx, " U = ",maxval(U_n_2D(1,nod_in_elem2D(1:num,n)))
           ! end if
            wdx = atmdata(i_xwind,i) - wdx !- U_n_2D(1,i) ! wind from data - ocean current ( x direction)
           ! if (i == 10000) then
           !   write(*,*) "wdx_cor = ", wdx
           ! end if
            wdy = atmdata(i_ywind,i) - wdy !- U_n_2D(2,i) ! wind from data - ocean current ( y direction)     \

            !SH Original part remove if it works
            !SH wdx = atmdata(i_xwind,i)! - U_n_2D(1,i) ! wind from data - ocean current ( x direction)
            !SH wdy = atmdata(i_ywind,i)! - U_n_2D(2,i) ! wind from data - ocean current ( y direction)
         endif
         wndm = SQRT( wdx * wdx + wdy * wdy )
!      call ncar_ocean_fluxes (wndm, atmdata(i_tair,i), ts, q, qs, z, avail, &
!                              cd, ch, ce, ustar, bstar       )
         ! TF (temperature) should be interpolated on elem + Kelvin  (NCEP/NCAR is in K.)
         ! If no temperature calculation, then constant T
         if (type_task>2) then
!        if u use elem2D
!            zst = sum(w_cv(1:4,i)*TF(1,elem2D_nodes(:,i)))+273.15
            zst = TF(1,i)+273.15
         else
            zst = T_const+273.15
         endif

         q_sat = 0.98 * 640380. / rhoa * EXP( -5107.4 / zst )

         call core_coeff_2z(2.0_wp, 10.0_wp, zst, atmdata(i_tair,i), &
                           q_sat, atmdata(i_humi,i), wndm, Cd, Ch, Ce, t_zu, q_zu)

!         call core_coeff_2z(2.0_wp, 10.0_wp, 237.15_wp+5.0_wp, 237.15_wp+5.0_wp, &
!                           q_sat, 0.008_wp, wndm, Cd, Ch, Ce, t_zu, q_zu)

!         Cd = (0.6 + 0.07 * wdx ) /1000

         rtmp = rhoa * wndm * Cd
!         taum(i)    = rtmp * wndm  ! stress module (if we need to save it)

!!!===========================================  changing of taux and tauy here ====================
         taux_node(i) = rtmp * wdx  ! new taux (stress along x)
         tauy_node(i) = rtmp * wdy  ! new tauy (stress along y)
!         taux(i) = rhoa * wdx * wdx * (0.6 + 0.07 * wdx ) /1000
!         tauy(i) = 0.0
!!!===========================================  changing of taux and tauy here ====================

         zevap = MAX( 0.e0, rhoa    *Ce*( q_sat - q_zu ) * wndm )   ! Evaporation
         zqsb  = rhoa*cpa*Ch*( zst - t_zu ) * wndm     ! Sensible Heat

         if ( nm_net_flux == 0 ) then
            ! NCEP case
            zqlw  = (  atmdata(i_qlw,i) - Stef * zst*zst*zst*zst  )    ! Long  Wave
         else
            ! CoastDat case with Net flux
            zqlw  = atmdata(i_qlw,i)
         endif

         zqla  = Lv * zevap   ! Latent Heat
!!!================================ data for use in ocean part ================================
         ! Downward Non Solar heat flux over the ocean
         qns(i) = zqlw - zqsb - zqla
         ! zqlw   = downward longwave heat over the ocean
         ! - zqsb = downward sensible heat over the ocean
         ! - zqla = downward latent   heat over the ocean
         ! qns    = downward non solar heat over the ocean
         emp(i) = zevap - atmdata(i_prec,i)/nm_prec_coef
         qsr(i) = (1.0_wp - albo ) * atmdata(i_qsr,i)
         mslp(i)= atmdata(i_mslp,i)
!!!================================ data for use in ocean part ================================
      end do
!!$OMP END DO
      ! convert tau to elements
!!$OMP DO
      do i = 1, elem2D
         taux(i) = sum(w_cv(1:4,i)*taux_node(elem2D_nodes(:,i)))
         tauy(i) = sum(w_cv(1:4,i)*tauy_node(elem2D_nodes(:,i)))
      end do
!!$OMP END DO
!!$OMP END PARALLEL

   END SUBROUTINE apply_atm_fluxes


   SUBROUTINE core_coeff_2z(zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd, Ch, Ce, T_zu, q_zu)
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  core_coeff_2z  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004).
      !!
      !! ** Method  :   I N E R T I A L   D I S S I P A T I O N   M E T H O D
      !!      Momentum, Latent and sensible heat exchange coefficients
      !!      Caution: this procedure should only be used in cases when air
      !!      temperature (T_air) and air specific humidity (q_air) are at 2m
      !!      whereas wind (dU) is at 10m.
      !!
      !! References :   Large & Yeager, 2004 : ???
      !! code was adopted from NEMO 3.3.1
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      real(wp)            :: dU10        ! dU                             [m/s]
      real(wp)            :: dT          ! air/sea temperature difference   [K]
      real(wp)            :: dq          ! air/sea humidity difference      [K]
      real(wp)            :: Cd_n10      ! 10m neutral drag coefficient
      real(wp)            :: Ce_n10      ! 10m neutral latent coefficient
      real(wp)            :: Ch_n10      ! 10m neutral sensible coefficient
      real(wp)            :: sqrt_Cd_n10 ! root square of Cd_n10
      real(wp)            :: sqrt_Cd     ! root square of Cd
      real(wp)            :: T_vpot      ! virtual potential temperature    [K]
      real(wp)            :: T_star      ! turbulent scale of tem. fluct.
      real(wp)            :: q_star      ! turbulent humidity of temp. fluct.
      real(wp)            :: U_star      ! turb. scale of velocity fluct.
      real(wp)            :: L           ! Monin-Obukov length              [m]
      real(wp)            :: zeta_u      ! stability parameter at height zu
      real(wp)            :: zeta_t      ! stability parameter at height zt
      real(wp)            :: U_n10       ! neutral wind velocity at 10m     [m]
      real(wp)            :: xlogt , xct , zpsi_hu , zpsi_ht , zpsi_m
      real(wp)            :: stab        ! 1st guess stability test integer
      !!
      real(wp), intent(in)   :: &
         zt,      &     ! height for T_zt and q_zt                   [m]
         zu             ! height for dU                              [m]
      real(wp), intent(in)   ::  &
         sst,      &     ! sea surface temperature              [Kelvin]
         T_zt,     &     ! potential air temperature            [Kelvin]
         q_sat,    &     ! sea surface specific humidity         [kg/kg]
         q_zt,     &     ! specific air humidity                 [kg/kg]
         dU              ! relative wind module |U(zu)-U(0)|       [m/s]
      real(wp), intent(out)  ::  &
         Cd,       &     ! transfer coefficient for momentum         (tau)
         Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
         Ce,       &     ! transfert coefficient for evaporation   (Q_lat)
         T_zu,     &     ! air temp. shifted at zu                     [K]
         q_zu            ! spec. hum.  shifted at zu               [kg/kg]

      integer :: j_itt
      integer,  parameter :: nb_itt = 5   ! number of itterations
      real(wp), parameter ::                        &
         grav   = 9.8,      &  ! gravity
         kappa  = 0.4          ! von Karman's constant
      !!----------------------------------------------------------------------
      !!  * Start

      !! Initial air/sea differences
      dU10 = max(0.5_wp, dU)      !  we don't want to fall under 0.5 m/s
      dT = T_zt - sst
      dq = q_zt - q_sat

      !! Neutral Drag Coefficient :
      stab = 0.5 + sign(0.5_wp,dT)                 ! stab = 1  if dT > 0  -> STABLE
      Cd_n10  = 1E-3*( 2.7/dU10 + 0.142 + dU10/13.09 )
      sqrt_Cd_n10 = sqrt(Cd_n10)
      Ce_n10  = 1E-3*( 34.6 * sqrt_Cd_n10 )
      Ch_n10  = 1E-3*sqrt_Cd_n10*(18*stab + 32.7*(1 - stab))

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = Cd_n10 ;  Ce = Ce_n10 ;  Ch = Ch_n10 ;  sqrt_Cd = sqrt(Cd)

      !! Initializing z_u values with z_t values :
      T_zu = T_zt ;  q_zu = q_zt

      !!  * Now starting iteration loop
      do j_itt=1, nb_itt
         dT = T_zu - sst ;  dq = q_zu - q_sat ! Updating air/sea differences
         T_vpot = T_zu*(1. + 0.608*q_zu)    ! Updating virtual potential temperature at zu
         U_star = sqrt_Cd*dU10                ! Updating turbulent scales :   (L & Y eq. (7))
         T_star  = Ch/sqrt_Cd*dT              !
         q_star  = Ce/sqrt_Cd*dq              !
         !!
         L = (U_star*U_star) &                ! Estimate the Monin-Obukov length at height zu
              & / (kappa*grav/T_vpot*(T_star*(1.+0.608*q_zu) + 0.608*T_zu*q_star))
         !! Stability parameters :
         zeta_u  = zu/L  ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
         zeta_t  = zt/L  ;  zeta_t = sign( min(abs(zeta_t),10.0), zeta_t )
         zpsi_hu = psi_h(zeta_u)
         zpsi_ht = psi_h(zeta_t)
         zpsi_m  = psi_m(zeta_u)
         !!
         !! Shifting the wind speed to 10m and neutral stability : (L & Y eq.(9a))
!        U_n10 = dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u)))
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is now imposed.
         U_n10 = max(0.25 , dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m)))
         !!
         !! Shifting temperature and humidity at zu :          (L & Y eq. (9b-9c))
!        T_zu = T_zt - T_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
         T_zu = T_zt - T_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
!        q_zu = q_zt - q_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
         q_zu = q_zt - q_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
         !!
         !! q_zu cannot have a negative value : forcing 0
         stab = 0.5 + sign(0.5_wp,q_zu) ;  q_zu = stab*q_zu
         !!
         !! Updating the neutral 10m transfer coefficients :
         Cd_n10  = 1E-3 * (2.7/U_n10 + 0.142 + U_n10/13.09)    ! L & Y eq. (6a)
         sqrt_Cd_n10 = sqrt(Cd_n10)
         Ce_n10  = 1E-3 * (34.6 * sqrt_Cd_n10)                 ! L & Y eq. (6b)
         stab    = 0.5 + sign(0.5_wp,zeta_u)
         Ch_n10  = 1E-3*sqrt_Cd_n10*(18.*stab + 32.7*(1-stab)) ! L & Y eq. (6c-6d)
         !!
         !!
         !! Shifting the neutral 10m transfer coefficients to (zu,zeta_u) :
!        xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u))
         xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m)
         Cd = Cd_n10/(xct*xct) ; sqrt_Cd = sqrt(Cd)
         !!
!        xlogt = log(zu/10.) - psi_h(zeta_u)
         xlogt = log(zu/10.) - zpsi_hu
         !!
         xct = 1. + Ch_n10*xlogt/kappa/sqrt_Cd_n10
         Ch  = Ch_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
         xct = 1. + Ce_n10*xlogt/kappa/sqrt_Cd_n10
         Ce  = Ce_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
         !!
      end do
      !!

      !
   END SUBROUTINE core_coeff_2z

   FUNCTION psi_h( zta )
      !! Psis, L & Y eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      real(wp)             :: X2
      real(wp)             :: X
      real(wp)             :: stabit
      !
      real(wp), intent(in) ::   zta
      real(wp)             ::   psi_h
      !-------------------------------------------------------------------------------

      X2 = sqrt(abs(1. - 16.*zta))  ;  X2 = max(X2 , 1.) ;  X  = sqrt(X2)
      stabit    = 0.5 + sign(0.5_wp,zta)
      psi_h = -5.*zta*stabit  &                                       ! Stable
         &    + (1. - stabit)*(2.*log( (1. + X2)/2. ))                 ! Unstable
   END FUNCTION psi_h

   FUNCTION psi_m( zta )
   !! Psis, L & Y eq. (8c), (8d), (8e)
   !-------------------------------------------------------------------------------
      real(wp)             :: X2
      real(wp)             :: X
      real(wp)             :: stabit
      !!
      real(wp), intent(in) ::   zta

      real(wp), parameter  :: pi = 3.141592653589793_wp
      real(wp)             :: psi_m
      !-------------------------------------------------------------------------------

      X2 = sqrt(abs(1. - 16.*zta))  ;  X2 = max(X2 , 1.0) ;  X  = sqrt(X2)
      stabit    = 0.5 + sign(0.5_wp,zta)
      psi_m = -5.*zta*stabit  &                                                          ! Stable
         &    + (1. - stabit)*(2*log((1. + X)/2) + log((1. + X2)/2) - 2*atan(X) + pi/2)  ! Unstable

      !
   END FUNCTION psi_m

   SUBROUTINE getcoeffld(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE getcoeffld ***
      !!
      !! ** Purpose : read fields from files, interpolate on model mesh and prepare interpolation coefficients
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
!      integer,intent(in)   :: idate ! initialization date
!      integer,intent(in)   :: isec ! initialization seconds
      real(wp),intent(in)   :: rdate ! initialization date

      integer              :: iost !I/O status
      integer              :: ncid      ! netcdf file id
      ! ID dimensions and variables:
      integer              :: id_data
      integer              :: nf_start(4)
      integer              :: nf_edges(4)
!      integer              :: zero_year,yyyy,mm,dd
!      character(len = 256) :: att_string ! attribute
      integer              :: fld_idx, i,j,ii, ip1, jp1, extrp
      integer              :: sbc_alloc, itot

      real(wp)             :: denom, x1, x2, y1, y2, x, y
      real(wp)             :: now_date

      real(wp), allocatable, dimension(:,:)  :: sbcdata1,sbcdata2
      real(wp)             :: data1,data2
      real(wp)             :: delta_t   ! time(t_indx) - time(t_indx+1)

      integer              :: elnodes(4) !4 nodes from one element
      integer              :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )

      ALLOCATE( sbcdata1(nc_Nlon,nc_Nlat), sbcdata2(nc_Nlon,nc_Nlat),&
                &      STAT=sbc_alloc )
!                data1(elem2D),data2(elem2D), &
      if( sbc_alloc /= 0 )   STOP 'getcoeffld: failed to allocate arrays'

      ! find time index in files
      now_date = rdate
      call binarysearch(nc_Ntime,nc_time,now_date,t_indx)
      if ( t_indx < nc_Ntime ) then
         t_indx_p1 = t_indx + 1
         delta_t = nc_time(t_indx_p1) - nc_time(t_indx)
      else ! NO extrapolation to future
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
         WRITE(*,*) '     WARNING:  time lager than last time step in forcing file,'
         WRITE(*,*) '        last time step will be used as a constant field'
         one_field = .true.
      end if

      if ( now_date < nc_time(1) ) then  ! NO extrapolation back in time
         t_indx = 1
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
      end if
      itot = i_totfl
            if (nm_calc_flux==0) then
         itot = 8
            end if

      do fld_idx = 1, itot
         !open file sbc_flfi
         iost = nf_open(sbc_flfi(fld_idx)%file_name,NF_NOWRITE,ncid)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         ! get variable id
         iost = nf_inq_varid(ncid, sbc_flfi(fld_idx)%var_name, id_data)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         !read data from file
         nf_start(1)=1
         nf_edges(1)=nc_Nlon
         nf_start(2)=1
         nf_edges(2)=nc_Nlat
         nf_start(3)=t_indx
         nf_edges(3)=1
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata1)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         ! read next time step in file (check for +1 done before)
         nf_start(3)=t_indx_p1
         nf_edges(3)=1
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata2)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         !flip data in case of lat from -90 to 90
!!!! WARNING
         if ( flip_lat == 1 ) then
             sbcdata1=sbcdata1(:,nc_Nlat:1:-1)
             sbcdata2=sbcdata2(:,nc_Nlat:1:-1)
         end if


         ! bilinear space interpolation, and time interpolation ,
         ! data is assumed to be sampled on a regular grid
!!$OMP PARALLEL
!!$OMP DO
         do ii = 1, myDim_nod2D+eDim_nod2D
            i = bilin_indx_i(ii)
            j = bilin_indx_j(ii)
            ip1 = i + 1
            jp1 = j + 1

!!WARNING !! get coordinates of element (mean of nodes in elem)
!            elnodes = elem2D_nodes(:,ii) !! 4 nodes in element
!            numnodes = 4
!            if(elnodes(1)==elnodes(4)) numnodes = 3  !! set to 3 if we have triangle
!            y = sum(coord_nod2d(2,elnodes(1:numnodes)))/dble(numnodes) /rad
!            x = sum(coord_nod2d(1,elnodes(1:numnodes)))/dble(numnodes) /rad
!  use these lines in case if we use sbc on nodes
            x  = coord_nod2d(1,ii)/rad
            y  = coord_nod2d(2,ii)/rad
            extrp = 0
            if ( i == 0 ) then
               i   = nc_Nlon
               ip1 = i
               extrp = extrp + 1
            end if
            if ( i == -1 ) then
               i   = 1
               ip1 = i
               extrp = extrp + 1
            end if
            if ( j == 0 ) then
               j   = nc_Nlat
               jp1 = j
               extrp = extrp + 2
            end if
            if ( j == -1 ) then
               j   = 1
               jp1 = j
               extrp = extrp + 2
            end if

            x1 = nc_lon(i)
            x2 = nc_lon(ip1)
            y1 = nc_lat(j)
            y2 = nc_lat(jp1)

            if ( extrp == 0 ) then
            ! if point inside forcing domain
               denom = (x2 - x1)*(y2 - y1)
               data1 = ( sbcdata1(i,j)   * (x2-x)*(y2-y)   + sbcdata1(ip1,j)    * (x-x1)*(y2-y) + &
                     sbcdata1(i,jp1) * (x2-x)*(y-y1)   + sbcdata1(ip1, jp1) * (x-x1)*(y-y1)     ) / denom
               data2 = ( sbcdata2(i,j)   * (x2-x)*(y2-y)   + sbcdata2(ip1,j)    * (x-x1)*(y2-y) + &
                     sbcdata2(i,jp1) * (x2-x)*(y-y1)   + sbcdata2(ip1, jp1) * (x-x1)*(y-y1)     ) / denom
            else if ( extrp == 1 ) then !  "extrapolation" in x direction
               denom = (y2 - y1)
               data1 = ( sbcdata1(i,j)   * (y2-y)   + sbcdata1(ip1, jp1) * (y-y1) ) / denom
               data2 = ( sbcdata2(i,j)   * (y2-y)   + sbcdata2(ip1, jp1) * (y-y1) ) / denom
            else if ( extrp == 2 ) then !  "extrapolation" in y direction
               denom = (x2 - x1)
               data1 = ( sbcdata1(i,j)   * (x2-x)   + sbcdata1(ip1, jp1) * (x-x1) ) / denom
               data2 = ( sbcdata2(i,j)   * (x2-x)   + sbcdata2(ip1, jp1) * (x-x1) ) / denom
            else if ( extrp == 3 ) then !  "extrapolation" in x and y direction
               data1 = sbcdata1(i,j)
               data2 = sbcdata2(i,j)
            end if
            ! calculate new coefficients for interpolations
            coef_a(fld_idx, ii) = ( data2 - data1 ) / delta_t !( nc_time(t_indx+1) - nc_time(t_indx) )
            coef_b(fld_idx, ii) = data1 - coef_a(fld_idx, ii) * nc_time(t_indx)

         end do !ii
!!$OMP END DO
!!$OMP END PARALLEL
         iost = nf_close(ncid)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)

      end do !fld_idx

      DEALLOCATE( sbcdata1, sbcdata2 )


   END SUBROUTINE getcoeffld

   SUBROUTINE data_timeinterp(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE data_timeinterp ***
      !!
      !! ** Purpose : interpolation of fields(interpolated on model grid) from IN files in time
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
!      integer,intent(in) :: ndate ! date
!      integer,intent(in) :: nsec  ! seconds
      real(wp),intent(in)    :: rdate  ! seconds

     ! assign data from interpolation to taux and tauy
      integer            :: fld_idx, i,j,ii
      real(wp)           :: now_date

!      now_date = ndate+nsec/86400.0_wp
      now_date = rdate
!!$OMP PARALLEL
!!$OMP DO
      do i = 1, myDim_nod2D+eDim_nod2D
         do fld_idx = 1, i_totfl

            atmdata(fld_idx,i) = now_date * coef_a(fld_idx,i) + coef_b(fld_idx,i)

         end do !fld_idx
      end do !elem2D
!!$OMP END DO
!!$OMP END PARALLEL
   END SUBROUTINE data_timeinterp

   SUBROUTINE sbc_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ini ***
      !!
      !! ** Purpose : inizialization of ocean forcing
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      integer            :: idate ! initialization date
      integer            :: isec  ! initialization seconds
      integer            :: iost  ! I/O status
      integer            :: sbc_alloc                   !: allocation status

      real(wp)           :: tx, ty

      integer            :: node_size


      namelist/nam_sbc/ nm_sbc, nm_tauwind, nm_xwind0, nm_ywind0, nm_sbc_ftype, &
                        nm_xwind_file, nm_ywind_file, nm_humi_file, nm_qsr_file, &
                        nm_qlw_file, nm_tair_file, nm_prec_file, nm_humi0, nm_qsr0, &
                        nm_qlw0, nm_tair, nm_prec,nm_xwind_var,nm_ywind_var,nm_humi_var, &
                        nm_qsr_var, nm_qlw_var, nm_tair_var, nm_prec_var, &
                        nm_mslp_var, nm_mslp_file, nm_mslp, &
                        nm_cloud_var, nm_cloud_file, nm_cloud, &
                        depth_swr, nm_prec_coef, nm_net_flux, nm_calc_flux, &
                        nm_nc_secstep, nm_nc_iyear, nm_nc_imm, nm_nc_idd

      write(*,*) "Start: Ocean forcing inizialization."

      node_size = myDim_nod2D+eDim_nod2D

      idate = time_jd0
      isec  = 86400.0 *(time_jd0 - idate)
! file for debuging
!      open(unit=unit_deb,file='obc_deb.dat')

      ! OPEN and read namelist for SBC
      open( unit=nm_sbc_unit, file='namelist_bc.nml', form='formatted', access='sequential', status='old', iostat=iost )
      if( iost == 0 ) then
         WRITE(*,*) '     file   : ', 'namelist_bc.nml',' open ok'
      else
         WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist_bc.nml',' ; iostat=',iost
         STOP 'ERROR: --> sbc_ini'
      endif
      READ( nm_sbc_unit, nml=nam_sbc, iostat=iost )
      close( nm_sbc_unit )
      if( nm_sbc == -1 ) then
         ! IF module not in use
         return
      endif
      write(*,*) "Start: Ocean forcing inizialization."
      write(*,*) "   Surface boundary conditions parameters:"
      write(*,*) "      nm_sbc        = ", nm_sbc,"   ! 1= constant, 2=from file, -1=module not in use"
      write(*,*) "      nm_tauwind    = ", nm_tauwind ," ! if 1 = wind stress, 2 = wind 10 m  .    "
      write(*,*) "      nm_sbc_ftype  = ", nm_sbc_ftype ," ! input file type: 1 = ASCII, 2 = netcdf   "
      if (nm_sbc==1) then
         write(*,*) "      nm_xwind0     = ", nm_xwind0 ," ! constant 10m. wind value in i-direction/stress, if nm_sbc=1   "
         write(*,*) "      nm_ywind0     = ", nm_ywind0 ," ! constant 10m. wind value in j-direction, if nm_sbc=1  "
         write(*,*) "      nm_humi0      = ", nm_humi0 ," ! constant humidity "
         write(*,*) "      nm_qsr0       = ", nm_qsr0 ," ! constant solar heat  "
         write(*,*) "      nm_qlw0       = ", nm_qlw0 ," ! constant Long wave "
         write(*,*) "      nm_tair       = ", nm_tair," ! constant 2m air temperature "
         write(*,*) "      nm_prec       = ", nm_prec ," ! constant total precipitation "
         write(*,*) " WARNING:: only wind constants in use, rest set to 0 !!!    "
      else
         write(*,*) "      nm_xwind_file = ", trim(nm_xwind_file) ," ! name of file with winds, if nm_sbc=2 "
         write(*,*) "      nm_ywind_file = ", trim(nm_ywind_file) ," ! name of file with winds, if nm_sbc=2 "
         write(*,*) "      nm_humi_file  = ", trim(nm_humi_file) ," ! name of file with humidity "
         write(*,*) "      nm_qsr_file   = ", trim(nm_qsr_file) ," ! name of file with solar heat "
         write(*,*) "      nm_qlw_file   = ", trim(nm_qlw_file) ," ! name of file with Long wave "
         write(*,*) "      nm_tair_file  = ", trim(nm_tair_file) ," ! name of file with 2m air temperature "
         write(*,*) "      nm_prec_file  = ", trim(nm_prec_file) ," ! name of file with total precipitation "
         write(*,*) "      nm_mslp_file  = ", trim(nm_mslp_file)," !air_pressure_at_sea_level "
         write(*,*) "      nm_cloud_file  = ", trim(nm_cloud_file)," !clouds "
         write(*,*) "      nm_xwind_var  = ", trim(nm_xwind_var) ," ! name of variable in file with wind "
         write(*,*) "      nm_ywind_var  = ", trim(nm_ywind_var) ," ! name of variable in file with wind "
         write(*,*) "      nm_humi_var   = ", trim(nm_humi_var) ," ! name of variable in file with humidity  "
         write(*,*) "      nm_qsr_var    = ", trim(nm_qsr_var) ," ! name of variable in file with solar heat "
         write(*,*) "      nm_qlw_var    = ", trim(nm_qlw_var) ," ! name of variable in file with Long wave "
         write(*,*) "      nm_tair_var   = ", trim(nm_tair_var) ," ! name of variable in file with 2m air temperature "
         write(*,*) "      nm_prec_var   = ", trim(nm_prec_var) ," ! name of variable in file with total precipitation  "
         write(*,*) "      nm_mslp_var   = ", trim(nm_mslp_var) ," ! name of variable in file with air_pressure_at_sea_level "
         write(*,*) "      nm_cloud_var   = ", trim(nm_cloud_var) ," ! name of variable in file with clouds "
      endif
      write(*,*) "      depth_swr     = ", depth_swr ," !  depth of swr penetration"
      write(*,*) "      nm_net_flux   = ", nm_net_flux," !  key for downward longwave heat over the ocean: 0 - , 1 - Net"
      write(*,*) "      nm_prec_coef  = ", nm_prec_coef," !  precipitation will be devide by this constant(3600-CoastDat,1-NCEP)"
      write(*,*) "      nm_nc_secstep = ", nm_nc_secstep ," !  time units coef (86400 CoastDat, 24 NCEP)"
      write(*,*) "      nm_nc_iyear   = ", nm_nc_iyear ," ! initial year of time axis in netCDF (1948 like CoastDat,1800 NCEP)"
      write(*,*) "      nm_nc_imm     = ", nm_nc_imm ," ! initial month of time axis in netCDF "
      write(*,*) "      nm_nc_idd     = ", nm_nc_idd ," ! initial day of time axis in netCDF "
      write(*,*) "      nm_calc_flux  = ", nm_calc_flux ,"! 1= calculate I0,... (GOTM, 0= use atm. model fluxes"
      ALLOCATE( coef_a(i_totfl,node_size), coef_b(i_totfl,node_size), &
              & atmdata(i_totfl,node_size), &
                   &      STAT=sbc_alloc )
      if( sbc_alloc /= 0 )   STOP 'sbc_ini: failed to allocate arrays'



      ALLOCATE( bilin_indx_i(node_size),bilin_indx_j(node_size), &
              & qns(node_size), emp(node_size), qsr(node_size),  &
                   &      STAT=sbc_alloc )
!              & qns_2(nod2D), emp_2(nod2D), qsr_2(nod2D),  &
!              & taux_node_2(nod2D), tauy_node_2(nod2D),  &

      if( sbc_alloc /= 0 )   STOP 'sbc_ini: failed to allocate arrays'

! IF constant wind/stress
      if( nm_sbc == 1 ) then
         STOP "ERROR: sbc_ini: nm_sbc == 1, constant wind/stress is under developing "

         if( nm_tauwind == 1 ) then
            write(*,*) "WARNING:: nm_sbc == 1 and  nm_tauwind == 1, "
            write(*,*) "           to calculate wind at 10m will be used simple Large and Pond (1981) _revers_ parametrization."
            call stress2wind_scalar(nm_xwind0, nm_ywind0, tx, ty)
            taux = nm_xwind0
            tauy = nm_ywind0
            taux_node = nm_xwind0
            tauy_node = nm_ywind0

         else
            write(*,*) "WARNING:: nm_sbc == 1 and  nm_tauwind == 2, "
            write(*,*) "           to calculate wind stress will be used simple Large and Pond (1981) parametrization."
!            STOP "ERROR: sbc_ini"
            call wind2stress_scalar(nm_xwind0, nm_ywind0, tx, ty)
            taux = tx
            tauy = ty
            taux_node = tx
            tauy_node = ty
            windx = nm_xwind0
            windy = nm_ywind0
         endif
      endif


! ========================== ASCII PART start===============================================================
! OPEN ASCII file with winds if needed and check number of nodes/elements (elem2D)
      if( nm_sbc == 2 .AND. nm_sbc_ftype == 1 ) then
         if( nm_tauwind == 2 ) then
            write(*,*) "This function not implemented yet.nm_sbc == 2  nm_sbc_ftype == 1 nm_tauwind == 2 "
            STOP "ERROR: sbc_ini"
         endif

         call ascii_sbc_ini(idate,isec)

      endif
! ========================== NC PART start===============================================================
      if( nm_sbc == 2 .AND. nm_sbc_ftype == 2 ) then
         call nc_sbc_ini(idate,isec)
      endif


      write(*,*) "DONE:  Ocean forcing inizialization."
   END SUBROUTINE sbc_ini

   SUBROUTINE data2tau_interp(ndate,nsec)
      IMPLICIT NONE
      integer,intent(in) :: ndate ! date
      integer,intent(in) :: nsec  ! seconds
     ! assign data from interpolation to taux and tauy
      integer           :: i                    ! loop variables

      if( nm_tauwind == 1 ) then
         do i = 1, nod2D
            taux_node(i) = (ndate+nsec/86400.0_wp) * coef_a(i_xwind,i) + coef_b(i_xwind,i)
            tauy_node(i) = (ndate+nsec/86400.0_wp) * coef_a(i_ywind,i) + coef_b(i_ywind,i)
         enddo
         do i = 1, elem2D
            taux(i) = sum(w_cv(1:4,i)*taux_node(elem2D_nodes(:,i)))
            tauy(i) = sum(w_cv(1:4,i)*tauy_node(elem2D_nodes(:,i)))
         end do
      else
         do i = 1, nod2D
            STOP "not implemented data2tau_interp"
!            taux(i) = wind2stress_scalar((ndate+nsec/86400.0_wp) * coef_a(i_xwind,i) + coef_b(i_xwind,i))
!            tauy(i) = wind2stress_scalar((ndate+nsec/86400.0_wp) * coef_a(i_ywind,i) + coef_b(i_ywind,i))
         enddo
      endif
   END SUBROUTINE data2tau_interp

   SUBROUTINE data2tau_const
     ! assign data from file to taux and tauy
      IMPLICIT NONE
      integer           :: i                    ! loop variables

      if( nm_tauwind == 1 ) then
!         taux = datawx_f
!         tauy = datawy_f
         taux_node = datawx_f
         tauy_node = datawy_f
         do i = 1, elem2D
            taux(i) = sum(w_cv(1:4,i)*taux_node(elem2D_nodes(:,i)))
            tauy(i) = sum(w_cv(1:4,i)*tauy_node(elem2D_nodes(:,i)))
         end do
      else
         do i = 1, nod2D
            STOP "not implemented data2tau_interp"
!            taux(i) = wind2stress_scalar(datawx_f(i))
!            tauy(i) = wind2stress_scalar(datawy_f(i))
         enddo
      endif
   END SUBROUTINE data2tau_const

   SUBROUTINE sbc_do
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_do ***
      !!
      !! ** Purpose : provide at each time-step: wind stress, ...
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer            :: rdate ! date
      integer            :: rsec  ! seconds

      if( nm_sbc == -1 ) then
         ! IF module not in use
         return
      endif

!     currently not in use
!      if( nm_sbc == 1 ) then
         ! IF constant values were used
         !assuming that taux and tauy were not changed during previews time step after initialization
!         return
!      endif

      rdate = time_jd
      rsec  = 86400.0 *(time_jd - rdate)

      if( nm_sbc_ftype == 1 ) then
         if( one_field ) then
            ! IF one time step in IN files
            !assuming that taux and tauy were not changed during previews time step after initialization
            return
         endif

         ! ASCII case
         if( rdate+rsec/86400.0_wp > datewx_s + secwx_s/86400.0_wp ) then
            call sbc_do_ascii(rdate,rsec)
         end if

        !put new (interpolated on a new date) values to taux and tauy

         call data2tau_interp(rdate,rsec)

      else
         ! NetCDF case

         if( .not. one_field ) then
            ! IF more field available
            if( time_jd > nc_time(t_indx_p1) ) then
               ! get new coefficients for time interpolation on model grid for all data
!write(*,*) 'ini0'
               call getcoeffld(time_jd)
            endif
         endif
         ! interpolate in time
         call data_timeinterp(time_jd)
         ! change model ocean parameters
         if (nm_calc_flux==1) then
!write(*,*) 'ini01'
            call apply_atm_fluxes_gotm(time_jd)
         else
!write(*,*) 'ini1'
            call apply_atm_fluxes
         end if


      endif




   END SUBROUTINE sbc_do


   SUBROUTINE sbc_do_ascii(rdate,rsec)
   ! check if it is time to change coeficients (ASCII case)
      integer,intent(in) :: rdate ! date
      integer,intent(in) :: rsec  ! seconds

      integer :: iost, iost2 ! I/O status

         !if rdate later than second date in forcing
         ! read next forcing steps while rdate will not be early than second date in forcing or no more data in file
         iost = 0
         iost2= 0
         do while ( (rdate+rsec/86400.0_wp) >= (datewx_s + secwx_s/86400.0_wp) .and. iost == 0 .and. iost2 == 0 )
            datawx_f = datawx_s
            datawy_f = datawy_s
            datewx_f = datewx_s
            datewy_f = datewy_s
            secwx_f  = secwx_s
            secwy_f  = secwy_s

            call read_sbc_ascii(xwind_file_unit,nod2D,datawx_s,datewx_s,secwx_s,iost)
            call read_sbc_ascii(ywind_file_unit,nod2D,datawy_s,datewy_s,secwy_s,iost2)

         enddo
         ! check if reach end of file than use previews field as a constant
         if( (iost == END_OF_FILE) .or. (iost2 == END_OF_FILE) ) then
!            WRITE(*,*) '     WARNING: sbc_ini: current date lager than last time step in forcing file,'
!            WRITE(*,*) '        last time step will be used as a constant field'
            one_field = .true.
            coef_b(i_xwind,:) = datawx_f
            coef_a(i_xwind,:) = 0.0_wp
            coef_b(i_ywind,:) = datawy_f
            coef_a(i_ywind,:) = 0.0_wp
         else
            if( iost /= 0 ) call err_call(iost,nm_xwind_file)
         endif
         if( .not. one_field ) then
         ! prepare coeficients for liniar interpolations
            ! check if dates as not descending
            if( datewx_s+secwx_s/86400.0_wp <= datewx_f + secwx_f/86400.0_wp ) then
               WRITE(*,*) 'ERROR sbc_do: descending dates in file: ', nm_xwind_file
               STOP 'ERROR: sbc_do'
            endif
            if( datewy_s+secwy_s/86400.0_wp <= datewy_f + secwy_f/86400.0_wp ) then
               WRITE(*,*) 'ERROR sbc_do: descending dates in file: ', nm_ywind_file
               STOP 'ERROR: sbc_do'
            endif

            ! calculate new coeficients for interpolations
            coef_a(i_xwind,:) = ( datawx_s - datawx_f )/( datewx_s +&
                                  secwx_s/86400.0_wp - datewx_f - secwx_f/86400.0_wp )
            coef_b(i_xwind,:) = datawx_f - coef_a(i_xwind,:) * ( datewx_f +&
                                secwx_f/86400.0_wp )
            coef_a(i_ywind,:) = ( datawy_s - datawy_f )/( datewy_s + &
                                secwy_s/86400.0_wp - datewy_f - secwy_f/86400.0_wp )
            coef_b(i_ywind,:) = datawy_f - coef_a(i_ywind,:) * ( datewy_f + &
                                secwy_f/86400.0_wp )
         endif

   END SUBROUTINE sbc_do_ascii

   SUBROUTINE read_sbc_ascii(funit,nod,dataf,date,sec,ierr)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_do ***
      !!
      !! ** Purpose : read 2d field (vector)
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      integer,intent(in)    :: funit
      integer,intent(in)    :: nod

      integer, intent(out)  :: ierr
      real(wp), intent(out) :: dataf(:) ! data

      integer, intent(out) :: date ! date
      integer, intent(out) :: sec  ! seconds


      integer           :: yy,mm,dd,hh,min,ss   ! time variables
      character         :: c1,c2,c3,c4
      character(len=20) :: cbuf
      integer           :: i                    ! loop variables
      integer           :: aloc_stat            ! loop variables

      real(wp),allocatable,dimension(:) :: tmp_data(:) ! data

!   integer   :: i,j,k  ! loop variables
      allocate(tmp_data(nod),stat=aloc_stat)
      if (aloc_stat /= 0) STOP 'read_sbc_ascii: Error allocating memory (tmp_data)'
      ierr = 0
      read(funit,'(A19)',ERR=100,END=110) cbuf
      read(cbuf,900,ERR=100,END=110) yy,c1,mm,c2,dd,hh,c3,min,c4,ss
      do i = 1, nod
         read(funit,*,ERR=100,END=110) tmp_data(i)
      end do
      ! if all ok, then rewrite output
      sec = hh * 3600 + min * 60 + ss
      date = julday(yy,mm,dd) - 0.5
      dataf = tmp_data
      deallocate(tmp_data)

      return

100 ierr = READ_ERROR
      deallocate(tmp_data)

      return

110 ierr = END_OF_FILE
      deallocate(tmp_data)

      return

900 format(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)


   END SUBROUTINE read_sbc_ascii


   SUBROUTINE wind2stress_scalar(U10x,U10y,tx,ty)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  wind2stress ***
      !!
      !! ** Purpose : convert wind 10m to wind stress (return scalar)
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      real(wp), intent(in)  :: U10x ! wind X direction
      real(wp), intent(in)  :: U10y ! wind Y direction
      real(wp), intent(out) :: tx   ! tau X direction
      real(wp), intent(out) :: ty   ! tau X direction

      real(wp)             :: U10

!     Large and Pond (1981), J. Phys. Oceanog., 11, 324-336.
!     Tau = Cd * ro_air * U10^2
!     U10    = wind speed at 10 m above the sea surface
!     ro_air =  1.22 kg m-3
!     Cd     = dimensionless drag coefficient :
!         1.2 x 10^-3 for 4 < U10 < 11 m s-1
!         10^-3 (0.49 + 0.065 U10)     for     11 < U10 < 25 m s-1
      U10 = sqrt( U10x*U10x + U10y*U10y )
      if ( U10 <= 11 ) then
         tx = 1.2d-3 * 1.22 * U10 * U10x
         ty = 1.2d-3 * 1.22 * U10 * U10y
      else
         tx = 1.d-3 * (0.49 + 0.065 * U10) * 1.22 * U10 * U10x
         ty = 1.d-3 * (0.49 + 0.065 * U10) * 1.22 * U10 * U10y
      endif

   END SUBROUTINE wind2stress_scalar

   SUBROUTINE stress2wind_scalar(tx, ty, U10x, U10y)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  wind2stress ***
      !!
      !! ** Purpose : convert wind stress to wind 10m (return scalar)
      !! ** Method  :  it is very bad way to do it, sorry
      !! ** Action  :
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      real(wp), intent(out)  :: U10x ! wind X direction
      real(wp), intent(out)  :: U10y ! wind Y direction
      real(wp), intent(in) :: tx   ! tau X direction
      real(wp), intent(in) :: ty   ! tau X direction

!     Tau = Cd * ro_air * U10^2
!     U10    = wind speed at 10 m above the sea surface
!     ro_air =  1.22 kg m-3
!     Cd     = dimensionless drag coefficient : 1.2 x 10^-3
!

      U10y = sqrt(ty / (1.22 * 1.2d-3 * sqrt(tx*tx/(ty*ty)+1) ))
      U10x = sqrt(tx / (1.22 * 1.2d-3 * sqrt(ty*ty/(tx*tx)+1) ))


   END SUBROUTINE stress2wind_scalar

   SUBROUTINE err_call(iost,fname)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  err_call ***
      !!
      !! ** Purpose : call Error
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      integer, intent(in)            :: iost
      character(len=256), intent(in) :: fname
      write(*,*) 'ERROR: I/O status=',iost,' file= ',fname
      STOP 'ERROR:  stop'


   END SUBROUTINE err_call

   FUNCTION julday(yyyy,mm,dd)

   IMPLICIT NONE
      integer, INTENT(IN) :: mm, dd, yyyy
      integer             :: julday
! In this routine julday returns the Julian Day Number that begins at noon of the calendar
!    date specified by month mm , day dd , and year yyyy , all integer variables. Positive year
!    signifies A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D. (from Num. Rec.)
      integer, PARAMETER  :: IGREG=15+31*(10+12*1582)
! Gregorian Calendar adopted Oct. 15, 1582.
      integer             :: ja,jm,jy

      jy = yyyy
      if (jy == 0) STOP 'julday: there is no year zero'
      if (jy < 0) jy=jy+1
      if (mm > 2) then
         jm=mm+1
      else
         jy=jy-1
         jm=mm+13
      endif
      julday=int(365.25_wp*jy)+int(30.6001_wp*jm)+dd+1720995
!Test whether to change to Gregorian Calendar.
      if (dd+31*(mm+12*yyyy) >= IGREG) then
         ja=int(0.01*jy)
         julday=julday+2-ja+int(0.25_wp*ja)
      end if
   END FUNCTION julday


   SUBROUTINE calendar_date(julian,yyyy,mm,dd)

!  Converts a Julian day to a calendar date (year, month and day). Numerical Recipes
   IMPLICIT NONE
!
      integer,intent(in) :: julian
      integer            :: yyyy,mm,dd

      integer, parameter :: IGREG=2299161
      integer            :: ja,jb,jc,jd,je
      real(wp)           :: x
!
!-----------------------------------------------------------------------
      if (julian >= IGREG ) then
         x = ((julian-1867216)-0.25)/36524.25
         ja = julian+1+int(x)-int(0.25*x)
      else
         ja = julian
      end if

      jb = ja+1524
      jc = int(6680 + ((jb-2439870)-122.1)/365.25)
      jd = int(365*jc+(0.25*jc))
      je = int((jb-jd)/30.6001)

      dd = jb-jd-int(30.6001*je)
      mm = je-1
      if (mm > 12) mm = mm-12
      yyyy = jc - 4715
      if (mm > 2) yyyy = yyyy-1
      if (yyyy <= 0) yyyy = yyyy-1

      return
   END SUBROUTINE calendar_date

   SUBROUTINE ascii_sbc_ini(idate,isec)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ascii_ini ***
      !!
      !! ** Purpose : inizialization of ocean forcing from ASCII file
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      integer,intent(in) :: idate ! initialization date
      integer,intent(in) :: isec ! initialization seconds
      integer :: iost, iost2 ! I/O status
      integer :: tmpi ! temporary integer
      integer :: i    ! loop variable
      integer :: sbc_alloc

! for X direction
         ALLOCATE( datawx_f(nod2D), datawx_s(nod2D), &
               & datawy_f(nod2D), datawy_s(nod2D), &
                     &      STAT=sbc_alloc )
         if( sbc_alloc /= 0 )   STOP 'sbc_ini: failed to allocate arrays'

         open( unit=xwind_file_unit, FILE=nm_xwind_file, status='old', iostat=iost )
         if( iost == 0 ) then
            WRITE(*,*) '     file   : ', trim(nm_xwind_file),' open ok'
         else
            WRITE(*,*) 'ERROR: --> bad opening file   : ', trim(nm_xwind_file),' ; iostat=',iost
            STOP 'ERROR: --> sbc_ini'
         endif

         READ( unit=xwind_file_unit, FMT=*,iostat=iost) tmpi
         if( tmpi /= nod2D ) then
            WRITE(*,*) 'ERROR: --> file   : ', trim(nm_xwind_file),' ; wrong nod2D'
            STOP 'ERROR: --> sbc_ini'
         endif
! for Y direction

         open( unit=ywind_file_unit, FILE=nm_ywind_file, status='old', iostat=iost )
         if( iost == 0 ) then
            WRITE(*,*) '     file   : ', trim(nm_ywind_file),' open ok'
         else
            WRITE(*,*) 'ERROR: --> bad opening file   : ', trim(nm_ywind_file),' ; iostat=',iost
            STOP 'ERROR: --> sbc_ini'
         endif

         READ( unit=ywind_file_unit, FMT=*,iostat=iost) tmpi
         if( tmpi /= nod2D ) then
            WRITE(*,*) 'ERROR: --> file   : ', trim(nm_ywind_file),' ; wrong nod2D'
            STOP 'ERROR: --> sbc_ini'
         endif

!read first vectors
         call read_sbc_ascii(xwind_file_unit,nod2D,datawx_f,datewx_f,secwx_f,iost)
         if( iost /= 0 ) call err_call(iost,nm_xwind_file)
         call read_sbc_ascii(ywind_file_unit,nod2D,datawy_f,datewy_f,secwy_f,iost)
         if( iost /= 0 ) call err_call(iost,nm_ywind_file)
!try to read second time step in file, if it does not exist than use first one as not changed winf during simulation
         call read_sbc_ascii(xwind_file_unit,nod2D,datawx_s,datewx_s,secwx_s,iost)
         if( iost == END_OF_FILE ) then
            WRITE(*,*) '     WARNING: sbc_ini: found only one time step in wind file'
            one_field = .true.
         else
            if( iost /= 0 ) call err_call(iost,nm_xwind_file)
         endif
         if( datewx_s+secwx_s/86400.0_wp <= datewx_f + secwx_f/86400.0_wp ) then
            WRITE(*,*) 'ERROR sbc_ini: descending dates in file: ', nm_xwind_file
            STOP 'ERROR: sbc_ini'
         endif

         call read_sbc_ascii(ywind_file_unit,nod2D,datawy_s,datewy_s,secwy_s,iost)
         if( iost == END_OF_FILE ) then
            WRITE(*,*) '     WARNING: sbc_ini: found only one time step in wind file'
            one_field = .true.
         else
            if( iost /= 0 ) call err_call(iost,nm_ywind_file)
         endif
         if( datewy_s+secwy_s/86400.0_wp <= datewy_f + secwy_f/86400.0_wp ) then
            WRITE(*,*) 'ERROR sbc_ini: descending dates in file: ', nm_ywind_file
            STOP 'ERROR: sbc_ini'
         endif


         if( one_field ) then
! IF constant wind/stress field from file
               call data2tau_const
         else

! find proper time in file according to ini date,make first interpolation
! if proper time not found use last or first known data
            if( idate+isec/86400.0_wp <= datewx_f + secwx_f/86400.0_wp ) then
            ! if initial date early than first date in forcing
               coef_b(i_xwind,:) = datawx_f
               coef_a(i_xwind,:) = 0.0_wp
               coef_b(i_ywind,:) = datawy_f
               coef_a(i_ywind,:) = 0.0_wp

               call data2tau_interp(idate,isec)
               WRITE(*,*) '     WARNING: sbc_ini: initial time step early than first time step in forcing file,'
               WRITE(*,*) '        first time step will be used as a constant field'
! rewind data files and put first date from file to initial date, same with data
               rewind( unit=xwind_file_unit, iostat=iost )
               read( unit=xwind_file_unit, FMT=*,iostat=iost) tmpi
               call read_sbc_ascii(xwind_file_unit,nod2D,datawx_s,datewx_s,secwx_s,iost)
               datawx_f=datawx_s
               datewx_f=idate
               datewx_f=isec
               rewind( unit=ywind_file_unit, iostat=iost )
               read( unit=ywind_file_unit, FMT=*,iostat=iost) tmpi
               call read_sbc_ascii(ywind_file_unit,nod2D,datawy_s,datewy_s,secwy_s,iost)
               datawy_f=datawy_s
               datewy_f=idate
               datewy_f=isec



            endif
            if( idate+isec/86400.0_wp > datewx_f + secwx_f/86400.0_wp ) then
            ! if initial date later than first date in forcing
            ! reread forcing while initial date will not be early than second date in forcing or no more data in file
               iost = 0
               iost2= 0
               do while ( (idate+isec/86400.0_wp) > (datewx_s + secwx_s/86400.0_wp) .and. iost == 0 .and. iost2 == 0 )
                  datawx_f = datawx_s
                  datawy_f = datawy_s
                  datewx_f = datewx_s
                  datewy_f = datewy_s
                  secwx_f  = secwx_s
                  secwy_f  = secwy_s

                  call read_sbc_ascii(xwind_file_unit,nod2D,datawx_s,datewx_s,secwx_s,iost)
                  call read_sbc_ascii(ywind_file_unit,nod2D,datawy_s,datewy_s,secwy_s,iost2)
               enddo
               if( (iost == END_OF_FILE) .or. (iost2 == END_OF_FILE) ) then
                  WRITE(*,*) '     WARNING: sbc_ini: initial time step lager than last time step in forcing file,'
                  WRITE(*,*) '        last time step will be used as a constant field'
                  one_field = .true.
               else
                  if( iost /= 0 ) call err_call(iost,nm_xwind_file)
               endif
               if( .not. one_field ) then
         ! prepare coefficients for linear interpolations
                  if( datewx_s+secwx_s/86400.0_wp <= datewx_f + secwx_f/86400.0_wp ) then
                     WRITE(*,*) 'ERROR sbc_ini: descending dates in file: ', nm_xwind_file
                     STOP 'ERROR: sbc_ini'
                  endif
                  if( datewy_s+secwy_s/86400.0_wp <= datewy_f + secwy_f/86400.0_wp ) then
                     WRITE(*,*) 'ERROR sbc_ini: descending dates in file: ', nm_ywind_file
                     STOP 'ERROR: sbc_ini'
                  endif

                  coef_a(i_xwind,:) = ( datawx_s - datawx_f )/( datewx_s +&
 secwx_s/86400.0_wp - datewx_f - secwx_f/86400.0_wp )
                  coef_b(i_xwind,:) = datawx_f - coef_a(i_xwind,:) * ( datewx_f +&
 secwx_f/86400.0_wp )
                  coef_a(i_ywind,:) = ( datawy_s - datawy_f )/( datewy_s +&
 secwy_s/86400.0_wp - datewy_f - secwy_f/86400.0_wp )
                  coef_b(i_ywind,:) = datawy_f - coef_a(i_ywind,:) * ( datewy_f +&
 secwy_f/86400.0_wp )

                  call data2tau_interp(idate,isec)
               endif

            endif

         endif


   END SUBROUTINE ascii_sbc_ini

   SUBROUTINE sbc_end

      IMPLICIT NONE

      if( nm_sbc == 2 .AND. nm_sbc_ftype == 1 ) then
         DEALLOCATE( datawx_f, datawx_s, datawy_f, datawy_s )
      endif

      if( nm_sbc == 2 .AND. nm_sbc_ftype == 2 ) then
         DEALLOCATE( nc_lon, nc_lat, nc_time)
      endif

      DEALLOCATE( coef_a, coef_b, atmdata, &
                  &  bilin_indx_i, bilin_indx_j,  &
                  &  qns, emp, qsr)




   END SUBROUTINE sbc_end

   SUBROUTINE check_nferr(iost,fname)
   IMPLICIT NONE
      character(len=256), intent(in) :: fname
      integer, intent(in) :: iost

      if (iost .ne. NF_NOERR) then
         write(*,*) 'ERROR: I/O status= "',trim(nf_strerror(iost)),'";',iost,' file= ',fname
         STOP 'ERROR: stop'
      endif
   END SUBROUTINE

   SUBROUTINE binarysearch(length, array, value, ind)!, delta)
      ! Given an array and a value, returns the index of the element that
      ! is closest to, but less than, the given value.
      ! Uses a binary search algorithm.
      ! "delta" is the tolerance used to determine if two values are equal
      ! if ( abs(x1 - x2) <= delta) then
      !    assume x1 = x2
      ! endif
      !org. source from: https://github.com/cfinch/Shocksolution_Examples/blob/master/FORTRAN/BilinearInterpolation/interpolation.f90

      IMPLICIT NONE
      integer,  intent(in) :: length
      real(wp), dimension(length), intent(in) :: array
      real(wp), intent(in) :: value
   !   real, intent(in), optional :: delta

   !   integer :: binarysearch
      integer, intent(out) :: ind

      integer :: left, middle, right
      real(wp):: d

   !   if (present(delta) .eqv. .true.) then
   !      d = delta
   !   else
      d = 1e-9
   !   endif
      left = 1
      right = length
      do
         if (left > right) then
            exit
         endif
         middle = nint((left+right) / 2.0)
         if ( abs(array(middle) - value) <= d) then
            ind = middle
            return
         else if (array(middle) > value) then
            right = middle - 1
         else
            left = middle + 1
         end if
      end do
      ind = right

   END SUBROUTINE binarysearch

   SUBROUTINE apply_atm_fluxes_gotm(gdate)
      !!----------------------------------------------------------------------
      !! ** Purpose : Change model variables according to atm fluxes
      !! source of original code: GOTM
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      real(wp),intent(in)  :: gdate  ! date

      integer             :: rdate ! date (days)
      integer             :: rsec  ! seconds
      integer             :: yearday, yyyy,dd,mm

      integer             :: i
      integer             :: hum_method = 4
      integer             :: back_radiation_method = 1
!      real(wp)            :: wndm    ! delta of wind module and ocean curent module
      real(wp)            :: wdx,wdy ! delta of wind x/y and ocean curent x/y
      real(wp)            :: tx,ty ! wind stress x/y
      real(wp)            :: zst     ! surface temperature in Kelvin

      real(wp)            :: qa ! specific humidity (kg/kg)
      real(wp)            :: ea ! actual water vapor pressure in Pascal
      real(wp)            :: es ! saturation vapor pressure
      real(wp)            :: qs ! saturation specific humidity
      real(wp)            :: rhoa ! air density

      real(wp)            :: qe  ! turbulent sensible heat flux (W/m2)
      real(wp)            :: qh  ! turbulent latent heat flux (W/m2)
      real(wp)            :: qb  ! long-wave back radiation (W/m2)
      real(wp)            :: evap  ! evaporation/condensation (m/s)
      real(wp)            :: zenith_angle
      real(wp)            :: hh
      real(wp)            :: I_0

      rdate = gdate
      rsec  = 86400.0 *(gdate - rdate)
      hh = rsec*(1.0_wp/3600)

      call calendar_date(rdate,yyyy,dd,mm)
      yearday = rdate-julday(yyyy,1,1)+1-0.5

      do i = 1,nod2D
         windx(i) = atmdata(i_xwind,i)
         windy(i) = atmdata(i_ywind,i)
         ! IF you going to use winds on nodes (nod2d) U_n and V_n should be changed to Unode and Vnode
         if (type_task>1) then
            wdx = atmdata(i_xwind,i) - Unode(1,i) ! wind from data - ocean current ( x direction)
            wdy = atmdata(i_ywind,i) - Vnode(1,i) ! wind from data - ocean current ( y direction)
         else
! WARNING:: surface curent U and V should be used here, while it should be interpolate by compute_el2nodes_2D_2fields from fv_utilit
!           so if u use type_task = 1 and winds you are welcome to add some code
            wdx = atmdata(i_xwind,i)! - U_n_2D(1,i) ! wind from data - ocean current ( x direction)
            wdy = atmdata(i_ywind,i)! - U_n_2D(2,i) ! wind from data - ocean current ( y direction)
         endif
!         wndm = SQRT( wdx * wdx + wdy * wdy )

         ! TF (temperature) should be interpolated on elem + Kelvin  (NCEP/NCAR is in K.)
         ! If no temperature calculation, then constant T
         if (type_task>2) then
!        if u use elem2D
!            zst = sum(w_cv(1:4,i)*TF(1,elem2D_nodes(:,i)))+273.15
            zst = TF(1,i)+273.15
         else
            zst = T_const+273.15
         endif

  !subroutine humidity(hum_method,hum,              airp,             tw, ta,               qa,qs,rhoa,ea,es)
         call humidity(hum_method,atmdata(i_humi,i),atmdata(i_mslp,i),zst,atmdata(i_tair,i),qa,qs,rhoa,ea,es)

  !subroutine back_radiation(method,               dlat,                tw, ta,               cloud,ea,qa,qb)
         call back_radiation(back_radiation_method,coord_nod2d(2,i)/rad,zst,atmdata(i_tair,i),atmdata(i_cloud,i),ea,qa,qb)

  !subroutine fairall(sst,airt,             u10,v10,precip,                        qs,qa,rhoa,evap,taux,tauy,qe,qh)
         call fairall(zst,atmdata(i_tair,i),wdx,wdy,atmdata(i_prec,i)/nm_prec_coef,qs,qa,rhoa,evap,tx,ty,qe,qh)

         zenith_angle = solar_zenith_angle(yearday,hh,coord_nod2d(1,i)/rad,coord_nod2d(2,i)/rad)

!         function short_wave_radiation(zenith_angle,yday,dlon,dlat,cloud)

         I_0 = short_wave_radiation(zenith_angle,yearday,coord_nod2d(1,i)/rad,coord_nod2d(2,i)/rad,atmdata(i_cloud,i))


!!!===========================================  changing of taux and tauy here ====================
         taux_node(i) = tx
         tauy_node(i) = ty
!!!===========================================  changing of taux and tauy here ====================
         emp(i) = evap - atmdata(i_prec,i)/nm_prec_coef
         qns(i) = qb+qe+qh
         qsr(i) = I_0
!!!================================ data for use in ocean part ================================
         ! Downward Non Solar heat flux over the ocean
!         qns(i) = zqlw - zqsb - zqla
         ! zqlw   = downward longwave heat over the ocean
         ! - zqsb = downward sensible heat over the ocean
         ! - zqla = downward latent   heat over the ocean
         ! qns    = downward non solar heat over the ocean
!         emp(i) = zevap - atmdata(i_prec,i)/nm_prec_coef
!         qsr(i) = (1.0_wp - albo ) * atmdata(i_qsr,i)
!         mslp(i)= atmdata(i_mslp,i)
!!!================================ data for use in ocean part ================================
      end do
      ! convert tau to elements
      do i = 1, elem2D
         taux(i) = sum(w_cv(1:4,i)*taux_node(elem2D_nodes(:,i)))
         tauy(i) = sum(w_cv(1:4,i)*tauy_node(elem2D_nodes(:,i)))
      end do

   END SUBROUTINE apply_atm_fluxes_gotm

!-----------------------------------------------------------------------
! This subroutine taken from GOTM src, used to compare with our old fluxes,
!   original code from (see description of code)
!BOP
!
! !ROUTINE: Heat and momentum fluxes according to Fairall et al.
!
! !INTERFACE:
   subroutine fairall(sst,airt,u10,v10,precip,qs,qa,rhoa,evap,taux,tauy,qe,qh)
!
! !DESCRIPTION:
!  The surface momentum flux vector, $(\tau_x^s,\tau_y^s)$,
!  in [N\,m$^{-2}$],
!  the latent heat flux, $Q_e$,
!  and the sensible heat flux, $Q_h$, both in [W\,m$^{-2}$]
!  are calculated here according to the \cite{Fairalletal96a} bulk
!  formulae, which are build on the Liu-Katsaros-Businger
!  (\cite{Liuetal79}) method.
!  Cool skin and warm layer effects are considered according to the
!  suggestions of \cite{Fairalletal96b}.
!
!  The air temperature {\tt airt} and the sea surface temperature
!  {\tt sst} may be given in Kelvin or Celsius:
!  if they are $>$ 100 - Kelvin is assumed.
!
!  This piece of code has been adapted from the COARE code originally
!  written by David Rutgers and Frank Bradley - see
!  http://www.coaps.fsu.edu/COARE/flux\_algor/flux.html.
!
! !USES:
!   use airsea_variables, only: const06,rgas,rho_0,g,rho_0,kappa
!   use airsea_variables, only: qs,qa,rhoa
!   use airsea_variables, only: cpa,cpw
!   use airsea, only: rain_impact,calc_evaporation
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(wp), intent(in)                :: sst,airt,u10,v10,precip
   real(wp), intent(in)                :: qa ! specific humidity (kg/kg)
   real(wp), intent(in)                :: qs ! saturation specific humidity
   real(wp), intent(in)                :: rhoa ! air density
!
! !INPUT/OUTPUT PARAMETERS:
   real(wp), intent(out)               :: evap
!
! !OUTPUT PARAMETERS:
   real(wp), intent(out)               :: taux,tauy,qe,qh
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips
!
! !DEFINED PARAMETERS:
!  Fairall LKB roughness Reynolds number to Von Karman
   real(wp),parameter        :: fdg = 1.0          ! non-dimensional

!  Beta parameter evaluated from Fairall low windspeed turbulence data.
   real(wp),parameter        :: beta = 1.2         ! non-dimensional

!  Zabl      Height (m) of atmospheric boundary layer.
   real(wp),parameter        :: Zabl = 600.0       ! in [m]

   real(wp), parameter       :: r3 = 1.0/3.0
!
!  Liu et al. (1979) look-up table coefficients to compute roughness
!  Reynolds number for temperature (rt) and moisture (rq) as a
!  function of wind Reynolds number (rr):
!
!       rt = Liu_a(:,1) * Rr   ** Liu_b(:,1)    temperature
!       rq = Liu_a(:,2) * Rr   ** Liu_b(:,2)    moisture
!
   real(wp),parameter, dimension(8,2) :: Liu_a = reshape ( &
                 (/ 0.177,  1.376,    1.026,      1.625,   &
                    4.661, 34.904, 1667.190, 588000.0,     &
                    0.292,  1.808,    1.393,      1.956,   &
                    4.994, 30.709, 1448.680, 298000.0 /),  &
                 (/ 8, 2 /) )

   real(wp),parameter, dimension(8,2) :: Liu_b = reshape ( &
                 (/  0.0,    0.929, -0.599, -1.018,        &
                    -1.475, -2.067, -2.907, -3.935,        &
                     0.0,    0.826, -0.528, -0.870,        &
                    -1.297, -1.845, -2.682, -3.616 /),     &
                 (/ 8, 2 /) )

   real(wp),parameter, dimension(9) :: Liu_Rr =            &
                 (/    0.0,  0.11,   0.825,   3.0,         &
                      10.0, 30.0,  100.0,   300.0,         &
                    1000.0 /)
!
!  Height (m) of surface air temperature measurement.
   real(wp), parameter       ::  zt= 2.0
!  Height (m) of surface air humidity measurement
   real(wp), parameter       ::  zq= 2.0
!  Height (m) of surface winds measurement
   real(wp), parameter       ::  zw= 10.0
   integer,  parameter       :: itermax = 20

   real(wp), parameter       :: wgust=0.0


!   real(wp)                    :: psi

   real(wp), parameter         :: cpa=1008.
   real(wp), parameter         :: cpw=3985.
   real(wp), parameter         :: emiss=0.97
   real(wp), parameter         :: bolz=5.67e-8
   real(wp), parameter         :: kelvin=273.16
   real(wp), parameter         :: const06=0.62198
   real(wp), parameter         :: rgas = 287.1    !
   real(wp), parameter         :: g = 9.81        ! [m/s2]
   real(wp), parameter         :: rho_0 = 1025.   ! [kg/m3]
   real(wp), parameter         :: kappa = 0.41    ! von Karman
   logical                     :: rain_impact = .true.
   logical                     :: calc_evaporation =.true.

! !LOCAL VARIABLES:
   real(wp)                  :: tmp,cff,wgus
   real(wp)                  :: L
   real(wp)                  :: Cd
   real(wp)                  :: ta,ta_k,tw,tw_k
   integer                   :: ier,iter,k
   real(wp)                  :: vis_air
   real(wp)                  :: tpsi,qpsi,wpsi,ZWoL,oL,ZToL,ZQoL,ZoW,ZoT, ZoQ
   real(wp)                  :: Wstar,Tstar, Qstar, delQ, delT, rr,rt,rq
   real(wp)                  :: TVstar,Bf, upvel,delw,Wspeed, w
   real(wp)                  :: ri,cd_rain
   real(wp)                  :: x1,x2,x3
   real(wp)                  :: x
   real(wp)                  :: rainfall
   real(wp), parameter       :: eps=1.0e-12
!EOP
!-----------------------------------------------------------------------
!BOC
   evap = 0.0_wp
   w = sqrt(u10*u10+v10*v10)

   if (sst .lt. 100.) then
      tw  = sst
      tw_k= sst+kelvin
   else
      tw  = sst-kelvin
      tw_k= sst
   end if

   if (airt .lt. 100.) then
      ta_k  = airt + kelvin
      ta = airt
   else
      ta  = airt - kelvin
      ta_k = airt
   end if

!
!  Initialize.
!
   qe   = 0.0_wp
   qh   = 0.0_wp
   taux = 0.0_wp
   tauy = 0.0_wp
   delw=sqrt(w*w+wgust*wgust)
   if (delw .ne. 0.0) then
!-----------------------------------------------------------------------
!     Compute Monin-Obukhov similarity parameters for wind (Wstar),
!     heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------

!     Kinematic viscosity of dry air (m2/s), Andreas (1989).
      vis_air=1.326e-5*(1.0+ta*(6.542e-3+ta*(8.301e-6-4.84e-9*ta)))

!     Compute latent heat of vaporization (J/kg) at sea surface
      L = (2.501-0.00237*tw)*1.e6
!
!     Assume that wind is measured relative to sea surface and include
!     gustiness.
!     Initialize.
      ier = 0
      delq=qa-qs
      delt=ta-tw

!     Initial guesses for Monin-Obukhov similarity scales.
      ZWoL=0.0
      ZoW=0.0005
      Wstar=0.04*delw
      Tstar=0.04*delt
      Qstar=0.04*delq
      TVstar=Tstar*(1.0+0.61*qa)+0.61*ta_k*Qstar

!     Compute Richardson number.
      ri=g*zw*(delt+0.61*ta_k*delq)/(ta_k*delw*delw)

!     Fairall computes turbulent fluxes only if Ri< 0.25
      if ( ri .le. 0.25) then
!        Iterate until convergence or when IER is negative.  It usually
!        converges within four iterations.
         do iter=1,itermax
            if ( ier .ge. 0 ) then
!              Compute Monin-Obukhov stability parameter, Z/L.
               oL=g*kappa*TVstar/(ta_k*(1.0+0.61*qa)*Wstar*Wstar)
               ZWoL=zw*oL
               ZToL=zt*oL
               ZQoL=zq*oL

!              Evaluate stability functions at Z/L.
               wpsi=psi(1,ZWoL)
               tpsi=psi(2,ZToL)
               qpsi=psi(2,ZQoL)

!              Compute wind scaling parameters, Wstar.
               ZoW=0.011*Wstar*Wstar/g+0.11*vis_air/Wstar
               Wstar=delw*kappa/(log(zw/ZoW)-wpsi)

!              Computes roughness Reynolds number for wind (Rr), heat (Rt),
!              and moisture (Rq). Use Liu et al. (1976) look-up table to
!              compute "Rt" and "Rq" as function of "Rr".
               rr=ZoW*Wstar/vis_air
               if ((rr .ge. 0.0).and.(rr .lt. 1000.0)) then
                  do k=1,8
                     if ((liu_rr(k).le.rr).and.(rr .lt. liu_rr(k+1))) then
                        rt=liu_a(k,1)*rr**liu_b(k,1)
                        rq=liu_a(k,2)*rr**liu_b(k,2)
                     end if
                  end do

!                Compute heat and moisture scaling parameters,
!                Tstar and Qstar.
                  cff=vis_air/Wstar
                  ZoT=rt*cff
                  ZoQ=rq*cff
                  cff=kappa*fdg
                  Tstar=(delt)*cff/(log(zt/ZoT)-tpsi)
                  Qstar=(delq)*cff/(log(zq/ZoQ)-qpsi)

!                 Compute gustiness in wind speed.
                  TVstar=Tstar*(1.0+0.61*qa)+0.61*ta_k*Qstar
                  bf=-g/ta_k*Wstar*TVstar
                  if (bf .gt. 0) then
                     wgus=beta*(bf*Zabl)**r3
                  else
                     wgus=0.0_wp
                  end if
                  delw=sqrt(w*w+wgus*wgus)
               else
                  ier = -2
               end if
            end if
         end do

!        Compute transfer coefficients for momentun (Cd), heat (Ch),
!        and moisture (Ce).
         if (ier .ge. 0) then
            Wspeed=sqrt(w*w+wgus*wgus)
            Cd=Wstar*Wstar/(Wspeed*Wspeed)

!           Compute turbulent sensible heat flux (W/m2), qe.
!           out of ocean is negative
            qe=cpa*rhoa*Wstar*Tstar

!           compute sensible heatflux due to rain fall
            if (rain_impact) then
!              units of qs and qa - should be kg/kg
               rainfall=precip * 1000. ! (convert from m/s to kg/m2/s)
               x1 = 2.11e-5*(ta_k/kelvin)**1.94
               x2 = 0.02411*(1.0+ta*(3.309e-3-1.44e-6*ta))/(rhoa*cpa)
               x3 = qa * L /(rgas * ta_K * ta_K)
               cd_rain = 1.0/(1.0+const06*(x3*L*x1)/(cpa*x2))
               cd_rain = cd_rain*cpw*((tw-ta) + (qs-qa)*L/cpa)
               qe = qe - rainfall * cd_rain
            end if

!           Compute turbulent latent heat flux (W/m2), qh.
            qh=L*rhoa*Wstar*Qstar

!           Compute Webb correction (Webb effect) to latent heat flux
            upvel=-1.61*Wstar*Qstar-(1.0+1.61*qa)*Wstar*Tstar/ta_k
            qh=qh-rhoa*L*upvel*qa

!           calculation of evaporation/condensation in m/s
            if (rain_impact .and. calc_evaporation) then
               evap = rhoa/rho_0*Wstar*Qstar
            else
               evap = 0.0_wp
            end if

!           Compute wind stress components (N/m2), Tau.
            cff=rhoa*Cd*Wspeed
            taux=(cff*u10)
            tauy=(cff*v10)

!           Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!           according to Caldwell and Elliott (1971, JPO)
            if ( rain_impact ) then
               tmp  = 0.85d0 * rainfall
               taux  = taux + tmp * u10
               tauy  = tauy + tmp * v10
            end if

         end if ! ier >0
      end if ! Ri < 0.25
   end if  !delw != 0.0

   return
   end subroutine fairall
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the humidity \label{sec:humidity}
!
! !INTERFACE:
   subroutine humidity(hum_method,hum,airp,tw,ta,qa,qs,rhoa,ea,es)
!
! !DESCRIPTION:
!
! This routine calculated the saturation vapour pressure at SST and at
! air temperature, as well as the saturation specific humidty and the
! specific humidity. For the latter, four methods are implemented,
! and the method has to be chosen in the namelist file {\tt airsea.nml}
! as parameter {\tt hum\_method}, see \sect{sec:init-air-sea} for details.
!
!
! !USES:
!   use airsea_variables, only: kelvin,const06,rgas
!   use airsea_variables, only: es,ea,qs,qa,rhoa
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: hum_method
   real(wp), intent(in)                :: hum,airp,tw,ta

   real(wp), intent(out)               :: qa ! specific humidity (kg/kg)
   real(wp), intent(out)               :: ea ! actual water vapor pressure in Pascal
   real(wp), intent(out)               :: es ! saturation vapor pressure
   real(wp), intent(out)               :: qs ! saturation specific humidity
   real(wp), intent(out)               :: rhoa ! air density
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips, Hans Burchard & Karsten Bolding
!
! !DEFINED PARAMETERS:
   real(wp), parameter         :: kelvin=273.16
   real(wp), parameter         :: const06=0.62198
   real(wp), parameter         :: rgas = 287.1
!  Note shift of indices for coefficients compared to Lowe (1977, J. Appl. Met.)
   real(wp), parameter       :: a1=6.107799961
   real(wp), parameter       :: a2=4.436518521e-1
   real(wp), parameter       :: a3=1.428945805e-2
   real(wp), parameter       :: a4=2.650648471e-4
   real(wp), parameter       :: a5=3.031240396e-6
   real(wp), parameter       :: a6=2.034080948e-8
   real(wp), parameter       :: a7=6.136820929e-11
!
! !LOCAL VARIABLES:
   real(wp)        :: rh,twet,twet_k,dew,dew_k
!EOP
!-----------------------------------------------------------------------
!BOC
!  saturation vapor pressure - using SST
   es = a1 +tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*(a6+tw*a7)))))
   es = es * 100.0 ! Conversion millibar --> Pascal

!  correction for seawater, following Kraus 1972
!  correcting for salt water assuming 98% RH
   es=0.98 * es
!  saturation specific humidity
   qs = const06*es/(airp-0.377*es)

!  must be also calcuated for airtemperature, depending on humidity input
!  see ../ncdf/ncdf_meteo.F90 for defined constants
   select case (hum_method)
      case (1) ! relative humidity in % given
         rh = 0.01 * hum
!        saturation vapor pressure at that air temperature
         ea = a1 +ta*(a2+ta*(a3+ta*(a4+ta*(a5+ta*(a6+ta*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
!        get actual vapor pressure
         ea = rh * ea
!        convert to specific humidity
         qa = const06*ea/(airp-0.377*ea)
      case (2)  ! Specific humidity from wet bulb temperature
!        calculate the SVP at wet bulb temp then
!        use the psychrometer formula to get vapour pressure
!        See Smithsonian Met tables 6th Edition pg 366 eqn 3
!        Make sure this is in degC
         if (hum .lt. 100 ) then
            twet_k=hum + kelvin
            twet=hum
         else
            twet=hum - kelvin
            twet_k=hum
         end if
!        saturation vapor pressure at wet bulb temperature
         ea = a1 +twet*(a2+twet*(a3+twet*(a4+twet*(a5+twet*(a6+twet*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
!        actual vapor pressure
         ea = ea - 6.6e-4*(1+1.15e-3*twet)*airp*(ta-twet)
!        specific humidity in kg/kg
         qa = const06*ea/(airp-0.377*ea)
      case (3)  ! Specific humidity from dew point temperature
         if (hum .lt. 100.) then
            dew = hum
            dew_k = hum + kelvin
         else
            dew = hum - kelvin
            dew_k = hum
         end if
         ea = a1 +dew*(a2+dew*(a3+dew*(a4+dew*(a5+dew*(a6+dew*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
         qa = const06*ea/(airp-0.377*ea)
      case (4)
!        specific humidity in kg/kg is given
         qa = hum
!        actual water vapor pressure in Pascal
         ea = qa *airp/(const06+0.378*qa)
      case default
         stop 'not a valid hum_method bulk_fluxes()'
   end select

   rhoa = airp/(rgas*(ta+kelvin)*(1.0+const06*qa))

   return
   end subroutine humidity
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------


   function psi(iflag, ZoL)
!=======================================================================
!                                                                      !
!  This function evaluates the stability function, PSI, for wind       !
!  speed (iflag=1) or for air temperature and moisture (iflag=2)       !
!  profiles as function of the stability parameter, ZoL (z/L).         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!=======================================================================
!
      real(wp)  :: psi
!
!  Imported variable declarations.
!
   integer,  intent(in) :: iflag
   real(wp), intent(in) :: ZoL
!
!  Local variable declarations.
!
   real(wp), parameter :: r3 = 1.0/3.0
   real(wp), parameter :: sqr3 = 1.7320508
   real(wp), parameter :: pi=3.141592653589
   real(wp)            :: Fw, chic, chik, psic, psik

!  Initialize for the zero "ZoL" case.
!
   psi=0.0
!
!  Unstable conditions.
!
   if (ZoL .lt. 0.0) then
      chik=(1.0-16.0*ZoL)**0.25
      if (iflag .eq. 1) then
         psik=2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik*chik))-   &
              2.0*ATAN(chik)+ 0.5*pi
      else if (iflag .eq. 2) then
            psik=2.0*LOG(0.5*(1.0+chik*chik))
      end if
!
!  For very unstable conditions, use free-convection (Fairall).
!
      chic=(1.0-12.87*ZoL)**r3
      psic=1.5*LOG(r3*(1.0+chic+chic*chic))-                    &
         sqr3*ATAN((1.0+2.0*chic)/sqr3)+ pi/sqr3
!
!  Match Kansas and free-convection forms with weighting Fw.
!
      Fw=1.0/(1.0+ZoL*ZoL)
      psi=Fw*psik+(1.0-Fw)*psic
!
!  Stable conditions.
!
   else if (ZoL .gt. 0.0) then
      psi=-4.7*ZoL
   end if

   return
   end function psi

!EOC
!-----------------------------------------------------------------------
!Copyright (C) 2007 - Adolf Stips
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the long-wave back radiation \label{sec:back-rad}
!
! !INTERFACE:
   subroutine back_radiation(method,dlat,tw,ta,cloud,ea,qa,qb)
!
! !DESCRIPTION:
!
! Here, the long-wave back radiation is calculated by means of one out
! of four methods, which depend on the value given to the parameter
! {\tt method}:
! {\tt method}=1: \cite{Clarketal74},
! {\tt method}=2: \cite{HastenrathLamb78},
! {\tt method}=3: \cite{Bignamietal95},
! {\tt method}=4: \cite{BerliandBerliand52}.
! It should ne noted that the latitude must here be given in degrees.
!
! !USES:
!   use airsea_variables, only: emiss,bolz
!   use airsea_variables, only: ea,qa
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method
   real(wp), intent(in)                :: dlat,tw,ta,cloud
   real(wp), intent(in)                :: qa ! specific humidity (kg/kg)
   real(wp), intent(in)                :: ea ! actual water vapor pressure in Pascal
!
! !OUTPUT PARAMETERS:
   real(wp), intent(out)               :: qb
!
! !REVISION HISTORY:
!  Original author(s): Adols Stips, Hans Burchard & Karsten Bolding
!   real(wp), public, parameter         :: cpw=3985.
   real(wp), parameter         :: emiss=0.97
   real(wp), parameter         :: bolz=5.67e-8

! !LOCAL VARIABLES:

   integer, parameter   :: clark=1      ! Clark et al, 1974
   integer, parameter   :: hastenrath=2 ! Hastenrath and Lamb, 1978
   integer, parameter   :: bignami=3    ! Bignami et al., 1995 - Medsea
   integer, parameter   :: berliand=4   ! Berliand and Berliand, 1952 - ROMS


   real(wp), parameter, dimension(91)  :: cloud_correction_factor = (/ &
     0.497202,     0.501885,     0.506568,     0.511250,     0.515933, &
     0.520616,     0.525299,     0.529982,     0.534665,     0.539348, &
     0.544031,     0.548714,     0.553397,     0.558080,     0.562763, &
     0.567446,     0.572129,     0.576812,     0.581495,     0.586178, &
     0.590861,     0.595544,     0.600227,     0.604910,     0.609593, &
     0.614276,     0.618959,     0.623641,     0.628324,     0.633007, &
     0.637690,     0.642373,     0.647056,     0.651739,     0.656422, &
     0.661105,     0.665788,     0.670471,     0.675154,     0.679837, &
     0.684520,     0.689203,     0.693886,     0.698569,     0.703252, &
     0.707935,     0.712618,     0.717301,     0.721984,     0.726667, &
     0.731350,     0.736032,     0.740715,     0.745398,     0.750081, &
     0.754764,     0.759447,     0.764130,     0.768813,     0.773496, &
     0.778179,     0.782862,     0.787545,     0.792228,     0.796911, &
     0.801594,     0.806277,     0.810960,     0.815643,     0.820326, &
     0.825009,     0.829692,     0.834375,     0.839058,     0.843741, &
     0.848423,     0.853106,     0.857789,     0.862472,     0.867155, &
     0.871838,     0.876521,     0.881204,     0.885887,     0.890570, &
     0.895253,     0.899936,     0.904619,     0.909302,     0.913985, &
     0.918668 /)

   real(wp)                  :: ccf
   real(wp)                  :: x1,x2,x3
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  calculate cloud correction factor,fortran counts from 1 !
   ccf= cloud_correction_factor(nint(abs(dlat))+1)

   select case(method)
      case(clark)
!        Clark et al. (1974) formula.
!        unit of ea is Pascal, must hPa
!        Black body defect term, clouds, water vapor correction
         x1=(1.0-ccf*cloud*cloud)*(tw**4)
         x2=(0.39-0.05*sqrt(ea*0.01))
!        temperature jump term
         x3=4.0*(tw**3)*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case(hastenrath) ! qa in g(water)/kg(wet air)
!        Hastenrath and Lamb (1978) formula.
         x1=(1.0-ccf*cloud*cloud)*(tw**4)
         x2=(0.39-0.056*sqrt(1000.0*qa))
         x3=4.0*(tw**3)*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case(bignami)
!        Bignami et al. (1995) formula (Med Sea).
!        unit of ea is Pascal, must hPa
         ccf=0.1762
         x1=(1.0+ccf*cloud*cloud)*ta**4
         x2=(0.653+0.00535*(ea*0.01))
         x3= emiss*(tw**4)
         qb=-bolz*(-x1*x2+x3)
      case(berliand)
!        Berliand & Berliand (1952) formula (ROMS).
         x1=(1.0-0.6823*cloud*cloud)*ta**4
         x2=(0.39-0.05*sqrt(0.01*ea))
         x3=4.0*ta**3*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case default
   end select

   return
   end subroutine back_radiation
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the solar zenith angle \label{sec:swr}
!
! !INTERFACE:
   function solar_zenith_angle(yday,hh,dlon,dlat)
!
! !DESCRIPTION:
!  This subroutine calculates the solar zenith angle as being used both
!  in albedo_water() and short_wave_radiation(). The result is in degrees.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: yday
   real(wp), intent(in)                :: hh
   real(wp), intent(in)                :: dlon,dlat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   real(wp), parameter       :: pi=3.14159265358979323846
   real(wp), parameter       :: deg2rad=pi/180.
   real(wp), parameter       :: rad2deg=180./pi

   real(wp)                  :: rlon,rlat
   real(wp)                  :: yrdays
   real(wp)                  :: th0,th02,th03,sundec
   real(wp)                  :: thsun,coszen

   real(wp)                  :: solar_zenith_angle
!EOP
!-----------------------------------------------------------------------
!BOC
!  from now on everything in radians
   rlon = deg2rad*dlon
   rlat = deg2rad*dlat

   yrdays=365.25

   th0 = 2.*pi*yday/yrdays
   th02 = 2.*th0
   th03 = 3.*th0
!  sun declination :
   sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)         &
           - 0.006758*cos(th02) + 0.000907*sin(th02)                 &
           - 0.002697*cos(th03) + 0.001480*sin(th03)
!  sun hour angle :
   thsun = (hh-12.)*15.*deg2rad + rlon

!  cosine of the solar zenith angle :
   coszen =sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)
   if (coszen .lt. 0.0_wp) coszen = 0.0_wp

   solar_zenith_angle = rad2deg*acos(coszen)

   return
   end function solar_zenith_angle
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the short--wave radiation \label{sec:swr}
!
! !INTERFACE:
   function short_wave_radiation(zenith_angle,yday,dlon,dlat,cloud)
!
! !DESCRIPTION:
!  This subroutine calculates the short--wave net radiation based on
!  solar zenith angle, year day, longitude, latitude, and fractional cloud cover.
!  No corrections for albedo - must be done by calls to albedo\_water() and
!  if ice is included albedo\_ice().
!  The basic formula for the short-wave radiation at the surface, $Q_s$,
!  has been taken from \cite{RosatiMiyacoda88}, who adapted the work
!  of \cite{Reed77} and \cite{SimpsonPaulson99}:
!
!  \begin{equation}
!  Q_s=Q_{tot} (1-0.62 C + 0.0019 \beta) (1-\alpha),
!  \end{equation}
!
!  with the total radiation reaching the surface under clear skies,
!  $Q_{tot}$, the fractional cloud cover, $C$, the solar noon altitude,
!  $\beta$, and the albedo, $\alpha$.
!  This piece of code has been taken the MOM-I (Modular Ocean Model)
!  version at the INGV (Istituto Nazionale di Geofisica e Vulcanologia,
!  see {\tt http://www.bo.ingv.it/}).
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(wp), intent(in)                :: zenith_angle
   integer, intent(in)                 :: yday
   real(wp), intent(in)                :: dlon,dlat
   real(wp), intent(in)                :: cloud
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   real(wp), parameter       :: pi=3.14159265358979323846
   real(wp), parameter       :: deg2rad=pi/180.
   real(wp), parameter       :: rad2deg=180./pi

   real(wp), parameter       :: solar=1350.
   real(wp), parameter       :: eclips=23.439*deg2rad
   real(wp), parameter       :: tau=0.7
   real(wp), parameter       :: aozone=0.09

   real(wp)                  :: coszen,sunbet
   real(wp)                  :: qatten,qzer,qdir,qdiff,qtot,qshort
   real(wp)                  :: rlon,rlat,eqnx
   real(wp)                  :: yrdays

   real(wp)                  :: short_wave_radiation
!EOP
!-----------------------------------------------------------------------
!BOC
   coszen = cos(deg2rad*zenith_angle)
   if (coszen .le. 0.0) then
      coszen = 0.0
      qatten = 0.0
   else
      qatten = tau**(1.0_wp/coszen)
   end if

   qzer  = coszen * solar
   qdir  = qzer * qatten
   qdiff = ((1.0_wp-aozone)*qzer - qdir) * 0.5
   qtot  =  qdir + qdiff

!  from now on everything in radians
   rlon = deg2rad*dlon
   rlat = deg2rad*dlat

   yrdays=365.
   eqnx = (yday-81.)/yrdays*2.*pi
!  sin of the solar noon altitude in radians :
   sunbet=sin(rlat)*sin(eclips*sin(eqnx))+cos(rlat)*cos(eclips*sin(eqnx))
!  solar noon altitude in degrees :
   sunbet = asin(sunbet)*rad2deg

!  radiation as from Reed(1977), Simpson and Paulson(1979)
!  calculates SHORT WAVE FLUX ( watt/m*m )
!  Rosati,Miyakoda 1988 ; eq. 3.8
!  clouds from COADS perpetual data set
!#if 1
   qshort  = qtot*(1-0.62*cloud + .0019*sunbet)
   if(qshort .gt. qtot ) then
      qshort  = qtot
   end if
!#else
!  original implementation
!   if(cloud .lt. 0.3) then
!      qshort  = qtot
!   else
!      qshort  = qtot*(1-0.62*cloud + 0.0019*sunbet)
!   endif
!#endif
   short_wave_radiation = qshort

   return
   end function short_wave_radiation
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

END MODULE fv_sbc
