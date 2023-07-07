module module_cfradial_output
  use module_configuration, only : RKIND
  implicit none
  private
  public :: cfradial_type

  integer, parameter :: MAX_RAYS = 85000 ! 40000 ! 85000
  integer, parameter :: MAX_NGATES = 2500
  integer, parameter, public :: MAX_SWEEPS = 230 ! 15 ! 230

  type config_output_type
      logical :: output = .TRUE.
      integer :: varid = -99999
  end type config_output_type

  ! runtime-configurable output fields:
  type(config_output_type) :: config_zhh
  type(config_output_type) :: config_zhh_attenuated
  type(config_output_type) :: config_stdv_zhh
  type(config_output_type) :: config_zhh_noise_added
  type(config_output_type) :: config_zhh_cloud
  type(config_output_type) :: config_zhh_rain
  type(config_output_type) :: config_zhh_ice
  type(config_output_type) :: config_zhh_snow
  type(config_output_type) :: config_zhh_graupel
  type(config_output_type) :: config_zvv
  type(config_output_type) :: config_zvv_attenuated
  type(config_output_type) :: config_stdv_zvv
  type(config_output_type) :: config_zvv_noise_added
  type(config_output_type) :: config_zvh
  type(config_output_type) :: config_zvh_attenuated
  type(config_output_type) :: config_rhohv
  type(config_output_type) :: config_zdr
  type(config_output_type) :: config_zdr_attenuated
  type(config_output_type) :: config_ldrh
  type(config_output_type) :: config_adp
  type(config_output_type) :: config_kdp
  type(config_output_type) :: config_phidp
  type(config_output_type) :: config_phidp_wdp    
  type(config_output_type) :: config_dphase
  type(config_output_type) :: config_ah
  type(config_output_type) :: config_av
  type(config_output_type) :: config_vel
  type(config_output_type) :: config_vel_noise_added
  type(config_output_type) :: config_stdv_vel
  type(config_output_type) :: config_velunf
  type(config_output_type) :: config_velunf_fvel
  type(config_output_type) :: config_velunf_noise_added
  type(config_output_type) :: config_velunf_fvel_noise_added ! Added by B. Klotz (10/11/2021)
  type(config_output_type) :: config_veltru
  type(config_output_type) :: config_veltru_fvel
  type(config_output_type) :: config_veltru_noise_added
  type(config_output_type) :: config_veltru_fvel_noise_added ! Added by B. Klotz (10/11/2021)
  type(config_output_type) :: config_veltru_folded
  type(config_output_type) :: config_veltru_folded_noise_added
  type(config_output_type) :: config_veltru_folded_fvel ! Added by B. Klotz (10/11/2021)
  type(config_output_type) :: config_veltru_folded_fvel_noise_added ! Added by B. Klotz (10/11/2021)
  type(config_output_type) :: config_dvh ! Added by B. Klotz (10/5/2021)
  type(config_output_type) :: config_dopp ! Added by B. Klotz (10/5/2021)
  type(config_output_type) :: config_fvel ! Added by B. Klotz (10/5/2021)
  type(config_output_type) :: config_swtot
  type(config_output_type) :: config_swh
  type(config_output_type) :: config_sws
  type(config_output_type) :: config_swt
  type(config_output_type) :: config_swv
  type(config_output_type) :: config_snr_zhh
  type(config_output_type) :: config_snr_zvv
  type(config_output_type) :: config_x
  type(config_output_type) :: config_y
  type(config_output_type) :: config_z
  type(config_output_type) :: config_WRF_TEMPERATURE
  type(config_output_type) :: config_WRF_HT
  type(config_output_type) :: config_WRF_RHO_D
  type(config_output_type) :: config_WRF_QCLOUD
  type(config_output_type) :: config_WRF_QRAIN
  type(config_output_type) :: config_WRF_QICE
  type(config_output_type) :: config_WRF_QSNOW
  type(config_output_type) :: config_WRF_QGRAUP
  type(config_output_type) :: config_WRF_QNRAIN
  type(config_output_type) :: config_WRF_QNICE
  type(config_output_type) :: config_WRF_QNSNOW
  type(config_output_type) :: config_WRF_QNGRAUPEL
  type(config_output_type) :: config_WRF_U
  type(config_output_type) :: config_WRF_V
  type(config_output_type) :: config_WRF_W

  type, public :: volume_field_type
      character(len=256) :: field_name
      character(len=256) :: field_description
      character(len=256) :: field_units
      character(len=8) :: source
      integer :: varid
      real(kind=RKIND), allocatable, dimension(:,:) :: data
      type(volume_field_type), pointer :: next => NULL()
  end type volume_field_type

  type(volume_field_type), target  :: root

  type, public :: volume_type
      integer :: number
      integer :: nrays
      integer :: ngates
      integer :: sweep
      real(kind=RKIND)    :: beamwidth_h
      real(kind=RKIND)    :: beamwidth_v
      real(kind=RKIND)    :: fold_limit_lower
      real(kind=RKIND)    :: fold_limit_upper
      integer,           dimension(MAX_SWEEPS) :: sweep_number
      character(len=32), dimension(MAX_SWEEPS) :: sweep_mode
      real(kind=RKIND),              dimension(MAX_SWEEPS) :: fixed_angle
      integer,           dimension(MAX_SWEEPS) :: sweep_start_ray_index
      integer,           dimension(MAX_SWEEPS) :: sweep_end_ray_index
      character(len=19) :: time_coverage_start
      character(len=19) :: time_coverage_end
      character(len=19) :: reference_time
      integer           :: start_time_integer_seconds
      real(kind=RKIND), dimension(MAX_RAYS)  :: time
      real(kind=RKIND), dimension(MAX_NGATES) :: range
      real(kind=RKIND), dimension(MAX_RAYS)  :: latitude
      real(kind=RKIND), dimension(MAX_RAYS)  :: longitude
      real(kind=RKIND), dimension(MAX_RAYS)  :: altitude
      real(kind=RKIND), dimension(MAX_RAYS)  :: azimuth
      real(kind=RKIND), dimension(MAX_RAYS)  :: elevation
      integer, dimension(MAX_RAYS)  :: georefs_applied
      real(kind=RKIND), dimension(MAX_RAYS)  :: heading
      real(kind=RKIND), dimension(MAX_RAYS)  :: roll
      real(kind=RKIND), dimension(MAX_RAYS)  :: pitch
      real(kind=RKIND), dimension(MAX_RAYS)  :: drift
      real(kind=RKIND), dimension(MAX_RAYS)  :: rotation
      real(kind=RKIND), dimension(MAX_RAYS)  :: tilt

      ! "instrument_parameters" meta_group
      real(kind=RKIND), dimension(MAX_RAYS)  :: nyquist_velocity
      real(kind=RKIND), dimension(MAX_RAYS)  :: unambiguous_range

      ! "platform_velocity" meta_group
      real(kind=RKIND), dimension(MAX_RAYS)  :: eastward_velocity
      real(kind=RKIND), dimension(MAX_RAYS)  :: northward_velocity
      real(kind=RKIND), dimension(MAX_RAYS)  :: vertical_velocity
      real(kind=RKIND), dimension(MAX_RAYS)  :: eastward_wind
      real(kind=RKIND), dimension(MAX_RAYS)  :: northward_wind
      real(kind=RKIND), dimension(MAX_RAYS)  :: vertical_wind
      real(kind=RKIND), dimension(MAX_RAYS)  :: heading_rate
      real(kind=RKIND), dimension(MAX_RAYS)  :: roll_rate
      real(kind=RKIND), dimension(MAX_RAYS)  :: pitch_rate

      ! "geometry_correction" meta_group
      real(kind=RKIND) :: azimuth_correction = 0.0
      real(kind=RKIND) :: elevation_correction = 0.0
      real(kind=RKIND) :: range_correction = 0.0
      real(kind=RKIND) :: longitude_correction = 0.0
      real(kind=RKIND) :: latitude_correction = 0.0
      real(kind=RKIND) :: pressure_altitude_correction = 0.0
      real(kind=RKIND) :: radar_altitude_correction = 0.0
      real(kind=RKIND) :: eastward_ground_speed_correction = 0.0
      real(kind=RKIND) :: northward_ground_speed_correction = 0.0
      real(kind=RKIND) :: vertical_velocity_correction = 0.0
      real(kind=RKIND) :: heading_correction = 0.0
      real(kind=RKIND) :: roll_correction = 0.0
      real(kind=RKIND) :: pitch_correction = 0.0
      real(kind=RKIND) :: drift_correction = 0.0
      real(kind=RKIND) :: rotation_correction = 0.0
      real(kind=RKIND) :: tilt_correction = 0.0

      ! other variables needed for each beam, but not necessarily output in the CFRadial files.
      real(kind=RKIND), dimension(MAX_RAYS)  :: aircraft_xgrid
      real(kind=RKIND), dimension(MAX_RAYS)  :: aircraft_ygrid
      real(kind=RKIND), dimension(MAX_RAYS)  :: aircraft_z

    contains
      procedure :: initialize_new_field => volume_initialize_new_field
      procedure :: point_to => volume_point_to
      procedure :: point_to_data => volume_point_to_data
      procedure :: new_array => volume_new_array
      procedure :: next_point => volume_next_point
      procedure :: final => volume_destruct
      procedure :: radial_velocity
      procedure :: radial_velocity_fvel
      procedure :: fold_radial_velocity
      procedure :: compute_snr
      procedure :: compute_ips
      procedure :: add_zhh_noise
      procedure :: compute_phidp
      procedure :: attenuate

  end type volume_type

  interface volume_type
      module procedure volume_initialize
  end interface volume_type

  type cfradial_type
      integer :: ncid = -999999
      integer :: volume_number_varid
      integer :: platform_type_varid
      integer :: instrument_type_varid
      integer :: primary_axis_varid
      integer :: time_coverage_start_varid
      integer :: time_coverage_end_varid
      integer :: time_varid
      integer :: range_varid
      integer :: latitude_varid
      integer :: longitude_varid
      integer :: altitude_varid
      integer :: altitude_agl_varid
      integer :: sweep_number_varid
      integer :: sweep_mode_varid
      integer :: fixed_angle_varid
      integer :: sweep_start_ray_index_varid
      integer :: sweep_end_ray_index_varid
      integer :: azimuth_varid
      integer :: elevation_varid
      integer :: georefs_applied_varid
      integer :: heading_varid
      integer :: roll_varid
      integer :: pitch_varid
      integer :: drift_varid
      integer :: rotation_varid
      integer :: tilt_varid

      ! "radar_parameters" meta_group
      integer :: radar_beam_width_h_varid
      integer :: radar_beam_width_v_varid

      ! "instrument_parameters" meta_group
      integer :: nyquist_velocity_varid
      integer :: unambiguous_range_varid

      ! "geometry_correction" meta_group
      integer :: azimuth_correction_varid
      integer :: elevation_correction_varid
      integer :: range_correction_varid
      integer :: longitude_correction_varid
      integer :: latitude_correction_varid
      integer :: pressure_altitude_correction_varid
      integer :: radar_altitude_correction_varid
      integer :: eastward_ground_speed_correction_varid
      integer :: northward_ground_speed_correction_varid
      integer :: vertical_velocity_correction_varid
      integer :: heading_correction_varid
      integer :: roll_correction_varid
      integer :: pitch_correction_varid
      integer :: drift_correction_varid
      integer :: rotation_correction_varid
      integer :: tilt_correction_varid

      ! "platform_velocity" meta_group
      integer :: eastward_velocity_varid
      integer :: northward_velocity_varid
      integer :: vertical_velocity_varid
      integer :: eastward_wind_varid
      integer :: northward_wind_varid
      integer :: vertical_wind_varid
      integer :: heading_rate_varid
      integer :: roll_rate_varid
      integer :: pitch_rate_varid

    contains
      procedure :: open => cfradial_open
      procedure :: close => cfradial_close
      procedure :: prepare_metadata => cfradial_prepare_metadata
      procedure :: add_to_output_metadata => cfradial_add_to_output_metadata
      procedure :: write_volume
      procedure :: final => cfradial_destruct
  end type cfradial_type

contains

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine cfradial_open ( self, &
       &                     flnm, namelist_file, &
       &                     meters_between_gates, meters_to_center_of_first_gate, &
       &                     nrays, max_gates, sweep, time_coverage_start, &
       &                     fold_limit_lower, fold_limit_upper )
    use netcdf, only : NF90_CLOBBER
    use netcdf, only : NF90_64BIT_OFFSET
    use netcdf, only : NF90_GLOBAL
    use netcdf, only : NF90_HDF5
    use netcdf, only : NF90_INT
    use netcdf, only : NF90_CHAR
    use netcdf, only : NF90_BYTE
    use netcdf, only : NF90_FILL_BYTE
    use netcdf, only : NF90_SHORT
    use netcdf, only : NF90_FLOAT
    use netcdf, only : NF90_DOUBLE
    use netcdf, only : nf90_create
    use netcdf, only : nf90_def_dim
    use netcdf, only : nf90_def_var
    use netcdf, only : nf90_put_att
    use netcdf, only : nf90_enddef
    use kwm_date_utilities, only : geth_newdate
    class (cfradial_type), intent(inout) :: self
    character(len=*),      intent(in) :: flnm
    character(len=*),      intent(in) :: namelist_file
    real(kind=RKIND),                  intent(in) :: meters_between_gates
    real(kind=RKIND),                  intent(in) :: meters_to_center_of_first_gate
    integer,               intent(in) :: nrays
    integer,               intent(in) :: max_gates
    integer,               intent(in) :: sweep
    character(len=19),     intent(in) :: time_coverage_start
    real(kind=RKIND),                  intent(in) :: fold_limit_lower
    real(kind=RKIND),                  intent(in) :: fold_limit_upper
    integer :: ierr
    integer :: time_dimid
    integer :: range_dimid
    ! integer :: n_points_dimid
    integer :: sweep_dimid
    integer :: string_length_32_dimid
    integer :: string_length_64_dimid

    character(len=256) :: history_string
    character(len=8) :: history_date_string
    character(len=10) :: history_time_string
    character(len=5) :: history_zone_string
    integer, dimension(8) :: date_values_array
    integer, dimension(0) :: dum0

    integer :: namelist_unit

    namelist/config_output/ config_zhh,            &
         &                  config_zhh_attenuated, &
         &                  config_stdv_zhh,       &
         &                  config_zhh_noise_added,&
         &                  config_zhh_cloud,      &
         &                  config_zhh_rain,       &
         &                  config_zhh_ice,        &
         &                  config_zhh_snow,       &
         &                  config_zhh_graupel,    &
         &                  config_zvv,            &
         &                  config_zvv_attenuated, &
         &                  config_zvv_noise_added,&
         &                  config_zvh,            &
         &                  config_zvh_attenuated, &
         &                  config_rhohv,          &
         &                  config_zdr,            &
         &                  config_zdr_attenuated, &
         &                  config_ldrh,           &
         &                  config_Adp,            &
         &                  config_Kdp,            &
         &                  config_PHIdp,          &
         &                  config_PHIdp_wDP,      &
         &                  config_DPhase,         &
         &                  config_Ah,             &
         &                  config_Av,             &
         &                  config_Av,             &
         &                  config_vel,            &
         &                  config_vel_noise_added, &
         &                  config_stdv_vel,       &
         &                  config_velunf,         &
         &                  config_velunf_fvel,    &
         &                  config_velunf_noise_added,  &
         &                  config_velunf_fvel_noise_added,  &
         &                  config_veltru,         &
         &                  config_veltru_fvel,    &
         &                  config_veltru_noise_added,         &
         &                  config_veltru_fvel_noise_added,         &
         &                  config_veltru_folded,  &
         &                  config_veltru_folded_noise_added,  &
         &                  config_veltru_folded_fvel,  &
         &                  config_veltru_folded_fvel_noise_added,  &
         &                  config_dvh,            &
         &                  config_dopp,           &
         &                  config_fvel,           &
         &                  config_swtot,          &
         &                  config_swh,            &
         &                  config_sws,            &
         &                  config_swt,            &
         &                  config_swv,            &
         &                  config_snr_zhh,        &
         &                  config_snr_zvv,        &
         &                  config_x,              &
         &                  config_y,              &
         &                  config_z,              &
         &                  config_WRF_TEMPERATURE,&
         &                  config_WRF_HT,         &
         &                  config_WRF_RHO_D,      &
         &                  config_WRF_QCLOUD,     &
         &                  config_WRF_QRAIN,      &
         &                  config_WRF_QICE,       &
         &                  config_WRF_QSNOW,      &
         &                  config_WRF_QGRAUP,     &
         &                  config_WRF_QNRAIN,     &
         &                  config_WRF_QNICE,      &
         &                  config_WRF_QNSNOW,     &
         &                  config_WRF_QNGRAUPEL,  &
         &                  config_WRF_U,          &
         &                  config_WRF_V,          &
         &                  config_WRF_W

    open(newunit=namelist_unit, form='formatted', action='read', file=trim(namelist_file) )
    read(namelist_unit, config_output)
    close(namelist_unit)
    ! write(*,config_output)

    ! ierr = nf90_create ( trim(flnm), NF90_CLOBBER .or. NF90_64BIT_OFFSET, self%ncid )
    ! call error_handler(ierr, "Problem opening file " //trim(flnm) // " for output")

    ierr = nf90_create ( trim(flnm), IOR(NF90_CLOBBER,NF90_HDF5), self%ncid )
    call error_handler(ierr, "Problem opening file " //trim(flnm) // " for output")

    !
    ! "The time coordinate indicates the number of rays in the file."
    !

    ierr = nf90_def_dim ( self%ncid, "time", nrays, time_dimid)
    call error_handler(ierr, "Problem defining time dimension")

    !
    ! "The range coordinate indicates the maximum number of gates for any ray in the file."
    !

    ierr = nf90_def_dim ( self%ncid, "range", max_gates, range_dimid )
    call error_handler(ierr, "Problem defining dimension 'range'")

    !
    ! "The n_points dimension indicates the total number of gates
    ! stored in all of the rays. It is equal to the sum of ray_n_gates
    ! over all rays."
    !
    ! "n_points is required if the number of gates varies by ray. It
    ! must not be included if the number of gates is fixed."
    !

    ! For simplicity for starters, I'll assume number of gates is fixes, so we 
    ! dispense with n_points dimension
    ! ierr = nf90_def_dim ( self%ncid, "n_points", nrays*max_gates, n_points_dimid )
    ! call error_handler ( ierr , "Problem defining n_points dimension" )

    ierr = nf90_def_dim ( self%ncid, "sweep", sweep, sweep_dimid)
    call error_handler ( ierr , "Problem defining sweep dimension")

    ierr = nf90_def_dim ( self%ncid, "string_length_64", 64, string_length_64_dimid )
    call error_handler ( ierr , "Problem defining string_length_64 dimension")

    ierr = nf90_def_dim ( self%ncid, "string_length_32", 32, string_length_32_dimid )
    call error_handler ( ierr , "Problem defining string_length_32 dimension")

    !
    ! Global Attributes
    !

    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "Conventions", &
         "CF/Radial instrument_parameters radar_parameters platform_velocity geometry_correction" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "version", "1.4" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "title", "WRF-derived simulated APAR" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "institution", "NCAR/MMM and NCAR/EOL" )
    ! ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "references", "George Bryan CM1" )
    ! ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "source", "CM1_as_APAR" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "references", "WRF and CRSIM to CFRadial" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "source", "WRF_as_APAR" )

    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "scan_name", "" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "scan_id", 0 )

    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "platform_is_mobile", "true" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "n_gates_vary", "false" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "ray_times_increase", "true" )

    call date_and_time ( VALUES=date_values_array )

    write(history_string, '("Created ",I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2,"Z")') &
         date_values_array(1), date_values_array(2), date_values_array(3), &
         date_values_array(5), date_values_array(6), date_values_array(7)
    ! Correct for time zone.
    call geth_newdate ( history_string(9:27), history_string(9:27), date_values_array(4)*(-60) )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "history", trim(history_string) )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "comment", "APAR simulated by WRF" )
    ierr = nf90_put_att ( self%ncid, NF90_GLOBAL, "instrument_name", "APAR simulated by WRF" )

    !
    ! Define variables.
    !

    !
    ! Global variables.
    !

    ierr = nf90_def_var ( self%ncid, "volume_number", NF90_INT, varid=self%volume_number_varid)
    ierr = nf90_put_att ( self%ncid, self%volume_number_varid, "long_name", "data_volume_index_number")
    ierr = nf90_put_att ( self%ncid, self%volume_number_varid, "units", "")
    ierr = nf90_put_att ( self%ncid, self%volume_number_varid, "_FillValue", -9999)

    ierr = nf90_def_var ( self%ncid, "platform_type", NF90_CHAR, dimids=(/string_length_32_dimid/), varid=self%platform_type_varid)
    ierr = nf90_put_att ( self%ncid, self%platform_type_varid, "long_name", "platform_type")
    ierr = nf90_put_att ( self%ncid, self%platform_type_varid, "options", "fixed, vehicle, ship, aircraft, aircraft_fore,&
         & aircraft_aft, aircraft_tail, aircraft_belly, aircraft_roof, aircraft_nose, satellite_orbit, satellite_geostat" )

    ierr = nf90_def_var ( self%ncid, "instrument_type", NF90_CHAR, dimids=(/string_length_32_dimid/), varid=self%instrument_type_varid)
    ierr = nf90_put_att ( self%ncid, self%instrument_type_varid, "long_name", "type_of_instrument")
    ierr = nf90_put_att ( self%ncid, self%instrument_type_varid, "options", "radar, lidar" )

    ierr = nf90_def_var ( self%ncid, "primary_axis", NF90_CHAR, dimids=(/string_length_32_dimid/), varid=self%primary_axis_varid)
    ierr = nf90_put_att ( self%ncid, self%primary_axis_varid, "long_name", "primary_axis_of_rotation")
    ierr = nf90_put_att ( self%ncid, self%primary_axis_varid, "options", "axis_z, axis_y, axis_x, axis_z_prime, axis_y_prime, axis_x_prime" )


    ierr = nf90_def_var ( self%ncid, "time_coverage_start", NF90_CHAR, dimids=(/string_length_32_dimid/), varid=self%time_coverage_start_varid)
    ierr = nf90_put_att ( self%ncid, self%time_coverage_start_varid, "long_name", "data_volume_start_time_utc")
    ierr = nf90_put_att ( self%ncid, self%time_coverage_start_varid, "comment", "ray times are relative to start time in secs" )

    ierr = nf90_def_var ( self%ncid, "time_coverage_end", NF90_CHAR, dimids=(/string_length_32_dimid/), varid=self%time_coverage_end_varid)
    ierr = nf90_put_att ( self%ncid, self%time_coverage_end_varid, "long_name", "data_volume_end_time_utc")


    !
    ! Coordinate variables.
    !

    ierr = nf90_def_var ( self%ncid, "time", NF90_DOUBLE, dimids=(/time_dimid/), varid=self%time_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%time_varid, "standard_name", "time" )
    ierr = nf90_put_att ( self%ncid, self%time_varid, "long_name", "time_in_seconds_since_volume_start" )
    ierr = nf90_put_att ( self%ncid, self%time_varid, "units", "seconds since "//time_coverage_start//"Z" )


    ierr = nf90_def_var ( self%ncid, "range", NF90_FLOAT, dimids=(/range_dimid/), varid=self%range_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%range_varid, "standard_name", "projection_range_coordinate" )
    ierr = nf90_put_att ( self%ncid, self%range_varid, "long_name", "range_to_measurement_volume" )
    ierr = nf90_put_att ( self%ncid, self%range_varid, "units", "meters" )
    ierr = nf90_put_att ( self%ncid, self%range_varid, "spacing_is_constant", "true" )
    ierr = nf90_put_att ( self%ncid, self%range_varid, "meters_to_center_of_first_gate", real(meters_to_center_of_first_gate,4) )
    ierr = nf90_put_att ( self%ncid, self%range_varid, "meters_between_gates", real(meters_between_gates,4) )
    ierr = nf90_put_att ( self%ncid, self%range_varid, "axis", "radial_range_coordinate" )

    !
    ! Location variables.
    !

    ierr = nf90_def_var ( self%ncid, "latitude", NF90_DOUBLE, dimids=(/time_dimid/), varid=self%latitude_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%latitude_varid, "long_name", "latitude")
    ierr = nf90_put_att ( self%ncid, self%latitude_varid, "units", "degrees_north")
    ! ierr = nf90_put_att ( self%ncid, self%latitude_varid, "_FillValue", -9999.D0)
    ierr = nf90_put_att ( self%ncid, self%latitude_varid, "_FillValue", -9999.0_8)

    ierr = nf90_def_var ( self%ncid, "longitude", NF90_DOUBLE, dimids=(/time_dimid/), varid=self%longitude_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%longitude_varid, "long_name", "longitude")
    ierr = nf90_put_att ( self%ncid, self%longitude_varid, "units", "degrees_east")
    ! ierr = nf90_put_att ( self%ncid, self%longitude_varid, "_FillValue", -9999.D0)
    ierr = nf90_put_att ( self%ncid, self%longitude_varid, "_FillValue", -9999.0_8)

    ierr = nf90_def_var ( self%ncid, "altitude", NF90_DOUBLE, dimids=(/time_dimid/), varid=self%altitude_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%altitude_varid, "long_name", "altitude")
    ierr = nf90_put_att ( self%ncid, self%altitude_varid, "units", "meters")
    ierr = nf90_put_att ( self%ncid, self%altitude_varid, "positive", "up")
    ierr = nf90_put_att ( self%ncid, self%altitude_varid, "_FillValue", -9999.0_8)

    ierr = nf90_def_var ( self%ncid, "altitude_agl", NF90_DOUBLE, dimids=(/time_dimid/), varid=self%altitude_agl_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%altitude_agl_varid, "long_name", "altitude_above_ground_level")
    ierr = nf90_put_att ( self%ncid, self%altitude_agl_varid, "units", "meters")
    ierr = nf90_put_att ( self%ncid, self%altitude_agl_varid, "positive", "up")
    ierr = nf90_put_att ( self%ncid, self%altitude_agl_varid, "_FillValue", -9999.0_8)

    !
    ! Sweep variables.
    !

    ierr = nf90_def_var ( self%ncid, "sweep_number", NF90_INT, dimids=(/sweep_dimid/), varid=self%sweep_number_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%sweep_number_varid, "long_name", "sweep_index_number_0_based" )
    ierr = nf90_put_att ( self%ncid, self%sweep_number_varid, "units", "" )
    ierr = nf90_put_att ( self%ncid, self%sweep_number_varid, "_FillValue", -9999 )
    
    ierr = nf90_def_var ( self%ncid, "sweep_mode", NF90_CHAR, dimids=(/string_length_32_dimid,sweep_dimid/), varid=self%sweep_mode_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%sweep_mode_varid, "long_name", "scan_mode_for_sweep")
    ierr = nf90_put_att ( self%ncid, self%sweep_mode_varid, "options", "sector, coplane, rhi, vertical_pointing, idle,&
         & azimuth_surveillance, elevation_surveillance, sunscan, pointing,&
         & calibration, manual_ppi, manual_rhi, sunscan_rhi" )


    ierr = nf90_def_var ( self%ncid, "fixed_angle", NF90_FLOAT, dimids=(/sweep_dimid/), varid=self%fixed_angle_varid , &
         &                deflate_level=2, shuffle=.true. )
    ! ierr = nf90_put_att ( self%ncid, self%fixed_angle_varid, "long_name", "ray_target_fixed_angle" )
    ierr = nf90_put_att ( self%ncid, self%fixed_angle_varid, "long_name", "target_fixed_angle" )
    ierr = nf90_put_att ( self%ncid, self%fixed_angle_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%fixed_angle_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "sweep_start_ray_index", NF90_INT, dimids=(/sweep_dimid/), varid=self%sweep_start_ray_index_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%sweep_start_ray_index_varid, "long_name", "index_of_first_ray_in_sweep" )
    ierr = nf90_put_att ( self%ncid, self%sweep_start_ray_index_varid, "units", "" )
    ierr = nf90_put_att ( self%ncid, self%sweep_start_ray_index_varid, "_FillValue", -9999 )

    ierr = nf90_def_var ( self%ncid, "sweep_end_ray_index", NF90_INT, dimids=(/sweep_dimid/), varid=self%sweep_end_ray_index_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%sweep_end_ray_index_varid, "long_name", "index_of_last_ray_in_sweep" )
    ierr = nf90_put_att ( self%ncid, self%sweep_end_ray_index_varid, "units", "" )
    ierr = nf90_put_att ( self%ncid, self%sweep_end_ray_index_varid, "_FillValue", -9999 )

    !
    ! Sensor pointing variables.
    !

    ierr = nf90_def_var ( self%ncid, "azimuth", NF90_FLOAT, dimids=(/time_dimid/), varid=self%azimuth_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%azimuth_varid, "standard_name", "ray_azimuth_angle" )
    ierr = nf90_put_att ( self%ncid, self%azimuth_varid, "long_name", "azimuth_angle_from_true_north" )
    ierr = nf90_put_att ( self%ncid, self%azimuth_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%azimuth_varid, "axis", "radial_azimuth_coordinate" )
    ierr = nf90_put_att ( self%ncid, self%azimuth_varid, "_FillValue", -9999._4 )

    ierr = nf90_def_var ( self%ncid, "elevation", NF90_FLOAT, dimids=(/time_dimid/), varid=self%elevation_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%elevation_varid, "standard_name", "ray_elevation_angle" )
    ierr = nf90_put_att ( self%ncid, self%elevation_varid, "long_name", "elevation_angle_from_horizontal_plane" )
    ierr = nf90_put_att ( self%ncid, self%elevation_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%elevation_varid, "axis", "radial_elevation_coordinate" )
    ierr = nf90_put_att ( self%ncid, self%elevation_varid, "positive", "up" )
    ierr = nf90_put_att ( self%ncid, self%elevation_varid, "_FillValue", -9999._4 )

    ierr = nf90_def_var ( self%ncid, "georefs_applied", NF90_BYTE, dimids=(/time_dimid/), varid=self%georefs_applied_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%georefs_applied_varid, "standard_name", "georefs_applied" )
    ierr = nf90_put_att ( self%ncid, self%georefs_applied_varid, "long_name", "georefs_applied" )
    ierr = nf90_put_att ( self%ncid, self%georefs_applied_varid, "units", "" )
    ierr = nf90_put_att ( self%ncid, self%georefs_applied_varid, "_FillValue", NF90_FILL_BYTE )

    !
    ! Moving platform georeference variables.
    !

    ierr = nf90_def_var ( self%ncid, "heading", NF90_FLOAT, dimids=(/time_dimid/), varid=self%heading_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%heading_varid, "long_name", "platform_heading_angle" )
    ierr = nf90_put_att ( self%ncid, self%heading_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%heading_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "roll", NF90_FLOAT, dimids=(/time_dimid/), varid=self%roll_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%roll_varid, "long_name", "platform_roll_angle" )
    ierr = nf90_put_att ( self%ncid, self%roll_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%roll_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "pitch", NF90_FLOAT, dimids=(/time_dimid/), varid=self%pitch_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%pitch_varid, "long_name", "platform_pitch_angle" )
    ierr = nf90_put_att ( self%ncid, self%pitch_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%pitch_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "drift", NF90_FLOAT, dimids=(/time_dimid/), varid=self%drift_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%drift_varid, "long_name", "platform_drift_angle" )
    ierr = nf90_put_att ( self%ncid, self%drift_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%drift_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "rotation", NF90_FLOAT, dimids=(/time_dimid/), varid=self%rotation_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%rotation_varid, "long_name", "ray_rotation_angle_relative_to_platform" )
    ierr = nf90_put_att ( self%ncid, self%rotation_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%rotation_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "tilt", NF90_FLOAT, dimids=(/time_dimid/), varid=self%tilt_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%tilt_varid, "long_name", "ray_tilt_angle_relative_to_platform" )
    ierr = nf90_put_att ( self%ncid, self%tilt_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%tilt_varid, "_FillValue", -9999.0_4 )

    !
    !  Instrument_parameters meta-group
    !

    ierr = nf90_def_var ( self%ncid, "nyquist_velocity", NF90_FLOAT, dimids=(/time_dimid/), varid=self%nyquist_velocity_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%nyquist_velocity_varid, "meta_group", "instrument_parameters" )
    ierr = nf90_put_att ( self%ncid, self%nyquist_velocity_varid, "long_name", "unambiguous_doppler_velocity" )
    ierr = nf90_put_att ( self%ncid, self%nyquist_velocity_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%nyquist_velocity_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "unambiguous_range", NF90_FLOAT, dimids=(/time_dimid/), varid=self%unambiguous_range_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%unambiguous_range_varid, "meta_group", "instrument_parameters" )
    ierr = nf90_put_att ( self%ncid, self%unambiguous_range_varid, "long_name", "nyquist_velocity" )
    ierr = nf90_put_att ( self%ncid, self%unambiguous_range_varid, "units", "meters" )
    ierr = nf90_put_att ( self%ncid, self%unambiguous_range_varid, "_FillValue", -9999.0_4 )

    !
    !  geometry_correction meta_group
    !

    ierr = nf90_def_var ( self%ncid, "azimuth_correction", NF90_FLOAT, dimids=dum0, varid=self%azimuth_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%azimuth_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%azimuth_correction_varid, "long_name", "azimuth_angle_correction" )
    ierr = nf90_put_att ( self%ncid, self%azimuth_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%azimuth_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "elevation_correction", NF90_FLOAT, dimids=dum0, varid=self%elevation_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%elevation_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%elevation_correction_varid, "long_name", "ray_elevation_angle_correction" )
    ierr = nf90_put_att ( self%ncid, self%elevation_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%elevation_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "range_correction", NF90_FLOAT, dimids=dum0, varid=self%range_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%range_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%range_correction_varid, "long_name", "range_to_center_of_measurement_volume_correction" )
    ierr = nf90_put_att ( self%ncid, self%range_correction_varid, "units", "meters" )
    ierr = nf90_put_att ( self%ncid, self%range_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "longitude_correction", NF90_FLOAT, dimids=dum0, varid=self%longitude_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%longitude_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%longitude_correction_varid, "long_name", "longitude_correction" )
    ierr = nf90_put_att ( self%ncid, self%longitude_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%longitude_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "latitude_correction", NF90_FLOAT, dimids=dum0, varid=self%latitude_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%latitude_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%latitude_correction_varid, "long_name", "latitude_correction" )
    ierr = nf90_put_att ( self%ncid, self%latitude_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%latitude_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "pressure_altitude_correction", NF90_FLOAT, dimids=dum0, varid=self%pressure_altitude_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%pressure_altitude_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%pressure_altitude_correction_varid, "long_name", "pressure_altitude_correction" )
    ierr = nf90_put_att ( self%ncid, self%pressure_altitude_correction_varid, "units", "meters" )
    ierr = nf90_put_att ( self%ncid, self%pressure_altitude_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "radar_altitude_correction", NF90_FLOAT, dimids=dum0, varid=self%radar_altitude_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%radar_altitude_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%radar_altitude_correction_varid, "long_name", "radar_altitude_correction" )
    ierr = nf90_put_att ( self%ncid, self%radar_altitude_correction_varid, "units", "meters" )
    ierr = nf90_put_att ( self%ncid, self%radar_altitude_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "eastward_ground_speed_correction", NF90_FLOAT, dimids=dum0, varid=self%eastward_ground_speed_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%eastward_ground_speed_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%eastward_ground_speed_correction_varid, "long_name", "eastward_ground_speed_correction" )
    ierr = nf90_put_att ( self%ncid, self%eastward_ground_speed_correction_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%eastward_ground_speed_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "northward_ground_speed_correction", NF90_FLOAT, dimids=dum0, varid=self%northward_ground_speed_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%northward_ground_speed_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%northward_ground_speed_correction_varid, "long_name", "northward_ground_speed_correction" )
    ierr = nf90_put_att ( self%ncid, self%northward_ground_speed_correction_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%northward_ground_speed_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "vertical_velocity_correction", NF90_FLOAT, dimids=dum0, varid=self%vertical_velocity_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_correction_varid, "long_name", "platform_vertical_velocity_correction" )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_correction_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "heading_correction", NF90_FLOAT, dimids=dum0, varid=self%heading_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%heading_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%heading_correction_varid, "long_name", "platform_heading_angle_correction" )
    ierr = nf90_put_att ( self%ncid, self%heading_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%heading_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "roll_correction", NF90_FLOAT, dimids=dum0, varid=self%roll_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%roll_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%roll_correction_varid, "long_name", "platform_roll_angle_correction" )
    ierr = nf90_put_att ( self%ncid, self%roll_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%roll_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "pitch_correction", NF90_FLOAT, dimids=dum0, varid=self%pitch_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%pitch_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%pitch_correction_varid, "long_name", "platform_pitch_angle_correction" )
    ierr = nf90_put_att ( self%ncid, self%pitch_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%pitch_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "drift_correction", NF90_FLOAT, dimids=dum0, varid=self%drift_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%drift_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%drift_correction_varid, "long_name", "platform_drift_angle_correction" )
    ierr = nf90_put_att ( self%ncid, self%drift_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%drift_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "rotation_correction", NF90_FLOAT, dimids=dum0, varid=self%rotation_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%rotation_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%rotation_correction_varid, "long_name", "ray_rotation_angle_relative_to_platform_correction" )
    ierr = nf90_put_att ( self%ncid, self%rotation_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%rotation_correction_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "tilt_correction", NF90_FLOAT, dimids=dum0, varid=self%tilt_correction_varid )
    ierr = nf90_put_att ( self%ncid, self%tilt_correction_varid, "meta_group", "geometry_correction" )
    ierr = nf90_put_att ( self%ncid, self%tilt_correction_varid, "long_name", "ray_tilt_angle_relative_to_platform_correction" )
    ierr = nf90_put_att ( self%ncid, self%tilt_correction_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%tilt_correction_varid, "_FillValue", -9999.0_4 )

    !
    ! Moving platform platform-velocity variables.
    !

    ierr = nf90_def_var ( self%ncid, "eastward_velocity", NF90_FLOAT, dimids=(/time_dimid/), varid=self%eastward_velocity_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%eastward_velocity_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%eastward_velocity_varid, "long_name", "platform_eastward_velocity" )
    ierr = nf90_put_att ( self%ncid, self%eastward_velocity_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%eastward_velocity_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "northward_velocity", NF90_FLOAT, dimids=(/time_dimid/), varid=self%northward_velocity_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%northward_velocity_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%northward_velocity_varid, "long_name", "platform_northward_velocity" )
    ierr = nf90_put_att ( self%ncid, self%northward_velocity_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%northward_velocity_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "vertical_velocity", NF90_FLOAT, dimids=(/time_dimid/), varid=self%vertical_velocity_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_varid, "long_name", "platform_vertical_velocity" )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%vertical_velocity_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "eastward_wind", NF90_FLOAT, dimids=(/time_dimid/), varid=self%eastward_wind_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%eastward_wind_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%eastward_wind_varid, "long_name", "eastward_wind" )
    ierr = nf90_put_att ( self%ncid, self%eastward_wind_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%eastward_wind_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "northward_wind", NF90_FLOAT, dimids=(/time_dimid/), varid=self%northward_wind_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%northward_wind_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%northward_wind_varid, "long_name", "northward_wind" )
    ierr = nf90_put_att ( self%ncid, self%northward_wind_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%northward_wind_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "vertical_wind", NF90_FLOAT, dimids=(/time_dimid/), varid=self%vertical_wind_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%vertical_wind_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%vertical_wind_varid, "long_name", "upward_air_velocity" )
    ierr = nf90_put_att ( self%ncid, self%vertical_wind_varid, "units", "meters per second" )
    ierr = nf90_put_att ( self%ncid, self%vertical_wind_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "heading_rate", NF90_FLOAT, dimids=(/time_dimid/), varid=self%heading_rate_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%heading_rate_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%heading_rate_varid, "long_name", "platform_heading_angle_rate_of_change" )
    ierr = nf90_put_att ( self%ncid, self%heading_rate_varid, "units", "degrees per second" )
    ierr = nf90_put_att ( self%ncid, self%heading_rate_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "roll_rate", NF90_FLOAT, dimids=(/time_dimid/), varid=self%roll_rate_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%roll_rate_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%roll_rate_varid, "long_name", "platform_roll_angle_rate_of_change" )
    ierr = nf90_put_att ( self%ncid, self%roll_rate_varid, "units", "degrees per second" )
    ierr = nf90_put_att ( self%ncid, self%roll_rate_varid, "_FillValue", -9999.0_4 )

    ierr = nf90_def_var ( self%ncid, "pitch_rate", NF90_FLOAT, dimids=(/time_dimid/), varid=self%pitch_rate_varid, &
         &                deflate_level=2, shuffle=.true. )
    ierr = nf90_put_att ( self%ncid, self%pitch_rate_varid, "meta_group", "platform_velocity" )
    ierr = nf90_put_att ( self%ncid, self%pitch_rate_varid, "long_name", "platform_pitch_angle_rate_of_change" )
    ierr = nf90_put_att ( self%ncid, self%pitch_rate_varid, "units", "degrees per second" )
    ierr = nf90_put_att ( self%ncid, self%pitch_rate_varid, "_FillValue", -9999.0_4 )

    !
    !  Radar parameters
    !

    ierr = nf90_def_var ( self%ncid, "radar_beam_width_h", NF90_FLOAT, varid=self%radar_beam_width_h_varid )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_h_varid, "long_name", "half_power_radar_beam_width_h_channel" )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_h_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_h_varid, "_FillValue", -9999.0_4 )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_h_varid, "meta_group" , "radar_parameters" )

    ierr = nf90_def_var ( self%ncid, "radar_beam_width_v", NF90_FLOAT, varid=self%radar_beam_width_v_varid )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_v_varid, "long_name", "half_power_radar_beam_width_v_channel" )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_v_varid, "units", "degrees" )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_v_varid, "FillValue", -9999.0_4 )
    ierr = nf90_put_att ( self%ncid, self%radar_beam_width_v_varid, "meta_group", "radar_parameters" )

    !
    ! Moments variables
    !

    ! Unfolded radial velocity
    if (config_velunf%output) then
        ierr = nf90_def_var ( self%ncid, "VELUNF", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_velunf%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_velunf%varid, "standard_name", "unfolded_radial_velocity_of_scatterers_away_from_instrument")
        ierr = nf90_put_att ( self%ncid, config_velunf%varid, "long_name", "unfolded_radial_velocity")
        ierr = nf90_put_att ( self%ncid, config_velunf%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_velunf%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_velunf%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_velunf%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    ! Unfolded radial velocity with fall velocity included
    !!! Added by B. Klotz (10/11/2021)
    if (config_velunf_fvel%output) then
       ierr = nf90_def_var ( self%ncid, "VELUNF_FVEL", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_velunf_fvel%varid, &
         &                deflate_level=2, shuffle=.true. )
       ierr = nf90_put_att ( self%ncid, config_velunf_fvel%varid, "standard_name", "unfolded_radial_velocity_of_scatterers_away_from_instrument_with_fall_velocity")
       ierr = nf90_put_att ( self%ncid, config_velunf_fvel%varid, "long_name", "unfolded_radial_velocity_with_fvel")
       ierr = nf90_put_att ( self%ncid, config_velunf_fvel%varid, "units", "meters per second")
       ierr = nf90_put_att ( self%ncid, config_velunf_fvel%varid, "field_folds", "false")
       ierr = nf90_put_att ( self%ncid, config_velunf_fvel%varid, "_FillValue", -9999.0_4 )
       ierr = nf90_put_att ( self%ncid, config_velunf_fvel%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_velunf_noise_added%output) then
        ierr = nf90_def_var ( self%ncid, "VELUNF_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_velunf_noise_added%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_velunf_noise_added%varid, "standard_name", "unfolded_radial_velocity_noise_added")
        ierr = nf90_put_att ( self%ncid, config_velunf_noise_added%varid, "long_name", "unfolded_radial_velocity_noise_added")
        ierr = nf90_put_att ( self%ncid, config_velunf_noise_added%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_velunf_noise_added%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_velunf_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_velunf_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    !!! Added by B. Klotz (10/11/2021)
    if (config_velunf_fvel_noise_added%output) then
        ierr = nf90_def_var ( self%ncid, "VELUNF_FVEL_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_velunf_fvel_noise_added%varid, &
             &     deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_velunf_fvel_noise_added%varid, "standard_name", "unfolded_radial_velocity_with_fall_velocity_and_noise_added")
        ierr = nf90_put_att ( self%ncid, config_velunf_fvel_noise_added%varid, "long_name", "unfolded_radial_velocity_with_fvel_noise_added")
        ierr = nf90_put_att ( self%ncid, config_velunf_fvel_noise_added%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_velunf_fvel_noise_added%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_velunf_fvel_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_velunf_fvel_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    ! Folded radial velocity
    if (config_vel%output) then
        ierr = nf90_def_var ( self%ncid, "VEL", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_vel%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_vel%varid, "standard_name", "radial_velocity_of_scatterers_away_from_instrument")
        ierr = nf90_put_att ( self%ncid, config_vel%varid, "long_name", "radial_velocity")
        ierr = nf90_put_att ( self%ncid, config_vel%varid, "units", "meters per second")
        if ( fold_limit_lower > -1.E4 .and. fold_limit_upper < 1.E4 ) then
            ierr = nf90_put_att ( self%ncid, config_vel%varid, "field_folds", "true")
            ierr = nf90_put_att ( self%ncid, config_vel%varid, "fold_limit_lower", fold_limit_lower)
            ierr = nf90_put_att ( self%ncid, config_vel%varid, "fold_limit_upper", fold_limit_upper)
        else
            ierr = nf90_put_att ( self%ncid, config_vel%varid, "field_folds", "false")
        endif
        ierr = nf90_put_att ( self%ncid, config_vel%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_vel%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_vel_noise_added%output) then
        ierr = nf90_def_var ( self%ncid, "VEL_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_vel_noise_added%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "standard_name", "radial_velocity_noise_added")
        ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "long_name", "radial_velocity_noise_added")
        ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "units", "meters per second")
        if ( fold_limit_lower > -1.E4 .and. fold_limit_upper < 1.E4 ) then
            ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "field_folds", "true")
            ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "fold_limit_lower", fold_limit_lower)
            ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "fold_limit_upper", fold_limit_upper)
        else
            ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "field_folds", "false")
        endif
        ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_vel_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_stdv_vel%output) then
        ierr = nf90_def_var ( self%ncid, "STDV_VEL", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_stdv_vel%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_stdv_vel%varid, "standard_name", "standard deviation of radial velocity")
        ierr = nf90_put_att ( self%ncid, config_stdv_vel%varid, "long_name", "standard deviation of radial velocity")
        ierr = nf90_put_att ( self%ncid, config_stdv_vel%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_stdv_vel%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_stdv_vel%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_veltru%output) then
        ! "True" radial velocity (i.e., without aircraft motion considered)
        ierr = nf90_def_var ( self%ncid, "VELTRU", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_veltru%varid, "standard_name", "radial_velocity_not_considering_aircraft_motion")
        ierr = nf90_put_att ( self%ncid, config_veltru%varid, "long_name", "true_radial_velocity")
        ierr = nf90_put_att ( self%ncid, config_veltru%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_veltru%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_veltru%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_veltru%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_veltru%output) then
    ! "True" radial velocity (i.e., without aircraft motion considered), fall velocity included
       ierr = nf90_def_var ( self%ncid, "VELTRU_FVEL", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru_fvel%varid, &
         &                deflate_level=2, shuffle=.true. )
       ierr = nf90_put_att ( self%ncid, config_veltru_fvel%varid, "standard_name", "radial_velocity_not_considering_aircraft_motion_with_fall_velocity")
       ierr = nf90_put_att ( self%ncid, config_veltru_fvel%varid, "long_name", "true_radial_velocity_with_fvel")
       ierr = nf90_put_att ( self%ncid, config_veltru_fvel%varid, "units", "meters per second")
       ierr = nf90_put_att ( self%ncid, config_veltru_fvel%varid, "field_folds", "false")
       ierr = nf90_put_att ( self%ncid, config_veltru_fvel%varid, "_FillValue", -9999.0_4 )
       ierr = nf90_put_att ( self%ncid, config_veltru_fvel%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_veltru_noise_added%output) then
        ! "True" radial velocity (i.e., without aircraft motion considered)
        ierr = nf90_def_var ( self%ncid, "VELTRU_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru_noise_added%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_veltru_noise_added%varid, "standard_name", "radial_velocity_noise_added_not_considering_aircraft_motion")
        ierr = nf90_put_att ( self%ncid, config_veltru_noise_added%varid, "long_name", "true_radial_velocity_noise_added")
        ierr = nf90_put_att ( self%ncid, config_veltru_noise_added%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_veltru_noise_added%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_veltru_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_veltru_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    !!! Added by B. Klotz (10/11/2021)
    if (config_veltru_fvel_noise_added%output) then
    ! "True" radial velocity (i.e., without aircraft motion considered)
        ierr = nf90_def_var ( self%ncid, "VELTRU_FVEL_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru_fvel_noise_added%varid, &
         &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_veltru_fvel_noise_added%varid, "standard_name", &
             &   "radial_velocity_with_fall_velocity_and_noise_added_not_considering_aircraft_motion")
        ierr = nf90_put_att ( self%ncid, config_veltru_fvel_noise_added%varid, "long_name", "true_radial_velocity__with_fvel_noise_added")
        ierr = nf90_put_att ( self%ncid, config_veltru_fvel_noise_added%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_veltru_fvel_noise_added%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_veltru_fvel_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_veltru_fvel_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_veltru_folded%output) then
        ierr = nf90_def_var ( self%ncid, "VELTRUF", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru_folded%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded%varid, "standard_name", "folded_radial_velocity_not_considering_aircraft_motion")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded%varid, "long_name", "folded_true_radial_velocity")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_veltru_folded_noise_added%output) then
        ierr = nf90_def_var ( self%ncid, "VELTRUF_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru_folded_noise_added%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_noise_added%varid, "standard_name", "folded_noise_added_radial_velocity_not_considering_aircraft_motion")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_noise_added%varid, "long_name", "folded_noise_added_true_radial_velocity")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_noise_added%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_noise_added%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    !!! Added by B. Klotz (10/11/2021)
    if (config_veltru_folded_fvel%output) then
        ierr = nf90_def_var ( self%ncid, "VELTRUF_FVEL", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru_folded_fvel%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel%varid, "standard_name", "folded_radial_velocity_with_fall_velocity_not_considering_aircraft_motion")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel%varid, "long_name", "folded_true_radial_velocity_with_fvel")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    !!! Added by B. Klotz (10/11/2021)
    if (config_veltru_folded_fvel_noise_added%output) then
        ierr = nf90_def_var ( self%ncid, "VELTRUF_FVEL_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_veltru_folded_fvel_noise_added%varid, &
            &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel_noise_added%varid, "standard_name", &
            &       "folded_noise_added_radial_velocity_with_fall_velocity_not_considering_aircraft_motion")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel_noise_added%varid, "long_name", "folded_noise_added_true_radial_velocity")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel_noise_added%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel_noise_added%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_veltru_folded_fvel_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

!!!! Added by B. Klotz (10/05/2021)
    !if (config_dvh%output) then
         ! Mean Doppler velocity from CRSIM - w is considered 0
    !    ierr = nf90_def_var ( self%ncid, "DVh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_dvh%varid, &
    !         &                deflate_level=2, shuffle=.true. )
    !    ierr = nf90_put_att ( self%ncid, config_dvh%varid, "standard_name", "Mean Doppler velocity considering hydrometeor fall velocity")
    !    ierr = nf90_put_att ( self%ncid, config_dvh%varid, "long_name", "mean_doppler_velocity")
    !    ierr = nf90_put_att ( self%ncid, config_dvh%varid, "units", "meters per second")
    !    ierr = nf90_put_att ( self%ncid, config_dvh%varid, "field_folds", "false")
    !    ierr = nf90_put_att ( self%ncid, config_dvh%varid, "_FillValue", -9999.0_4 )
    !    ierr = nf90_put_att ( self%ncid, config_dvh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    !endif

    !if (config_dopp%output) then
     ! Doppler velocity from CRSIM - w is not set to 0
    !    ierr = nf90_def_var ( self%ncid, "Dopp", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_dopp%varid, &
    !         &                deflate_level=2, shuffle=.true. )
    !    ierr = nf90_put_att ( self%ncid, config_dopp%varid, "standard_name", "Doppler velocity considering hydrometeor fall velocity")
    !    ierr = nf90_put_att ( self%ncid, config_dopp%varid, "long_name", "doppler_velocity")
    !    ierr = nf90_put_att ( self%ncid, config_dopp%varid, "units", "meters per second")
    !    ierr = nf90_put_att ( self%ncid, config_dopp%varid, "field_folds", "false")
    !    ierr = nf90_put_att ( self%ncid, config_dopp%varid, "_FillValue", -9999.0_4 )
    !    ierr = nf90_put_att ( self%ncid, config_dopp%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    !endif

    if (config_fvel%output) then
    ! Hydrometeor fall speed from CRSIM
        ierr = nf90_def_var ( self%ncid, "fvel", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_fvel%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_fvel%varid, "standard_name", "Hydrometeor fall velocity")
        ierr = nf90_put_att ( self%ncid, config_fvel%varid, "long_name", "fall_velocity")
        ierr = nf90_put_att ( self%ncid, config_fvel%varid, "units", "meters per second")
        ierr = nf90_put_att ( self%ncid, config_fvel%varid, "field_folds", "false")
        ierr = nf90_put_att ( self%ncid, config_fvel%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_fvel%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
!!!! End addition by B. Klotz

    if (config_zhh%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh%varid, "standard_name", "Zhh")
        ierr = nf90_put_att ( self%ncid, config_zhh%varid, "long_name", "Reflectivity at hh polarization")
        ierr = nf90_put_att ( self%ncid, config_zhh%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

 !   if (config_zhh%output) then
 !       ierr = nf90_def_var ( self%ncid, "DBZ", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh%varid, &
 !            &                deflate_level=2, shuffle=.true. )
 !       ierr = nf90_put_att ( self%ncid, config_zhh%varid, "standard_name", "DBZ")
 !       ierr = nf90_put_att ( self%ncid, config_zhh%varid, "long_name", "Reflectivity at hh polarization")
 !       ierr = nf90_put_att ( self%ncid, config_zhh%varid, "units", "dBZ")
 !       ierr = nf90_put_att ( self%ncid, config_zhh%varid, "_FillValue", -9999.0_4 )
 !       ierr = nf90_put_att ( self%ncid, config_zhh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
 !   endif
    
    if (config_zhh_attenuated%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh_attenuated", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh_Attenuated%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh_Attenuated%varid, "standard_name", "Zhh_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zhh_Attenuated%varid, "long_name", "Zhh_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zhh_Attenuated%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh_Attenuated%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh_Attenuated%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_stdv_zhh%output) then
        ierr = nf90_def_var ( self%ncid, "STDV_Zhh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_stdv_zhh%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_stdv_zhh%varid, "standard_name", "STDV_Zhh")
        ierr = nf90_put_att ( self%ncid, config_stdv_zhh%varid, "long_name", "Standard Deviation of Zhh")
        ierr = nf90_put_att ( self%ncid, config_stdv_zhh%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_stdv_zhh%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_stdv_zhh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_zhh_noise_added%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh_noise_added%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh_noise_added%varid, "standard_name", "Zhh_noise_added")
        ierr = nf90_put_att ( self%ncid, config_zhh_noise_added%varid, "long_name", "Zhh_attenuated with added noise")
        ierr = nf90_put_att ( self%ncid, config_zhh_noise_added%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_zhh_cloud%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh_cloud", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh_cloud%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh_cloud%varid, "standard_name", "Zhh cloud")
        ierr = nf90_put_att ( self%ncid, config_zhh_cloud%varid, "long_name", "Zhh contribution from cloud")
        ierr = nf90_put_att ( self%ncid, config_zhh_cloud%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh_cloud%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh_cloud%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zhh_rain%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh_rain", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh_rain%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh_rain%varid, "standard_name", "Zhh rain")
        ierr = nf90_put_att ( self%ncid, config_zhh_rain%varid, "long_name", "Zhh contribution from rain")
        ierr = nf90_put_att ( self%ncid, config_zhh_rain%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh_rain%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh_rain%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zhh_ice%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh_ice", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh_ice%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh_ice%varid, "standard_name", "Zhh ice")
        ierr = nf90_put_att ( self%ncid, config_zhh_ice%varid, "long_name", "Zhh contribution from ice")
        ierr = nf90_put_att ( self%ncid, config_zhh_ice%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh_ice%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh_ice%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zhh_snow%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh_snow", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh_snow%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh_snow%varid, "standard_name", "Zhh snow")
        ierr = nf90_put_att ( self%ncid, config_zhh_snow%varid, "long_name", "Zhh contribution from snow")
        ierr = nf90_put_att ( self%ncid, config_zhh_snow%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh_snow%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh_snow%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zhh_graupel%output) then
        ierr = nf90_def_var ( self%ncid, "Zhh_graupel", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zhh_graupel%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zhh_graupel%varid, "standard_name", "Zhh graupel")
        ierr = nf90_put_att ( self%ncid, config_zhh_graupel%varid, "long_name", "Zhh contribution from graupel")
        ierr = nf90_put_att ( self%ncid, config_zhh_graupel%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zhh_graupel%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zhh_graupel%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zvv%output) then
        ierr = nf90_def_var ( self%ncid, "Zvv", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zvv%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zvv%varid, "standard_name", "Zvv")
        ierr = nf90_put_att ( self%ncid, config_zvv%varid, "long_name", "Reflectivity at vv polarization")
        ierr = nf90_put_att ( self%ncid, config_zvv%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zvv%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zvv%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zvv_attenuated%output) then
        ierr = nf90_def_var ( self%ncid, "Zvv_attenuated", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zvv_Attenuated%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zvv_Attenuated%varid, "standard_name", "Zvv_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zvv_Attenuated%varid, "long_name", "Zvv_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zvv_Attenuated%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zvv_Attenuated%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zvv_Attenuated%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_stdv_zvv%output) then
        ierr = nf90_def_var ( self%ncid, "STDV_Zvv", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_stdv_zvv%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_stdv_zvv%varid, "standard_name", "STDV_Zvv")
        ierr = nf90_put_att ( self%ncid, config_stdv_zvv%varid, "long_name", "Standard Deviation of Zvv")
        ierr = nf90_put_att ( self%ncid, config_stdv_zvv%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_stdv_zvv%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_stdv_zvv%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_zvv_noise_added%output) then
        ierr = nf90_def_var ( self%ncid, "Zvv_noise_added", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zvv_noise_added%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zvv_noise_added%varid, "standard_name", "Zvv_noise_added")
        ierr = nf90_put_att ( self%ncid, config_zvv_noise_added%varid, "long_name", "Zvv_attenuated with added noise")
        ierr = nf90_put_att ( self%ncid, config_zvv_noise_added%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zvv_noise_added%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zvv_noise_added%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_zvh%output) then
        ierr = nf90_def_var ( self%ncid, "Zvh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zvh%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zvh%varid, "standard_name", "Zvh")
        ierr = nf90_put_att ( self%ncid, config_zvh%varid, "long_name", "Reflectivity at vh polarization")
        ierr = nf90_put_att ( self%ncid, config_zvh%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zvh%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zvh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zvh_attenuated%output) then
        ierr = nf90_def_var ( self%ncid, "Zvh_attenuated", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zvh_attenuated%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zvh_attenuated%varid, "standard_name", "Zvh_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zvh_attenuated%varid, "long_name", "Zvh_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zvh_attenuated%varid, "units", "dBZ")
        ierr = nf90_put_att ( self%ncid, config_zvh_attenuated%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zvh_attenuated%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_rhohv%output) then
        ierr = nf90_def_var ( self%ncid, "RHOhv", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_rhohv%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_rhohv%varid, "standard_name", "RHOhv")
        ierr = nf90_put_att ( self%ncid, config_rhohv%varid, "long_name", "Cross-correlation coefficient")
        ierr = nf90_put_att ( self%ncid, config_rhohv%varid, "units", "-")
        ierr = nf90_put_att ( self%ncid, config_rhohv%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_rhohv%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zdr%output) then
        ierr = nf90_def_var ( self%ncid, "Zdr", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zdr%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zdr%varid, "standard_name", "Zdr")
        ierr = nf90_put_att ( self%ncid, config_zdr%varid, "long_name", "Differential reflectivity")
        ierr = nf90_put_att ( self%ncid, config_zdr%varid, "units", "dB")
        ierr = nf90_put_att ( self%ncid, config_zdr%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zdr%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_zdr_attenuated%output) then
        ierr = nf90_def_var ( self%ncid, "Zdr_attenuated", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_zdr_attenuated%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_zdr_attenuated%varid, "standard_name", "Zdr_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zdr_attenuated%varid, "long_name", "Zdr_attenuated")
        ierr = nf90_put_att ( self%ncid, config_zdr_attenuated%varid, "units", "dB")
        ierr = nf90_put_att ( self%ncid, config_zdr_attenuated%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_zdr_attenuated%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_ldrh%output) then
        ierr = nf90_def_var ( self%ncid, "LDRh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_ldrh%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_ldrh%varid, "standard_name", "LDRh")
        ierr = nf90_put_att ( self%ncid, config_ldrh%varid, "long_name", "log_linear_deolarization_ratio_h")
        ierr = nf90_put_att ( self%ncid, config_ldrh%varid, "units", "dB")
        ierr = nf90_put_att ( self%ncid, config_ldrh%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_ldrh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_ah%output) then
        ierr = nf90_def_var ( self%ncid, "Ah", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_ah%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_ah%varid, "standard_name", "Ah")
        ierr = nf90_put_att ( self%ncid, config_ah%varid, "long_name", "Specific horizontal attenuation")
        ierr = nf90_put_att ( self%ncid, config_ah%varid, "units", "dB/km")
        ierr = nf90_put_att ( self%ncid, config_ah%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_ah%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_av%output) then
        ierr = nf90_def_var ( self%ncid, "Av", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_av%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_av%varid, "standard_name", "Av")
        ierr = nf90_put_att ( self%ncid, config_av%varid, "long_name", "Specific vertical attenuation")
        ierr = nf90_put_att ( self%ncid, config_av%varid, "units", "dB/km")
        ierr = nf90_put_att ( self%ncid, config_av%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_av%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_Adp%output) then
        ierr = nf90_def_var ( self%ncid, "Adp", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_Adp%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_Adp%varid, "standard_name", "Adp")
        ierr = nf90_put_att ( self%ncid, config_Adp%varid, "long_name", "Differential attenuation")
        ierr = nf90_put_att ( self%ncid, config_Adp%varid, "units", "dB/km")
        ierr = nf90_put_att ( self%ncid, config_Adp%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_Adp%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_Kdp%output) then
        ierr = nf90_def_var ( self%ncid, "Kdp", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_kdp%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_kdp%varid, "standard_name", "Kdp")
        ierr = nf90_put_att ( self%ncid, config_kdp%varid, "long_name", "Specific differential phase")
        ierr = nf90_put_att ( self%ncid, config_kdp%varid, "units", "deg/km")
        ierr = nf90_put_att ( self%ncid, config_kdp%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_kdp%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_DPhase%output) then
        ierr = nf90_def_var ( self%ncid, "DPhase", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_dphase%varid, &
         &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_dphase%varid, "standard_name", "DPhase")
        ierr = nf90_put_att ( self%ncid, config_dphase%varid, "long_name", "Differential backscatter phase")
        ierr = nf90_put_att ( self%ncid, config_dphase%varid, "units", "deg")
        ierr = nf90_put_att ( self%ncid, config_dphase%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_dphase%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_PHIdp%output) then
        ierr = nf90_def_var ( self%ncid, "PHIdp", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_phidp%varid, &
         &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_phidp%varid, "standard_name", "PHIdp")
        ierr = nf90_put_att ( self%ncid, config_phidp%varid, "long_name", "Total differential phase")
        ierr = nf90_put_att ( self%ncid, config_phidp%varid, "units", "deg")
        ierr = nf90_put_att ( self%ncid, config_phidp%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_phidp%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_PHIdp_wDP%output) then
        ierr = nf90_def_var ( self%ncid, "PHIdp_DPh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_phidp_wdp%varid, &
         &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_phidp_wdp%varid, "standard_name", "PHIdp_DPh")
        ierr = nf90_put_att ( self%ncid, config_phidp_wdp%varid, "long_name", "Total differential phase with backscatter phase")
        ierr = nf90_put_att ( self%ncid, config_phidp_wdp%varid, "units", "deg")
        ierr = nf90_put_att ( self%ncid, config_phidp_wdp%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_phidp_wdp%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_swtot%output) then
        ierr = nf90_def_var ( self%ncid, "SWtot", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_swtot%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_swtot%varid, "standard_name", "SWtot")
        ierr = nf90_put_att ( self%ncid, config_swtot%varid, "long_name", "Total spectrum width")
        ierr = nf90_put_att ( self%ncid, config_swtot%varid, "units", "m/s")
        ierr = nf90_put_att ( self%ncid, config_swtot%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_swtot%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_swh%output) then
        ierr = nf90_def_var ( self%ncid, "SWh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_swh%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_swh%varid, "standard_name", "SWh")
        ierr = nf90_put_att ( self%ncid, config_swh%varid, "long_name", "Spectrum width due to hydrometeors")
        ierr = nf90_put_att ( self%ncid, config_swh%varid, "units", "m/s")
        ierr = nf90_put_att ( self%ncid, config_swh%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_swh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_sws%output) then
        ierr = nf90_def_var ( self%ncid, "SWs", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_sws%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_sws%varid, "standard_name", "SWs")
        ierr = nf90_put_att ( self%ncid, config_sws%varid, "long_name", "Spectrum width due to wind shear")
        ierr = nf90_put_att ( self%ncid, config_sws%varid, "units", "m/s")
        ierr = nf90_put_att ( self%ncid, config_sws%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_sws%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_swt%output) then
        ierr = nf90_def_var ( self%ncid, "SWt", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_swt%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_swt%varid, "standard_name", "SWt")
        ierr = nf90_put_att ( self%ncid, config_swt%varid, "long_name", "Spectrum width due to turbulence")
        ierr = nf90_put_att ( self%ncid, config_swt%varid, "units", "m/s")
        ierr = nf90_put_att ( self%ncid, config_swt%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_swt%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_swv%output) then
        ierr = nf90_def_var ( self%ncid, "SWv", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_swv%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_swv%varid, "standard_name", "SWv")
        ierr = nf90_put_att ( self%ncid, config_swv%varid, "long_name", "Spectrum width due to cross wind")
        ierr = nf90_put_att ( self%ncid, config_swv%varid, "units", "m/s")
        ierr = nf90_put_att ( self%ncid, config_swv%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_swv%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_snr_zhh%output) then
        ierr = nf90_def_var ( self%ncid, "SNR_zhh", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_snr_zhh%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_snr_zhh%varid, "standard_name", "SNR_zhh")
        ierr = nf90_put_att ( self%ncid, config_snr_zhh%varid, "long_name", "signal-to-noise ratio for Zhh")
        ierr = nf90_put_att ( self%ncid, config_snr_zhh%varid, "units", "dB")
        ierr = nf90_put_att ( self%ncid, config_snr_zhh%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_snr_zhh%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif
    
    if (config_snr_zvv%output) then
        ierr = nf90_def_var ( self%ncid, "SNR_zvv", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_snr_zvv%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_snr_zvv%varid, "standard_name", "SNR_zvv")
        ierr = nf90_put_att ( self%ncid, config_snr_zvv%varid, "long_name", "signal-to-noise ratio for Zvv")
        ierr = nf90_put_att ( self%ncid, config_snr_zvv%varid, "units", "dB")
        ierr = nf90_put_att ( self%ncid, config_snr_zvv%varid, "_FillValue", -9999.0_4 )
        ierr = nf90_put_att ( self%ncid, config_snr_zvv%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_x%output) then
        ierr = nf90_def_var ( self%ncid, "x", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_x%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_x%varid, "standard_name", "x")
        ierr = nf90_put_att ( self%ncid, config_x%varid, "long_name", "x")
        ierr = nf90_put_att ( self%ncid, config_x%varid, "units", "-")
        ierr = nf90_put_att ( self%ncid, config_x%varid, "_FillValue", -1.E36_4 )
        ierr = nf90_put_att ( self%ncid, config_x%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_y%output) then
        ierr = nf90_def_var ( self%ncid, "y", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_y%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_y%varid, "standard_name", "y")
        ierr = nf90_put_att ( self%ncid, config_y%varid, "long_name", "y")
        ierr = nf90_put_att ( self%ncid, config_y%varid, "units", "-")
        ierr = nf90_put_att ( self%ncid, config_y%varid, "_FillValue", -1.E36_4 )
        ierr = nf90_put_att ( self%ncid, config_y%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    if (config_z%output) then
        ierr = nf90_def_var ( self%ncid, "z", NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=config_z%varid, &
             &                deflate_level=2, shuffle=.true. )
        ierr = nf90_put_att ( self%ncid, config_z%varid, "standard_name", "z")
        ierr = nf90_put_att ( self%ncid, config_z%varid, "long_name", "z")
        ierr = nf90_put_att ( self%ncid, config_z%varid, "units", "-")
        ierr = nf90_put_att ( self%ncid, config_z%varid, "_FillValue", -1.E36_4 )
        ierr = nf90_put_att ( self%ncid, config_z%varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    endif

    ierr = nf90_enddef ( self%ncid )
    call error_handler(ierr, "(CFRADIAL_OPEN) Problem enddef")

  end subroutine cfradial_open

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine cfradial_prepare_metadata(self)
    implicit none
    class(cfradial_type), intent(inout) :: self
    type(volume_field_type), pointer :: ptr

    ptr => root
    do while ( associated ( ptr%next ) )
        ptr => ptr%next
        if ( ptr%source == "COMPUTE" ) then
            ! skip
        elseif ( ptr%source == "WRF" ) then
            ptr%varid = -99999
            if ( ptr%field_name == "TEMPERATURE" .and. .not. config_WRF_TEMPERATURE%output ) cycle
            if ( ptr%field_name == "HT"          .and. .not. config_WRF_HT%output ) cycle
            if ( ptr%field_name == "RHO_D"       .and. .not. config_WRF_RHO_D%output ) cycle
            if ( ptr%field_name == "U"           .and. .not. config_WRF_U%output ) cycle
            if ( ptr%field_name == "V"           .and. .not. config_WRF_V%output ) cycle
            if ( ptr%field_name == "W"           .and. .not. config_WRF_W%output ) cycle
            if ( ptr%field_name == "QCLOUD"      .and. .not. config_WRF_QCLOUD%output ) cycle
            if ( ptr%field_name == "QRAIN"       .and. .not. config_WRF_QRAIN%output ) cycle
            if ( ptr%field_name == "QICE"        .and. .not. config_WRF_QICE%output ) cycle
            if ( ptr%field_name == "QSNOW"       .and. .not. config_WRF_QSNOW%output ) cycle
            if ( ptr%field_name == "QGRAUP"      .and. .not. config_WRF_QGRAUP%output ) cycle
            if ( ptr%field_name == "QNRAIN"      .and. .not. config_WRF_QNRAIN%output ) cycle
            if ( ptr%field_name == "QNICE"       .and. .not. config_WRF_QNICE%output ) cycle
            if ( ptr%field_name == "QNSNOW"      .and. .not. config_WRF_QNSNOW%output ) cycle
            if ( ptr%field_name == "QNGRAUPEL"   .and. .not. config_WRF_QNGRAUPEL%output ) cycle
            call self%add_to_output_metadata("WRF_"//ptr%field_name, ptr%field_description, ptr%field_units, ptr%varid )
        else
            call self%add_to_output_metadata(ptr%field_name, ptr%field_description, ptr%field_units, ptr%varid )
        endif
    enddo

  end subroutine cfradial_prepare_metadata

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine cfradial_add_to_output_metadata ( self, field_name, field_description, field_units, varid )
    use netcdf, only : NF90_FLOAT
    use netcdf, only : nf90_redef
    use netcdf, only : nf90_inq_dimid
    use netcdf, only : nf90_def_var
    use netcdf, only : nf90_put_att
    use netcdf, only : nf90_enddef
    implicit none
    class (cfradial_type), intent(in) :: self
    character(len=*), intent(in) :: field_name, field_description, field_units
    integer, intent(out) :: varid
    integer :: ierr, range_dimid, time_dimid
    character(len=33), parameter :: here = "CFRADIAL_ADD_TO_OUTPUT_METADATA: "

    ierr = nf90_inq_dimid ( self%ncid, "range", range_dimid )
    call error_handler(ierr, "Problem get range_dimid")

    ierr = nf90_inq_dimid ( self%ncid, "time", time_dimid )
    call error_handler(ierr, "Problem get time_dimid")

    ierr = nf90_redef(self%ncid)
    call error_handler(ierr, "Problem redef")

    ierr = nf90_def_var ( self%ncid, trim(field_name), NF90_FLOAT, dimids=(/range_dimid,time_dimid/), varid=varid, &
         &                deflate_level=2, shuffle=.true. )
    call error_handler(ierr, here//"Problem define variable '"//trim(field_name)//"'")
    ierr = nf90_put_att ( self%ncid, varid, "standard_name", field_name)
    call error_handler(ierr, here//"Problem put attribute standard_name for variable '"//trim(field_name)//"'")
    ierr = nf90_put_att ( self%ncid, varid, "long_name", trim(field_description) )
    call error_handler(ierr, here//"Problem put attribute long_name for variable '"//trim(field_name)//"'")
    ierr = nf90_put_att ( self%ncid, varid, "units", trim(field_units) )
    call error_handler(ierr, here//"Problem put attribute units for variable '"//trim(field_name)//"'")
    ierr = nf90_put_att ( self%ncid, varid, "_FillValue", -9999.0_4 )
    call error_handler(ierr, here//"Problem put attribute  _FillValue for variable '"//trim(field_name)//"'")
    ierr = nf90_put_att ( self%ncid, varid, "coordinates", "elevation azimuth range heading roll pitch rotation tilt" )
    call error_handler(ierr, here//"Problem put attribute coordinates for variable '"//trim(field_name)//"'")


    ierr = nf90_enddef ( self%ncid )
    call error_handler(ierr, "Problem enddef")

  end subroutine cfradial_add_to_output_metadata

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine write_volume ( self , volume , primary_axis )
    use netcdf, only : nf90_put_var
    use kwm_date_utilities, only : geth_newdate
    implicit none
    class (cfradial_type), intent(inout) :: self
    type(volume_type) :: volume
    character(len=*), intent(in) :: primary_axis
    integer :: ierr
    type(volume_field_type), pointer :: ptr
    real(kind=RKIND), pointer, dimension(:,:) :: arr
    ! print*, 'BEGIN write volume'

    ierr = nf90_put_var(self%ncid, self%volume_number_varid, volume%number)
    call error_handler(ierr, "Problem put volume_number")

    ierr = nf90_put_var(self%ncid, self%platform_type_varid, "aircraft_fore")
    call error_handler(ierr, "Problem put platform_type")

    ierr = nf90_put_var(self%ncid, self%instrument_type_varid, "radar")
    call error_handler(ierr, "Problem put instrument_type")

    select case (primary_axis)
    case default
        write(*,'("CFRADIAL:  WRITE_VOLUME: Primary axis: ", A)') trim(primary_axis)
        stop "CFRADIAL:  WRITE_VOLUME: Unrecognized primary axis"
    case ("X")
        ierr = nf90_put_var(self%ncid, self%primary_axis_varid, "axis_x")
        call error_handler(ierr, "Problem put primary_axis")
    case ("X-Prime")
        ierr = nf90_put_var(self%ncid, self%primary_axis_varid, "axis_x_prime")
        call error_handler(ierr, "Problem put primary_axis")
    case ("Y")
        ierr = nf90_put_var(self%ncid, self%primary_axis_varid, "axis_y")
        call error_handler(ierr, "Problem put primary_axis")
    case ("Y-Prime")
        ierr = nf90_put_var(self%ncid, self%primary_axis_varid, "axis_y_prime")
        call error_handler(ierr, "Problem put primary_axis")
    case ("Z")
        ierr = nf90_put_var(self%ncid, self%primary_axis_varid, "axis_z")
        call error_handler(ierr, "Problem put primary_axis")
    case ("Z-Prime")
        ierr = nf90_put_var(self%ncid, self%primary_axis_varid, "axis_z_prime")
        call error_handler(ierr, "Problem put primary_axis")
    end select
        

    ierr = nf90_put_var(self%ncid, self%time_coverage_start_varid, volume%time_coverage_start//"Z" )
    call error_handler(ierr, "Problem put time_coverage_start")

    ! Compute <time_coverage_end> based on the time of the last beam in the volume.
    ! call geth_newdate ( volume%time_coverage_end, volume%reference_time, volume%start_time_integer_seconds + int(volume%time(volume%nrays)) )
    call geth_newdate ( volume%time_coverage_end, volume%reference_time, int(volume%time(volume%nrays)) )
    print*, "time_coverage_end = ", volume%reference_time,  " + ", int(volume%time(volume%nrays)), " = ", volume%time_coverage_end

    ierr = nf90_put_var(self%ncid, self%time_coverage_end_varid, volume%time_coverage_end//"Z" )
    call error_handler(ierr, "Problem put time_coverage_end")

    ierr = nf90_put_var(self%ncid, self%time_varid, volume%time(1:volume%nrays))
    call error_handler(ierr, "Problem put time")

    ierr = nf90_put_var(self%ncid, self%range_varid, volume%range(1:volume%ngates))
    call error_handler(ierr, "Problem put range")

    ierr = nf90_put_var(self%ncid, self%latitude_varid, volume%latitude(1:volume%nrays))
    call error_handler(ierr, "Problem put latitude")

    ierr = nf90_put_var(self%ncid, self%longitude_varid, volume%longitude(1:volume%nrays))
    call error_handler(ierr, "Problem put longitude")

    ierr = nf90_put_var(self%ncid, self%altitude_varid, volume%altitude(1:volume%nrays))
    call error_handler(ierr, "Problem put altitude")

    ! Put altitude_agl == altitude for our simulated no-terrain cases.
    ierr = nf90_put_var(self%ncid, self%altitude_agl_varid, volume%altitude(1:volume%nrays))
    call error_handler(ierr, "Problem put altitude_agl")

    ierr = nf90_put_var(self%ncid, self%sweep_number_varid, volume%sweep_number(1:volume%sweep))
    call error_handler(ierr, "Problem put sweep_number")

    ierr = nf90_put_var(self%ncid, self%sweep_mode_varid, volume%sweep_mode(1:volume%sweep))
    call error_handler(ierr, "Problem put sweep_mode")

    ierr = nf90_put_var(self%ncid, self%fixed_angle_varid, volume%fixed_angle(1:volume%sweep))
    call error_handler(ierr, "Problem put fixed_angle")

    ierr = nf90_put_var(self%ncid, self%sweep_start_ray_index_varid, volume%sweep_start_ray_index(1:volume%sweep))
    call error_handler(ierr, "Problem put sweep_start_ray_index")

    ierr = nf90_put_var(self%ncid, self%sweep_end_ray_index_varid, volume%sweep_end_ray_index(1:volume%sweep))
    call error_handler(ierr, "Problem put sweep_end_ray_index")

    ierr = nf90_put_var(self%ncid, self%azimuth_varid, volume%azimuth(1:volume%nrays))
    call error_handler(ierr, "Problem put azimuth")

    ierr = nf90_put_var(self%ncid, self%elevation_varid, volume%elevation(1:volume%nrays))
    call error_handler(ierr, "Problem put elevation")

    ierr = nf90_put_var(self%ncid, self%georefs_applied_varid, volume%georefs_applied(1:volume%nrays))
    call error_handler(ierr, "Problem put georefs_applied")

    ierr = nf90_put_var(self%ncid, self%heading_varid, volume%heading(1:volume%nrays))
    call error_handler(ierr, "Problem put heading")

    ierr = nf90_put_var(self%ncid, self%roll_varid, volume%roll(1:volume%nrays))
    call error_handler(ierr, "Problem put roll")

    ierr = nf90_put_var(self%ncid, self%pitch_varid, volume%pitch(1:volume%nrays))
    call error_handler(ierr, "Problem put pitch")

    ierr = nf90_put_var(self%ncid, self%drift_varid, volume%drift(1:volume%nrays))
    call error_handler(ierr, "Problem put drift")

    ierr = nf90_put_var(self%ncid, self%rotation_varid, volume%rotation(1:volume%nrays))
    call error_handler(ierr, "Problem put rotation")

    ierr = nf90_put_var(self%ncid, self%tilt_varid, volume%tilt(1:volume%nrays))
    call error_handler(ierr, "Problem put tilt")

    ierr = nf90_put_var(self%ncid, self%nyquist_velocity_varid, volume%nyquist_velocity(1:volume%nrays))
    call error_handler(ierr, "Problem put nyquist_velocity")

    ierr = nf90_put_var(self%ncid, self%unambiguous_range_varid, volume%unambiguous_range(1:volume%nrays))
    call error_handler(ierr, "Problem put unambiguous_range")

    ! geometry_correction meta_group:

    ierr = nf90_put_var(self%ncid, self%azimuth_correction_varid, volume%azimuth_correction)
    call error_handler(ierr, "Problem put azimuth_correction")

    ierr = nf90_put_var(self%ncid, self%elevation_correction_varid, volume%elevation_correction)
    call error_handler(ierr, "Problem put elevation_correction")

    ierr = nf90_put_var(self%ncid, self%range_correction_varid, volume%range_correction)
    call error_handler(ierr, "Problem put range_correction")

    ierr = nf90_put_var(self%ncid, self%longitude_correction_varid, volume%longitude_correction)
    call error_handler(ierr, "Problem put longitude_correction")

    ierr = nf90_put_var(self%ncid, self%latitude_correction_varid, volume%latitude_correction)
    call error_handler(ierr, "Problem put latitude_correction")

    ierr = nf90_put_var(self%ncid, self%pressure_altitude_correction_varid, volume%pressure_altitude_correction)
    call error_handler(ierr, "Problem put pressure_altitude_correction")

    ierr = nf90_put_var(self%ncid, self%radar_altitude_correction_varid, volume%radar_altitude_correction)
    call error_handler(ierr, "Problem put radar_altitude_correction")

    ierr = nf90_put_var(self%ncid, self%eastward_ground_speed_correction_varid, volume%eastward_ground_speed_correction)
    call error_handler(ierr, "Problem put eastward_ground_speed_correction")

    ierr = nf90_put_var(self%ncid, self%northward_ground_speed_correction_varid, volume%northward_ground_speed_correction)
    call error_handler(ierr, "Problem put northward_ground_speed_correction")

    ierr = nf90_put_var(self%ncid, self%vertical_velocity_correction_varid, volume%vertical_velocity_correction)
    call error_handler(ierr, "Problem put vertical_velocity_correction")

    ierr = nf90_put_var(self%ncid, self%heading_correction_varid, volume%heading_correction)
    call error_handler(ierr, "Problem put heading_correction")

    ierr = nf90_put_var(self%ncid, self%roll_correction_varid, volume%roll_correction)
    call error_handler(ierr, "Problem put roll_correction")

    ierr = nf90_put_var(self%ncid, self%pitch_correction_varid, volume%pitch_correction)
    call error_handler(ierr, "Problem put pitch_correction")

    ierr = nf90_put_var(self%ncid, self%drift_correction_varid, volume%drift_correction)
    call error_handler(ierr, "Problem put drift_correction")

    ierr = nf90_put_var(self%ncid, self%rotation_correction_varid, volume%rotation_correction)
    call error_handler(ierr, "Problem put rotation_correction")

    ierr = nf90_put_var(self%ncid, self%tilt_correction_varid, volume%tilt_correction)
    call error_handler(ierr, "Problem put tilt_correction")

    !

    ierr = nf90_put_var(self%ncid, self%eastward_velocity_varid, volume%eastward_velocity(1:volume%nrays))
    call error_handler(ierr, "Problem put eastward_velocity")

    ierr = nf90_put_var(self%ncid, self%northward_velocity_varid, volume%northward_velocity(1:volume%nrays))
    call error_handler(ierr, "Problem put northward_velocity")

    ierr = nf90_put_var(self%ncid, self%vertical_velocity_varid, volume%vertical_velocity(1:volume%nrays))
    call error_handler(ierr, "Problem put vertical_velocity")

    ierr = nf90_put_var(self%ncid, self%eastward_wind_varid, volume%eastward_wind(1:volume%nrays))
    call error_handler(ierr, "Problem put eastward_wind")

    ierr = nf90_put_var(self%ncid, self%northward_wind_varid, volume%northward_wind(1:volume%nrays))
    call error_handler(ierr, "Problem put northward_wind")

    ierr = nf90_put_var(self%ncid, self%vertical_wind_varid, volume%vertical_wind(1:volume%nrays))
    call error_handler(ierr, "Problem put vertical_wind")

    ierr = nf90_put_var(self%ncid, self%heading_rate_varid, volume%heading_rate(1:volume%nrays))
    call error_handler(ierr, "Problem put heading_rate")

    ierr = nf90_put_var(self%ncid, self%roll_rate_varid, volume%roll_rate(1:volume%nrays))
    call error_handler(ierr, "Problem put roll_rate")

    ierr = nf90_put_var(self%ncid, self%pitch_rate_varid, volume%pitch_rate(1:volume%nrays))
    call error_handler(ierr, "Problem put pitch_rate")

    ierr = nf90_put_var ( self%ncid, self%radar_beam_width_h_varid, volume%beamwidth_h )
    call error_handler(ierr, "Problem put radar_beam_width_h" )

    ierr = nf90_put_var ( self%ncid, self%radar_beam_width_v_varid, volume%beamwidth_v )
    call error_handler(ierr, "Problem put radar_beam_width_v" )

    if (config_x%output) then
        arr => volume%point_to_data("VX")
        ierr = nf90_put_var(self%ncid, config_x%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put x")
        nullify(arr)
    endif

    if (config_y%output) then
        arr => volume%point_to_data("VY")
        ierr = nf90_put_var(self%ncid, config_y%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put y")
        nullify(arr)
    endif

    if (config_z%output) then
        arr => volume%point_to_data("VZ")
        ierr = nf90_put_var(self%ncid, config_z%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put z")
        nullify(arr)
    endif

    if (config_vel%output) then
        arr => volume%point_to_data("VEL")
        ierr = nf90_put_var(self%ncid, config_vel%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VEL")
        nullify(arr)
    endif
    
    if (config_vel_noise_added%output) then
        arr => volume%point_to_data("VEL_noise_added")
        ierr = nf90_put_var(self%ncid, config_vel_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VEL_noise_added")
        nullify(arr)
    endif

    if (config_stdv_vel%output) then
        arr => volume%point_to_data("STDV_VEL")
        ierr = nf90_put_var(self%ncid, config_stdv_vel%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put STDV_VEL")
        nullify(arr)
    endif
    
!    if (config_stdv_vel%output) then
!        arr => volume%point_to_data("STDV_VEL")
!        ierr = nf90_put_var(self%ncid, config_stdv_vel%varid, arr(1:volume%ngates,1:volume%nrays))
!        call error_handler(ierr, "Problem put STDV_VEL")
!        nullify(arr)
!    endif
    
    if (config_velunf%output) then
        arr => volume%point_to_data("VELUNF")
        ierr = nf90_put_var(self%ncid, config_velunf%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELUNF")
        nullify(arr)
    endif

    if (config_velunf_fvel%output) then
        arr => volume%point_to_data("VELUNF_FVEL")
        ierr = nf90_put_var(self%ncid, config_velunf_fvel%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELUNF_FVEL")
        nullify(arr)
    endif
    
    if (config_velunf_noise_added%output) then
        arr => volume%point_to_data("VELUNF_noise_added")
        ierr = nf90_put_var(self%ncid, config_velunf_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELUNF_noise_added")
        nullify(arr)
    endif

!!! Added by B. Klotz (10/11/2021)
    if (config_velunf_fvel_noise_added%output) then
        arr => volume%point_to_data("VELUNF_FVEL_noise_added")
        ierr = nf90_put_var(self%ncid, config_velunf_fvel_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELUNF_FVEL_noise_added")
        nullify(arr)
    endif

    if (config_veltru%output) then
        arr => volume%point_to_data("VELTRU")
        ierr = nf90_put_var(self%ncid, config_veltru%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRU")
        nullify(arr)
    endif

    if (config_veltru_fvel%output) then
        arr => volume%point_to_data("VELTRU_FVEL")
        ierr = nf90_put_var(self%ncid, config_veltru_fvel%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRU_FVEL")
        nullify(arr)
    endif

    if (config_veltru_noise_added%output) then
        arr => volume%point_to_data("VELTRU_noise_added")
        ierr = nf90_put_var(self%ncid, config_veltru_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRU_noise_added")
        nullify(arr)
    endif

!!! Added by B. Klotz (10/11/2021)
    if (config_veltru_fvel_noise_added%output) then
        arr => volume%point_to_data("VELTRU_FVEL_noise_added")
        ierr = nf90_put_var(self%ncid, config_veltru_fvel_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRU_FVEL_noise_added")
        nullify(arr)
    endif
    
    if (config_veltru_folded%output) then
        arr => volume%point_to_data("VELTRUF")
        ierr = nf90_put_var(self%ncid, config_veltru_folded%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRUF")
        nullify(arr)
    endif

    if (config_veltru_folded_noise_added%output) then
        arr => volume%point_to_data("VELTRUF_noise_added")
        ierr = nf90_put_var(self%ncid, config_veltru_folded_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRUF_noise_added")
        nullify(arr)
    endif

!!! Added by B. Klotz (10/11/2021)
    if (config_veltru_folded_fvel%output) then
        arr => volume%point_to_data("VELTRUF_FVEL")
        ierr = nf90_put_var(self%ncid, config_veltru_folded_fvel%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRUF_FVEL")
        nullify(arr)
    endif

    if (config_veltru_folded_fvel_noise_added%output) then
        arr => volume%point_to_data("VELTRUF_FVEL_noise_added")
        ierr = nf90_put_var(self%ncid, config_veltru_folded_fvel_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put VELTRUF_FVEL_noise_added")
        nullify(arr)
    endif

!!!! Added by B. Klotz (10/05/2021)
    !if (config_dvh%output) then
    !    arr => volume%point_to_data("DVh")
    !    ierr = nf90_put_var(self%ncid, config_dvh%varid, arr(1:volume%ngates,1:volume%nrays))
    !    call error_handler(ierr, "Problem put DVh")
    !    nullify(arr)
    !endif

    !if (config_dopp%output) then
    !    arr => volume%point_to_data("Dopp")
    !    ierr = nf90_put_var(self%ncid, config_dopp%varid, arr(1:volume%ngates,1:volume%nrays))
    !    call error_handler(ierr, "Problem put Dopp")
    !    nullify(arr)
    !endif

    if (config_fvel%output) then
        arr => volume%point_to_data("fvel")
        ierr = nf90_put_var(self%ncid, config_fvel%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put fvel")
        nullify(arr)
    endif
!!!! End addition by B. Klotz
    
    if (config_zhh%output) then
        arr => volume%point_to_data("Zhh")
        ierr = nf90_put_var(self%ncid, config_zhh%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh")
        nullify(arr)
    endif
    
    if (config_zhh_Attenuated%output) then
        arr => volume%point_to_data("Zhh_attenuated")
        ierr = nf90_put_var(self%ncid, config_zhh_Attenuated%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh_attenuated")
        nullify(arr)
    endif

    if (config_stdv_Zhh%output) then
        arr => volume%point_to_data("STDV_Zhh")
        ierr = nf90_put_var(self%ncid, config_stdv_zhh%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put STDV_Zhh")
        nullify(arr)
    endif

    if (config_stdv_Zvv%output) then
        arr => volume%point_to_data("STDV_Zvv")
        ierr = nf90_put_var(self%ncid, config_stdv_zvv%varid, real(arr(1:volume%ngates,1:volume%nrays),4))
        call error_handler(ierr, "Problem put STDV_Zvv")
        nullify(arr)
    endif
    
    if (config_Zhh_noise_added%output) then
        arr => volume%point_to_data("Zhh_noise_added")
        ierr = nf90_put_var(self%ncid, config_zhh_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh_noise_added")
        nullify(arr)
    endif

    if (config_zhh_cloud%output) then
        arr => volume%point_to_data("Zhh_cloud")
        ierr = nf90_put_var(self%ncid, config_zhh_cloud%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh_cloud")
        nullify(arr)
    endif

    if (config_zhh_rain%output) then
        arr => volume%point_to_data("Zhh_rain")
        ierr = nf90_put_var(self%ncid, config_zhh_rain%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh_rain")
        nullify(arr)
    endif

    if (config_zhh_ice%output) then
        arr => volume%point_to_data("Zhh_ice")
        ierr = nf90_put_var(self%ncid, config_zhh_ice%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh_ice")
        nullify(arr)
    endif

    if (config_zhh_snow%output) then
        arr => volume%point_to_data("Zhh_snow")
        ierr = nf90_put_var(self%ncid, config_zhh_snow%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh_snow")
        nullify(arr)
    endif

    if (config_zhh_graupel%output) then
        arr => volume%point_to_data("Zhh_graupel")
        ierr = nf90_put_var(self%ncid, config_zhh_graupel%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zhh_graupel")
        nullify(arr)
    endif

    if (config_zvv%output) then
        arr => volume%point_to_data("Zvv")
        ierr = nf90_put_var(self%ncid, config_zvv%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zvv")
        nullify(arr)
    endif

    if (config_zvv_attenuated%output) then
        arr => volume%point_to_data("Zvv_attenuated")
        ierr = nf90_put_var(self%ncid, config_zvv_Attenuated%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zvv_attenuated")
        nullify(arr)
    endif

    if (config_Zvv_noise_added%output) then
        arr => volume%point_to_data("Zvv_noise_added")
        ierr = nf90_put_var(self%ncid, config_zvv_noise_added%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zvv_noise_added")
        nullify(arr)
    endif

    if (config_zvh%output) then
        arr => volume%point_to_data("Zvh")
        ierr = nf90_put_var(self%ncid, config_zvh%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zvh")
        nullify(arr)
    endif

    if (config_zvh_attenuated%output) then
        arr => volume%point_to_data("Zvh_attenuated")
        ierr = nf90_put_var(self%ncid, config_zvh_attenuated%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zvh_attenuated")
        nullify(arr)
    endif

    if (config_rhohv%output) then
        arr => volume%point_to_data("RHOhv")
        ierr = nf90_put_var(self%ncid, config_rhohv%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put RHOhv")
        nullify(arr)
    endif

    if (config_zdr%output) then
        arr => volume%point_to_data("Zdr")
        ierr = nf90_put_var(self%ncid, config_zdr%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zdr")
        nullify(arr)
    endif

    if (config_zdr_attenuated%output) then
        arr => volume%point_to_data("Zdr_attenuated")
        ierr = nf90_put_var(self%ncid, config_zdr_attenuated%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Zdr_attenuated")
        nullify(arr)
    endif

    if (config_ldrh%output) then
        arr => volume%point_to_data("LDRh")
        ierr = nf90_put_var(self%ncid, config_ldrh%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put LDRh")
        nullify(arr)
    endif

    if (config_ah%output) then
        arr => volume%point_to_data("Ah")
        ierr = nf90_put_var(self%ncid, config_ah%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Ah")
        nullify(arr)
    endif

    if (config_av%output) then
        arr => volume%point_to_data("Av")
        ierr = nf90_put_var(self%ncid, config_av%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Av")
        nullify(arr)
    endif

    if (config_adp%output) then
        arr => volume%point_to_data("Adp")
        ierr = nf90_put_var(self%ncid, config_Adp%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Adp")
        nullify(arr)
    endif

    if (config_kdp%output) then
        arr => volume%point_to_data("Kdp")
        ierr = nf90_put_var(self%ncid, config_kdp%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put Kdp")
        nullify(arr)
    endif

    if (config_dphase%output) then
        arr => volume%point_to_data("DPhase")
        ierr = nf90_put_var(self%ncid, config_dphase%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put DPhase")
        nullify(arr)
    endif

    if (config_phidp%output) then
        arr => volume%point_to_data("PHIdp")
        ierr = nf90_put_var(self%ncid, config_phidp%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put PHIdp")
        nullify(arr)
    endif

    if (config_phidp_wdp%output) then
        arr => volume%point_to_data("PHIdp_DPh")
        ierr = nf90_put_var(self%ncid, config_phidp_wdp%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put PHIdp_DPh")
        nullify(arr)
    endif

    if (config_swtot%output) then
        arr => volume%point_to_data("SWtot")
        ierr = nf90_put_var(self%ncid, config_swtot%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put SWtot")
        nullify(arr)
    endif

    if (config_sws%output) then
        arr => volume%point_to_data("SWs")
        ierr = nf90_put_var(self%ncid, config_sws%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put SWs")
        nullify(arr)
    endif

    if (config_swh%output) then
        arr => volume%point_to_data("SWh")
        ierr = nf90_put_var(self%ncid, config_swh%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put SWh")
        nullify(arr)
    endif

    if (config_swv%output) then
        arr => volume%point_to_data("SWv")
        ierr = nf90_put_var(self%ncid, config_swv%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put SWv")
        nullify(arr)
    endif

    if (config_swt%output) then
        arr => volume%point_to_data("SWt")
        ierr = nf90_put_var(self%ncid, config_swt%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put SWt")
        nullify(arr)
    endif

    if (config_snr_zhh%output) then
        arr => volume%point_to_data("SNR_zhh")
        ierr = nf90_put_var(self%ncid, config_SNR_zhh%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put SNR_zhh")
        nullify(arr)
    endif

    if (config_snr_zvv%output) then
        arr => volume%point_to_data("SNR_zvv")
        ierr = nf90_put_var(self%ncid, config_SNR_zvv%varid, arr(1:volume%ngates,1:volume%nrays))
        call error_handler(ierr, "Problem put SNR_zvv")
        nullify(arr)
    endif
    
    ptr => root
    do while ( associated ( ptr%next ) ) 
        ptr => ptr%next

        if ( ptr%source == "COMPUTE" ) then
            ! skip
        else
            if (ptr%varid > -99998) then
                ierr = nf90_put_var(self%ncid, ptr%varid, ptr%data)
                call error_handler(ierr, "Problem put field "//trim(ptr%field_name) )
            endif
        endif
    enddo

    ! print*, 'END   write volume'
  end subroutine write_volume

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine cfradial_close(self)
    use netcdf, only : nf90_close
    implicit none
    class (cfradial_type), intent(inout) :: self
    integer :: ierr
    ierr = nf90_close(self%ncid)
    call error_handler(ierr, "Problem closing file for output")
  end subroutine cfradial_close

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine cfradial_destruct(self)
    implicit none
    class (cfradial_type), intent(inout) :: self
    write(*,'("CFRADIAL_DESTRUCT: address: ",I20)') loc(self)
  end subroutine cfradial_destruct

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  function volume_initialize ( number, reference_time, current_time_seconds, &
       &                       scan_max_range_in_meters, scan_meters_between_gates, &
       &                       scan_fold_limit_lower, scan_fold_limit_upper, &
       &                       ncid1, list1 )
    use kwm_date_utilities, only : geth_newdate
    implicit none
    type(volume_type) :: volume_initialize
    integer,          intent(in) :: number
    character(len=*), intent(in) :: reference_time
    integer,          intent(in) :: current_time_seconds
    real(kind=RKIND),             intent(in) :: scan_max_range_in_meters
    real(kind=RKIND),             intent(in) :: scan_meters_between_gates
    real(kind=RKIND),             intent(in) :: scan_fold_limit_lower
    real(kind=RKIND),             intent(in) :: scan_fold_limit_upper

    integer,                             intent(in) :: ncid1
    character(len=256), dimension(:), intent(in) :: list1

    real(kind=RKIND), parameter :: badval = -1.E36
    type(volume_field_type), pointer :: ptr
    integer :: i

    ! write(*,'("VOLUME_INITIALIZE: address: ",I20)') loc(volume_initialize)
    volume_initialize%number = number
    volume_initialize%ngates = int ( ceiling ( scan_max_range_in_meters / scan_meters_between_gates ) )
    volume_initialize%fold_limit_lower = scan_fold_limit_lower
    volume_initialize%fold_limit_upper = scan_fold_limit_upper

    volume_initialize%reference_time = reference_time
    call geth_newdate ( volume_initialize%time_coverage_start , volume_initialize%reference_time, current_time_seconds )
    print*, 'volume_initialize:  time_coverage_start = ', volume_initialize%reference_time, " + ", current_time_seconds, " = ", volume_initialize%time_coverage_start
    volume_initialize%start_time_integer_seconds = current_time_seconds

    volume_initialize%nrays = 0
    volume_initialize%sweep = 0
    volume_initialize%sweep_number = 0
    volume_initialize%sweep_mode = " "
    volume_initialize%fixed_angle = -9999.
    volume_initialize%sweep_start_ray_index = 0
    volume_initialize%sweep_end_ray_index = 0
    volume_initialize%time_coverage_end = " "
    volume_initialize%time = badval
    volume_initialize%range = badval
    volume_initialize%latitude = badval
    volume_initialize%longitude = badval
    volume_initialize%altitude = badval
    volume_initialize%azimuth = badval
    volume_initialize%elevation = badval
    volume_initialize%georefs_applied = 0
    volume_initialize%heading = badval
    volume_initialize%roll = badval
    volume_initialize%pitch = badval
    volume_initialize%drift = badval
    volume_initialize%rotation = badval
    volume_initialize%tilt = badval

    ! And clear our anything that might be in our linked list.
    ptr => root
    do while ( associated ( ptr%next ) ) 
        ptr => ptr%next
        deallocate(ptr%data)
        root%next => ptr%next
        nullify(ptr)
        ptr => root
    enddo

    do i = 1, size(list1)
        ! Exclude fields:
        if ( list1(i) == "T" ) cycle
        if ( list1(i) == "P" ) cycle
        if ( list1(i) == "PB" ) cycle
        if ( list1(i) == "PH" ) cycle
        if ( list1(i) == "PHB" ) cycle
        if ( list1(i) == "CLDFRA" ) cycle
        ! if ( list1(i) == "QNGRAUPEL" ) cycle
        ! if ( list1(i) == "QNRAIN" ) cycle
        ! if ( list1(i) == "QNSNOW" ) cycle
        ! if ( list1(i) == "QNICE" ) cycle
        ! if ( list1(i) == "QGRAUP" ) cycle
        ! if ( list1(i) == "QSNOW" ) cycle
        ! if ( list1(i) == "QICE" ) cycle
        ! if ( list1(i) == "QRAIN" ) cycle
        ! if ( list1(i) == "QCLOUD" ) cycle
        if ( list1(i) == "QVAPOR" ) cycle
        if ( list1(i) == "P_HYD" ) cycle

        call volume_initialize%initialize_new_field ( ncid1  , list1(i) )
    enddo

    ! Data needed for calls to CRSIM:
    call volume_initialize%initialize_new_field ( -9999  , "TEMPERATURE", source="WRF", description = "Temperature",         units = "K" )
    call volume_initialize%initialize_new_field ( -9999  , "RHO_D", source="WRF", description = "rho_d",         units = "kg/m^3" )

    ! call volume%initialize_new_field ( -9999  , "PRESSURE",    source="WRF", description = "Pressure",            units = "Pa" )
    call volume_initialize%initialize_new_field ( -9999  , "HT",      source="WRF", &
         &                             description = "Geopotential Height", units = "m" )

  end function volume_initialize

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine volume_initialize_new_field(self, ncid, name, source, description, units)
    use netcdf, only : NF90_NOERR
    use netcdf, only : NF90_GLOBAL
    use netcdf, only : nf90_inq_varid
    use netcdf, only : nf90_get_att
    implicit none
    class(volume_type), intent(inout) :: self
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    character(len=*), optional, intent(in) :: source
    character(len=*), optional, intent(in) :: description
    character(len=*), optional, intent(in) :: units
    type(volume_field_type), pointer :: ptr
    integer :: local_varid, ierr
    character(len=256) :: title_string
    character(len=256) :: input_file_string

    ptr => root
    do while ( associated ( ptr%next ) ) 
        ptr => ptr%next
        if ( trim(ptr%field_name) == name ) then
            write(*,'("Field ", A, " already associated")') trim(name)
            stop "module_cfradial_output:volume_initialize_new_field"
        endif
    enddo

    allocate(ptr%next)
    ptr => ptr%next
    ptr%field_name = name

    ! Determine what type of file this is coming from.

    ptr%source = ""
    if ( present ( source ) ) then
        ptr%source = source
    else
        ! WRF output has a "TITLE" attribute:
        ierr = nf90_get_att(ncid, NF90_GLOBAL, "TITLE", title_string)
        if (ierr == NF90_NOERR) then
            ptr%source = "WRF"
        endif

        ! CRSIM output has a "WRF_input_file" attribute:
        ierr = nf90_get_att(ncid, NF90_GLOBAL, "WRF_input_file",input_file_string)
        if (ierr == NF90_NOERR) then
            ptr%source = "CRSIM"
        endif
    endif

    if ( .not. present ( source ) ) then
        ierr = nf90_inq_varid(ncid, trim(ptr%field_name), local_varid)
        call error_handler(ierr, "Problem inq varid for variable '"//trim(ptr%field_name)//"'")
    endif

    if ( present ( description ) ) then
        ptr%field_description = description
    else
        ierr = nf90_get_att(ncid, local_varid, "description", ptr%field_description)
        call error_handler(ierr, "Problem get att description")
    endif

    if ( present ( units ) ) then
        ptr%field_units = units
    else
        ierr = nf90_get_att(ncid, local_varid, "units", ptr%field_units)
        call error_handler(ierr, "Problem get att units")
    endif

  end subroutine volume_initialize_new_field

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  type (volume_field_type) function volume_point_to(self, name) result (ptr)
    implicit none
    class(volume_type), intent(inout) :: self
    character(len=*), intent(in) :: name
    pointer ptr
    ptr => root
    do while ( associated ( ptr%next ) ) 
        ptr => ptr%next
        if ( trim(ptr%field_name) == name ) then
            return
        endif
    enddo
    write(*,'("MODULE_CFRADIAL_OUTPUT:VOLUME_POINT_TO:  Field ", A, " not found.")') trim(name)
    stop
  end function volume_point_to

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  function volume_point_to_data(self, name) result (data)
    !
    ! PURPOSE:
    !
    !       Return a pointer to the %data array from the entry of the 
    !       requested name in the <volume_field_type> linked list.
    !
    implicit none
    class(volume_type), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(kind=RKIND), dimension(:,:), pointer :: data
    type(volume_field_type), pointer :: ptr
    ptr => root
    do while ( associated ( ptr%next ) ) 
        ptr => ptr%next
        if ( trim(ptr%field_name) == name ) then
            data => ptr%data
            return
        endif
    enddo

    write(*,'("MODULE_CFRADIAL_OUTPUT:VOLUME_POINT_TO_DATA:  Field ", A, " not found.")') trim(name)
    stop
  end function volume_point_to_data

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  function volume_new_array(self, name) result (data)
    !
    ! PURPOSE:
    !
    !      Create a new entry in our <volume_field_type> linked list,
    !      allocating the %data array for the entry, but not initializing 
    !      values in the array.  Returns a pointer to the %data array.
    !
    implicit none
    class(volume_type), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(kind=RKIND), dimension(:,:), pointer :: data
    type(volume_field_type), pointer :: ptr

    ! 
    !  Scoot to the end of our linked list, checking to make sure that
    !  we don't already have an entry of the requested name already defined.
    !

    ptr => root
    do while ( associated ( ptr%next ) ) 
        ptr => ptr%next
        if ( trim(ptr%field_name) == name ) then
            write(*,'("field ''",A,"'' already exists....  Stop.")') trim(name)
            stop
        endif
    enddo

    !
    !  We're at the end of our linked list.  Add a new entry to the end.
    !

    allocate(ptr%next)
    ptr => ptr%next
    ptr%field_name = name
    ptr%source = "COMPUTE"

    !
    !  Allocate the %data portion of the entry, and return a pointer
    !  to that %data array.
    !

    allocate ( ptr%data ( self%ngates , self%nrays ) )
    data => ptr%data

  end function volume_new_array

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine volume_destruct(self)
    implicit none
    class (volume_type), intent(inout) :: self
    write(*,'("VOLUME_DESTRUCT: address: ",I20)') loc(self)
  end subroutine volume_destruct

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  type (volume_field_type) function volume_next_point ( self , inptr ) result ( ptr )
    !
    ! PURPOSE:
    !
    !     Return a pointer to the next entry in the <volume_field_type> linked list.
    !     If the current pointer is a null pointer, point to the first entry in the list.
    !     If the next entry cannot be found (i.e., we're already at the last entry),
    !     return a null pointer.
    !
    implicit none
    class(volume_type), intent(inout) :: self
    type(volume_field_type), pointer :: inptr
    pointer ptr
    
    !
    ! If the input pointer is not associated, return a null pointer.
    !

    if (.not. associated (inptr) ) then
        ptr => root%next
        return
    endif

    !
    ! If we can see the next entry, return a pointer to that entry,
    ! nullifying the input pointer.
    !

    if ( associated ( inptr%next ) ) then
        ptr => inptr%next
        nullify(inptr)
        return
    endif

    !
    ! If there's no next entry, return a null pointer.
    !
    
    nullify(ptr)

  end function volume_next_point

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine radial_velocity ( self, grid_dx )
    !
    ! PURPOSE:
    !
    !        Compute the (unfolded) radial velocity field of a full volume of sweeps and beams.
    !        This routine implements equation 26 from Lee et al. (1994).
    !
    ! INPUT:
    !
    !        self
    !               The volume_type class, with the following fields already filled in:
    !
    !                   nrays
    !                   ngates
    !                   range
    !                   elevation
    !                   heading
    !                   drift
    !                   pitch
    !                   eastward_velocity
    !                   northward_velocity
    !                   aircraft_xgrid
    !                   aircraft_ygrid
    !                   aircraft_z
    !
    !               Also, the following must already be defined and filled in as entries 
    !               in the <volume_field_type> linked list:
    !
    !                   U    :: the west-east component of the model wind field at the gate location.
    !                   V    :: the south-north component of the model wind field at the gate location.
    !                   W    :: the vertical component of the model wind field at the gate location.
    !                   VX   :: the x-coordinate of the gate location, in grid units
    !                   VY   :: the y-coordinate of the gate location, in grid units
    !                   VZ   :: the z-coordinate of the gate location, in grid units
    !
    !        grid_dx
    !                The grid spacing [meters] of the model grid.
    !
    ! OUTPUT:
    ! 
    !        self
    !               The volume_type class, with the "VELUNF" and "VELTRU" fields in the
    !               <volume_field_type> linked list newly defined and filled in:
    !
    !                  VELUNF  :: the (unfolded) radial velocity at the gate location.
    !
    !                  VELTRU  :: radial velocity not consdering aircraft motion.
    !
    ! SIDE EFFECTS:
    !
    !       The field "VELUNF" is created in the <volume_field_type> linked list, and filled in with the
    !       (unfolded) radial velocity.
    !
    !       The field "VELTRU" is created in the <volume_field_type> linked list, and filled in with the
    !       true (i.e., not considering aircraft motion) radial velocity.
    !
    !       In case of error in finding or creating entries in the <volume_field_type> linked list, 
    !       the program will print an error message and stop.
    !
    
    use module_llxy, only : RAD_PER_DEG
    implicit none
    class (volume_type) :: self
    real(kind=RKIND), intent(in) :: grid_dx
    real(kind=RKIND), pointer, dimension(:,:) :: uptr   ! Pointer to the  "U"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vptr   ! Pointer to the  "V"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: wptr   ! Pointer to the  "W"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vx     ! Pointer to the "VX"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vy     ! Pointer to the "VY"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vz     ! Pointer to the "VZ"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: velunf ! Pointer to the "VELUNF" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: veltru ! Pointer to the "VELTRU" field in the <volume_field_type> linked list

    real(kind=RKIND)    :: x     ! x position of a gate [m]
    real(kind=RKIND)    :: y     ! y position of a gate [m]
    real(kind=RKIND)    :: z     ! z position of a gate [m]
    real(kind=RKIND)    :: r     ! gate range, distance from radar (m)
    real(kind=RKIND)    :: u     ! west-east   component of the model wind field at a gate [m/s]
    real(kind=RKIND)    :: v     ! south-north component of the model wind field at a gate [m/s]
    real(kind=RKIND)    :: w     ! vertical    component of the model wind field at a gate [m/s]
    real(kind=RKIND)    :: phi   ! beam elevation angle [radians]
    real(kind=RKIND)    :: VG    ! aircraft horizontal ground speed [m/s]
    real(kind=RKIND)    :: WG    ! aircraft vertical ground speed [m/s]
    real(kind=RKIND)    :: H     ! aircraft heading [radians]
    real(kind=RKIND)    :: T     ! aircraft track [radians] = heading + drift
    real(kind=RKIND)    :: P     ! aircraft pitch [radians]
    real(kind=RKIND)    :: L     ! Distance [m] between aircraft INS and radar.  Only longitudinal distance considered?
    real(kind=RKIND)    :: dHdt  ! Rate of change of heading [radians/s]
    real(kind=RKIND)    :: dPdt  ! Rate of change of pitch [radias/s]
    real(kind=RKIND)    :: vt    ! Terminal fallspeed of the scatterer [m/s; positive downwards]
    integer :: iray  ! Loop control variable, looping over the beams in a volume
    integer :: igate ! loop control variable, looping over the gates in a beam
    integer :: k

    uptr => self%point_to_data("U")
    vptr => self%point_to_data("V")
    wptr => self%point_to_data("W")
    vx   => self%point_to_data("VX")
    vy   => self%point_to_data("VY")
    vz   => self%point_to_data("VZ")

    velunf => self%new_array("VELUNF")
    velunf = -9999 ! Default (i.e., no-data) value which will get overwritten below if we have the data.

    veltru => self%new_array("VELTRU")
    veltru = -9999 ! Default (i.e., no-data) value which will get overwritten below if we have the data.

    do iray = 1, self%nrays

        phi = self%elevation(iray) * RAD_PER_DEG
        T   = ( self%heading(iray) + self%drift(iray) ) * RAD_PER_DEG
        P   = self%pitch(iray) * RAD_PER_DEG
        H   = self%heading(iray) * RAD_PER_DEG
        VG = sqrt(self%eastward_velocity(iray)**2 + self%northward_velocity(iray)**2)
        WG  = self%vertical_velocity(iray)
        L   = 0.0 ! Maybe an aircraft parameter.
        dHdt = 0.0 ! until we compute it for the aircraft.
        dPdt = 0.0 ! until we compute it for the aircraft.

        GATE_LOOP : do igate = 1, self%ngates

            if ( uptr(igate,iray) < -999 ) cycle GATE_LOOP
            if ( vptr(igate,iray) < -999 ) cycle GATE_LOOP
            if ( wptr(igate,iray) < -999 ) cycle GATE_LOOP

            u = uptr(igate,iray)
            v = vptr(igate,iray)
            w = wptr(igate,iray)

            ! Compute gate location (relative to aircraft) as an offset from the aircraft in meters.
            r = self%range(igate)
            x = grid_dx * ( vx(igate,iray) - self%aircraft_xgrid(iray))
            y = grid_dx * ( vy(igate,iray) - self%aircraft_ygrid(iray))
            z =             vz(igate,iray) - self%aircraft_z(iray)
            vt = 0.0 ! Until we start considering this.

            velunf(igate,iray) = &
                 &  (1./r)*(u*x+v*y)+(w-vt-WG)*sin(phi)-(VG/r)*(x*sin(T)+y*cos(T))  &
                 &  + (L/r)*(x*((1.+cos(P))*cos(H)*(dHdt)-sin(P)*sin(H)*(dPdt))  &
                 &         - y*((1.+cos(P))*sin(H)*(dHdt)+sin(P)*cos(H)*(dPdt))  &
                 &         + z*cos(P)*dPdt)
            
            if (abs(velunf(igate,iray)) > 500) then
                print*, 'iray = ', iray
                print*, 'eastward_velocity(iray) = ', self%eastward_velocity(iray)
                print*, 'northward_velocity(iray) = ', self%northward_velocity(iray)
                print*, 'igate = ', igate
                print*, 'u = ', u
                print*, 'v = ', v
                print*, 'w = ', w
                print*, 'r = ', r
                print*, 'x = ', x
                print*, 'y = ', y
                print*, 'z = ', z
                print*, 'T = ', T
                print*, 'VG = ', VG
                print*, 'WG = ', WG
                print*, 'vt = ', vt
                stop "VELUNF"
            endif

            veltru(igate,iray) = &
                 &  (1./r)*(u*x+v*y)+(w-vt-0.0)*sin(phi)-(0.0/r)*(x*sin(T)+y*cos(T))  &
                 &  + (L/r)*(x*((1.+cos(P))*cos(H)*(dHdt)-sin(P)*sin(H)*(dPdt))  &
                 &         - y*((1.+cos(P))*sin(H)*(dHdt)+sin(P)*cos(H)*(dPdt))  &
                 &         + z*cos(P)*dPdt)
                 

        enddo GATE_LOOP

    enddo

  end subroutine radial_velocity

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine radial_velocity_fvel ( self, grid_dx )
    ! Added by B. Klotz (10/7/2021)
    ! PURPOSE:
    !
    !        Compute the (unfolded) radial velocity field of a full volume of sweeps and beams.
    !        This routine implements equation 26 from Lee et al. (1994). And it includes
    !        The average fall velocity within each gate from the crsim module
    !
    ! INPUT:
    !
    !        self
    !               The volume_type class, with the following fields already filled in:
    !
    !                   nrays
    !                   ngates
    !                   range
    !                   elevation
    !                   heading
    !                   drift
    !                   pitch
    !                   eastward_velocity
    !                   northward_velocity
    !                   aircraft_xgrid
    !                   aircraft_ygrid
    !                   aircraft_z
    !
    !               Also, the following must already be defined and filled in as entries
    !               in the <volume_field_type> linked list:
    !
    !                   U    :: the west-east component of the model wind field at the gate location.
    !                   V    :: the south-north component of the model wind field at the gate location.
    !                   W    :: the vertical component of the model wind field at the gate location.
    !                   VX   :: the x-coordinate of the gate location, in grid units
    !                   VY   :: the y-coordinate of the gate location, in grid units
    !                   VZ   :: the z-coordinate of the gate location, in grid units
    !
    !        grid_dx
    !                The grid spacing [meters] of the model grid.
    !
    ! OUTPUT:
    !
    !        self
    !               The volume_type class, with the "VELUNF" and "VELTRU" fields in the
    !               <volume_field_type> linked list newly defined and filled in:
    !
    !                  VELUNF  :: the (unfolded) radial velocity at the gate location.
    !
    !                  VELTRU  :: radial velocity not consdering aircraft motion.
    !
    ! SIDE EFFECTS:
    !
    !       The field "VELUNF" is created in the <volume_field_type> linked list, and filled in with the
    !       (unfolded) radial velocity.
    !
    !       The field "VELTRU" is created in the <volume_field_type> linked list, and filled in with the
    !       true (i.e., not considering aircraft motion) radial velocity.
    !
    !       In case of error in finding or creating entries in the <volume_field_type> linked list,
    !       the program will print an error message and stop.
    !

    ! B. Klotz (10/8/2021)
    ! Added variables for inclusion of fvel
    
    use module_llxy, only : RAD_PER_DEG
    implicit none
    class (volume_type) :: self
    real(kind=RKIND), intent(in) :: grid_dx
    real(kind=RKIND), pointer, dimension(:,:) :: uptr   ! Pointer to the  "U"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vptr   ! Pointer to the  "V"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: wptr   ! Pointer to the  "W"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vx     ! Pointer to the "VX"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vy     ! Pointer to the "VY"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: vz     ! Pointer to the "VZ"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: fvel     ! Pointer to the "fvel"     field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: velunf_fvel ! Pointer to the "VELUNF" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: veltru_fvel ! Pointer to the "VELTRU" field in the <volume_field_type> linked list

    real(kind=RKIND)    :: x     ! x position of a gate [m]
    real(kind=RKIND)    :: y     ! y position of a gate [m]
    real(kind=RKIND)    :: z     ! z position of a gate [m]
    real(kind=RKIND)    :: r     ! gate range, distance from radar (m)
    real(kind=RKIND)    :: u     ! west-east   component of the model wind field at a gate [m/s]
    real(kind=RKIND)    :: v     ! south-north component of the model wind field at a gate [m/s]
    real(kind=RKIND)    :: w     ! vertical    component of the model wind field at a gate [m/s]
    real(kind=RKIND)    :: phi   ! beam elevation angle [radians]
    real(kind=RKIND)    :: VG    ! aircraft horizontal ground speed [m/s]
    real(kind=RKIND)    :: WG    ! aircraft vertical ground speed [m/s]
    real(kind=RKIND)    :: H     ! aircraft heading [radians]
    real(kind=RKIND)    :: T     ! aircraft track [radians] = heading + drift
    real(kind=RKIND)    :: P     ! aircraft pitch [radians]
    real(kind=RKIND)    :: L     ! Distance [m] between aircraft INS and radar.  Only longitudinal distance considered?
    real(kind=RKIND)    :: dHdt  ! Rate of change of heading [radians/s]
    real(kind=RKIND)    :: dPdt  ! Rate of change of pitch [radias/s]
    real(kind=RKIND)    :: vt    ! Terminal fallspeed of the scatterer [m/s; positive downwards]
    integer :: iray  ! Loop control variable, looping over the beams in a volume
    integer :: igate ! loop control variable, looping over the gates in a beam
    integer :: k

    uptr => self%point_to_data("U")
    vptr => self%point_to_data("V")
    wptr => self%point_to_data("W")
    vx   => self%point_to_data("VX")
    vy   => self%point_to_data("VY")
    vz   => self%point_to_data("VZ")
    fvel => self%point_to_data("fvel")

    velunf_fvel => self%new_array("VELUNF_FVEL")
    velunf_fvel = -9999 ! Default (i.e., no-data) value which will get overwritten below if we have the data.

    veltru_fvel => self%new_array("VELTRU_FVEL")
    veltru_fvel = -9999 ! Default (i.e., no-data) value which will get overwritten below if we have the data.

    do iray = 1, self%nrays

        phi = self%elevation(iray) * RAD_PER_DEG
        T   = ( self%heading(iray) + self%drift(iray) ) * RAD_PER_DEG
        P   = self%pitch(iray) * RAD_PER_DEG
        H   = self%heading(iray) * RAD_PER_DEG
        VG = sqrt(self%eastward_velocity(iray)**2 + self%northward_velocity(iray)**2)
        WG  = self%vertical_velocity(iray)
        L   = 0.0 ! Maybe an aircraft parameter.
        dHdt = 0.0 ! until we compute it for the aircraft.
        dPdt = 0.0 ! until we compute it for the aircraft.

        GATE_LOOP : do igate = 1, self%ngates

            if ( uptr(igate,iray) < -999 ) cycle GATE_LOOP
            if ( vptr(igate,iray) < -999 ) cycle GATE_LOOP
            if ( wptr(igate,iray) < -999 ) cycle GATE_LOOP

            u = uptr(igate,iray)
            v = vptr(igate,iray)
            w = wptr(igate,iray)
            vt = fvel(igate,iray)

            ! Compute gate location (relative to aircraft) as an offset from the aircraft in meters.
            r = self%range(igate)
            x = grid_dx * ( vx(igate,iray) - self%aircraft_xgrid(iray))
            y = grid_dx * ( vy(igate,iray) - self%aircraft_ygrid(iray))
            z =             vz(igate,iray) - self%aircraft_z(iray)
            !vt = 0.0 ! Until we start considering this.

            if (vt >= 0) then
               velunf_fvel(igate,iray) = &
                 &  (1./r)*(u*x+v*y)+(w-vt-WG)*sin(phi)-(VG/r)*(x*sin(T)+y*cos(T))  &
                 &  + (L/r)*(x*((1.+cos(P))*cos(H)*(dHdt)-sin(P)*sin(H)*(dPdt))  &
                 &         - y*((1.+cos(P))*sin(H)*(dHdt)+sin(P)*cos(H)*(dPdt))  &
                 &         + z*cos(P)*dPdt)
            else
               velunf_fvel(igate,iray) = -9999.0
            endif

            !print *, 'ray, gate, FVEL, VELUNF_FVEL: ',iray,igate,vt,velunf_fvel(igate,iray)
            !if (abs(velunf_fvel(igate,iray)) > 500) then
            !    print*, 'iray = ', iray
            !!    print*, 'eastward_velocity(iray) = ', self%eastward_velocity(iray)
            !    print*, 'northward_velocity(iray) = ', self%northward_velocity(iray)
            !    print*, 'igate = ', igate
            !    print*, 'u = ', u
            !    print*, 'v = ', v
            !    print*, 'w = ', w
            !    print*, 'r = ', r
            !    print*, 'x = ', x
            !    print*, 'y = ', y
            !    print*, 'z = ', z
            !    print*, 'T = ', T
            !    print*, 'VG = ', VG
            !    print*, 'WG = ', WG
            !    print*, 'vt = ', vt
            !    stop "VELUNF_FVEL"
            !endif

            if (vt >= 0) then
               veltru_fvel(igate,iray) = &
                 &  (1./r)*(u*x+v*y)+(w-vt-0.0)*sin(phi)-(0.0/r)*(x*sin(T)+y*cos(T))  &
                 &  + (L/r)*(x*((1.+cos(P))*cos(H)*(dHdt)-sin(P)*sin(H)*(dPdt))  &
                 &         - y*((1.+cos(P))*sin(H)*(dHdt)+sin(P)*cos(H)*(dPdt))  &
                 &         + z*cos(P)*dPdt)
            else
                veltru_fvel(igate,iray) = -9999.0
            endif

        enddo GATE_LOOP

    enddo

  end subroutine radial_velocity_fvel

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine fold_radial_velocity(self, snr_mask_threshold)
    implicit none
    class (volume_type) :: self
    real(kind=RKIND), intent(In) :: snr_mask_threshold

    ! B. Klotz (10/8/2021)
    ! Added variables for inclusion of fvel

    ! 
    ! Velocity folding
    !

    real(kind=RKIND), pointer, dimension(:,:) :: vel    ! Pointer to the "VEL" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: velunf ! Pointer to the "VELUNF" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: velunf_fvel
    real(kind=RKIND), pointer, dimension(:,:) :: veltru
    real(kind=RKIND), pointer, dimension(:,:) :: veltru_fvel  
    real(kind=RKIND), pointer, dimension(:,:) :: veltruf
    real(kind=RKIND), pointer, dimension(:,:) :: veltruf_fvel  
    real(kind=RKIND), pointer, dimension(:,:) :: stdv_vel
    real(kind=RKIND), pointer, dimension(:,:) :: snr_zhh
    real(kind=RKIND), pointer, dimension(:,:) :: vel_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: velunf_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: velunf_fvel_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: veltru_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: veltru_fvel_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: veltruf_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: veltruf_fvel_noise_added
    integer :: iray
    integer :: igate
    integer :: k
    real(kind=RKIND) :: rn_draw

    if ( self%fold_limit_lower < -1.E4 ) return
    if ( self%fold_limit_upper >  1.E4 ) return

    velunf => self%point_to_data("VELUNF")
    veltru => self%point_to_data("VELTRU")
    velunf_fvel => self%point_to_data("VELUNF_FVEL")
    veltru_fvel => self%point_to_data("VELTRU_FVEL")
    snr_zhh    => self%point_to_data("SNR_zhh")
    
    vel    => self%new_array("VEL")
    vel = -9999. ! Default (i.e., no-data) value which will get overwritten below if we have the data
    
    velunf_noise_added => self%new_array("VELUNF_noise_added")
    velunf_noise_added = -9999.
      
    velunf_fvel_noise_added => self%new_array("VELUNF_FVEL_noise_added")
    velunf_fvel_noise_added = -9999.
      
    veltru_noise_added => self%new_array("VELTRU_noise_added")
    veltru_noise_added = -9999.

    veltru_fvel_noise_added => self%new_array("VELTRU_FVEL_noise_added")
    veltru_fvel_noise_added = -9999.

    veltruf_noise_added => self%new_array("VELTRUF_noise_added")
    veltruf_noise_added = -9999.

    veltruf_fvel_noise_added => self%new_array("VELTRUF_FVEL_noise_added")
    veltruf_fvel_noise_added = -9999.
      
    vel_noise_added => self%new_array("VEL_noise_added")
    vel_noise_added = -9999.

    ! Apply noise
    stdv_vel => self%point_to_data("STDV_VEL")
    do iray = 1, self%nrays
        do igate = 1, self%ngates
            if (velunf(igate,iray) > -9999) then
                if ( (stdv_vel(igate,iray) > -9999) .and. (SNR_zhh(igate,iray) > snr_mask_threshold) ) then
                    rn_draw = random_normal()
                    velunf_noise_added(igate,iray) = velunf(igate,iray) + stdv_vel(igate,iray) * rn_draw
                    veltru_noise_added(igate,iray) = veltru(igate,iray) + stdv_vel(igate,iray) * rn_draw
                    velunf_fvel_noise_added(igate,iray) = velunf_fvel(igate,iray) + stdv_vel(igate,iray) * rn_draw
                    veltru_fvel_noise_added(igate,iray) = veltru_fvel(igate,iray) + stdv_vel(igate,iray) * rn_draw
                endif
            endif
        enddo
    enddo

    veltruf => self%new_array("VELTRUF")
    veltruf = -9999.

    veltruf_fvel => self%new_array("VELTRUF_FVEL")
    veltruf_fvel = -9999.

    RAY_LOOP : do iray = 1, self%nrays
        self%nyquist_velocity(iray)   = self%fold_limit_upper
        ! self%unambiguous_range(iray) = 71379.16
        GATE_LOOP : do igate = 1, self%ngates

            if ( abs(velunf_noise_added(igate,iray)) > 9998 ) cycle GATE_LOOP

            vel(igate,iray) = velunf(igate,iray)
            vel_noise_added(igate,iray) = velunf_noise_added(igate,iray)! velunf(igate,iray)

            do k = 1, 10
                if ( vel(igate,iray) < self%fold_limit_upper ) exit
                vel(igate,iray) = self%fold_limit_lower + (vel(igate,iray)-self%fold_limit_upper)
                if (k==10) then
                    print*, 'velunf = ', velunf(igate,iray)
                    stop "A) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( vel(igate,iray) > self%fold_limit_lower ) exit
                vel(igate,iray) = self%fold_limit_upper + (vel(igate,iray)-self%fold_limit_lower)
                if (k==10) then
                    print*, 'velunf = ', velunf(igate,iray)
                    stop "B) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( vel_noise_added(igate,iray) < self%fold_limit_upper ) exit
                vel_noise_added(igate,iray) = self%fold_limit_lower + (vel_noise_added(igate,iray)-self%fold_limit_upper)
                if (k==10) then
                    print*, 'velunf_noise_added = ', velunf_noise_added(igate,iray)
                    stop "C) Fold?"
                endif
            enddo
            
            do k = 1, 10
                if ( vel_noise_added(igate,iray) > self%fold_limit_lower ) exit
                vel_noise_added(igate,iray) = self%fold_limit_upper + (vel_noise_added(igate,iray)-self%fold_limit_lower)
                if (k==10) then
                    print*, 'velunf_noise_added = ', velunf_noise_added(igate,iray)
                    stop "D) Fold?"
                endif
            enddo
            
            ! And do the folding for the "true" vel field and the one includeing fvel.
            veltruf(igate,iray) = veltru(igate,iray)
            veltruf_noise_added(igate,iray) = veltru_noise_added(igate,iray)

            veltruf_fvel(igate,iray) = veltru_fvel(igate,iray)
            veltruf_fvel_noise_added(igate,iray) = veltru_noise_added(igate,iray)

            do k = 1, 10
                if ( veltruf(igate,iray) < self%fold_limit_upper ) exit
                veltruf(igate,iray) = self%fold_limit_lower + (veltruf(igate,iray)-self%fold_limit_upper)
                if (k==10) then
                    print*, 'veltru = ', veltru(igate,iray)
                    stop "C) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( veltruf(igate,iray) > self%fold_limit_lower ) exit
                veltruf(igate,iray) = self%fold_limit_upper + (veltruf(igate,iray)-self%fold_limit_lower)
                if (k==10) then
                    print*, 'veltru = ', veltru(igate,iray)
                    stop "D) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( veltruf_noise_added(igate,iray) < self%fold_limit_upper ) exit
                veltruf_noise_added(igate,iray) = self%fold_limit_lower + (veltruf_noise_added(igate,iray)-self%fold_limit_upper)
                if (k==10) then
                    print*, 'veltru_noise_added = ', veltru_noise_added(igate,iray)
                    stop "E) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( veltruf_noise_added(igate,iray) > self%fold_limit_lower ) exit
                veltruf_noise_added(igate,iray) = self%fold_limit_upper + (veltruf_noise_added(igate,iray)-self%fold_limit_lower)
                if (k==10) then
                    print*, 'veltru_noise_added = ', veltru_noise_added(igate,iray)
                    stop "F) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( veltruf_fvel(igate,iray) < self%fold_limit_upper ) exit
                veltruf_fvel(igate,iray) = self%fold_limit_lower + (veltruf_fvel(igate,iray)-self%fold_limit_upper)
                if (k==10) then
                    print*, 'veltru_fvel = ', veltru_fvel(igate,iray)
                    stop "G) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( veltruf_fvel(igate,iray) > self%fold_limit_lower ) exit
                veltruf_fvel(igate,iray) = self%fold_limit_upper + (veltruf_fvel(igate,iray)-self%fold_limit_lower)
                if (k==10) then
                    print*, 'veltru_fvel = ', veltru_fvel(igate,iray)
                    stop "H) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( veltruf_fvel_noise_added(igate,iray) < self%fold_limit_upper ) exit
                veltruf_fvel_noise_added(igate,iray) = self%fold_limit_lower + (veltruf_fvel_noise_added(igate,iray)-self%fold_limit_upper)
                if (k==10) then
                    print*, 'veltru_fvel_noise_added = ', veltru_fvel_noise_added(igate,iray)
                    stop "I) Fold?"
                endif
            enddo

            do k = 1, 10
                if ( veltruf_fvel_noise_added(igate,iray) > self%fold_limit_lower ) exit
                veltruf_fvel_noise_added(igate,iray) = self%fold_limit_upper + (veltruf_fvel_noise_added(igate,iray)-self%fold_limit_lower)
                if (k==10) then
                    print*, 'veltru_fvel_noise_added = ', veltru_fvel_noise_added(igate,iray)
                    stop "J) Fold?"
                endif
            enddo
            
        enddo GATE_LOOP
    enddo RAY_LOOP

  end subroutine fold_radial_velocity

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine attenuate(self, meters_between_gates)
    implicit none
    class (volume_type) :: self
    real(kind=RKIND), intent(in) :: meters_between_gates

    !
    ! Apply crsim-computed attenuation to reflectivity fields.
    !

    integer :: iray
    integer :: igate
    real(kind=RKIND)    :: sumattH
    real(kind=RKIND)    :: sumattV
    real(kind=RKIND)    :: sumattD

    real(kind=RKIND), pointer, dimension(:,:) :: zhh    ! Pointer to the "Zhh" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: ah     ! Pointer to the "Ah" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: zhh_attenuated  ! Pointer to new field in the <volume_field_type> linked list

    real(kind=RKIND), pointer, dimension(:,:) :: zvv    ! Pointer to the "Zvv" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: av     ! Pointer to the "Av" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: zvv_attenuated  ! Pointer to new field in the <volume_field_type> linked list

    real(kind=RKIND), pointer, dimension(:,:) :: zvh    ! Pointer to the "Zvh" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: zvh_attenuated  ! Pointer to new field in the <volume_field_type> linked list

    real(kind=RKIND), pointer, dimension(:,:) :: zdr    ! Pointer to the "Zdr" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: adp    ! Pointer to the "Adp" field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: zdr_attenuated   ! Pointer to new field in the <volume_field_type> linked list


    zhh => self%point_to_data("Zhh")
    ah  => self%point_to_data("Ah")
    zhh_attenuated    => self%new_array("Zhh_attenuated")
    zhh_attenuated = -9999.0

    zvv => self%point_to_data("Zvv")
    av  => self%point_to_data("Av")
    zvv_attenuated    => self%new_array("Zvv_attenuated")
    zvv_attenuated = -9999.0

    zvh => self%point_to_data("Zvh")
    zvh_attenuated    => self%new_array("Zvh_attenuated")
    zvh_attenuated = -9999.0

    zdr => self%point_to_data("Zdr")
    adp => self%point_to_data("Av")
    zdr_attenuated    => self%new_array("Zdr_attenuated")
    zdr_attenuated = -9999.0

    RAY_LOOP : do iray = 1, self%nrays

        !  For each ray, integrate the attenuation term (already interpolated to the rays/gates), 
        !  along the ray.
        
        ! Assuming 1-way attenuation, so it's doubled for Zhh and Zvv, 
        ! and Ah and Av both applied separately for Zvh.

        sumattH = 0.0
        sumattV = 0.0
        sumattD = 0.0
        GATE_LOOP : do igate = 1, self%ngates
            if (zhh(igate,iray) > -999) then
                zhh_attenuated(igate,iray) = zhh(igate,iray) - 2*sumattH
            endif
            if (zvv(igate,iray) > -999) then
                zvv_attenuated(igate,iray) = zvv(igate,iray) - 2*sumattV
            endif
            if (zvh(igate,iray) > -999) then
                zvh_attenuated(igate,iray) = zvh(igate,iray) - sumattV - sumattH
            endif
            if (zdr(igate,iray) > -999) then
                zdr_attenuated(igate,iray) = zdr(igate,iray) - 2*sumattD
            endif
            if (Ah(igate,iray) > -999) then
                sumattH = sumattH + (Ah(igate,iray)*1.E-3*meters_between_gates)
            endif
            if (Av(igate,iray) > -999) then
                sumattV = sumattV + (Av(igate,iray)*1.E-3*meters_between_gates)
            endif
            if (Adp(igate,iray) > -999) then
                sumattD = sumattD + (Adp(igate,iray)*1.E-3*meters_between_gates)
            endif
        enddo GATE_LOOP
    enddo RAY_LOOP
    
  end subroutine attenuate

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine compute_snr(self)
    ! Signal-to-Noise Ratio
    implicit none
    class (volume_type) :: self
    integer :: iray
    integer :: igate
    real(kind=RKIND), pointer, dimension(:,:) :: SNR_zhh  ! Pointer to new field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: SNR_zvv  ! Pointer to new field in the <volume_field_type> linked list
    real(kind=RKIND), pointer, dimension(:,:) :: Zhh_attenuated
    real(kind=RKIND), pointer, dimension(:,:) :: Zvv_attenuated
    real(kind=RKIND), parameter :: MDS1KM = -30.0 ! Minimum detectable signal (dBZ) at 1km

    SNR_zhh => self%new_array("SNR_zhh")
    SNR_zvv => self%new_array("SNR_zvv")
    
    Zhh_attenuated => self%point_to_data("Zhh_attenuated")
    Zvv_attenuated => self%point_to_data("Zvv_attenuated")

    do iray = 1, self%nrays
        do igate = 1, self%ngates
            if (Zhh_attenuated(igate,iray) > -9999) then
                SNR_zhh(igate,iray) = Zhh_attenuated(igate,iray) - ( MDS1KM + 20.0*log10( 1.E-3 * self%range(igate) ) )
            else
                SNR_zhh(igate,iray) = -9999.0
            endif
            if (Zvv_attenuated(igate,iray) > -9999) then
                SNR_zvv(igate,iray) = Zvv_attenuated(igate,iray) - ( MDS1KM + 20.0*log10( 1.E-3 * self%range(igate) ) )
            else
                SNR_zvv(igate,iray) = -9999.0
            endif
        enddo
    enddo

  end subroutine compute_snr

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine compute_ips(self, revisits_per_acquisition_time, pulse_repetition_time, radar_frequency_gigahertz)
    !
    ! Independent Pulse Sampling
    !
    implicit none
    class (volume_type)                       :: self
    integer,                       intent(in) :: revisits_per_acquisition_time
    real(kind=RKIND),              intent(in) :: pulse_repetition_time
    real(kind=RKIND),              intent(in) :: radar_frequency_gigahertz
    real(kind=RKIND), pointer, dimension(:,:) :: SWtot
    real(kind=RKIND), pointer, dimension(:,:) :: SNR_zhh
    real(kind=RKIND), pointer, dimension(:,:) :: SNR_zvv
    real(kind=RKIND), pointer, dimension(:,:) :: stdv_zhh
    real(kind=RKIND), pointer, dimension(:,:) :: stdv_zvv
    real(kind=RKIND), pointer, dimension(:,:) :: stdv_vel
    real(kind=RKIND)                          :: rho
    real(kind=RKIND)                          :: lambda
    integer                                   :: iray
    integer                                   :: igate
    integer                                   :: L
    real(kind=RKIND)                          :: K
    real(kind=RKIND)                          :: YU3_zhh
    real(kind=RKIND)                          :: YU3_zvv
    real(kind=RKIND)                          :: YU4_zhh

    real(kind=RKIND)                          :: SNR_zhh_lin
    real(kind=RKIND)                          :: SNR_zvv_lin
    
    integer         , parameter               :: m = 1
    real(kind=RKIND), parameter               :: pi = 4.0*atan(1.0_RKIND)
    

    lambda = 2.9979E8 / (radar_frequency_gigahertz * 1.E9)
    
    SWtot => self%point_to_data("SWtot")
    SNR_zhh => self%point_to_data("SNR_zhh")
    SNR_zvv => self%point_to_data("SNR_zvv")
    stdv_zhh => self%new_array("STDV_Zhh")
    stdv_zvv => self%new_array("STDV_Zvv")
    stdv_vel => self%new_array("STDV_VEL")
    stdv_zhh = -9999.
    stdv_zvv = -9999.
    stdv_vel = -9999.
    L = revisits_per_acquisition_time

    do iray = 1, self%nrays
        do igate = 1, self%ngates
            if (SWtot(igate,iray) > -9998) then
                
                !
                !  The SNR array as computed earlier is in dB. Eq. (3) from Yu2007
                !  needs linear units:
                !
    
                SNR_zhh_lin = 10**(SNR_zhh(igate,iray)*0.1)
                SNR_zvv_lin = 10**(SNR_zvv(igate,iray)*0.1)

                !
                ! OK, Scott gave me this equation for rho but where does it come from?
                !
                
                rho = exp(-8.0*(pi*SWtot(igate,iray)*m*pulse_repetition_time/lambda)**2)

                !
                ! Equations (3) and (4) from Yu et al. 2007
                !
                
                K = (lambda*lambda)/(32*(pi*pulse_repetition_time*rho)**2)
                
                YU3_zhh = 1.0/(2.0*L) * ( (1.0+rho**2) + (1.0/SNR_zhh_lin)**2 + (2.0/SNR_zhh_lin) )
                YU3_zvv = 1.0/(2.0*L) * ( (1.0+rho**2) + (1.0/SNR_zvv_lin)**2 + (2.0/SNR_zvv_lin) )

                YU4_zhh = (K/real(L,kind=RKIND)) * ( (1.0-rho**2) + (1.0/SNR_zhh_lin)**2 + (2.0/SNR_zhh_lin) )

                !
                !  And compute a standard deviation from YU3
                !
                
                stdv_Zhh(igate,iray) = 10*log10(1.0+sqrt(YU3_zhh))
                stdv_Zvv(igate,iray) = 10*log10(1.0+sqrt(YU3_zvv))
                stdv_Vel(igate,iray) = sqrt(YU4_zhh)
                
            endif
        enddo
    enddo

  end subroutine compute_ips

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

  subroutine add_zhh_noise(self, snr_mask_threshold)
    implicit none
    class (volume_type) :: self
    real(kind=RKIND) :: snr_mask_threshold
    real(kind=RKIND), pointer, dimension(:,:) :: Zhh_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: STDV_Zhh
    real(kind=RKIND), pointer, dimension(:,:) :: SNR_Zhh
    real(kind=RKIND), pointer, dimension(:,:) :: Zhh_attenuated
    real(kind=RKIND), pointer, dimension(:,:) :: Zvv_noise_added
    real(kind=RKIND), pointer, dimension(:,:) :: STDV_Zvv
    real(kind=RKIND), pointer, dimension(:,:) :: SNR_Zvv
    real(kind=RKIND), pointer, dimension(:,:) :: Zvv_attenuated
    integer :: igate
    integer :: iray
    
    Zhh_noise_added => self%new_array("Zhh_noise_added")
    Zhh_attenuated => self%point_to_data("Zhh_attenuated")
    STDV_Zhh => self%point_to_data("STDV_Zhh")
    SNR_Zhh => self%point_to_data("SNR_zhh")

    Zvv_noise_added => self%new_array("Zvv_noise_added")
    Zvv_attenuated => self%point_to_data("Zvv_attenuated")
    STDV_Zvv => self%point_to_data("STDV_Zvv")
    SNR_Zvv => self%point_to_data("SNR_zvv")
    
    do iray = 1, self%nrays
        do igate = 1, self%ngates
            if ((STDV_Zhh(igate,iray) > -9998).and.(SNR_zhh(igate,iray) > snr_mask_threshold)) then
                Zhh_noise_added(igate,iray) = Zhh_attenuated(igate,iray) + STDV_zhh(igate,iray) * random_normal()
            else
                Zhh_noise_added(igate,iray) = -9999.
            endif
            if ((STDV_Zvv(igate,iray) > -9998).and.(SNR_zvv(igate,iray) > snr_mask_threshold)) then
                Zvv_noise_added(igate,iray) = Zvv_attenuated(igate,iray) + STDV_zvv(igate,iray) * random_normal()
            else
                Zvv_noise_added(igate,iray) = -9999.
            endif
        enddo
    enddo
    
  end subroutine add_zhh_noise

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

 subroutine compute_phidp(self, snr_mask_threshold)
  implicit none
  class (volume_type) :: self
  real(kind=RKIND), intent(In) :: snr_mask_threshold

  integer :: iray
  integer :: igate
  integer :: id1, ii
  real(kind=RKIND) :: tmpsum = 0.0
  real(kind=RKIND) :: xarr(3) = (/ 0.0, 0.0, 0.0 /)
  real(kind=RKIND) :: tarr(3) = (/ 0.0, 0.0, 0.0 /)
  real(kind=RKIND) :: rng_diff
  real(kind=RKIND) :: missing_val
  real(kind=RKIND), pointer, dimension(:,:) :: zhh    ! Pointer to the "Zhh" field in the <volume_field_type> linked list
  real(kind=RKIND), pointer, dimension(:,:) :: kdp    ! Pointer to the "Kdp" field in the <volume_field_type> linked list
  real(kind=RKIND), pointer, dimension(:,:) :: dphase   ! Pointer to the "DPhase" field in the <volume_field_type> linked list
  real(kind=RKIND), pointer, dimension(:,:) :: phidp ! Pointer to new field in the <volume_field_type> linked list
  real(kind=RKIND), pointer, dimension(:,:) :: phidp_wdp ! Pointer to new field in the <volume_field_type> linked list 
  real(kind=RKIND), pointer, dimension(:,:) :: SNR_Zhh
  real(kind=RKIND), pointer, dimension(:,:) :: STDV_Zhh
  

  kdp => self%point_to_data("Kdp")
  dphase => self%point_to_data("DPhase")
  SNR_Zhh => self%point_to_data("SNR_zhh")
  STDV_Zhh => self%point_to_data("STDV_Zhh")
  zhh => self%point_to_data("Zhh")
  phidp => self%new_array("PHIdp")
  phidp_wdp  => self%new_array("PHIdp_DPh")
  phidp = -9999.0
  phidp_wdp = -9999.0
  missing_val = -9999.0    

  RAY_LOOP : do iray = 1, self%nrays
     !tmpsum = 0.0    
     GATE_LOOP : do igate = 9, self%ngates

        ! From Eq. (1) in Tromel et al. 2013 (JAMC)
        if ((STDV_Zhh(igate,iray) > -9998).and.(SNR_Zhh(igate,iray) > snr_mask_threshold)) then

           if (zhh(igate,iray) < 20.) then
               !ivals = (/ igate-8,igate-7,igate-6,igate-5,igate-4,igate-3,igate-2,igate-1,igate /)
               id1 = 8
           else if (zhh(igate,iray) >= 20. .and. zhh(igate,iray) < 35.) then
               !ivals = (/ igate-4,igate-3,igate-2,igate-1,igate /)
               id1 = 4
           else if (zhh(igate,iray) >= 35.) then
               !ivals = (/ igate-2,igate-1,igate /)
               id1 = 2
           else
               !ivals = (/ igate-2,igate-1,igate /)
               id1 = 2
           endif

           rng_diff = ( self%range(igate) - self%range(igate-id1) ) * 1.E-3 ! in km
           !print *, rng_diff, missing_val
           
           if (igate .eq. 9) then ! .and. kdp(igate,iray) > -1000) then
              xarr = (/ ( 2. * rng_diff ) * kdp(igate,iray), dphase(igate,iray), dble(0.0) /)
              xarr = dble(xarr)

               if (xarr(1) > dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = SUM(xarr)
                  phidp(igate,iray) = xarr(1)
               else if (xarr(1) > dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) < dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(1) + xarr(2)
                  phidp(igate,iray) = xarr(1)
               else if (xarr(1) > dble(-9000.0) .and. xarr(2) < dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(1) + xarr(3)
                  phidp(igate,iray) = xarr(1)
               else if (xarr(1) > dble(-9000.0) .and. xarr(2) < dble(-9000.0) .and. xarr(3) < dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(1)
                  phidp(igate,iray) = xarr(1)
               else if (xarr(1) < dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(2) + xarr(3)
                  phidp(igate,iray) = missing_val
               else if (xarr(1) < dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) < dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(2)
                  phidp(igate,iray) = missing_val
               else if (xarr(1) < dble(-9000.0) .and. xarr(2) < dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = missing_val
                  phidp(igate,iray) = missing_val
               else
                  phidp_wdp(igate,iray) = missing_val
                  phidp(igate,iray) = missing_val
               endif
             
              !print *, SUM(xarr, MASK=DMOD(xarr,missing_val)==0)
              !phidp(igate,iray) = SUM(xarr, MASK=DMOD(xarr, missing_val)==1)
               print *, xarr, phidp(igate,iray)
               !call sleep(3)
           else if (igate .gt. 9) then ! .and. kdp(igate,iray) > -1000) then
               !print *, iray, igate
              xarr = (/ ( 2. * rng_diff ) * ( kdp(igate,iray) ), phidp_wdp(igate-1, iray), dphase(igate,iray) /)
              tarr = (/ ( 2. * rng_diff ) * ( kdp(igate,iray) ), phidp(igate-1, iray), dble(0.0) /)
              xarr = dble(xarr)
              tarr = dble(tarr)
               if (xarr(1) > dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = SUM(xarr)
               else if (xarr(1) > dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) < dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(1) + xarr(2)
               else if (xarr(1) > dble(-9000.0) .and. xarr(2) < dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(1) + xarr(3)
               else if (xarr(1) > dble(-9000.0) .and. xarr(2) < dble(-9000.0) .and. xarr(3) < dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(1)
               else if (xarr(1) < dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(2) + xarr(3)
               else if (xarr(1) < dble(-9000.0) .and. xarr(2) > dble(-9000.0) .and. xarr(3) < dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(2)
               else if (xarr(1) < dble(-9000.0) .and. xarr(2) < dble(-9000.0) .and. xarr(3) > dble(-9000.0)) then
                  phidp_wdp(igate,iray) = xarr(3)
               else
                  phidp_wdp(igate,iray) = missing_val
               endif

               if (tarr(1) > dble(-9000.0) .and. tarr(2) > dble(-9000.0)) then
                  phidp(igate,iray) = tarr(1) + tarr(2)
               else if (tarr(1) > dble(-9000.0) .and. tarr(2) < dble(-9000.0)) then
                  phidp(igate,iray) = tarr(1)
               else if (tarr(1) < dble(-9000.0) .and. tarr(2) > dble(-9000.0)) then
                  phidp(igate,iray) = tarr(2)
               else if (tarr(1) < dble(-9000.0) .and. tarr(2) < dble(-9000.0)) then
                  phidp(igate,iray) = missing_val
               else
                  phidp(igate,iray) = missing_val
               endif

               print *, iray, igate
               print *, xarr, phidp_wdp(igate,iray)
               print *, tarr, phidp(igate,iray)

               !call sleep(2)
        
           else
              phidp(igate,iray) = missing_val
              phidp_wdp(igate,iray) = missing_val
           endif
           !phidp(igate,iray) = ( 2.0 * ( 1.E-3 * self%range(igate) ) * ( kdp(igate,iray) - kdp(1,iray) ) ) + dphase(igate,iray)
           !if (phidp(igate,iray) > 180.) then
           !   phidp(igate,iray) = phidp(igate,iray)-360.
           !endif
        else
           phidp(igate,iray) = missing_val
           phidp_wdp(igate,iray) = missing_val
        endif

     enddo GATE_LOOP
  enddo RAY_LOOP

 end subroutine compute_phidp

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------


  FUNCTION random_normal() RESULT (ran_norm)

    ! Adapted from the following Fortran 77 code
    !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
    !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

    !  The function random_normal() returns a normally distributed pseudo-random
    !  number with zero mean and unit variance.   This version uses the default
    !  uniform random number generator which is in your fortran library.

    !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
    !  and J.F. Monahan augmented with quadratic bounding curves.

    !  Fortran 90 version by Alan Miller (alan @ mel.dms.csiro.au)

    IMPLICIT NONE
    REAL :: ran_norm

    !     Local variables
    REAL, PARAMETER :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,   &
         half = 0.5, r1 = 0.27597, r2 = 0.27846
    REAL            :: u, v, x, y, q

    !     Generate P = (u,v) uniform in rectangle enclosing acceptance region

    DO
        CALL RANDOM_NUMBER(u)
        CALL RANDOM_NUMBER(v)
        v = 1.7156 * (v - half)

        !     Evaluate the quadratic form
        x = u - s
        y = ABS(v) - t
        q = x**2 + y*(a*y - b*x)

        !     Accept P if inside inner ellipse
        IF (q < r1) EXIT
        !     Reject P if outside outer ellipse
        IF (q > r2) CYCLE
        !     Reject P if outside acceptance region
        IF (v**2 < -4.0*LOG(u)*u**2) EXIT
    END DO

    !     Return ratio of P's coordinates as the normal deviate
    ran_norm = v/u
    RETURN

  END FUNCTION random_normal

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
  
  subroutine error_handler(ierr,message)
    use netcdf, only : NF90_NOERR
    use netcdf, only : nf90_strerror
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message
    if ( ierr == NF90_NOERR ) return
    write(*,'("NETCDF >> ", A)') nf90_strerror(ierr)
    write(*,'("MODULE_CFRADIAL_OUTPUT: ",A)') message
    stop
  end subroutine error_handler

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

end module module_cfradial_output
