#include "Args.h"

#include <cstdarg>
#include <iostream>
#include <fstream>
//#include <regex>

bool Args::parseArgs(int argc, char *argv[]) {

  std::unordered_set<std::string> keywordSet ( { } );

  tdrp_override_t override;
  TDRP_init_override(&override);

  char tmp_str[BUFSIZ];

  // Fill override from command line arguments

  for (int i = 1; i < argc; i++) {
    const char *arg = argv[i] + 1;
    std::unordered_set<std::string>::const_iterator found = keywordSet.find(arg);
    if ( found != keywordSet.end() ) {
      sprintf(tmp_str, "%s = %s;", arg, argv[++i]);
      TDRP_add_override(&override, tmp_str);      
    } else {
      // This one has specific meaning to TDRP.
      if ( ! strcmp(argv[i], "-params" ) )
	i++;
    }
  }

  // Now parse the config file
  
  char *paramsPath;
  
  if (params.loadFromArgs(argc, argv, override.list, &paramsPath))
    badparms("Problem with command line args");

  TDRP_free_override(&override);

  // and finally process the arguments
  return params_to_hash();
}


// Fill the given configHash with the data from the params object

#define CONFIG_INSERT_INT(X) configHash->insert(X, params.X.toInt())
#define CONFIG_INSERT_FLOAT(X) configHash->insert(X, params.X.toFloat())
#define CONFIG_INSERT_STR(X) configHash->insert(X, params.X)

bool params_to_hash(QHash<QString, QString> *configHash) {

  // string arguments

  CONFIG_INSERT_STRING("adjust_background");
  CONFIG_INSERT_STRING("allow_background_missing_values");
  CONFIG_INSERT_STRING("allow_negative_angles");
  CONFIG_INSERT_STRING("array_order");
  CONFIG_INSERT_STRING("bg_interpolation");
  CONFIG_INSERT_STRING("bkgd_obs_interpolation");
  CONFIG_INSERT_STRING("data_directory");
  CONFIG_INSERT_STRING("debug_bgU_nc");
  CONFIG_INSERT_STRING("debug_bgU_overwrite");
  CONFIG_INSERT_STRING("fractl_nc_file");
  CONFIG_INSERT_STRING("horizontal_radar_appx");
  CONFIG_INSERT_STRING("load_background");
  CONFIG_INSERT_STRING("load_bg_coefficients");
  CONFIG_INSERT_STRING("mask_reflectivity");
  CONFIG_INSERT_STRING("mode");
  CONFIG_INSERT_STRING("output_asi");
  CONFIG_INSERT_STRING("output_COAMPS");
  CONFIG_INSERT_STRING("output_directory");
  CONFIG_INSERT_STRING("output_mish");
  CONFIG_INSERT_STRING("output_netcdf");
  CONFIG_INSERT_STRING("output_qc");
  CONFIG_INSERT_STRING("output_txt");
  CONFIG_INSERT_STRING("preprocess_obs");
  CONFIG_INSERT_STRING("projection");
  CONFIG_INSERT_STRING("qr_variable");
  CONFIG_INSERT_STRING("radar_dbz");
  CONFIG_INSERT_STRING("radar_sw");
  CONFIG_INSERT_STRING("radar_vel");
  CONFIG_INSERT_STRING("ref_state");
  CONFIG_INSERT_STRING("save_mish");
  CONFIG_INSERT_STRING("use_fractl_errors");

  
  // int arguments

  CONFIG_INSERT_INT("bkgd_kd_num_neighbors");
  CONFIG_INSERT_INT("debug_kd");
  CONFIG_INSERT_INT("debug_kd_step");
  CONFIG_INSERT_INT("dynamic_stride");
  CONFIG_INSERT_INT("num_iterations");
  CONFIG_INSERT_INT("radar_skip");
  CONFIG_INSERT_INT("radar_stride");
  CONFIG_INSERT_INT("ref_time");
  CONFIG_INSERT_INT("spline_approximation");

  // float arguments

  CONFIG_INSERT_FLOAT("aeri_qv_error");
  CONFIG_INSERT_FLOAT("aeri_rhoa_error");
  CONFIG_INSERT_FLOAT("aeri_rhou_error");
  CONFIG_INSERT_FLOAT("aeri_rhov_error");
  CONFIG_INSERT_FLOAT("aeri_rhow_error");
  CONFIG_INSERT_FLOAT("aeri_tempk_error");
  CONFIG_INSERT_FLOAT("amv_rhou_error");
  CONFIG_INSERT_FLOAT("amv_rhov_error");
  CONFIG_INSERT_FLOAT("ascat_rhou_error");
  CONFIG_INSERT_FLOAT("ascat_rhov_error");
  CONFIG_INSERT_FLOAT("bg_obs_error");
  CONFIG_INSERT_FLOAT("bg_interpolation_error");    
  CONFIG_INSERT_FLOAT("bg_qr_error");
  CONFIG_INSERT_FLOAT("bg_qv_error");
  CONFIG_INSERT_FLOAT("bg_rhoa_error");
  CONFIG_INSERT_FLOAT("bg_rhou_error");
  CONFIG_INSERT_FLOAT("bg_rhov_error");
  CONFIG_INSERT_FLOAT("bg_rhow_error");
  CONFIG_INSERT_FLOAT("bg_tempk_error");
  CONFIG_INSERT_FLOAT("bkgd_kd_max_distance");
  CONFIG_INSERT_FLOAT("dbz_pseudow_weight");
  CONFIG_INSERT_FLOAT("dropsonde_qv_error");
  CONFIG_INSERT_FLOAT("dropsonde_rhoa_error");
  CONFIG_INSERT_FLOAT("dropsonde_rhou_error");
  CONFIG_INSERT_FLOAT("dropsonde_rhov_error");
  CONFIG_INSERT_FLOAT("dropsonde_rhow_error");
  CONFIG_INSERT_FLOAT("dropsonde_tempk_error");
  CONFIG_INSERT_FLOAT("flightlevel_qv_error");
  CONFIG_INSERT_FLOAT("flightlevel_rhoa_error");
  CONFIG_INSERT_FLOAT("flightlevel_rhou_error");
  CONFIG_INSERT_FLOAT("flightlevel_rhov_error");
  CONFIG_INSERT_FLOAT("flightlevel_rhow_error");
  CONFIG_INSERT_FLOAT("flightlevel_tempk_error");
  CONFIG_INSERT_FLOAT("i_background_roi");
  CONFIG_INSERT_FLOAT("i_filter_length");
  CONFIG_INSERT_FLOAT("i_incr");
  CONFIG_INSERT_FLOAT("i_max");
  CONFIG_INSERT_FLOAT("i_max_wavenumber");
  CONFIG_INSERT_FLOAT("i_max_wavenumber_qr");
  CONFIG_INSERT_FLOAT("i_max_wavenumber_qv");
  CONFIG_INSERT_FLOAT("i_max_wavenumber_rhoa");
  CONFIG_INSERT_FLOAT("i_max_wavenumber_rhou");
  CONFIG_INSERT_FLOAT("i_max_wavenumber_rhov");
  CONFIG_INSERT_FLOAT("i_max_wavenumber_rhow");
  CONFIG_INSERT_FLOAT("i_max_wavenumber_tempk");
  CONFIG_INSERT_FLOAT("i_min");
  CONFIG_INSERT_FLOAT("i_reflectivity_roi");
  CONFIG_INSERT_FLOAT("i_spline_cutoff");
  CONFIG_INSERT_FLOAT("insitu_qv_error");
  CONFIG_INSERT_FLOAT("insitu_rhoa_error");
  CONFIG_INSERT_FLOAT("insitu_rhou_error");
  CONFIG_INSERT_FLOAT("insitu_rhov_error");
  CONFIG_INSERT_FLOAT("insitu_rhow_error");
  CONFIG_INSERT_FLOAT("insitu_tempk_error");
  CONFIG_INSERT_FLOAT("j_background_roi");
  CONFIG_INSERT_FLOAT("j_filter_length");
  CONFIG_INSERT_FLOAT("j_incr");
  CONFIG_INSERT_FLOAT("j_max");
  CONFIG_INSERT_FLOAT("j_max_wavenumber");
  CONFIG_INSERT_FLOAT("j_max_wavenumber_qr");
  CONFIG_INSERT_FLOAT("j_max_wavenumber_qv");
  CONFIG_INSERT_FLOAT("j_max_wavenumber_rhoa");
  CONFIG_INSERT_FLOAT("j_max_wavenumber_rhou");
  CONFIG_INSERT_FLOAT("j_max_wavenumber_rhov");
  CONFIG_INSERT_FLOAT("j_max_wavenumber_rhow");
  CONFIG_INSERT_FLOAT("j_max_wavenumber_tempk");
  CONFIG_INSERT_FLOAT("j_min");
  CONFIG_INSERT_FLOAT("j_reflectivity_roi");
  CONFIG_INSERT_FLOAT("j_spline_cutoff");
  CONFIG_INSERT_FLOAT("k_filter_length");
  CONFIG_INSERT_FLOAT("k_incr");
  CONFIG_INSERT_FLOAT("k_max");
  CONFIG_INSERT_FLOAT("k_max_wavenumber");
  CONFIG_INSERT_FLOAT("k_max_wavenumber_qr");
  CONFIG_INSERT_FLOAT("k_max_wavenumber_qv");
  CONFIG_INSERT_FLOAT("k_max_wavenumber_rhoa");
  CONFIG_INSERT_FLOAT("k_max_wavenumber_rhou");
  CONFIG_INSERT_FLOAT("k_max_wavenumber_rhov");
  CONFIG_INSERT_FLOAT("k_max_wavenumber_rhow");
  CONFIG_INSERT_FLOAT("k_max_wavenumber_tempk");
  CONFIG_INSERT_FLOAT("k_min");
  CONFIG_INSERT_FLOAT("k_reflectivity_roi");
  CONFIG_INSERT_FLOAT("k_spline_cutoff");
  CONFIG_INSERT_FLOAT("lidar_min_error");
  CONFIG_INSERT_FLOAT("lidar_power_error");
  CONFIG_INSERT_FLOAT("lidar_sw_error");
  CONFIG_INSERT_FLOAT("max_radar_elevation");
  CONFIG_INSERT_FLOAT("mc_weight");
  CONFIG_INSERT_FLOAT("melting_zone_width");
  CONFIG_INSERT_FLOAT("mesonet_qv_error");
  CONFIG_INSERT_FLOAT("mesonet_rhoa_error");
  CONFIG_INSERT_FLOAT("mesonet_rhou_error");
  CONFIG_INSERT_FLOAT("mesonet_rhov_error");
  CONFIG_INSERT_FLOAT("mesonet_rhow_error");
  CONFIG_INSERT_FLOAT("mesonet_tempk_error");
  CONFIG_INSERT_FLOAT("mixed_phase_dbz");
  CONFIG_INSERT_FLOAT("mtp_rhoa_error");
  CONFIG_INSERT_FLOAT("mtp_tempk_error");
  CONFIG_INSERT_FLOAT("output_latlon_increment");
  CONFIG_INSERT_FLOAT("output_pressure_increment");
  CONFIG_INSERT_FLOAT("qscat_rhou_error");
  CONFIG_INSERT_FLOAT("qscat_rhov_error");
  CONFIG_INSERT_FLOAT("radar_fallspeed_error");
  CONFIG_INSERT_FLOAT("radar_min_error");
  CONFIG_INSERT_FLOAT("radar_sw_error");
  CONFIG_INSERT_FLOAT("rain_dbz");
  CONFIG_INSERT_FLOAT("ref_lat");
  CONFIG_INSERT_FLOAT("ref_lon");
  CONFIG_INSERT_FLOAT("sfmr_windspeed_error");

}
