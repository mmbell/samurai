#include "Args.h"

#include <cstdarg>
#include <iostream>
#include <fstream>

#include <QVariant>

//#include <regex>

bool Args::parseArgs(int argc, char *argv[])
{

  bool ret_val = true;
  
  std::unordered_set<std::string> keywordSet ( { } );

  tdrp_override_t override;
  TDRP_init_override(&override);

  char tmp_str[BUFSIZ];

  // Fill override from command line arguments

  for (int i = 1; i < argc; i++) {
    if ( !strcmp(argv[i], "--") ||
	 !strcmp(argv[i], "-h") ||
	 !strcmp(argv[i], "-help") ||
	 !strcmp(argv[i], "-man")) {
      params.usage(std::cout);
      exit (0);
    }
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
  
  if (params.loadFromArgs(argc, argv, override.list, &paramsPath)) {
    std::cout << "Problem with command line args" << std::endl;
    ret_val = false;
  }

  TDRP_free_override(&override);

  return ret_val;
}


// Fill the given configHash with the data from the params object
#define STR_HELPER(_X_) QString(#_X_)
#define STR(_X_) STR_HELPER(_X_)

#define CONFIG_INSERT_BOOL(_X_) configHash->insert(#_X_, params._X_ ? "true" : "false")
#define CONFIG_INSERT_INT(_X_) configHash->insert(#_X_, QVariant(params._X_).toString())
#define CONFIG_INSERT_FLOAT(_X_) configHash->insert(#_X_, QVariant(params._X_).toString())
#define CONFIG_INSERT_STR(_X_) configHash->insert(#_X_, params._X_)
#define CONFIG_INSERT_FLOAT_ARRAY(_X_, _iter_) \
  configHash->insert(STR(_X_) + "_" + QVariant(iter).toString(), QVariant(params._X_[_iter_]).toString())

#define CONFIG_INSERT_MAP_VALUE(_X_, _table_) configHash->insert(#_X_, _table_[params._X_])

// This needs to match bkgd_interp_method_t in paramdef.samurai

QString interp_map[] = { "none", "spline", "kd_tree", "fractl" };

// This needs to match mode_t in paramdef.samurai

QString mode_map[] = { "XYZ", "RTZ" };

// This needs to match projection_t in paramdef.samurai

QString projection_map[] = { "PROJ_LAMBERT_CONFORMAL_CONIC",
				 "PROJ_TRANSVERSE_MERCATOR_EXACT" };

bool Args::paramsToHash(QHash<QString, QString> *configHash) {

  // string arguments

  CONFIG_INSERT_BOOL(adjust_background);
  CONFIG_INSERT_BOOL(allow_background_missing_values);
  CONFIG_INSERT_BOOL(allow_negative_angles);
  CONFIG_INSERT_STR(array_order);
  CONFIG_INSERT_STR(bg_interpolation);
  CONFIG_INSERT_MAP_VALUE(bkgd_obs_interpolation, interp_map);
  CONFIG_INSERT_STR(data_directory);
  CONFIG_INSERT_STR(debug_bgu_nc);
  CONFIG_INSERT_STR(debug_bgu_overwrite);
  CONFIG_INSERT_STR(fractl_nc_file);
  CONFIG_INSERT_BOOL(horizontal_radar_appx);
  CONFIG_INSERT_BOOL(load_background);
  CONFIG_INSERT_BOOL(load_bg_coefficients);
  CONFIG_INSERT_FLOAT(mask_reflectivity);
  CONFIG_INSERT_MAP_VALUE(mode, mode_map);
  CONFIG_INSERT_BOOL(output_asi);
  CONFIG_INSERT_BOOL(output_COAMPS);
  CONFIG_INSERT_STR(output_directory);
  CONFIG_INSERT_BOOL(output_mish);
  CONFIG_INSERT_BOOL(output_netcdf);
  CONFIG_INSERT_BOOL(output_qc);
  CONFIG_INSERT_BOOL(output_txt);
  CONFIG_INSERT_BOOL(preprocess_obs);
  CONFIG_INSERT_MAP_VALUE(projection, projection_map);
  CONFIG_INSERT_STR(qr_variable);
  CONFIG_INSERT_STR(radar_dbz);
  CONFIG_INSERT_STR(radar_sw);
  CONFIG_INSERT_STR(radar_vel);
  CONFIG_INSERT_STR(ref_state);
  CONFIG_INSERT_BOOL(save_mish);
  CONFIG_INSERT_BOOL(use_fractl_errors);
  
  // int arguments

  CONFIG_INSERT_INT(bkgd_kd_num_neighbors);
  CONFIG_INSERT_INT(debug_kd);
  CONFIG_INSERT_INT(debug_kd_step);
  CONFIG_INSERT_INT(dynamic_stride);
  CONFIG_INSERT_INT(num_iterations);
  CONFIG_INSERT_INT(radar_skip);
  CONFIG_INSERT_INT(radar_stride);
  CONFIG_INSERT_INT(ref_time);

  CONFIG_INSERT_INT(i_qr_bcL);
  CONFIG_INSERT_INT(i_qr_bcR);
  CONFIG_INSERT_INT(i_qv_bcL);
  CONFIG_INSERT_INT(i_qv_bcR);
  CONFIG_INSERT_INT(i_rhoa_bcL);
  CONFIG_INSERT_INT(i_rhoa_bcR);
  CONFIG_INSERT_INT(i_rhou_bcL);
  CONFIG_INSERT_INT(i_rhou_bcR);
  CONFIG_INSERT_INT(i_rhov_bcL);
  CONFIG_INSERT_INT(i_rhov_bcR);
  CONFIG_INSERT_INT(i_rhow_bcL);
  CONFIG_INSERT_INT(i_rhow_bcR);
  CONFIG_INSERT_INT(i_tempk_bcL);
  CONFIG_INSERT_INT(i_tempk_bcR);
  CONFIG_INSERT_INT(j_qr_bcL);
  CONFIG_INSERT_INT(j_qr_bcR);
  CONFIG_INSERT_INT(j_qv_bcL);
  CONFIG_INSERT_INT(j_qv_bcR);
  CONFIG_INSERT_INT(j_rhoa_bcL);
  CONFIG_INSERT_INT(j_rhoa_bcR);
  CONFIG_INSERT_INT(j_rhou_bcL);
  CONFIG_INSERT_INT(j_rhou_bcR);
  CONFIG_INSERT_INT(j_rhov_bcL);
  CONFIG_INSERT_INT(j_rhov_bcR);
  CONFIG_INSERT_INT(j_rhow_bcL);
  CONFIG_INSERT_INT(j_rhow_bcR);
  CONFIG_INSERT_INT(j_tempk_bcL);
  CONFIG_INSERT_INT(j_tempk_bcR);
  CONFIG_INSERT_INT(k_qr_bcL);
  CONFIG_INSERT_INT(k_qr_bcR);
  CONFIG_INSERT_INT(k_qv_bcL);
  CONFIG_INSERT_INT(k_qv_bcR);
  CONFIG_INSERT_INT(k_rhoa_bcL);
  CONFIG_INSERT_INT(k_rhoa_bcR);
  CONFIG_INSERT_INT(k_rhou_bcL);
  CONFIG_INSERT_INT(k_rhou_bcR);
  CONFIG_INSERT_INT(k_rhov_bcL);
  CONFIG_INSERT_INT(k_rhov_bcR);
  CONFIG_INSERT_INT(k_rhow_bcL);
  CONFIG_INSERT_INT(k_rhow_bcR);
  CONFIG_INSERT_INT(k_tempk_bcL);
  CONFIG_INSERT_INT(k_tempk_bcR);
  
  // float arguments

  CONFIG_INSERT_FLOAT(aeri_qv_error);
  CONFIG_INSERT_FLOAT(aeri_rhoa_error);
  CONFIG_INSERT_FLOAT(aeri_rhou_error);
  CONFIG_INSERT_FLOAT(aeri_rhov_error);
  CONFIG_INSERT_FLOAT(aeri_rhow_error);
  CONFIG_INSERT_FLOAT(aeri_tempk_error);
  CONFIG_INSERT_FLOAT(amv_rhou_error);
  CONFIG_INSERT_FLOAT(amv_rhov_error);
  CONFIG_INSERT_FLOAT(ascat_rhou_error);
  CONFIG_INSERT_FLOAT(ascat_rhov_error);
  CONFIG_INSERT_FLOAT(bg_obs_error);
  CONFIG_INSERT_FLOAT(bg_interpolation_error);
  CONFIG_INSERT_FLOAT(bkgd_kd_max_distance);
  CONFIG_INSERT_FLOAT(dbz_pseudow_weight);
  CONFIG_INSERT_FLOAT(dropsonde_qv_error);
  CONFIG_INSERT_FLOAT(dropsonde_rhoa_error);
  CONFIG_INSERT_FLOAT(dropsonde_rhou_error);
  CONFIG_INSERT_FLOAT(dropsonde_rhov_error);
  CONFIG_INSERT_FLOAT(dropsonde_rhow_error);
  CONFIG_INSERT_FLOAT(dropsonde_tempk_error);
  CONFIG_INSERT_FLOAT(flightlevel_qv_error);
  CONFIG_INSERT_FLOAT(flightlevel_rhoa_error);
  CONFIG_INSERT_FLOAT(flightlevel_rhou_error);
  CONFIG_INSERT_FLOAT(flightlevel_rhov_error);
  CONFIG_INSERT_FLOAT(flightlevel_rhow_error);
  CONFIG_INSERT_FLOAT(flightlevel_tempk_error);
  CONFIG_INSERT_FLOAT(i_background_roi);
  // CONFIG_INSERT_FLOAT(i_filter_length);
  CONFIG_INSERT_FLOAT(i_incr);
  CONFIG_INSERT_FLOAT(i_max);
#if 0  
  CONFIG_INSERT_FLOAT(i_max_wavenumber);
  CONFIG_INSERT_FLOAT(i_max_wavenumber_qr);
  CONFIG_INSERT_FLOAT(i_max_wavenumber_qv);
  CONFIG_INSERT_FLOAT(i_max_wavenumber_rhoa);
  CONFIG_INSERT_FLOAT(i_max_wavenumber_rhou);
  CONFIG_INSERT_FLOAT(i_max_wavenumber_rhov);
  CONFIG_INSERT_FLOAT(i_max_wavenumber_rhow);
  CONFIG_INSERT_FLOAT(i_max_wavenumber_tempk);
#endif
  CONFIG_INSERT_FLOAT(i_min);
  CONFIG_INSERT_FLOAT(i_reflectivity_roi);
  // CONFIG_INSERT_FLOAT(i_spline_cutoff);
  CONFIG_INSERT_FLOAT(insitu_qv_error);
  CONFIG_INSERT_FLOAT(insitu_rhoa_error);
  CONFIG_INSERT_FLOAT(insitu_rhou_error);
  CONFIG_INSERT_FLOAT(insitu_rhov_error);
  CONFIG_INSERT_FLOAT(insitu_rhow_error);
  CONFIG_INSERT_FLOAT(insitu_tempk_error);
  CONFIG_INSERT_FLOAT(j_background_roi);
  // CONFIG_INSERT_FLOAT(j_filter_length);
  CONFIG_INSERT_FLOAT(j_incr);
  CONFIG_INSERT_FLOAT(j_max);
#if 0  
  CONFIG_INSERT_FLOAT(j_max_wavenumber);
  CONFIG_INSERT_FLOAT(j_max_wavenumber_qr);
  CONFIG_INSERT_FLOAT(j_max_wavenumber_qv);
  CONFIG_INSERT_FLOAT(j_max_wavenumber_rhoa);
  CONFIG_INSERT_FLOAT(j_max_wavenumber_rhou);
  CONFIG_INSERT_FLOAT(j_max_wavenumber_rhov);
  CONFIG_INSERT_FLOAT(j_max_wavenumber_rhow);
  CONFIG_INSERT_FLOAT(j_max_wavenumber_tempk);
#endif
  CONFIG_INSERT_FLOAT(j_min);
  CONFIG_INSERT_FLOAT(j_reflectivity_roi);
  // CONFIG_INSERT_FLOAT(j_spline_cutoff);
  // CONFIG_INSERT_FLOAT(k_filter_length);
  CONFIG_INSERT_FLOAT(k_incr);
  CONFIG_INSERT_FLOAT(k_max);
#if 0
  CONFIG_INSERT_FLOAT(k_max_wavenumber);
  CONFIG_INSERT_FLOAT(k_max_wavenumber_qr);
  CONFIG_INSERT_FLOAT(k_max_wavenumber_qv);
  CONFIG_INSERT_FLOAT(k_max_wavenumber_rhoa);
  CONFIG_INSERT_FLOAT(k_max_wavenumber_rhou);
  CONFIG_INSERT_FLOAT(k_max_wavenumber_rhov);
  CONFIG_INSERT_FLOAT(k_max_wavenumber_rhow);
  CONFIG_INSERT_FLOAT(k_max_wavenumber_tempk);
#endif  
  CONFIG_INSERT_FLOAT(k_min);
  CONFIG_INSERT_FLOAT(k_reflectivity_roi);
  // CONFIG_INSERT_FLOAT(k_spline_cutoff);
  CONFIG_INSERT_FLOAT(lidar_min_error);
  CONFIG_INSERT_FLOAT(lidar_power_error);
  CONFIG_INSERT_FLOAT(lidar_sw_error);
  CONFIG_INSERT_FLOAT(max_radar_elevation);
  // CONFIG_INSERT_FLOAT(mc_weight);
  CONFIG_INSERT_FLOAT(melting_zone_width);
  CONFIG_INSERT_FLOAT(mesonet_qv_error);
  CONFIG_INSERT_FLOAT(mesonet_rhoa_error);
  CONFIG_INSERT_FLOAT(mesonet_rhou_error);
  CONFIG_INSERT_FLOAT(mesonet_rhov_error);
  CONFIG_INSERT_FLOAT(mesonet_rhow_error);
  CONFIG_INSERT_FLOAT(mesonet_tempk_error);
  CONFIG_INSERT_FLOAT(mixed_phase_dbz);
  CONFIG_INSERT_FLOAT(mtp_rhoa_error);
  CONFIG_INSERT_FLOAT(mtp_tempk_error);
  CONFIG_INSERT_FLOAT(output_latlon_increment);
  CONFIG_INSERT_FLOAT(output_pressure_increment);
  CONFIG_INSERT_FLOAT(qscat_rhou_error);
  CONFIG_INSERT_FLOAT(qscat_rhov_error);
  CONFIG_INSERT_FLOAT(radar_fallspeed_error);
  CONFIG_INSERT_FLOAT(radar_min_error);
  CONFIG_INSERT_FLOAT(radar_sw_error);
  CONFIG_INSERT_FLOAT(rain_dbz);
  CONFIG_INSERT_FLOAT(sfmr_windspeed_error);

  for (int iter = 1; iter <= params.num_iterations; iter++) {
    CONFIG_INSERT_FLOAT_ARRAY(_bg_qr_error, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_bg_rhoa_error, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_bg_rhou_error, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_bg_rhov_error, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_bg_rhow_error, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_bg_tempk_error, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_bg_qv_error, iter);

    CONFIG_INSERT_FLOAT_ARRAY(_i_filter_length, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_j_filter_length, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_k_filter_length, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_i_spline_cutoff, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_j_spline_cutoff, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_k_spline_cutoff, iter);    
    CONFIG_INSERT_FLOAT_ARRAY(_i_max_wavenumber, iter);
    CONFIG_INSERT_FLOAT_ARRAY(_j_max_wavenumber, iter);
  }
  return true;
}
