/***************************************************/
/* CONSTANTS */
/***************************************************/
#define IDENT_LEN 4
#define PROJ_NAME_LEN 20
#define FLIGHT_NUM_LEN 8
#define FAC_NAME_LEN 8
#define RAD_NAME_LEN 8
#define MAX_NUM_PARMS 40 
#define PARM_NAME_LEN 8 
#define PARM_DESC_LEN 40
#define PARM_UNIT_LEN 8
#define THRESHOLD_FLD_LEN 8
#define MAX_BEAMS 1300 
#define MAX_GATES 1500
#define MAX_KEYS 8
#define SIGN16 0x8000
# define SHORTS_TO_INTS(x) ((((x)-1)>>1)+1)
# define INTS_TO_BYTES(x) ((x)<<2)
/***************************************************/
/* STRUCTURES */
/***************************************************/
// Required pragma to support 64-bit operation on 32-bit files
#pragma pack(push)
#pragma pack(4)
struct key_table_info {
    int offset;
    int size;
    int type;
};

struct sswb_info {
    /* parameters from the first version */
    int last_used;		/* Unix time */
    int start_time;
    int stop_time;
    int sizeof_file;
    int compression_flag;
    int volume_time_stamp;	/* to reference current volume */
    int num_params;		/* number of parameters */
	
    /* end of first version parameters */
	
    char radar_name[8];
	
    double d_start_time;
    double d_stop_time;
    /*
     * "last_used" is an age off indicator where > 0 implies Unix time
     * of the last access and
     * 0 implies this sweep should not be aged off
     */
    int version_num;
    int num_key_tables;
    int status;
    int place_holder[7];
    struct key_table_info key_table[MAX_KEYS];
    /*
     * offset and key info to a table containing key value such as
     * the rot. angle and the offset to the corresponding ray
     * in the disk file
     */

};
struct vold_info {

   short ver_num;
   short vol_num;
   int max_bytes;
   char proj_name[PROJ_NAME_LEN];
   short year;
   short mon;
   short day;
   short hour;
   short min;
   short sec;
   char flight_num[FLIGHT_NUM_LEN];
   char gen_fac_name[FAC_NAME_LEN];
   short gen_year;
   short gen_mon;
   short gen_day;
   short num_sensors;

};
struct radd_info {
   
   char rad_name[RAD_NAME_LEN];
   float rad_constant;
   float peak_pow;
   float noise_pos;
   float rec_gain;
   float ant_gain;
   float sys_gain;
   float horiz_beam_width;
   float vert_beam_width;
   short rad_type;
   short scan_mode;
   float scan_rate;
   float start_ang;
   float stop_ang;
   short num_param_desc;
   short num_desc;
   short compress_flag;
   short data_reduc_flag;
   float data_reduc_parm1;
   float data_reduc_parm2;
   float radar_lon;
   float radar_lat;
   float radar_alt;
   float unambig_vel;
   float unambig_range;
   short num_freq;
   short num_pulse;
   float freq1;
   float freq2;
   float freq3;
   float freq4;
   float freq5;
   float pulse_per1;
   float pulse_per2;
   float pulse_per3;
   float pulse_per4;
   float pulse_per5;

};

struct cfac_info {
   
   float c_azimuth;
   float c_elevation;
   float c_range_delay;
   float c_rad_lon;
   float c_rad_lat;
   float c_alt_msl;
   float c_alt_agl;
   float c_ew_grspeed;
   float c_ns_grspeed;
   float c_vert_vel;
   float c_head;
   float c_roll;
   float c_pitch;
   float c_drift;
   float c_rotang;
   float c_tiltang;

};

struct parm_info {

   char parm_name[PARM_NAME_LEN];
   char parm_desc[PARM_DESC_LEN];
   char parm_unit[PARM_UNIT_LEN];
   short ipp;
   short trans_freq;
   float rec_bandwidth;
   short pulse_width;
   short polarization;
   short num_samples;
   short parm_type;
   char threshold_fld[THRESHOLD_FLD_LEN];
   float threshold_val;
   float scale_fac;
   float offset_fac;
   int baddata_flag;
    
};

struct celv_info {

   int total_gates;
   float gate_spacing[MAX_GATES];

};

struct swib_info {

   char rad_name[RAD_NAME_LEN];
   int sweep_num;
   int num_rays;
   float start_ang;
   float stop_ang;
   float fixed_ang;
   int filter_flag;
   
};

struct ryib_info {

   int sweep_num;
   int julian_day;
   short hour;
   short min;
   short sec;
   short msec;
   float azimuth;
   float elevation;
   float peak_power;
   float scan_rate;
   int ray_status;

};

struct asib_info {

   float lon;
   float lat;
   float alt_msl;
   float alt_agl;
   float ew_gspeed;
   float ns_gspeed;
   float vert_vel;
   float head;
   float roll;
   float pitch;
   float drift;
   float rot_ang;
   float tilt_ang;
   float ew_horiz_wind;
   float ns_horiz_wind;
   float vert_wind;
   float head_change;
   float pitch_change;

};

struct rdat_info {

   char parm_name[PARM_NAME_LEN];
   float data[MAX_GATES];

};

struct radar_angles {
    /* all angles are in radians
     */
    float azimuth;
    float elevation;
    float x;
    float y;
    float z;
    float psi;
    float rotation_angle;
    float tilt;
};

struct rot_table_entry {
    float rotation_angle;
    int offset;
    int size;
};

struct rktb_info {
    float angle2ndx;
    int ndx_que_size;
    int first_key_offset;
    int angle_table_offset;
    int num_rays;
};
#pragma pack(pop)
//#pragma options align=reset
