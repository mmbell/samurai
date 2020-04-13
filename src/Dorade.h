/*
 *  Dorade.h
 *  tcvar
 *
 *  Copyright 2008 Michael Bell
 *  All rights reserved.
 *
 */

#ifndef DORADE_H
#define DORADE_H

#include "read_dorade.h"
#include "datetime.h"
#include <cmath>
#include <string>

class Dorade
{
	
public:
	// Constructor / Destructor
	Dorade();
	Dorade(const std::string& swpFilename);
	~Dorade();
	
	// I/O
	bool readSwpfile();
	bool readSwpfile(const std::string& refname, const std::string& velname, const std::string& swname);
	bool writeSwpfile();
	bool writeSwpfile(const std::string& suffix);
	bool writeDoradefile(const std::string& doradeFilename);
	std::string getFilename();
	std::string getRadarname();
	void setFilename(const std::string& newname);
	
	// Data Retrieval
	float getAzimuth(int& ray);
	float getElevation(int& ray);
	float* getReflectivity(int& ray);
	float* getRadialVelocity(int& ray);	
	float* getSpectrumWidth(int& ray);
	float* getRayData(int &ray, const std::string& field);
	int getNumRays();
	int getNumGates();
	float* getGateSpacing();
	float getRadarLat();
	float getRadarLon();
	float getRadarAlt();
	float getRadarAltMSL();
	float getRadarLat(const int& ray);
	float getRadarLon(const int& ray);
	float getRadarAlt(const int& ray);
	float getRadarAltMSL(const int& ray);
	datetime getRayTime(int& ray);
	float getFLwind_u(const int& ray);
	float getFLwind_v(const int& ray);
	float getHeading(const int& ray);
	float getBeamwidthDeg();
	
	// Editing
	bool copyField(const std::string& oldFieldName,const std::string& newFieldName, 
				   const std::string& newFieldDesc,const std::string& newFieldUnits);
	//bool deleteField(const std::string& fldname);
	
private:
	bool swap_bytes;
	struct sswb_info *ssptr;
	struct vold_info *vptr;
	struct radd_info *rptr;
	struct celv_info *cptr;
	struct cfac_info *cfptr;
	struct parm_info *pptr;
	struct swib_info *sptr;
	struct ryib_info *ryptr;
	struct asib_info *aptr; 
	struct rdat_info **dptr;
	struct radar_angles *rangles;
	struct rktb_info *rkptr;
	struct rot_table_entry *rtptr;
	unsigned int rktb_size;
	int volumeTime;
	short int volumeDate;
	void swapVolHeader();
	void swapRadarHeader();
	bool machineBigEndian();
	int short swap2(char *ov);
	int swap4(char *ov);
	int refIndex;
	int velIndex;
	int swIndex;
	std::string ref_fld;
	std::string vel_fld;
	std::string sw_fld;	
	std::string filename;
	/*double isnanf(double x)   {  return  (((*(int *)&(x) & 0x7f800000L) == 0x7f800000L) && \
										  ((*(int *)&(x) & 0x007fffffL) != 0x00000000L)); }; */
	double RADIANS(double x)  { return ((x)*0.017453292); };
	double FMOD360(double x)  { return (fmod((double)((x)+720.), (double)360.)); };
    double DEGREES(double x)  { return ((x)*57.29577951); };
	
	/* PROTOTYPES */
	void sweepread(const char swp_fname[],struct sswb_info *ssptr, struct vold_info *vptr,
				   struct radd_info *rptr,struct celv_info *cptr,
				   struct cfac_info *cfptr,struct parm_info *pptr,
				   struct swib_info *sptr,struct ryib_info *ryptr,
				   struct asib_info *aptr, struct rdat_info **dptr);
	void sweepwrite(const char swp_fname[],struct sswb_info *ssptr, struct vold_info *vptr,
				   struct radd_info *rptr,struct celv_info *cptr,
				   struct cfac_info *cfptr,struct parm_info *pptr,
				   struct swib_info *sptr,struct ryib_info *ryptr,
				   struct asib_info *aptr, struct rdat_info **dptr, int doradeflag);
	void read_sswb(FILE *fp,struct sswb_info *ssptr);
	void read_vold(FILE *fp,struct vold_info *vptr);
	void read_radd(FILE *fp,struct radd_info *rptr);
	void read_cfac(FILE *fp,struct cfac_info *cptr);
	void read_parm(FILE *fp,struct parm_info *pptr);
	void read_celv(FILE *fp,struct celv_info *cptr,int desc_length);
	void read_swib(FILE *fp,struct swib_info *sptr);
	void read_ryib(FILE *fp,struct ryib_info *rptr);
	void read_asib(FILE *fp,struct asib_info *aptr);
	void read_rktb(FILE *fp);
	void read_rdat(FILE *fp,int fld_num,
				   int desc_len,int *match,short parm_type,
				   int beam_count,int total_gates,
				   short compression,int baddata_flag,
				   float scale_fac, struct rdat_info *rdat);
	//void read_ch_arr(FILE *fp,int arrsize);
	void read_sh_arr(FILE *fp,int arrsize,int beam_count,
					 int total_gates,short compression,
					 int baddata_flag,float scale_fac,
					 struct rdat_info *rdat);
	//void read_lg_arr(FILE *fp,int arrsize);
	//void read_fl_arr(FILE *fp,int arrsize);
	void get_field(struct parm_info parm[],int num_desc,int *fld_num);
	int read_int(FILE *fp);
	void skip_bytes(FILE *fp,int numskip);
	int dd_hrd16_uncompressx(short *ss,short *dd, int flag,
							 int *empty_run,int wmax ,int beam_count);
	int dd_compress(unsigned short *src, unsigned short *dst, unsigned short flag, int n );
	void calcAirborneAngles(struct asib_info *asib, struct cfac_info *cfac, struct ryib_info* ra,
							struct radar_angles *angles);
};

#endif
