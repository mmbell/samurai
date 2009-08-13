/*
 *  Dorade.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell
 *  All rights reserved.
 *
 */

#ifndef DORADE_H
#define DORADE_H


#include "read_dorade.h"
#include <QDir>
#include <QDateTime>
#include <cmath>
#include "precision.h"


class Dorade
{
	
public:
	Dorade(const QString& swpFilename);
	~Dorade();
	bool readSwpfile();
	float getAzimuth(int& ray);
	float getElevation(int& ray);
	float* getReflectivity(int& ray);
	float* getRadialVelocity(int& ray);	
	float* getSpectrumWidth(int& ray);
	int getNumRays();
	int getNumGates();
	float* getGateSpacing();
	float getRadarLat();
	float getRadarLon();
	float getRadarAlt();
	QDateTime getRayTime(int& ray);
	
private:
	bool swap_bytes;	
	struct vold_info *vptr;
	struct radd_info *rptr;
	struct celv_info *cptr;
	struct cfac_info *cfptr;
	struct parm_info *pptr;
	struct swib_info *sptr;
	struct ryib_info *ryptr;
	struct asib_info *aptr; 
	struct rdat_info **dptr;
	long int volumeTime;
	short int volumeDate;
	void swapVolHeader();
	void swapRadarHeader();
	bool machineBigEndian();
	int short swap2(char *ov);
	int long swap4(char *ov);
	int refIndex;
	int velIndex;
	int swIndex;
	QString filename;
	double isnanf(double x)   {  return  (((*(long *)&(x) & 0x7f800000L) == 0x7f800000L) && \
										  ((*(long *)&(x) & 0x007fffffL) != 0x00000000L)); };
	double RADIANS(double x)  { return ((x)*0.017453292); };
	double FMOD360(double x)  { return (fmod((double)((x)+720.), (double)360.)); };
    double DEGREES(double x)  { return ((x)*57.29577951); };
	
	/* PROTOTYPES */
	void sweepread(const char swp_fname[],struct vold_info *vptr,
				   struct radd_info *rptr,struct celv_info *cptr,
				   struct cfac_info *cfptr,struct parm_info *pptr,
				   struct swib_info *sptr,struct ryib_info *ryptr,
				   struct asib_info *aptr, struct rdat_info **dptr);
	void read_vold(FILE *fp,struct vold_info *vptr);
	void read_radd(FILE *fp,struct radd_info *rptr);
	void read_cfac(FILE *fp,struct cfac_info *cptr);
	void read_parm(FILE *fp,struct parm_info *pptr);
	void read_celv(FILE *fp,struct celv_info *cptr,int desc_length);
	void read_swib(FILE *fp,struct swib_info *sptr);
	void read_ryib(FILE *fp,struct ryib_info *rptr);
	void read_asib(FILE *fp,struct asib_info *aptr);
	void read_rdat(FILE *fp,int fld_num,
				   int desc_len,int *match,short parm_type,
				   int beam_count,long total_gates,
				   short compression,long baddata_flag,
				   float scale_fac, struct rdat_info *rdat);
	void read_ch_arr(FILE *fp,int arrsize);
	void read_sh_arr(FILE *fp,int arrsize,int beam_count,
					 long total_gates,short compression,
					 long baddata_flag,float scale_fac,
					 struct rdat_info *rdat);
	void read_lg_arr(FILE *fp,int arrsize);
	void read_fl_arr(FILE *fp,int arrsize);
	void get_field(struct parm_info parm[],int num_desc,int *fld_num);
	long read_long(FILE *fp);
	void skip_bytes(FILE *fp,int numskip);
	int dd_hrd16_uncompressx(short *ss,short *dd, int flag,
							 int *empty_run,int wmax ,int beam_count);
	void calcAirborneAngles(struct asib_info *asib, struct cfac_info *cfac, struct ryib_info* ra);
};

#endif
