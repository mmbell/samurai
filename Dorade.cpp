/*
 *  Dorade.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell
 *  All rights reserved.
 *
 */

#include "Dorade.h"

Dorade::Dorade(const QString& swpFilename) 
{
	swap_bytes = false;
	vptr = new vold_info;
	rptr = new radd_info;
	cptr = new celv_info;
	cfptr = new cfac_info;
	pptr = new parm_info[20];
	sptr = new swib_info;
	ryptr = new ryib_info[500];
	aptr = new asib_info[500];
	dptr = new rdat_info*[20];
	for (int i=0; i<20; i++) {
		dptr[i] = new rdat_info[500];
	}
	// struct rdat_info *dptr = new rdat_info[500];
	
	
	filename = swpFilename;
}

Dorade::~Dorade()
{
	
	delete vptr;
	delete rptr;
	delete cptr;
	delete cfptr;
	delete[] pptr;
	delete sptr;
	delete[] ryptr;
	delete[] aptr;
	for (int i=0; i<20; i++) {
		delete[] dptr[i];
	}
	delete[] dptr;
  	
}

bool Dorade::readSwpfile()
{

	// Check the byte order
	if (!machineBigEndian()) {
		swap_bytes = true;
	}
	
	// Read in a dorade file
	QByteArray bytes  = filename.toAscii();
	const char* ccfilename = bytes.data();
	sweepread(ccfilename, vptr, rptr, cptr,
			  cfptr, pptr, sptr, ryptr, aptr, dptr);
	
	/* ref_fld = "DZ";
	vel_fld = "VG";
	sw_fld = "SW"; */
	
	return true;
}

float Dorade::getAzimuth(int& ray)
{
	if (ray<(sptr->num_rays)) {
		return ryptr[ray].azimuth;
	} else {
		return -999;
	}
}

float Dorade::getElevation(int& ray)
{
	if (ray<(sptr->num_rays)) {
		return ryptr[ray].elevation;
	} else {
		return -999;
	}
	
}

float* Dorade::getReflectivity(int& ray)
{
	if (ray<(sptr->num_rays)) {
		return dptr[refIndex][ray].data;
	} else {
		return 0;
	}
}

float* Dorade::getRadialVelocity(int& ray)
{
	if (ray<(sptr->num_rays)) {
		return dptr[velIndex][ray].data ;
	} else {
		return 0;
	}
}

float* Dorade::getSpectrumWidth(int &ray)
{
	if (ray<(sptr->num_rays)) {
		return dptr[swIndex][ray].data;
	} else {
		return 0;
	}
}

int Dorade::getNumRays()
{
	return sptr->num_rays;
}

int Dorade::getNumGates()
{
	return cptr->total_gates;
}

float* Dorade::getGateSpacing()
{
	return cptr->gate_spacing;
}

float Dorade::getRadarAlt()
{
	return aptr->alt_agl;
}

float Dorade::getRadarLat()
{
	return aptr->lat;
}

float Dorade::getRadarLon()
{
	return aptr->lon;
}

QDateTime Dorade::getRayTime(int& ray)
{
	int year = vptr->year;
	int jd = ryptr[ray].julian_day;
	//int month = vptr->mon;
	//int day = vptr->day;
	QDate date(year,1,1);
	date = date.addDays(jd-1);
	int hour = ryptr[ray].hour;
	int min = ryptr[ray].min;
	int sec = ryptr[ray].sec;
	int msec = ryptr[ray].msec;
	//QDate date = QDate(year, month, day);
	QTime time = QTime(hour, min, sec, msec);
	return QDateTime(date, time, Qt::UTC);
}
	
bool Dorade::machineBigEndian(){
	
    union {
		unsigned char byte[4];
		int val;
    }word;
	
    word.val = 0;
    word.byte[3] = 0x01;
    
    return word.val == 1;
}

int short Dorade::swap2(char *ov)		/* swap integer*2 */
{
    union {
		int short newval;
		char nv[2];
    }u;
    u.nv[1] = *ov++; u.nv[0] = *ov++;
    return(u.newval);
}

int long Dorade::swap4(char *ov )		/* swap integer*4 */
{
	union {
		int long newval;
		char nv[4];
	}u;
	
	u.nv[3] = *ov++; u.nv[2] = *ov++; u.nv[1] = *ov++; u.nv[0] = *ov++;
	
	return(u.newval);
}
/**************************************************/
void Dorade::sweepread(const char swp_fname[],struct vold_info *vptr,
			   struct radd_info *rptr,struct celv_info *cptr,
			   struct cfac_info *cfptr,struct parm_info *pptr,
			   struct swib_info *sptr,struct ryib_info *ryptr,
			   struct asib_info *aptr, struct rdat_info **dptr)
{
	
	/****************************************************/
	/* THIS SUBROUTINE READS THE DESCRIPTORS FROM A     */
	/* SWEEP FILE & PASSES THE VOLD, RADD, CFAC, PARM,  */
	/* SWIB BACK TO THE FORTRAN PROGRAM                 */
	/* swp_fname=name of sweep file                     */
	/* *vptr=pointer to vold descriptor                 */
	/* *rptr=pointer to radd descriptor                 */
	/* *cptr=pointer to cfac descriptor                 */
	/* *pptr=pointer to parm descriptor                 */
	/* *sptr=pointer to swib descriptor                 */
	/* *ryptr=pointer to ryib descriptor                */
	/* *aptr=pointer to asib descriptor                 */
	/* *dptr=pointer to rdat descriptor                 */
	/* tot_gates=total number of gates                  */
	/* arr=array to hold gate spacing                   */
	/* *sptr=pointer of swib descriptor                 */
	/****************************************************/
	
	FILE *fp;
	int i=0;
	char identifier[IDENT_LEN];
	int desc_len;
	int fld_num=0;
	int match;
	int found=FALSE;
	int beam=-1;
	i = 0;
	/* ADD NULL CHARACTER TO END OF STRING
	for (unsigned int i=0;i<strlen(swp_fname);i++) {
		if (isspace(swp_fname[i])) {break;}
		else {len++;}
	}
	swp_fname[len]='\0'; */
	
	/* READ THE SWEEP FILE */
	/* OPEN THE SWEEP FILE */
	if ( (fp = fopen(swp_fname,"rb"))==NULL) {
		printf("Can't open %s\n",swp_fname);
		return;
	}
	
	while ( !feof(fp) ) {
		
		/* READ THE DESCRIPTOR IDENTIFIER */
		if ( (fread(identifier,sizeof(char),IDENT_LEN,fp)) != IDENT_LEN) {
			printf("sweep file read error..can't read identifier\n");
			exit(-1);
		}
		/*printf ("reading %s\n",identifier);*/
		
		/* READ THE DESCRIPTOR LENGTH */
		desc_len=read_long(fp);
		
		if ( (strncmp(identifier,"VOLD",IDENT_LEN)) == 0) {
			/* READ THE VOLUME DESCRIPTOR */
			/*printf ("reading vold\n");*/
			read_vold(fp,vptr);
			
		} else if ( (strncmp(identifier,"RADD",IDENT_LEN)) == 0) {
			/* READ THE RADAR DESCRIPTOR */
			/*printf ("reading radd\n");*/
			read_radd(fp,rptr);
			if (rptr->num_param_desc > MAX_NUM_PARMS) {
				printf ("WARNING: NUMBER OF PARAMETERS ");
				printf ("GREATER THAN MAX_NUM_PARMS\n");
				printf ("NUMBER OF PARAMETERS: %d\n",rptr->num_param_desc);
				printf ("MAX_NUM_PARMS: %d\n",MAX_NUM_PARMS);
				printf ("SEE READ_DORADE.H\n");
			}
			
			
		} else if ( (strncmp(identifier,"CFAC",IDENT_LEN)) == 0) {
			/* READ THE CFAC DESCRIPTOR */
			/*printf ("reading cfac\n");*/
			read_cfac(fp,cfptr);
			
		} else if ( (strncmp(identifier,"PARM",IDENT_LEN)) == 0) {
			/* READ THE PARAMETER DESCRIPTOR */
			/*printf ("reading parm\n");*/
			read_parm(fp,pptr);
			*pptr++;
			
		} else if ( (strncmp(identifier,"CELV",IDENT_LEN)) == 0) {
			/*printf ("reading CELV\n");*/
			/* CHECK & MAKE SURE NUMBER OF GATES DOESN'T
            EXCEED LIMIT */
			read_celv(fp,cptr,desc_len);
			if (cptr->total_gates>MAX_GATES) {
				printf ("WARNING: NUMBER OF GATES ");
				printf ("GREATER THAN MAX_GATES\n");
				printf ("SEE READ_DORADE.H\n");
				printf ("NUMBER OF GATES: %d\n",(int)cptr->total_gates);
				printf ("MAX_GATES: %d\n",MAX_GATES);
			}
			
		} else if ( (strncmp(identifier,"SWIB",IDENT_LEN)) == 0) {
			/* READ THE SWEEP INFO DESCRIPTOR */
			/*printf ("reading swib\n");*/
			read_swib(fp,sptr);
			if (sptr->num_rays > MAX_BEAMS) {
				printf ("WARNING: NUMBER OF BEAMS ");
				printf ("GREATER THAN MAX_BEAMS\n");
				printf ("SEE READ_DORADE.H\n");
				printf ("NUMBER OF BEAMS: %d\n",(int)sptr->num_rays);
				printf ("MAX_BEAMS: %d\n",MAX_BEAMS);
			}
			
			
		} else if ( (strncmp(identifier,"RYIB",IDENT_LEN)) == 0) {
			/* READ THE RAY INFO DESCRIPTOR */
			/*printf ("reading ryib\n");*/
			beam++;
			read_ryib(fp,&ryptr[beam]);
			fld_num=0;
			/* GO BACK TO FIRST PARAMETER DESCRIPTOR */
			for (i=0;i<rptr->num_param_desc;i++) {*pptr--;}
			/* DID WE FIND THE FIELD??
			if (beam==2 && found==FALSE) {
				printf ("%s NOT FOUND..\n",fld_name);
				printf ("VALID FIELD NAMES:\n");
				for (i=0;i<rptr->num_param_desc;i++) {
					printf ("%s\n",pptr->parm_name);
					*pptr++;
				}
				printf ("EXITING..\n");
				exit(-1);
			}
			*/
			
		} else if ( (strncmp(identifier,"ASIB",IDENT_LEN)) == 0) {
			/* READ THE PLATFORM INFO DESCRIPTOR */
			/*printf ("reading asib\n");*/
			read_asib(fp,&aptr[beam]);
			//*aptr++;
			
		} else if ( (strncmp(identifier,"RDAT",IDENT_LEN)) == 0) {
			/* READ THE DATA DESCRIPTOR */
			// printf ("reading rdat %s %d:\n", fld_name, desc_len);
			
			match=FALSE;
			read_rdat(fp,fld_num,desc_len,&match,
					  pptr->parm_type,beam,cptr->total_gates,
					  rptr->compress_flag,pptr->baddata_flag,
					  pptr->scale_fac,&dptr[fld_num][beam]);
			fld_num++;
			*pptr++;
						
			if (match==TRUE) {
				// *dptr++;
				found=match;
			}
			
		} else if ( (strncmp(identifier,"NULL",IDENT_LEN)) == 0) {
			break;
			
		} else {
			skip_bytes(fp,desc_len-(IDENT_LEN+sizeof(long)));
		} /* endif */

	} /* endwhile */

	fclose(fp);
	
	if (rptr->scan_mode == 9) {
		// Airborne data, need to calculate ground relative azimuth and elevation
		for (int i=0; i < sptr->num_rays; i++) {
			calcAirborneAngles(&aptr[i], cfptr, &ryptr[i]);
		}
	}
}
/**************************************************/
void Dorade::read_vold(FILE *fp,struct vold_info *vptr) 
{
	
	/* READ THE VOLUME DESCRIPTOR */
	if ( fread ((char *)vptr,sizeof (struct vold_info),1,fp) !=1 ) {
		puts("ERROR READING VOLUME DESCRIPTOR\n");
		exit(-1);
	} /* endif */

}
/**************************************************/
void Dorade::read_radd(FILE *fp,struct radd_info *rptr) 
{
	
	/* READ THE RADAR DESCRIPTOR */
	if ( fread ((char *)rptr,sizeof (struct radd_info),1,fp) !=1 ) {
		puts("ERROR READING RADAR DESCRIPTOR\n");
		exit(-1);
	} /* endif */

}
/**************************************************/
void Dorade::read_cfac(FILE *fp,struct cfac_info *cptr) 
{
	
	/* READ THE CFAC DESCRIPTOR */
	if ( fread ((char *)cptr,sizeof (struct cfac_info),1,fp) !=1 ) {
		puts("ERROR READING CFAC DESCRIPTOR\n");
		exit(-1);
	} /* endif */

}
/**************************************************/
void Dorade::read_parm(FILE *fp,struct parm_info *pptr) 
{
	/* READ THE PARMAMETER DESCRIPTOR */
	if ( fread ((char *)pptr,sizeof (struct parm_info),1,fp) !=1 ) {
		puts("ERROR READING PARAMETER DESCRIPTOR\n");
		exit(-1);
	} /* endif */

   pptr->parm_name[PARM_NAME_LEN-1]='\0';
   pptr->parm_desc[PARM_DESC_LEN-1]='\0';
   pptr->parm_unit[PARM_UNIT_LEN-1]='\0';

}
/***************************************************/
void Dorade::read_celv(FILE *fp,struct celv_info *cptr,int desc_len)
{
	
	int skip;
	
	/* TOTAL GATES */
	cptr->total_gates=read_long(fp);
	
	/* ALLOCATE THE ARRAY */
	/*
	 cptr->gate_spacing=calloc(cptr->total_gates,sizeof(float));
	 if (!cptr->gate_spacing) {
		 printf ("doublelocation error..aborting..\n");
		 exit(1);
	 } 
	 */
	
	/* GATE SPACING */
	if ( (fread(cptr->gate_spacing,sizeof(float),cptr->total_gates,fp))
		 != (unsigned int)cptr->total_gates) {
		puts("ERROR READING CELV DESCRIPTOR\n");
	}
	
	skip=desc_len-(sizeof(float)*cptr->total_gates+12);
	skip_bytes(fp,skip);
	
}
/**************************************************/
void Dorade::read_swib(FILE *fp,struct swib_info *sptr) 
{
	
	/* READ THE SWEEP INFO DESCRIPTOR */
	if ( fread ((char *)sptr,sizeof (struct swib_info),1,fp) !=1 ) {
		puts("ERROR READING SWIB DESCRIPTOR\n");
		exit(-1);
	} /* endif */

}
/**************************************************/
void Dorade::read_ryib(FILE *fp,struct ryib_info *rptr) 
{
	
	/* READ THE RAY INFO DESCRIPTOR */
	if ( fread ((char *)rptr,sizeof (struct ryib_info),1,fp) !=1 ) {
		puts("ERROR READING RYIB DESCRIPTOR\n");
		exit(-1);
	} /* endif */

}
/**************************************************/
void Dorade::read_asib(FILE *fp,struct asib_info *aptr) 
{
	
	/* READ THE PLATFORM INFO DESCRIPTOR */
	if ( fread ((char *)aptr,sizeof (struct asib_info),1,fp) !=1 ) {
		puts("ERROR READING ASIB DESCRIPTOR\n");
		exit(-1);
	} /* endif */

}
/***************************************************/
void Dorade::read_rdat(FILE *fp,int fld_num,
               int desc_len,int *match,short parm_type,
			   int beam_count,long total_gates,
               short compression,long baddata_flag,
			   float scale_fac,struct rdat_info *rdat)
{
	
	/* fp=pointer to sweep file
	*  fld_name=user supplied field name
	*  fld_num=index of parameter descriptor
	*  desc_len=length of RDAT descriptor
	*  match=flag to indicate field match found
	*  parm_type=data type of parameter
	*  dptr=pointer to rdat structure
	*/ 
	
	int strsize,datasize,arrsize;
	strsize = datasize = arrsize = 0;
	char tempname[PARM_NAME_LEN];
	//QString ref_fld = "DBZ";
	QString ref_fld = "DZ";
	QString vel_fld = "VG";
	//QString sw_fld = "SW";
	QString sw_fld = "SPEC_WDT";
	memset(rdat->parm_name,' ',PARM_NAME_LEN);
	memset(tempname,' ',PARM_NAME_LEN);
	
	/* READ THE PARAMETER NAEM */
	if ( (fread(tempname,sizeof(char),PARM_NAME_LEN,fp))
         != PARM_NAME_LEN)
	{printf("sweep file read error..can't read parameter name\n");}
	
	/* CALCULATE LENGTH OF TEMPNAME */
	for (unsigned int strsize=0;strsize<strlen(tempname);strsize++) {
		if (isspace(tempname[strsize])) {break;}
	}
		
	/* FIND THE CORRECT FIELD */
	/* Modified to read all fields, but record ref, vel, and sw indices - MB */
	QString fld_name(tempname);
	if (fld_name.size() > 8) fld_name.resize(8);
	if (fld_name.trimmed() == ref_fld) {
		refIndex = fld_num;
	} else if (fld_name.trimmed() == vel_fld) {
		velIndex = fld_num;
	} else if (fld_name.trimmed() == sw_fld) {
		swIndex = fld_num;
	}
	
	*match=TRUE;
	/* CALCULATE SIZE OF DATA */
	strncpy(rdat->parm_name,tempname,strsize);
	if (parm_type==1) {datasize=sizeof(char);}
	else if (parm_type==2) {datasize=sizeof(short);}
	else if (parm_type==3) {datasize=sizeof(long);}
	else if (parm_type==4) {datasize=sizeof(float);}
	/* SIZE OF ARRAY */
	arrsize=(desc_len-(IDENT_LEN+sizeof(long)
					   +PARM_NAME_LEN))/datasize;
	/* READ IN THE DATA */
	if (datasize==1) {
		read_ch_arr(fp,arrsize);
	} else if (datasize==2) {
		read_sh_arr(fp,arrsize,beam_count,total_gates,compression,
					baddata_flag,scale_fac,rdat);
	} else if (datasize==3) {
		read_lg_arr(fp,arrsize);
	} else if (datasize==4) {
		read_fl_arr(fp,arrsize);
	} /* endif */

	//} else {
	//	skip_bytes(fp,desc_len-(IDENT_LEN+sizeof(long)+PARM_NAME_LEN));
	//}

}

/***************************************************/
void Dorade::read_sh_arr(FILE *fp,int arrsize,int beam_count,
				 long total_gates,short compression,
				 long baddata_flag,float scale_fac,
                 struct rdat_info *rdat) {
	
	static short arr_com[MAX_GATES], arr_uncom[MAX_GATES];
	int i,num;
	int empty_run=0;
	
	if (compression==TRUE) {
		/* READ A RAY OF DATA */
		if ( (fread(arr_com,sizeof(short),arrsize,fp)) != (unsigned int)arrsize)
		{printf("sweep file read error..can't read data\n");}
		
		/* UNCOMPRESS A RAY OF DATA */
		num=dd_hrd16_uncompressx(arr_com,arr_uncom,baddata_flag,
								 &empty_run,total_gates,beam_count);
		
	} else {
		/* READ A RAY OF DATA */
		if ( (fread(arr_uncom,sizeof(short),arrsize,fp)) != (unsigned int)arrsize)
		{printf("sweep file read error..can't read data\n");}
	} /* endif */

   /* SCALE THE DATA */
   for (i=0;i<total_gates;i++) {
	   if (arr_uncom[i] != baddata_flag) {
		   rdat->data[i]=(float)arr_uncom[i]/scale_fac;
	   } else {
		   rdat->data[i]=(float)arr_uncom[i];
	   }
   }

}

/* Not implemented */
void Dorade::read_ch_arr(FILE *fp,int arrsize) {
}
void Dorade::read_lg_arr(FILE *fp,int arrsize) {
}
void Dorade::read_fl_arr(FILE *fp,int arrsize) {
}

/***************************************************/
void Dorade::get_field(struct parm_info parm[],int num_desc,int *fld_num)
{	

	/* GET THE NAME OF THE DESIRED FIELD */
	printf ("Please choose number of the desired field:\n");
	
	for (int i=0;i<num_desc;i++) {
		printf ("%2d.  %s\n",i+1,parm[i].parm_name);
	}
	
	scanf("%d",fld_num);
	
}
/***************************************************/
long Dorade::read_long(FILE *fp)
{
	
	long temp;
	
	if ( (fread(&temp,sizeof(long),1,fp)) != 1) {
		printf("sweep file read error..\n");
	}
	return temp;
	
}
/***************************************************/
void Dorade::skip_bytes(FILE *fp,int numskip)
{
	
	/* SKIP TO THE RIGHT BYTE! */
	
	if (fseek(fp,numskip,SEEK_CUR)) {
		printf("Seek Error..aborting..\n");
		exit(1);
	}
	
}
/***************************************************/
/* Dick Oye's decompression routine */
int Dorade::dd_hrd16_uncompressx(short *ss,short * dd,
						 int flag,int *empty_run,
						 int wmax ,int beam_count)

//  short *ss, *dd;
//  int flag, *empty_run, wmax;
//  int beam_count;
{
    /*
     * routine to unpacks actual data assuming MIT/HRD compression where:
     * ss points to the first 16-bit run-length code for the compressed data
     * dd points to the destination for the unpacked data
     * flag is the missing data flag for this dataset that is inserted
     *     to fill runs of missing data.
     * empty_run pointer into which the number of missing 16-bit words
     *    can be stored. This will be 0 if the last run contained data.
# wmax indicate the maximum number of 16-bit words the routine can
     *    unpack. This should stop runaways.
     */
    int n, mark, wcount=0;
	
    while(*ss != 1) {           /* 1 is an end of compression flag */
        n = *ss & 0x7fff;       /* nab the 16-bit word count */
        if(wcount+n > wmax) {
            printf("Uncompress failure %d %d %d at %d\n"
                   , wcount, n, wmax, beam_count);
            mark = 0;
            break;
        }
        else {
            wcount += n;                /* keep a running tally */
        }
        if( *ss & 0x8000 ) {    /* high order bit set implies data! */
            *empty_run = 0;
            ss++;
            for(; n--;) {
                *dd++ = *ss++;
            }
        }
else {                  /* otherwise fill with flags */
*empty_run = n;
ss++;
for(; n--;) {
	*dd++ = flag;
}
}
    }
return(wcount);
}

void Dorade::calcAirborneAngles(struct asib_info *asib, struct cfac_info *cfac, struct ryib_info* ra)
{
    /* compute the true azimuth, elevation, etc. from platform
     * parameters using Testud's equations with their different
     * definitions of rotation angle, etc.
     */
    double R, H, P, D, T, theta_a, tau_a, lambda_t;
    double sinP, cosP, sinT, cosT, sinD, cosD;
    double sin_theta_rc, cos_theta_rc, sin_tau_a, cos_tau_a;
    double xsubt, ysubt, zsubt;
	
    double d;
 	/*
     * see Wen-Chau Lee's paper
     * "Mapping of the Airborne Doppler Radar Data"
     */
    d = asib->roll;
    R = isnanf(d) ? 0 : RADIANS(d +cfac->c_roll);
	
    d = asib->pitch;
    P = isnanf(d) ? 0 : RADIANS(d +cfac->c_pitch);
	
    d = asib->head;
    H = isnanf(d) ? 0 : RADIANS(d +cfac->c_head);
	
    d = asib->drift;
    D = isnanf(d) ? 0 : RADIANS(d +cfac->c_drift);
	
    sinP = sin(P);
    cosP = cos(P);
    sinD = sin(D);
    cosD = cos(D);
    
    T = H + D;
		
    sinT = sin(T);
    cosT = cos(T);
	
    d = asib->rot_ang;
    theta_a = isnanf(d) ? 0 : RADIANS(d +cfac->c_rotang);
    
    d = asib->tilt_ang;
    tau_a = isnanf(d) ? 0 : RADIANS(d +cfac->c_tiltang);
    sin_tau_a = sin(tau_a);
    cos_tau_a = cos(tau_a);
    sin_theta_rc = sin(theta_a + R); /* roll corrected rotation angle */
    cos_theta_rc = cos(theta_a + R); /* roll corrected rotation angle */
	
    xsubt = (cos_theta_rc * sinD * cos_tau_a * sinP
					 + cosD * sin_theta_rc * cos_tau_a
					 -sinD * cosP * sin_tau_a);
    ysubt = (-cos_theta_rc * cosD * cos_tau_a * sinP
					 + sinD * sin_theta_rc * cos_tau_a
					 + cosP * cosD * sin_tau_a);
    zsubt = (cosP * cos_tau_a * cos_theta_rc
					 + sinP * sin_tau_a);
	
    //ra->rotation_angle = theta_t = atan2(xsubt, zsubt);
    //ra->tilt = tau_t = asin(ysubt);
    lambda_t = atan2(xsubt, ysubt);
	double azrad = fmod(lambda_t + T, 6.283185307);
	double elrad = asin(zsubt);
    ra->azimuth = FMOD360(DEGREES(azrad)+360.);
    ra->elevation = DEGREES(elrad);    
	return;
	
}

