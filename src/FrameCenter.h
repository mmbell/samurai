/*
 *  FrameCenter.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef FRAMECENTER
#define FRAMECENTER

#include "datetime.h"

class FrameCenter
{
		
public:
	FrameCenter();
	FrameCenter(datetime& time, float lat, float lon, float u, float v);
	~FrameCenter();
	
	void printString();
	
	float getLat() const;
	void setLat(const float& lat);
	
	float getLon() const;
	void setLon(const float& lon);
	
	datetime getTime() const;
	void setTime(const datetime& time);
	
	float getUmean() const;
	void setUmean(const float& u);
	
	float getVmean() const;
	void setVmean(const float& v);
	
private:
	
	float latitude;
	float longitude;
	datetime time;
	float Um;
	float Vm;		
	
};

#endif

