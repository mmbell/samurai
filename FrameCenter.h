/*
 *  FrameCenter.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef FRAMECENTER
#define FRAMECENTER

#include <QDateTime>

class FrameCenter
{
		
public:
	FrameCenter();
	FrameCenter(QDateTime& time, float& lat, float& lon, float& u, float& v);
	~FrameCenter();
	
	void printString();
	
	float getLat() const;
	void setLat(const float& lat);
	
	float getLon() const;
	void setLon(const float& lon);
	
	QDateTime getTime() const;
	void setTime(const QDateTime& time);
	
	float getUmean() const;
	void setUmean(const float& u);
	
	float getVmean() const;
	void setVmean(const float& v);
	
private:
	
	float latitude;
	float longitude;
	QDateTime time;
	float Um;
	float Vm;		
	
};

#endif

