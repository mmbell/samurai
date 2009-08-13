/*
 *  TCcenter.h
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#ifndef TCCENTER
#define TCCENTER

#include <QDateTime>

class TCcenter
{
		
public:
	TCcenter();
	TCcenter(QDateTime& time, float& lat, float& lon, float& u, float& v);
	~TCcenter();
	
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

