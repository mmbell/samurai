/*
 *  TCcenter.cpp
 *  samurai
 *
 *  Copyright 2008 Michael Bell. All rights reserved.
 *
 */

#include "TCcenter.h"

TCcenter::TCcenter()
{
	latitude = -999;
	longitude = -999;
	Um = -999;
	Vm = -999;
	time = QDateTime();
}

TCcenter::TCcenter(QDateTime& t, float& lat, float& lon, float& u, float& v)
{
	time = t;
	latitude = lat;
	longitude = lon;
	Um = u;
	Vm = v;
}

TCcenter::~TCcenter()
{
}

float TCcenter::getLat() const
{
	return latitude;
}

void TCcenter::setLat(const float& lat)
{
	latitude = lat;
}


float TCcenter::getLon() const
{
	return longitude;
}

void TCcenter::setLon(const float& lon)
{
	longitude = lon;
}

QDateTime TCcenter::getTime() const
{
	return time;
}

void TCcenter::setTime(const QDateTime& obTime)
{
	time = QDateTime(obTime);
}

float  TCcenter::getUmean() const
{
	return Um;
}

void  TCcenter::setUmean(const float& u)
{
	Um = u;
}

float  TCcenter::getVmean() const
{
	return Vm;
}

void  TCcenter::setVmean(const float& v)
{
	Vm = v;
}
