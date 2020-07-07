#ifndef SAMURAI_INCLUDE_DATETIME_H_
#define SAMURAI_INCLUDE_DATETIME_H_

#include <sstream>
#include "date/tz.h"
#include <iomanip>
#include <ctime>
#include <iostream>
//#include <ostream>
//#include <istream>

// Using std::chrono::system_clock, assuming time zones aren't used:
using datetime = std::chrono::time_point<std::chrono::system_clock>;

datetime  ParseDate(const char *string, const char *fmt);
datetime  ParseDate(std::string string_in, const char *fmt);
datetime  ParseTime(const char *string, const char *fmt);
datetime  ParseTime(std::string string_in, const char *fmt);
//std::string TimeString(datetime& in);
std::string PrintTime(datetime in);
std::string PrintDate(datetime in);

int64_t Time(datetime in);
int64_t Date(datetime in);



#endif // SAMURAI_INCLUDE_DATETIME_H_
