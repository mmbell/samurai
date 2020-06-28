#include "datetime.h"
#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>
#include <time.h>
#include <locale>
#include <iomanip>

datetime  ParseDate(std::string in, const char *fmt) {
  return ParseDate(in.c_str(), fmt);
}

datetime  ParseDate(const char *string, const char *fmt) {
   char buffer[512];
  std::tm tm = {};
  std::stringstream ss(string);
  //ss >> std::get_time(&tm, fmt);
  if (strptime(string, fmt, &tm) == NULL) {
      std::cout << "ERROR PARSING TIME; quitting." <<  " ( string = " << string << ", fmt = " << fmt << " ) " << std::endl;
      exit(1);
  }
//  std::cout << "------- ParseDate : " << string << " & " << fmt << " ---------- " << std::endl;
//  std::cout << "DEBUG: tm.tm_hour = " << tm.tm_hour << std::endl;
//  std::cout << "DEBUG: tm.tm_min  = " << tm.tm_min << std::endl;
//  std::cout << "DEBUG: tm.tm_sec  = " << tm.tm_sec << std::endl;
  std::time_t tt = timegm(&tm);
  //std::cout << "DEBUG: tt = " << tt << " put_time : " << std::put_time(&tm, "%c %Z") << " gmtime : " << std::gmtime(&tt) << std::endl;
  strftime(buffer, 512, "%c %Z", &tm);
//  std::cout << "DEBUG: tt = " << tt << " put_time : " << buffer << " gmtime : " << std::gmtime(&tt) << std::endl;
  std::gmtime(&tt);
 // auto tt2 = std::mktime(&tt);
  //auto tp = std::chrono::system_clock::from_time_t(std::mktime(&tm));
  auto tp = std::chrono::system_clock::from_time_t(tt);
  { 
     using namespace std::chrono;
     using namespace date;
//     std::cout << "ParseDate tp<3> = " << tp <<  "   ( In: " << string << " )" << std::endl;
  }
  return tp;
}


datetime  ParseTime(const char *string, const char *fmt) {
   std::string epochstring = "19700101" + std::string(string);
   std::string epochfmt = "%Y%m%d" + std::string(fmt);
  char buffer[512];
  std::tm tm = {};
  std::stringstream ss(epochstring);
  if (strptime(epochstring.c_str(), epochfmt.c_str(), &tm) == NULL) {
      std::cout << "ERROR PARSING TIME; quitting." <<  " ( epochstring = " << epochstring << ", epochfmt = " << epochfmt << " ) " << std::endl;
      exit(1);
  }
  std::time_t tt = timegm(&tm);
  strftime(buffer, 512, "%c %Z", &tm);
  std::gmtime(&tt);
  auto tp = std::chrono::system_clock::from_time_t(tt);
  { 
     using namespace std::chrono;
     using namespace date;
//     std::cout << "ParseDate tp<3> = " << tp <<  "   ( In: " << string << " )" << std::endl;
  }
  return tp;
}

int64_t Time(datetime in) {
  datetime start_of_day = date::floor<date::days>(in);
  auto diff = in - start_of_day;
  return std::chrono::duration_cast<std::chrono::seconds>(diff).count();
}

int64_t Date(datetime in) {
   auto tse = in.time_since_epoch();
  return std::chrono::duration_cast<std::chrono::seconds>(tse).count();
}



std::string PrintTime(datetime in) {
   using namespace std::chrono;
   using namespace date;

    datetime time = date::floor<date::days>(in);

    auto diff = in - time;
    //int seconds = std::chrono::duration_cast<std::chrono::seconds>(diff).count();

    std::time_t tt = std::chrono::system_clock::to_time_t(diff + std::chrono::time_point<std::chrono::system_clock>{});
    //std::cout << "DEBUG(PrintTime) tt = " << tt << std::endl;

    auto tm = std::gmtime(&tt);

   std::stringstream ss;
   char buffer[512];
   strftime(buffer, 512, "%H:%M:%S", tm);
   //ss << std::put_time(tm, "%H:%M:%S");
   std::string buf(buffer);
   //buf = buffer;
   return buf;
}

std::string PrintDate(datetime in) {
   using namespace std::chrono;
   using namespace date;

//   std::cout << "DEBUG: PrintTime => " << in << std::endl;

   std::stringstream ss;
   ss << in;
   std::string buf;
   buf = ss.str();
   return buf;
}

