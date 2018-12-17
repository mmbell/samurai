#include <iostream>
#include "timers.h"

// -------------------------- Timer Helpers -----------------------------

namespace Timers {
  // These come from the original RadarWind code

  void addDeltaTime(struct timeval * ptva,
		    double * psum)
  {
    struct timeval tvb;
    if (gettimeofday( &tvb, NULL) != 0)
      return;
    double deltaSec = tvb.tv_sec - ptva->tv_sec
      + 1.e-6 * (tvb.tv_usec - ptva->tv_usec);
    ptva->tv_sec = tvb.tv_sec;
    ptva->tv_usec = tvb.tv_usec;
    if (psum != NULL) (*psum) += deltaSec;
  }

  void printRunTime(const std::string &str,
			       struct timeval * ptva)
  {
    struct timeval tvb;
    if (gettimeofday( &tvb, NULL) != 0)
      return;
    double deltaSec = tvb.tv_sec - ptva->tv_sec
      + 1.e-6 * (tvb.tv_usec - ptva->tv_usec);
    std::cout << "runTime: " << str << ": " << deltaSec << std::endl;
    ptva->tv_sec = tvb.tv_sec;
    ptva->tv_usec = tvb.tv_usec;
  }

}
