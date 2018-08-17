#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
#include <string>

namespace Timers {
  void addDeltaTime(struct timeval * ptva, double * psum);
  void printRunTime(const std::string &str, struct timeval * ptva);
}

#define START_TIMER(s) struct timeval timer_##s; Timers::addDeltaTime(&timer_##s, NULL);
#define PRINT_TIMER(msg, s) Timers::printRunTime(msg, &timer_##s);

#endif
