#ifndef DMLC_TIMER_H_
#define DMLC_TIMER_H_
/*!
 * \file timer.h
 * \brief cross platform timer for timing
 * \author Tianqi Chen
 */
#include <time.h>
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif
//#include "./logging.h"

namespace dmlc {
/*!
 * \brief return time in seconds
 */
inline double GetTime(void) {
  // TODO: use c++11 chrono when c++11 was available
  #ifdef __MACH__
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  CHECK(clock_get_time(cclock, &mts) == 0) << "failed to get time";
  mach_port_deallocate(mach_task_self(), cclock);
  return static_cast<double>(mts.tv_sec) + static_cast<double>(mts.tv_nsec) * 1e-9;
  #else
  #if defined(__unix__) || defined(__linux__)
  timespec ts;
  CHECK(clock_gettime(CLOCK_REALTIME, &ts) == 0) << "failed to get time";
  return static_cast<double>(ts.tv_sec) + static_cast<double>(ts.tv_nsec) * 1e-9;
  #else
  return static_cast<double>(time(NULL));
  #endif
  #endif
}
}  // namespace dmlc
#endif