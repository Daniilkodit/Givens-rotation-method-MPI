#include "header.h"

double
  get_full_time ()
{
  struct timeval t;
  gettimeofday (&t, NULL);
  return t.tv_sec + t.tv_usec / 1e6;
}

double
  get_cpu_time ()
{
  struct rusage t;
  getrusage (RUSAGE_THREAD, &t);
  return t.ru_utime.tv_sec + t.ru_utime.tv_usec / 1e6;
}
