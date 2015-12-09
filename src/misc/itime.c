#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
void
itime_(
 iarray )
  int iarray[];
{
  struct tm* tm;		/* date/time structure */
  time_t clock;			/* current time in seconds  */

  clock = time(0);
  tm = localtime(&clock);
  iarray[0] = tm->tm_hour;	/* hour */
  iarray[1] = tm->tm_min;	/* minute */
  iarray[2] = tm->tm_sec;	/* second */
}
