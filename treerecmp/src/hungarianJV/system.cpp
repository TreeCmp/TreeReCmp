
/* THIS CODE IS SOURCED FROM http://www.assignmentproblems.com/LAPJV.htm*/

// System dependent routines
// File: system.cpp

#include "hungarianJV/system.h"
#include <stdlib.h>

namespace tools {

void seedRandom(unsigned int seed)
// seed for random number generator.
{
  srand(seed);
  return;   
}
	      
double randomX(void)
// random number between 0.0 and 1.0 (uncluded).
{
  double rrr;
  
  rrr = (double) rand() / (double) RAND_MAX;
  return rrr;
}
 
double seconds()
// cpu time in seconds since start of run.
{
  double secs;
   
  secs = (double)(clock() / 1000.0);
  return(secs);
}


} // end of namespace