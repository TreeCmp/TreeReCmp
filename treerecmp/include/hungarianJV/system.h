/* THIS CODE IS SOURCED FROM http://www.assignmentproblems.com/LAPJV.htm*/
// System dependent definitions and functions
// File: system.h

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>     // for seconds()
namespace tools {

/*************** FUNCTIONS  *******************/

extern void seedRandom(unsigned int seed);
	      
extern double randomX(void);

extern double seconds();


} // end of namespace