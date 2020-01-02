#include <time.h>

double MPI_Wtime(void);

double clock_time_()
{
  return MPI_Wtime();
//    return(((double)clock())/((double)CLOCKS_PER_SEC));
}
