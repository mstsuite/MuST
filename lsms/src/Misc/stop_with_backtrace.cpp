#ifdef HAS_BACKTRACE
#include <execinfo.h>
#endif
#include <stdio.h>
#include <stdlib.h>

extern "C"
{
#ifdef HAS_BACKTRACE
void stop_with_backtrace_(const char *name, int len)
{
  void *array[20];
  size_t num;
  char **strings;

  printf("STOP called in %.*s\n",len,name);

  num = backtrace(array,10);
  strings = backtrace_symbols(array,num);

  printf("Obtained %zd stack frames.\n",num);
  for(size_t i=0; i<num; i++)
    printf("  %s\n",strings[i]);

  free(strings);

  exit(1);
}
#else
void stop_with_backtrace_(const char *name, int len)
{ 
  printf("STOP called in %.*s\n",len,name);
  printf("backtrace not available.\n");
  exit(1);
}
#endif

}
