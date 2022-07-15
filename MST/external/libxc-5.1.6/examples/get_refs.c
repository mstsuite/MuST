#include <stdio.h>

#include <xc.h>

int main()
{
  xc_func_type func;

  xc_func_init(&func, XC_GGA_X_B88, XC_UNPOLARIZED);

  printf("'%s', defined in the reference(s):\n", func.info->name);
  for(int i=0;i<5;i++)
    if(func.info->refs[i]!=NULL)
      printf("%s (DOI %s)\n",func.info->refs[i]->ref,func.info->refs[i]->doi);

  xc_func_end(&func);
}
