#ifndef LSMS_ACCELERATOR_HPP
#define LSMS_ACCELERATOR_HPP

void acceleratorInitialize(int sz, int nthreads);
void acceleratorFinalize(void);
void acceleratorPrint(void);

extern "C"
{
void accelerator_initialize_(int *sz);
void accelerator_finalize_(void);
}

#endif
