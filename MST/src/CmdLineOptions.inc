!  =============================================================================
!  Table of command line options and their corresponding keys
!  =============================================================================
   n = n + 1; Options(1,n)='-wp';   Options(2,n)='--write-procs';    Options(3,n)='--write_procs'
              OptionKeys(n)='Number of Writing Processes'
              OptionValues(n) = '0'
!
   n = n + 1; Options(1,n)='-rp';   Options(2,n)='--read-procs';     Options(3,n)='--read_procs'
              OptionKeys(n)='Number of Reading Processes'
              OptionValues(n) = '0'
!
   n = n + 1; Options(1,n)='-app';  Options(2,n)='--atoms-per-proc'; Options(3,n)='--atoms_per_proc'
              OptionKeys(n)='Number of Atoms Per Process'
              OptionValues(n) = '0'
!
   n = n + 1; Options(1,n)='-cpu';  Options(2,n)='--cpu-only'; Options(3,n)='--cpu_only'
              OptionKeys(n)='Run on CPU without Acceleration'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-abm';  Options(2,n)='--acc-bigmat'; Options(3,n)='--acc_bigmat'
              OptionKeys(n)='Accelerate KKR Matrix Calculation'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-tl';  Options(2,n)='--timing-lsms'; Options(3,n)='--timing_lsms'
              OptionKeys(n)='Timing the LSMS Matrix Calculations'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-pb';  Options(2,n)='--print-blocking'; Options(3,n)='--print_blocking'
              OptionKeys(n)='Print Blocking Details in LSMS Matrix Inverse'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-nb';  Options(2,n)='--no-blocking'; Options(3,n)='--no_blocking'
              OptionKeys(n)='The LSMS Matrix Inverse is Performed without Blocking'
              OptionValues(n) = '-1'
