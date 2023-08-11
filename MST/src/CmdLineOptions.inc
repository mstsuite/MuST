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
