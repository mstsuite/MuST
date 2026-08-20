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
   n = n + 1; Options(1,n)='-cbm';  Options(2,n)='--cpu-bigmat'; Options(3,n)='--cpu_bigmat'
              OptionKeys(n)='KKR Matrix Calculation on CPU'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-agm';  Options(2,n)='--acc-gijmat'; Options(3,n)='--acc_gijmat'
              OptionKeys(n)='Accelerate Gij Matrix Calculation'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-cgm';  Options(2,n)='--cpu-gijmat'; Options(3,n)='--cpu_gijmat'
              OptionKeys(n)='Gij Matrix Calculation on CPU'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-tl';  Options(2,n)='--timing-lsms'; Options(3,n)='--timing_lsms'
              OptionKeys(n)='Timing the LSMS Matrix Calculations'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-nbc';  Options(2,n)='--no-blocking-constr'; Options(3,n)='--no_blocking_constr'
              OptionKeys(n)='Constructing the LSMS Matrix without using blocks'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-pb';  Options(2,n)='--print-blocking'; Options(3,n)='--print_blocking'
              OptionKeys(n)='Print Blocking Details in LSMS Matrix Inverse'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-nbi';  Options(2,n)='--no-blocking-inverse'; Options(3,n)='--no_blocking_inverse'
              OptionKeys(n)='The LSMS Matrix Inverse is Performed without Blocking'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-re';  Options(2,n)='--report-energy'; Options(3,n)='--report_energy'
              OptionKeys(n)='Report GPU Energy Consuption'
              OptionValues(n) = '-1'
!
   n = n + 1; Options(1,n)='-ptm';  Options(2,n)='--print-tau'; Options(3,n)='--print_tau'
              OptionKeys(n)='Print the Tau Matrix'
              OptionValues(n) = '0'
!
!  ===================================================================
!  Emulation precision parameter and contour zone parameters are defined as follows:
!  (Assuming N = number of energy contour zones)
!      precision param 1 is applied to energy contour zone 1, for which Re[energy] starts from the bottom
!                           of the valence band up to energy contour zone param 1
!      precision param 2 is applied to energy contour zone 2, which is right next to zone 1 and ends where
!                           Im[energy] = energy contour zone param 2
!      precision param 3 is applied to energy contour zone 3, which is right next to zone 2 and ends where
!                           Im[energy] = energy contour zone param 3
!      ...,...
!      precision param N is applied to energy contour zone N, which is right next to zone N-1 and ends where
!                           Im[energy] = energy contour zone param N = 0.0.
!  The constrain: energy contour zone param 2 > energy contour zone param 3 > ... 
!                                             > energy contour zone param N = 0.0
!  ===================================================================

   n = n + 1; Options(1,n)='-em';  Options(2,n)='--emul-scheme'; Options(3,n)='--emul_scheme'
              OptionKeys(n)='Emulation Scheme'
              OptionValues(n) = '0'  ! other possible values are: 1 for Ozaki 1, 2 for Ozaki 2
!
   n = n + 1; Options(1,n)='-emp';  Options(2,n)='--emul-param'; Options(3,n)='--emul_param'
              OptionKeys(n)='Emulation Precision Parameters'
              OptionValues(n) = '12,14,16'
!
   n = n + 1; Options(1,n)='-cz';  Options(2,n)='--contour-zone'; Options(3,n)='--contour_zone'
              OptionKeys(n)='Energy Contour Zone Parameters'
              OptionValues(n) = '0.000,0.010'
