#!/usr/bin/env perl

# Copyright (C) 2017 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

use Data::Dumper;
use FindBin;                 # locate this script
use lib "$FindBin::Bin/.";   # use current directory

use maple2c_work;
use maple2c_derivatives;

die "Usage: $0 srcdir functional max_order (optional args) (simplify)\n"
    if ($#ARGV < 2);

$config{"srcdir"}     = $ARGV[0];
$config{"functional"} = $ARGV[1];
$config{"max_order"}  = $ARGV[2];
$config{"simplify"}   = ($#ARGV >= 3 && $ARGV[3] eq "yes") ? 1 : 0;
$config{"prefix"}     = "";

$config{"simplify_begin"} = ($config{'simplify'} == 1) ? "simplify(" : "";
$config{"simplify_end"}   = ($config{'simplify'} == 1) ? ", symbolic)" : "";

$config{"replace"} = [];

# find out where the maple file resides
if(-f $config{"functional"}.".mpl"){
  $config{"mathfile"} = $config{"functional"}.".mpl";
}elsif(-f $config{"srcdir"}."/maple/".$config{"functional"}.".mpl"){
  $config{"mathfile"} = $config{"srcdir"}."/maple/".$config{"functional"}.".mpl";
}else{
  $temp = $config{"functional"};
  $temp =~ s/hyb_//;
  $temp =~ s/_.*$//;

  if(-f $config{"srcdir"}."/maple/".$temp."_exc/".$config{"functional"}.".mpl"){
    $config{"mathfile"} =  $config{"srcdir"}."/maple/".$temp."_exc/".$config{"functional"}.".mpl";
  }elsif(-f $config{"srcdir"}."/maple/".$temp."_vxc/".$config{"functional"}.".mpl"){
    $config{"mathfile"} =  $config{"srcdir"}."/maple/".$temp."_vxc/".$config{"functional"}.".mpl";
  }
}

open my $in, '<', $config{"mathfile"} or die "File $mathfile does not exist\n";

# Find out the type of functional
while($_ = <$in>){
  if(/^\(\* type:\s(\S*)\s/){
    $config{"functype"} = $1;
  };
  if(/^\(\* replace:\s*"([^"]*)"\s*->\s*"([^"]*)"/){
    push @{$config{"replace"}}, "$1";
    push @{$config{"replace"}}, "$2";
  };
  if(/^\(\* prefix:/){
    while( ($_ = <$in>) && ! /^\*\)/ ){
      $config{"prefix"} .= $_;
    }
  }
}
close($in);

my %commands = (
  "lda_exc"   => \&work_lda_exc,
  "lda_vxc"   => \&work_lda_vxc,
  "gga_exc"   => \&work_gga_exc,
  "gga_vxc"   => \&work_gga_vxc,
  "mgga_exc"  => \&work_mgga_exc,
  "mgga_vxc"  => \&work_mgga_vxc,
    );

if ($commands{$config{"functype"}}) {
  $commands{$config{"functype"}}->();
} else {
  die "No such type: ".$config{"functype"}."\n";
}

#####################################################################
# Returns the variables of a LDA and constructs all spin variants of
# the derivatives
sub lda_var_der() {
  # these are the variables that the functional depends on
  my @variables = ("rho_0_", "rho_1_");

  # the definition of the derivatives that libxc transmits to the calling program
  my @partials = (
    ["zk"],
    ["vrho"],
    ["v2rho2"],
    ["v3rho3"],
    ["v4rho4"]
      );

  my @derivatives = ();
  for(my $order=0; $order< $config{"max_order"}+1; $order++){
    $derivatives[$order] = [];
    foreach my $der (@{$partials[$order]}){
      push @{$derivatives[$order]}, maple2c_derivatives($der, "lda");
    }
  }

  return (\@variables, \@derivatives);
}


#####################################################################
# Returns the variables of a GGA and constructs all spin variants of
# the derivatives
sub gga_var_der() {
  # these are the variables that the functional depends on
  my @variables = ("rho_0_", "rho_1_", "sigma_0_", "sigma_1_", "sigma_2_");

  # the definition of the derivatives that libxc transmits to the calling program
  my @partials = (
    ["zk"],
    ["vrho", "vsigma"],
    ["v2rho2", "v2rhosigma", "v2sigma2"],
    ["v3rho3", "v3rho2sigma", "v3rhosigma2", "v3sigma3"],
    ["v4rho4", "v4rho3sigma", "v4rho2sigma2", "v4rhosigma3", "v4sigma4"]
      );

  my @derivatives = ();
  for(my $order=0; $order< $config{"max_order"}+1; $order++){
    $derivatives[$order] = [];
    foreach my $der (@{$partials[$order]}){
      push @{$derivatives[$order]}, maple2c_derivatives($der, "gga");
    }
  }

  return (\@variables, \@derivatives);
}


#####################################################################
# Returns the variables of a MGGA and constructs all spin variants of
# the derivatives
sub mgga_var_der() {
  # these are the variables that the functional depends on
  my @variables = ("rho_0_", "rho_1_", "sigma_0_", "sigma_1_", "sigma_2_", "lapl_0_", "lapl_1_", "tau_0_", "tau_1_");

  # the definition of the derivatives that libxc transmits to the calling program
  my @partials = (
    ["zk"],
    ["vrho", "vsigma", "vlapl", "vtau"],
    ["v2rho2", "v2rhosigma", "v2rholapl", "v2rhotau", "v2sigma2",
     "v2sigmalapl", "v2sigmatau", "v2lapl2", "v2lapltau", "v2tau2"],
    ["v3rho3", "v3rho2sigma", "v3rho2lapl", "v3rho2tau", "v3rhosigma2",
     "v3rhosigmalapl", "v3rhosigmatau", "v3rholapl2", "v3rholapltau",
     "v3rhotau2", "v3sigma3", "v3sigma2lapl", "v3sigma2tau", "v3sigmalapl2",
     "v3sigmalapltau", "v3sigmatau2", "v3lapl3", "v3lapl2tau", "v3lapltau2",
     "v3tau3"],
    ["v4rho4", "v4rho3sigma", "v4rho3lapl", "v4rho3tau", "v4rho2sigma2",
     "v4rho2sigmalapl", "v4rho2sigmatau", "v4rho2lapl2", "v4rho2lapltau",
     "v4rho2tau2", "v4rhosigma3", "v4rhosigma2lapl", "v4rhosigma2tau",
     "v4rhosigmalapl2", "v4rhosigmalapltau", "v4rhosigmatau2",
     "v4rholapl3", "v4rholapl2tau", "v4rholapltau2", "v4rhotau3",
     "v4sigma4", "v4sigma3lapl", "v4sigma3tau", "v4sigma2lapl2",
     "v4sigma2lapltau", "v4sigma2tau2", "v4sigmalapl3", "v4sigmalapl2tau",
     "v4sigmalapltau2", "v4sigmatau3", "v4lapl4", "v4lapl3tau",
     "v4lapl2tau2", "v4lapltau3", "v4tau4",
    ]
      );

  my @derivatives = ();
  for(my $order=0; $order< $config{"max_order"}+1; $order++){
    $derivatives[$order] = [];
    foreach my $der (@{$partials[$order]}){
      push @{$derivatives[$order]}, maple2c_derivatives($der, "mgga");
    }
  }

  return (\@variables, \@derivatives);
}

#####################################################################
# This separates the derivatives of vxc into derivatives of vxc_0
# and vxc_1. All other derivatives (e.g. vsigma) are ignored.
sub filter_vxc_derivatives {
  my @all_derivatives = @{$_[0]};

  my @derivatives  = ();
  my @derivatives1 = ();
  my @derivatives2 = ();

  for(my $order=0; $order < $#all_derivatives; $order++){
    $derivatives [$order] = [];
    $derivatives1[$order] = [];
    $derivatives2[$order] = [];
    foreach my $der (@{$all_derivatives[$order + 1]}){
      if(${$der}[0][0] > 0){
        ${$der}[0][0]--;
        push @{$derivatives1[$order]}, $der;
        push @{$derivatives [$order]}, $der;
      }elsif(${$der}[0][1] > 0){
        ${$der}[0][1]--;
        push @{$derivatives2[$order]}, $der;
        push @{$derivatives [$order]}, $der;
      }
    }
  }

  return (\@derivatives, \@derivatives1, \@derivatives2);
}

#####################################################################
# sort by character and then by number
sub sort_alphanumerically {
  $a =~ /([^_]+)_([0-9]+)_/;
  my ($name1, $num1) = ($1, $2);

  $b =~ /([^_]+)_([0-9]+)_/;
  my ($name2, $num2) = ($1, $2);

  return $name1 cmp $name2 || $num1 <=> $num2;
}


#####################################################################
# Process a LDA functional for the energy
sub work_lda_exc {
  my $variables, $derivatives;

  ($variables, $derivatives) = lda_var_der();

  # get arguments of the functions
  $input_args  = "const double *rho";
  $output_args = ", double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, )";

  my ($der_def, @out_c) =
      maple2c_create_derivatives($variables, $derivatives, "mf");
  my $out_c = join(", ", @out_c);
  $out_c = ", $out_c" if ($out_c ne "");

  # we join all the pieces
  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1) -> \\
  $config{'simplify_begin'} \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1)) \\
  $config{'simplify_end'}:

(* mf is energy per unit volume *)
mf   := (r0, r1) -> eval(dens(r0, r1)*mzk(r0, r1)):

\$include <util.mpl>
";
  my $maple_zk = " zk_0_ = mzk(".join(", ", @{$variables}).")";

  # we build 2 variants of the functional, for unpolarized, and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
\n",
      );

  maple2c_run($variables, $derivatives, \@variants, 0, $input_args, $output_args);
}

#####################################################################
# Process a LDA functional for the potential
sub work_lda_vxc {
  my $variables, $all_derivatives;

  ($variables, $all_derivatives) = lda_var_der();
  my $derivatives, $derivatives1, $derivatives2;
  ($derivatives, $derivatives1, $derivatives2) = filter_vxc_derivatives($all_derivatives);

  # get arguments of the functions
  $input_args  = "const double *rho";
  $output_args = "LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, )";

  # we obtain the missing pieces for maple
  # unpolarized calculation
  my ($der_def_unpol, @out_c_unpol) =
      maple2c_create_derivatives($variables, $derivatives1, "mf0", "unpol");
  my $out_c_unpol = join(", ", @out_c_unpol);
  $out_c_unpol = ", ".$out_c_unpol if (! $out_c_unpol =~ /^ *$/);

  # polarized calculation
  my ($der_def_pol, @out_c_pol1) =
      maple2c_create_derivatives($variables, $derivatives1, "mf0", "pol");
  my ($der_def_pol2, @out_c_pol2) =
      maple2c_create_derivatives($variables, $derivatives2, "mf1", "pol");

  $der_def_pol .= $der_def_pol2;

  push(@out_c_pol1, @out_c_pol2);
  my $out_c_pol = join(", ", sort sort_alphanumerically @out_c_pol1);
  $out_c_pol = ", ".$out_c_pol if (! $out_c_pol =~ /^ *$/);

  # we join all the pieces
  my $maple_code1 = "
(* mf is the up potential *)
mzk   := (r0, r1) -> $config{'simplify_begin'} f(r_ws(dens(r0, r1)), zeta(r0, r1)) $config{'simplify_end'}:
mf0   := (r0, r1) -> eval(mzk(r0, r1)):
mf1   := (r0, r1) -> eval(mzk(r1, r0)):

\$include <util.mpl>
";

  my $maple_vrho0 = "vrho_0_ = mf0(".join(", ", @{$variables}).")";
  my $maple_vrho1 = "vrho_1_ = mf1(".join(", ", @{$variables}).")";

  # we build 2 variants of the functional, for unpolarized, and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:

$der_def_unpol

$maple_code1
C([$maple_vrho0$out_c_unpol], optimize, deducetypes=false):
",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):

$der_def_pol

$maple_code1
C([$maple_vrho0, $maple_vrho1$out_c_pol], optimize, deducetypes=false):
"
      );

  maple2c_run($variables, $derivatives, \@variants, 1, $input_args, $output_args);
}


#####################################################################
# Process a GGA functional for the energy
sub work_gga_exc {
  my $variables, $derivatives;

  ($variables, $derivatives) = gga_var_der();

  # get arguments of the functions
  $input_args  = "const double *rho, const double *sigma";
  $output_args = ", double *zk GGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, )";

  my ($der_def, @out_c) =
      maple2c_create_derivatives($variables, $derivatives, "mf");
  my $out_c = join(", ", @out_c);
  $out_c = ", $out_c" if ($out_c ne "");

  # we join all the pieces
  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1, s0, s1, s2) -> \\
  $config{'simplify_begin'} \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2)) \\
  $config{'simplify_end'}:

(* mf is energy per unit volume *)
mf   := (r0, r1, s0, s1, s2) -> eval(dens(r0, r1)*mzk(r0, r1, s0, s1, s2)):

\$include <util.mpl>
";
  my $maple_zk = "zk_0_ = mzk(".join(", ", @{$variables}).")";

  # we build 2 variants of the functional, for unpolarized and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
\n",
      );

  maple2c_run($variables, $derivatives, \@variants, 0, $input_args, $output_args);
}

#####################################################################
# Process a GGA functional for the potential
sub work_gga_vxc {
  my $variables, $all_derivatives;

  ($variables, $all_derivatives) = gga_var_der();
  my $derivatives, $derivatives1, $derivatives2;
  ($derivatives, $derivatives1, $derivatives2) = filter_vxc_derivatives($all_derivatives);

  # get arguments of the functions
  $input_args  = "const double *rho, const double *sigma";
  $output_args = "GGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, )";

  # we obtain the missing pieces for maple
  # unpolarized calculation
  my ($der_def_unpol, @out_c_unpol) =
      maple2c_create_derivatives($variables, $derivatives1, "mf0", "unpol");
  my $out_c_unpol = join(", ", @out_c_unpol);
  $out_c_unpol = ", ".$out_c_unpol if (! $out_c_unpol =~ /^ *$/);

  # polarized calculation
  my ($der_def_pol, @out_c_pol1) =
      maple2c_create_derivatives($variables, $derivatives1, "mf0", "pol");
  my ($der_def_pol2, @out_c_pol2) =
      maple2c_create_derivatives($variables, $derivatives2, "mf1", "pol");

  $der_def_pol .= $der_def_pol2;

  push(@out_c_pol1, @out_c_pol2);
  my $out_c_pol = join(", ", sort sort_alphanumerically @out_c_pol1);
  $out_c_pol = ", ".$out_c_pol if (! $out_c_pol =~ /^ *$/);

  # we join all the pieces
  my $maple_code1 = "
(* mf is the up potential *)
mzk   := (r0, r1, s0, s1, s2) -> \\
  $config{'simplify_begin'} \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2)) \\
  $config{'simplify_end'}:
mf0   := (r0, r1, s0, s1, s2) -> eval(mzk(r0, r1, s0, s1, s2)):
mf1   := (r0, r1, s0, s1, s2) -> eval(mzk(r1, r0, s2, s1, s0)):

\$include <util.mpl>
";

  my $maple_vrho0 = "vrho_0_ = mf0(".join(", ", @{$variables}).")";
  my $maple_vrho1 = "vrho_1_ = mf1(".join(", ", @{$variables}).")";

  # we build 2 variants of the functional, for unpolarized and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):

$der_def_unpol

$maple_code1
C([$maple_vrho0$out_c_unpol], optimize, deducetypes=false):
",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):

$der_def_pol

$maple_code1
C([$maple_vrho0, $maple_vrho1$out_c_pol], optimize, deducetypes=false):
"
      );

  maple2c_run(variables, $derivatives, \@variants, 1, $input_args, $output_args);
}

#####################################################################
# Process a MGGA functional for the energy
sub work_mgga_exc {
  my $variables, $derivatives;

  ($variables, $derivatives) = mgga_var_der();

  # get arguments of the functions
  $input_args  = "const double *rho, const double *sigma, const double *lapl, const double *tau";
  $output_args = ", double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, )";

  my ($der_def, @out_c) =
      maple2c_create_derivatives($variables, $derivatives, "mf");
  my $out_c = join(", ", @out_c);
  $out_c = ", $out_c" if ($out_c ne "");

  # we join all the pieces
  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> \\
  $config{'simplify_begin'} \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)) \\
  $config{'simplify_end'}:

(* mf is energy per unit volume *)
mf   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(dens(r0, r1)*mzk(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)):

\$include <util.mpl>
";
  my $maple_zk = " zk_0_ = mzk(".join(", ", @{$variables}).")";

  # we build 2 variants of the functional, for unpolarized and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
\n",
      );

  maple2c_run($variables, $derivatives, \@variants, 0, $input_args, $output_args);
}

#####################################################################
# Process a MGGA functional for the potential
sub work_mgga_vxc {
  my $variables, $all_derivatives;

  ($variables, $all_derivatives) = mgga_var_der();
  my $derivatives, $derivatives1, $derivatives2;
  ($derivatives, $derivatives1, $derivatives2) = filter_vxc_derivatives($all_derivatives);

  # get arguments of the functions
  $input_args  = "const double *rho, const double *sigma, const double *lapl, const double *tau";
  $output_args = "MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, )";

  # we obtain the missing pieces for maple
  # unpolarized calculation
  my ($der_def_unpol, @out_c_unpol) =
      maple2c_create_derivatives($variables, $derivatives1, "mf0", "unpol");
  my $out_c_unpol = join(", ", @out_c_unpol);
  $out_c_unpol = ", ".$out_c_unpol if (! $out_c_unpol =~ /^ *$/);

  # polarized calculation
  my ($der_def_pol, @out_c_pol1) =
      maple2c_create_derivatives($variables, $derivatives1, "mf0", "pol");
  my ($der_def_pol2, @out_c_pol2) =
      maple2c_create_derivatives($variables, $derivatives2, "mf1", "pol");

  $der_def_pol .= $der_def_pol2;

  push(@out_c_pol1, @out_c_pol2);
  my $out_c_pol = join(", ", sort sort_alphanumerically @out_c_pol1);
  $out_c_pol = ", ".$out_c_pol if (! $out_c_pol =~ /^ *$/);

  # we join all the pieces
  my $maple_code1 = "
(* mf is the up potential *)
mzk   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> \\
  $config{'simplify_begin'} \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)) \\
  $config{'simplify_end'}:
mf0   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(mzk(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)):
mf1   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(mzk(r1, r0, s2, s1, s0, l1, l0, tau1, tau0)):

\$include <util.mpl>
";

  my $maple_vrho0 = "vrho_0_ = mf0(".join(", ", @{$variables}).")";
  my $maple_vrho1 = "vrho_1_ = mf1(".join(", ", @{$variables}).")";

  # we build 2 variants of the functional, for unpolarized and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

$der_def_unpol

$maple_code1
C([$maple_vrho0$out_c_unpol], optimize, deducetypes=false):
",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

$der_def_pol

$maple_code1
C([$maple_vrho0, $maple_vrho1$out_c_pol], optimize, deducetypes=false):
"
      );

  maple2c_run($variables, $derivatives, \@variants, 1, $input_args, $output_args);
}
