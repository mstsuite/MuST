#!/usr/bin/env perl

use maple2c_derivatives;
use Data::Dumper;
use Text::Wrap qw(wrap);

local $Text::Wrap::columns = 120;

my @replace = (
  '\[n0, s0, u0\]' , "",
  '\[n1, s2, u1\]' , "",
  '\[n0, n1, s0, s1, s2, u0, u1, ked1, ked2\]' , "",
  'Derivative' , "",

  '\[1, 0, 0\]\[ked1\]' , "ked1_vrho[0]",
  '\[0, 1, 0\]\[ked1\]' , "ked1_vsigma[0]",
  '\[0, 0, 1\]\[ked1\]' , "ked1_vlapl[0]",

  '\[2, 0, 0\]\[ked1\]' , "ked1_v2rho2[0]",
  '\[1, 1, 0\]\[ked1\]' , "ked1_v2rhosigma[0]",
  '\[1, 0, 1\]\[ked1\]' , "ked1_v2rholapl[0]",
  '\[0, 2, 0\]\[ked1\]' , "ked1_v2sigma2[0]",
  '\[0, 1, 1\]\[ked1\]' , "ked1_v2sigmalapl[0]",
  '\[0, 0, 2\]\[ked1\]' , "ked1_v2lapl2[0]",

  '\[3, 0, 0\]\[ked1\]' , "ked1_v3rho3[0]",
  '\[2, 1, 0\]\[ked1\]' , "ked1_v3rho2sigma[0]",
  '\[2, 0, 1\]\[ked1\]' , "ked1_v3rho2lapl[0]",
  '\[1, 2, 0\]\[ked1\]' , "ked1_v3rhosigma2[0]",
  '\[1, 1, 1\]\[ked1\]' , "ked1_v3rhosigmalapl[0]",
  '\[1, 0, 2\]\[ked1\]' , "ked1_v3rholapl2[0]",
  '\[0, 3, 0\]\[ked1\]' , "ked1_v3sigma3[0]",
  '\[0, 2, 1\]\[ked1\]' , "ked1_v3sigma2lapl[0]",
  '\[0, 1, 2\]\[ked1\]' , "ked1_v3sigmalapl2[0]",
  '\[0, 0, 3\]\[ked1\]' , "ked1_v3lapl3[0]",

  '\[4, 0, 0\]\[ked1\]' , "ked1_v4rho4[0]",
  '\[3, 1, 0\]\[ked1\]' , "ked1_v4rho3sigma[0]",
  '\[3, 0, 1\]\[ked1\]' , "ked1_v4rho3lapl[0]",
  '\[2, 2, 0\]\[ked1\]' , "ked1_v4rho2sigma2[0]",
  '\[2, 1, 1\]\[ked1\]' , "ked1_v4rho2sigmalapl[0]",
  '\[2, 0, 2\]\[ked1\]' , "ked1_v4rho2lapl2[0]",
  '\[1, 3, 0\]\[ked1\]' , "ked1_v4rhosigma3[0]",
  '\[1, 2, 1\]\[ked1\]' , "ked1_v4rhosigma2lapl[0]",
  '\[1, 1, 2\]\[ked1\]' , "ked1_v4rhosigmalapl2[0]",
  '\[1, 0, 3\]\[ked1\]' , "ked1_v4rholapl3[0]",
  '\[0, 4, 0\]\[ked1\]' , "ked1_v4sigma4[0]",
  '\[0, 3, 1\]\[ked1\]' , "ked1_v4sigma3lapl[0]",
  '\[0, 2, 2\]\[ked1\]' , "ked1_v4sigma2lapl2[0]",
  '\[0, 1, 3\]\[ked1\]' , "ked1_v4sigmalapl3[0]",
  '\[0, 0, 4\]\[ked1\]' , "ked1_v4lapl4[0]",
    );

my @ders = (
  "vrho", "vsigma", "vlapl", "vtau",

  "v2rho2", "v2rhosigma", "v2rholapl", "v2rhotau", "v2sigma2",
  "v2sigmalapl", "v2sigmatau", "v2lapl2", "v2lapltau", "v2tau2",

  "v3rho3", "v3rho2sigma", "v3rho2lapl", "v3rho2tau",
  "v3rhosigma2", "v3rhosigmalapl", "v3rhosigmatau", "v3rholapl2",
  "v3rholapltau", "v3rhotau2", "v3sigma3", "v3sigma2lapl",
  "v3sigma2tau", "v3sigmalapl2", "v3sigmalapltau", "v3sigmatau2",
  "v3lapl3", "v3lapl2tau", "v3lapltau2", "v3tau3",

  "v4rho4", "v4rho3sigma", "v4rho3lapl", "v4rho3tau",
  "v4rho2sigma2", "v4rho2sigmalapl", "v4rho2sigmatau", "v4rho2lapl2",
  "v4rho2lapltau", "v4rho2tau2", "v4rhosigma3", "v4rhosigma2lapl",
  "v4rhosigma2tau", "v4rhosigmalapl2", "v4rhosigmalapltau", "v4rhosigmatau2",
  "v4rholapl3", "v4rholapl2tau", "v4rholapltau2", "v4rhotau3",
  "v4sigma4", "v4sigma3lapl", "v4sigma3tau", "v4sigma2lapl2",
  "v4sigma2lapltau", "v4sigma2tau2", "v4sigmalapl3", "v4sigmalapl2tau",
  "v4sigmalapltau2", "v4sigmatau3", "v4lapl4", "v4lapl3tau",
  "v4lapl2tau2", "v4lapltau3", "v4tau4"
);
my @ders_def = ();

foreach $der (@ders){
  push @ders_def, maple2c_derivatives($der, "mgga");
}

# create replaces. E.g.
#  '\[1, 0, 0, 0, 0, 0, 0, 0, 0\]\[mgga\]' , "mgga_vrho[0]",
foreach $der (@ders_def){
  ${$der}[1] =~ s/_(\d+)_/[$1]/;
  push @replace, '\['.join(", ", @{${$der}[0]}).'\]\[mgga\]';
  push @replace, "mgga_".${$der}[1];
}

my %all_out = ();
foreach $der (@ders_def){
  # let us do one order each time
  next if(! (${$der}[1] =~ /^v4/));

  $mder = "";
  $mder = $mder." {n0, ".${$der}[0][0]."}";
  $mder = $mder.",{n1, ".${$der}[0][1]."}";
  $mder = $mder.",{s0, ".${$der}[0][2]."}";
  $mder = $mder.",{s1, ".${$der}[0][3]."}";
  $mder = $mder.",{s2, ".${$der}[0][4]."}";
  $mder = $mder.",{u0, ".${$der}[0][5]."}";
  $mder = $mder.",{u1, ".${$der}[0][6]."}";
  $mder = $mder.",{t0, ".${$der}[0][7]."}";
  $mder = $mder.",{t1, ".${$der}[0][8]."}";

  open(OUT, ">/tmp/math.m");
  print OUT "Print[ToString[FullSimplify[D[mgga[n0, n1, s0, s1, s2, u0, u1, ked1[n0, s0, u0], ked2[n1, s2, u1]], $mder]], FormatType -> InputForm]]";
  close(OUT);

  $out = `math -script /tmp/math.m`;
  chomp($out);

  for(my $j=0; $j<$#replace; $j+=2){
    my ($from, $to) = ($replace[$j], $replace[$j+1]);
    $out =~ s/$from/$to/g;

    $from =~ s/ked1/ked2/g;
    $to   =~ s/ked1/ked2/g;
    $out  =~ s/$from/$to/g;
  }

  # let us get rid of the powers
  $out =~ s/([A-Za-z][A-Za-z_\[\]0-9]*)\^2/$1*$1/g;
  $out =~ s/([A-Za-z][A-Za-z_\[\]0-9]*)\^3/$1*$1*$1/g;
  $out =~ s/([A-Za-z][A-Za-z_\[\]0-9]*)\^4/$1*$1*$1*$1/g;

  # we order the lines in a more useful way
  ${$der}[1] =~ /(.\d?)(.*?)\[(\d+)\]/;

  $out = ${$der}[1]." = ".$out.";\n";

  if($3 eq "0"){
    $out1 .= wrap('', '  ', $out);
  }else{
    $out2 .= wrap('  ', '    ', $out);;
  }
}

# convert tabs to spaces
$out1 =~ s/\t/        /g;
$out2 =~ s/\t/        /g;
print "$out1\nif(func->nspin == XC_POLARIZED){\n$out2}\n";
