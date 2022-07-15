#!/usr/bin/env perl

# Copyright (C) 2006-2007 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

if(@ARGV < 2) {
    print STDERR "Usage: get_funcs.pl srcdir builddir\n";
    exit(1);
}

$srcdir = shift;
$builddir = shift;

my @funcs = ("hyb_lda", "lda", "hyb_gga", "gga", "hyb_mgga", "mgga");
my %all_ids;

open(DOCS, ">$builddir/libxc_docs.txt") or die("Could not open '$builddir/libxc_docs.txt.'\n");

$s0 = ""; $s3 = ""; $s4 = ""; $s5 = ""; $s6 = "";
$publiclist = ""; $xclist = ""; $fxclist = ""; $xcf90list = ""; $xcfclist = "";

foreach $func (@funcs){
  undef %deflist_f;
  undef %deflist_c;
  undef %num;

  read_file($srcdir, $func);

  $s1 = ""; $s2 = "";
  foreach $key (sort { $a <=> $b } keys %deflist_f) {
    $s = sprintf "%s %-30s %3s  /*%-70s*/\n", "#define ",
      $deflist_f{$key}, $key, $deflist_c{$key};

    if($key < 100000){
      $s0 .= $s;
    }else{
      $s6 .= $s;
    }

    $t = $deflist_f{$key};
    $t =~ s/XC_(.*)/\L$1/;

    $s4 .= ",\n" if($s4);
    $s4 .= sprintf "{\"%s\", %d}", $t, $key;

    if($key < 100000){
      $s3 .= sprintf "  %s %-30s = %3s  ! %s\n", "integer, parameter ::",
        $deflist_f{$key}, $key, $deflist_c{$key};

      $s5 .= sprintf "  %s %-30s = %3s  ! %s\n", "integer(c_int), parameter, public ::",
      $deflist_f{$key}, $key, $deflist_c{$key};
    }

    $s1 .= "extern xc_func_info_type xc_func_info_$t;\n";
    $s2 .= "  &xc_func_info_$t,\n";
  }

  open(OUT, ">$builddir/funcs_$func.c") or die("Could not open '$builddir/funcs_$func.c'.\n");
  print OUT <<EOF
#include "util.h"

$s1

const xc_func_info_type *xc_${func}_known_funct[] = {
$s2  NULL
};
EOF
    ;
  close OUT;
}

close DOCS;

open(OUT, ">$builddir/funcs_key.c") or die("Could not open '$builddir/funcs_key.c'.\n");
print OUT <<EOF
#include "util.h"

xc_functional_key_t xc_functional_keys[] = {
$s4,
{"", -1}
};
EOF
;

open(OUT, ">$builddir/xc_funcs.h") or die("Could not open '$builddir/xc_funcs.h'.\n");
print OUT $s0;
close OUT;

open(OUT, ">$builddir/xc_funcs_worker.h") or die("Could not open '$builddir/xc_funcs_worker.h'.\n");
print OUT $s6;
close OUT;

open(OUT, ">$builddir/libxc_inc.f90") or die("Could not open '$builddir/libxc_incs.f90'.\n");
print OUT <<EOF
$s5
EOF
  ;
close OUT;

sub read_file() {
  my ($dir, $type) = @_;
  $type =~ s/(.*)/\L$1/;

  my $TYPE = $type;
  $TYPE =~ s/(.*)/\U$1/;

  # we remove the hyb from the filenames
  $xc_info_exe = "$builddir/xc-info";

  opendir(DIR, "$dir/") || die "cannot opendir '$dir': $!";
  while($_ = readdir(DIR)){
    $ftype = ${type};
    $ftype =~ s/^hyb_//;
    next if(!/(?:hyb_)?${ftype}_.*\.c$/);

    $file = $_;
    open(IN, "<$dir/$_") or die("Could not open '$dir/$_'.\n");
    while($_=<IN>){
      if(/#define\s+(XC_(${TYPE})_\S+)\s+(\S+)\s+\/\*(.*?)\s*\*\//){
        $deflist_f{$3} = $1;
        $deflist_c{$3} = $4;
        $num{$1} = $3;

        # check if ID is already in use
        if ( $all_ids{$3} ){
          printf stderr "Error: ID $2 repeated in\n  $1\n  $all_ids{$2}\n";
          exit 1;
        }else{
          $all_ids{$3} = $1;
        }
      }

      if(/^(const |)xc_func_info_type xc_func_info_${type}/){
        $infostr = "";
        while($_=<IN>){
          if(/([^}])*};/) {
            $infostr .= $1;
            last;
          }
          # remove C comments
          $_ =~ s|/\*[^\*]*\*/||;
          # remove braces
          $_ =~ s/[{}]//g;
          chomp($_);
          $infostr .= $_;
        }
        @infos = split('"', $infostr);
        @infos0 = split(',', $infos[0]);
        $infos0[0] =~ s/^\s*//;
        $infos0[1] =~ s/^\s*//;
        @infos2 = split(',', $infos[2]);

        for($ii = 0; $ii <= $#infos2; $ii++) {
          # remove leading spaces
          $infos2[$ii] =~ s/^\s*//;
        }

        print DOCS "Number         : $num{$infos0[0]}\n";
        print DOCS "File           : $file\n";
        print DOCS "Codename       : $infos0[0]\n";
        print DOCS "Kind           : $infos0[1]\n";
        print DOCS "Description 1  : $infos[1]\n";
        $deflist_c{$num{$infos0[0]}} =~ s/^\s*//;
        $deflist_c{$num{$infos0[0]}} =~ s/\s*$//;
        if($deflist_c{$num{$infos0[0]}} ne $infos[1]) {
          print DOCS "Description 2  : $deflist_c{$num{$infos0[0]}}\n";
        }
        #infos2[0] will be blank
        print DOCS "Family         : $infos2[1]\n";

        if(-e "$xc_info_exe" && -x "$xc_info_exe") {
          $xc_info = `$xc_info_exe $num{$infos0[0]}`;
          @refs = split('\n', $xc_info);
          if($refs[4] =~ /Reference\(s\)/) {
            print DOCS "References     : ";
            print DOCS $refs[5] . "\n";
            $ref_start = 6;
          } else {
            print DOCS $refs[4] . "\n";
            print DOCS "References     : ";
            print DOCS $refs[7] . "\n";
            $ref_start = 8;
          }
          for($ii = $ref_start; $ii <= $#refs; $ii++) {
            print DOCS "                 " . $refs[$ii] . "\n";
          }
        } else {
          # print only the names of the variables in references.c
          print DOCS "References     :";
          for($ii = 2; $ii <= 6; $ii++) {
            if($infos2[$ii] ne "NULL") {
              $infos2[$ii] =~ s/&xc_ref_//;
              print DOCS " $infos2[$ii]";
            }
          }
          print DOCS "\n";
        }

        if(($infos2[7] =~ /XC_FLAGS_(.)D/) != 1) {
          print STDERR $infos2[7], "\n";
          print STDERR "$infos0[0]: Must set exactly one dimensionality flag.\n";
          exit(1);
        }
        print DOCS "Dimensionality : $1\n";

        print DOCS "Quantities     : ";
        @quantities = ($infos2[7] =~ /XC_FLAGS_HAVE_(.XC)/g);
        print DOCS join(" ", @quantities) . "\n";

        $infos2[7] =~ s/XC_FLAGS_.D//;
        $infos2[7] =~ s/XC_FLAGS_HAVE_.XC//g;
        $infos2[7] =~ s/\|//g;
        $infos2[7] =~ s/^\s*//;
        $infos2[7] =~ s/^s*$//;

        print DOCS "Other flags    : $infos2[7]\n";

        open(IN, "<$srcdir/$file");
        chomp(my @lines = <IN>);
        close(IN);

        $shortname = lc(substr($infos0[0], 3));
        @lines = grep {/xc_${shortname}_set_params\(xc_func_type/} @lines;
        $set_params = shift @lines;

        if($set_params ne "") {
          if($set_params !~ /void/) {
            $set_params = "void $set_params";
          }
          print DOCS $set_params . "\n";
        }

        print DOCS "min dens       : $infos2[8]\n";
        print DOCS "min grad       : $infos2[9]\n";
        print DOCS "min tau        : $infos2[10]\n";
        print DOCS "min zeta       : $infos2[11]\n";
        print DOCS "init           : $infos2[12]\n";
        # apparently, the end is always NULL
        print DOCS "end            : $infos2[13]\n";
        print DOCS "work lda       : $infos2[14]\n";
        print DOCS "work gga       : $infos2[15]\n";
        print DOCS "work mgga      : $infos2[16]\n";
        print DOCS "----------------------------\n";

        if($num{$infos0[0]} eq "") {
          print STDERR "ERROR: missing number\n";
          print STDERR $infos0[0], "\n";
          exit(1);
        }

        if($deflist_f{$num{$infos0[0]}} ne $infos0[0]) {
          print STDERR $deflist_f{$num{$infos0[0]}} . " " . $infos0[0] . "\n";
          print STDERR "Mismatch of names.\n";
          exit(1);
        }
      }
    }
    close(IN);
  }
  closedir DIR;
}
