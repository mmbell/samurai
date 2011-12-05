#!/usr/bin/perl -w
require 'met_formulas.pl';

open (IN, "0912-allfields.txt") or die;
while(<IN>) {
  @line = split;
#  if ($line[0] == 2.0) {
  if ($line[0] !~ /Height/) {
    # Calculate hstar
      $z = $line[0]/.25;
      $r = int($line[1]);
    $rho = $line[5];
    #$rhoprime = $rho - 1.4683*exp(-0.00013592*2000);
    $temp = $line[6];
    $press = $temp*287*$rho/100;
    $qsat = w_sat($temp,$press);
#    print "TPQ: $temp, $press, $qsat\n";
    $hstar = 1005.7*$temp + 2.5e3*$qsat + 9.81*$line[0]*1000;
#    if (($line[0] == 0) or ($line[1] == 0)) {
	$psi[$z][$r] = 0;
#    } else {
#	$psi[$z][$r] = -$line[3]*($line[1]*1000)*$rho*250. + $psi[$z-1][$r];
##	print "$z $r $psi[$z][$r]\n";
#    }
#    if (($z >0) and ($r>0)) {
#    	$ucheck = -($psi[$z][$r] - $psi[$z-1][$r])/(250.*$r*1000*$rho);
##	print "U: $z $r $ucheck\n";
#    }     
#    printf "%13e %13e %13e %13e %13e %13e %13e\n", $line[1], $line[2], $line[3], $line[4], $hstar, $line[7], $rho;
    printf "%13e %13e %13e %13e %13e %13e %13e\n", $line[0]*1000, $line[1], $line[2], $psi[$z][$r], $hstar, $line[7], $rho; 
  }     
}
