#!/usr/bin/perl -w
use POSIX;
require 'met_formulas.pl';
$pi = acos(-1.);

open (IN, "0912-allfields.txt") or die;
while(<IN>) {
    @line = split;
    if ($line[0] !~ /Height/) {
	$z = $line[0]*1000;
	$r = int($line[1]);
	$vt = $line[2];
	$vr = $line[3];
	$wspd = sqrt($vt*$vt + $vr*$vr);
	$wdir = atan2(-$vr,-$vt)*180/$pi;
	$w = $line[4];
	$rho = $line[5];
    #$rhoprime = $rho - 1.4683*exp(-0.00013592*2000);                                                                                              
	$temp = $line[6];
	$qv = $line[7];
	$press = $temp*287*$rho/100;
	$e = ($press*$qv/1000)/($qv/1000 + 0.622);
	#print "$e\n";
	$dewp = dewpoint($e)-273.15;
        $temp -= 273.15;
	$latrad = 25. * $pi/180.0;
	$fac_lat = 111.13209 - 0.56605 * cos(2.0 * $latrad)
	  + 0.00012 * cos(4.0 * $latrad) - 0.000002 * cos(6.0 * $latrad);
	$fac_lon = 111.41513 * cos($latrad)
	  - 0.09455 * cos(3.0 * $latrad) + 0.00012 * cos(5.0 * $latrad);
	$oblon = -90.+$r/$fac_lon;
	$oblat = 25.;

	#print " 120000 $oblat $oblon -999 -999 -999 -999 $z $press $wdir $wspd $temp $dewp -999 -999 -999 $w\n";
        print " 120000 $oblat $oblon -999 -999 -999 -999 $z $press -999 -999 $temp $dewp -999 -999 -999 $w\n";
    }
}
#ob.setLat(lineparts[1].toFloat());
#ob.setLon(lineparts[2].toFloat());
#ob.setAltitude(lineparts[7].toFloat());
#ob.setPressure(lineparts[8].toFloat());
#ob.setTemperature(lineparts[11].toFloat() + 273.15);
#ob.setDewpoint(lineparts[12].toFloat() + 273.15);
#ob.setWindDirection(lineparts[9].toFloat());
#ob.setWindSpeed(lineparts[10].toFloat());
#ob.setVerticalVelocity(lineparts[16].toFloat());
#ob.setObType(MetObs::flightlevel);
