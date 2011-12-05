#!/usr/bin/perl -w
use POSIX;

# samurai_out2bg.pl
# By Michael Bell
#

# Modify the following variables
$tclat = 18.0;
$tclon = 150.456451;
$time = 1221350400;
$tcu = -4.0;
$tcv = 0.0;

## Datasets (Input files)
if ($#ARGV < 0) {
  print "Usage: out2bg.pl <samurai_XYZ_*.out>\n";
  exit();
} else {
  $infile = $ARGV[0];
}

# Output files
$outfile = "samurai_newBackground.in";
open(OUT,">$outfile") or die;

# Constants
$pi = acos(-1);
$dtor=$pi/180.;
$f_lat = $tclat;
$coriolisf = 2*(7.292)*sin($f_lat*$dtor);
print "Coriolis parameter: $coriolisf (10-5s-1)\n";

# Read the samurai file
&read_dotout;

print "Processing...\n";
for($y=1;$y<=$maxy-1;$y++) {
  for($x=1;$x<=$maxx-1;$x++) {
    for($z=0;$z<=$maxz;$z++) {
      # Define the box
      $xdist = $xfirst + $x*$xdelta;
      $ydist = $yfirst + $y*$ydelta;
      $heightkm = $zfirst + $z*$zdelta;
      $latrad = $tclat * $pi/180.0;
      $fac_lat = 111.13209 - 0.56605 * cos(2.0 * $latrad)
	+ 0.00012 * cos(4.0 * $latrad) - 0.000002 * cos(6.0 * $latrad);
      $fac_lon = 111.41513 * cos($latrad)
	- 0.09455 * cos(3.0 * $latrad) + 0.00012 * cos(5.0 * $latrad);
      $lat = $tclat + $ydist/$fac_lat;
      $lon = $tclon + $xdist/$fac_lon;
      #@fields = qw(X Y Z rhoE u v w vort div qv' rho' T' P' h' udx udy udz vdx vdy vdz wdx wdy wdz rhowdz MC);
      #                      0 1 2 3 4    5   6   7    8  9  10 11  12  13  14  15  16  17  18  19  20     21
      $u = $fielddata[$x][$y][$z][1]+$tcu;
      $v = $fielddata[$x][$y][$z][2]+$tcv;
      $w = $fielddata[$x][$y][$z][3];
      $qvprime = $fielddata[$x][$y][$z][6];
      $rhoaprime = $fielddata[$x][$y][$z][7];
      $tprime = $fielddata[$x][$y][$z][8];
      $alt = $heightkm * 1000;
      if ($alt == 0) { $alt = 10; }
      $rhobar = getReferenceVariable(1, $alt);
      $rhoa = $rhobar + $rhoaprime/100.;
      $qbar = bhypInvTransform(getReferenceVariable(0, $alt));
      $qv = $qbar + $qvprime;
      if ($qv < 0) { $qv = 0; }
      $temp = $tprime + getReferenceVariable(4, $alt);
      printf OUT ("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", $time, $lat, $lon, $alt, $u, $v, $w, $temp, $qv, $rhoa);
    }
  }
}

sub read_dotout {
  open(IN, "$infile") or die;
  $x = $y = $z = -32768;
  $xdelta = $ydelta = $zdelta = 0;
  while(<IN>) {
    if (!/X/) {
      @line=split;
      if ($x == -32768) {
	$minx = $line[0];
	$miny = $line[1];
	$minz = $line[2];
	$x = $y = $z = 0;
      } else {
	if ($line[2] == $minz) {
	  $z = 0;
	  if ($line[1] == $miny) {
	    $y = 0;
	    if ($line[0] == $minx) {
	      $x = 0;
	    } else {
	      $x++;
	      if ($xdelta == 0) {
		$xdelta = $line[0] - $minx;
	      }
	    }
	  } else {
	    $y++;
	    if ($ydelta == 0) {
	      $ydelta = $line[1] - $miny;
	    }
	  }
	} else {
	  $z++;
	  if ($zdelta == 0) {
	    $zdelta = $line[2] - $minz;
	  }
	}
      }
      for $f (3 .. 24) {
	$fielddata[$x][$y][$z][$f-3] = $line[$f];
      }
    }
  }
  $maxx = $x;
  $maxy = $y;
  $maxz = $z;
  $xfirst = $minx;
  $yfirst = $miny;
  $zfirst = $minz;
  
}

sub getReferenceVariable
{
  my $refVariable = $_[0];
  my $heightm = $_[1];

  $qvbhypcoeff[0] = 9.4826;
  $qvbhypcoeff[1] = -0.0026721;
  $qvbhypcoeff[2] = 2.8312e-07;
  $qvbhypcoeff[3] = -1.3217e-11;
  $qvbhypcoeff[4] = 2.2749e-16;
  
  $rhoacoeff[0] = 1.1439;
  $rhoacoeff[1] = -0.00010117;
  $rhoacoeff[2] = 3.2486e-09;
  $rhoacoeff[3] = -3.4898e-14;
  $rhoacoeff[4] = -2.6925e-19;
  
  $dpdzcoeff[0] = -11.432;
  $dpdzcoeff[1] = 0.0010633;
  $dpdzcoeff[2] = -4.0545e-08;
  $dpdzcoeff[3] =  7.9634e-13;
  $dpdzcoeff[4] = -5.8778e-18;
    
  my $qvbhypref = 0;
  my $rhoaref = 1;
  my $rhoref = 2;
  my $href = 3;
  my $tempref = 4;
  my $pressref = 5;

  if ($refVariable == $qvbhypref) {
    my $qvbhyp = 0.;
    for ($i = 0; $i < 5; $i++) {
      $power = $heightm**$i; 
      $qvbhyp += $qvbhypcoeff[$i] * $power;
    }
    if ($qvbhyp < 0.) { $qvbhyp = 0.; }
    return $qvbhyp;
  } elsif ($refVariable == $rhoaref) {
    my $rhoa = 0.;
    for ($i = 0; $i < 5; $i++) {
      $power = $heightm**$i; 
      $rhoa += $rhoacoeff[$i] * $power;
    }
    return $rhoa;
  } elsif ($refVariable == $rhoref) {
    my $rho = 0.;
    my $qvbhyp = 0.;
    my $rhoa = 0.;
    for ($i = 0; $i < 5; $i++) {
      $power = $heightm**$i; 
      $rhoa += $rhoacoeff[$i] * $power;
      $qvbhyp += $qvbhypcoeff[$i] * $power;
    }
    if ($qvbhyp < 0.) { $qvbhyp = 0.; }
    my $qv = bhypInvTransform($qvbhyp);
    $rho = $rhoa*$qv/1000. + $rhoa;
    return $rho;
  } elsif (($refVariable == $href) or ($refVariable = $tempref) or ($refVariable == $pressref)) {
    # Integrate hydrostatic equation to get pressure and/or solve for T or h
    my $press = 0.;
    my $temp = 0.;
    my $rho = 0.;
    my $qvbhyp = 0.;
    my $rhoa = 0.;
    for ($i = 0; $i < 5; $i++) {
      $power = $heightm**$i;
      $power1 = $heightm**($i+1);
      $press += $dpdzcoeff[$i] * $power1 / ($i+1);
      $rhoa += $rhoacoeff[$i] * $power;
      $qvbhyp += $qvbhypcoeff[$i] * $power;
    }
    if ($qvbhyp < 0.) { $qvbhyp = 0.; }
    my $qv = bhypInvTransform($qvbhyp);
    $rho = $rhoa*$qv/1000. + $rhoa;
    $press += 101510.0;
    $temp = $press/(286.9*$rhoa + 461.5*$rhoa*$qv/1000.);
    my $h = 1005.7*$temp + 9.81*$heightm + 2.5e3*$qv;
    if ($refVariable == $href) {
      return $h;
    } elsif ($refVariable == $tempref) {
      return $temp;
    } elsif ($refVariable == $pressref) {
      return $press;
    } 
  }
  
  return 0;
}

sub bhypTransform
{
  
  my $qv = $_[0];
  my $qvbhyp = 0.5*(($qv + 1.e-7) - 1.e-14/($qv + 1.e-7));
  return $qvbhyp;
	
}

sub bhypInvTransform
{
  my $qvbhyp = $_[0];
  my $qv = 0.;
  if ($qvbhyp > 0) {
    $qv = sqrt($qvbhyp*$qvbhyp + 1.e-14) + $qvbhyp - 1.e-7;
  }
  return $qv;
}
