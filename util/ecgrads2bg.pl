#!/usr/bin/perl -w
require 'met_formulas.pl';
use POSIX;
# Use wgrib to parse the ECMWF analysis for a SAMURAI background field
# Template for commands
#$wgrib = "wgrib -s $infile | grep \":$var:$press mb:\" | wgrib -i -text $infile -o $outfile";
if ($#ARGV < 0) {
  print "Usage: ecmwfgrads2bg.pl <dtg>\n";
  exit();
}
$dtg = $ARGV[0];
# Modify the center of the background field
$tclat = 18.00;
$tclon = 135.00;
$maxradius = 2000;

$hr = substr($dtg,-2);
$year = substr($dtg,0,4);
$mon = substr($dtg,4,2);
$day = substr($dtg,6,2);
$sec = "00";
$min = "00";
# Convert time to seconds and record as time
$platform = `uname`;
$linuxstr = $year."-".$mon."-".$day." ".$hr.":".$min.":".$sec;
$macstr = $mon.$day.$hr.$min.$year.".".$sec;
if ($platform =~ /Linux/) {
       $inttime = `date -u -d '$linuxstr' +%s`;
print "datestr=$linuxstr and inttime=$inttime\n";
} elsif ($platform =~ /Darwin/) {
       $inttime = `date -j -u +%s "$macstr"`;
print "datestr=$macstr and inttime=$inttime\n";
} else {
       print "Not sure about $platform, trying Linux date command...\n";
       $inttime = `date -u -d '$linuxstr' +%s`;
}
chomp($inttime);
print "Reading $dtg ($inttime Unix time)\n";

$ecout = $dtg."_samuraibg.txt";
open (EC, ">$ecout") or die;
$fac_lon = 111.41*cos($tclat/57.296);
$fac_lat = 111.13;
print "Getting data within $maxradius km of $tclat, $tclon\n";
#@pressfields = qw(Z U V W T Q);
#@sfcfields = qw(Z 10U 10V SP 2T 2D);
@presslevels = qw(1000 975 950 925 900 850 800 750 700 650 600 550 500 450 400 350 300 250 200 150 100 70 50);
@pvars = ("geop", "uwnd", "vwnd", "wwnd", "tmpk", "shum");
@sfcvars = ("null", "usfc", "vsfc", "pmsl", "tsfc", "tdsf");

#Start with the surface
$zdim = 1;
for $n (1 .. $#sfcvars) {
  $var = $sfcvars[$n];
  $gradsfile = $dtg."_".$var.".dat";
  open(IN, $gradsfile) or die("Couldn't open $gradsfile :$!\n");
  binmode IN;
  $xdim = 841;
  $ydim = 401;
  $zdim = 1;
  $numvar = 1;
  $tdim = 1;
  $bufsize = 4*$xdim*$ydim;
  for $t (0 .. $tdim-1) {
    print "Reading $var...\n";
    for $z (0 .. $zdim-1) {
      read(IN,$buf,$bufsize);
      #$swapbuf = reverse $buf;
      @levarr = unpack("f*",$buf);
      #@levarr = reverse @levarr;
      #print "@levarr\n";
      $count=0;
      for $y (0 .. $ydim-1) {
	for $x (0 .. $xdim-1) {
	  $cart_data[$n][$x][$y][$z]=$levarr[$count];
	  $count++;
	}
      }
      #print "$n, 100, 100, 0 = $cart_data[$n][100][100][0]\n";
    }
  }
  close(IN);
  
}

# Adjust the surface fields to relevant variables
for $i (0 .. $xdim-1) {
  for $j (0 .. $ydim-1) {
      $lat = -20+$j*0.25;
      $lon =  60+$i*0.25;
      $ydist = ($lat-$tclat)*$fac_lat;
      $xdist = ($lon-$tclon)*$fac_lon;
      $dist = sqrt($xdist**2 + $ydist**2);
      if ($dist < $maxradius) {
	  $cart_data[0][$i][$j][0] = 10.;
	  $temp = $cart_data[4][$i][$j][0];
	  $dewp = $cart_data[5][$i][$j][0];
	  $press = $cart_data[3][$i][$j][0]/100.;
	  $vp = e_from_dp($dewp);
          $rhoa = 100*($press-$vp)/($temp*287);
	  $qv = 1000*.622*$vp/($press-$vp);
	  if ($qv < 0) {
	    # Problem!
	    print "Negative Qv at $i, $j!\n";
	    print "$temp, $dewp, $press, $rhoa, $vp, $qv\n";
	    exit();
	  }
	  $cart_data[3][$i][$j][0] = 0.;
	  #$cart_data[4][$i][$j][0] = 1005.7*$temp + 9.80665*$cart_data[0][$i][$j][0] + 2.5e3*$qv;
	  $cart_data[5][$i][$j][0] = $qv;
	  $cart_data[6][$i][$j][0] = $rhoa;
	  printf EC ("%d\t%e\t%e\t",$inttime,$lat,$lon);
	  for $f (0 .. $#pvars+1) {
	      printf EC ("%e\t",$cart_data[$f][$i][$j][0]);
	  }
	  print EC "\n";
      }
  }
}

#Now move up the grid
for $n (0 .. $#pvars) {
  $var = $pvars[$n];
  $gradsfile = $dtg."_".$var.".dat";
  open(IN, $gradsfile) or die("Couldn't open $gradsfile :$!\n");
  binmode IN;
  $xdim = 841;
  $ydim = 401;
  $zdim = 26;
  $numvar = 1;
  $tdim = 1;
  $bufsize = 4*$xdim*$ydim;
  for $t (0 .. $tdim-1) {
    print "Reading $var...\n";
    for $z (1 .. $zdim) {
      read(IN,$buf,$bufsize);
      #$swapbuf = reverse $buf;
      @levarr = unpack("f*",$buf);
      #@levarr = reverse @levarr;
      #print "\tUnpacking level $z...\n";
      $count=0;
      for $y (0 .. $ydim-1) {
	for $x (0 .. $xdim-1) {
	  $cart_data[$n][$x][$y][$z]=$levarr[$count];
	  $count++;
	}
      }
      #print "$n, 100, 100, $z = $cart_data[$n][100][100][$z]\n";
      
    }
  }
  close(IN);
}

  # Adjust the fields to relevant variables
for $p (0 .. $#presslevels) {
  $z = $p+1;	    
  for $i (0 .. $xdim-1) {
      for $j (0 .. $ydim-1) {
	$lat = -20+$j*0.25;
	$lon =  60+$i*0.25;
	$ydist = ($lat-$tclat)*$fac_lat;
	$xdist = ($lon-$tclon)*$fac_lon;
	$dist = sqrt($xdist**2 + $ydist**2);
	if ($dist < $maxradius) {	  

	  #$cart_data[0][$i][$j][$z] = $cart_data[0][$i][$j][$z]/9.80665;
	  $temp = $cart_data[4][$i][$j][$z];
	  $press = $presslevels[$p];
	  $cart_data[3][$i][$j][$z] /= -(9.80665*$rhoa);
	  $qv = $cart_data[5][$i][$j][$z];
          $vp = $press * $qv/(0.622 + $qv);
          $rhoa = 100*($press-$vp)/($temp*287);
	  $qv *= 1000;
	  #$cart_data[4][$i][$j][$z] = 1005.7*$temp + 9.80665*$cart_data[0][$i][$j][$z] + 2.5e3*$qv;
	  $cart_data[5][$i][$j][$z] = $qv;
	  $cart_data[6][$i][$j][$z] = $rhoa;

	  printf EC ("%d\t%e\t%e\t",$inttime,$lat,$lon);
	  for $f (0 .. $#pvars+1) {
	    printf EC ("%e\t",$cart_data[$f][$i][$j][$z]);
	    #print "$f, $i, $j, $z defined...\n";
	  }
	  print EC "\n";
	}
      }
    }
}

close(EC);

# Now write the columns for easy vertical interpolation
$bg = $dtg."_Background.in";
open (BG, ">$bg") or die;
for $i (0 .. $xdim-1) {
  for $j (0 .. $ydim-1) {
    $lat = -20+$j*0.25;
    $lon =  60+$i*0.25;
    $ydist = ($lat-$tclat)*$fac_lat;
    $xdist = ($lon-$tclon)*$fac_lon;
    $dist = sqrt($xdist**2 + $ydist**2);
    if ($dist < $maxradius) {
      for $k (0 .. $#presslevels+1) {
	printf BG ("%d\t%e\t%e\t",$inttime,$lat,$lon);
	for $f (0 .. $#pvars+1) {
	  printf BG ("%e\t",$cart_data[$f][$i][$j][$k]);
	  if ($cart_data[$f][$i][$j][$k] eq "") {
	    print "$f, $i, $j, $k undefined!\n";
	    $stop = <STDIN>;
	  }
	}
	print BG "\n";
      }
    }
  }
}
