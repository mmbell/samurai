#!/usr/bin/perl

# Constants

$_A_="5.0065";
$_B_="19.83923";
$C_P="1005.7";
$R_D="287.";
$E_3="6.1078";
$T_3="273.15";
$EPSILON="0.622";

sub w_sat {
    my $t = shift;
    my $p = shift;
#/*
## Saturation mixing ratio in g/kg at t (deg. K) and p (mb)
##/
	my $e = e_sw ($t);

	return (1000.0 * $EPSILON * $e / ($p - $e));
}


sub e_sw {

    my $t = shift;

#/*
# * Saturation vapor pressure in mb at temperature t (in deg. K)
# */

	my $index;
	my $frac;
#/*
# * Interpolate from the table of vapor pressures if ETBL_MIN < t < ETBL_MAX
# */
#	if (ETBL_MIN < t && t < ETBL_MAX)
#	{
#		index = (t - ETBL_MIN) / ETBL_STEP;
#		frac = fmod (t, ETBL_STEP) / ETBL_STEP;
#		return ((1.0 - frac) * E_tbl[index] + frac * E_tbl[index+1]);
#	}
#/*
## Otherwise, use the exact (but slow) formula
##/
#	else
		return ($E_3 * exp ($_A_ * log ($T_3 / $t)) * 
			exp (($_A_ + $_B_) * (1 - $T_3 / $t)));
}


sub theta_dry {

    my $t = shift;
    my $p = shift;
#/*
## Potential temperature (deg. K) of dry air at t (deg. K) and p (mb)
## 
##			      (R /C )
##	                        d  p
##	theta = t * (1000 / p)
##/
    my ($u, $diff, $diff2, $diff3);

    $u = 1000. / $p;
#/*
## Use Taylor series expansion about 700mb for pressures down to 500mb
##/
    if ($p > 500.) {
	$diff = $u - 1000. / 700.;
	$diff2 = $diff * $diff;
	$diff3 = $diff * $diff2;
	
	return ($t * (1.10714599 +
		     0.22116497 * $diff + 
		     -0.05531763 * $diff2 + 
		     0.02213146 * $diff3 + 
		     -0.01051376 * $diff * $diff3 + 
		     0.00546766 * $diff2 * $diff3 +
		     -0.00300743 * $diff3 * $diff3));
    }
#
# Use Taylor series expansion about 350mb for pressures down to 250mb
#/
    elsif ($p > 250.) {
	$diff = $u - 1000. / 350.;
	$diff2 = $diff * $diff;
	$diff3 = $diff * $diff2;
	
	return ($t * (1.34930719 +
		     0.13476972 * $diff +
		     -0.01685425 * $diff2 +
		     0.00337152 * $diff3 + 
		     -0.00080084 * $diff * $diff3 + 
		     0.00020824 * $diff2 * $diff3 +
		     -0.00005727 * $diff3 * $diff3));
	
    }
#
# Use Taylor series expansion about 175mb for pressures down to 125mb
#/
    elsif ($p > 125.) {
	$diff = $u - 1000. / 175.;
	$diff2 = $diff * $diff;
	$diff3 = $diff * $diff2;
	
	return ($t * (1.64443524 +
		     0.08212365 * $diff +
		     -0.00513518 * $diff2 +
		     0.00051362 * $diff3 + 
		     -0.00006100 * $diff * $diff3 + 
		     0.00000793 * $diff2 * $diff3 + 
		     -0.00000109 * $diff3 * $diff3));
    }
#
# Otherwise, use the exact form
#/
    else {
	return ($t * pow ($u, .28537338));
    }
}
	
sub lcl_temp {
    my $temp = shift;
    my $dp = shift;
#
# Calculate the temperature (deg. K) of the lifting condensation 
# level from the surface temp (deg. K) and dewpoint (deg. K)
#/

    if ($dp > $temp) {
	return (0.0);	# error */
    }
    return (1. / (1. / ($dp - 56.) + log ($temp / $dp) / 800.) + 56.);
}


sub lcl_pres {
    my $dp = shift;
    my $temp = shift;
    my $pres = shift;
#
# Find the pressure (mb) of the lifting condensation level from the 
# surface dewpoint (deg. K), temperature (deg. K), and pressure (mb).
#
# Adapted from Carl Mohr's RSANAL program.
#/

    my $w = w_sat ($dp, $pres);	# mixing ratio#/
    my $theta = theta_dry ($temp, $pres);
    my $plcl = $pres;
    my $test = 1.0;
    
    if ($dp > $temp) {
	return (0.0);	# error#/
    }
    while ($test > 0.01 || $test < -0.01) {
	$test = t_mr ($plcl, $w) - theta_to_t ($theta, $plcl);
	$plcl *= pow (2.0, (0.02 * $test));
    }

    return ($plcl);
}


sub t_mr {

    my $p = shift;
    my $w = shift;
#
# Calculate the temperature (deg. K) of saturated air at
# pressure p (in mb) and with mixing ratio w (in g/kg)
#
# The vapor pressure e is calculated using Note 2 of the Herzegh memo.
# The formula for the temperature was derived from Note 6 by substituting
# a two term Taylor series expansion for the ln(T_3/T) term and solving for
# T.
#/

    my $e = ($p * $w) / (1000. * $EPSILON + $w);
	
    return ($T_3 * $_A_ / 
	    ($_A_ - $_B_ + sqrt (square ($_B_) + 2 * $_A_ * log ($E_3 / $e))));
}


sub theta_to_t {

    my $theta = shift;
    my $p = shift;
#
# Temperature of dry air (deg. K) at potential temperature (deg. K) and p (mb)
#/

    return ($theta * pow ((1000.0 / $p), -$R_D / $C_P));
}



sub theta_e{
    my $t = shift; 
    my $dp = shift;
    my $p = shift;
#
# Equivalent potential temperature at temperature t (K), dewpoint dp (K),
# and pressure p (mb)
#/

    my $w = w_sat ($dp, $p);
#    my $theta = theta_dry ($t, $p);
    # Take into account mixing ratio
    my $theta = $t * pow ((1000.0 / $p),(.28537338*(1-0.00028*$w)));

    my $t_l = lcl_temp ($t, $dp);
    if ($t_l != 0) {
	return ($theta * exp ((3.376 / $t_l - 0.00254) * $w * (1 + 0.00081 * $w)));
    } else {
	return (-32767.);
    }
}


sub t_sat {
    my $ept = shift;
    my $p = shift;
#
# Temperature of saturated air (deg. K) with equivalent potential
# temperature ept (deg. K) and pressure p (mb)
#/

    my $t_s = $T_3;
    my $delta = 60.0;
    my $x;

    $x = $ept - theta_e ($t_s, $t_s, $p);

    while ($x > 0.01 || $x < -0.01) {
	$t_s += $x > 0.0 ? $delta : -$delta;
	$delta /= 2.0;
	if ($delta == 0.0) {
	    $delta = 60.0;
	}
	$x = $ept - theta_e ($t_s, $t_s, $p);
    }
    
    return ($t_s);
}



sub dewpoint {
    my $e = shift;
#
# Calculate the dewpoint from the vapor pressure
#/
    my $u = log ($e / $E_3);

    return (237.3 * $u / (17.2694 - $u) + $T_3);
}


sub e_from_dp {
    my $dp = shift;
#
# Calculate the vapor pressure from the dewpoint
# (This uses Note 5 from the Herzegh memo, inverted to represent e in
#  terms of the dewpoint)
#/

    return ($E_3 * exp (17.2694 * (1.0 - 237.3 / ($dp - $T_3 + 237.3))));
}



sub t_v {
    my $t = shift;
    my $p = shift;
    my $e = shift;
#
# * Calculate the virtual temperature from the temperature (K), pressure (mb), 
# and vapor pressure (mb)
#/

    return ($t / (1 - ($e / $p) * (1 - $EPSILON)));
}




sub t_wet {
    my $t = shift;
    my $p = shift;
    my $rh = shift;
#
# Calculate wet bulb temperature from temperature (K), pressure (mb), and 
# relative humidity (%).
# 
# This is an iterative solution using notes 3 and 6 of the Herzegh
# memo of 18 March 88 and tweaking t_wet until we're close (This seems
# like a poor algorithm, but I don't have time to make it better)
#/

    my ($e_sat, $e, $e_left, $e_right, $e_guess);
    my ($tw_left, $tw_right, $tw_guess);
#
# Find our vapor pressure
#/
    $e_sat = e_sw ($t);
    $e = $rh * 0.01 * $e_sat;
#
# Initialize
#/
    $tw_left = $t - 80; # Assume t_w is within 80 deg. of t */
    $e_left = e_sw ($tw_left) - $p * ($t - $tw_left) * 0.00066 * 
	(0.6859 + 0.00115 * $tw_left);

    $tw_right = $t;
    $e_right = $e_sat;

    while (1)
    {
	$tw_guess = $tw_left + ($e - $e_left) / ($e_right - $e_left) * 
	    ($tw_right - $tw_left);
	
	$e_guess = e_sw ($tw_guess) - $p * ($t - $tw_guess) * 0.00066 * 
	    (0.6859 + 0.00115 * $tw_guess);

    #
    # Quit when the vapor pressure for our guess is within 0.01 mb of the 
    # actual
    #/
	if (fabs ($e_guess - $e) < 0.01) {	
	    last;
	}

	if ($e_guess < $e)
	{
	    $tw_left = $tw_guess;
	    $e_left = $e_guess;
	} else {
	    $tw_right = $tw_guess;
	    $e_right = $e_guess;
	}
    }
    
    return ($tw_guess);
}




sub square {
    my $x = shift;
    
    return ($x * $x);
}


sub ten_to_the {
    my $x = shift;

    return (pow (10.0, $x));
}
