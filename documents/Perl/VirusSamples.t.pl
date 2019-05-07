#!/usr/bin/perl

use strict;
use warnings;

BEGIN 
{
    unshift @INC, '../modules';
}

use Jls::VirusSamples;

my $founder = Jls::VirusSamples->new;

#use JSON::JSON;
#my $json = JSON->new->allow_blessed->convert_blessed->pretty( 1 );

# open( my $o_f, ">check.out" );
# print $o_f "@$birth_ages\n";
# close( $o_f );

use Statistics::R;
my $R = Statistics::R->new();

$R->run( 'set.seed( 31415 )' );

my $EPS = 0.005;
my $RELATIVE_ERROR = 1.0e-06;

my $m = 2; # Poisson mean number of offspring
my $mu = 10.0; # mean time to birth
my $sigma = 0.1; # st dev time to birth

my $q = Jls::Virus::extinction_probability( $m );
my $q0 = 0.2031878;
die unless abs( 1.0 - $q / $q0 ) < $RELATIVE_ERROR;

my $m0 = Jls::Virus::harris_sevastyanov( $m, $RELATIVE_ERROR );
die unless abs( 1.0 - (1.0 - $q) * $m / $m0 ) < $RELATIVE_ERROR;

# Checks mean of R Poisson random number generator.

my $n = 1000;

my $cmd_s = <<EOF;
    array <- rpois( n = $n, lambda = $m ) 
    mean( array )
EOF

my $mean = Jls::RIO::o_perl( $R->run( $cmd_s ) ); # mean Poisson $m = 2.
my $mean_0 = 1.994;
die unless abs( 1 - $mean / $mean_0 ) < $EPS;

# Checks mean of R gamma random number generator.

my $shape = ($mu / $sigma) * ($mu / $sigma);
my $lambda = $mu / $sigma / $sigma;

$cmd_s = <<EOF;
    array <- rgamma( n = $n, shape = $shape, rate = $lambda ) 
    mean( array )
EOF

$mean = Jls::RIO::o_perl( $R->run( $cmd_s ) ); 
# mean gamma $shape / $lambda = = 10000 / 1000 = 10.
$mean_0 = 10.00549;
die unless abs( 1 - $mean / $mean_0 ) < $EPS;

# Checks random birth times specifically for Virus package.

$R->run( 'set.seed( 31415 )' );

$m = 10; # Poisson mean number of offspring
$mu = 2.0; # mean time to birth
$sigma = 0.2; # st dev time to birth

# deterministic $mu = 1.0 because ! defined $mu (1.0) ! defined $sigma (0.0)
my $birth_ages = Jls::Virus::_random_birth_ages( $R, $m ); 

die unless @$birth_ages == 15;
for (my $i = 0; $i != @$birth_ages; $i++) { die unless $birth_ages->[ $i ] == 1.0; }

$R->run( 'set.seed( 31415 )' );

# deterministic $mu = 2.0 because ! defined $sigma (0.0)
$birth_ages = Jls::Virus::_random_birth_ages( $R, $m, $mu ); 

die unless @$birth_ages == 15;
for (my $i = 0; $i != @$birth_ages; $i++) { die unless $birth_ages->[ $i ] == $mu; }

# random gamma( $mu, $sigma ) in ascending order
$birth_ages = Jls::Virus::_random_birth_ages( $R, $m, $mu, $sigma ); 

my $birth_ages_0 =
[
    1.722800,
    1.769459,
    1.771818,
    1.854591,
    1.868233,
    1.918559, # 5
    1.939816,
    1.971192,
    2.031599,
    2.280907,
    2.345833,
    2.447765
];

die unless @$birth_ages == @$birth_ages_0;
for (my $i = 0; $i != @$birth_ages_0; $i++) 
{ 
    die unless $birth_ages->[ $i ] == $birth_ages_0->[ $i ]; 
}

# Tests the Virus object.

$R->run( 'set.seed( 31415 )' );

$m = 10; # Poisson mean number of offspring
$mu = 2.0; # mean time to birth
$sigma = 0.2; # st dev time to birth

die unless ! defined $founder->get_parent();
die unless $founder->get_time() == 0.0;
die unless $founder->stringify() eq "0\tundef\t0";

my $daughters;
my @time0 = ( 3.0, 1.0, 2.0 ); 

$daughters = $founder->daughters( \@time0 );
@time0 = sort { $a <=> $b } @time0;

die unless @$daughters == 3;
for (my $i = 0; $i != @$daughters; $i++) 
{
    $daughters->[ $i ]->set_samples( $i );
    die unless $daughters->[ $i ]->{ 'time' } == $time0[ $i ]; 
    die unless $daughters->[ $i ]->stringify() eq $time0[ $i ] . "\t0\t" . $i; 
}

$daughters = $founder->random_daughters( $R, $m ); 

die unless @$daughters == 15;
for (my $i = 0; $i != @$daughters; $i++) 
{ 
    die unless $daughters->[ $i ]->{ 'time' } == 1.0; 
}

$R->run( 'set.seed( 31415 )' );

$daughters = $founder->random_daughters( $R, $m, $mu ); 

die unless @$daughters == 15;
for (my $i = 0; $i != @$daughters; $i++) 
{ 
    die unless $daughters->[ $i ]->{ 'time' } == $mu; 
}

$daughters = $founder->random_daughters( $R, $m, $mu, $sigma ); 

# Prints Virus.

die unless @$daughters == @$birth_ages_0;
for (my $i = 0; $i != @$birth_ages_0; $i++) 
{ 
    die unless $daughters->[ $i ]->get_time() == $birth_ages_0->[ $i ]; 
}

my $times = [ 0, 2, 4, 6 ];
my $granddaughters = $daughters->[ 5 ]->daughters( $times );

die unless @$granddaughters == @$times;
for (my $i = 0; $i != @$times; $i++) 
{
    $granddaughters->[ $i ]->set_samples( $i );
    die unless 
        $granddaughters->[ $i ]->get_parent() == 
        $daughters->[ 5 ];
    die unless 
        $granddaughters->[ $i ]->get_time() == 
        $times->[ $i ] + $daughters->[ 5 ]->get_time();
    die unless $granddaughters->[ $i ]->stringify() eq 
        $granddaughters->[ $i ]->get_time() . "\t" . $daughters->[ 5 ]->get_time() . "\t" . $i; 
}

die unless $founder->add_samples(  0 )->get_samples() ==  0;
die unless $founder->add_samples(  3 )->get_samples() ==  3;
die unless $founder->add_samples( -2 )->get_samples() ==  1;

die unless $founder->set_samples(  0 )->get_samples() ==  0;
die unless $founder->set_samples(  3 )->get_samples() ==  3;
die unless $founder->set_samples( -2 )->get_samples() == -2;

# Tests get_samples and add_samples.



1