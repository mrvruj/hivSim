#!/usr/bin/perl

use strict;
use warnings;

package Jls::VirusUtil;

use Exporter;
our @ISA = 'Exporter';
our @EXPORT = 'cmp_time random_birth_ages poisson_generating_function extinction_probability harris_sevastyanov';

# Compares Virus on 'time'.

sub cmp_time
{
    my $a_ = shift;
    my $b_ = shift;
    
    return $a_->{ 'time' } <=> $b_->{ 'time' }
}

# Returns in ascending order some simulated birth times from a mother.
# Requires an R session.

sub random_birth_ages
{
    my $R_ = shift; # R object
    my $r0_ = shift; # Poisson mean number of offspring
    my $mu_ = shift; # (1.0) mean time to birth
    my $sigma_ = shift; # (0.0) st dev time to birth

    $mu_ = 1.0 unless defined $mu_;
    $sigma_ = 0.0 unless defined $sigma_;

    use Jls::RIO;
        
# Poisson Random Variate

my $cmd_s = <<EOF;
    array <- rpois( n = 1, lambda = $r0_ ) 
    array
EOF

    my $n = Jls::RIO::o_perl( $R_->run( $cmd_s ) ); # number of offspring

    if ($sigma_ == 0.0) 
    { 
        my @offspring = ($mu_) x $n; 
        return \@offspring;
    }
    
    my $shape = ($mu_ / $sigma_) * ($mu_ / $sigma_);
    my $lambda = $mu_ / $sigma_ / $sigma_;
    
# Gamma Random Variates

$cmd_s = <<EOF;
    array <- rgamma( n = $n, shape = $shape, rate = $lambda ) 
    array
EOF

    my $birth_ages = Jls::RIO::o_perl( $R_->run( $cmd_s ) ); # offspring birthtimes
    
    if ($n == 0) { return []; }
    elsif ($n == 1) { return [ $birth_ages ]; }
    
    my @sorted_birth_ages = sort { $a <=> $b } @$birth_ages;
    
    return \@sorted_birth_ages;
}

# Returns the Poisson generating function.

sub poisson_generating_function 
# $x_ # argument
# $r0_ # Poisson mean number of offspring = reproductive number
{ 
    my $x_ = shift;
    my $r0_ = shift; # Poisson mean number of offspring

    return exp( $r0_ * ($x_ - 1.0) ); 
}

# Returns the extinction probability for a Poisson branching process.

sub extinction_probability
# $r0_ = shift; # Poisson mean number of offspring
# $relative_tolerance_ # (1.0e-06)
{
    my $r0_ = shift; # Poisson mean number of offspring
    my $relative_tolerance_ = shift;

    my $RELATIVE_TOLERANCE = 1.0e-06;
    $relative_tolerance_ = $RELATIVE_TOLERANCE
        unless defined $relative_tolerance_;

    die "Jls::Virus::extinction_probability : $r0_ <= 1.0" if $r0_ <= 1.0;
    
    my $q0 = 0.0;
    my $q1;
    
    while (1)
    {
        $q1 = poisson_generating_function( $q0, $r0_ );
        last if (abs( 1.0 - $q0 / $q1 ) < $RELATIVE_TOLERANCE);
        $q0 = $q1;
    }
    
    return $q1;
}
    
# The Harris-Sevastyanov transformation
#     Makes a subcritical branching process immortal.

sub harris_sevastyanov
# $r0_ # Poisson mean number of offspring
# $relative_tolerance_ # (1.0e-06)
{
    my $r0_ = shift; # Poisson mean number of offspring
    my $relative_tolerance_ = shift;

    return (1.0 - extinction_probability( $r0_, $relative_tolerance_ )) * $r0_;
}
    
1;
