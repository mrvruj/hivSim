#!/usr/bin/perl

use strict;
use warnings;

package Jls::VirusSamples;

use Exporter;
our @ISA = ( 'Exporter', 'Jls::Virus' );
our @EXPORT = 'new initialize stringify get_samples set_samples add_samples';

use base qw( Jls::Virus );

# Counts sampled descendants of the base class Virus.

sub new # VirusSamples
{
    my $proto_ = shift;
    my $class_ = ref( $proto_ ) || $proto_;
    my $self = bless {}, $class_;
    
    $self->initialize( $proto_, @_ );
    return $self;
}

sub initialize
{
    my $self = shift;
    my $proto_ = shift;

    $self->SUPER::initialize( $proto_, @_ );
    $self->{ 'samples' } = 0;

    return $self;
}

sub get_samples { return shift()->{ 'samples' }; }

sub stringify
{
    my $self = shift;

    my $s = $self->{ 'samples' };
    $s .= "\t" . $self->SUPER::stringify();

    return $s;
}

sub set_samples
# $samples_
{
    my $self = shift;
    my $samples_ = shift;

    $self->{ 'samples' } = $samples_;

    return $self;
}

# Adds to the count of VirusSamples. 
#     Negative increments are permitted, but module dies on negative sample totals.

sub add_samples
# $samples_
{
    my $self = shift;
    my $samples_ = shift; 
    
    $self->{ 'samples' } += $samples_;
    die 'VirusSamples->add_samples : samples < 0' if ($self->{ 'samples' } < 0);
    
    return $self;
}
   
1;
