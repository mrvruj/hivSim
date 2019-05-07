#!/usr/bin/perl

use strict;
use warnings;

package Jls::Virus;

use Exporter;
our @ISA = 'Exporter';
our @EXPORT = 'new daughters get_time random_daughters cmp_time';

# The founder is Virus::new and has 'parent' => undef.
# The founder's descendants are $virus->new and have 'parent' => $virus.

sub new 
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
    my $time_ = shift; # time of birth
    
    $self->{ 'parent' } = ref( $proto_ ) ? $proto_ : undef;
    $time_ = 0.0 unless defined $time_;
    
    $self->{ 'time' } = $time_;
}

sub get_parent { return shift()->{ 'parent' }; }
sub get_time { return shift()->{ 'time' }; }

sub stringify
{
    my $self = shift;
    
    my $s = '';
    
    $s .= $self->get_time();
    
    my $parent = $self->get_parent();
    my $parent_time = defined $parent ? $parent->get_time() : 'undef';
    $s .= "\t" . $parent_time;
    
    return $s;
}
    
# Returns unblessed Virus for use in JSON, e.g.

sub TO_JSON 
{
    return { %{ shift() } };
}

# Compares Virus on 'time'.

sub cmp_time
{
    use Jls::VirusUtil;
    return Jls::VirusUtil::cmp_time( shift, shift );
}

# Returns simulated \@daughters of a Virus, sorted on 'time'.

sub daughters
# $self 
# \@birth_ages_ 
{
    my $self = shift; # mother giving birth
    my $birth_ages_ = shift; # mother's age at births

    my @birth_ages = sort { $a <=> $b } @$birth_ages_;
    my @daughters = ();
    
    foreach my $birth_age (@birth_ages)
    {
        push
        ( 
            @daughters, 
            $self->new( $self->{ 'time' } + $birth_age )
        );
    }
    
    return \@daughters;
}

# Returns simulated \@daughters of a Virus, sorted on 'time'.

sub random_daughters
# $R_ # R object
# $r0_ # Poisson mean number of offspring = reproductive number
# $mu_ # mean time to birth
# $sigma_ # st dev time to birth
{
    my $self = shift; # parent giving birth
    my $R_ = shift; # R object
    my $r0_ = shift; # Poisson mean number of offspring
    my $mu_ = shift; # mean time to birth
    my $sigma_ = shift; # st dev time to birth

    my $birth_ages = _random_birth_ages( $R_, $r0_, $mu_, $sigma_ );

    return $self->daughters( $birth_ages );
}

# Returns in ascending order some simulated birth times from a mother.
# Requires an R session.

sub _random_birth_ages
{
    use Jls::VirusUtil;
    return Jls::VirusUtil::random_birth_ages( shift, shift, shift, shift );
}

sub _poisson_generating_function 
# $x_ # argument
# $r0_ # Poisson mean number of offspring = reproductive number
{ 
    use Jls::VirusUtil;
    return Jls::VirusUtil::poisson_generating_function ( shift, shift );
}

# Returns the extinction probability for a Poisson branching process.

sub extinction_probability
# $r0_ = shift; # Poisson mean number of offspring
# $relative_tolerance_ # (1.0e-06)
{
    use Jls::VirusUtil;
    return Jls::VirusUtil::extinction_probability ( shift, shift );
}
    
# The Harris-Sevastyanov transformation
#     Makes a subcritical branching process immortal.

sub harris_sevastyanov
# $r0_ # Poisson mean number of offspring
# $relative_tolerance_ # (1.0e-06)
{
    use Jls::VirusUtil;
    return Jls::VirusUtil::harris_sevastyanov ( shift, shift );
}
    
1;
