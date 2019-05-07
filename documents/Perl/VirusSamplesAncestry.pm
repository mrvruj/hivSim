#!/usr/bin/perl

use strict;
use warnings;

package Jls::VirusSamplesAncestry;

use Exporter;
our @ISA = 'Exporter';
our @EXPORT = 'new TO_JSON stringify to_population_count add_samples sample_spectrum';

# Stores the parameters for a random ancestry of Virus-es.
# Performs the Harris-Sevastyanov transformation if an immortal process is preferred.

sub new 
# $class_ = ref ($proto_) || $proto_;
# $m_ = shift; # Poisson mean number of offspring
# $mu_ = shift; # (1.0) mean time to birth
# $sigma_ = shift; # (0.0) st dev time to birth
# $is_immortal_ = shift; # (0) ? immortalize Poisson branching process (H-S transform) ?
# $is_with_replacement_ # (0) ? sample with replacement ?
{
    my $proto_ = shift;
    my $class_ = ref ($proto_) || $proto_;
    
    my $r0_ = shift; # Poisson mean number of offspring = reproductive number
    my $mu_ = shift; # (1.0) mean time to birth
    my $sigma_ = shift; # (0.0) st dev time to birth
    my $is_immortal_ = shift; # (0) ? the immortalized Poisson branching process ?
    my $is_with_replacement_ = shift; # (0) ? sample with replacement ?
    
    $mu_ = 1.0 unless defined $mu_;
    $sigma_ = 0.0 unless defined $sigma_;
    $is_immortal_ = 0 unless defined $is_immortal_;
    $is_with_replacement_ = 0 unless defined $is_with_replacement_;

    die "Jls::VirusAncestryBase::new : $r0_ <= 1.0" if $r0_ <= 1.0;
    
    use Jls::VirusUtil;
    my $r0 = $is_immortal_ ? Jls::VirusUtil::harris_sevastyanov( $r0_ ) : $r0_;
            
    my $self = 
    {
        'ancestors' => [],
        'descendants' => [],

        'r0' => $r0,
        'mu' => $mu_,
        'sigma' => $sigma_,
        'is_immortal' => $is_immortal_,
        'is_with_replacement' => $is_with_replacement_,

        'shuffles' => [],
        'sampled' => 0,
        'sample_spectrum' => [],
    };
    
    return bless $self, $class_;
}

sub get_r0 { return shift()->{ 'r0' }; }
sub get_mu { return shift()->{ 'mu' }; }
sub get_sigma { return shift()->{ 'sigma' }; }
sub get_is_immortal { return shift()->{ 'is_immortal' }; }
sub get_is_with_replacement { return shift()->{ 'is_with_replacement' }; }
sub get_ancestors { return shift()->{ 'ancestors' }; }
sub get_descendants { return shift()->{ 'descendants' }; }
sub get_sampled { return shift()->{ 'sampled' }; }
sub get_sample_spectrum { return shift()->{ 'sample_spectrum' }; }

sub get_last_birth_time {  return shift()->{ 'ancestors' }->[ -1 ]->get_time(); }
sub get_next_birth_time {  return shift()->{ 'descendants' }->[ 0 ]->get_time(); }

sub TO_JSON { return { %{ shift() } }; } # Returns unblessed VirusAncestry for use in JSON, e.g.

sub stringify
{
    my $self = shift;    

    my $s = '';
    
    foreach my $virus (@{ $self->get_ancestors() })
    {
        if ($s ne '') { $s .= "\n"; }
        $s .= $virus->stringify();
    }

    $s .= "\n";

    foreach my $virus (@{ $self->get_descendants() })
    {
        $s .= "\n"; 
        $s .= $virus->stringify();
    }
    
    return $s;
}
    
# Simulates the next generation of ancestor and descendant Virus-es
#     until \@descendants exceed $population_count_, and
# Prepares also for sampling with a random shuffle of the \@descendants.

sub to_population_count 
# $self
# $R_ # R object
# $population_count_ # terminate simulation when descendants reach this count 
{
    my $self = shift;
    my $R_ = shift;
    my $population_count_ = shift;

    if (@{ $self->get_descendants() } < $population_count_)
    {
        $self->{ 'sampled' } = 0;
    }
    
    while (@{ $self->get_descendants() } < $population_count_) 
    {
        $self->_next_generation( $R_ );
    }

    $self->_init_shuffles( $R_ );

    return $self;
}    

# Returns an unsorted incremental sample of indices from $descendants Virus-es.

sub add_samples
# $self
# $samples_ # number of samples added
{
    my $self = shift;
    my $samples_ = shift; # number of samples added

    die unless 0 <= $samples_;

    if ($samples_ == 0) { return $self; }
    
    my $sampled = $self->{ 'sampled' };

    $self->{ 'sampled' } += $samples_; # live viruses added to samples
    die unless $self->{ 'sampled' } <= @{ $self->{ 'descendants' } };
    
    $self->_update( $sampled );
    $self->_calc_sample_spectrum();

    return $self;
}

# Calculates the sample spectrum of ancestors (excluding unsampled ancestors).

sub _calc_sample_spectrum
# $self
{
    my $self = shift;

    my $sample_spectrum = []; 
    
    my $ancestors = $self->{ 'ancestors' };
    my $descendants = $self->{ 'descendants' };

    foreach my $ancestor (@$ancestors)
    {
        if (defined $ancestor->get_parent()) { $sample_spectrum->[ $ancestor->get_samples() ]++; }
    }

    foreach my $descendant (@$descendants)
    {
        $sample_spectrum->[ $descendant->get_samples() ]++; 
    }
    
    for (my $i = 0; $i <= $self->get_sampled(); $i++)
    {
        $sample_spectrum->[ $i ] = 0 unless defined $sample_spectrum->[ $i ];
    }
    
    $self->{ 'sample_spectrum' } = $sample_spectrum;
    
    return $self;
}

# Simulates the next generation of ancestor and descendant Virus-es.
# Initializes with $founder in \@descendants.
# Updates by 
#    (1) moving earliest descendant into the \@ancestors, and 
#    (2) adding daughters to \@descendants.

sub _next_generation 
# $self
# $R_ # R object
{
    my $self = shift;
    my $R_ = shift;

    my $ancestors = $self->{ 'ancestors' };
    my $descendants = $self->{ 'descendants' };

    my $r0 = $self->{ 'r0' };
    my $mu = $self->{ 'mu' };
    my $sigma = $self->{ 'sigma' };
    my $is_immortal = $self->{ 'is_immortal' };

    my $founder;
    
    use Jls::VirusSamples;
    
    if (! @$ancestors) 
    {            
        push( @$descendants, Jls::VirusSamples->new() ); # founder
    }

    my $current = shift @$descendants;
    
    if (! $is_immortal && ! defined $current) # Renews extinct process.
    {
        $self->{ 'ancestors' } = [];
        $self->{ 'descendants' } = [];
        
        push( @$descendants, Jls::VirusSamples->new() ); # founder
        $current = shift @$descendants;
    }
            
    my $daughters = $current->random_daughters( $R_, $r0, $mu, $sigma );
    
    while ($is_immortal && ! @$daughters) 
    {
        $daughters = $current->random_daughters( $R_, $r0, $mu, $sigma );
    }
    
    push( @$ancestors, $current );

    if (! @$daughters) { return $self; }
   
    if (! @$descendants)
    {
        $self->{ 'descendants' } = $daughters;
        return $self;
    }
    
    my $i;
    
    for ($i = @$descendants; 
         $i > 0 && 
             $daughters->[ 0 ]->get_time() < 
             $self->{ 'descendants' }->[ $i - 1 ]->get_time();
         $i-- ) {}
    
    # $i - 1 can be kept in the @$descendants.
    
    if ($i < @$descendants)
    {
        my @reorder = splice( @$descendants, $i );

        use Jls::Array;
        my $tail = 
            Jls::Array::merge( \@reorder, $daughters, \&Jls::VirusUtil::cmp_time );
            
        push( @$descendants,  @$tail )
    }
    else { push( @$descendants,  @$daughters ); }

    $self->{ 'ancestors' } = $ancestors;
    $self->{ 'descendants' } = $descendants;
    
    return $self;
}

# Initializes a the random shuffle of length $self->get_descendants().

sub _init_shuffles
# $self
# $R_ # R object
{
    my $self = shift;
    my $R_ = shift; 

    my @shuffles = ();
    
    my $descendants = $self->get_descendants();

    if ($self->get_is_with_replacement()) # Initializes samples with replacement.
    {
        foreach my $j ( 0...(@$descendants - 1) )
        {
            push( @shuffles, int( rand( @$descendants ) ) );
        }
    }        
    else # Initializes samples without replacement in shuffled order.
    {
        use List::Util qw(shuffle);
        @shuffles = shuffle( 0...(@$descendants - 1) );
    }

    $self->{ 'shuffles' } = \@shuffles;

    return $self;
}

# Returns Virus ancestors_ with count of sampled $descendants.

sub _update
{
    my $self = shift;
    my $sampled_ = shift; 
    
    my $ancestors = $self->{ 'ancestors' };
    my $descendants = $self->{ 'descendants' };
    my $shuffles = $self->{ 'shuffles' };

    for (my $i = $sampled_; $i != $self->{ 'sampled' }; $i++)
    {
        my $sample = $shuffles->[ $i ];
        my $ancestor = $descendants->[ $sample ];
        
        do 
        {
            $ancestor->add_samples( 1 );
            $ancestor = $ancestor->{ 'parent' };
        }
        while (defined $ancestor);
    }
}

1;
