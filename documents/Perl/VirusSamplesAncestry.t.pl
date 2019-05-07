#!/usr/bin/perl

use strict;
use warnings;

BEGIN
{
    unshift @INC, '../modules';
}

use Jls::VirusSamplesAncestry;

use Statistics::R;

my $R = Statistics::R->new();

my $population_count = 10;
my $samples = 5; # Tests sampling.

# input

my $r0 = 2; # Poisson mean number of offspring
my $mu = 10.0; # mean time to birth
my $sigma = 0.1; # st dev time to birth
my $immortal = 0;
my $is_with_replacement = 1;

# test

$R->run( 'set.seed( 31415 )' ); # R random seed
srand( 31415 ); # Perl random seed

# Checks simulation of Virus ancestry, sampling with replacement.

my $virus_ancestry = Jls::VirusSamplesAncestry->new( $r0, $mu, $sigma, 0, $is_with_replacement );

my $shuffles0 = '1 9 2 8 10 10 2 6 1 7 4';

my $stringify0 = <<'END';
0	0	undef
0	9.845249	0
0	9.861224	0
0	9.886034	0
0	9.957134	0
0	9.974201	0
0	19.709181	9.886034

0	19.788834	9.861224
0	19.798822	9.861224
0	19.82462	9.861224
0	19.947186	9.957134
0	19.982708	9.861224
0	20.061111	9.974201
0	20.079493	9.861224
0	20.079911	9.974201
0	29.623042	19.709181
0	29.669188	19.709181
0	29.674494	19.709181
END

my $sample_spectrum0 = '';

$virus_ancestry->to_population_count( $R, $population_count );
die unless "@{ $virus_ancestry->{ 'shuffles' } }" eq $shuffles0;
die unless $virus_ancestry->stringify() . "\n" eq $stringify0;
die unless $virus_ancestry->get_last_birth_time() == 19.709181;
die unless $virus_ancestry->get_next_birth_time() == 19.788834;
die unless $virus_ancestry->get_sampled() == 0;
die unless "@{ $virus_ancestry->get_sample_spectrum() }" eq $sample_spectrum0;

$stringify0 = <<'END';
5	0	undef
0	9.845249	0
2	9.861224	0
3	9.886034	0
0	9.957134	0
0	9.974201	0
3	19.709181	9.886034

0	19.788834	9.861224
1	19.798822	9.861224
1	19.82462	9.861224
0	19.947186	9.957134
0	19.982708	9.861224
0	20.061111	9.974201
0	20.079493	9.861224
0	20.079911	9.974201
1	29.623042	19.709181
1	29.669188	19.709181
1	29.674494	19.709181
END

$sample_spectrum0 = '9 5 1 2 0 0';

$virus_ancestry->add_samples( $samples );
die unless "@{ $virus_ancestry->{ 'shuffles' } }" eq $shuffles0;
die unless $virus_ancestry->stringify() . "\n" eq $stringify0;
die unless $virus_ancestry->get_last_birth_time() == 19.709181;
die unless $virus_ancestry->get_next_birth_time() == 19.788834;
die unless $virus_ancestry->get_sampled() == $samples;
die unless "@{ $virus_ancestry->get_sample_spectrum() }" eq $sample_spectrum0;

$stringify0 = <<'END';
10	0	undef
0	9.845249	0
5	9.861224	0
4	9.886034	0
0	9.957134	0
1	9.974201	0
4	19.709181	9.886034

0	19.788834	9.861224
2	19.798822	9.861224
2	19.82462	9.861224
0	19.947186	9.957134
0	19.982708	9.861224
0	20.061111	9.974201
1	20.079493	9.861224
1	20.079911	9.974201
1	29.623042	19.709181
1	29.669188	19.709181
2	29.674494	19.709181
END

$sample_spectrum0 = '6 5 3 0 2 1 0 0 0 0 0';

$virus_ancestry->add_samples( $samples );
die unless "@{ $virus_ancestry->{ 'shuffles' } }" eq $shuffles0;
die unless $virus_ancestry->stringify() . "\n" eq $stringify0;
die unless $virus_ancestry->get_last_birth_time() == 19.709181;
die unless $virus_ancestry->get_next_birth_time() == 19.788834;
die unless $virus_ancestry->get_sampled() == 2 * $samples;
die unless "@{ $virus_ancestry->get_sample_spectrum() }" eq $sample_spectrum0;

exit( 0 );

1
